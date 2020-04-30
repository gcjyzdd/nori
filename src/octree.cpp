#include "nori/octree.h"

NORI_NAMESPACE_BEGIN

inline BoundingBox3f createBBox(const Point3f& p1, const Point3f& p2,
                                const Point3f& p3) {
  BoundingBox3f bb(p1, p2);
  bb.expandBy(p3);
  return bb;
}

inline std::vector<BoundingBox3f> splitBBox(const BoundingBox3f& box) {
  auto c = box.getCenter();
  auto e = box.getExtents() / 2.F;
  auto x = e(0);
  auto y = e(1);
  auto z = e(2);
  return {BoundingBox3f(c, box.max),
          BoundingBox3f(c + Point3f(-x, 0, 0), c + Point3f(0, y, z)),
          BoundingBox3f(c + Point3f(-x, -y, 0), c + Point3f(0, 0, z)),
          BoundingBox3f(c + Point3f(0, -y, 0), c + Point3f(x, 0, z)),
          BoundingBox3f(c + Point3f(0, 0, -z), c + Point3f(x, y, 0)),
          BoundingBox3f(c + Point3f(-x, 0, -z), c + Point3f(0, y, 0)),
          BoundingBox3f(c + Point3f(-x, -y, -z), c + Point3f(0, 0, 0)),
          BoundingBox3f(c + Point3f(0, -y, -z), c + Point3f(x, 0, 0))};
}

static OctreeNode* buildNode(Octree* tree, const BoundingBox3f& bbox,
                             const std::vector<uint32_t>& indices) {
  if (indices.size() == 0) return nullptr;

  if (indices.size() <= N_LEAF) {
    auto ret = new OctreeNode(tree);
    ret->mBbox = bbox;
    for (auto idx : indices) {
      ret->mIndices.push_back(idx);
    }

    return ret;
  }

  std::vector<std::vector<uint32_t>> triList(N_LEAF);

  const auto& pos = tree->getMeshPtr()->getVertexPositions();
  const auto& face = tree->getMeshPtr()->getIndices();
  auto boxes = splitBBox(bbox);
  // detect if a triangle overlaps a bbox
  for (auto idx : indices) {
    auto box = createBBox(pos.col(face(0, idx)), pos.col(face(1, idx)),
                          pos.col(face(2, idx)));

    for (size_t j = 0; j < N_LEAF; ++j) {
      if (boxes[j].overlaps(box)) triList[j].push_back(idx);
    }
  }

  OctreeNode* node = new OctreeNode(tree);
  node->mBbox = bbox;
  for (size_t i = 0; i < N_LEAF; ++i)
    node->mChildren[i] = buildNode(tree, boxes[i], triList[i]);

  return node;
}

struct BBoxIts {
  BBoxIts(uint32_t i, float f) : index{i}, near{f} {}

  bool operator<(const BBoxIts& b) const { return near < b.near; }

  uint32_t index;
  float near;
};

void OctreeNode::traverse(Ray3f& ray, Intersection& its, uint32_t& f,
                          bool& found, bool shadowRay) const {
  auto mesh = mRoot->getMeshPtr();

  for (const auto& idx : mIndices) {
    float u, v, t;
    if (mesh->rayIntersect(idx, ray, u, v, t)) {
      if (shadowRay) return;
      ray.maxt = its.t = t;
      its.uv = Point2f(u, v);
      its.mesh = mesh;
      f = idx;
      found = true;
    }
  }

  if (mIndices.size() == 0) {
    float near, far;
    std::vector<BBoxIts> dist;
    dist.reserve(N_LEAF);
    for (uint32_t i = 0; i < N_LEAF; ++i) {
      if (mChildren[i] && mChildren[i]->mBbox.rayIntersect(ray, near, far)) {
        dist.emplace_back(i, near);
      }
    }
    if (dist.size() > 1) std::sort(dist.begin(), dist.end());
    for (const auto& it : dist) {
      mChildren[it.index]->traverse(ray, its, f, found, shadowRay);
      if (found) return;
    }
  }
}

Octree::Octree(Mesh* mesh) : mMesh{mesh} {}

void Octree::build() {
  auto bbox = mMesh->getBoundingBox();
  std::vector<uint32_t> indices(mMesh->getIndices().cols());
  for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;
  mRootNode = buildNode(this, bbox, indices);
}

bool Octree::rayIntersect(Ray3f& ray, Intersection& its, uint32_t& f,
                          bool shadowRay) {
  bool found = false;
  if (mRootNode->mBbox.rayIntersect(ray))
    mRootNode->traverse(ray, its, f, found, shadowRay);
  return found;
}

NORI_NAMESPACE_END
