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
    for (int i = 0; i < indices.size(); ++i) {
      ret->mIndices[i] = indices[i];
    }

    return ret;
  }

  std::vector<std::vector<uint32_t>> triList(N_LEAF);

  const auto& pos = tree->getMeshPtr()->getVertexPositions();
  const auto& face = tree->getMeshPtr()->getIndices();
  // detect if a triangle overlaps a bbox
  for (size_t i = 0; i < indices.size(); ++i) {
    auto box = createBBox(pos.col(face(0, i)), pos.col(face(1, i)),
                          pos.col(face(2, i)));
    auto boxes = splitBBox(bbox);
    for (size_t j = 0; j < N_LEAF; ++j) {
      if (boxes[j].overlaps(box)) triList[j].push_back(indices[j]);
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

void Octree::rayIntersect(const Ray3f& ray, float& u, float& v, float& t) {}

NORI_NAMESPACE_END
