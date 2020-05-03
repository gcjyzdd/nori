#include "nori/octree.h"
#include <atomic>
#include <chrono>

NORI_NAMESPACE_BEGIN

std::atomic<int> numIter{0};
uint32_t OctreeNode::mCounter = 0;

inline BoundingBox3f createBBox(const Point3f& p1, const Point3f& p2,
                                const Point3f& p3) {
  BoundingBox3f bb(p1);
  bb.expandBy(p2);
  bb.expandBy(p3);
  return bb;
}

inline std::array<BoundingBox3f, NUM_NODE> splitBBox(const BoundingBox3f& box) {
  Point3f c = box.getCenter();
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

static OctreeNode* buildNode(Octree* tree, uint32_t depth,
                             const BoundingBox3f& bbox,
                             const std::vector<uint32_t>& indices) {
  if (indices.size() == 0) return nullptr;

  if (indices.size() <= N_LEAF || depth >= MAX_TREE_DEPTH) {
    auto ret = new OctreeNode(tree, depth);
    ret->mBbox = bbox;
    for (auto idx : indices) {
      ret->mIndices.push_back(idx);
    }

    return ret;
  }

  std::vector<std::vector<uint32_t>> triList(NUM_NODE);

  const auto& pos = tree->getMeshPtr()->getVertexPositions();
  const auto& face = tree->getMeshPtr()->getIndices();
  auto boxes = splitBBox(bbox);
  // detect if a triangle overlaps a bbox
  for (auto idx : indices) {
    auto box = createBBox(pos.col(face(0, idx)), pos.col(face(1, idx)),
                          pos.col(face(2, idx)));

    for (size_t j = 0; j < NUM_NODE; ++j) {
      if (boxes[j].overlaps(box)) triList[j].push_back(idx);
    }
  }

  OctreeNode* node = new OctreeNode(tree, depth);
  node->mBbox = bbox;
  for (size_t i = 0; i < NUM_NODE; ++i)
    node->mChildren[i] = buildNode(tree, depth + 1, boxes[i], triList[i]);

  return node;
}

struct BBoxIts {
  BBoxIts(uint32_t i, float f) : index{i}, near{f} {}

  bool operator<(const BBoxIts& b) const { return near < b.near; }

  uint32_t index;
  float near;
};

#define USE_SORT true

void OctreeNode::traverse(Ray3f& ray, Intersection& its, uint32_t& f,
                          bool& found, bool shadowRay) const {
  ++numIter;
  auto mesh = mRoot->getMeshPtr();

  for (const auto& idx : mIndices) {
    float u, v, t;
    if (mesh->rayIntersect(idx, ray, u, v, t)) {
      if (shadowRay) {
        found = true;
        return;
      }
      if (ray.maxt < t) continue;
      ray.maxt = its.t = t;
      its.uv = Point2f(u, v);
      its.mesh = mesh;
      f = idx;
      found = true;
    }
  }

  if (mIndices.size() == 0) {
#ifdef USE_SORT
    float near, far;
    std::vector<BBoxIts> dist;
    dist.reserve(NUM_NODE);
    for (uint32_t i = 0; i < NUM_NODE; ++i) {
      if (mChildren[i] && mChildren[i]->mBbox.rayIntersect(ray, near, far)) {
        dist.emplace_back(i, near);
      }
    }
    if (dist.size() > 1) std::sort(dist.begin(), dist.end());
    for (const auto& it : dist) {
      mChildren[it.index]->traverse(ray, its, f, found, shadowRay);
      if (found) return;
    }
#else
    for (auto c : mChildren) {
      if (c && c->mBbox.rayIntersect(ray))
        c->traverse(ray, its, f, found, shadowRay);
    }
#endif
  }
}

void OctreeNode::accept(NodeVisitor& visitor) const {
  visitor.accumulateInfo(*this);
  for (auto child : mChildren) {
    if (child) child->accept(visitor);
  }
}

void NodeVisitor::accumulateInfo(const OctreeNode& node) {
  if (node.mDepth > depth) depth = node.mDepth;
  if (node.mIndices.size() == 0) {
    innerNodes += 1;
  } else {
    leafNodes += 1;
    numTri += node.mIndices.size();
  }
}

void NodeVisitor::print() {
  cout << "Max depth of nodes: " << depth << "\n";
  cout << "Number of interior nodes: " << innerNodes << "\n";
  cout << "Number of leaf nodes: " << leafNodes << "\n";
  cout << "Average number of triangles per leaf node: "
       << (float)numTri / leafNodes << "\n";
}

Octree::Octree(Mesh* mesh) : mMesh{mesh} {}

Octree::~Octree() { cout << "numIter = " << numIter << "\n"; }

void Octree::build() {
  auto bbox = mMesh->getBoundingBox();
  const float margin = 0.1F;
  const auto delta = Point3f(margin, margin, margin);
  auto augBBox = BoundingBox3f(bbox.min - delta, bbox.max + delta);

  std::vector<uint32_t> indices(mMesh->getIndices().cols());
  for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;

  auto begin = std::chrono::steady_clock::now();

  mRootNode = buildNode(this, 0, augBBox, indices);

  auto end = std::chrono::steady_clock::now();

  std::cout << "Construction time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     begin)
                   .count()
            << "[ms]" << std::endl;
  printState();
}

bool Octree::rayIntersect(Ray3f& ray, Intersection& its, uint32_t& f,
                          bool shadowRay) const {
  bool found = false;
  if (mRootNode->mBbox.rayIntersect(ray))
    mRootNode->traverse(ray, its, f, found, shadowRay);
  return found;
}

void Octree::printState() {
  NodeVisitor v;
  mRootNode->accept(v);
  v.print();
}

NORI_NAMESPACE_END
