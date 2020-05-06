#include "nori/octree.h"
#include <Eigen/Geometry>
#include <atomic>
#include <chrono>

namespace {
/**
 * Merge two matrices
 *
 * */
template <typename M>
inline void mergeMatrix(M& A, const M& B) {
  if (B.size() == 0) return;
  M C(A.rows(), A.cols() + B.cols());
  C << A, B;
  A = C;
}
}  // namespace

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

OctreeNode* buildNode(Octree* tree, uint32_t depth, const BoundingBox3f& bbox,
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

  const auto& pos = tree->getVertexPositions();
  const auto& face = tree->getIndices();
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
  BBoxIts(uint32_t i, float n) : index{i}, near{n} {}

  bool operator<(const BBoxIts& b) const { return near < b.near; }

  uint32_t index;
  float near;
};

#define USE_SORT true

void OctreeNode::traverse(Ray3f& ray, Intersection& its, uint32_t& f,
                          bool& found, bool shadowRay) const {
  ++numIter;

  const auto& ids = mRoot->getIdentities();
  for (const auto& idx : mIndices) {
    float u, v, t;
    auto objId = ids(idx);
    auto mesh = mRoot->m_mesh_map.at(objId);
    if (mRoot->rayIntersect(idx, ray, u, v, t)) {
      if (shadowRay) {
        found = true;
        return;
      }
      if (ray.maxt < t) continue;
      ray.maxt = its.t = t;
      its.uv = Point2f(u, v);
      its.mesh = mesh;
      f = idx - mRoot->m_base_index[objId];
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
      if (its.t < it.near) continue;
      mChildren[it.index]->traverse(ray, its, f, found, shadowRay);
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

Octree::Octree(Mesh* mesh, uint32_t id)
    : m_bbox{mesh->getBoundingBox()},
      m_V{mesh->getVertexPositions()},
      m_F{mesh->getIndices()} {
  m_ID.resize(1, m_F.cols());
  m_ID.array() = id;
  m_mesh_map[id] = mesh;
  m_base_index[id] = 0;
}

Octree::~Octree() { cout << "numIter = " << numIter << "\n"; }

void Octree::build() {
  auto bbox = m_bbox;
  const float margin = 0.1F;
  const auto delta = Point3f(margin, margin, margin);
  auto augBBox = BoundingBox3f(bbox.min - delta, bbox.max + delta);

  std::vector<uint32_t> indices(m_F.cols());
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

void Octree::addMesh(Mesh* newMesh, uint32_t id) {
  m_mesh_map[id] = newMesh;
  m_base_index[id] = m_F.cols();

  auto c1 = m_V.cols();

  mergeMatrix(m_V, newMesh->getVertexPositions());

  MatrixXu F2 = newMesh->getIndices().array() + (MatrixXu::Scalar)c1;
  mergeMatrix(m_F, F2);

  MatrixXu ids(1, newMesh->getIndices().cols());
  ids.array() = id;
  mergeMatrix(m_ID, ids);

  m_bbox.expandBy(newMesh->getBoundingBox());
}

void Octree::printState() {
  NodeVisitor v;
  mRootNode->accept(v);
  v.print();
}

bool Octree::rayIntersect(uint32_t index, const Ray3f& ray, float& u, float& v,
                          float& t) const {
  uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
  const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

  /* Find vectors for two edges sharing v[0] */
  Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

  /* Begin calculating determinant - also used to calculate U parameter */
  Vector3f pvec = ray.d.cross(edge2);

  /* If determinant is near zero, ray lies in plane of triangle */
  float det = edge1.dot(pvec);

  if (det > -1e-8f && det < 1e-8f) return false;
  float inv_det = 1.0f / det;

  /* Calculate distance from v[0] to ray origin */
  Vector3f tvec = ray.o - p0;

  /* Calculate U parameter and test bounds */
  u = tvec.dot(pvec) * inv_det;
  if (u < 0.0 || u > 1.0) return false;

  /* Prepare to test V parameter */
  Vector3f qvec = tvec.cross(edge1);

  /* Calculate V parameter and test bounds */
  v = ray.d.dot(qvec) * inv_det;
  if (v < 0.0 || u + v > 1.0) return false;

  /* Ray intersects triangle -> compute t */
  t = edge2.dot(qvec) * inv_det;

  return t >= ray.mint && t <= ray.maxt;
}

NORI_NAMESPACE_END
