#pragma once
#include <nori/mesh.h>
#include <array>
#include <unordered_map>
#include <vector>

NORI_NAMESPACE_BEGIN

static const unsigned int NUM_NODE = 8;
static const unsigned int N_LEAF = 8;  // set it to 4 causes issues! Investigate
static const unsigned int MAX_TREE_DEPTH = 10;

class Octree;
class NodeVisitor;
class OctreeNode {
 public:
  OctreeNode(Octree* tree, uint32_t d = 1)
      : mRoot{tree}, mDepth{d}, id{mCounter} {
    ++mCounter;
    memset(&mChildren[0], 0, NUM_NODE * sizeof(OctreeNode*));
  }

  void traverse(Ray3f& ray, Intersection& its, uint32_t& f, bool& found,
                bool shadowRay) const;

  void accept(NodeVisitor& visitor) const;

 public:
  Octree* mRoot{nullptr};
  std::array<OctreeNode*, NUM_NODE> mChildren;
  // uint32_t mIndices[N_LEAF];
  std::vector<uint32_t> mIndices;

  // define bbox of the node
  // uint32_t mNodeIdx;  // calculate bbox based on index
  BoundingBox3f mBbox;

  uint32_t mDepth{0};
  uint32_t id{0};

  static uint32_t mCounter;
};

static OctreeNode* buildNode(Octree* tree, uint32_t depth,
                             const BoundingBox3f& bbox,
                             const std::vector<uint32_t>& indices);

class NodeVisitor {
 public:
  NodeVisitor() {}

  void accumulateInfo(const OctreeNode& node);

  void print();

  uint32_t depth = 0;
  uint32_t innerNodes = 0;
  uint32_t leafNodes = 0;
  uint32_t numTri = 0;
};

class Octree {
 public:
  Octree(Mesh* mesh, uint32_t id);
  ~Octree();
  void build();

  bool rayIntersect(Ray3f& ray, Intersection& its, uint32_t& f,
                    bool shadowRay) const;

  void addMesh(Mesh* mesh, uint32_t id);

  const BoundingBox3f& getBoundingBox() const { return m_bbox; }

  /// Return a pointer to the vertex positions
  const MatrixXf& getVertexPositions() const { return m_V; }

  /// Return a pointer to the vertex normals (or \c nullptr if there are none)
  const MatrixXf& getVertexNormals() const { return m_N; }

  /// Return a pointer to the texture coordinates (or \c nullptr if there are
  /// none)
  const MatrixXf& getVertexTexCoords() const { return m_UV; }

  /// Return a pointer to the triangle vertex index list
  const MatrixXu& getIndices() const { return m_F; }

  const MatrixXu& getIdentities() const { return m_ID; }

  void printState();

  bool rayIntersect(uint32_t index, const Ray3f& ray, float& u, float& v,
                    float& t) const;

  std::unordered_map<uint32_t, Mesh*> m_mesh_map;
  std::unordered_map<uint32_t, uint32_t> m_base_index;

 private:
  OctreeNode* mRootNode;

  BoundingBox3f m_bbox;  ///< Bounding box of the mesh
  MatrixXf m_V;          ///< Vertex positions
  MatrixXf m_N;          ///< Vertex normals
  MatrixXf m_UV;         ///< Vertex texture coordinates
  MatrixXu m_F;          ///< Faces
  MatrixXu m_ID;         ///< Face ID's
};

NORI_NAMESPACE_END
