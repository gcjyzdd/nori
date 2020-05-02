#pragma once
#include <nori/mesh.h>
#include <array>
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
  Octree(Mesh* mesh);
  ~Octree();
  void build();

  bool rayIntersect(Ray3f& ray, Intersection& its, uint32_t& f,
                    bool shadowRay) const;

  Mesh* getMeshPtr() { return mMesh; }

  void printState();

 private:
  Mesh* mMesh;
  OctreeNode* mRootNode;
};

NORI_NAMESPACE_END
