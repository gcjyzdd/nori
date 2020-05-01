#pragma once
#include <nori/mesh.h>
#include <vector>

NORI_NAMESPACE_BEGIN

static const unsigned int NUM_NODE = 8;
static const unsigned int N_LEAF = 8;  // set it to 4 causes issues! Investigate
static const unsigned int MAX_TREE_DEPTH = 8;

class Octree;
class OctreeNode {
 public:
  OctreeNode(Octree* tree, uint32_t d = 1)
      : mRoot{tree}, mDepth{d}, id{mCounter} {
    ++mCounter;
  }

  void traverse(Ray3f& ray, Intersection& its, uint32_t& f, bool& found,
                bool shadowRay) const;

 public:
  Octree* mRoot{nullptr};
  OctreeNode* mChildren[NUM_NODE];
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

class Octree {
 public:
  Octree(Mesh* mesh);
  ~Octree();
  void build();

  bool rayIntersect(Ray3f& ray, Intersection& its, uint32_t& f,
                    bool shadowRay) const;

  Mesh* getMeshPtr() { return mMesh; }

 private:
  Mesh* mMesh;
  OctreeNode* mRootNode;
};

NORI_NAMESPACE_END
