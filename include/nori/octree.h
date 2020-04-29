#pragma once
#include <nori/mesh.h>
#include <vector>

NORI_NAMESPACE_BEGIN

static const unsigned int NUM_NODE = 8;
static const unsigned int N_LEAF = 8;

class Octree;
class OctreeNode {
 public:
  OctreeNode(Octree* tree) : mRoot{tree} {}

 public:
  Octree* mRoot;
  OctreeNode* mChildren[NUM_NODE];
  uint32_t mIndices[N_LEAF];

  // define bbox of the node
  // uint32_t mNodeIdx;  // calculate bbox based on index
  BoundingBox3f mBbox;
};

static OctreeNode* buildNode(Octree* tree, const BoundingBox3f& bbox,
                             const std::vector<uint32_t>& indices);

class Octree {
 public:
  Octree(Mesh* mesh);
  void build();

  void rayIntersect(const Ray3f& ray, float& u, float& v, float& t);

  Mesh* getMeshPtr() { return mMesh; }

 private:
  Mesh* mMesh;
  OctreeNode* mRootNode;
};

NORI_NAMESPACE_END
