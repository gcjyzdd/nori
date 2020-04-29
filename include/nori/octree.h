#pragma once
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

static const unsigned int NUM_NODE = 8;
static const unsigned int N_LEAF = 8;

class Octree;
class OctreeNode {
 public:
  OctreeNode(Octree* tree) : mRoot{tree} {}

  static OctreeNode* build();

 private:
  Octree* mRoot;
  OctreeNode* mChildren[NUM_NODE];
  uint32_t mIndices[N_LEAF];

  // define bbox of the node
  Point3f mBoxMin;
  Point3f mBoxMax;
};

class Octree {
 public:
  Octree(Mesh* mesh);
  void build();

  void rayIntersect(const Ray3f& ray, float& u, float& v, float& t);

 private:
  Mesh* mMesh;
  OctreeNode* mRootNode;
};

NORI_NAMESPACE_END
