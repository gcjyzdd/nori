#include "nori/octree.h"

NORI_NAMESPACE_BEGIN

OctreeNode* OctreeNode::build() { return nullptr; }

Octree::Octree(Mesh* mesh) : mMesh{mesh} {}

void Octree::build() {}

void Octree::rayIntersect(const Ray3f& ray, float& u, float& v, float& t) {}

NORI_NAMESPACE_END
