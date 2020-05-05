/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/warp.h>
#include <pcg32.h>
#include <Eigen/Geometry>

namespace {
pcg32 rng;
const static int MAGIC_NUMBER = 50;
}  // namespace

NORI_NAMESPACE_BEGIN

Mesh::Mesh() {}

Mesh::~Mesh() {
  delete m_bsdf;
  delete m_emitter;
}

void Mesh::activate() {
  if (!m_bsdf) {
    /* If no material was assigned, instantiate a diffuse BRDF */
    m_bsdf = static_cast<BSDF *>(
        NoriObjectFactory::createInstance("diffuse", PropertyList()));
  }
  uint32_t n = getTriangleCount();
  m_dpdf.reserve(n);
  m_area.reserve(n);
  for (uint32_t i = 0; i < n; ++i) {
    float a = surfaceArea(i);
    m_area.push_back(a);
    m_dpdf.append(a);
    m_reciprocal_area += a;
  }
  m_dpdf.normalize();
  m_reciprocal_area = 1.F / m_reciprocal_area;
  for (int i = 0; i < MAGIC_NUMBER; ++i) rng.nextFloat();
}

float Mesh::surfaceArea(uint32_t index) const {
  uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

  const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

  return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v,
                        float &t) const {
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

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
  BoundingBox3f result(m_V.col(m_F(0, index)));
  result.expandBy(m_V.col(m_F(1, index)));
  result.expandBy(m_V.col(m_F(2, index)));
  return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
  return (1.0f / 3.0f) * (m_V.col(m_F(0, index)) + m_V.col(m_F(1, index)) +
                          m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
  auto t = obj->getClassType();
  switch (t) {
    case EBSDF:
      if (m_bsdf)
        throw NoriException("Mesh: tried to register multiple BSDF instances!");
      m_bsdf = static_cast<BSDF *>(obj);
      break;

    case EEmitter: {
      Emitter *emitter = static_cast<Emitter *>(obj);
      if (m_emitter)
        throw NoriException(
            "Mesh: tried to register multiple Emitter instances!");
      m_emitter = emitter;
    } break;

    default:
      throw NoriException("Mesh::addChild(<%s>) is not supported!",
                          classTypeName(obj->getClassType()));
  }
}

std::string Mesh::toString() const {
  // cout << "center: " << m_V.array().rowwise().mean().transpose() << "\n";
  return tfm::format(
      "Mesh[\n"
      "  name = \"%s\",\n"
      "  vertexCount = %i,\n"
      "  triangleCount = %i,\n"
      "  bsdf = %s,\n"
      "  emitter = %s\n"
      "]",
      m_name, m_V.cols(), m_F.cols(),
      m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
      m_emitter ? indent(m_emitter->toString()) : std::string("null"));
}

void Mesh::unitSquare2tri(const Point2f &sample, Point3f &p, Point3f &n,
                          float &f) const {
  auto idx = m_dpdf.sample(rng.nextFloat());
  f = m_reciprocal_area;
  auto c = m_F.col(idx);

  float t = std::sqrt(sample(0));
  float alpha = 1.F - t;
  float beta = sample(1) * t;

  Point3f p0 = m_V.col(c(0));
  Point3f p1 = m_V.col(c(1));
  Point3f p2 = m_V.col(c(2));

  p = alpha * p0 + beta * p1 + (1.F - alpha - beta) * p2;
  if (m_N.size() > 0) {
    n = alpha * m_N.col(c(0)) + beta * m_N.col(c(1)) +
        (1.F - alpha - beta) * m_N.col(c(2));
  } else {
    n = (p1 - p0).cross(p2 - p0);
  }
  n.normalize();
}

void Mesh::sample(EmitterQueryRecord &record, const Point2f &sample) {
  unitSquare2tri(sample, record.p, record.normal, record.pdf);
  record.color = m_emitter->eval(record);
}

std::string Intersection::toString() const {
  if (!mesh) return "Intersection[invalid]";

  return tfm::format(
      "Intersection[\n"
      "  p = %s,\n"
      "  t = %f,\n"
      "  uv = %s,\n"
      "  shFrame = %s,\n"
      "  geoFrame = %s,\n"
      "  mesh = %s\n"
      "]",
      p.toString(), t, uv.toString(), indent(shFrame.toString()),
      indent(geoFrame.toString()),
      mesh ? mesh->toString() : std::string("null"));
}

NORI_NAMESPACE_END
