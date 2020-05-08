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

#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

namespace {
float rationalFun(const nori::Vector3f &wv, const nori::Vector3f &wh,
                  float alpha) {
  float c = wv.dot(wh) * wv(2);
  if (c <= 0) return 0;

  float b = 1.F / (alpha * std::sqrt(1.F - wv(2) * wv(2)) / wv(2));
  if (b < 1.6F) {
    float b2 = b * b;
    return (3.535F * b + 2.181F * b2) / (1.F + 2.276F * b + 2.577F * b2);
  } else {
    return 1.F;
  }
}
}  // namespace

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
 public:
  Microfacet(const PropertyList &propList) {
    /* RMS surface roughness */
    m_alpha = propList.getFloat("alpha", 0.1f);

    /* Interior IOR (default: BK7 borosilicate optical glass) */
    m_intIOR = propList.getFloat("intIOR", 1.5046f);

    /* Exterior IOR (default: air) */
    m_extIOR = propList.getFloat("extIOR", 1.000277f);

    /* Albedo of the diffuse base material (a.k.a "kd") */
    m_kd = propList.getColor("kd", Color3f(0.5f));

    /* To ensure energy conservation, we must scale the
       specular component by 1-kd.

       While that is not a particularly realistic model of what
       happens in reality, this will greatly simplify the
       implementation. Please see the course staff if you're
       interested in implementing a more realistic version
       of this BRDF. */
    m_ks = 1 - m_kd.maxCoeff();
  }

  /// Evaluate the BRDF for the given pair of directions
  Color3f eval(const BSDFQueryRecord &bRec) const {
    auto &wi = bRec.wi;
    auto &wo = bRec.wo;
    Vector3f wh = wi + wo;
    wh = wh / (wh.dot(wh));

    float D = Warp::squareToBeckmannPdf(wh, m_alpha);
    float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
    float num = m_ks * D * F * rationalFun(wi, wh, m_alpha) *
                rationalFun(wo, wh, m_alpha);
    float den = 4.F * wi(2) * wo(2) * wh(2);
    return m_kd / M_PI + num / den;
  }

  /// Evaluate the sampling density of \ref sample() wrt. solid angles
  float pdf(const BSDFQueryRecord &bRec) const {
    auto &wi = bRec.wi;
    auto &wo = bRec.wo;
    Vector3f wh = wi + wo;
    float Jh = 0.25F / (wh.dot(wo));
    return m_ks * Warp::squareToBeckmannPdf(wh, m_alpha) * Jh +
           (1.F - m_ks) * wo(2) / M_PI;
  }

  /// Sample the BRDF
  Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
    float epsilon = _sample(0);
    m_dpdf.sampleReuse(epsilon);
    if (_sample(0) < m_ks) {  // specular
      bRec.wo = Warp::squareToUniformHemisphere(_sample);
    } else {
      bRec.wo = Warp::squareToCosineHemisphere(_sample);
    }

    // Note: Once you have implemented the part that computes the scattered
    // direction, the last part of this function should simply return the
    // BRDF value divided by the solid angle density and multiplied by the
    // cosine factor from the reflection equation, i.e.
    return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
  }

  bool isDiffuse() const {
    /* While microfacet BRDFs are not perfectly diffuse, they can be
       handled by sampling techniques for diffuse/non-specular materials,
       hence we return true here */
    return true;
  }

  std::string toString() const {
    return tfm::format(
        "Microfacet[\n"
        "  alpha = %f,\n"
        "  intIOR = %f,\n"
        "  extIOR = %f,\n"
        "  kd = %s,\n"
        "  ks = %f\n"
        "]",
        m_alpha, m_intIOR, m_extIOR, m_kd.toString(), m_ks);
  }

 private:
  float m_alpha;
  float m_intIOR, m_extIOR;
  float m_ks;
  DiscretePDF m_dpdf;
  Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
