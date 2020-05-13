#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/mesh.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <pcg32.h>

namespace {
const float eps = 1e-3F;
const float PT = 0.95F;  // Russian Roulette Path Termination
const int MIN_PATH = 3;
}  // namespace

NORI_NAMESPACE_BEGIN

class MisIntegrator : public Integrator {
 public:
  MisIntegrator(const PropertyList &props) { /* No parameters this time */
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
    if (!scene->hasEmitter()) {
      return Color3f(0.0F);
    }

    Color3f throughput(1.F);
    Color3f emission(0.F);
    float pdf = 1.0F;
    float eta = 1.F;
    int paths = 0;
    bool lastSpecular = true;
    float pLight = 0;
    float pBrdf = 1;
    bool first = true;
    Ray3f ray1(ray);
    Intersection prevIts;
    Color3f prevAlbedo;
    while (true) {
      /* Find the surface that is visible in the requested direction */
      Intersection its;
      if (!scene->rayIntersect(ray1, its)) {
        throughput = Color3f(0.0F);
        break;
      }

      if (its.mesh->isEmitter()) {
        EmitterQueryRecord rec;
        if (lastSpecular || first) {
          emission +=
              its.mesh->getEmitter()->eval(rec).array() * throughput.array();
        } else {
          rec.mesh = its.mesh;
          rec.triIndex = its.face;
          rec.pdf = scene->pdfEmitter(rec);
          auto lc = its.mesh->getEmitter()->eval(rec);

          float ny = its.shFrame.n.dot(-ray1.d);
          if (ny > 0) {
            float r2 = its.t * its.t;

            float sa = ny / (rec.pdf * r2);
            pLight = 1.F / sa;
            emission +=
                throughput.array() * lc.array() * pBrdf / (pBrdf + pLight);
          }
        }
      }
      first = false;

      auto bsdf = its.mesh->getBSDF();
      if (bsdf->isDiffuse()) {
        // sample a light source
        EmitterQueryRecord rec;
        scene->sampleEmitter(rec, sampler->next1D(), sampler->next2D());
        Point3f seg = rec.p - its.p;
        float segLen = seg.norm();
        Point3f wo = seg.normalized();
        if (its.shFrame.n.dot(wo) >= eps) {
          Ray3f shadowRay(its.p, wo);
          Intersection shadowIts;
          if (!scene->rayIntersect(shadowRay, shadowIts) ||
              (segLen - shadowIts.t) <= eps) {
            float nx = its.shFrame.n.dot(wo);
            float ny = shadowIts.shFrame.n.dot(-wo);
            float r2 = segLen * segLen;
            float gxy = nx * ny / r2;
            if (gxy > 0) {
              BSDFQueryRecord bsdfQuery(its.shFrame.toLocal(-ray1.d),
                                        its.shFrame.toLocal(wo),
                                        EMeasure::ESolidAngle);
              auto albedo = bsdf->eval(bsdfQuery);
              float sa = ny / (rec.pdf * r2);
              pLight = 1.F / sa;
              pBrdf = bsdf->pdf(bsdfQuery);
              Color3f Lr = albedo.array() * rec.color.array() * gxy;
              float wl = pLight / (pBrdf + pLight);
              emission += throughput.array() * Lr.array() / rec.pdf * wl;
            }
          }
        }
        lastSpecular = false;
      } else {
        lastSpecular = true;
      }

      BSDFQueryRecord bsdfQuery(its.shFrame.toLocal(-ray1.d));
      auto c = bsdf->sample(bsdfQuery, sampler->next2D());
      pBrdf = bsdf->pdf(bsdfQuery);
      prevAlbedo = bsdf->eval(bsdfQuery);

      eta *= bsdfQuery.eta;
      float pt = paths <= MIN_PATH
                     ? 1.F
                     : std::min(0.99F, throughput.maxCoeff() * eta * eta);

      if (sampler->next1D() < pt) {
        throughput = throughput.array() * c.array() / pt;
        ray1 = Ray3f(its.p, its.shFrame.toWorld(bsdfQuery.wo));
      } else {
        throughput = Color3f(0.F);
        break;
      }
      ++paths;
      prevIts = its;
    }

    /* Return the component-wise absolute
       value of the shading normal as a color */

    return emission;
  }

  std::string toString() const {
    return "MultipleImportanceSamplingIntegrator[]";
  }
};

NORI_REGISTER_CLASS(MisIntegrator, "path_mis");
NORI_NAMESPACE_END
