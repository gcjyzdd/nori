#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <pcg32.h>

namespace {
const float eps = 1e-3F;
const float PT = 0.95F;  // Russian Roulette Path Termination
const int MIN_PATH = 3;
}  // namespace

NORI_NAMESPACE_BEGIN

class EmitterSamplingIntegrator : public Integrator {
 public:
  EmitterSamplingIntegrator(
      const PropertyList &props) { /* No parameters this time */
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
    if (!scene->hasEmitter()) {
      return Color3f(0.0F);
    }

    Color3f color(1.F);
    Color3f emission(0.F);
    float pdf = 1.0F;
    float eta = 1.F;
    int paths = 0;
    Ray3f ray1(ray);
    while (true) {
      /* Find the surface that is visible in the requested direction */
      Intersection its;
      if (!scene->rayIntersect(ray1, its)) {
        color = Color3f(0.0F);
        break;
      }

      if (ray.rowIdx == 362 && ray.columnIdx == 284) {
        // cout << "hello\n";
      }

      auto bsdf = its.mesh->getBSDF();
      if (bsdf->isDiffuse()) {
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
            float gxy = its.shFrame.n.dot(wo) * shadowIts.shFrame.n.dot(-wo) /
                        (segLen * segLen);
            BSDFQueryRecord bsdfQuery(its.shFrame.toLocal(-ray.d),
                                      its.shFrame.toLocal(wo),
                                      EMeasure::ESolidAngle);
            auto albedo = bsdf->eval(bsdfQuery);
            Color3f Lr = albedo.array() * rec.color.array() * gxy;
            emission += color.array() * Lr.array() / (rec.pdf);
          }
        }
      }

      BSDFQueryRecord bsdfQuery(its.shFrame.toLocal(-ray1.d));
      auto c = bsdf->sample(bsdfQuery, sampler->next2D());

      eta *= bsdfQuery.eta;
      float pt =
          paths <= MIN_PATH
              ? 1.F
              : std::min(0.99F, (emission + color).maxCoeff() * eta * eta);

      if (sampler->next1D() < pt) {
        color = color.array() * c.array() / pt;
        ray1 = Ray3f(its.p, its.shFrame.toWorld(bsdfQuery.wo));
      } else {
        color = Color3f(0.F);
        break;
      }
      ++paths;
    }

    /* Return the component-wise absolute
       value of the shading normal as a color */

    return emission;
  }

  std::string toString() const { return "MaterialSamplingIntegrator[]"; }
};

NORI_REGISTER_CLASS(EmitterSamplingIntegrator, "path_ems");
NORI_NAMESPACE_END