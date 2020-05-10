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

class MaterialSamplingIntegrator : public Integrator {
 public:
  MaterialSamplingIntegrator(
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

      EmitterQueryRecord rec;
      if (its.mesh->isEmitter()) {
        emission += its.mesh->getEmitter()->eval(rec).array() * color.array();
      }

      if (ray.rowIdx == 362 && ray.columnIdx == 284) {
        // cout << "hello\n";
      }

      auto bsdf = its.mesh->getBSDF();

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

NORI_REGISTER_CLASS(MaterialSamplingIntegrator, "path_mats");
NORI_NAMESPACE_END