#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <pcg32.h>

namespace {
pcg32 rng;
const float eps = 1e-3F;
}  // namespace

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
 public:
  WhittedIntegrator(const PropertyList &props) { /* No parameters this time */
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
    /* Find the surface that is visible in the requested direction */
    Intersection its;
    if (!scene->rayIntersect(ray, its)) return Color3f(0.0F);

    if (!scene->hasEmitter()) return Color3f(0);

    Color3f color(0.F);

    float cosThetaX = its.shFrame.n.dot(-ray.d);
    EmitterQueryRecord rec;
    if (its.mesh->isEmitter())
      color += its.mesh->getEmitter()->eval(rec) * cosThetaX / (its.t * its.t);
    for (int i = 0; i < NUM_SAMPLE; ++i) {
      scene->sampleEmitter(rec, rng.nextFloat(),
                           Point2f(rng.nextFloat(), rng.nextFloat()));
      Point3f seg = rec.p - its.p;
      float segLen = seg.norm();
      Point3f wo = seg.normalized();
      Ray3f shadowRay(its.p, wo);
      Intersection shadowIts;
      if (!scene->rayIntersect(shadowRay, shadowIts) ||
          (segLen - shadowIts.t) > eps)
        continue;
      float cosThetaY = shadowIts.shFrame.n.dot(-shadowRay.d);
      if (cosThetaY < 0) continue; // backface

      BSDFQueryRecord bsdfQuery(its.shFrame.toLocal(-ray.d),
                                its.shFrame.toLocal(wo), EMeasure::ESolidAngle);
      auto albedo = its.mesh->getBSDF()->eval(bsdfQuery);
      color += albedo.array() * rec.color.array() * cosThetaX * cosThetaY /
               (rec.pdf * shadowIts.t * shadowIts.t);
      if (color.prod() < 0) {
        cout << "error";
      }
    }

    /* Return the component-wise absolute
       value of the shading normal as a color */

    return color / NUM_SAMPLE;
  }

  std::string toString() const { return "WhittedIntegrator[]"; }

  enum { NUM_SAMPLE = 1 };
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END