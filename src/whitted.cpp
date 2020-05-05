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

    EmitterQueryRecord rec;
    if (its.mesh->isEmitter()) {
      color += its.mesh->getEmitter()->eval(rec);
    }

    for (int i = 0; i < NUM_SAMPLE; ++i) {
      scene->sampleEmitter(rec, rng.nextFloat(),
                           Point2f(rng.nextFloat(), rng.nextFloat()));
      Point3f seg = rec.p - its.p;
      float segLen = seg.norm();
      Point3f wo = seg.normalized();
      if (its.shFrame.n.dot(wo) <= eps)
        continue;  // light source below the surface
      Ray3f shadowRay(its.p, wo);
      Intersection shadowIts;
      if (!scene->rayIntersect(shadowRay, shadowIts) ||
          (segLen - shadowIts.t) > eps)
        continue;

      float gxy = its.shFrame.n.dot(wo) * shadowIts.shFrame.n.dot(-wo) /
                  (segLen * segLen);
      if (gxy < 0) {
        cout << "error: r = " << ray.rowIdx << ", c = " << ray.columnIdx
             << ", cos th = " << its.shFrame.n.dot(wo) << "\n";
      }
      BSDFQueryRecord bsdfQuery(its.shFrame.toLocal(-ray.d),
                                its.shFrame.toLocal(wo), EMeasure::ESolidAngle);
      auto albedo = its.mesh->getBSDF()->eval(bsdfQuery);
      Color3f Lr = albedo.array() * rec.color.array() * gxy;
      color += Lr / (rec.pdf);
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