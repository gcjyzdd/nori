#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <pcg32.h>

namespace {
const float eps = 1e-3F;
const float pt = 0.95F;  // Russian Roulette Path Termination
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

    if (ray.rowIdx == 387 && ray.columnIdx == 296) {
      cout << "hello\n";
    }

    auto bsdf = its.mesh->getBSDF();
    if (bsdf->isDiffuse()) {
      scene->sampleEmitter(rec, sampler->next1D(), sampler->next2D());
      Point3f seg = rec.p - its.p;
      float segLen = seg.norm();
      Point3f wo = seg.normalized();
      if (its.shFrame.n.dot(wo) <= eps)
        return color;  // light source below the surface
      Ray3f shadowRay(its.p, wo);
      Intersection shadowIts;
      if (!scene->rayIntersect(shadowRay, shadowIts) ||
          (segLen - shadowIts.t) > eps)
        return color;

      float gxy = its.shFrame.n.dot(wo) * shadowIts.shFrame.n.dot(-wo) /
                  (segLen * segLen);
      if (gxy < 0) {
        return color; // backface
        cout << "error: r = " << ray.rowIdx << ", c = " << ray.columnIdx
             << ", cos thx = " << its.shFrame.n.dot(wo)
             << ", cos thy = " << shadowIts.shFrame.n.dot(-wo) << "\n";
      }
      BSDFQueryRecord bsdfQuery(its.shFrame.toLocal(-ray.d),
                                its.shFrame.toLocal(wo), EMeasure::ESolidAngle);
      auto albedo = bsdf->eval(bsdfQuery);
      Color3f Lr = albedo.array() * rec.color.array() * gxy;
      color += Lr / (rec.pdf);
    } else {  // specular
      BSDFQueryRecord bsdfQuery(its.shFrame.toLocal(-ray.d));
      auto c = bsdf->sample(bsdfQuery, sampler->next2D());

      if (sampler->next1D() < pt) {
        Ray3f newRay(its.p, its.shFrame.toWorld(bsdfQuery.wo));
        return c.array() * Li(scene, sampler, newRay) / pt;
      } else {
        return color;
      }
    }

    /* Return the component-wise absolute
       value of the shading normal as a color */

    return color;
  }

  std::string toString() const { return "WhittedIntegrator[]"; }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END