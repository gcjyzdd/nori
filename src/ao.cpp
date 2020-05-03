#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN

pcg32 rng;

class AoIntegrator : public Integrator {
 public:
  AoIntegrator(const PropertyList &props) { /* No parameters this time */
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
    /* Find the surface that is visible in the requested direction */
    Intersection its;
    if (!scene->rayIntersect(ray, its)) return Color3f(0.0f);

    Point2f sample;
    Color3f color(0.F);
    for (int i = 0; i < NUM_SAMPLE; ++i) {
      sample = Point2f(rng.nextFloat(), rng.nextFloat());
      auto d = Warp::squareToCosineHemisphere(sample);
      auto f = Warp::squareToCosineHemispherePdf(d);

      auto localDir = its.shFrame.toWorld(d);
      Ray3f shadowRay(its.p, localDir);
      Intersection shadowIts;
      if (scene->getAccel()->rayIntersect(shadowRay, shadowIts, true)) continue;

      color += std::max(localDir.dot(its.shFrame.n), 0.F) / (M_PI * f);
    }

    return color / NUM_SAMPLE;
  }

  std::string toString() const { return "AoIntegrator[]"; }

 private:
  enum { NUM_SAMPLE = 2 };
};

NORI_REGISTER_CLASS(AoIntegrator, "ao");
NORI_NAMESPACE_END