#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

static const float OneOver4Pi2 = 0.25F / (M_PI * M_PI);

class SimpleIntegrator : public Integrator {
 public:
  explicit SimpleIntegrator(const PropertyList &props)
      : mPos{props.getPoint("position")}, mColor{props.getColor("energy")} {}

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
    /* Find the surface that is visible in the requested direction */
    Intersection its;
    if (!scene->rayIntersect(ray, its)) return Color3f(0.0F);

    Point3f hit = its.p;
    Point3f d = mPos - hit;
    Point3f nd = d.normalized();
    Ray3f shadowRay(hit, nd);
    Intersection shadowIts;

    if (scene->getAccel()->rayIntersect(shadowRay, shadowIts, true))
      // if (scene->rayIntersect(shadowRay, shadowIts))
      return Color3f(0.0F);  // light is occluded

    /* Return the component-wise absolute
       value of the shading normal as a color */
    float cosTheta = std::max(its.shFrame.toLocal(nd)(2), 0.F);
    return mColor * (cosTheta * OneOver4Pi2) / d.dot(d);
  }

  std::string toString() const { return "SimpleIntegrator[]"; }

  Point3f mPos;
  Color3f mColor;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END