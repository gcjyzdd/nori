#pragma once

#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
 public:
  explicit AreaLight(const PropertyList& props)
      : mColor{props.getColor("radiance")} {}

  void sample(EmitterQueryRecord& record, const Point2f& sample) const override;

  Color3f eval(const EmitterQueryRecord& record) const override {
    return mColor;
  }

  std::string toString() const { return "AreaLight[]"; }

 private:
  Color3f mColor;
};

void AreaLight::sample(EmitterQueryRecord& record,
                       const Point2f& sample) const {}

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
