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

#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

class Mesh;

struct EmitterQueryRecord {
  /// Outgoing direction (in the local frame)
  Point3f p;

  Point3f normal;

  float pdf;

  Color3f color;

  const Mesh* mesh{nullptr};

  uint32_t triIndex{0};
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
 public:
  virtual void sample(EmitterQueryRecord& record,
                      const Point2f& sample) const = 0;

  virtual Color3f eval(const EmitterQueryRecord& record) const = 0;

  /**
   * \brief Return the type of object (i.e. Mesh/Emitter/etc.)
   * provided by this instance
   * */
  EClassType getClassType() const { return EEmitter; }
};

NORI_NAMESPACE_END
