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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

inline float toTent(float x) { return x < 0.5F ? -1.F+std::sqrt(2.F*x) : 1-std::sqrt(2.F*(1.F-x)); }

inline float toTentPdf(float x) { return 1.F - std::abs(x); }

Point2f Warp::squareToTent(const Point2f &sample) {
    return Point2f(toTent(sample(0)), toTent(sample(1)));
}

float Warp::squareToTentPdf(const Point2f &p) {
    return toTentPdf(p(0)) * toTentPdf(p(1));
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r = std::sqrt(sample(0));
    return Point2f(r * std::cos(2 * M_PI * sample(1)),
                 r * std::sin(2 * M_PI * sample(1)));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {  
    return p.norm() >=1.F?0:1.F / M_PI;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float th = std::acos(2 * sample(0)-1);
    float phi = 2.F * M_PI * sample(1);
    float s = std::sin(th);
    return Vector3f(s * std::cos(phi), s * std::sin(phi), std::cos(th));
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return 0.25F / M_PI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float th = std::acos(sample(0));
    float phi = 2.F * M_PI * sample(1);
    float s = std::sin(th);
    return Vector3f(s * std::cos(phi), s * std::sin(phi), std::cos(th));
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return v(2) >=0 ? 0.5F / M_PI : 0;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    float th = std::asin(std::sqrt(sample(0)));
    float phi = 2.F * M_PI * sample(1);
    float s = std::sin(th);
    return Vector3f(s * std::cos(phi), s * std::sin(phi), std::cos(th));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    return v(2) >= 0 ? v(2) * 1.F / M_PI : 0;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
}

NORI_NAMESPACE_END
