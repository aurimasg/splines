
#pragma once


#include "Base.h"


/**
 * Represents a point in a two-dimensional plane defined as two floating point
 * values X and Y.
 *
 * This struct can be used as a vector and as a point.
 *
 * @see FloatPoint3
 * @see FloatPoint4
 */
struct FloatPoint final {

    /**
     * Constructs point with X and Y values set to zero.
     */
    constexpr FloatPoint() {
    }


    /**
     * Constructs point with given X and Y values.
     */
    constexpr FloatPoint(const double x, const double y);

public:

    /**
     * Returns squared distance from this point to another point.
     */
    double DistanceToSquared(const FloatPoint &point) const;


    /**
     * Returns length of vector defined by X and Y components of this point.
     */
    double Length() const;


    /**
     * Returns squared length of vector defined by X and Y components of this
     * point.
     */
    double LengthSquared() const;


    /**
     * Returns normalized version of this vector. If this vector has length of
     * zero, vector with both components set to zero will be returned.
     */
    FloatPoint UnitVector() const;


    /**
     * Returns dot product of this vector and a given vector.
     */
    double Dot(const FloatPoint &point) const;

public:

    double X = 0;
    double Y = 0;
};


static inline void operator+=(FloatPoint &v, const FloatPoint p) {
    v.X += p.X;
    v.Y += p.Y;
}


static inline FloatPoint operator+(const FloatPoint p0, const FloatPoint p1) {
    return FloatPoint(p0.X + p1.X, p0.Y + p1.Y);
}


static constexpr FloatPoint operator-(const FloatPoint p0,
    const FloatPoint p1)
{
    return FloatPoint(p0.X - p1.X, p0.Y - p1.Y);
}


static inline FloatPoint operator*(const FloatPoint p, const double s) {
    return FloatPoint(p.X * s, p.Y * s);
}


static inline FloatPoint operator*(const double s, const FloatPoint p) {
    return FloatPoint(p.X * s, p.Y * s);
}


static inline FloatPoint operator/(const FloatPoint p, const double s) {
    return FloatPoint(p.X / s, p.Y / s);
}


inline constexpr FloatPoint::FloatPoint(const double x, const double y)
:   X(x),
    Y(y)
{
}


inline double FloatPoint::DistanceToSquared(const FloatPoint &point) const {
    return (point - *this).LengthSquared();
}


inline double FloatPoint::Length() const {
    return Sqrt(LengthSquared());
}


inline double FloatPoint::LengthSquared() const {
    return X * X + Y * Y;
}


inline double FloatPoint::Dot(const FloatPoint &point) const {
    return (X * point.X) + (Y * point.Y);
}
