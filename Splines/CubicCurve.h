
#pragma once


#include "FloatPoint.h"


/**
 * Represents a two dimensional cubic curve.
 *
 * CubicCurve is defined by start point, end point and two control points.
 */
struct CubicCurve final {

    /**
     * Constructs cubic curve with given points.
     *
     * Constructs curve with coordinates specified as four points.
     *
     * @param p0 Start point of curve.
     *
     * @param p1 First control point of curve.
     *
     * @param p2 Second control point of curve.
     *
     * @param p3 End point of curve.
     */
    constexpr CubicCurve(const FloatPoint p0, const FloatPoint p1,
        const FloatPoint p2, const FloatPoint p3);

public:

    /**
     * Start point of the curve.
     */
    FloatPoint P0;


    /**
     * First control point of the curve.
     */
    FloatPoint P1;


    /**
     * Second control point of the curve.
     */
    FloatPoint P2;


    /**
     * End point.
     */
    FloatPoint P3;
};


inline constexpr CubicCurve::CubicCurve(const FloatPoint p0,
    const FloatPoint p1, const FloatPoint p2, const FloatPoint p3)
:   P0(p0),
    P1(p1),
    P2(p2),
    P3(p3)
{
}
