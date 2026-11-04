
#pragma once


#include "FloatPoint.h"


/**
 * Cubic curve with implied start and end points at (0, 0) and (1, 0)
 * respectively.
 */
struct UnitCubicCurve final {
    UnitCubicCurve(const FloatPoint p0, const FloatPoint p1);

    const FloatPoint P0 = FloatPoint(1.0 / 3.0, 0.0);
    const FloatPoint P1 = FloatPoint(2.0 / 3.0, 0.0);

private:

    UnitCubicCurve() = delete;
};


inline UnitCubicCurve::UnitCubicCurve(const FloatPoint p0, const FloatPoint p1)
:   P0(p0),
    P1(p1)
{
}
