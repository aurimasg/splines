
#pragma once


#include "CurveArray.h"
#include "FloatPoint.h"


/**
 * Control points of the two ends of an open S9 spline.
 */
struct S9SplineEndControls final {
    FloatPoint StartControl;
    FloatPoint EndControl;
    bool HasStartControl = false;
    bool HasEndControl = false;
};


extern void S9SolveAcyclicExplicit(const U32 count,
    const FloatPoint *const knots, const S9SplineEndControls &controls,
    CurveArray<CubicCurve> &curves) noexcept;


extern void S9SolveAcyclic(const U32 count, const FloatPoint *const knots,
    CurveArray<CubicCurve> &curves) noexcept;


extern void S9SolveCyclic(const U32 count, const FloatPoint *const knots,
    CurveArray<CubicCurve> &curves) noexcept;
