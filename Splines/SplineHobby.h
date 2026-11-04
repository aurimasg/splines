
#pragma once


#include "Base.h"
#include "CurveArray.h"
#include "FloatPoint.h"


/**
 * Implements Hobby spline.
 *
 * @param count A number of knots in knot array. If this value is less than 1,
 * this function does nothing.
 *
 * @param knots Pointer to knots. If `count` is greater than 0, this pointer
 * must not be nullptr.
 */
extern void SolveHobbyCyclic(const U32 count, const FloatPoint *knots,
    CurveArray<CubicCurve> &curves);


extern void SolveHobbyAcyclic(const U32 count, const FloatPoint *knots,
    CurveArray<CubicCurve> &curves, const double omega);
