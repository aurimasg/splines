
#pragma once


#include "Base.h"
#include "CurveArray.h"
#include "FloatPoint.h"


/**
 * Solves an open spline through the specified knots.
 *
 * Input requirements
 * The input must contain at least two knots with finite coordinates. The
 * distance between consecutive knots must be greater than `DBL_EPSILON`. The
 * function does not take ownership of the knot array and reads it only for
 * the duration of the call.
 *
 * @param count Number of elements in knots. The value must be at least two.
 *
 * @param knots Pointer to count input knots. The pointer must not be null.
 *
 * @param curves Receives one cubic curve for each pair of consecutive knots.
 * Its previous contents are replaced.
 */
extern void SolveCubicG2Acyclic(const U32 count, const FloatPoint *const knots,
    CurveArray<CubicCurve> &curves) noexcept;


/**
 * Solves a closed spline through the specified knots.
 *
 * Input requirements
 * The input must contain at least three knots with finite coordinates. The
 * distance between consecutive knots and between the final and first knots
 * must be greater than `DBL_EPSILON`. The function does not take ownership of
 * the knot array and reads it only for the duration of the call.
 *
 * @param count Number of elements in knots. The value must be at least three.
 *
 * @param knots Pointer to count input knots. The pointer must not be null.
 *
 * @param curves Receives count cubic curves. The final curve ends at the
 * first knot. Its previous contents are replaced.
 */
extern void SolveCubicG2Cyclic(const U32 count, const FloatPoint *const knots,
    CurveArray<CubicCurve> &curves) noexcept;
