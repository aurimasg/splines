
#pragma once


#include "Base.h"
#include "CurveArray.h"
#include "FloatPoint.h"


/**
 * Solves a closed cubic spline through the input knots.
 *
 * Continuity
 *
 * Equal parameter intervals give C2 continuity at all joins, including
 * the join between the last and first curves.
 *
 * @param count Number of input knots. Fewer than three leave curves unchanged.
 *
 * @param knots Borrowed array of count finite knot positions. It must be
 * non-null when count is positive and valid throughout this call.
 * No ownership is transferred and no pointer is kept.
 *
 * @param curves Receives the solved spline when count is at least three.
 * Previous contents are replaced and storage is reused. No reference is
 * kept.
 */
extern void SolveCubicC2Cyclic(const U32 count,
    const FloatPoint *knots, CurveArray<CubicCurve> &curves);


/**
 * Solves an open cubic spline through the input knots.
 *
 * Continuity
 *
 * Equal parameter intervals give C2 continuity at all interior joins.
 * The second derivative is zero at both endpoints.
 *
 * @param count Number of input knots. Fewer than three leave curves unchanged.
 *
 * @param knots Borrowed array of count finite knot positions. It must be
 * non-null when count is positive and valid throughout this call.
 * No ownership is transferred and no pointer is kept.
 *
 * @param curves Receives the solved spline when count is at least three.
 * Previous contents are replaced and storage is reused. No reference is
 * kept.
 */
extern void SolveCubicC2Acyclic(const U32 count,
    const FloatPoint *knots, CurveArray<CubicCurve> &curves);
