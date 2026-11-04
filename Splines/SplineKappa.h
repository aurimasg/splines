
#pragma once


#include "Base.h"
#include "CurveArray.h"
#include "FloatPoint.h"


/**
 * Solves a closed Kappa spline as connected quadratic curves.
 *
 * @param knots Array of count finite knot positions. The function reads it
 * during the call, takes no ownership, and keeps no reference. A null
 * pointer produces an empty result.
 *
 * @param count Number of knots. At least three are required. A smaller count
 * produces an empty result.
 *
 * @param curves Destination whose contents are replaced. Storage is reused.
 * The result is empty if the iteration produces nonfinite coordinates.
 */
extern void SolveKappaCyclic(const FloatPoint *knots, const U32 count,
    CurveArray<QuadraticCurve> &curves);


/**
 * Solves an open Kappa spline as connected quadratic curves.
 *
 * @param knots Array of count finite knot positions. The function reads it
 * during the call, takes no ownership, and keeps no reference. A null
 * pointer produces an empty result.
 *
 * @param count Number of knots. At least four are required. A smaller count
 * produces an empty result.
 *
 * @param curves Destination whose contents are replaced. Storage is reused.
 * The result is empty if the iteration produces nonfinite coordinates.
 */
extern void SolveKappaAcyclic(const FloatPoint *knots, const U32 count,
    CurveArray<QuadraticCurve> &curves);
