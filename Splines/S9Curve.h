
#pragma once


#include "Base.h"
#include "CurveArray.h"


/*
 * Rounding rate K of the curvature law of S9Curve. The corner rate at
 * which the rounding term has half of its full weight.
 */
static constexpr double S9CornerRoundingRate = 0.5;


STATIC_ASSERT(S9CornerRoundingRate > 0.0);


static constexpr double S9CornerRounding(const double rate) noexcept {
    return rate / (rate + S9CornerRoundingRate);
}


/*
 * Default L of the rounding term.
 */
static constexpr double S9DefaultCornerRoundingLift = 1.5;


STATIC_ASSERT(S9DefaultCornerRoundingLift >= 0.0);


/**
 * Keeps the inputs of one S9 curve.
 */
struct S9CurveParameters final {

    /**
     * Angle alpha of the start in radians, counter-clockwise from the
     * chord.
     */
    double Alpha = 0.0;


    /**
     * Angle beta of the end in radians, clockwise from the chord.
     */
    double Beta = 0.0;


    /**
     * Signed curvature length of the start, from minus one to one. Magnitude
     * sets the corner rate of the start, one for an end without a corner.
     * Negative value reverses the start.
     */
    double StartLength = 1.0;


    /**
     * Signed curvature length of the end, in the same range and with the
     * same meaning.
     */
    double EndLength = 1.0;


    /**
     * L of the rounding term of the start corner, at least zero.
     */
    double StartRoundingLift = S9DefaultCornerRoundingLift;


    /**
     * L of the rounding term of the end corner, at least zero.
     */
    double EndRoundingLift = S9DefaultCornerRoundingLift;
};


/**
 * Keeps one S9 curve for an implicit unit chord and evaluates it.
 *
 * The chord runs from (0, 0) to (1, 0) with an upward positive Y axis,
 * and a positive curvature turns counter-clockwise. The turn rate per unit
 * of normalized arc length t, from zero through one, is
 *
 *     f(t) = StartLine * g0(t) + EndLine * g1(t)
 *     g0(t) = (1 - t) + c0(t)
 *     g1(t) = t + c1(1 - t)
 *     x = d * R
 *     c(d) = Q (p(x) + phi q(x))
 *     p(x) = 1 / (1 + x)^1.5
 *     q(x) = L x / (1 + x)^3
 *     phi = R / (R + K)
 */
struct S9Curve final {

    /**
     * Forward tangent angle in radians at normalized arc length zero,
     * measured counter-clockwise from the chord. It is the angle alpha of
     * the start, turned by PI for a reversed start.
     */
    double StartAngle = 0.0;


    /**
     * Forward tangent angle in radians at normalized arc length one,
     * measured counter-clockwise from the chord. It is minus the angle beta
     * of the end, turned by PI for a reversed end.
     */
    double EndAngle = 0.0;


    /**
     * Rate R of the start corner.
     */
    double StartR = 0.0;


    /**
     * Rate R of the end corner.
     */
    double EndR = 0.0;


    /**
     * L of the rounding term of the start corner, the factor of q in
     * the curvature law. The solve copies it from
     * S9CurveParameters::StartRoundingLift.
     */
    double StartL = S9DefaultCornerRoundingLift;


    /**
     * L of the rounding term of the end corner. The solve copies it
     * from S9CurveParameters::EndRoundingLift.
     */
    double EndL = S9DefaultCornerRoundingLift;


    /**
     * Net turn of the curve in radians. It is the difference of the two
     * angles, possibly plus or minus a full turn.
     */
    double Turn = 0.0;


    /**
     * Turn rate of the line part at the start. The solve writes it.
     */
    double StartLine = 0.0;


    /**
     * Turn rate of the line part at the end. The solve writes it.
     */
    double EndLine = 0.0;


    /**
     * Length of the curve for a chord of length one.
     */
    double Length = 1.0;


    /**
     * Solves the curve for given curvature lengths.
     *
     * @param parameters Angles, curvature lengths and rounding lifts of the
     * curve. Both angles must be finite, both lengths from minus one
     * to one, and both lifts at least zero. It is only read during the
     * call.
     *
     * Returns true when a curve that meets both tangent angles was found.
     */
    bool TrySolveWithLengths(const S9CurveParameters &parameters) noexcept;


    /**
     * Solves the curve for given lengths from a nearby solved curve.
     *
     * @param parameters Angles, curvature lengths and rounding lifts of the
     * curve, as for TrySolveWithLengths. It is only read during the call.
     *
     * @param nearby Solved curve of nearby angles. It is only read. It can be
     * this object.
     *
     * Returns true when a curve that meets both tangent angles was found.
     */
    bool TrySolveWithLengthsNear(const S9CurveParameters &parameters,
        const S9Curve &nearby) noexcept;


    /**
     * Approximates the solved curve with cubic curves from (0, 0) to (1, 0)
     * with an upward positive Y axis.
     *
     * Failure clears curves and leaves one point at (0, 0). The causes
     * are a curve that was not solved, an exhausted measurement limit, and an
     * exhausted cubic curve limit.
     *
     * @param maximumError Positive limit for the distance between the exact
     * curve and the cubic curves, in unit chord coordinates.
     *
     * @param curves Destination sequence. Its contents are replaced and its
     * storage is reused. No reference to it is kept.
     *
     * Returns true when a curve was written.
     */
    bool TryApproximateWithCubics(const double maximumError,
        CurveArray<CubicCurve> &curves) const noexcept;


    /**
     * Returns the turn J0 = 1/2 + R0 that the start basis makes over the
     * whole curve.
     */
    constexpr double StartBasisTurn() const noexcept;


    /**
     * Returns the turn J1 = 1/2 + R1 that the end basis makes over the whole
     * curve.
     */
    constexpr double EndBasisTurn() const noexcept;


    /**
     * Returns the start basis curvature g0(t) = (1 - t) + c0(t) at one
     * normalized arc length.
     *
     * @param t Normalized arc length from zero through one.
     */
    double StartBasis(const double t) const noexcept;


    /**
     * Returns the end basis curvature g1(t) = t + c1(1 - t) at one
     * normalized arc length.
     *
     * @param t Normalized arc length from zero through one.
     */
    double EndBasis(const double t) const noexcept;


    /**
     * Returns the integral of the start basis curvature from zero to t.
     *
     * @param t Normalized arc length from zero through one.
     */
    double StartBasisIntegral(const double t) const noexcept;


    /**
     * Returns the integral of the end basis curvature from zero to t.
     *
     * @param t Normalized arc length from zero through one.
     */
    double EndBasisIntegral(const double t) const noexcept;


    /**
     * Returns the turn rate f(t) of the curve at one normalized arc length.
     *
     * @param t Normalized arc length from zero through one.
     */
    double TurnRate(const double t) const noexcept;


    /**
     * Returns the signed curvature of the curve for a chord of length one.
     *
     * @param t Normalized arc length from zero through one.
     */
    double Curvature(const double t) const noexcept;


    /**
     * Returns the change of the tangent angle in radians from the start to t.
     *
     * @param t Normalized arc length from zero through one.
     */
    double TurnFromStart(const double t) const noexcept;


    /**
     * Returns the tangent angle at one normalized arc length, measured from
     * the chord. Needs a solved curve.
     *
     * @param t Normalized arc length from zero through one.
     */
    double TangentAngle(const double t) const noexcept;
};


inline constexpr double S9Curve::StartBasisTurn() const noexcept {
    return 0.5 + StartR;
}


inline constexpr double S9Curve::EndBasisTurn() const noexcept {
    return 0.5 + EndR;
}


inline double S9Curve::Curvature(const double t) const noexcept
{
    ASSERT(Length > 0.0);

    const double turnRate = TurnRate(t);

    return turnRate / Length;
}


inline double S9Curve::TangentAngle(const double t) const noexcept
{
    const double turn = TurnFromStart(t);

    return StartAngle + turn;
}
