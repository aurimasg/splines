
#pragma once


#include "Base.h"
#include "CurveArray.h"


/*
 * Rounding λ_c of the curvature law of S9Curve. The corner λ at
 * which the rounding term has half of its full weight.
 */
static constexpr double S9RoundingLambda = 0.5;


STATIC_ASSERT(S9RoundingLambda > 0.0);


static constexpr double S9CornerRounding(const double lambda) noexcept {
    return lambda / (lambda + S9RoundingLambda);
}


/*
 * Default η of the rounding term.
 */
static constexpr double S9RoundingEta = 1.5;


STATIC_ASSERT(S9RoundingEta >= 0.0);


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
     * Signed weight of the start, from minus one to one. Magnitude sets the
     * corner λ of the start, one for an end without a corner. Negative
     * value reverses the start.
     */
    double StartWeight = 1.0;


    /**
     * Signed weight of the end, in the same range and with the same
     * meaning.
     */
    double EndWeight = 1.0;


    /**
     * η of the rounding term of the start corner, at least zero.
     */
    double StartEta = S9RoundingEta;


    /**
     * η of the rounding term of the end corner, at least zero.
     */
    double EndEta = S9RoundingEta;
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
 *     x = λ r
 *     c(r) = Q (p(x) + φ q(x))
 *     p(x) = 1 / (1 + x)^1.5
 *     q(x) = η x / (1 + x)^3
 *     φ = λ / (λ + λ_c)
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
     * λ of the start corner.
     */
    double StartLambda = 0.0;


    /**
     * λ of the end corner.
     */
    double EndLambda = 0.0;


    /**
     * η of the rounding term of the start corner, the factor of q in
     * the curvature law. The solve copies it from
     * S9CurveParameters::StartEta.
     */
    double StartEta = S9RoundingEta;


    /**
     * η of the rounding term of the end corner. The solve copies it
     * from S9CurveParameters::EndEta.
     */
    double EndEta = S9RoundingEta;


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
     * Solves the curve for given weights.
     *
     * @param parameters Angles, weights and rounding η of the curve. Both
     * angles must be finite, both weights from minus one to one, and both
     * η at least zero. It is only read during the call.
     *
     * Returns true when a curve that meets both tangent angles was found.
     */
    bool TrySolveWithWeights(const S9CurveParameters &parameters) noexcept;


    /**
     * Solves the curve for given weights from a nearby solved curve.
     *
     * @param parameters Angles, weights and rounding η of the curve, as
     * for TrySolveWithWeights. It is only read during the call.
     *
     * @param nearby Solved curve of nearby angles. It is only read. It can be
     * this object.
     *
     * Returns true when a curve that meets both tangent angles was found.
     */
    bool TrySolveWithWeightsNear(const S9CurveParameters &parameters,
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
     * Returns the turn J0 = 1/2 + λ0 that the start basis makes over the
     * whole curve.
     */
    constexpr double StartBasisTurn() const noexcept;


    /**
     * Returns the turn J1 = 1/2 + λ1 that the end basis makes over the whole
     * curve.
     */
    constexpr double EndBasisTurn() const noexcept;


    /**
     * Returns the turn in radians that the start corner term carries over
     * the whole curve. It is StartLine times the start corner λ.
     */
    constexpr double StartCornerTurn() const noexcept;


    /**
     * Returns the turn in radians that the end corner term carries over the
     * whole curve. It is EndLine times the end corner λ.
     */
    constexpr double EndCornerTurn() const noexcept;


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
     * Returns the integral of the start basis curvature from zero to t,
     * which is G0(t) = t - t^2 / 2 plus the turn of the start corner term
     * over the distance t.
     *
     * @param t Normalized arc length from zero through one.
     */
    double StartBasisIntegral(const double t) const noexcept;


    /**
     * Returns the integral of the end basis curvature from zero to t,
     * which is G1(t) = t^2 / 2 plus the turn of the end corner term over
     * the width t from the start.
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
     * Returns the derivative of the turn rate with respect to normalized arc
     * length. The end corner is measured backward from the end, thus its
     * slope changes sign. Needs a solved curve.
     *
     * @param t Normalized arc length from zero through one.
     */
    double TurnRateSlope(const double t) const noexcept;


    /**
     * Returns the signed curvature of the curve for a chord of length one.
     *
     * @param t Normalized arc length from zero through one.
     */
    double Curvature(const double t) const noexcept;


    /**
     * Returns the derivative of the signed curvature with respect to arc
     * length, for a chord of length one. It is the turn rate slope over the
     * square of the length. Needs a solved curve.
     *
     * @param t Normalized arc length from zero through one.
     */
    double CurvatureSlope(const double t) const noexcept;


    /**
     * Returns the change of the tangent angle in radians from the start to
     * t, which is StartLine * G0(t) + EndLine * G1(t) in closed form.
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
    return 0.5 + StartLambda;
}


inline constexpr double S9Curve::EndBasisTurn() const noexcept {
    return 0.5 + EndLambda;
}


inline constexpr double S9Curve::StartCornerTurn() const noexcept {
    return StartLine * StartLambda;
}


inline constexpr double S9Curve::EndCornerTurn() const noexcept {
    return EndLine * EndLambda;
}


inline double S9Curve::Curvature(const double t) const noexcept
{
    ASSERT(Length > 0.0);

    const double turnRate = TurnRate(t);

    return turnRate / Length;
}


inline double S9Curve::CurvatureSlope(const double t) const noexcept
{
    ASSERT(Length > 0.0);

    const double slope = TurnRateSlope(t);

    return slope / (Length * Length);
}


inline double S9Curve::TangentAngle(const double t) const noexcept
{
    const double turn = TurnFromStart(t);

    return StartAngle + turn;
}
