
#include "AngleRange.h"
#include "SplineCubicG2.h"
#include "TridiagonalSolver.h"
#include "UnitCubicCurve.h"


// G2-continuous cubic spline.
//
// Based on https://spline.technology/paper1.pdf
//
//
// A quick overview of terminology used here. Please read this before reading
// code.
//
// • Alpha. An angle on the left of segment. From -PI to PI. Measured in
//   radians.
//
// • Beta. An angle on the right side of segment. From -PI to PI. Measured in
//   radians.
//
// • Knot. User-supplied point through which spline is interpolated.
//
// • Knot angle. Absolute angle at knot. Measured in radians.
//
// • Segment. Piece of spline between two knots.
//
// • Segment angle. Absolute angle of a line segment between two knots of the
//   spline.
//
//
// Overview of the solver.
//
// Solver has a few steps. First step is to find initial guess for knot angles.
// This is a rough approximation calculated locally. Next step is to clamp knot
// angles to the middle of the allowed range for each knot. Since this
// specific spline is made by series of segments with two angle parameters,
// after increasing angle beyond certain amount the curve at the beginning or
// end reaches infinite curvature and flips to the other side. This flip point
// is a limit for angle rotation.
//
// Next step is a series of sweeps over the knots. Each sweep solves the angle
// at every knot in turn with the ITP root finding method, with angles at the
// neighbouring knots fixed, see SolveKnotAngle. The sweeps stop when a sweep
// changes no knot angle by more than MaximumITPError.
//
// A final step is a few Newton style iterations to finish the solution. A
// tridiagonal Jacobian matrix is constructed with analytically calculated
// partial derivatives for each knot along (the formulas of base curve allows
// that) with error vector (curvature differences at each knot). For acyclic
// case the end knots follow their neighbours through the end rules, so only
// interior knots are unknowns and end rules fold into the first and last row.
// Once Jacobian matrix is solved, knot angles are adjusted accordingly. Once
// the system is strictly within bounds at each knot, it takes about 3-4
// iterations to achieve maximum error of only 1e-16 radians which is far
// beyond any visible differences compared to even higher precision. When the
// Newton phase diverges, which shows as a growing step or a knot outside its
// rotation range, or does not converge within its iteration limit, the solver
// restores the angles from before the phase, tightens the sweep tolerance by
// FallbackITPErrorFactor and continues sweeps. The total number of sweeps
// stays within MaximumITPIterationCount. The partial derivatives come from
// CalculateSegmentCurvatures, which evaluates the curvatures at both ends of
// one segment and their derivatives with respect to both segment angles at
// once.


/**
 * Returns angle of a vector, in radians. Positive angle is counter-clock-wise.
 * Input vector does not have to be normalized.
 *
 * @param v Vector.
 */
static inline double VectorAngle(const FloatPoint v) {
    return atan2(-v.Y, v.X);
}


/**
 * Normalized a given angle in radians to -PI..PI range.
 */
static double NormalizeAngle(const double radians)
{
    static constexpr double PI2 = M_PI + M_PI;
    static constexpr double PI2_R = 1.0 / PI2;

    const double r = radians * PI2_R;

    return PI2 * (r - round(r));
}


/**
 * Converts absolute knot angle to alpha angle relative to segment.
 */
static double AbsoluteToAlpha(const double knotAngle,
    const double segmentAngle)
{
    ASSERT(segmentAngle >= -M_PI);
    ASSERT(segmentAngle <= M_PI);

    return knotAngle - segmentAngle;
}


/**
 * Converts alpha knot angle relative to segment to absolute angle.
 */
static double AlphaToAbsolute(const double knotAngle,
    const double segmentAngle)
{
    ASSERT(segmentAngle >= -M_PI);
    ASSERT(segmentAngle <= M_PI);

    return knotAngle + segmentAngle;
}


/**
 * Converts absolute knot angle to beta angle relative to segment.
 */
static double AbsoluteToBeta(const double knotAngle,
    const double segmentAngle)
{
    ASSERT(segmentAngle >= -M_PI);
    ASSERT(segmentAngle <= M_PI);

    return segmentAngle - knotAngle;
}


/**
 * Converts beta knot angle relative to segment to absolute angle.
 */
static double BetaToAbsolute(const double knotAngle,
    const double segmentAngle)
{
    ASSERT(segmentAngle >= -M_PI);
    ASSERT(segmentAngle <= M_PI);

    return -knotAngle + segmentAngle;
}


/**
 * Rotate point (1, 0) around coordinate system origina by a given angle in
 * radians.
 *
 * @param radians Rotation angle.
 */
static FloatPoint RotatePositiveXAxisUnitVector(const double radians)
{
    const double s = sin(radians);
    const double c = cos(radians);

    return FloatPoint(c, -s);
}


/**
 * Rotates a given point around coordinate system origin.
 */
static FloatPoint RotatePointAroundOrigin(const FloatPoint v,
    const double radians)
{
    const double s = Sin(radians);
    const double c = Cos(radians);

    const double x = v.X * c + v.Y * s;
    const double y = -v.X * s + v.Y * c;

    return FloatPoint(x, y);
}


/**
 * Single segment data. Keeps two angles used as parameters to curve generator
 * and chord length.
 */
struct SplineSegment final {
    SplineSegment(const double a0, const double a1);

    UnitCubicCurve GetUnitCurve() const;

    const double A0 = 0;
    const double A1 = 0;
    const double ArmLength0 = 0;
    const double ArmLength1 = 0;

private:

    SplineSegment() = delete;
};


/**
 * Calculates arm length at one end of unit cubic. The result is the distance
 * from that end to its control point as a fraction of the chord.
 */
static double CubicArmLength(const double a, const double b)
{
    const double x = a - 0.3 * sin(2.0 * b - 0.4 * sin(2.0 * b));

    // The arm length is (5 / 12) (cos(x) - 0.2 cos(3x)).
    //
    // The identity cos(3x) = 4 cos(x)^3 - 3 cos(x) turns it into
    // cos(x) (2 - cos(x)^2) / 3, which needs one cosine instead of two.
    const double c = cos(x);

    return c * (2.0 - c * c) / 3.0;
}


inline SplineSegment::SplineSegment(const double a0, const double a1)
:   A0(a0),
    A1(a1),
    ArmLength0(CubicArmLength(a0, a1)),
    ArmLength1(CubicArmLength(a1, a0))
{
}


inline UnitCubicCurve SplineSegment::GetUnitCurve() const
{
    const FloatPoint p0 = RotatePositiveXAxisUnitVector(A0) * ArmLength0;
    const FloatPoint p1 = RotatePositiveXAxisUnitVector(A1) * ArmLength1;

    return UnitCubicCurve(p0, FloatPoint(1.0 - p1.X, p1.Y));
}


/**
 * Curvatures at both ends of one segment and their partial derivatives
 * with respect to both angles of a segment.
 */
struct SegmentCurvatures final {
    constexpr SegmentCurvatures(const double start, const double end,
        const double startByAlpha, const double startByBeta,
        const double endByAlpha, const double endByBeta) noexcept;

    double Start = 0;
    double End = 0;
    double StartByAlpha = 0;
    double StartByBeta = 0;
    double EndByAlpha = 0;
    double EndByBeta = 0;

private:

    SegmentCurvatures() = delete;
};


inline constexpr SegmentCurvatures::SegmentCurvatures(const double start,
    const double end, const double startByAlpha, const double startByBeta,
    const double endByAlpha, const double endByBeta) noexcept
:   Start(start),
    End(end),
    StartByAlpha(startByAlpha),
    StartByBeta(startByBeta),
    EndByAlpha(endByAlpha),
    EndByBeta(endByBeta)
{
}


/**
 * Calculates curvatures at both ends of one segment and their partial
 * derivatives with respect to two angles of a segment. The result uses the
 * same curvature expressions as CalculateScaledCurvatureDelta, so the Newton
 * iterations get the residual and the Jacobian from one place.
 *
 * @param alpha Angle at the start of the segment relative to the segment
 * direction, in radians.
 *
 * @param beta Angle at the end of the segment relative to the segment
 * direction, in radians.
 *
 * @param length Segment length. It must be greater than zero.
 */
static SegmentCurvatures CalculateSegmentCurvatures(const double alpha,
    const double beta, const double length) noexcept
{
    ASSERT(length > 0.0);

    // Solver turns the tangent at each knot until the curvature on both sides
    // of the knot is the same. A Newton step needs, for each knot, curvature
    // difference and how fast this difference changes with angles of the knot
    // and of its two neighbours. Both come from segments. The curvature at
    // the start of a segment depends on alpha (the knot at the start) and on
    // beta (the knot at the end), and so does the curvature at corresponding
    // end. This function calculates these two curvatures and their
    // derivatives for a single segment.
    //
    // The arm lengths come from CubicArmLength. Each arm depends on the
    // angle at its own end directly and on the opposite angle through a
    // turn term t -
    //
    //     arm = L(own - t(opposite))
    //     t(b) = 0.3 sin(2b - 0.4 sin(2b))
    //     L(u) = cos(u) (2 - cos(u)^2) / 3
    //
    // then with the derivatives
    //
    //     t'(b) = 0.3 cos(2b - 0.4 sin(2b)) (2 - 0.8 cos(2b))
    //     L'(u) = -sin(u) (2 - 3 cos(u)^2) / 3

    const double turnAlphaArgument = 2.0 * alpha - 0.4 * sin(2.0 * alpha);
    const double turnBetaArgument = 2.0 * beta - 0.4 * sin(2.0 * beta);
    const double turnAlpha = 0.3 * sin(turnAlphaArgument);
    const double turnBeta = 0.3 * sin(turnBetaArgument);
    const double turnAlphaSlope =
        0.3 * cos(turnAlphaArgument) * (2.0 - 0.8 * cos(2.0 * alpha));
    const double turnBetaSlope =
        0.3 * cos(turnBetaArgument) * (2.0 - 0.8 * cos(2.0 * beta));

    // Start arm and end arm with their slopes L'. Derivative of an arm by its
    // own angle is L'. Derivative by the opposite angle is -L' t', because
    // the opposite angle enters only through the turn term.
    const double startCosine = cos(alpha - turnBeta);
    const double endCosine = cos(beta - turnAlpha);
    const double startArm =
        startCosine * (2.0 - startCosine * startCosine) / 3.0;
    const double endArm = endCosine * (2.0 - endCosine * endCosine) / 3.0;
    const double startArmSlope = -sin(alpha - turnBeta) *
        (2.0 - 3.0 * startCosine * startCosine) / 3.0;
    const double endArmSlope = -sin(beta - turnAlpha) *
        (2.0 - 3.0 * endCosine * endCosine) / 3.0;

    // Curvature at the start of a cubic is proportional to the cross product
    // of its first two control polygon edges divided by the cube of the first
    // edge length. For a segment with arms A and B this gives
    //
    //     start curvature ~ (B sin(alpha + beta) - sin(alpha)) / (A^2 length)
    //     end curvature   ~ (A sin(alpha + beta) - sin(beta)) / (B^2 length)
    //
    // Common factor is then left out. It cancels in the curvature difference
    // at knot, which is what the solver reduces to zero.

    const double sineSum = sin(alpha + beta);
    const double cosineSum = cos(alpha + beta);
    const double startNumerator = endArm * sineSum - sin(alpha);
    const double endNumerator = startArm * sineSum - sin(beta);
    const double startScale = 1.0 / (startArm * startArm * length);
    const double endScale = 1.0 / (endArm * endArm * length);

    // Derivatives by the quotient rule. For n / (arm^2 length), derivative by
    // an angle is (n' - 2 n arm' / arm) / (arm^2 length), where n' and arm'
    // are derivatives by that angle.

    const double startQuotient =
        2.0 * startNumerator * startArmSlope / startArm;
    const double endQuotient = 2.0 * endNumerator * endArmSlope / endArm;

    const double start = startNumerator * startScale;
    const double end = endNumerator * endScale;
    const double startByAlpha = (endArm * cosineSum - cos(alpha) -
        endArmSlope * turnAlphaSlope * sineSum - startQuotient) * startScale;
    const double startByBeta = (endArm * cosineSum + endArmSlope * sineSum +
        startQuotient * turnBetaSlope) * startScale;
    const double endByAlpha = (startArm * cosineSum + startArmSlope * sineSum +
        endQuotient * turnAlphaSlope) * endScale;
    const double endByBeta = (startArm * cosineSum - cos(beta) -
        startArmSlope * turnBetaSlope * sineSum - endQuotient) * endScale;

    return SegmentCurvatures(start, end, startByAlpha, startByBeta,
        endByAlpha, endByBeta);
}


/**
 * One row of the Newton system of a knot. Jacobian elements for the angles at
 * the previous, center and the next knot, and the residual.
 */
struct KnotRow final {
    constexpr KnotRow(const double subdiagonal, const double diagonal,
        const double superdiagonal, const double residual) noexcept;

    double Subdiagonal = 0;
    double Diagonal = 0;
    double Superdiagonal = 0;
    double Residual = 0;

private:

    KnotRow() = delete;
};


inline constexpr KnotRow::KnotRow(const double subdiagonal,
    const double diagonal, const double superdiagonal,
    const double residual) noexcept
:   Subdiagonal(subdiagonal),
    Diagonal(diagonal),
    Superdiagonal(superdiagonal),
    Residual(residual)
{
}


/**
 * Calculates Newton row of a knot from the curvatures of the segment that
 * ends at the knot and the segment that starts at it. The residual is
 * curvature difference at the knot. The alpha angle of segment moves with the
 * knot at its start and the beta angle moves against the knot at its end,
 * which sets signs of Jacobian elements.
 *
 * @param previous Curvatures of the segment that ends at the knot.
 *
 * @param next Curvatures of the segment that starts at the knot.
 */
static KnotRow CalculateKnotRow(const SegmentCurvatures &previous,
    const SegmentCurvatures &next) noexcept
{
    const double subdiagonal = -previous.EndByAlpha;
    const double diagonal = next.StartByAlpha + previous.EndByBeta;
    const double superdiagonal = -next.StartByBeta;
    const double residual = next.Start - previous.End;

    return KnotRow(subdiagonal, diagonal, superdiagonal, residual);
}


/**
 * Calculates rotation limit axis given opposite angle.
 */
static double CalculateRotationLimitAxis(const double beta)
{
    const double beta2 = beta + beta;

    // Given opposite angle, calculate angle at which arm length becomes zero.
    // This will be the angle beyond which curve will flip.
    return 0.3 * sin(beta2 - (0.4 * sin(beta2)));
}


static AngleRange CalculateRotationRangeForSegmentPair(
    const double prevKnotAngle, const double nextKnotAngle,
    const double prevSegmentAngle, const double nextSegmentAngle)
{
    // Rotation limit axis for beta angle of the previous segment.
    const double betaLimit = BetaToAbsolute(CalculateRotationLimitAxis(
        AbsoluteToAlpha(prevKnotAngle, prevSegmentAngle)), prevSegmentAngle);

    // Rotation limit axis for alpha angle of the next segment.
    const double alphaLimit = AlphaToAbsolute(CalculateRotationLimitAxis(
        AbsoluteToBeta(nextKnotAngle, nextSegmentAngle)), nextSegmentAngle);

    const double axis = (betaLimit + alphaLimit) / 2.0;
    const double range = (M_PI - Abs(NormalizeAngle(betaLimit - alphaLimit))) / 2.0;

    return AngleRange(axis - range, axis + range);
}


/**
 * Calculates curvature delta at the shared point of two neighbouring
 * segments, scaled by denominators of both curvatures. Curvature at the start
 * of the next segment is a quotient n1 / d1 and the curvature at the end of
 * the previous segment is n0 / d0, see CalculateSegmentCurvatures. Result is
 * n1 d0 - n0 d1. Inside the rotation range of the center knot both d0 and d1
 * are positive, so result has the sign and the root of curvature delta, but
 * it has no poles at the ends of the range. This keeps the of the root finder
 * effective.
 *
 * @param x Absolute angle in radians at the previous knot.
 *
 * @param y Absolute angle in radians at the center knot. This is the knot for
 * which the curvature delta is calculated.
 *
 * @param z Absolute angle in radians at the next knot.
 *
 * @param a0 Absolute angle in radians of the previous segment, from -PI to PI.
 *
 * @param a1 Absolute angle in radians of the next segment, from -PI to PI.
 *
 * @param h0 Previous segment length. It must be greater than zero.
 *
 * @param h1 Next segment length. It must be greater than zero.
 */
static double CalculateScaledCurvatureDelta(const double x,
    const double y, const double z, const double a0,
    const double a1, const double h0, const double h1) noexcept
{
    ASSERT(h0 > 0.0);
    ASSERT(h1 > 0.0);

    // Alpha/beta angles for both segments.
    const double alpha0 = AbsoluteToAlpha(x, a0);
    const double beta0 = AbsoluteToBeta(y, a0);
    const double alpha1 = AbsoluteToAlpha(y, a1);
    const double beta1 = AbsoluteToBeta(z, a1);

    // Arm lengths.
    const double alphaArm0 = CubicArmLength(alpha0, beta0);
    const double betaArm0 = CubicArmLength(beta0, alpha0);
    const double alphaArm1 = CubicArmLength(alpha1, beta1);
    const double betaArm1 = CubicArmLength(beta1, alpha1);

    const double numerator1 =
        (betaArm1 * h1) * sin(alpha1 + beta1) - h1 * sin(alpha1);
    const double denominator1 = Square(alphaArm1 * h1);

    const double numerator0 =
        (alphaArm0 * h0) * sin(alpha0 + beta0) - h0 * sin(beta0);
    const double denominator0 = Square(betaArm0 * h0);

    return numerator1 * denominator0 - numerator0 * denominator1;
}


static double ClampKnotAngle(const double prevKnotAngle,
    const double currentKnotAngle, const double nextKnotAngle,
    const double prevSegmentAngle, const double nextSegmentAngle)
{
    const AngleRange absoluteRotationRange =
        CalculateRotationRangeForSegmentPair(prevKnotAngle,
            nextKnotAngle, prevSegmentAngle, nextSegmentAngle);

    const double min = absoluteRotationRange.MinAngle;
    const double max = absoluteRotationRange.MaxAngle;
    const double sum = min + max;

    return sum * 0.5;
}


/*
 * Number of times the root finder halves the rotation range in the worst
 * case. Angle it returns is within the range divided by 2^32 of the root, the
 same as after 32 bisection steps.
 */
static constexpr U32 ItpHalvingCount = 32U;


/*
 * Weight of the truncation step of the root finder, relative to the width
 * of the rotation range. See SolveKnotAngle.
 */
static constexpr double ItpTruncationWeight = 0.2;


/**
 * Solves the angle at the center knot at which the curvatures of the two
 * neighbouring segments match, with angles at the neighbouring knots fixed.
 * The scaled curvature delta decreases strictly over rotation range of the
 * knot and has poles of opposite sign at its ends, so the range brackets
 * exactly one root.
 *
 * The root is found with the ITP method of Oliveira and Takahashi (2020),
 * see https://en.wikipedia.org/wiki/ITP_method.
 *
 * @param prevKnotAngle Absolute angle in radians at the previous knot.
 *
 * @param nextKnotAngle Absolute angle in radians at the next knot.
 *
 * @param prevSegmentAngle Absolute angle in radians of the previous segment,
 * from -PI to PI.
 *
 * @param nextSegmentAngle Absolute angle in radians of the next segment, from
 * -PI to PI.
 *
 * @param prevSegmentLength Length of the previous segment. It must be greater
 * than zero.
 *
 * @param nextSegmentLength Length of the next segment. It must be greater than
 * zero.
 *
 * @return Absolute angle in radians at the center knot. If scaled curvature
 * delta does not change sign over the rotation range, the upper end of the
 * range is returned.
 */
static double SolveKnotAngle(const double prevKnotAngle,
    const double nextKnotAngle, const double prevSegmentAngle,
    const double nextSegmentAngle, const double prevSegmentLength,
    const double nextSegmentLength) noexcept
{
    ASSERT(prevSegmentLength > 0.0);
    ASSERT(nextSegmentLength > 0.0);

    const AngleRange absoluteRotationRange =
        CalculateRotationRangeForSegmentPair(prevKnotAngle,
            nextKnotAngle, prevSegmentAngle, nextSegmentAngle);

    double min = absoluteRotationRange.MinAngle;
    double max = absoluteRotationRange.MaxAngle;

    double minValue = CalculateScaledCurvatureDelta(prevKnotAngle, min,
        nextKnotAngle, prevSegmentAngle, nextSegmentAngle,
        prevSegmentLength, nextSegmentLength);

    double maxValue = CalculateScaledCurvatureDelta(prevKnotAngle, max,
        nextKnotAngle, prevSegmentAngle, nextSegmentAngle,
        prevSegmentLength, nextSegmentLength);

    ASSERT(NumberIsFinite(minValue));
    ASSERT(NumberIsFinite(maxValue));

    if (minValue == 0.0) {
        return min;
    }

    if (maxValue == 0.0) {
        return max;
    }

    // Without a sign change there is no root to find. Bisection would end
    // at the upper end of the range in this case.
    if ((minValue < 0.0) == (maxValue < 0.0)) {
        return max;
    }

    // The steps below expect a negative value at min and a positive value
    // at max. Flip the sign of all values otherwise.
    double sign = 1.0;

    if (minValue > 0.0) {
        sign = -1.0;
    }

    minValue *= sign;
    maxValue *= sign;

    const double range = max - min;
    const double epsilon =
        ldexp(range, -static_cast<int>(ItpHalvingCount + 1U));
    const double truncation = ItpTruncationWeight / range;
    const U32 maximumStepCount = ItpHalvingCount + 1U;

    for (U32 step = 0; step < maximumStepCount; step++) {
        const double width = max - min;

        if (width <= 2.0 * epsilon) {
            break;
        }

        const double middle = 0.5 * (min + max);

        // Regula falsi estimate.
        //
        // It is inside the bracket because the values at the ends have
        // opposite signs.
        const double secant =
            (maxValue * min - minValue * max) / (maxValue - minValue);

        // Truncation.
        //
        // Move estimate towards the middle by offset that shrinks with the
        // square of the bracket width, but never move it past the middle.
        const double offset = truncation * width * width;

        double direction = 1.0;

        if (secant > middle) {
            direction = -1.0;
        }

        double truncated = middle;

        if (offset <= Abs(middle - secant)) {
            truncated = secant + direction * offset;
        }

        // Projection.
        //
        // Keep the probe inside a window around the middle whose radius
        // halves with each step.
        const int remaining = static_cast<int>(maximumStepCount - step);

        double radius = ldexp(epsilon, remaining) - 0.5 * width;

        if (radius < 0.0) {
            radius = 0.0;
        }

        double probe = middle - direction * radius;

        if (Abs(truncated - middle) <= radius) {
            probe = truncated;
        }

        const double value = sign * CalculateScaledCurvatureDelta(
            prevKnotAngle, probe, nextKnotAngle, prevSegmentAngle,
            nextSegmentAngle, prevSegmentLength, nextSegmentLength);

        if (value > 0.0) {
            max = probe;
            maxValue = value;
        } else if (value < 0.0) {
            min = probe;
            minValue = value;
        } else {
            return probe;
        }
    }

    return 0.5 * (min + max);
}


static double AcyclicStartAngle(const double firstSegmentAngle,
    const double secondKnotAngle)
{
    const double beta = AbsoluteToBeta(secondKnotAngle, firstSegmentAngle);
    const double alpha = 0.5 * Sin(2.0 * beta);

    return AlphaToAbsolute(alpha, firstSegmentAngle);
}


static double AcyclicEndAngle(const double finalSegmentAngle,
    const double firstKnotAngle)
{
    const double alpha = AbsoluteToAlpha(firstKnotAngle, finalSegmentAngle);
    const double beta = 0.5 * Sin(2.0 * alpha);

    return BetaToAbsolute(beta, finalSegmentAngle);
}


static constexpr U32 MaximumCenteringIterationCount = 3;

static constexpr U32 MaximumITPIterationCount = 20;

static constexpr U32 MaximumNewtonIterationCount = 12;


/*
 * Sweep that changes no knot angle by more than this, in radians, ends the
 * sweeps and starts the Newton phase.
 */
static constexpr double MaximumITPError = 1e-2;


/*
 * Sweep tolerance is multiplied by this factor after failed Newton phase, see
 * SolveCubicG2Acyclic.
 */
static constexpr double FallbackITPErrorFactor = 0.01;

static constexpr double MaximumNewtonError = 1e-15;


void SolveCubicG2Acyclic(const U32 count, const FloatPoint * const knots,
    CurveArray<CubicCurve> &curves) noexcept
{
    const U32 segmentCount = count - 1;

    double *segmentLengths = AllocElements<double>(segmentCount);
    double *segmentAngles = AllocElements<double>(segmentCount);

    for (U32 i = 0; i < segmentCount; i++) {
        const FloatPoint current = knots[i];
        const FloatPoint next = knots[i + 1];
        const FloatPoint d = next - current;

        segmentLengths[i] = d.Length();

        segmentAngles[i] = VectorAngle(d);
    }

    const U32 n = count - 1;

    double *knotAngles = AllocElements<double>(count);

    for (U32 i = 1; i < n; i++) {
        const int prevIndex = i - 1;
        const int nextIndex = i + 1;

        const FloatPoint current = knots[i];

        const FloatPoint prevToCurrent = current - knots[prevIndex];
        const FloatPoint currentToNext = knots[nextIndex] - current;

        const double a0 = VectorAngle(prevToCurrent);
        const double a1 = VectorAngle(currentToNext);
        const double l0 = prevToCurrent.Length();
        const double l1 = currentToNext.Length();
        const double combinedLengths = l0 + l1;

        ASSERT(combinedLengths > DBL_EPSILON);

        const double bend = a1 - a0;

        knotAngles[i] = a0 + ((bend * l0) / combinedLengths);
    }

    knotAngles[0] = AlphaToAbsolute(0, segmentAngles[0]);

    knotAngles[n] = BetaToAbsolute(0, segmentAngles[segmentCount - 1]);

    // Centering iterations.
    for (U32 i = 0; i < MaximumCenteringIterationCount; i++) {
        const double angle0 = AcyclicStartAngle(segmentAngles[0], knotAngles[1]);

        knotAngles[0] = angle0;

        for (U32 j = 1; j < n; j++) {
            const double prevKnotAngle = knotAngles[j - 1];
            const double currentKnotAngle = knotAngles[j];
            const double nextKnotAngle = knotAngles[j + 1];
            const double prevSegmentAngle = segmentAngles[j - 1];
            const double nextSegmentAngle = segmentAngles[j];

            knotAngles[j] = ClampKnotAngle(prevKnotAngle, currentKnotAngle,
                nextKnotAngle, prevSegmentAngle, nextSegmentAngle);
        }

        const double angle1 = AcyclicEndAngle(segmentAngles[n - 1], knotAngles[n - 1]);

        knotAngles[n] = angle1;
    }

    // Tridiagonal Jacobian matrix elements and the residual vector of
    // interior knots. Two knots have no interior knot, so sweeps alone solve
    // them and run to the Newton tolerance.
    const U32 interiorCount = count - 2U;

    double *a = AllocElements<double>(count);
    double *b = AllocElements<double>(count);
    double *c = AllocElements<double>(count);
    double *x = AllocElements<double>(count);

    // Curvatures of all segments with their partial derivatives.
    SegmentCurvatures *segments =
        AllocElements<SegmentCurvatures>(segmentCount);

    // TKnot angles saved before a Newton phase, to restore them when the
    // phase fails.
    double *savedKnotAngles = AllocElements<double>(count);

    double tolerance = MaximumITPError;

    if (interiorCount == 0U) {
        tolerance = MaximumNewtonError;
    }

    U32 sweepCount = 0;
    bool converged = false;

    // Rounds of sweeps and a Newton phase. Sweeps stop when sweep changes no
    // knot angle by more than the tolerance. When the Newton phase then
    // fails, its start is restored, tolerance is tightened and sweeps
    // continue, as long as sweep limit allows.
    while (!converged) {
        double sweepError = 0;

        // Sweeps of ITP root finding, one knot at a time.
        do {
            const double angle0 =
                AcyclicStartAngle(segmentAngles[0], knotAngles[1]);

            double maxError = Abs(knotAngles[0] - angle0);

            knotAngles[0] = angle0;

            for (U32 j = 1; j < n; j++) {
                const double prevKnotAngle = knotAngles[j - 1];
                const double nextKnotAngle = knotAngles[j + 1];
                const double prevSegmentAngle = segmentAngles[j - 1];
                const double nextSegmentAngle = segmentAngles[j];
                const double prevSegmentLength = segmentLengths[j - 1];
                const double nextSegmentLength = segmentLengths[j];

                const double angle = SolveKnotAngle(prevKnotAngle,
                    nextKnotAngle, prevSegmentAngle, nextSegmentAngle,
                    prevSegmentLength, nextSegmentLength);

                maxError = Max(maxError, Abs(knotAngles[j] - angle));

                knotAngles[j] = angle;
            }

            const double angle1 =
                AcyclicEndAngle(segmentAngles[n - 1], knotAngles[n - 1]);

            maxError = Max(maxError, Abs(knotAngles[n] - angle1));

            knotAngles[n] = angle1;

            sweepCount++;
            sweepError = maxError;
        } while (sweepError >= tolerance and
            sweepCount < MaximumITPIterationCount);

        if (interiorCount == 0U) {
            break;
        }

        // Sweep sets the first knot from the second one before it solves the
        // second one, so the start rule is one sweep behind. Apply it again,
        // so the Newton phase starts from a state that satisfies both end rules.
        knotAngles[0] = AcyclicStartAngle(segmentAngles[0], knotAngles[1]);

        for (U32 j = 0; j < count; j++) {
            savedKnotAngles[j] = knotAngles[j];
        }

        // Newton iterations for interior knots.
        //
        //
        // End knots follow their neighbours through AcyclicStartAngle and
        // AcyclicEndAngle, so their derivatives fold into the first and the
        // last row of Jacobian -
        //
        //     dy0 / dy1 = -cos(2 beta0)
        //     dyn / dy(n - 1) = -cos(2 alpha(n - 1))
        //
        // where beta0 is the beta angle of the first segment and
        // alpha(n - 1) is the alpha angle of the last segment. A step that
        // grows, or a knot pushed outside its rotation range, means that the
        // phase diverges.

        double previousStep = M_PI;
        bool diverged = false;

        for (U32 i = 0; i < MaximumNewtonIterationCount; i++) {
            for (U32 j = 0; j < segmentCount; j++) {
                const double segmentAngle = segmentAngles[j];
                const double alpha =
                    AbsoluteToAlpha(knotAngles[j], segmentAngle);
                const double beta =
                    AbsoluteToBeta(knotAngles[j + 1], segmentAngle);

                segments[j] = CalculateSegmentCurvatures(alpha, beta,
                    segmentLengths[j]);
            }

            for (U32 j = 1; j < n; j++) {
                const KnotRow row =
                    CalculateKnotRow(segments[j - 1], segments[j]);

                a[j - 1] = row.Subdiagonal;
                b[j - 1] = row.Diagonal;
                c[j - 1] = row.Superdiagonal;
                x[j - 1] = row.Residual;
            }

            const double startBeta =
                AbsoluteToBeta(knotAngles[1], segmentAngles[0]);
            const double endAlpha =
                AbsoluteToAlpha(knotAngles[n - 1], segmentAngles[n - 1]);

            b[0] -= a[0] * cos(2.0 * startBeta);
            b[interiorCount - 1U] -=
                c[interiorCount - 1U] * cos(2.0 * endAlpha);

            // Tridiagonal solver needs at least two equations.
            if (interiorCount == 1U) {
                x[0] = x[0] / b[0];
            } else {
                SolveTridiagonal(interiorCount, a, b, c, x, x);
            }

            double maxStep = 0;

            for (U32 j = 1; j < n; j++) {
                maxStep = Max(maxStep, Abs(x[j - 1]));
            }

            if (maxStep > previousStep) {
                diverged = true;
                break;
            }

            previousStep = maxStep;

            // Update interior knot angles, then the end knots.
            for (U32 j = 1; j < n; j++) {
                knotAngles[j] -= x[j - 1];
            }

            knotAngles[0] = AcyclicStartAngle(segmentAngles[0], knotAngles[1]);
            knotAngles[n] =
                AcyclicEndAngle(segmentAngles[n - 1], knotAngles[n - 1]);

            for (U32 j = 1; j < n; j++) {
                const AngleRange range =
                    CalculateRotationRangeForSegmentPair(knotAngles[j - 1],
                        knotAngles[j + 1], segmentAngles[j - 1],
                        segmentAngles[j]);

                if (knotAngles[j] < range.MinAngle or
                    knotAngles[j] > range.MaxAngle) {
                    diverged = true;
                }
            }

            if (diverged) {
                break;
            }

            if (maxStep <= MaximumNewtonError) {
                converged = true;
                break;
            }
        }

        if (converged) {
            break;
        }

        // Newton phase shit the bed. Restore its start. More sweeps with a
        // tighter tolerance give the next phase a better start.
        for (U32 j = 0; j < count; j++) {
            knotAngles[j] = savedKnotAngles[j];
        }

        if (sweepCount >= MaximumITPIterationCount) {
            break;
        }

        tolerance *= FallbackITPErrorFactor;
    }

    // Produce spline.
    curves.Restart(knots[0]);

    for (U32 i = 0; i < n; i++) {
        const FloatPoint p0 = knots[i];
        const FloatPoint p1 = knots[i + 1];

        const double segmentAngle = segmentAngles[i];

        const double a0 = AbsoluteToAlpha(knotAngles[i], segmentAngle);
        const double a1 = AbsoluteToBeta(knotAngles[i + 1],
            segmentAngle);

        const SplineSegment segment(a0, a1);

        const double segmentLength = segmentLengths[i];
        const UnitCubicCurve curve = segment.GetUnitCurve();

        // Control points, unrotated.
        const FloatPoint c0u(curve.P0.X * segmentLength,
            curve.P0.Y * segmentLength);

        const FloatPoint c1u(curve.P1.X * segmentLength,
            curve.P1.Y * segmentLength);

        const FloatPoint rotatedFirstControl = RotatePointAroundOrigin(c0u,
            segmentAngle);
        const FloatPoint rotatedSecondControl = RotatePointAroundOrigin(c1u,
            segmentAngle);
        curves.CurveTo(p0 + rotatedFirstControl, p0 + rotatedSecondControl, p1);
    }

    Dealloc(savedKnotAngles);
    Dealloc(segments);
    Dealloc(x);
    Dealloc(c);
    Dealloc(b);
    Dealloc(a);
    Dealloc(knotAngles);
    Dealloc(segmentAngles);
    Dealloc(segmentLengths);
}


void SolveCubicG2Cyclic(const U32 count, const FloatPoint * const knots,
    CurveArray<CubicCurve> &curves) noexcept
{
    // See SolveCubicG2Acyclic for more comments.
    const U32 firstKnotIndex = 0;
    const U32 firstSegmentIndex = 0;
    const U32 finalKnotIndex = count - 1;
    const U32 finalSegmentIndex = count - 1;

    double *segmentLengths = AllocElements<double>(count);
    double *segmentAngles = AllocElements<double>(count);

    // Calculate segment properties for all segments except the last one.
    for (U32 i = 0; i < count - 1; i++) {
        const FloatPoint current = knots[i];
        const FloatPoint next = knots[i + 1];
        const FloatPoint d = next - current;

        segmentLengths[i] = d.Length();

        segmentAngles[i] = VectorAngle(d);
    }

    // Calculate segment properties for the last segment which goes from the
    // last knot to the first one.
    const FloatPoint current = knots[finalKnotIndex];
    const FloatPoint next = knots[firstKnotIndex];
    const FloatPoint d = next - current;

    segmentLengths[finalSegmentIndex] = d.Length();

    segmentAngles[finalSegmentIndex] = VectorAngle(d);

    // Calculate knot properties.
    double *knotAngles = AllocElements<double>(count);

    // For all knots except the first one and the last one.
    for (U32 i = 1; i < count - 1; i++) {
        const int prevIndex = i - 1;
        const int nextIndex = i + 1;

        const FloatPoint current = knots[i];

        const FloatPoint prevToCurrent = current - knots[prevIndex];
        const FloatPoint currentToNext = knots[nextIndex] - current;

        const double a0 = VectorAngle(prevToCurrent);
        const double a1 = VectorAngle(currentToNext);
        const double l0 = prevToCurrent.Length();
        const double l1 = currentToNext.Length();
        const double combinedLengths = l0 + l1;

        ASSERT(combinedLengths > DBL_EPSILON);

        const double bend = a1 - a0;

        knotAngles[i] = a0 + ((bend * l0) / combinedLengths);
    }

    // First and last knots.
    knotAngles[firstKnotIndex] = AlphaToAbsolute(0, segmentAngles[firstSegmentIndex]);

    knotAngles[finalKnotIndex] = BetaToAbsolute(0, segmentAngles[finalSegmentIndex]);

    // Centering iterations.
    for (U32 i = 0; i < MaximumCenteringIterationCount; i++) {
        knotAngles[firstKnotIndex] = ClampKnotAngle(
            knotAngles[finalKnotIndex], knotAngles[firstKnotIndex], knotAngles[1],
            segmentAngles[finalSegmentIndex], segmentAngles[firstSegmentIndex]);

        for (U32 j = 1; j < count - 1; j++) {
            const double prevKnotAngle = knotAngles[j - 1];
            const double currentKnotAngle = knotAngles[j];
            const double nextKnotAngle = knotAngles[j + 1];
            const double prevSegmentAngle = segmentAngles[j - 1];
            const double nextSegmentAngle = segmentAngles[j];

            knotAngles[j] = ClampKnotAngle(prevKnotAngle, currentKnotAngle,
                nextKnotAngle, prevSegmentAngle, nextSegmentAngle);
        }

        knotAngles[finalKnotIndex] = ClampKnotAngle(
            knotAngles[finalKnotIndex - 1], knotAngles[finalKnotIndex],
            knotAngles[firstKnotIndex],
            segmentAngles[finalSegmentIndex - 1],
            segmentAngles[finalSegmentIndex]);
    }

    double *a = AllocElements<double>(count);
    double *b = AllocElements<double>(count);
    double *c = AllocElements<double>(count);
    double *x = AllocElements<double>(count);

    SegmentCurvatures *segments = AllocElements<SegmentCurvatures>(count);

    double *savedKnotAngles = AllocElements<double>(count);

    double tolerance = MaximumITPError;
    U32 sweepCount = 0;
    bool converged = false;

    while (not converged) {
        double sweepError = 0;

        do {
            const double angle0 = SolveKnotAngle(knotAngles[finalKnotIndex],
                knotAngles[1], segmentAngles[finalSegmentIndex],
                segmentAngles[firstSegmentIndex],
                segmentLengths[finalSegmentIndex],
                segmentLengths[firstSegmentIndex]);

            double maxError = Abs(knotAngles[firstKnotIndex] - angle0);

            knotAngles[0] = angle0;

            for (U32 j = 1; j < count - 1; j++) {
                const double prevKnotAngle = knotAngles[j - 1];
                const double nextKnotAngle = knotAngles[j + 1];
                const double prevSegmentAngle = segmentAngles[j - 1];
                const double nextSegmentAngle = segmentAngles[j];
                const double prevSegmentLength = segmentLengths[j - 1];
                const double nextSegmentLength = segmentLengths[j];

                const double angle = SolveKnotAngle(prevKnotAngle,
                    nextKnotAngle, prevSegmentAngle, nextSegmentAngle,
                    prevSegmentLength, nextSegmentLength);

                maxError = Max(maxError, Abs(knotAngles[j] - angle));

                knotAngles[j] = angle;
            }

            const double angle1 = SolveKnotAngle(
                knotAngles[finalKnotIndex - 1], knotAngles[firstKnotIndex],
                segmentAngles[finalSegmentIndex - 1],
                segmentAngles[finalSegmentIndex],
                segmentLengths[finalSegmentIndex - 1],
                segmentLengths[finalSegmentIndex]);

            maxError = Max(maxError, Abs(knotAngles[finalKnotIndex] - angle1));

            knotAngles[finalKnotIndex] = angle1;

            sweepCount++;
            sweepError = maxError;
        } while (sweepError >= tolerance and
            sweepCount < MaximumITPIterationCount);

        for (U32 j = 0; j < count; j++) {
            savedKnotAngles[j] = knotAngles[j];
        }

        double previousStep = M_PI;
        bool diverged = false;

        for (U32 i = 0; i < MaximumNewtonIterationCount; i++) {
            for (U32 j = 0; j < finalSegmentIndex; j++) {
                const double segmentAngle = segmentAngles[j];
                const double alpha =
                    AbsoluteToAlpha(knotAngles[j], segmentAngle);
                const double beta =
                    AbsoluteToBeta(knotAngles[j + 1], segmentAngle);

                segments[j] = CalculateSegmentCurvatures(alpha, beta,
                    segmentLengths[j]);
            }

            const double finalSegmentAngle = segmentAngles[finalSegmentIndex];
            const double finalAlpha =
                AbsoluteToAlpha(knotAngles[finalKnotIndex], finalSegmentAngle);
            const double finalBeta =
                AbsoluteToBeta(knotAngles[firstKnotIndex], finalSegmentAngle);

            segments[finalSegmentIndex] = CalculateSegmentCurvatures(
                finalAlpha, finalBeta, segmentLengths[finalSegmentIndex]);

            const KnotRow firstRow = CalculateKnotRow(
                segments[finalSegmentIndex], segments[firstSegmentIndex]);

            a[0] = firstRow.Subdiagonal;
            b[0] = firstRow.Diagonal;
            c[0] = firstRow.Superdiagonal;
            x[0] = firstRow.Residual;

            for (U32 j = 1; j < count; j++) {
                const KnotRow row =
                    CalculateKnotRow(segments[j - 1], segments[j]);

                a[j] = row.Subdiagonal;
                b[j] = row.Diagonal;
                c[j] = row.Superdiagonal;
                x[j] = row.Residual;
            }

            SolveTridiagonalCyclic(count, a, b, c, x, x);

            double maxStep = 0;

            for (U32 j = 0; j < count; j++) {
                maxStep = Max(maxStep, Abs(x[j]));
            }

            if (maxStep > previousStep) {
                diverged = true;
                break;
            }

            previousStep = maxStep;

            for (U32 j = 0; j < count; j++) {
                knotAngles[j] -= x[j];
            }

            const AngleRange firstRange = CalculateRotationRangeForSegmentPair(
                knotAngles[finalKnotIndex], knotAngles[1],
                segmentAngles[finalSegmentIndex],
                segmentAngles[firstSegmentIndex]);

            if (knotAngles[0] < firstRange.MinAngle or
                knotAngles[0] > firstRange.MaxAngle) {
                diverged = true;
            }

            for (U32 j = 1; j < count - 1; j++) {
                const AngleRange range =
                    CalculateRotationRangeForSegmentPair(knotAngles[j - 1],
                        knotAngles[j + 1], segmentAngles[j - 1],
                        segmentAngles[j]);

                if (knotAngles[j] < range.MinAngle or
                    knotAngles[j] > range.MaxAngle) {
                    diverged = true;
                }
            }

            const AngleRange lastRange = CalculateRotationRangeForSegmentPair(
                knotAngles[finalKnotIndex - 1], knotAngles[firstKnotIndex],
                segmentAngles[finalSegmentIndex - 1],
                segmentAngles[finalSegmentIndex]);

            if (knotAngles[finalKnotIndex] < lastRange.MinAngle or
                knotAngles[finalKnotIndex] > lastRange.MaxAngle) {
                diverged = true;
            }

            if (diverged) {
                break;
            }

            if (maxStep <= MaximumNewtonError) {
                converged = true;
                break;
            }
        }

        if (converged) {
            break;
        }

        for (U32 j = 0; j < count; j++) {
            knotAngles[j] = savedKnotAngles[j];
        }

        if (sweepCount >= MaximumITPIterationCount) {
            break;
        }

        tolerance *= FallbackITPErrorFactor;
    }

    // Produce spline.
    curves.Restart(knots[firstKnotIndex]);

    for (U32 i = 0; i < count - 1; i++) {
        const FloatPoint p0 = knots[i];
        const FloatPoint p1 = knots[i + 1];

        const double segmentAngle = segmentAngles[i];

        const double a0 = AbsoluteToAlpha(knotAngles[i], segmentAngle);
        const double a1 = AbsoluteToBeta(knotAngles[i + 1],
            segmentAngle);

        const SplineSegment segment(a0, a1);

        const double segmentLength = segmentLengths[i];
        const UnitCubicCurve curve = segment.GetUnitCurve();

        // Control points, unrotated.
        const FloatPoint c0u(curve.P0.X * segmentLength,
            curve.P0.Y * segmentLength);

        const FloatPoint c1u(curve.P1.X * segmentLength,
            curve.P1.Y * segmentLength);

        const FloatPoint rotatedFirstControl = RotatePointAroundOrigin(c0u,
            segmentAngle);
        const FloatPoint rotatedSecondControl = RotatePointAroundOrigin(c1u,
            segmentAngle);
        curves.CurveTo(p0 + rotatedFirstControl, p0 + rotatedSecondControl, p1);
    }

    // Close.
    const FloatPoint p0 = knots[finalKnotIndex];
    const FloatPoint p1 = knots[firstKnotIndex];

    const double segmentAngle = segmentAngles[finalSegmentIndex];

    const double a0 = AbsoluteToAlpha(knotAngles[finalKnotIndex], segmentAngle);
    const double a1 = AbsoluteToBeta(knotAngles[firstKnotIndex],
        segmentAngle);

    const SplineSegment segment(a0, a1);

    const double segmentLength = segmentLengths[finalSegmentIndex];
    const UnitCubicCurve curve = segment.GetUnitCurve();

    // Control points, unrotated.
    const FloatPoint c0u(curve.P0.X * segmentLength,
        curve.P0.Y * segmentLength);

    const FloatPoint c1u(curve.P1.X * segmentLength,
        curve.P1.Y * segmentLength);

    const FloatPoint rotatedFirstControl = RotatePointAroundOrigin(c0u,
        segmentAngle);
    const FloatPoint rotatedSecondControl = RotatePointAroundOrigin(c1u,
        segmentAngle);
    curves.CurveTo(p0 + rotatedFirstControl, p0 + rotatedSecondControl, p1);

    Dealloc(savedKnotAngles);
    Dealloc(segments);
    Dealloc(x);
    Dealloc(c);
    Dealloc(b);
    Dealloc(a);
    Dealloc(knotAngles);
    Dealloc(segmentAngles);
    Dealloc(segmentLengths);
}
