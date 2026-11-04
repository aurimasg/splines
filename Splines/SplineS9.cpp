
#include "AngleRange.h"
#include "S9Curve.h"
#include "SplineS9.h"
#include "TridiagonalSolver.h"


// Conceptually this solver is similar to CubicG2, but uses non-polynomial
// curves as building blocks.


/**
 * Returns angle of a vector, in radians. Positive angle is counter-clock-wise.
 * Input vector does not have to be normalized.
 *
 * @param v Vector.
 */
static inline double S9VectorAngle(const FloatPoint v) {
    return atan2(-v.Y, v.X);
}


/**
 * Normalized a given angle in radians to -PI..PI range.
 */
static double S9NormalizeAngle(const double radians)
{
    static constexpr double PI2 = M_PI + M_PI;
    static constexpr double PI2_R = 1.0 / PI2;

    const double r = radians * PI2_R;

    return PI2 * (r - round(r));
}


/**
 * Converts absolute knot angle to alpha angle relative to segment.
 */
static double S9AbsoluteToAlpha(const double knotAngle,
    const double segmentAngle)
{
    ASSERT(segmentAngle >= -M_PI);
    ASSERT(segmentAngle <= M_PI);

    return knotAngle - segmentAngle;
}


/**
 * Converts alpha knot angle relative to segment to absolute angle.
 */
static double S9AlphaToAbsolute(const double knotAngle,
    const double segmentAngle)
{
    ASSERT(segmentAngle >= -M_PI);
    ASSERT(segmentAngle <= M_PI);

    return knotAngle + segmentAngle;
}


/**
 * Converts absolute knot angle to beta angle relative to segment.
 */
static double S9AbsoluteToBeta(const double knotAngle,
    const double segmentAngle)
{
    ASSERT(segmentAngle >= -M_PI);
    ASSERT(segmentAngle <= M_PI);

    return segmentAngle - knotAngle;
}


/**
 * Converts beta knot angle relative to segment to absolute angle.
 */
static double S9BetaToAbsolute(const double knotAngle,
    const double segmentAngle)
{
    ASSERT(segmentAngle >= -M_PI);
    ASSERT(segmentAngle <= M_PI);

    return -knotAngle + segmentAngle;
}


/**
 * Rotates a given point around coordinate system origin.
 */
static FloatPoint S9RotatePointAroundOrigin(const FloatPoint v,
    const double radians)
{
    const double s = Sin(radians);
    const double c = Cos(radians);

    const double x = v.X * c + v.Y * s;
    const double y = -v.X * s + v.Y * c;

    return FloatPoint(x, y);
}


static inline double S9InfiniteCurvaturePole(const double beta) noexcept
{
    static constexpr const double K = 0.3;

    const double beta2 = beta * 2.0;

    const double v = sinh(
        (1.0 + K) / (2.0 * K) * atan(K * sin(beta2) / (1.0 + K * cos(beta2))));

    return v / 2.0;
}


/**
 * Returns the rotation range of a knot. Absolute angles between curvature
 * pole of the segment before it and the pole of the segment after it, see
 * S9InfiniteCurvaturePole.
 */
static AngleRange S9CalculateRotationRangeForSegmentPair(
    const double prevKnotAngle, const double nextKnotAngle,
    const double prevSegmentAngle, const double nextSegmentAngle)
{
    // Rotation limit axis for beta angle of the previous segment.
    const double betaLimit = S9BetaToAbsolute(S9InfiniteCurvaturePole(
        S9AbsoluteToAlpha(prevKnotAngle, prevSegmentAngle)), prevSegmentAngle);

    // Rotation limit axis for alpha angle of the next segment.
    const double alphaLimit = S9AlphaToAbsolute(S9InfiniteCurvaturePole(
        S9AbsoluteToBeta(nextKnotAngle, nextSegmentAngle)), nextSegmentAngle);

    const double axis = (betaLimit + alphaLimit) / 2.0;
    const double range = (M_PI - Abs(S9NormalizeAngle(betaLimit - alphaLimit))) / 2.0;

    return AngleRange(axis - range, axis + range);
}


/**
 * Keeps the last solved S9 curve of one segment. The next solve of the
 * segment continues it, which is much faster, see S9TrySolveCurveNear.
 */
struct S9SegmentCurve final {

    // Last solved shape. It is valid only when IsKnown is true.
    S9Curve Shape;

    // True after the first successful solve of the segment.
    bool IsKnown = false;
};


/**
 * Returns the last solved shape of a segment, or null when there is none or
 * curve is nullptr.
 */
static inline const S9Curve *S9KnownShape(
    const S9SegmentCurve *const curve) noexcept
{
    if (curve == nullptr or !curve->IsKnown) {
        return nullptr;
    }

    return &curve->Shape;
}


/**
 * Keeps a solved shape as the last curve of a segment. A nullptr curve keeps
 * nothing.
 */
static inline void S9KeepShape(const S9Curve &shape,
    S9SegmentCurve *const curve) noexcept
{
    if (curve == nullptr) {
        return;
    }

    curve->Shape = shape;
    curve->IsKnown = true;
}


/**
 * Returns signed tangent length of one endpoint.
 *
 * @param alpha Finite angle of this endpoint in radians. At the right end,
 * use the clockwise angle from the chord.
 *
 * @param beta Finite angle of the opposite endpoint in radians, in the same
 * convention.
 */
static inline double S9TangentLength(const double alpha,
    const double beta) noexcept
{
    ASSERT(isfinite(alpha));
    ASSERT(isfinite(beta));

    static constexpr const double T = 0.9;

    const double axis = S9InfiniteCurvaturePole(beta);
    const double phase = alpha - axis;

    return tanh(T * cos(phase)) / tanh(T);
}


/**
 * Returns end weight of one signed tangent length.
 *
 * @param length Finite signed tangent length from minus one through one.
 */
static inline double S9EndWeight(const double length) noexcept
{
    ASSERT(isfinite(length));
    ASSERT(fabs(length) <= 1.0);

    return Max(pow(Abs(length), 2.8), 1.7e-10);
}


static inline double S9RoundingLift(const double length) noexcept
{
    ASSERT(isfinite(length));
    ASSERT(fabs(length) <= 1.0);

    return 2.4 * pow(cos((M_PI / 2.0) * length), 6.0);
}


static inline S9CurveParameters S9Lengths(const double alpha,
    const double beta) noexcept
{
    ASSERT(isfinite(alpha));
    ASSERT(isfinite(beta));

    const double startTangent = S9TangentLength(alpha, beta);
    const double endTangent = S9TangentLength(beta, alpha);
    const double startWeight = S9EndWeight(startTangent);
    const double endWeight = S9EndWeight(endTangent);

    // One for two flat ends, near two when either end is a corner.
    const double middle = 2.0 - startWeight * endWeight;

    ASSERT(middle >= 1.0);
    ASSERT(middle <= 2.0);

    S9CurveParameters parameters;

    parameters.Alpha = alpha;
    parameters.Beta = beta;
    parameters.StartLength = startWeight / middle;
    parameters.EndLength = endWeight / middle;
    parameters.StartRoundingLift = S9RoundingLift(parameters.StartLength);
    parameters.EndRoundingLift = S9RoundingLift(parameters.EndLength);

    ASSERT(parameters.StartLength > 0.0);
    ASSERT(parameters.StartLength <= 1.0);
    ASSERT(parameters.EndLength > 0.0);
    ASSERT(parameters.EndLength <= 1.0);

    if (startTangent < 0.0) {
        parameters.StartLength = -parameters.StartLength;
    }

    if (endTangent < 0.0) {
        parameters.EndLength = -parameters.EndLength;
    }

    return parameters;
}


static inline bool S9TrySolveCurve(const double alpha, const double beta,
    S9Curve &curve) noexcept
{
    ASSERT(isfinite(alpha));
    ASSERT(isfinite(beta));

    const S9CurveParameters parameters = S9Lengths(alpha, beta);

    return curve.TrySolveWithLengths(parameters);
}


static inline bool S9TrySolveCurveNear(const double alpha, const double beta,
    const S9Curve &nearby, S9Curve &curve) noexcept
{
    ASSERT(isfinite(alpha));
    ASSERT(isfinite(beta));
    ASSERT(nearby.Length > 0.0);

    const S9CurveParameters parameters = S9Lengths(alpha, beta);

    return curve.TrySolveWithLengthsNear(parameters, nearby);
}


static inline bool S9TryCalculateEndCurvatures(const double alpha,
    const double beta, const S9Curve *const nearby, S9Curve &curve,
    double &startCurvature, double &endCurvature) noexcept
{
    ASSERT(isfinite(alpha));
    ASSERT(isfinite(beta));

    S9Curve solved;
    bool isSolved = false;

    if (nearby != nullptr) {
        isSolved = S9TrySolveCurveNear(alpha, beta, *nearby, solved);
    } else {
        isSolved = S9TrySolveCurve(alpha, beta, solved);
    }

    if (!isSolved) {
        return false;
    }

    startCurvature = solved.Curvature(0.0);
    endCurvature = solved.Curvature(1.0);
    curve = solved;

    return true;
}


/**
 * Returns curvature delta at the center knot of two neighbouring segments.
 * Curvature at the start of the next segment minus curvature at the end of
 * the previous segment, each for its own chord length. The delta becomes very
 * large when a curvature length at the center knot goes to zero.
 *
 * The angles x, y and z are the absolute angles of the previous, center
 * and next knots, a0 and a1 those of the two segments, and h0 and h1 the
 * two chord lengths, greater than DBL_EPSILON. The optional last curves
 * of the two segments are continued when they are known, and both are
 * replaced when both solves succeed.
 */
static double S9CalculateCurvatureDelta(const double x,
    const double y, const double z, const double a0,
    const double a1, const double h0, const double h1,
    S9SegmentCurve *const prevCurve, S9SegmentCurve *const nextCurve)
{
    ASSERT(h0 > DBL_EPSILON);
    ASSERT(h1 > DBL_EPSILON);

    const double alpha0 = S9AbsoluteToAlpha(x, a0);
    const double beta0 = S9AbsoluteToBeta(y, a0);
    const double alpha1 = S9AbsoluteToAlpha(y, a1);
    const double beta1 = S9AbsoluteToBeta(z, a1);

    double startCurvature0 = 0.0;
    double unitCurvature0 = 0.0;
    double unitCurvature1 = 0.0;
    double endCurvature1 = 0.0;

    const S9Curve *const nearbyShape0 = S9KnownShape(prevCurve);
    const S9Curve *const nearbyShape1 = S9KnownShape(nextCurve);

    S9Curve shape0;
    S9Curve shape1;

    const bool isSolved0 = S9TryCalculateEndCurvatures(alpha0, beta0,
        nearbyShape0, shape0, startCurvature0, unitCurvature0);

    if (!isSolved0) {
        return NAN;
    }

    const bool isSolved1 = S9TryCalculateEndCurvatures(alpha1, beta1,
        nearbyShape1, shape1, unitCurvature1, endCurvature1);

    if (!isSolved1) {
        return NAN;
    }

    S9KeepShape(shape0, prevCurve);
    S9KeepShape(shape1, nextCurve);

    // Curvature is inversely proportional to the chord length.
    const double curvature0 = unitCurvature0 / h0;
    const double curvature1 = unitCurvature1 / h1;

    return curvature1 - curvature0;
}


/**
 * Returns the middle of rotation range of a knot, which is the start angle of
 * the solver for it.
 */
static double S9ClampKnotAngle(const double prevKnotAngle,
    const double currentKnotAngle, const double nextKnotAngle,
    const double prevSegmentAngle, const double nextSegmentAngle)
{
    const AngleRange absoluteRotationRange =
        S9CalculateRotationRangeForSegmentPair(prevKnotAngle,
            nextKnotAngle, prevSegmentAngle, nextSegmentAngle);

    const double min = absoluteRotationRange.MinAngle;
    const double max = absoluteRotationRange.MaxAngle;
    const double sum = min + max;

    return sum * 0.5;
}


/*
 * Distance, in radians, from an end of rotation range at which the solver
 * reads the sign of curvature delta. At the end itself curvature length is
 * zero, and rounding can make it negative, which reverses that end and can
 * change sign. The step is at most a quarter of the range.
 */
static constexpr double S9SignProbeStep = 1e-4;


/**
 * Finds the knot angle that makes the curvature delta zero with 32
 * bisection steps over the complete rotation range. It is the fallback of
 * S9SolveKnotAngle. The sign near the start of the range is read with
 * solves from nothing. The steps continue the last curves of the two
 * segments when they are given. The angles are absolute and in radians,
 * the lengths greater than DBL_EPSILON, and the last curves can be null.
 */
static double S9BisectionAngle(const double prevKnotAngle,
    const double nextKnotAngle, const double prevSegmentAngle,
    const double nextSegmentAngle, const double prevSegmentLength,
    const double nextSegmentLength, S9SegmentCurve *const prevCurve,
    S9SegmentCurve *const nextCurve)
{
    const AngleRange absoluteRotationRange =
        S9CalculateRotationRangeForSegmentPair(prevKnotAngle,
            nextKnotAngle, prevSegmentAngle, nextSegmentAngle);

    double min = absoluteRotationRange.MinAngle;
    double max = absoluteRotationRange.MaxAngle;

    const double sum = min + max;

    double angle = sum / 2.0;

    // Sign is read a small step inside the range, see S9SignProbeStep.
    const double range = max - min;
    const double probeStep = Min(S9SignProbeStep, 0.25 * range);
    const double probe = min + probeStep;
    const double v = S9CalculateCurvatureDelta(prevKnotAngle, probe,
        nextKnotAngle, prevSegmentAngle, nextSegmentAngle, prevSegmentLength,
        nextSegmentLength, nullptr, nullptr);

    const I8 sa = static_cast<I8>(copysign(1, v));

    for (U8 v = 0; v < 32; v++) {
        const double y = S9CalculateCurvatureDelta(prevKnotAngle, angle,
            nextKnotAngle, prevSegmentAngle, nextSegmentAngle,
            prevSegmentLength, nextSegmentLength, prevCurve, nextCurve);

        const I8 s = static_cast<I8>(copysign(1, y));

        if (sa == s) {
            min = angle;
        } else {
            max = angle;
        }

        const double sum = min + max;

        angle = sum / 2.0;
    }

    return angle;
}


/*
 * Limits of S9SolveKnotAngle. Tolerance is the bracket width at which the
 * root is found. It is at the noise floor of the knot angles - with smaller
 * width angles still move by up to 3e-13 radians from sweep to sweep because
 * of rounding in the S9 solves. The first step starts the search for a
 * bracket when the knot did not move in the previous sweep, and the smallest
 * step is its least value when it did.
 */
static constexpr double S9KnotAngleTolerance = 1e-13;
static constexpr U32 S9GreatestRootIterationCount = 64U;
static constexpr double S9FirstBracketStep = 1e-3;
static constexpr double S9SmallestBracketStep = 1e-9;
static constexpr double S9BracketStepGrowth = 4.0;


/**
 * Returns true if a finite curvature delta is on the low side of the root of
 * a knot, for sign of the deltas near the start of the rotation range.
 */
static inline bool S9DeltaIsBelowRoot(const double value,
    const bool lowIsNegative) noexcept
{
    const bool isNegative = signbit(value);

    return isNegative == lowIsNegative;
}


/**
 * Finds knot angle that makes curvature delta zero, starting at the current
 * angle of the knot.
 *
 * This function needs the sign of delta near the start of the rotation range.
 * The first call for a knot reads it with solves from nothing and keeps it in
 * lowSign, later calls use the stored sign. It then searches for a bracket of
 * the root with steps that start at twice the previous update of the knot and
 * grow 4 times each step, and refines the bracket with the Illinois form of
 * the false position method, with the middle of the bracket every 4th step.
 * All deltas after the sign continue the last curves of the two segments.
 */
static double S9SolveKnotAngle(const double prevKnotAngle,
    const double currentKnotAngle, const double nextKnotAngle,
    const double prevSegmentAngle, const double nextSegmentAngle,
    const double prevSegmentLength, const double nextSegmentLength,
    const double previousUpdate, S9SegmentCurve &prevCurve,
    S9SegmentCurve &nextCurve, I8 &lowSign) noexcept
{
    ASSERT(isfinite(currentKnotAngle));
    ASSERT(isfinite(previousUpdate));
    ASSERT(prevSegmentLength > DBL_EPSILON);
    ASSERT(nextSegmentLength > DBL_EPSILON);
    ASSERT(lowSign >= -1 and lowSign <= 1);

    const AngleRange absoluteRotationRange =
        S9CalculateRotationRangeForSegmentPair(prevKnotAngle,
            nextKnotAngle, prevSegmentAngle, nextSegmentAngle);

    const double min = absoluteRotationRange.MinAngle;
    const double max = absoluteRotationRange.MaxAngle;
    const double range = max - min;
    const double probeStep = Min(S9SignProbeStep, 0.25 * range);
    const double low = min + probeStep;
    const double high = max - probeStep;

    if (!(high > low)) {
        lowSign = 0;

        return S9BisectionAngle(prevKnotAngle, nextKnotAngle,
            prevSegmentAngle, nextSegmentAngle, prevSegmentLength,
            nextSegmentLength, &prevCurve, &nextCurve);
    }

    // Solves from nothing - the last curves are no good start at a pole.
    if (lowSign == 0) {
        const double lowValue = S9CalculateCurvatureDelta(prevKnotAngle, low,
            nextKnotAngle, prevSegmentAngle, nextSegmentAngle,
            prevSegmentLength, nextSegmentLength, nullptr, nullptr);

        if (!isfinite(lowValue)) {
            return S9BisectionAngle(prevKnotAngle, nextKnotAngle,
                prevSegmentAngle, nextSegmentAngle, prevSegmentLength,
                nextSegmentLength, &prevCurve, &nextCurve);
        }

        lowSign = 1;

        if (signbit(lowValue)) {
            lowSign = -1;
        }
    }

    const bool lowIsNegative = lowSign < 0;
    const double guess = Clamp(currentKnotAngle, low, high);
    const double guessValue = S9CalculateCurvatureDelta(prevKnotAngle, guess,
        nextKnotAngle, prevSegmentAngle, nextSegmentAngle, prevSegmentLength,
        nextSegmentLength, &prevCurve, &nextCurve);

    if (!isfinite(guessValue)) {
        lowSign = 0;

        return S9BisectionAngle(prevKnotAngle, nextKnotAngle,
            prevSegmentAngle, nextSegmentAngle, prevSegmentLength,
            nextSegmentLength, &prevCurve, &nextCurve);
    }

    // Search for a bracket - lower is below the root and upper above it.
    double firstStep = S9FirstBracketStep;

    if (previousUpdate != 0.0) {
        const double doubleUpdate = 2.0 * Abs(previousUpdate);

        firstStep = Max(doubleUpdate, S9SmallestBracketStep);
    }

    double step = Min(firstStep, 0.25 * range);
    double lower = guess;
    double lowerValue = guessValue;
    double upper = guess;
    double upperValue = guessValue;
    const bool guessIsBelow = S9DeltaIsBelowRoot(guessValue, lowIsNegative);
    bool hasBracket = false;

    for (U32 index = 0; index < S9GreatestRootIterationCount; index++) {
        double candidate = Max(upper - step, low);

        if (guessIsBelow) {
            candidate = Min(lower + step, high);
        }

        const double candidateValue = S9CalculateCurvatureDelta(prevKnotAngle,
            candidate, nextKnotAngle, prevSegmentAngle, nextSegmentAngle,
            prevSegmentLength, nextSegmentLength, &prevCurve, &nextCurve);

        if (!isfinite(candidateValue)) {
            break;
        }

        const bool isBelow = S9DeltaIsBelowRoot(candidateValue,
            lowIsNegative);

        if (isBelow) {
            lower = candidate;
            lowerValue = candidateValue;
        } else {
            upper = candidate;
            upperValue = candidateValue;
        }

        if (isBelow != guessIsBelow) {
            hasBracket = true;
            break;
        }

        const bool isAtEnd = candidate <= low or candidate >= high;

        if (isAtEnd) {
            break;
        }

        step *= S9BracketStepGrowth;
    }

    if (!hasBracket) {
        lowSign = 0;

        return S9BisectionAngle(prevKnotAngle, nextKnotAngle,
            prevSegmentAngle, nextSegmentAngle, prevSegmentLength,
            nextSegmentLength, &prevCurve, &nextCurve);
    }

    // Illinois steps.
    //
    // Value at an end that holds for two steps is halved.
    I8 lastSide = 0;

    for (U32 index = 0; index < S9GreatestRootIterationCount; index++) {
        const double width = upper - lower;

        if (width <= S9KnotAngleTolerance) {
            break;
        }

        const double middle = lower + 0.5 * width;
        double candidate = middle;
        const bool takesMiddle = (index % 4U) == 3U;
        const double valueSpan = upperValue - lowerValue;

        if (!takesMiddle and valueSpan != 0.0) {
            const double share = -lowerValue / valueSpan;
            const double linear = lower + share * width;

            if (linear > lower and linear < upper) {
                candidate = linear;
            }
        }

        const double candidateValue = S9CalculateCurvatureDelta(prevKnotAngle,
            candidate, nextKnotAngle, prevSegmentAngle, nextSegmentAngle,
            prevSegmentLength, nextSegmentLength, &prevCurve, &nextCurve);

        if (!isfinite(candidateValue)) {
            lowSign = 0;

            return S9BisectionAngle(prevKnotAngle, nextKnotAngle,
                prevSegmentAngle, nextSegmentAngle, prevSegmentLength,
                nextSegmentLength, &prevCurve, &nextCurve);
        }

        if (candidateValue == 0.0) {
            return candidate;
        }

        const bool isBelow = S9DeltaIsBelowRoot(candidateValue,
            lowIsNegative);

        if (isBelow) {
            lower = candidate;
            lowerValue = candidateValue;

            if (lastSide < 0) {
                upperValue *= 0.5;
            }

            lastSide = -1;
        } else {
            upper = candidate;
            upperValue = candidateValue;

            if (lastSide > 0) {
                lowerValue *= 0.5;
            }

            lastSide = 1;
        }
    }

    return lower + 0.5 * (upper - lower);
}


static double S9AcyclicStartAngle(const double firstSegmentAngle,
    const double secondKnotAngle)
{
    const double beta = S9AbsoluteToBeta(secondKnotAngle, firstSegmentAngle);
    const double alpha = 0.5 * Sin(2.0 * beta);

    return S9AlphaToAbsolute(alpha, firstSegmentAngle);
}


static double S9AcyclicEndAngle(const double finalSegmentAngle,
    const double firstKnotAngle)
{
    const double alpha = S9AbsoluteToAlpha(firstKnotAngle, finalSegmentAngle);
    const double beta = 0.5 * Sin(2.0 * alpha);

    return S9BetaToAbsolute(beta, finalSegmentAngle);
}


static constexpr U32 S9MaximumCenteringIterationCount = 3;

static constexpr U32 S9MaximumITPIterationCount = 20;


/*
 * Sweeps stop when no knot angle moves more than this many radians. Near the
 * solution each sweep makes the largest update only three to four times
 * smaller.
 */
static constexpr double S9MaximumITPError = 1e-11;


/*
 * Fraction of a knot angle update that the solver keeps when the update
 * goes in the opposite direction to the previous update of the same knot.
 */
static constexpr double S9OscillationDamping = 0.7;


/*
 * Distance between the S9 curve of a segment and its cubic curves, as a
 * fraction of the chord length. The conversion checks three points per
 * cubic curve, thus the distance can be a little larger.
 */
static constexpr double S9MaximumCubicError = 1e-5;


/**
 * Returns the next angle of a knot from the angle that S9SolveKnotAngle
 * found. When the update goes in the opposite direction to the previous
 * update of the same knot, only S9OscillationDamping of it is stored, so
 * that a knot does not move back and forth between two angles. An update
 * of half a turn or more is kept complete - a knot angle and the same
 * angle plus pi give the same curves, because S9Lengths changes the sign
 * of the length and the S9 curve then turns the heading back by PI.
 *
 * @param previousUpdate Update of the knot in the previous sweep, zero
 * before the first one. It receives the update that is applied.
 */
static double S9CalculateDampedKnotAngle(const double currentAngle,
    const double solvedAngle, double &previousUpdate) noexcept
{
    const double update = solvedAngle - currentAngle;
    const bool turnsBack = update * previousUpdate < 0.0;
    const bool isSmall = Abs(update) < 0.5 * M_PI;
    double appliedUpdate = update;

    if (turnsBack and isSmall) {
        appliedUpdate = update * S9OscillationDamping;
    }

    previousUpdate = appliedUpdate;

    return currentAngle + appliedUpdate;
}


static inline bool S9TryCreateCurve(const double alpha, const double beta,
    const double maximumError, S9Curve &curve,
    CurveArray<CubicCurve> &curves) noexcept
{
    ASSERT(isfinite(alpha));
    ASSERT(isfinite(beta));
    ASSERT(maximumError > 0.0);

    S9Curve solved;
    const bool isSolved = S9TrySolveCurve(alpha, beta, solved);

    if (!isSolved) {
        const FloatPoint startPoint;

        curves.Erase();
        curves.Restart(startPoint);

        return false;
    }

    const bool isApproximated = solved.TryApproximateWithCubics(maximumError,
        curves);

    if (!isApproximated) {
        return false;
    }

    curve = solved;

    return true;
}


static inline bool S9TryCreateCurveNear(const double alpha, const double beta,
    const double maximumError, const S9Curve &nearby, S9Curve &curve,
    CurveArray<CubicCurve> &curves) noexcept
{
    ASSERT(isfinite(alpha));
    ASSERT(isfinite(beta));
    ASSERT(maximumError > 0.0);
    ASSERT(nearby.Length > 0.0);

    S9Curve solved;
    const bool isSolved = S9TrySolveCurveNear(alpha, beta, nearby, solved);

    if (!isSolved) {
        const FloatPoint startPoint;

        curves.Erase();
        curves.Restart(startPoint);

        return false;
    }

    const bool isApproximated = solved.TryApproximateWithCubics(maximumError,
        curves);

    if (!isApproximated) {
        return false;
    }

    curve = solved;

    return true;
}


static bool S9AppendCubics(const double alpha, const double beta,
    const FloatPoint startPoint, const FloatPoint endPoint,
    const double segmentAngle, const double segmentLength,
    const S9SegmentCurve *const cache,
    CurveArray<CubicCurve> &pieces, CurveArray<CubicCurve> &curves) noexcept
{
    ASSERT(segmentLength > DBL_EPSILON);

    S9Curve shape;
    const S9Curve *const nearby = S9KnownShape(cache);
    bool isCreated = false;

    if (nearby != nullptr) {
        isCreated = S9TryCreateCurveNear(alpha, beta,
            S9MaximumCubicError, *nearby, shape, pieces);
    } else {
        isCreated = S9TryCreateCurve(alpha, beta,
            S9MaximumCubicError, shape, pieces);
    }

    if (!isCreated) {
        return false;
    }

    const U32 curveCount = pieces.GetCurveCount();
    ASSERT(curveCount > 0U);

    const U32 pointCount = pieces.GetPointCount();
    const U32 lastIndex = pointCount - 1U;

    for (U32 index = 1U; index < lastIndex; index++) {
        const FloatPoint point = pieces.GetPointAtIndex(index);
        const FloatPoint mirrored(point.X * segmentLength,
            -point.Y * segmentLength);
        const FloatPoint rotated = S9RotatePointAroundOrigin(mirrored,
            segmentAngle);

        pieces.SetPointAtIndex(index, startPoint + rotated);
    }

    // The spline must go through its knots exactly.
    pieces.SetPointAtIndex(lastIndex, endPoint);

    for (U32 index = 0U; index < curveCount; ++index) {
        const CubicCurve cubic = pieces.GetCurveAtIndex(index);

        curves.CurveTo(cubic.P1, cubic.P2, cubic.P3);
    }

    return true;
}


/*
 * All knots and segments of one spline, for S9SolveKnotSystem. The
 * arrays belong to the caller.
 */
struct S9KnotSystem final {

    // Number of knots.
    U32 Count = 0U;

    // Number of segments. One less than the knot count for an open spline,
    // the knot count for a closed spline.
    U32 SegmentCount = 0U;

    // True for a closed spline. Then the last segment ends at knot zero.
    bool IsCyclic = false;

    // True when an end control point gives the angle of the first knot of
    // an open spline. Then StartAngle holds it.
    bool HasStartAngle = false;

    // True when an end control point gives the angle of the last knot of
    // an open spline. Then EndAngle holds it.
    bool HasEndAngle = false;

    // Absolute angle of the first knot in radians. It is valid only when
    // HasStartAngle is true.
    double StartAngle = 0.0;

    // Absolute angle of the last knot in radians. It is valid only when
    // HasEndAngle is true.
    double EndAngle = 0.0;

    // Absolute angle of each segment in radians.
    const double *SegmentAngles = nullptr;

    // Chord length of each segment, greater than DBL_EPSILON.
    const double *SegmentLengths = nullptr;

    // Absolute angle of each knot in radians. The solver updates them.
    double *KnotAngles = nullptr;

    // Last solved curve of each segment. The solver updates them.
    S9SegmentCurve *SegmentCurves = nullptr;

    // True for each knot that is pinned to the middle of its rotation
    // range, see S9ClassifyKnots. The solver updates them.
    bool *Pinned = nullptr;
};


/*
 * Limits of S9SolveKnotSystem. Sweeps continue while the greatest update is
 * at least the start error of the coupled Newton steps. Step is accepted when
 * the greatest curvature delta goes down by at least the decrease fraction,
 * and no knot moves more than the range fraction of its distance to an end of
 * its rotation range. When that cap scales a step below the stall factor for
 * the stall count of steps in a row, root is very near a pole and the
 * iteration gives up.
 */
static constexpr U32 S9KnotSweepCount = 12U;
static constexpr double S9CoupledStartError = 3e-4;
static constexpr U32 S9CoupledGreatestIterationCount = 40U;
static constexpr U32 S9CoupledGreatestHalvingCount = 5U;
static constexpr double S9CoupledDifferenceStep = 1e-6;
static constexpr double S9CoupledDecrease = 1e-4;
static constexpr double S9CoupledRangeFraction = 0.25;
static constexpr double S9CoupledStallFactor = 0.2;
static constexpr U32 S9CoupledGreatestStallCount = 3U;


/*
 * Width of the rotation range, in radians, under which a knot is pinned to
 * the middle of its range, see S9ClassifyKnots. Both ends of such a knot
 * are within 0.3 degrees of their curvature poles.
 */
static constexpr double S9PinnedRangeWidth = 1e-2;


/*
 * Width of the rotation range, in radians, at which a pinned knot becomes
 * free again. It is wider than S9PinnedRangeWidth so that a knot whose
 * range moves about the pinning width does not change its mark at every
 * iteration.
 */
static constexpr double S9UnpinRangeWidth = 2e-2;


/*
 * Passes of S9CentrePinnedKnots. Each pass settles a chain of pinned
 * knots by about a factor of five.
 */
static constexpr U32 S9CentringPassCount = 3U;


static U32 S9UnknownCount(const S9KnotSystem &system) noexcept
{
    ASSERT(system.Count >= 2U);

    if (system.IsCyclic) {
        return system.Count;
    }

    return system.Count - 2U;
}


/**
 * Returns the knot of an unknown. The index must be less than S9UnknownCount.
 */
static U32 S9KnotOfUnknown(const S9KnotSystem &system,
    const U32 unknown) noexcept
{
    if (system.IsCyclic) {
        return unknown;
    }

    return unknown + 1U;
}


/**
 * Returns the segment that ends at a knot. For an open spline the knot must
 * not be zero.
 */
static U32 S9SegmentBeforeKnot(const S9KnotSystem &system,
    const U32 knot) noexcept
{
    if (knot == 0U) {
        ASSERT(system.IsCyclic);

        return system.Count - 1U;
    }

    return knot - 1U;
}


/**
 * Returns the knot at which a segment ends.
 */
static U32 S9EndKnotOfSegment(const S9KnotSystem &system,
    const U32 segment) noexcept
{
    ASSERT(segment < system.SegmentCount);

    const U32 next = segment + 1U;

    if (system.IsCyclic and next == system.Count) {
        return 0U;
    }

    return next;
}


/**
 * Returns the angle of the first knot of an open spline. The angle of its end
 * control point, or the end rule from its neighbour.
 */
static double S9StartKnotAngle(const S9KnotSystem &system,
    const double *const angles) noexcept
{
    ASSERT(!system.IsCyclic);
    ASSERT(angles != nullptr);

    if (system.HasStartAngle) {
        return system.StartAngle;
    }

    return S9AcyclicStartAngle(system.SegmentAngles[0], angles[1]);
}


/**
 * Returns the angle of the last knot of an open spline. The angle of its end
 * control point, or the end rule from its neighbour.
 */
static double S9EndKnotAngle(const S9KnotSystem &system,
    const double *const angles) noexcept
{
    ASSERT(!system.IsCyclic);
    ASSERT(angles != nullptr);

    if (system.HasEndAngle) {
        return system.EndAngle;
    }

    const U32 last = system.Count - 1U;

    return S9AcyclicEndAngle(system.SegmentAngles[last - 1U],
        angles[last - 1U]);
}


/**
 * Sets the angles of the two end knots of an open spline from their end
 * control points or the end rules. It does nothing for a closed spline.
 */
static void S9ApplyEndRules(const S9KnotSystem &system,
    double *const angles) noexcept
{
    ASSERT(angles != nullptr);

    if (system.IsCyclic) {
        return;
    }

    const U32 last = system.Count - 1U;

    angles[0] = S9StartKnotAngle(system, angles);
    angles[last] = S9EndKnotAngle(system, angles);
}


static void S9ApplyEndControls(const S9SplineEndControls &controls,
    const FloatPoint *const knots, S9KnotSystem &system) noexcept
{
    ASSERT(knots != nullptr);
    ASSERT(!system.IsCyclic);
    ASSERT(system.Count >= 2U);
    ASSERT(system.SegmentAngles != nullptr);

    const U32 last = system.Count - 1U;
    const U32 lastSegment = system.SegmentCount - 1U;

    if (controls.HasStartControl) {
        const FloatPoint handle = controls.StartControl - knots[0];
        const double distance = handle.Length();

        ASSERT(isfinite(distance));
        ASSERT(distance > DBL_EPSILON);

        if (isfinite(distance) and distance > DBL_EPSILON) {
            const double segmentAngle = system.SegmentAngles[0];
            const double handleAngle = S9VectorAngle(handle);
            const double alpha = S9NormalizeAngle(handleAngle - segmentAngle);

            system.HasStartAngle = true;
            system.StartAngle = S9AlphaToAbsolute(alpha, segmentAngle);
        }
    }

    if (controls.HasEndControl) {
        // Spline arrives at the last knot from its control point.
        const FloatPoint handle = knots[last] - controls.EndControl;
        const double distance = handle.Length();

        ASSERT(isfinite(distance));
        ASSERT(distance > DBL_EPSILON);

        if (isfinite(distance) and distance > DBL_EPSILON) {
            const double segmentAngle = system.SegmentAngles[lastSegment];
            const double handleAngle = S9VectorAngle(handle);
            const double beta = S9NormalizeAngle(segmentAngle - handleAngle);

            system.HasEndAngle = true;
            system.EndAngle = S9BetaToAbsolute(beta, segmentAngle);
        }
    }
}


/**
 * Solves one segment for a set of knot angles and returns its end
 * curvatures for its chord length.
 */
static bool S9SolveSegment(const S9KnotSystem &system,
    const double *const angles, const U32 segment,
    S9SegmentCurve &cache, const bool keep, double &startCurvature,
    double &endCurvature) noexcept
{
    ASSERT(angles != nullptr);
    ASSERT(segment < system.SegmentCount);

    const U32 startKnot = segment;
    const U32 endKnot = S9EndKnotOfSegment(system, segment);
    const double segmentAngle = system.SegmentAngles[segment];
    const double length = system.SegmentLengths[segment];

    ASSERT(length > DBL_EPSILON);

    const double alpha = S9AbsoluteToAlpha(angles[startKnot], segmentAngle);
    const double beta = S9AbsoluteToBeta(angles[endKnot], segmentAngle);
    const S9Curve *const nearby = S9KnownShape(&cache);

    S9Curve shape;

    double unitStart = 0.0;
    double unitEnd = 0.0;

    const bool isSolved = S9TryCalculateEndCurvatures(alpha, beta, nearby,
        shape, unitStart, unitEnd);

    if (!isSolved) {
        return false;
    }

    if (keep) {
        S9KeepShape(shape, &cache);
    }

    // Curvature is inversely proportional to the chord length.
    startCurvature = unitStart / length;
    endCurvature = unitEnd / length;

    return true;
}


/**
 * Solves every segment that has no stored curve yet, from nothing, and keeps
 * the shapes.
 */
static void S9SolveUnknownSegments(const S9KnotSystem &system) noexcept
{
    for (U32 segment = 0U; segment < system.SegmentCount; segment++) {
        S9SegmentCurve &cache = system.SegmentCurves[segment];

        if (cache.IsKnown) {
            continue;
        }

        double startCurvature = 0.0;
        double endCurvature = 0.0;

        // A failed solve leaves the segment unknown.
        S9SolveSegment(system, system.KnotAngles, segment, cache, true,
            startCurvature, endCurvature);
    }
}


/**
 * Solves all segments for a set of knot angles and calculates the
 * curvature delta of every unknown knot.
 */
static bool S9EvaluateKnotSystem(const S9KnotSystem &system,
    const double *const angles, const bool keep, double *const startCurvatures,
    double *const endCurvatures, double *const deltas) noexcept
{
    for (U32 segment = 0U; segment < system.SegmentCount; segment++) {
        const bool isSolved = S9SolveSegment(system, angles, segment,
            system.SegmentCurves[segment], keep, startCurvatures[segment],
            endCurvatures[segment]);

        if (!isSolved) {
            return false;
        }
    }

    const U32 unknownCount = S9UnknownCount(system);

    for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
        const U32 knot = S9KnotOfUnknown(system, unknown);
        const U32 before = S9SegmentBeforeKnot(system, knot);

        deltas[unknown] = startCurvatures[knot] - endCurvatures[before];
    }

    return true;
}


/**
 * Calculates rotation range of one unknown knot, a probe step inside its
 * ends, see S9SignProbeStep. Returns true when the range is not empty.
 */
static bool S9CalculateKnotRange(const S9KnotSystem &system,
    const double *const angles, const U32 knot, double &low,
    double &high) noexcept
{
    const U32 before = S9SegmentBeforeKnot(system, knot);
    const U32 after = knot;
    const U32 nextKnot = S9EndKnotOfSegment(system, after);
    const AngleRange range = S9CalculateRotationRangeForSegmentPair(
        angles[before], angles[nextKnot], system.SegmentAngles[before],
        system.SegmentAngles[after]);
    const double width = range.MaxAngle - range.MinAngle;
    const double probeStep = Min(S9SignProbeStep, 0.25 * width);

    low = range.MinAngle + probeStep;
    high = range.MaxAngle - probeStep;

    return high > low;
}



/**
 * Returns the width of the rotation range of a knot in radians, see
 * S9CalculateRotationRangeForSegmentPair.
 *
 * @param angles Non-null absolute knot angles in radians.
 *
 * @param knot Index of an unknown knot, less than the knot count.
 */
static double S9KnotRangeWidth(const S9KnotSystem &system,
    const double *const angles, const U32 knot) noexcept
{
    ASSERT(angles != nullptr);
    ASSERT(knot < system.Count);

    const U32 before = S9SegmentBeforeKnot(system, knot);
    const U32 nextKnot = S9EndKnotOfSegment(system, knot);
    const AngleRange range = S9CalculateRotationRangeForSegmentPair(
        angles[before], angles[nextKnot], system.SegmentAngles[before],
        system.SegmentAngles[knot]);

    return range.MaxAngle - range.MinAngle;
}


/**
 * Marks every unknown knot whose rotation range is narrower than
 * S9PinnedRangeWidth as pinned, and every other unknown knot as free, in
 * system.Pinned. Both ends of a pinned knot are so near their curvature
 * poles that the end weights sit at their floor, where the curvature
 * delta no longer changes with the knot angle and can have no root. Such a
 * knot stays in the middle of its range and follows its neighbours. A
 * pinned knot is freed only when its range is wider than
 * S9UnpinRangeWidth.
 *
 * @param angles Non-null absolute knot angles in radians.
 *
 * @return True when the mark of a knot changed.
 */
static bool S9ClassifyKnots(const S9KnotSystem &system,
    const double *const angles) noexcept
{
    ASSERT(angles != nullptr);
    ASSERT(system.Pinned != nullptr);

    const U32 unknownCount = S9UnknownCount(system);
    bool changed = false;

    for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
        const U32 knot = S9KnotOfUnknown(system, unknown);
        const double width = S9KnotRangeWidth(system, angles, knot);
        bool isPinned = width < S9PinnedRangeWidth;

        // A pinned knot stays pinned until its range widens clearly.
        if (system.Pinned[knot] and width < S9UnpinRangeWidth) {
            isPinned = true;
        }

        if (system.Pinned[knot] != isPinned) {
            changed = true;
            system.Pinned[knot] = isPinned;
        }
    }

    return changed;
}


/**
 * Moves every pinned knot to the middle of its rotation range, in knot
 * order, with the end rules applied before each pass and after the last
 * one. The middle of a pinned knot depends on the angles of its
 * neighbours, which can be pinned too, thus S9CentringPassCount passes
 * settle chains of pinned knots.
 *
 * @param angles Non-null absolute knot angles in radians. Pinned knots are
 * replaced.
 */
static void S9CentrePinnedKnots(const S9KnotSystem &system,
    double *const angles) noexcept
{
    ASSERT(angles != nullptr);
    ASSERT(system.Pinned != nullptr);

    const U32 unknownCount = S9UnknownCount(system);

    for (U32 pass = 0U; pass < S9CentringPassCount; pass++) {
        S9ApplyEndRules(system, angles);

        for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
            const U32 knot = S9KnotOfUnknown(system, unknown);

            if (!system.Pinned[knot]) {
                continue;
            }

            const U32 before = S9SegmentBeforeKnot(system, knot);
            const U32 nextKnot = S9EndKnotOfSegment(system, knot);

            angles[knot] = S9ClampKnotAngle(angles[before], angles[knot],
                angles[nextKnot], system.SegmentAngles[before],
                system.SegmentAngles[knot]);
        }
    }

    S9ApplyEndRules(system, angles);
}


/**
 * Clamps every free knot into its rotation range and applies the end
 * rules. A knot leaves its range when the ranges move with its
 * neighbours.
 *
 * @param angles Non-null absolute knot angles in radians. Free knots
 * outside their ranges are replaced.
 */
static void S9ClampFreeKnots(const S9KnotSystem &system,
    double *const angles) noexcept
{
    ASSERT(angles != nullptr);
    ASSERT(system.Pinned != nullptr);

    const U32 unknownCount = S9UnknownCount(system);

    for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
        const U32 knot = S9KnotOfUnknown(system, unknown);

        if (system.Pinned[knot]) {
            continue;
        }

        double low = 0.0;
        double high = 0.0;
        const bool hasRange = S9CalculateKnotRange(system, angles, knot,
            low, high);

        if (hasRange) {
            angles[knot] = Clamp(angles[knot], low, high);
        }
    }

    S9ApplyEndRules(system, angles);
}


/**
 * Returns the greatest magnitude of the curvature deltas of the free
 * knots. Pinned knots have no root to reach and are left out.
 *
 * @param unknownCount Number of unknown knots, see S9UnknownCount.
 *
 * @param deltas Non-null curvature delta of each unknown knot.
 */
static double S9FreeDeltaNorm(const S9KnotSystem &system,
    const U32 unknownCount, const double *const deltas) noexcept
{
    ASSERT(deltas != nullptr);
    ASSERT(system.Pinned != nullptr);

    double greatest = 0.0;

    for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
        const U32 knot = S9KnotOfUnknown(system, unknown);

        if (system.Pinned[knot]) {
            continue;
        }

        greatest = Max(greatest, Abs(deltas[unknown]));
    }

    return greatest;
}


/**
 * Writes the row of a pinned knot into the Newton system: the linearized
 * constraint that the knot angle equals the middle of its rotation range.
 * The middle is differentiated by finite differences with respect to the
 * unknowns before and after the knot and the knot itself, with the end
 * rules applied, because the end knots of an open spline follow their
 * neighbours.
 *
 * @param unknown Index of the pinned knot among the unknowns, less than
 * unknownCount.
 *
 * @param angles Non-null absolute knot angles in radians.
 *
 * @param trialAngles Non-null scratch space for all knot angles.
 *
 * @param subdiagonal Non-null elements left of the diagonal, one per
 * unknown. The element of the row is replaced.
 *
 * @param diagonal Non-null diagonal elements, one per unknown. The
 * element of the row is replaced.
 *
 * @param superdiagonal Non-null elements right of the diagonal, one per
 * unknown. The element of the row is replaced.
 *
 * @param rightSide Non-null right side of the system, one per unknown.
 * The element of the row is replaced.
 */
static void S9WritePinnedRow(const S9KnotSystem &system,
    const U32 unknownCount, const U32 unknown, const double *const angles,
    double *const trialAngles, double *const subdiagonal,
    double *const diagonal, double *const superdiagonal,
    double *const rightSide) noexcept
{
    ASSERT(unknown < unknownCount);
    ASSERT(angles != nullptr);
    ASSERT(trialAngles != nullptr);

    const U32 count = system.Count;
    const U32 knot = S9KnotOfUnknown(system, unknown);
    const U32 before = S9SegmentBeforeKnot(system, knot);
    const U32 nextKnot = S9EndKnotOfSegment(system, knot);
    const double middle = S9ClampKnotAngle(angles[before], angles[knot],
        angles[nextKnot], system.SegmentAngles[before],
        system.SegmentAngles[knot]);
    const bool hasPrevious = unknown > 0U or system.IsCyclic;
    const bool hasNext = unknown + 1U < unknownCount or system.IsCyclic;
    const U32 previous = (unknown + unknownCount - 1U) % unknownCount;
    const U32 next = (unknown + 1U) % unknownCount;
    const U32 columns[3] = {previous, unknown, next};
    const bool present[3] = {hasPrevious, true, hasNext};
    double slopes[3] = {0.0, 0.0, 0.0};

    for (U32 side = 0U; side < 3U; side++) {
        if (!present[side]) {
            continue;
        }

        const U32 columnKnot = S9KnotOfUnknown(system, columns[side]);

        for (U32 index = 0U; index < count; index++) {
            trialAngles[index] = angles[index];
        }

        trialAngles[columnKnot] += S9CoupledDifferenceStep;

        S9ApplyEndRules(system, trialAngles);

        const double moved = S9ClampKnotAngle(trialAngles[before],
            trialAngles[knot], trialAngles[nextKnot],
            system.SegmentAngles[before], system.SegmentAngles[knot]);

        slopes[side] = (moved - middle) / S9CoupledDifferenceStep;
    }

    diagonal[unknown] = 1.0 - slopes[1];
    rightSide[unknown] = middle - angles[knot];

    if (hasPrevious) {
        subdiagonal[unknown] = -slopes[0];
    }

    if (hasNext) {
        superdiagonal[unknown] = -slopes[2];
    }
}


/**
 * Solves the two segments at a knot for a set of knot angles, without keeping
 * the shapes. Receives start and end curvatures of the segment before the knot
 * and of the segment after it. Returns true when both segments were solved.
 */
static bool S9SolveSegmentsAtKnot(const S9KnotSystem &system,
    const double *const angles, const U32 knot, double &beforeStart,
    double &beforeEnd, double &afterStart, double &afterEnd) noexcept
{
    const U32 before = S9SegmentBeforeKnot(system, knot);
    const U32 after = knot;
    S9SegmentCurve beforeCache = system.SegmentCurves[before];
    S9SegmentCurve afterCache = system.SegmentCurves[after];
    const bool beforeSolved = S9SolveSegment(system, angles, before,
        beforeCache, false, beforeStart, beforeEnd);

    if (!beforeSolved) {
        return false;
    }

    return S9SolveSegment(system, angles, after, afterCache, false,
        afterStart, afterEnd);
}


/**
 * Makes one sweep over the knots in which each free knot takes one Newton
 * step on its curvature delta, with the slope from finite difference.
 * Pinned knots are moved to the middle of their ranges first, see
 * S9ClassifyKnots. Step is
 * limited to the range fraction of rotation range of the knot and the a angle
 * is clamped into that range.
 *
 * Returns the greatest update of a knot angle in radians, or a negative value
 * when a segment could not be solved or a slope was zero.
 */
static double S9SweepKnotNewtonSteps(const S9KnotSystem &system,
    double *const trial) noexcept
{
    ASSERT(trial != nullptr);

    const U32 count = system.Count;
    const U32 unknownCount = S9UnknownCount(system);
    double *const angles = system.KnotAngles;
    double greatestUpdate = 0.0;

    S9ClassifyKnots(system, angles);
    S9CentrePinnedKnots(system, angles);

    for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
        const U32 knot = S9KnotOfUnknown(system, unknown);

        // A pinned knot follows its neighbours.
        if (system.Pinned[knot]) {
            continue;
        }

        const U32 before = S9SegmentBeforeKnot(system, knot);
        const U32 after = knot;
        double low = 0.0;
        double high = 0.0;
        const bool hasRange = S9CalculateKnotRange(system, angles, knot,
            low, high);

        if (!hasRange) {
            return -1.0;
        }

        const double entering = angles[knot];
        const double current = Clamp(entering, low, high);

        angles[knot] = current;

        S9ApplyEndRules(system, angles);

        double beforeStart = 0.0;
        double beforeEnd = 0.0;
        double afterStart = 0.0;
        double afterEnd = 0.0;
        const bool beforeSolved = S9SolveSegment(system, angles, before,
            system.SegmentCurves[before], true, beforeStart, beforeEnd);

        if (!beforeSolved) {
            return -1.0;
        }

        const bool afterSolved = S9SolveSegment(system, angles, after,
            system.SegmentCurves[after], true, afterStart, afterEnd);

        if (!afterSolved) {
            return -1.0;
        }

        const double delta = afterStart - beforeEnd;

        // Finite difference goes towards the inside of the range.
        double difference = S9CoupledDifferenceStep;

        if (current + difference > high) {
            difference = -difference;
        }

        for (U32 index = 0U; index < count; index++) {
            trial[index] = angles[index];
        }

        trial[knot] = current + difference;

        S9ApplyEndRules(system, trial);

        double trialBeforeStart = 0.0;
        double trialBeforeEnd = 0.0;
        double trialAfterStart = 0.0;
        double trialAfterEnd = 0.0;
        const bool trialSolved = S9SolveSegmentsAtKnot(system, trial, knot,
            trialBeforeStart, trialBeforeEnd, trialAfterStart, trialAfterEnd);

        if (!trialSolved) {
            return -1.0;
        }

        const double trialDelta = trialAfterStart - trialBeforeEnd;
        const double slope = (trialDelta - delta) / difference;

        if (slope == 0.0 or !isfinite(slope)) {
            return -1.0;
        }

        const double cap = S9CoupledRangeFraction * (high - low);
        const double update = Clamp(-delta / slope, -cap, cap);
        const double next = Clamp(current + update, low, high);
        const double wanted = Abs(update);
        const double moved = Abs(next - current);

        // The update includes the clamp into the range.
        double magnitude = Abs(next - entering);

        // A knot held at an end of its range by an outward step counts the
        // step it wanted as its update, thus it does not look converged.
        if (moved == 0.0) {
            magnitude = Max(magnitude, wanted);
        }

        greatestUpdate = Max(greatestUpdate, magnitude);
        angles[knot] = next;

        S9ApplyEndRules(system, angles);
    }

    return greatestUpdate;
}


/**
 * Makes Newton steps on all knot angles at once until the angles do not move
 * more than S9MaximumITPError. Curvature delta of a knot depends only on the
 * angles of the knot and its two neighbours, thus the Jacobian is
 * tridiagonal, cyclic for a closed spline. Its columns come from finite
 * differences. Each step is scaled so that no knot moves more than the range
 * fraction of its distance to an end of its rotation range, then halved until
 * the greatest curvature delta goes down.
 *
 * Iteration gives up when a segment cannot be solved, when the step cannot be
 * found, when no halving makes deltas go down, or when the range cap scales
 * the steps down for several iterations in a row. The knot angles then keep
 * the last accepted iterate and the segment curves match them. Convergence is
 * judged on the Newton step before the cap and the clamp, thus a knot that
 * cannot move does not look converged. All segment curves must be known.
 *
 * Knots pinned by S9ClassifyKnots are not unknowns of the curvature
 * equations. They stay in the middle of their ranges through the rows of
 * S9WritePinnedRow, and the norm, the cap and the convergence test take
 * only the free knots.
 *
 * Returns true when iteration converged.
 */
static bool S9SolveKnotSystemCoupled(const S9KnotSystem &system) noexcept
{
    const U32 count = system.Count;
    const U32 segmentCount = system.SegmentCount;
    const U32 unknownCount = S9UnknownCount(system);

    if (unknownCount == 0U) {
        return true;
    }

    double *const startCurvatures = AllocElements<double>(segmentCount);
    double *const endCurvatures = AllocElements<double>(segmentCount);
    double *const trialStartCurvatures = AllocElements<double>(segmentCount);
    double *const trialEndCurvatures = AllocElements<double>(segmentCount);
    double *const deltas = AllocElements<double>(unknownCount);
    double *const trialDeltas = AllocElements<double>(unknownCount);
    double *const subdiagonal = AllocElements<double>(unknownCount);
    double *const diagonal = AllocElements<double>(unknownCount);
    double *const superdiagonal = AllocElements<double>(unknownCount);
    double *const rightSide = AllocElements<double>(unknownCount);
    double *const step = AllocElements<double>(unknownCount);
    double *const lows = AllocElements<double>(unknownCount);
    double *const highs = AllocElements<double>(unknownCount);
    double *const trialAngles = AllocElements<double>(count);
    double *const angles = system.KnotAngles;
    bool isComplete = false;
    U32 stallCount = 0U;

    S9ClassifyKnots(system, angles);
    S9ClampFreeKnots(system, angles);
    S9CentrePinnedKnots(system, angles);

    bool isEvaluated = S9EvaluateKnotSystem(system, angles, true,
        startCurvatures, endCurvatures, deltas);
    U32 iteration = 0U;

    while (isEvaluated and iteration < S9CoupledGreatestIterationCount) {
        iteration++;

        // The ranges move with the knots, thus the marks are refreshed.
        if (iteration > 1U) {
            const bool changed = S9ClassifyKnots(system, angles);

            if (changed) {
                S9ClampFreeKnots(system, angles);
                S9CentrePinnedKnots(system, angles);

                isEvaluated = S9EvaluateKnotSystem(system, angles, true,
                    startCurvatures, endCurvatures, deltas);

                if (!isEvaluated) {
                    break;
                }
            }
        }

        const double deltaNorm = S9FreeDeltaNorm(system, unknownCount,
            deltas);

        for (U32 index = 0U; index < unknownCount; index++) {
            subdiagonal[index] = 0.0;
            diagonal[index] = 0.0;
            superdiagonal[index] = 0.0;
        }

        // Finite difference columns, one unknown at a time.
        bool columnsAreValid = true;

        for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
            const U32 knot = S9KnotOfUnknown(system, unknown);
            const U32 before = S9SegmentBeforeKnot(system, knot);
            const U32 after = knot;
            const bool hasRange = S9CalculateKnotRange(system, angles, knot,
                lows[unknown], highs[unknown]);

            if (!hasRange) {
                columnsAreValid = false;
                break;
            }

            for (U32 index = 0U; index < count; index++) {
                trialAngles[index] = angles[index];
            }

            double difference = S9CoupledDifferenceStep;

            if (angles[knot] + difference > highs[unknown]) {
                difference = -difference;
            }

            trialAngles[knot] = angles[knot] + difference;

            S9ApplyEndRules(system, trialAngles);

            double beforeStart = 0.0;
            double beforeEnd = 0.0;
            double afterStart = 0.0;
            double afterEnd = 0.0;
            const bool isSolved = S9SolveSegmentsAtKnot(system, trialAngles,
                knot, beforeStart, beforeEnd, afterStart, afterEnd);

            if (!isSolved) {
                columnsAreValid = false;
                break;
            }

            const double scale = 1.0 / difference;
            const double beforeStartSlope =
                (beforeStart - startCurvatures[before]) * scale;
            const double beforeEndSlope =
                (beforeEnd - endCurvatures[before]) * scale;
            const double afterStartSlope =
                (afterStart - startCurvatures[after]) * scale;
            const double afterEndSlope = (afterEnd - endCurvatures[after]) * scale;

            // Delta of this knot: start of the segment after it minus end
            // of the segment before it.
            diagonal[unknown] = afterStartSlope - beforeEndSlope;

            // Delta of the previous knot: the start of the segment before
            // this knot changes.
            if (unknown > 0U) {
                superdiagonal[unknown - 1U] = beforeStartSlope;
            } else if (system.IsCyclic) {
                superdiagonal[unknownCount - 1U] = beforeStartSlope;
            }

            // Delta of the next knot: the end of the segment after this
            // knot changes.
            if (unknown + 1U < unknownCount) {
                subdiagonal[unknown + 1U] = -afterEndSlope;
            } else if (system.IsCyclic) {
                subdiagonal[0] = -afterEndSlope;
            }
        }

        if (!columnsAreValid) {
            break;
        }

        for (U32 index = 0U; index < unknownCount; index++) {
            rightSide[index] = -deltas[index];
        }

        // A pinned knot has no delta to zero. Its row keeps it in the
        // middle of its range as the neighbours move.
        for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
            const U32 knot = S9KnotOfUnknown(system, unknown);

            if (system.Pinned[knot]) {
                S9WritePinnedRow(system, unknownCount, unknown, angles,
                    trialAngles, subdiagonal, diagonal, superdiagonal,
                    rightSide);
            }
        }

        if (unknownCount == 1U) {
            if (diagonal[0] == 0.0) {
                break;
            }

            step[0] = rightSide[0] / diagonal[0];
        } else if (system.IsCyclic) {
            SolveTridiagonalCyclic(unknownCount, subdiagonal, diagonal,
                superdiagonal, rightSide, step);
        } else {
            SolveTridiagonal(unknownCount, subdiagonal, diagonal, superdiagonal,
                rightSide, step);
        }

        // Cap the step of the free knots. The whole step is scaled, thus
        // its direction is stored.
        double scaleFactor = 1.0;
        double greatestStep = 0.0;
        bool stepIsFinite = true;

        for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
            const U32 knot = S9KnotOfUnknown(system, unknown);
            const double value = step[unknown];

            if (!isfinite(value)) {
                stepIsFinite = false;
                break;
            }

            // A pinned knot follows, thus its step is not capped.
            if (system.Pinned[knot]) {
                continue;
            }

            const double stepMagnitude = Abs(value);

            greatestStep = Max(greatestStep, stepMagnitude);

            double room = highs[unknown] - angles[knot];

            if (value < 0.0) {
                room = angles[knot] - lows[unknown];
            }

            const double positiveRoom = Max(room, 0.0);
            const double allowed = S9CoupledRangeFraction * positiveRoom +
                S9MaximumITPError;
            const double magnitude = Abs(value);

            if (magnitude > allowed) {
                scaleFactor = Min(scaleFactor, allowed / magnitude);
            }
        }

        if (!stepIsFinite) {
            break;
        }

        if (scaleFactor < S9CoupledStallFactor) {
            stallCount++;
        } else {
            stallCount = 0U;
        }

        if (stallCount >= S9CoupledGreatestStallCount) {
            break;
        }

        // Line search on the greatest delta.
        double factor = scaleFactor;
        bool isAccepted = false;
        double greatestUpdate = 0.0;

        for (U32 halving = 0U; halving <= S9CoupledGreatestHalvingCount;
            halving++) {

            for (U32 index = 0U; index < count; index++) {
                trialAngles[index] = angles[index];
            }

            for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
                const U32 knot = S9KnotOfUnknown(system, unknown);
                const double moved = angles[knot] + factor * step[unknown];

                trialAngles[knot] = moved;

                if (!system.Pinned[knot]) {
                    trialAngles[knot] = Clamp(moved, lows[unknown],
                        highs[unknown]);
                }
            }

            // Pinned knots take their exact middles, not the linear ones.
            S9CentrePinnedKnots(system, trialAngles);

            greatestUpdate = 0.0;

            for (U32 index = 0U; index < count; index++) {
                const double magnitude = Abs(trialAngles[index] - angles[index]);

                greatestUpdate = Max(greatestUpdate, magnitude);
            }

            const bool trialEvaluated = S9EvaluateKnotSystem(system,
                trialAngles, true, trialStartCurvatures, trialEndCurvatures,
                trialDeltas);

            if (trialEvaluated) {
                const double trialNorm = S9FreeDeltaNorm(system,
                    unknownCount, trialDeltas);
                const bool decreases = trialNorm <=
                    deltaNorm * (1.0 - S9CoupledDecrease * factor);

                // The Newton step itself must be tiny, not the movement
                // after the cap and the clamp, or a held knot looks
                // converged.
                const bool isTiny = greatestStep < S9MaximumITPError;

                if (decreases or isTiny) {
                    isAccepted = true;
                    break;
                }
            }

            factor *= 0.5;
        }

        if (!isAccepted) {
            break;
        }

        for (U32 index = 0U; index < count; index++) {
            angles[index] = trialAngles[index];
        }

        for (U32 segment = 0U; segment < segmentCount; segment++) {
            startCurvatures[segment] = trialStartCurvatures[segment];
            endCurvatures[segment] = trialEndCurvatures[segment];
        }

        for (U32 unknown = 0U; unknown < unknownCount; unknown++) {
            deltas[unknown] = trialDeltas[unknown];
        }

        if (greatestStep < S9MaximumITPError) {
            isComplete = true;
            break;
        }
    }

    // The segment curves can hold a rejected trial. Match them to the
    // accepted angles.
    if (isEvaluated) {
        S9EvaluateKnotSystem(system, angles, true, startCurvatures,
            endCurvatures, deltas);
    }

    Dealloc(trialAngles);
    Dealloc(highs);
    Dealloc(lows);
    Dealloc(step);
    Dealloc(rightSide);
    Dealloc(superdiagonal);
    Dealloc(diagonal);
    Dealloc(subdiagonal);
    Dealloc(trialDeltas);
    Dealloc(deltas);
    Dealloc(trialEndCurvatures);
    Dealloc(trialStartCurvatures);
    Dealloc(endCurvatures);
    Dealloc(startCurvatures);

    return isComplete;
}


/**
 * Solves the knot angles of a spline so that the curvature is the same on
 * both sides of every knot. It solves every unknown segment, makes sweeps
 * of S9SweepKnotNewtonSteps until the greatest update of a sweep is less
 * than S9CoupledStartError, and then calls S9SolveKnotSystemCoupled.
 *
 * Returns true when the knot angles converged. False when a step failed.
 * Knot angles then hold the last accepted iterate, which is a better
 * start for the sweeps of S9SolveKnotAngle than the initial angles.
 */
static bool S9SolveKnotSystem(const S9KnotSystem &system) noexcept
{
    ASSERT(system.Count >= 2U);
    ASSERT(system.SegmentAngles != nullptr);
    ASSERT(system.SegmentLengths != nullptr);
    ASSERT(system.KnotAngles != nullptr);
    ASSERT(system.SegmentCurves != nullptr);

    S9SolveUnknownSegments(system);

    double *const trial = AllocElements<double>(system.Count);
    double greatestUpdate = 1e300;
    U32 sweep = 0U;
    bool isComplete = false;
    bool isFailed = false;

    while (sweep < S9KnotSweepCount and
        greatestUpdate >= S9CoupledStartError) {
        greatestUpdate = S9SweepKnotNewtonSteps(system, trial);
        sweep++;

        if (greatestUpdate < 0.0) {
            isFailed = true;
            break;
        }

        if (greatestUpdate < S9MaximumITPError) {
            isComplete = true;
            break;
        }
    }

    Dealloc(trial);

    if (isFailed) {
        return false;
    }

    if (isComplete) {
        return true;
    }

    return S9SolveKnotSystemCoupled(system);
}


void S9SolveAcyclicExplicit(const U32 count,
    const FloatPoint * const knots, const S9SplineEndControls &controls,
    CurveArray<CubicCurve> &curves) noexcept
{
    const U32 segmentCount = count - 1;
    const U32 segmentCurveCount = segmentCount;

    double *segmentLengths = AllocElements<double>(segmentCount);
    double *segmentAngles = AllocElements<double>(segmentCount);

    for (U32 i = 0; i < segmentCount; i++) {
        const FloatPoint current = knots[i];
        const FloatPoint next = knots[i + 1];
        const FloatPoint d = next - current;

        segmentLengths[i] = d.Length();

        segmentAngles[i] = S9VectorAngle(d);
    }

    const U32 n = count - 1;

    double *knotAngles = AllocElements<double>(count);

    S9SegmentCurve *const segmentCurves =
        AllocElementsInit<S9SegmentCurve>(segmentCurveCount);

    bool *const pinned = AllocElements<bool>(count);

    for (U32 i = 0; i < count; i++) {
        pinned[i] = false;
    }

    S9KnotSystem system;

    system.Count = count;
    system.SegmentCount = segmentCount;
    system.IsCyclic = false;
    system.SegmentAngles = segmentAngles;
    system.SegmentLengths = segmentLengths;
    system.KnotAngles = knotAngles;
    system.SegmentCurves = segmentCurves;
    system.Pinned = pinned;

    S9ApplyEndControls(controls, knots, system);

    for (U32 i = 1; i < n; i++) {
        const int prevIndex = i - 1;
        const int nextIndex = i + 1;

        const FloatPoint current = knots[i];

        const FloatPoint prevToCurrent = current - knots[prevIndex];
        const FloatPoint currentToNext = knots[nextIndex] - current;

        const double a0 = S9VectorAngle(prevToCurrent);
        const double a1 = S9VectorAngle(currentToNext);
        const double l0 = prevToCurrent.Length();
        const double l1 = currentToNext.Length();
        const double combinedLengths = l0 + l1;

        ASSERT(combinedLengths > DBL_EPSILON);

        const double bend = a1 - a0;

        knotAngles[i] = a0 + ((bend * l0) / combinedLengths);
    }

    knotAngles[0] = S9AlphaToAbsolute(0, segmentAngles[0]);

    knotAngles[n] = S9BetaToAbsolute(0, segmentAngles[segmentCount - 1]);

    for (U32 i = 0; i < S9MaximumCenteringIterationCount; i++) {
        knotAngles[0] = S9StartKnotAngle(system, knotAngles);

        for (U32 j = 1; j < n; j++) {
            const double prevKnotAngle = knotAngles[j - 1];
            const double currentKnotAngle = knotAngles[j];
            const double nextKnotAngle = knotAngles[j + 1];
            const double prevSegmentAngle = segmentAngles[j - 1];
            const double nextSegmentAngle = segmentAngles[j];

            knotAngles[j] = S9ClampKnotAngle(prevKnotAngle, currentKnotAngle,
                nextKnotAngle, prevSegmentAngle, nextSegmentAngle);
        }

        knotAngles[n] = S9EndKnotAngle(system, knotAngles);
    }

    // Updates of the previous sweep, for oscillation damping.
    double *const previousUpdates = AllocElements<double>(count);

    for (U32 i = 0; i < count; i++) {
        previousUpdates[i] = 0;
    }

    // The stored sign of the delta near the start of the rotation range of
    // each knot.
    I8 *const lowSigns = AllocElements<I8>(count);

    for (U32 i = 0; i < count; i++) {
        lowSigns[i] = 0;
    }

    const bool isConverged = S9SolveKnotSystem(system);

    // Sweeps of S9SolveKnotAngle when the Newton steps did not converge.
    for (U32 i = 0 ; i < S9MaximumITPIterationCount and !isConverged; i++) {
        S9ClassifyKnots(system, knotAngles);

        const double angle0 = S9StartKnotAngle(system, knotAngles);

        double maxError = Abs(knotAngles[0] - angle0);

        knotAngles[0] = angle0;

        for (U32 j = 1; j < n; j++) {
            const double prevKnotAngle = knotAngles[j - 1];
            const double nextKnotAngle = knotAngles[j + 1];
            const double prevSegmentAngle = segmentAngles[j - 1];
            const double nextSegmentAngle = segmentAngles[j];
            const double prevSegmentLength = segmentLengths[j - 1];
            const double nextSegmentLength = segmentLengths[j];

            double solvedAngle = 0.0;

            if (system.Pinned[j]) {
                solvedAngle = S9ClampKnotAngle(prevKnotAngle, knotAngles[j],
                    nextKnotAngle, prevSegmentAngle, nextSegmentAngle);
            } else {
                solvedAngle = S9SolveKnotAngle(prevKnotAngle,
                    knotAngles[j], nextKnotAngle, prevSegmentAngle,
                    nextSegmentAngle, prevSegmentLength, nextSegmentLength,
                    previousUpdates[j], segmentCurves[j - 1],
                    segmentCurves[j], lowSigns[j]);
            }

            const double angle = S9CalculateDampedKnotAngle(knotAngles[j],
                solvedAngle, previousUpdates[j]);

            maxError = Max(maxError, Abs(knotAngles[j] - angle));

            knotAngles[j] = angle;
        }

        const double angle1 = S9EndKnotAngle(system, knotAngles);

        maxError = Max(maxError, Abs(knotAngles[n] - angle1));

        knotAngles[n] = angle1;

        if (maxError < S9MaximumITPError) {
            break;
        }
    }

    Dealloc(lowSigns);
    Dealloc(previousUpdates);

    curves.Restart(knots[0]);

    CurveArray<CubicCurve> pieces;
    bool splineIsComplete = true;

    for (U32 i = 0; i < n; i++) {
        const FloatPoint p0 = knots[i];
        const FloatPoint p1 = knots[i + 1];

        const double segmentAngle = segmentAngles[i];

        const double a0 = S9AbsoluteToAlpha(knotAngles[i], segmentAngle);
        const double a1 = S9AbsoluteToBeta(knotAngles[i + 1],
            segmentAngle);

        const double segmentLength = segmentLengths[i];

        const bool appended = S9AppendCubics(a0, a1, p0, p1,
            segmentAngle, segmentLength, &segmentCurves[i], pieces, curves);

        if (!appended) {
            splineIsComplete = false;
            break;
        }
    }

    // Spline with a missing segment is not usable.
    if (!splineIsComplete) {
        curves.Erase();
    }

    Dealloc(pinned);
    Dealloc(segmentCurves);
    Dealloc(knotAngles);
    Dealloc(segmentAngles);
    Dealloc(segmentLengths);
}


void S9SolveAcyclic(const U32 count, const FloatPoint * const knots,
    CurveArray<CubicCurve> &curves) noexcept
{
    const S9SplineEndControls controls;

    S9SolveAcyclicExplicit(count, knots, controls, curves);
}


void S9SolveCyclic(const U32 count, const FloatPoint * const knots,
    CurveArray<CubicCurve> &curves) noexcept
{
    const U32 firstKnotIndex = 0;
    const U32 firstSegmentIndex = 0;
    const U32 finalKnotIndex = count - 1;
    const U32 finalSegmentIndex = count - 1;
    const U32 segmentCurveCount = count;

    double *segmentLengths = AllocElements<double>(count);
    double *segmentAngles = AllocElements<double>(count);

    for (U32 i = 0; i < count - 1; i++) {
        const FloatPoint current = knots[i];
        const FloatPoint next = knots[i + 1];
        const FloatPoint d = next - current;

        segmentLengths[i] = d.Length();

        segmentAngles[i] = S9VectorAngle(d);
    }

    // The last segment goes from the last knot to the first one.
    const FloatPoint current = knots[finalKnotIndex];
    const FloatPoint next = knots[firstKnotIndex];
    const FloatPoint d = next - current;

    segmentLengths[finalSegmentIndex] = d.Length();

    segmentAngles[finalSegmentIndex] = S9VectorAngle(d);

    double *knotAngles = AllocElements<double>(count);

    for (U32 i = 1; i < count - 1; i++) {
        const int prevIndex = i - 1;
        const int nextIndex = i + 1;

        const FloatPoint current = knots[i];

        const FloatPoint prevToCurrent = current - knots[prevIndex];
        const FloatPoint currentToNext = knots[nextIndex] - current;

        const double a0 = S9VectorAngle(prevToCurrent);
        const double a1 = S9VectorAngle(currentToNext);
        const double l0 = prevToCurrent.Length();
        const double l1 = currentToNext.Length();
        const double combinedLengths = l0 + l1;

        ASSERT(combinedLengths > DBL_EPSILON);

        const double bend = a1 - a0;

        knotAngles[i] = a0 + ((bend * l0) / combinedLengths);
    }

    knotAngles[firstKnotIndex] = S9AlphaToAbsolute(0, segmentAngles[firstSegmentIndex]);

    knotAngles[finalKnotIndex] = S9BetaToAbsolute(0, segmentAngles[finalSegmentIndex]);

    for (U32 i = 0; i < S9MaximumCenteringIterationCount; i++) {
        knotAngles[firstKnotIndex] = S9ClampKnotAngle(
            knotAngles[finalKnotIndex], knotAngles[firstKnotIndex], knotAngles[1],
            segmentAngles[finalSegmentIndex], segmentAngles[firstSegmentIndex]);

        for (U32 j = 1; j < count - 1; j++) {
            const double prevKnotAngle = knotAngles[j - 1];
            const double currentKnotAngle = knotAngles[j];
            const double nextKnotAngle = knotAngles[j + 1];
            const double prevSegmentAngle = segmentAngles[j - 1];
            const double nextSegmentAngle = segmentAngles[j];

            knotAngles[j] = S9ClampKnotAngle(prevKnotAngle, currentKnotAngle,
                nextKnotAngle, prevSegmentAngle, nextSegmentAngle);
        }

        knotAngles[finalKnotIndex] = S9ClampKnotAngle(
            knotAngles[finalKnotIndex - 1], knotAngles[finalKnotIndex],
            knotAngles[firstKnotIndex],
            segmentAngles[finalSegmentIndex - 1],
            segmentAngles[finalSegmentIndex]);
    }

    // Updates of the previous sweep, for oscillation damping.
    double *const previousUpdates = AllocElements<double>(count);

    for (U32 i = 0; i < count; i++) {
        previousUpdates[i] = 0;
    }

    // Last solved curve of each segment, and the kept sign of the delta near
    // the start of the rotation range of each knot.
    S9SegmentCurve *const segmentCurves =
        AllocElementsInit<S9SegmentCurve>(segmentCurveCount);
    I8 *const lowSigns = AllocElements<I8>(count);

    for (U32 i = 0; i < count; i++) {
        lowSigns[i] = 0;
    }

    bool *const pinned = AllocElements<bool>(count);

    for (U32 i = 0; i < count; i++) {
        pinned[i] = false;
    }

    S9KnotSystem system;

    system.Count = count;
    system.SegmentCount = count;
    system.IsCyclic = true;
    system.SegmentAngles = segmentAngles;
    system.SegmentLengths = segmentLengths;
    system.KnotAngles = knotAngles;
    system.SegmentCurves = segmentCurves;
    system.Pinned = pinned;

    const bool isConverged = S9SolveKnotSystem(system);

    // Sweeps of S9SolveKnotAngle when the Newton steps did not converge.
    for (U32 i = 0 ; i < S9MaximumITPIterationCount and !isConverged; i++) {
        S9ClassifyKnots(system, knotAngles);

        const double solvedAngle0 = S9SolveKnotAngle(knotAngles[finalKnotIndex],
            knotAngles[firstKnotIndex], knotAngles[1],
            segmentAngles[finalSegmentIndex], segmentAngles[firstSegmentIndex],
            segmentLengths[finalSegmentIndex], segmentLengths[firstSegmentIndex],
            previousUpdates[firstKnotIndex], segmentCurves[finalSegmentIndex],
            segmentCurves[firstSegmentIndex], lowSigns[firstKnotIndex]);

        const double angle0 = S9CalculateDampedKnotAngle(
            knotAngles[firstKnotIndex], solvedAngle0,
            previousUpdates[firstKnotIndex]);

        double maxError = Abs(knotAngles[firstKnotIndex] - angle0);

        knotAngles[0] = angle0;

        for (U32 j = 1; j < count - 1; j++) {
            const double prevKnotAngle = knotAngles[j - 1];
            const double nextKnotAngle = knotAngles[j + 1];
            const double prevSegmentAngle = segmentAngles[j - 1];
            const double nextSegmentAngle = segmentAngles[j];
            const double prevSegmentLength = segmentLengths[j - 1];
            const double nextSegmentLength = segmentLengths[j];

            double solvedAngle = 0.0;

            if (system.Pinned[j]) {
                solvedAngle = S9ClampKnotAngle(prevKnotAngle, knotAngles[j],
                    nextKnotAngle, prevSegmentAngle, nextSegmentAngle);
            } else {
                solvedAngle = S9SolveKnotAngle(prevKnotAngle,
                    knotAngles[j], nextKnotAngle, prevSegmentAngle,
                    nextSegmentAngle, prevSegmentLength, nextSegmentLength,
                    previousUpdates[j], segmentCurves[j - 1],
                    segmentCurves[j], lowSigns[j]);
            }

            const double angle = S9CalculateDampedKnotAngle(knotAngles[j],
                solvedAngle, previousUpdates[j]);

            maxError = Max(maxError, Abs(knotAngles[j] - angle));

            knotAngles[j] = angle;
        }

        const double solvedAngle1 = S9SolveKnotAngle(
            knotAngles[finalKnotIndex - 1], knotAngles[finalKnotIndex],
            knotAngles[firstKnotIndex], segmentAngles[finalSegmentIndex - 1],
            segmentAngles[finalSegmentIndex],
            segmentLengths[finalSegmentIndex - 1],
            segmentLengths[finalSegmentIndex], previousUpdates[finalKnotIndex],
            segmentCurves[finalSegmentIndex - 1],
            segmentCurves[finalSegmentIndex], lowSigns[finalKnotIndex]);

        const double angle1 = S9CalculateDampedKnotAngle(
            knotAngles[finalKnotIndex], solvedAngle1,
            previousUpdates[finalKnotIndex]);

        maxError = Max(maxError, Abs(knotAngles[finalKnotIndex] - angle1));

        knotAngles[finalKnotIndex] = angle1;

        if (maxError < S9MaximumITPError) {
            break;
        }
    }

    Dealloc(lowSigns);
    Dealloc(previousUpdates);

    curves.Restart(knots[firstKnotIndex]);

    CurveArray<CubicCurve> pieces;
    bool splineIsComplete = true;

    for (U32 i = 0; i < count - 1; i++) {
        const FloatPoint p0 = knots[i];
        const FloatPoint p1 = knots[i + 1];

        const double segmentAngle = segmentAngles[i];

        const double a0 = S9AbsoluteToAlpha(knotAngles[i], segmentAngle);
        const double a1 = S9AbsoluteToBeta(knotAngles[i + 1],
            segmentAngle);

        const double segmentLength = segmentLengths[i];

        const bool appended = S9AppendCubics(a0, a1, p0, p1,
            segmentAngle, segmentLength, &segmentCurves[i], pieces, curves);

        if (!appended) {
            splineIsComplete = false;
            break;
        }
    }

    // Close.
    if (splineIsComplete) {
        const FloatPoint p0 = knots[finalKnotIndex];
        const FloatPoint p1 = knots[firstKnotIndex];

        const double segmentAngle = segmentAngles[finalSegmentIndex];

        const double a0 = S9AbsoluteToAlpha(knotAngles[finalKnotIndex],
            segmentAngle);
        const double a1 = S9AbsoluteToBeta(knotAngles[firstKnotIndex],
            segmentAngle);

        const double segmentLength = segmentLengths[finalSegmentIndex];

        splineIsComplete = S9AppendCubics(a0, a1, p0, p1,
            segmentAngle, segmentLength, &segmentCurves[finalSegmentIndex],
            pieces, curves);
    }

    // A spline with a missing segment is not usable.
    if (!splineIsComplete) {
        curves.Erase();
    }

    Dealloc(pinned);
    Dealloc(segmentCurves);
    Dealloc(knotAngles);
    Dealloc(segmentAngles);
    Dealloc(segmentLengths);
}
