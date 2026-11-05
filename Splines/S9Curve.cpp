
#include "GaussLegendre.h"
#include "S9Curve.h"


/*
 * The corner law has closed form integrals, thus the scalar solver uses
 * exact tangent angle. The unit corner term at the normalized arc length r
 * from its end is
 *
 *     x = λ r
 *     c(r) = Q (p(x) + φ q(x))
 *     p(x) = 1 / (1 + x)^1.5
 *     q(x) = η x / (1 + x)^3
 *     φ = λ / (λ + λ_c)
 *
 * Scale Q makes the integral over the curve equal to λ. Both terms
 * and their integrals are polynomials in u = 1 / sqrt(1 + x), therefore the
 * angle law needs no quadrature. The term decreases from its end when
 * η φ is at most 1.5.
 */


/**
 * Returns the mean of the unit corner integrand, corner term without its
 * scale Q, over an interval of scaled distances -
 *
 *     u1 = 1 / sqrt(1 + x1), u2 = 1 / sqrt(1 + x2)
 *     mean = u1^2 u2^2 (2 / (u1 + u2) + η φ / 2 (x1 u1^2 + x2 u2^2))
 *
 * @param nearDistance Scaled distance x1 of the end nearer the corner, at
 * least zero.
 *
 * @param farDistance Scaled distance x2 of the far end, at least
 * nearDistance.
 *
 * @param rounding Weight φ of the rounding term, from zero through one.
 *
 * @param eta η of the rounding term, at least zero.
 *
 * @return Positive mean.
 */
static inline double S9UnitCornerMean(const double nearDistance,
    const double farDistance, const double rounding,
    const double eta) noexcept
{
    ASSERT(nearDistance >= 0.0);
    ASSERT(farDistance >= nearDistance);
    ASSERT(rounding >= 0.0);
    ASSERT(rounding <= 1.0);
    ASSERT(eta >= 0.0);

    const double nearRoot = 1.0 / sqrt(1.0 + nearDistance);
    const double farRoot = 1.0 / sqrt(1.0 + farDistance);
    const double product = nearRoot * farRoot;
    const double nearSquare = nearRoot * nearRoot;
    const double farSquare = farRoot * farRoot;
    const double squareProduct = product * product;
    const double poleMean = 2.0 / (nearRoot + farRoot);
    const double roundingMean = 0.5 * eta * rounding *
        (nearDistance * nearSquare + farDistance * farSquare);
    const double mean = squareProduct * (poleMean + roundingMean);
    ASSERT(mean > 0.0);

    return mean;
}


/**
 * Returns the scale Q of the unit corner term for one corner λ.
 *
 * @param lambda Corner λ, at least zero. Zero gives zero, an end without
 * a corner.
 *
 * @param eta η of the rounding term, at least zero.
 */
static inline double S9UnitCornerScale(const double lambda,
    const double eta) noexcept
{
    ASSERT(lambda >= 0.0);
    ASSERT(eta >= 0.0);

    const double rounding = S9CornerRounding(lambda);
    const double mean = S9UnitCornerMean(0.0, lambda, rounding, eta);

    return lambda / mean;
}


/**
 * Returns the unit corner term at one distance from its end, for a scale
 * that the caller computed once with S9UnitCornerScale, because the
 * scale needs a square root and a division.
 *
 * @param lambda Corner λ, at least zero.
 *
 * @param eta η of the rounding term, at least zero.
 *
 * @param scale Value of S9UnitCornerScale for λ and η.
 *
 * @param distance Normalized arc length r from the end of the corner, from
 * zero to one.
 *
 * @return Turn rate of the term, at least zero.
 */
static inline double S9UnitCornerValueWithScale(const double lambda,
    const double eta, const double scale, const double distance) noexcept
{
    ASSERT(lambda >= 0.0);
    ASSERT(eta >= 0.0);
    ASSERT(scale >= 0.0);
    ASSERT(distance >= 0.0);
    ASSERT(distance <= 1.0);

    const double scaled = distance * lambda;
    const double rounding = S9CornerRounding(lambda);
    const double unitRoot = 1.0 / sqrt(1.0 + scaled);
    const double pole = unitRoot * unitRoot * unitRoot;
    const double roundingTerm = eta * rounding * scaled * pole * pole;

    return scale * (pole + roundingTerm);
}


/**
 * Returns the unit corner term at one distance from its end, see
 * S9UnitCornerValueWithScale.
 */
static inline double S9UnitCornerValue(const double lambda,
    const double eta, const double distance) noexcept
{
    const double scale = S9UnitCornerScale(lambda, eta);

    return S9UnitCornerValueWithScale(lambda, eta, scale, distance);
}


/**
 * Returns the derivative of the unit corner term with respect to the
 * distance from its end. The derivative is
 *
 *     c'(r) = Q λ (-1.5 + η φ (1 - 2 x) / (1 + x)^1.5) / (1 + x)^2.5
 *
 * At the end the slope is Q λ (η φ - 1.5). It is negative at every
 * distance when η φ is at most 1.5, and past half a pole distance for
 * every η.
 *
 * @param lambda Corner λ, at least zero.
 *
 * @param eta η of the rounding term, at least zero.
 *
 * @param distance Normalized arc length r from the end of the corner, from
 * zero to one.
 */
static inline double S9UnitCornerSlope(const double lambda,
    const double eta, const double distance) noexcept
{
    ASSERT(lambda >= 0.0);
    ASSERT(eta >= 0.0);
    ASSERT(distance >= 0.0);
    ASSERT(distance <= 1.0);

    const double scale = S9UnitCornerScale(lambda, eta);
    const double scaled = distance * lambda;
    const double rounding = S9CornerRounding(lambda);
    const double unitRoot = 1.0 / sqrt(1.0 + scaled);
    const double unitSquare = unitRoot * unitRoot;
    const double unitCube = unitSquare * unitRoot;
    const double unitFifth = unitSquare * unitCube;
    const double roundingSlope = eta * rounding * unitCube *
        (1.0 - 2.0 * scaled);

    return scale * lambda * unitFifth * (-1.5 + roundingSlope);
}


static inline double S9UnitCornerTurnFromEndWithScale(const double lambda,
    const double eta, const double scale, const double distance) noexcept
{
    ASSERT(lambda >= 0.0);
    ASSERT(eta >= 0.0);
    ASSERT(scale >= 0.0);
    ASSERT(distance >= 0.0);
    ASSERT(distance <= 1.0);

    const double scaledDistance = distance * lambda;
    const double rounding = S9CornerRounding(lambda);
    const double mean = S9UnitCornerMean(0.0, scaledDistance, rounding, eta);

    return scale * distance * mean;
}


static inline double S9UnitCornerTurnFromEnd(const double lambda,
    const double eta, const double distance) noexcept
{
    const double scale = S9UnitCornerScale(lambda, eta);

    return S9UnitCornerTurnFromEndWithScale(lambda, eta, scale, distance);
}


static inline double S9UnitCornerTurnToFarEndWithScale(const double lambda,
    const double eta, const double scale, const double width) noexcept
{
    ASSERT(lambda >= 0.0);
    ASSERT(eta >= 0.0);
    ASSERT(scale >= 0.0);
    ASSERT(width >= 0.0);
    ASSERT(width <= 1.0);

    const double nearDistance = (1.0 - width) * lambda;
    const double rounding = S9CornerRounding(lambda);
    const double mean = S9UnitCornerMean(nearDistance, lambda, rounding, eta);

    return scale * width * mean;
}


static inline double S9UnitCornerTurnToFarEnd(const double lambda,
    const double eta, const double width) noexcept
{
    const double scale = S9UnitCornerScale(lambda, eta);

    return S9UnitCornerTurnToFarEndWithScale(lambda, eta, scale, width);
}


double S9Curve::StartBasis(const double t) const noexcept
{
    ASSERT(t >= 0.0);
    ASSERT(t <= 1.0);

    const double corner = S9UnitCornerValue(StartLambda, StartEta, t);

    return (1.0 - t) + corner;
}


double S9Curve::EndBasis(const double t) const noexcept
{
    ASSERT(t >= 0.0);
    ASSERT(t <= 1.0);

    const double corner = S9UnitCornerValue(EndLambda, EndEta, 1.0 - t);

    return t + corner;
}


double S9Curve::StartBasisIntegral(const double t) const noexcept
{
    ASSERT(t >= 0.0);
    ASSERT(t <= 1.0);

    const double corner = S9UnitCornerTurnFromEnd(StartLambda, StartEta,
        t);

    return t - 0.5 * t * t + corner;
}


double S9Curve::EndBasisIntegral(const double t) const noexcept
{
    ASSERT(t >= 0.0);
    ASSERT(t <= 1.0);

    const double corner = S9UnitCornerTurnToFarEnd(EndLambda, EndEta, t);

    return 0.5 * t * t + corner;
}


double S9Curve::TurnRate(const double t) const noexcept
{
    const double startBasis = StartBasis(t);
    const double endBasis = EndBasis(t);

    return StartLine * startBasis + EndLine * endBasis;
}


double S9Curve::TurnRateSlope(const double t) const noexcept
{
    ASSERT(t >= 0.0);
    ASSERT(t <= 1.0);

    const double startCorner = S9UnitCornerSlope(StartLambda, StartEta,
        t);
    const double endCorner = S9UnitCornerSlope(EndLambda, EndEta,
        1.0 - t);
    const double startSlope = -1.0 + startCorner;
    const double endSlope = 1.0 - endCorner;

    return StartLine * startSlope + EndLine * endSlope;
}


double S9Curve::TurnFromStart(const double t) const noexcept
{
    const double startIntegral = StartBasisIntegral(t);
    const double endIntegral = EndBasisIntegral(t);

    return StartLine * startIntegral + EndLine * endIntegral;
}


/*
 * Smallest magnitude of a weight. A smaller magnitude is raised to it, and
 * the corner is then numerically infinite with a pole distance of about
 * 1e-8.
 */
static constexpr double S9SmallestWeight = 1e-8;


/*
 * Limits of the quadrature grid.
 */
static constexpr U32 S9BasePanelCount = 8U;
static constexpr U32 S9CoarseBasePanelCount = 4U;
static constexpr U32 S9GreatestGradedEdgeCount = 30U;
static constexpr U32 S9GreatestEdgeCount = S9BasePanelCount + 1U +
    2U * S9GreatestGradedEdgeCount;
static constexpr U32 S9GreatestNodeCount =
    (S9GreatestEdgeCount - 1U) * 16U;
static constexpr U32 S9GreatestCoarseNodeCount =
    (S9CoarseBasePanelCount + 2U * S9GreatestGradedEdgeCount) * 8U;

/*
 * Limits of the root search.
 */
static constexpr U32 S9ScanSampleCount = 1440U;
static constexpr double S9ScanReach = 6.0 * Pi;
static constexpr double S9RootTolerance = 1e-12;
static constexpr U32 S9RefineCount = 64U;
static constexpr double S9SameTurnTolerance = 1e-9;
static constexpr double S9SmallestChordLength = 1e-9;

static constexpr U32 S9ScanRestartInterval = 32U;


/*
 * Limits of cubic approximation.
 */
static constexpr double S9GreatestPanelTurn = 0.25;
static constexpr double S9IntegrationAllowance = 1e-12;
static constexpr U32 S9GreatestDepth = 34U;
static constexpr U32 S9GreatestCubicCount = 2048U;
static constexpr U32 S9GreatestPanelCount = 60000U;


STATIC_ASSERT(S9SmallestWeight > 0.0);
STATIC_ASSERT(S9SmallestWeight < 1.0);
STATIC_ASSERT(S9SmallestWeight * 134217728.0 > 0.5);
STATIC_ASSERT(S9ScanSampleCount > 2U);
STATIC_ASSERT(S9ScanSampleCount % 2U == 0U);
STATIC_ASSERT(S9ScanRestartInterval > 0U);
STATIC_ASSERT(S9ScanRestartInterval <= S9ScanSampleCount);
STATIC_ASSERT(S9GreatestCoarseNodeCount <= S9GreatestNodeCount);
STATIC_ASSERT(S9RootTolerance > 0.0);
STATIC_ASSERT(S9GreatestPanelTurn > 0.0);
STATIC_ASSERT(S9CoarseBasePanelCount <= S9BasePanelCount);



static constexpr U32 S9GreatestNearStepCount = 12U;
static constexpr double S9NearReach = 2.0;
static constexpr double S9GreatestNearAngleChange = 0.5;


static double S9ReduceAngle(const double angle) noexcept
{
    ASSERT(isfinite(angle));

    static constexpr double fullTurn = 2.0 * Pi;

    return remainder(angle, fullTurn);
}


static double S9Heading(const double forwardAngle,
    const double weight) noexcept
{
    ASSERT(isfinite(forwardAngle));
    ASSERT(isfinite(weight));

    double heading = forwardAngle;

    if (weight < 0.0) {
        heading += Pi;
    }

    return S9ReduceAngle(heading);
}


static double S9WeightLambda(const double weight) noexcept
{
    ASSERT(isfinite(weight));

    const double magnitude = Abs(weight);
    const double limited = Clamp(magnitude, S9SmallestWeight, 1.0);
    const double lambda = (1.0 - limited) / limited;

    ASSERT(lambda >= 0.0);
    ASSERT(isfinite(lambda));

    return lambda;
}


/**
 * Sets the angles, the corner λ, the rounding η and the net turn
 * of a shape from its parameters.
 */
static void S9PrepareShape(const S9CurveParameters &parameters,
    S9Curve &shape) noexcept
{
    ASSERT(isfinite(parameters.Alpha));
    ASSERT(isfinite(parameters.Beta));
    ASSERT(parameters.StartWeight >= -1.0);
    ASSERT(parameters.StartWeight <= 1.0);
    ASSERT(parameters.EndWeight >= -1.0);
    ASSERT(parameters.EndWeight <= 1.0);
    ASSERT(isfinite(parameters.StartEta));
    ASSERT(isfinite(parameters.EndEta));
    ASSERT(parameters.StartEta >= 0.0);
    ASSERT(parameters.EndEta >= 0.0);

    shape.StartAngle = S9Heading(parameters.Alpha, parameters.StartWeight);
    shape.EndAngle = S9Heading(-parameters.Beta, parameters.EndWeight);
    shape.StartLambda = S9WeightLambda(parameters.StartWeight);
    shape.EndLambda = S9WeightLambda(parameters.EndWeight);
    shape.StartEta = parameters.StartEta;
    shape.EndEta = parameters.EndEta;
    shape.Turn = S9ReduceAngle(shape.EndAngle - shape.StartAngle);
}


/**
 * Returns raw net turn of a shape.
 */
static inline constexpr double S9RawTurn(const S9Curve &shape) noexcept
{
    return shape.EndAngle - shape.StartAngle;
}


static void S9SetLines(const double share,
    S9Curve &shape) noexcept
{
    ASSERT(isfinite(share));
    ASSERT(isfinite(shape.Turn));
    ASSERT(shape.StartLambda >= 0.0);
    ASSERT(shape.EndLambda >= 0.0);

    const double startBasisTurn = shape.StartBasisTurn();
    const double endBasisTurn = shape.EndBasisTurn();

    shape.StartLine = share / startBasisTurn;
    shape.EndLine = (shape.Turn - share) / endBasisTurn;
}


/**
 * Returns start share of the net turn of a solved shape.
 */
static inline constexpr double S9StartShare(const S9Curve &shape) noexcept
{
    const double startBasisTurn = shape.StartBasisTurn();

    return shape.StartLine * startBasisTurn;
}


/**
 * Quadrature grid of one shape.
 */
struct S9Grid final {
    S9Grid() noexcept;

    // Number of valid nodes.
    U32 Count = 0U;

    // True when StartRate and EndRate are valid.
    bool HasRates = false;

    // Quadrature weight of each node.
    double Weight[S9GreatestNodeCount];

    // Start basis integral G0 at each node.
    double StartIntegral[S9GreatestNodeCount];

    // End basis integral G1 at each node.
    double EndIntegral[S9GreatestNodeCount];

    // Rate of change of the tangent angle at each node with the start
    // share - StartIntegral / J0 - EndIntegral / J1.
    double ShareRate[S9GreatestNodeCount];

    // Start basis curvature g0 at each node. Valid when HasRates is true.
    double StartRate[S9GreatestNodeCount];

    // End basis curvature g1 at each node. Valid when HasRates is true.
    double EndRate[S9GreatestNodeCount];
};


inline S9Grid::S9Grid() noexcept {
}


static inline constexpr void S9AddGradedEdges(const double lambda,
    const bool fromEnd, double (&edges)[S9GreatestEdgeCount],
    U32 &count) noexcept
{
    ASSERT(lambda >= 0.0);
    ASSERT(count <= S9GreatestEdgeCount);

    if (lambda <= 0.0) {
        return;
    }

    double distance = 0.125 / lambda;
    U32 added = 0U;

    while (distance < 0.5 and added < S9GreatestGradedEdgeCount and
        count < S9GreatestEdgeCount) {
        double edge = distance;

        if (fromEnd) {
            edge = 1.0 - distance;
        }

        edges[count] = edge;
        count++;
        added++;
        distance *= 2.0;
    }
}


static void S9BuildGrid(const S9Curve &shape,
    const U32 basePanelCount, const U32 nodesPerPanel, const bool withRates,
    S9Grid &grid) noexcept
{
    ASSERT(basePanelCount >= 1U);
    ASSERT(basePanelCount <= S9BasePanelCount);
    ASSERT(nodesPerPanel == 8U or nodesPerPanel == 16U);
    ASSERT(shape.StartLambda >= 0.0);
    ASSERT(shape.EndLambda >= 0.0);
    ASSERT(shape.StartEta >= 0.0);
    ASSERT(shape.EndEta >= 0.0);

    double edges[S9GreatestEdgeCount];
    U32 edgeCount = 0U;

    for (U32 index = 0U; index <= basePanelCount; index++) {
        edges[edgeCount] = double(index) / double(basePanelCount);
        edgeCount++;
    }

    S9AddGradedEdges(shape.StartLambda, false, edges, edgeCount);
    S9AddGradedEdges(shape.EndLambda, true, edges, edgeCount);

    // Insertion sort. The list has at most S9GreatestEdgeCount edges.
    for (U32 index = 1U; index < edgeCount; index++) {
        const double edge = edges[index];
        U32 position = index;

        while (position > 0U and edges[position - 1U] > edge) {
            edges[position] = edges[position - 1U];
            position--;
        }

        edges[position] = edge;
    }

    const double *nodeX = GaussLegendreX_16;
    const double *nodeScale = GaussLegendreScale_16;

    if (nodesPerPanel == 8U) {
        nodeX = GaussLegendreX_8;
        nodeScale = GaussLegendreScale_8;
    }

    const S9Curve local = shape;
    const double startLambda = local.StartLambda;
    const double endLambda = local.EndLambda;
    const double startEta = local.StartEta;
    const double endEta = local.EndEta;
    const double startScale = S9UnitCornerScale(startLambda, startEta);
    const double endScale = S9UnitCornerScale(endLambda, endEta);
    const double startBasisTurn = local.StartBasisTurn();
    const double endBasisTurn = local.EndBasisTurn();
    U32 count = 0U;

    for (U32 index = 1U; index < edgeCount; index++) {
        const double begin = edges[index - 1U];
        const double end = edges[index];
        const double width = end - begin;

        // A duplicate edge gives an empty panel, which adds nothing.
        if (width <= 0.0) {
            continue;
        }

        for (U32 node = 0U; node < nodesPerPanel; node++) {
            const U32 slot = count;

            ASSERT(slot < S9GreatestNodeCount);

            const double t = begin + width * nodeX[node];
            const double startCorner = S9UnitCornerTurnFromEndWithScale(
                startLambda, startEta, startScale, t);
            const double endCorner = S9UnitCornerTurnToFarEndWithScale(
                endLambda, endEta, endScale, t);
            const double startIntegral = t - 0.5 * t * t + startCorner;
            const double endIntegral = 0.5 * t * t + endCorner;

            grid.Weight[slot] = width * nodeScale[node];
            grid.StartIntegral[slot] = startIntegral;
            grid.EndIntegral[slot] = endIntegral;
            grid.ShareRate[slot] = startIntegral / startBasisTurn -
                endIntegral / endBasisTurn;
            count++;
        }
    }

    ASSERT(count > 0U);

    grid.Count = count;
    grid.HasRates = withRates;

    if (!withRates) {
        return;
    }

    // Second pass for the basis curvatures, with the same nodes.
    U32 slot = 0U;

    for (U32 index = 1U; index < edgeCount; index++) {
        const double begin = edges[index - 1U];
        const double end = edges[index];
        const double width = end - begin;

        if (width <= 0.0) {
            continue;
        }

        for (U32 node = 0U; node < nodesPerPanel; node++) {
            ASSERT(slot < count);

            const double t = begin + width * nodeX[node];
            const double startCorner = S9UnitCornerValueWithScale(
                startLambda, startEta, startScale, t);
            const double endCorner = S9UnitCornerValueWithScale(endLambda,
                endEta, endScale, 1.0 - t);

            grid.StartRate[slot] = (1.0 - t) + startCorner;
            grid.EndRate[slot] = t + endCorner;
            slot++;
        }
    }

    ASSERT(slot == count);
}


/*
 * End point of the curve in units of the curve length, and the derivative
 * of that end point with respect to the start share.
 */
struct S9EndPoint final {
    FloatPoint Point;
    FloatPoint Slope;
};


static S9EndPoint S9MeasureEndPoint(
    const S9Grid &grid, S9Curve &shape,
    const double share) noexcept
{
    ASSERT(isfinite(share));
    ASSERT(grid.Count > 0U);

    S9SetLines(share, shape);

    const double startAngle = shape.StartAngle;
    const double startLine = shape.StartLine;
    const double endLine = shape.EndLine;
    double horizontal = 0.0;
    double vertical = 0.0;
    double horizontalSlope = 0.0;
    double verticalSlope = 0.0;

    for (U32 index = 0U; index < grid.Count; index++) {
        const double startIntegral = grid.StartIntegral[index];
        const double endIntegral = grid.EndIntegral[index];
        const double angle = startAngle + startLine * startIntegral +
            endLine * endIntegral;

        // Tangent angle changes with the share at this rate.
        const double rate = grid.ShareRate[index];
        const double weight = grid.Weight[index];
        const double cosine = cos(angle);
        const double sine = sin(angle);

        horizontal += weight * cosine;
        vertical += weight * sine;
        horizontalSlope -= weight * sine * rate;
        verticalSlope += weight * cosine * rate;
    }

    S9EndPoint result;

    result.Point = FloatPoint(horizontal, vertical);
    result.Slope = FloatPoint(horizontalSlope, verticalSlope);

    return result;
}


static double S9MeasureTurning(const S9Grid &grid,
    S9Curve &shape, const double share) noexcept
{
    ASSERT(isfinite(share));
    ASSERT(grid.Count > 0U);
    ASSERT(grid.HasRates);

    S9SetLines(share, shape);

    const double startLine = shape.StartLine;
    const double endLine = shape.EndLine;
    double turning = 0.0;

    for (U32 index = 0U; index < grid.Count; index++) {
        const double turnRate = startLine * grid.StartRate[index] +
            endLine * grid.EndRate[index];

        turning += grid.Weight[index] * Abs(turnRate);
    }

    ASSERT(turning >= 0.0);

    return turning;
}


/**
 * Sets the length of a shape from the end point of its solved share,
 * whose distance from the origin must be greater than
 * S9SmallestChordLength.
 */
static void S9CompleteShape(const FloatPoint endPoint,
    S9Curve &shape) noexcept
{
    const double distance = Hypot(endPoint.X, endPoint.Y);

    ASSERT(distance > S9SmallestChordLength);

    shape.Length = 1.0 / distance;

    ASSERT(shape.Length > 0.0);
    ASSERT(isfinite(shape.Length));
}


/**
 * Tells whether two finite vertical end point parts bracket a root.
 */
static inline constexpr bool S9PartsBracketRoot(const double first,
    const double second) noexcept
{
    ASSERT(first == first);
    ASSERT(second == second);

    const bool firstIsNegative = first < 0.0;
    const bool secondIsNegative = second < 0.0;

    return firstIsNegative != secondIsNegative;
}


/**
 * Tells whether a measured end point is a root.
 */
static bool S9EndPointIsRoot(const FloatPoint endPoint) noexcept
{
    const bool isOnChordLine = isfinite(endPoint.Y) and
        Abs(endPoint.Y) <= S9RootTolerance;

    return isOnChordLine and endPoint.X > 0.0;
}


/*
 * The best root of a scan so far. The best root is the one with the least
 * total turning.
 */
struct S9RootChoice final {
    double Share = 0.0;
    double Distance = 0.0;
    double Turning = 0.0;
    FloatPoint EndPoint;
    bool Found = false;
};


static void S9OfferRoot(const S9Grid &fine,
    S9Curve &shape, const double share, const FloatPoint endPoint,
    S9RootChoice &choice) noexcept
{
    ASSERT(fine.Count > 0U);
    ASSERT(isfinite(share));

    if (!(endPoint.X > 0.0)) {
        return;
    }

    const double turning = S9MeasureTurning(fine, shape, share);

    if (!choice.Found or turning < choice.Turning) {
        choice.Share = share;
        choice.Distance = Hypot(endPoint.X, endPoint.Y);
        choice.Turning = turning;
        choice.EndPoint = endPoint;
        choice.Found = true;
    }
}


static bool S9TryRefineRoot(const S9Grid &grid,
    S9Curve &shape, const double low, const double high,
    const double lowPart, const double highPart, double &root,
    FloatPoint &endPoint) noexcept
{
    ASSERT(grid.Count > 0U);
    ASSERT(isfinite(low));
    ASSERT(isfinite(high));
    ASSERT(S9PartsBracketRoot(lowPart, highPart));

    if (Abs(lowPart) <= S9RootTolerance) {
        const S9EndPoint measured = S9MeasureEndPoint(grid,
            shape, low);

        root = low;
        endPoint = measured.Point;

        return true;
    }

    if (Abs(highPart) <= S9RootTolerance) {
        const S9EndPoint measured = S9MeasureEndPoint(grid,
            shape, high);

        root = high;
        endPoint = measured.Point;

        return true;
    }

    double lowShare = low;
    double highShare = high;
    double lowValue = lowPart;
    double highValue = highPart;
    double guess = lowShare + 0.5 * (highShare - lowShare);
    const double valueSpan = lowValue - highValue;

    if (valueSpan != 0.0) {
        const double fraction = lowValue / valueSpan;

        if (fraction > 0.02 and fraction < 0.98) {
            guess = lowShare + (highShare - lowShare) * fraction;
        }
    }

    for (U32 index = 0U; index < S9RefineCount; index++) {
        if (guess == lowShare or guess == highShare) {
            break;
        }

        const S9EndPoint measured = S9MeasureEndPoint(grid,
            shape, guess);
        const double value = measured.Point.Y;

        if (!isfinite(value)) {
            break;
        }

        if (Abs(value) <= S9RootTolerance) {
            root = guess;
            endPoint = measured.Point;

            return true;
        }

        const bool sameSideAsLow = (value < 0.0) == (lowValue < 0.0);

        if (sameSideAsLow) {
            lowShare = guess;
            lowValue = value;
        } else {
            highShare = guess;
            highValue = value;
        }

        // Newton step from this measurement, kept inside the bracket.
        const double slope = measured.Slope.Y;
        double next = lowShare + 0.5 * (highShare - lowShare);

        if (slope != 0.0 and isfinite(slope)) {
            const double newton = guess - value / slope;

            if (newton > lowShare and newton < highShare) {
                next = newton;
            }
        }

        guess = next;
    }

    return false;
}


static void S9RefineBracket(const S9Grid &fine,
    S9Curve &shape, const double low, const double high,
    const FloatPoint lowPoint, const FloatPoint highPoint,
    S9RootChoice &choice) noexcept
{
    ASSERT(fine.Count > 0U);
    ASSERT(isfinite(low));
    ASSERT(isfinite(high));

    const bool isFinite = isfinite(lowPoint.Y) and isfinite(highPoint.Y);

    if (!isFinite) {
        return;
    }

    const bool brackets = S9PartsBracketRoot(lowPoint.Y, highPoint.Y);

    if (!brackets) {
        return;
    }

    double root = 0.0;
    FloatPoint endPoint;
    const bool refined = S9TryRefineRoot(fine, shape, low, high,
        lowPoint.Y, highPoint.Y, root, endPoint);

    if (refined) {
        S9OfferRoot(fine, shape, root, endPoint, choice);
    }
}


static void S9ExamineBracket(const S9Grid &fine,
    S9Curve &shape, const double low, const double high,
    const double step, S9RootChoice &choice) noexcept
{
    ASSERT(fine.Count > 0U);
    ASSERT(isfinite(low));
    ASSERT(high > low);
    ASSERT(step > 0.0);

    const S9EndPoint fineLow = S9MeasureEndPoint(fine, shape,
        low);
    const S9EndPoint fineHigh = S9MeasureEndPoint(fine, shape,
        high);
    const bool lowIsRoot = S9EndPointIsRoot(fineLow.Point);
    const bool highIsRoot = S9EndPointIsRoot(fineHigh.Point);

    if (lowIsRoot) {
        S9OfferRoot(fine, shape, low, fineLow.Point, choice);
    }

    if (highIsRoot) {
        S9OfferRoot(fine, shape, high, fineHigh.Point, choice);
    }

    if (lowIsRoot or highIsRoot) {
        return;
    }

    const bool isFinite = isfinite(fineLow.Point.Y) and
        isfinite(fineHigh.Point.Y);

    if (!isFinite) {
        return;
    }

    const bool brackets = S9PartsBracketRoot(fineLow.Point.Y,
        fineHigh.Point.Y);

    if (brackets) {
        S9RefineBracket(fine, shape, low, high, fineLow.Point,
            fineHigh.Point, choice);

        return;
    }

    // The root is near one end. It is in the interval on the other side
    // of that end, which the coarse grid did not see as a bracket.
    const double lowerLow = low - step;
    const S9EndPoint fineLowerLow = S9MeasureEndPoint(fine,
        shape, lowerLow);
    const bool lowerBrackets = isfinite(fineLowerLow.Point.Y) and
        S9PartsBracketRoot(fineLowerLow.Point.Y, fineHigh.Point.Y);

    if (lowerBrackets) {
        S9RefineBracket(fine, shape, lowerLow, high,
            fineLowerLow.Point, fineHigh.Point, choice);

        return;
    }

    const double upperHigh = high + step;
    const S9EndPoint fineUpperHigh = S9MeasureEndPoint(fine,
        shape, upperHigh);

    S9RefineBracket(fine, shape, low, upperHigh, fineLow.Point,
        fineUpperHigh.Point, choice);
}


static void S9SampleShareScan(const S9Grid &grid,
    S9Curve &shape,
    FloatPoint (&points)[S9ScanSampleCount]) noexcept
{
    ASSERT(grid.Count > 0U);
    ASSERT(grid.Count <= S9GreatestCoarseNodeCount);
    ASSERT(isfinite(shape.StartAngle));
    ASSERT(isfinite(shape.Turn));

    const double startBasisTurn = shape.StartBasisTurn();
    const double endBasisTurn = shape.EndBasisTurn();
    const double firstShare = 0.5 * shape.Turn - S9ScanReach;
    const double step = 2.0 * S9ScanReach /
        double(S9ScanSampleCount - 1U);

    // Every entry below grid. Count is written before it is read.
    double cosines[S9GreatestCoarseNodeCount];
    double sines[S9GreatestCoarseNodeCount];
    double stepCosines[S9GreatestCoarseNodeCount];
    double stepSines[S9GreatestCoarseNodeCount];
    double weightMagnitude = 0.0;
    double greatestIntegralMagnitude = 0.0;

    // The angle of a node is linear in the share, therefore a share step adds a
    // fixed delta per node. Angle addition identities update the sin
    // and cos of each node with the sin and cos of its
    // delta, instead of calling sin() and cos() for every share.
    for (U32 index = 0U; index < grid.Count; index++) {
        const double startIntegral = grid.StartIntegral[index] /
            startBasisTurn;
        const double endIntegral = grid.EndIntegral[index] / endBasisTurn;
        const double rate = startIntegral - endIntegral;
        const double increment = step * rate;
        const double absoluteStartIntegral = Abs(startIntegral);
        const double absoluteEndIntegral = Abs(endIntegral);
        const double integralMagnitude = absoluteStartIntegral +
            absoluteEndIntegral;
        const double absoluteWeight = Abs(grid.Weight[index]);

        stepCosines[index] = cos(increment);
        stepSines[index] = sin(increment);
        weightMagnitude += absoluteWeight;
        greatestIntegralMagnitude = Max(greatestIntegralMagnitude,
            integralMagnitude);
    }

    const double lastShare = firstShare +
        step * double(S9ScanSampleCount - 1U);
    const double absoluteFirstShare = Abs(firstShare);
    const double absoluteLastShare = Abs(lastShare);
    const double greatestShare = Max(absoluteFirstShare, absoluteLastShare);
    const double absoluteStartAngle = Abs(shape.StartAngle);
    const double absoluteTurn = Abs(shape.Turn);
    const double phaseScale = 1.0 + absoluteStartAngle +
        (absoluteTurn + 2.0 * greatestShare) * greatestIntegralMagnitude;
    const double weightScale = Max(1.0, weightMagnitude);
    const double signAllowance = 64.0 * DBL_EPSILON *
        double(grid.Count + S9ScanRestartInterval + 8U) *
        phaseScale * weightScale;

    for (U32 sample = 0U; sample < S9ScanSampleCount; sample++) {
        const double share = firstShare + step * double(sample);
        double horizontal = 0.0;
        double vertical = 0.0;

        if (sample % S9ScanRestartInterval == 0U) {
            S9SetLines(share, shape);

            for (U32 index = 0U; index < grid.Count; index++) {
                const double angle = shape.StartAngle +
                    shape.StartLine * grid.StartIntegral[index] +
                    shape.EndLine * grid.EndIntegral[index];
                const double cosine = cos(angle);
                const double sine = sin(angle);
                const double weight = grid.Weight[index];

                cosines[index] = cosine;
                sines[index] = sine;
                horizontal += weight * cosine;
                vertical += weight * sine;
            }
        } else {
            for (U32 index = 0U; index < grid.Count; index++) {
                const double cosine = cosines[index] * stepCosines[index] -
                    sines[index] * stepSines[index];
                const double sine = sines[index] * stepCosines[index] +
                    cosines[index] * stepSines[index];
                const double weight = grid.Weight[index];

                cosines[index] = cosine;
                sines[index] = sine;
                horizontal += weight * cosine;
                vertical += weight * sine;
            }
        }

        const double absoluteHorizontal = Abs(horizontal);
        const double absoluteVertical = Abs(vertical);

        if (!isfinite(signAllowance) or !isfinite(horizontal) or
            !isfinite(vertical) or absoluteHorizontal <= signAllowance or
            absoluteVertical <= signAllowance) {
            const S9EndPoint direct = S9MeasureEndPoint(grid,
                shape, share);

            points[sample] = direct.Point;
        } else {
            points[sample] = FloatPoint(horizontal, vertical);
        }
    }
}


static void S9ScanShare(const S9Grid &coarse,
    const S9Grid &fine, S9Curve &shape,
    S9RootChoice &choice) noexcept
{
    ASSERT(coarse.Count > 0U);
    ASSERT(fine.Count > 0U);
    ASSERT(isfinite(shape.Turn));

    const double center = 0.5 * shape.Turn;
    const double step = 2.0 * S9ScanReach /
        double(S9ScanSampleCount - 1U);
    double previousShare = center - S9ScanReach;
    FloatPoint points[S9ScanSampleCount];

    S9SampleShareScan(coarse, shape, points);

    FloatPoint previous = points[0U];

    for (U32 sample = 1U; sample < S9ScanSampleCount; sample++) {
        const double share = center - S9ScanReach + step * double(sample);
        const FloatPoint current = points[sample];
        const bool isFinite = isfinite(previous.Y) and isfinite(current.Y);
        const bool brackets = isFinite and
            S9PartsBracketRoot(previous.Y, current.Y);

        if (brackets) {
            // A bracket of the rotated samples is checked with direct
            // measurements.
            const S9EndPoint directPrevious = S9MeasureEndPoint(
                coarse, shape, previousShare);
            const S9EndPoint directCurrent = S9MeasureEndPoint(
                coarse, shape, share);
            const bool directIsFinite = isfinite(directPrevious.Point.Y) and
                isfinite(directCurrent.Point.Y);
            const bool facesForward = directPrevious.Point.X > 0.0 or
                directCurrent.Point.X > 0.0;
            const bool directBrackets = directIsFinite and
                S9PartsBracketRoot(directPrevious.Point.Y,
                    directCurrent.Point.Y);

            if (facesForward and directBrackets) {
                S9ExamineBracket(fine, shape, previousShare, share,
                    step, choice);
            }
        }

        previousShare = share;
        previous = current;
    }
}


/*
 * Limits of the predicted root.
 */
static constexpr const double S9PredictedAngleLimit = 2.1;
static constexpr const double S9PredictedReach = 3.0;
static constexpr const U32 S9PredictedStepCount = 16U;


static bool S9TryPredictedRoot(const S9Grid &fine,
    S9Curve &shape, double &share, FloatPoint &endPoint) noexcept
{
    ASSERT(fine.Count > 0U);
    ASSERT(isfinite(shape.Turn));

    const double startBasisTurn = shape.StartBasisTurn();
    const double endBasisTurn = shape.EndBasisTurn();
    double startSum = 0.0;
    double endSum = 0.0;

    for (U32 index = 0U; index < fine.Count; index++) {
        startSum += fine.Weight[index] * fine.StartIntegral[index];
        endSum += fine.Weight[index] * fine.EndIntegral[index];
    }

    const double rate = startSum / startBasisTurn - endSum / endBasisTurn;
    double predicted = 0.5 * shape.Turn;

    if (rate != 0.0) {
        predicted = -(shape.StartAngle + shape.Turn * endSum / endBasisTurn) /
            rate;
    }

    if (!isfinite(predicted)) {
        return false;
    }

    double current = predicted;

    for (U32 index = 0U; index < S9PredictedStepCount; index++) {
        const S9EndPoint measured = S9MeasureEndPoint(fine, shape,
            current);
        const double value = measured.Point.Y;

        if (!isfinite(value)) {
            return false;
        }

        if (Abs(value) <= S9RootTolerance) {
            const double distance = Hypot(measured.Point.X, measured.Point.Y);
            const bool isAhead = measured.Point.X > 0.0 and
                distance > S9SmallestChordLength;
            const double offset = Abs(current - 0.5 * shape.Turn);
            const bool isCentral = offset <= S9PredictedReach;

            if (isAhead and isCentral) {
                share = current;
                endPoint = measured.Point;

                return true;
            }

            return false;
        }

        const double slope = measured.Slope.Y;

        if (slope == 0.0 or !isfinite(slope)) {
            return false;
        }

        const double next = current - value / slope;
        const double travel = Abs(next - predicted);

        if (!isfinite(next) or travel > S9PredictedReach) {
            return false;
        }

        current = next;
    }

    return false;
}


bool S9Curve::TrySolveWithWeights(
    const S9CurveParameters &parameters) noexcept
{
    ASSERT(isfinite(parameters.Alpha));
    ASSERT(isfinite(parameters.Beta));
    ASSERT(parameters.StartWeight >= -1.0);
    ASSERT(parameters.StartWeight <= 1.0);
    ASSERT(parameters.EndWeight >= -1.0);
    ASSERT(parameters.EndWeight <= 1.0);
    ASSERT(parameters.StartEta >= 0.0);
    ASSERT(parameters.EndEta >= 0.0);

    S9Curve candidate;

    S9PrepareShape(parameters, candidate);

    S9Grid coarse;
    S9Grid fine;

    S9BuildGrid(candidate, S9CoarseBasePanelCount, 8U, false, coarse);
    S9BuildGrid(candidate, S9BasePanelCount, 16U, true, fine);

    static constexpr double fullTurn = 2.0 * Pi;
    const double principalTurn = candidate.Turn;
    const double rawTurn = S9RawTurn(candidate);

    const bool headingsAreModerate =
        Abs(candidate.StartAngle) <= S9PredictedAngleLimit and
        Abs(candidate.EndAngle) <= S9PredictedAngleLimit;

    if (headingsAreModerate) {
        candidate.Turn = rawTurn;

        double predictedShare = 0.0;
        FloatPoint predictedEnd;
        const bool isPredicted = S9TryPredictedRoot(fine, candidate,
            predictedShare, predictedEnd);

        if (isPredicted) {
            S9SetLines(predictedShare, candidate);
            S9CompleteShape(predictedEnd, candidate);
            *this = candidate;

            return true;
        }

        candidate.Turn = principalTurn;
    }

    const double turns[4] = {
        rawTurn,
        principalTurn,
        principalTurn - fullTurn,
        principalTurn + fullTurn
    };

    S9RootChoice choice;
    double bestTurn = rawTurn;

    for (U32 index = 0U; index < 4U and !choice.Found; index++) {
        const double turn = turns[index];
        const double difference = Abs(turn - rawTurn);

        if (index > 0U and difference <= S9SameTurnTolerance) {
            continue;
        }

        candidate.Turn = turn;
        S9ScanShare(coarse, fine, candidate, choice);
        bestTurn = turn;
    }

    if (!choice.Found or choice.Distance <= S9SmallestChordLength) {
        return false;
    }

    candidate.Turn = bestTurn;
    S9SetLines(choice.Share, candidate);
    S9CompleteShape(choice.EndPoint, candidate);
    *this = candidate;

    return true;
}


bool S9Curve::TrySolveWithWeightsNear(const S9CurveParameters &parameters,
    const S9Curve &nearby) noexcept
{
    ASSERT(isfinite(parameters.Alpha));
    ASSERT(isfinite(parameters.Beta));
    ASSERT(parameters.StartWeight >= -1.0);
    ASSERT(parameters.StartWeight <= 1.0);
    ASSERT(parameters.EndWeight >= -1.0);
    ASSERT(parameters.EndWeight <= 1.0);
    ASSERT(parameters.StartEta >= 0.0);
    ASSERT(parameters.EndEta >= 0.0);
    ASSERT(nearby.Length > 0.0);

    S9Curve candidate;

    S9PrepareShape(parameters, candidate);

    const double startChange = candidate.StartAngle - nearby.StartAngle;
    const double endChange = candidate.EndAngle - nearby.EndAngle;
    const double rawTurn = S9RawTurn(candidate);
    const double principalDifference = Abs(nearby.Turn - candidate.Turn);
    const double rawDifference = Abs(nearby.Turn - rawTurn);
    double branchTurn = candidate.Turn;
    double branchDifference = principalDifference;

    if (rawDifference < principalDifference) {
        branchTurn = rawTurn;
        branchDifference = rawDifference;
    }

    const bool anglesAreNear =
        Abs(startChange) <= S9GreatestNearAngleChange and
        Abs(endChange) <= S9GreatestNearAngleChange and
        branchDifference < Pi;
    const double nearbyShare = S9StartShare(nearby);
    bool isSolved = false;

    if (anglesAreNear and isfinite(nearbyShare)) {
        S9Grid fine;

        S9BuildGrid(candidate, S9BasePanelCount, 16U, false, fine);
        candidate.Turn = branchTurn;

        double share = nearbyShare;

        for (U32 index = 0U; index < S9GreatestNearStepCount; index++) {
            const S9EndPoint measured = S9MeasureEndPoint(fine,
                candidate, share);
            const double value = measured.Point.Y;

            if (!isfinite(value)) {
                break;
            }

            if (Abs(value) <= S9RootTolerance) {
                const double distance = Hypot(measured.Point.X,
                    measured.Point.Y);

                if (measured.Point.X > 0.0 and
                    distance > S9SmallestChordLength) {
                    S9CompleteShape(measured.Point, candidate);
                    isSolved = true;
                }

                break;
            }

            const double slope = measured.Slope.Y;

            if (slope == 0.0 or !isfinite(slope)) {
                break;
            }

            const double next = share - value / slope;
            const double distance = Abs(next - nearbyShare);

            if (!isfinite(next) or distance > S9NearReach) {
                break;
            }

            share = next;
        }
    }

    if (isSolved) {
        *this = candidate;

        return true;
    }

    return TrySolveWithWeights(parameters);
}


struct S9Approximation final {
    const S9Curve *Shape = nullptr;
    CurveArray<CubicCurve> *Curves = nullptr;
    double MaximumError = 0.0;
    U32 RemainingMeasurements = 0U;
    U32 RemainingCubics = 0U;
    bool Failed = false;
};


static FloatPoint S9MeasurePanel(const S9Curve &shape,
    const double begin, const double width) noexcept
{
    ASSERT(begin >= 0.0);
    ASSERT(width > 0.0);
    ASSERT(begin + width <= 1.0 + DBL_EPSILON);

    double horizontal = 0.0;
    double vertical = 0.0;

    for (U32 index = 0U; index < 16U; index++) {
        const double t = Min(begin + width * GaussLegendreX_16[index], 1.0);
        const double angle = shape.TangentAngle(t);
        const double weight = width * GaussLegendreScale_16[index];

        horizontal += weight * cos(angle);
        vertical += weight * sin(angle);
    }

    return FloatPoint(horizontal * shape.Length, vertical * shape.Length);
}


static void S9AddCubic(S9Approximation &approximation,
    const FloatPoint startPoint, const FloatPoint endPoint,
    const double startAngle, const double endAngle,
    const double arcLength) noexcept
{
    ASSERT(approximation.Curves != nullptr);
    ASSERT(approximation.RemainingCubics > 0U);
    ASSERT(arcLength > 0.0);

    const double handle = arcLength / 3.0;
    const FloatPoint startDirection(cos(startAngle), sin(startAngle));
    const FloatPoint endDirection(cos(endAngle), sin(endAngle));
    const FloatPoint firstControl = startPoint + startDirection * handle;
    const FloatPoint secondControl = endPoint - endDirection * handle;

    approximation.Curves->CurveTo(firstControl, secondControl, endPoint);
    approximation.RemainingCubics--;
}


static void S9ApproximateRange(S9Approximation &approximation,
    const double begin, const double end, const U32 depth,
    const FloatPoint startPoint, FloatPoint &endPoint) noexcept
{
    ASSERT(approximation.Shape != nullptr);
    ASSERT(begin >= 0.0);
    ASSERT(end > begin);
    ASSERT(end <= 1.0 + DBL_EPSILON);
    ASSERT(depth <= S9GreatestDepth);

    const S9Curve &shape = *approximation.Shape;
    const double width = end - begin;
    const double middle = begin + 0.5 * width;
    const double quarterWidth = 0.25 * width;
    const double startAngle = shape.TangentAngle(begin);

    FloatPoint piecePoint = startPoint;

    FloatPoint samplePoints[3] = {
        startPoint,
        startPoint,
        startPoint
    };

    double sampleAngles[3] = {
        startAngle,
        startAngle,
        startAngle
    };

    double variation = 0.0;
    double pieceAngle = startAngle;

    for (U32 index = 0U; index < 4U; index++) {
        const double pieceBegin = begin + quarterWidth * double(index);
        const double pieceEnd = Min(pieceBegin + quarterWidth, 1.0);
        const FloatPoint displacement = S9MeasurePanel(shape,
            pieceBegin, quarterWidth);
        const double nextAngle = shape.TangentAngle(pieceEnd);

        piecePoint += displacement;
        variation += Abs(nextAngle - pieceAngle);
        pieceAngle = nextAngle;

        if (index < 3U) {
            samplePoints[index] = piecePoint;
            sampleAngles[index] = pieceAngle;
        }
    }

    const double endAngle = pieceAngle;
    const double splitTurn = endAngle - startAngle;
    const FloatPoint splitDisplacement = piecePoint - startPoint;
    const double turnMagnitude = Max(variation, Abs(splitTurn));
    bool needsSplit = false;

    if (turnMagnitude > S9GreatestPanelTurn) {
        needsSplit = true;
    }

    if (!needsSplit) {
        const FloatPoint wholeDisplacement = S9MeasurePanel(shape,
            begin, width);
        const FloatPoint difference = wholeDisplacement - splitDisplacement;
        const double error = Hypot(difference.X, difference.Y);
        const double share = S9IntegrationAllowance * width;
        const double scale = width + Abs(startPoint.X) + Abs(startPoint.Y) +
            Abs(piecePoint.X) + Abs(piecePoint.Y);
        const double floor = 32.0 * DBL_EPSILON * scale;
        const double panelAllowance = Max(share, floor);

        if (error > panelAllowance) {
            needsSplit = true;
        }
    }

    // The cubic curve is compared at three interior points, because it can
    // pass through the middle point and still leave the curve.
    if (!needsSplit) {
        const double arcLength = width * shape.Length;
        const double handle = arcLength / 3.0;
        const FloatPoint startDirection(cos(startAngle), sin(startAngle));
        const FloatPoint endDirection(cos(endAngle), sin(endAngle));
        const FloatPoint firstControl = startPoint + startDirection * handle;
        const FloatPoint secondControl = piecePoint - endDirection * handle;
        const double parameters[3] = {0.25, 0.5, 0.75};

        for (U32 index = 0U; index < 3U; index++) {
            const double t = parameters[index];
            const double s = 1.0 - t;
            const double startShare = s * s * s;
            const double firstShare = 3.0 * s * s * t;
            const double secondShare = 3.0 * s * t * t;
            const double endShare = t * t * t;
            const FloatPoint onCubic = startShare * startPoint +
                firstShare * firstControl + secondShare * secondControl +
                endShare * piecePoint;
            const double sampleAngle = sampleAngles[index];
            const FloatPoint across(-sin(sampleAngle), cos(sampleAngle));
            const FloatPoint cubicError = onCubic - samplePoints[index];
            const double crossDistance = cubicError.X * across.X +
                cubicError.Y * across.Y;

            if (Abs(crossDistance) > approximation.MaximumError) {
                needsSplit = true;
            }
        }
    }

    static constexpr U32 spent = 5U;

    if (approximation.RemainingMeasurements > spent) {
        approximation.RemainingMeasurements -= spent;
    } else {
        approximation.RemainingMeasurements = 0U;
        approximation.Failed = true;
    }

    if (needsSplit and depth < S9GreatestDepth and !approximation.Failed) {
        FloatPoint innerPoint;

        S9ApproximateRange(approximation, begin, middle, depth + 1U,
            startPoint, innerPoint);
        S9ApproximateRange(approximation, middle, end, depth + 1U,
            innerPoint, endPoint);

        return;
    }

    endPoint = piecePoint;

    if (approximation.Failed) {
        return;
    }

    if (approximation.RemainingCubics == 0U) {
        approximation.Failed = true;

        return;
    }

    const double arcLength = width * shape.Length;

    S9AddCubic(approximation, startPoint, endPoint, startAngle, endAngle,
        arcLength);
}


static bool S9TryApproximateSolvedShape(const S9Curve &solved,
    const double maximumError, CurveArray<CubicCurve> &curves) noexcept
{
    ASSERT(maximumError > 0.0);
    ASSERT(solved.Length > 0.0);

    const FloatPoint startPoint;
    S9Approximation approximation;

    approximation.Shape = &solved;
    approximation.Curves = &curves;
    approximation.MaximumError = maximumError;
    approximation.RemainingMeasurements = S9GreatestPanelCount;
    approximation.RemainingCubics = S9GreatestCubicCount;

    curves.Restart(startPoint);

    FloatPoint finalPoint;

    S9ApproximateRange(approximation, 0.0, 1.0, 0U, startPoint, finalPoint);

    const U32 curveCount = curves.GetCurveCount();
    if (approximation.Failed or curveCount == 0U) {
        curves.Erase();
        curves.Restart(startPoint);

        return false;
    }

    ASSERT(curveCount > 0U);

    const U32 pointCount = curves.GetPointCount();
    const U32 lastIndex = pointCount - 1U;
    const FloatPoint lastPoint(1.0, 0.0);
    const FloatPoint firstControl = curves.GetPointAtIndex(1U);
    const FloatPoint secondControl = curves.GetPointAtIndex(lastIndex - 1U);
    const FloatPoint firstEdge = firstControl - startPoint;
    const FloatPoint secondEdge = secondControl - lastPoint;
    const double firstHandle = Hypot(firstEdge.X, firstEdge.Y);
    const double secondHandle = Hypot(secondEdge.X, secondEdge.Y);
    const FloatPoint startDirection(cos(solved.StartAngle),
        sin(solved.StartAngle));
    const FloatPoint endDirection(cos(solved.EndAngle), sin(solved.EndAngle));

    const FloatPoint correctedFirstControl = startPoint +
        startDirection * firstHandle;
    const FloatPoint correctedSecondControl = lastPoint -
        endDirection * secondHandle;
    curves.SetPointAtIndex(1U, correctedFirstControl);
    curves.SetPointAtIndex(lastIndex - 1U, correctedSecondControl);
    curves.SetPointAtIndex(lastIndex, lastPoint);

    return true;
}


bool S9Curve::TryApproximateWithCubics(const double maximumError,
    CurveArray<CubicCurve> &curves) const noexcept
{
    ASSERT(maximumError > 0.0);

    const FloatPoint startPoint;
    const bool isSolved = Length > 0.0;

    ASSERT(isSolved);
    if (!isSolved) {
        curves.Erase();
        curves.Restart(startPoint);

        return false;
    }

    return S9TryApproximateSolvedShape(*this, maximumError, curves);
}
