
#include "SplineKappa.h"
#include "TridiagonalSolver.h"


/**
 * Maximum number of solver iterations. Should not be less than 20. Solver
 * stops earlier if the last iteration did not result it large enough change
 * from previous iteration.
 */
static constexpr U32 MaximumSolverIterations = 60;


/**
 * Solver stops once none of the points moved more than certain distance since
 * last iteration. This constant defines squared distance for at least one
 * point to move for solver to continue.
 *
 * Lower value means higher output precision. Should not be larger than 1e-16.
 */
static constexpr double SolverStopEpsilon = 1e-8;


/**
 * Find a single quadratic root within interval [0..1].
 */
static double FindQuadraticRoot(const double a, const double b, const double c)
{
    const double delta = b * b - 4.0 * a * c;

    if (delta < 0.0) {
        return 0.0;
    }

    if (delta > 0.0) {
        const double d = Sqrt(delta);
        const double q = -0.5 * (b + (b < 0.0 ? -d : d));
        const double rv0 = q / a;
        const double rv1 = c / q;

        if (FuzzyIsEqual(rv0, rv1)) {
            return rv0;
        }

        if (rv0 >= 0.0 and rv0 <= 1.0) {
            return rv0;
        }

        return rv1;
    }

    if (a != 0) {
        return -0.5 * b / a;
    }

    // Normally, should not happen.

    return 0.5;
}


static double KappaCubeRoot(const double x) {
    return cbrt(x);
}


/**
 * Find a single cubic root within interval [0..1].
 */
static double FindCubicRoot(const double coe0, const double coe1, const double coe2,
    const double coe3)
{
    if (FuzzyIsZero(coe0)) {
        return FindQuadraticRoot(coe1, coe2, coe3);
    }

    const double inva = 1.0 / coe0;

    const double a = coe1 * inva;
    const double b = coe2 * inva;
    const double c = coe3 * inva;

    const double Q = (a * a - b * 3.0) / 9.0;
    const double R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;

    const double R2 = R * R;
    const double Q3 = Q * Q * Q;
    const double R2subQ3 = R2 - Q3;
    const double adiv3 = a / 3.0;

    if (R2subQ3 < 0.0) {
        // If Q and R are real (always true when a, b, c are real) and R2 < Q3,
        // then the cubic equation has three real roots.
        const double theta = Acos(Clamp(R / Sqrt(Q3), -1.0, 1.0));
        const double negative2RootQ = -2.0 * Sqrt(Q);

        const double x0 = negative2RootQ * Cos(theta / 3.0) - adiv3;
        const double x1 = negative2RootQ * Cos((theta + 2.0 * M_PI) / 3.0) - adiv3;
        const double x2 = negative2RootQ * Cos((theta - 2.0 * M_PI) / 3.0) - adiv3;

        if (x0 >= 0.0 and x0 <= 1.0) {
            return x0;
        }

        if (x1 >= 0.0 and x1 <= 1.0) {
            return x1;
        }

        if (x2 >= 0.0 and x2 <= 1.0) {
            return x2;
        }
    }

    double A = Abs(R) + Sqrt(R2subQ3);

    A = KappaCubeRoot(A);

    if (R > 0.0) {
        A = -A;
    }

    if (A != 0.0) {
        A += Q / A;
    }

    return A - adiv3;
}


static double MaxParam(const FloatPoint c0, const FloatPoint cp,
    const FloatPoint c2)
{
    const double a = (c2 - c0).Dot(c2 - c0);
    const double b = 3.0 * (c2 - c0).Dot(c0 - cp);
    const double c = (3.0 * c0 - 2.0 * cp - c2).Dot(c0 - cp);
    const double d = -(c0 - cp).Dot(c0 - cp);

    return FindCubicRoot(a, b, c, d);
}


static double CalculateTriangleArea(const FloatPoint A, const FloatPoint B,
    const FloatPoint C)
{
    return Abs((A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y))) / 2.0;
}


static double Lambda(const FloatPoint c0, const FloatPoint c1,
    const FloatPoint c3, const FloatPoint c4)
{
    const double area_c0c1c3 = CalculateTriangleArea(c0, c1, c3);
    const double sqrt_area_c0c1c3 = sqrt(area_c0c1c3);
    const double area_c1c3c4 = CalculateTriangleArea(c1, c3, c4);
    const double sqrt_area_c1c3c4 = sqrt(area_c1c3c4);

    return sqrt_area_c0c1c3 / (sqrt_area_c0c1c3 + sqrt_area_c1c3c4);
}


/**
 * Writes all solved quadratics if all their coordinates are finite.
 * Leaves the destination empty when a coordinate is not finite.
 *
 * @param quadratics Non-null array of count solved quadratics. The array is
 * read during this call. The function does not keep a reference or take
 * ownership.
 *
 * @param count Number of quadratics to inspect, at least one.
 *
 * @param curves Empty destination. Its storage is reused. Points are copied
 * only when all quadratics pass the finite-coordinate check.
 */
static void KappaWriteFiniteCurves(const QuadraticCurve *const quadratics,
    const U32 count, CurveArray<QuadraticCurve> &curves) noexcept
{
    ASSERT(quadratics != nullptr);
    ASSERT(count > 0U);
    ASSERT(curves.GetPointCount() == 0U);

    for (U32 index = 0U; index < count; index++) {
        const QuadraticCurve &curve = quadratics[index];
        const bool startIsFinite = NumberIsFinite(curve.P0.X) and
            NumberIsFinite(curve.P0.Y);
        const bool controlIsFinite = NumberIsFinite(curve.P1.X) and
            NumberIsFinite(curve.P1.Y);
        const bool endIsFinite = NumberIsFinite(curve.P2.X) and
            NumberIsFinite(curve.P2.Y);

        if (not startIsFinite or not controlIsFinite or not endIsFinite) {
            return;
        }
    }

    curves.Restart(quadratics[0].P0);

    for (U32 index = 0U; index < count; index++) {
        curves.CurveTo(quadratics[index].P1, quadratics[index].P2);
    }
}


void SolveKappaCyclic(const FloatPoint *knots, const U32 count,
    CurveArray<QuadraticCurve> &curves)
{
    curves.Erase();

    if (count < 3U or knots == nullptr) {
        return;
    }

    const U32 n = count;

    QuadraticCurve *quads = AllocElementsInit<QuadraticCurve>(n);

    quads[0].P0 = (knots[0] + knots[n - 1]) / 2.0;
    quads[0].P1 = knots[0];
    quads[0].P2 = (knots[0] + knots[1]) / 2.0;

    for (U32 i = 1; i < n - 1; i++) {
        const FloatPoint mk0 = knots[i] + knots[i - 1];
        const FloatPoint mk1 = knots[i] + knots[i + 1];

        quads[i].P0 = mk0 / 2.0;
        quads[i].P1 = knots[i];
        quads[i].P2 = mk1 / 2.0;
    }

    quads[n - 1].P0 = (knots[n - 1] + knots[n - 2]) / 2.0;
    quads[n - 1].P1 = knots[n - 1];
    quads[n - 1].P2 = (knots[n - 1] + knots[0]) / 2.0;

    double *m = AllocElements<double>(n);
    double *l = AllocElements<double>(n);

    double *A = AllocElements<double>(n);
    double *B = AllocElements<double>(n);
    double *C = AllocElements<double>(n);

    double *X = AllocElements<double>(n);
    double *Y = AllocElements<double>(n);

    double *bx = AllocElements<double>(n);
    double *by = AllocElements<double>(n);

    for (U32 i = 0; i < n; i++) {
        X[i] = knots[i].X;
        Y[i] = knots[i].Y;
    }

    double maxError = 0.0;

    U32 it = MaximumSolverIterations;

    do {
        m[0] = MaxParam(quads[0].P0, knots[0], quads[0].P2);
        l[0] = Lambda(quads[0].P0, quads[0].P1, quads[1].P1, quads[1].P2);

        for (U32 i = 1; i < n - 1; i++) {
            m[i] = MaxParam(quads[i].P0, knots[i], quads[i].P2);
            l[i] = Lambda(quads[i].P0, quads[i].P1, quads[i + 1].P1, quads[i + 1].P2);
        }

        m[n - 1] = MaxParam(quads[n - 1].P0, knots[n - 1], quads[n - 1].P2);
        l[n - 1] = Lambda(quads[n - 1].P0, quads[n - 1].P1, quads[0].P1, quads[0].P2);

        // Init matrix diagonals.
        A[0] = (1.0 - l[n - 1]) * (1.0 - m[0]) * (1.0 - m[0]);
        B[0] = l[n - 1] * (1.0 - m[0]) * (1.0 - m[0]) + 2.0 * (1.0 - m[0]) * m[0] + (1.0 - l[0]) * m[0] * m[0];
        C[0] = l[0] * m[0] * m[0];

        for (U32 i = 1; i < n - 1; i++) {
            A[i] = (1.0 - l[i - 1]) * (1.0 - m[i]) * (1.0 - m[i]);
            B[i] = l[i - 1] * (1.0 - m[i]) * (1.0 - m[i]) + 2.0 * (1.0 - m[i]) * m[i] + (1.0 - l[i]) * m[i] * m[i];
            C[i] = l[i] * m[i] * m[i];
        }

        A[n - 1] = (1.0 - l[n - 2]) * (1.0 - m[n - 1]) * (1.0 - m[n - 1]);
        B[n - 1] = l[n - 2] * (1.0 - m[n - 1]) * (1.0 - m[n - 1]) + 2.0 * (1.0 - m[n - 1]) * m[n - 1] + (1.0 - l[n - 1]) * m[n - 1] * m[n - 1];
        C[n - 1] = l[n - 1] * m[n - 1] * m[n - 1];

        SolveTridiagonalCyclic(n, &A[0], &B[0], &C[0], &X[0], &bx[0]);
        SolveTridiagonalCyclic(n, &A[0], &B[0], &C[0], &Y[0], &by[0]);

        // Update control points.
        for (U32 i = 0; i < n; i++) {
            quads[i].P1 = FloatPoint(bx[i], by[i]);
        }

        // Update end points.
        const FloatPoint p0 = (1 - l[0]) * quads[0].P1 + l[0] * quads[1].P1;

        maxError = quads[1].P0.DistanceToSquared(p0);

        quads[1].P0 = p0;
        quads[0].P2 = p0;

        for (U32 i = 1; i < n - 1; i++) {
            const FloatPoint p = (1 - l[i]) * quads[i].P1 + l[i] * quads[i + 1].P1;

            maxError = Max(maxError, quads[i + 1].P0.DistanceToSquared(p));

            quads[i + 1].P0 = p;
            quads[i].P2 = p;
        }

        const FloatPoint p2 = (1 - l[n - 1]) * quads[n - 1].P1 + l[n - 1] * quads[0].P1;

        maxError = Max(maxError, quads[0].P0.DistanceToSquared(p2));

        quads[0].P0 = p2;
        quads[n - 1].P2 = p2;

        it--;
    } while (maxError > SolverStopEpsilon and it > 0);

    KappaWriteFiniteCurves(quads, n, curves);

    Dealloc(by);
    Dealloc(bx);
    Dealloc(Y);
    Dealloc(X);
    Dealloc(C);
    Dealloc(B);
    Dealloc(A);
    Dealloc(l);
    Dealloc(m);
    Dealloc(quads);
}


void SolveKappaAcyclic(const FloatPoint *knots, const U32 count,
    CurveArray<QuadraticCurve> &curves)
{
    curves.Erase();

    if (count < 4U or knots == nullptr) {
        return;
    }

    const U32 n = count - 2;

    QuadraticCurve *quads = AllocElementsInit<QuadraticCurve>(n);

    for (U32 i = 1; i <= n; i++) {
        const FloatPoint mk0 = knots[i] + knots[i - 1];
        const FloatPoint mk1 = knots[i] + knots[i + 1];

        quads[i - 1].P0 = mk0 / 2.0;
        quads[i - 1].P1 = knots[i];
        quads[i - 1].P2 = mk1 / 2.0;
    }

    quads[0].P0 = knots[0];
    quads[n - 1].P2 = knots[n + 1];

    double *t = AllocElements<double>(n);
    double *ld = AllocElements<double>(n - 1);

    double *A = AllocElements<double>(n);
    double *B = AllocElements<double>(n);
    double *C = AllocElements<double>(n);

    double *X = AllocElements<double>(n);
    double *Y = AllocElements<double>(n);

    double *bx = AllocElements<double>(n);
    double *by = AllocElements<double>(n);

    for (U32 i = 0; i < n; i++) {
        X[i] = knots[i + 1].X;
        Y[i] = knots[i + 1].Y;
    }

    for (U32 iter = 0; iter < 50; iter++) {
        for (U32 i = 0; i < n; i++) {
            t[i] = MaxParam(quads[i].P0, knots[i + 1], quads[i].P2);
        }

        for (int i = 0; i < n - 1; i++) {
            ld[i] = Lambda(quads[i].P0, quads[i].P1, quads[i + 1].P1, quads[i + 1].P2);
        }

        A[0] = 0.0;
        B[0] = 2.0 * (1 - t[0]) * t[0] + (1.0 - ld[0]) * t[0] * t[0];
        C[0] = ld[0] * t[0] * t[0];

        for (U32 i = 1; i < n - 1; i++) {
            A[i] = (1.0 - ld[i - 1]) * (1.0 - t[i]) * (1.0 - t[i]);
            B[i] = ld[i - 1] * (1.0 - t[i]) * (1.0 - t[i]) + 2.0 * (1.0 - t[i]) * t[i] + (1.0 - ld[i]) * t[i] * t[i];
            C[i] = ld[i] * t[i] * t[i];
        }

        A[n - 1] = (1.0 - ld[n - 2]) * (1.0 - t[n - 1]) * (1.0 - t[n - 1]);
        B[n - 1] = ld[n - 2] * (1.0 - t[n - 1]) * (1.0 - t[n - 1]) + 2.0 * (1.0 - t[n - 1]) * t[n - 1];
        C[n - 1] = 0.0;

        const FloatPoint v0 = knots[1] - (1 - t[0]) * (1 - t[0]) * knots[0];
        const FloatPoint v1 = knots[n] - t[n - 1] * t[n - 1] * knots[n + 1];

        X[0] = v0.X;
        Y[0] = v0.Y;
        X[n - 1] = v1.X;
        Y[n - 1] = v1.Y;

        SolveTridiagonal(n, &A[0], &B[0], &C[0], &X[0], &bx[0]);
        SolveTridiagonal(n, &A[0], &B[0], &C[0], &Y[0], &by[0]);

        for (U32 i = 0; i < n; i++) {
            quads[i].P1 = FloatPoint(bx[i], by[i]);
        }

        for (U32 i = 0; i < n - 1; i++) {
            quads[i + 1].P0 = quads[i].P2 = (1 - ld[i]) * quads[i].P1 + ld[i] * quads[i + 1].P1;
        }
    }

    KappaWriteFiniteCurves(quads, n, curves);

    Dealloc(by);
    Dealloc(bx);
    Dealloc(Y);
    Dealloc(X);
    Dealloc(C);
    Dealloc(B);
    Dealloc(A);
    Dealloc(ld);
    Dealloc(t);
    Dealloc(quads);
}
