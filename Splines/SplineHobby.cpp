
#include "DeduplicatePoints.h"
#include "SplineHobby.h"
#include "TridiagonalSolver.h"


static double Velocity(const double a, const double b)
{
    const double sa = Sin(a);
    const double sb = Sin(b);
    const double ca = Cos(a);
    const double cb = Cos(b);
    const double s5 = Sqrt(5.0);
    const double num = 4.0 + Sqrt(8.0) * (sa - sb / 16.0) * (sb - sa / 16.0) * (ca - cb);
    const double den = 2.0 + (s5 - 1.0) * ca + (3.0 - s5) * cb;

    return num / den;
}


static FloatPoint Rotate(const double x, const double y, const double sin, const double cos) {
    return FloatPoint(x * cos - y * sin, x * sin + y * cos);
}


static FloatPoint RotateAngle(const double x, const double y, const double alpha) {
    return Rotate(x, y, Sin(alpha), Cos(alpha));
}


static void DoSolveHobbyCyclic(const U32 count, const FloatPoint *knots,
    CurveArray<CubicCurve> &curves)
{
    ASSERT(count >= 3);
    ASSERT(knots != nullptr);

    const U32 n = count;

    double *D = AllocElements<double>(n);
    double *dx = AllocElements<double>(n);
    double *dy = AllocElements<double>(n);

    for (U32 i = 0; i < n; i++) {
        const U32 j = (i + 1) % n;

        dx[i] = knots[j].X - knots[i].X;
        dy[i] = knots[j].Y - knots[i].Y;
        D[i] = Sqrt(dx[i] * dx[i] + dy[i] * dy[i]);
    }

    double *gamma = AllocElements<double>(n);

    for (U32 i = 0; i < n; i++) {
        const U32 k = (i + n - 1) % n;
        const double sin = dy[k] / D[k];
        const double cos = dx[k] / D[k];
        const FloatPoint pp = Rotate(dx[i], dy[i], -sin, cos);

        gamma[i] = Atan2(pp.Y, pp.X);
    }

    double *a = AllocElements<double>(n);
    double *b = AllocElements<double>(n);
    double *c = AllocElements<double>(n);
    double *d = AllocElements<double>(n);

    for (U32 i = 0; i < n; i++) {
        const U32 j = (i + 1) % n;
        const U32 k = (i + n - 1) % n;

        a[i] = 1.0 / D[k];
        b[i] = (2.0 * D[k] + 2.0 * D[i]) / (D[k] * D[i]);
        c[i] = 1.0 / D[i];
        d[i] = -(2.0 * gamma[i] * D[i] + gamma[j] * D[k]) / (D[k] * D[i]);
    }

    double *alpha = AllocElements<double>(n);

    SolveTridiagonalCyclic(n, a, b, c, d, alpha);

    double *beta = AllocElements<double>(n);

    for (U32 i = 0; i < n; i++) {
        const U32 j = (i + 1) % n;

        beta[i] = -gamma[j] - alpha[j];
    }


    curves.Restart(knots[0]);

    for (U32 i = 0; i < n; i++) {
        const U32 j = (i + 1) % n;
        const double a = Velocity(alpha[i], beta[i]) * D[i] / 3;
        const double b = Velocity(beta[i], alpha[i]) * D[i] / 3;

        const FloatPoint xy0 = RotateAngle(dx[i], dy[i], alpha[i]).UnitVector();
        const FloatPoint xy1 = RotateAngle(dx[i], dy[i], -beta[i]).UnitVector();

        curves.CurveTo(knots[i] + (xy0 * a), knots[j] - (xy1 * b), knots[j]);
    }

    Dealloc(beta);
    Dealloc(alpha);
    Dealloc(d);
    Dealloc(c);
    Dealloc(b);
    Dealloc(a);
    Dealloc(gamma);
    Dealloc(dy);
    Dealloc(dx);
    Dealloc(D);
}


static void DoSolveHobbyAcyclic(const U32 count, const FloatPoint *knots,
    CurveArray<CubicCurve> &curves, const double omega)
{
    ASSERT(count >= 3);
    ASSERT(knots != nullptr);

    const U32 n = count - 1;

    double *D = new double[n];
    double *dx = new double[n];
    double *dy = new double[n];

    for (U32 i = 0; i < n; i++) {
        dx[i] = knots[i+1].X-knots[i].X;
        dy[i] = knots[i+1].Y-knots[i].Y;

        D[i] = Sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
    }

    double *gamma = new double[n+1];

    for (U32 i = 1; i < n; i++) {
        const double sin = dy[i-1] / D[i-1];
        const double cos = dx[i-1] / D[i-1];

        const FloatPoint pp = Rotate(dx[i], dy[i], -sin, cos);

        gamma[i] = Atan2(pp.Y, pp.X);
    }

    gamma[n] = 0;

    double *a = new double[n+1];
    double *b = new double[n+1];
    double *c = new double[n+1];
    double *d = new double[n+1];

    a[0] = 0;
    b[0] = 2 + omega;
    c[0] = 2 * omega + 1;
    d[0] = -c[0] * gamma[1];

    for (U32 i = 1; i < n; i++) {
        a[i] = 1 / D[i-1];
        b[i] = (2*D[i-1]+2*D[i])/(D[i-1]*D[i]);
        c[i] = 1 / D[i];
        d[i] = -(2*gamma[i]*D[i]+gamma[i+1]*D[i-1])/(D[i-1]*D[i]);
    }

    a[n] = 2 * omega + 1;
    b[n] = 2 + omega;
    c[n] = 0;
    d[n] = 0;

    double *alpha = new double[n+1];

    SolveTridiagonal(n + 1, a, b, c, d, alpha);

    double *beta = new double[n];

    for (U32 i = 0; i < n - 1; i++) {
        beta[i] = -gamma[i+1]-alpha[i+1];
    }

    beta[n-1] = -alpha[n];

    curves.Restart(knots[0]);

    for (U32 i = 0; i < n; i++) {
        const double a = Velocity(alpha[i], beta[i]) * D[i] / 3;
        const double b = Velocity(beta[i], alpha[i]) * D[i] / 3;

        const FloatPoint xy0 = RotateAngle(dx[i], dy[i], alpha[i]).UnitVector();
        const FloatPoint xy1 = RotateAngle(dx[i], dy[i], -beta[i]).UnitVector();

        curves.CurveTo(knots[i] + (xy0 * a), knots[i + 1] - (xy1 * b), knots[i + 1]);
    }

    delete [] beta;
    delete [] alpha;
    delete [] d;
    delete [] c;
    delete [] b;
    delete [] a;
    delete [] gamma;
    delete [] dy;
    delete [] dx;
    delete [] D;
}


void SolveHobbyCyclic(const U32 count, const FloatPoint *knots,
    CurveArray<CubicCurve> &curves)
{
    if (count > 0) {
        ASSERT(knots != nullptr);

        if (count == 1) {
            //l[0] = FloatPoint(0, 0);
            //r[0] = FloatPoint(0, 0);
        } else if (count == 2) {
            //l[0] = FloatPoint(0, 0);
            //r[0] = FloatPoint(0, 0);
            //l[1] = FloatPoint(0, 0);
            //r[1] = FloatPoint(0, 0);
        } else {
            DoSolveHobbyCyclic(count, knots, curves);
        }
    }
}


void SolveHobbyAcyclic(const U32 count, const FloatPoint *knots,
    CurveArray<CubicCurve> &curves, const double omega)
{
    if (count > 0) {
        ASSERT(knots != nullptr);

        FloatPoint *d = new FloatPoint[count];

        const U32 n = DeduplicatePoints(count, knots, d);

        if (n == 1) {
            curves.Restart(d[0]);
        } else if (n == 2) {
            curves.Restart(d[0]);

            const FloatPoint cp0 = InterpolateLinear(d[0], d[1], 1.0 / 3.0);
            const FloatPoint cp1 = InterpolateLinear(d[0], d[1], 2.0 / 3.0);

            curves.CurveTo(cp0, cp1, d[1]);
        } else if (n > 2) {
            DoSolveHobbyAcyclic(n, d, curves, omega);
        }

        delete [] d;
    }
}
