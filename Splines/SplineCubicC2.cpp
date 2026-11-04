
#include "SplineCubicC2.h"
#include "TridiagonalSolver.h"


static void DoSolveCubicC2Cyclic(const U32 count,
    const FloatPoint *knots, CurveArray<CubicCurve> &curves)
{
    ASSERT(count >= 3);
    ASSERT(knots != nullptr);

    const U32 n = count;

    double *ax = AllocElements<double>(n);
    double *bx = AllocElements<double>(n);
    double *cx = AllocElements<double>(n);
    double *dx = AllocElements<double>(n);
    double *ay = AllocElements<double>(n);
    double *by = AllocElements<double>(n);
    double *cy = AllocElements<double>(n);
    double *dy = AllocElements<double>(n);

    ax[0] = 1.0;
    bx[0] = 4.0;
    cx[0] = 1.0;
    dx[0] = 4.0 * knots[0].X + 2.0 * knots[1].X;

    ay[0] = 1.0;
    by[0] = 4.0;
    cy[0] = 1.0;
    dy[0] = 4.0 * knots[0].Y + 2.0 * knots[1].Y;

    for (U32 i = 1; i < n - 1; i++) {
        ax[i] = 1.0;
        bx[i] = 4.0;
        cx[i] = 1.0;
        dx[i] = 4.0 * knots[i].X + 2.0 * knots[i + 1].X;

        ay[i] = 1.0;
        by[i] = 4.0;
        cy[i] = 1.0;
        dy[i] = 4.0 * knots[i].Y + 2.0 * knots[i + 1].Y;
    }

    const U32 r = n - 1;

    ax[r] = 1.0;
    bx[r] = 4.0;
    cx[r] = 1.0;
    dx[r] = 4 * knots[r].X + 2.0 * knots[0].X;

    ay[r] = 1.0;
    by[r] = 4.0;
    cy[r] = 1.0;
    dy[r] = 4.0 * knots[r].Y + 2.0 * knots[0].Y;

    SolveTridiagonalCyclic(n, ax, bx, cx, dx, dx);
    SolveTridiagonalCyclic(n, ay, by, cy, dy, dy);

    curves.Restart(knots[0]);

    for (U32 i = 0; i < r; i++) {
        const FloatPoint cp0(dx[i], dy[i]);
        const FloatPoint cp1(2.0 * knots[i + 1].X - dx[i + 1],
            2.0 *knots[i + 1].Y - dy[i + 1]);

        curves.CurveTo(cp0, cp1, knots[i + 1]);
    }

    const FloatPoint cp0(dx[r], dy[r]);
    const FloatPoint cp1(2.0 * knots[0].X - dx[0],
        2.0 * knots[0].Y - dy[0]);

    curves.CurveTo(cp0, cp1, knots[0]);

    Dealloc(dy);
    Dealloc(cy);
    Dealloc(by);
    Dealloc(ay);
    Dealloc(dx);
    Dealloc(cx);
    Dealloc(bx);
    Dealloc(ax);
}


static void DoSolveCubicC2Acyclic(const U32 count,
    const FloatPoint *knots, CurveArray<CubicCurve> &curves)
{
    const U32 n = count - 1;

    double *ax = AllocElements<double>(n);
    double *bx = AllocElements<double>(n);
    double *cx = AllocElements<double>(n);
    double *dx = AllocElements<double>(n);
    double *ay = AllocElements<double>(n);
    double *by = AllocElements<double>(n);
    double *cy = AllocElements<double>(n);
    double *dy = AllocElements<double>(n);

    ax[0] = 0.0;
    bx[0] = 2.0;
    cx[0] = 1.0;
    dx[0] = knots[0].X + 2.0 * knots[1].X;

    ay[0] = 0.0;
    by[0] = 2.0;
    cy[0] = 1.0;
    dy[0] = knots[0].Y + 2.0 * knots[1].Y;

    for (U32 i = 1; i < n-1; i++) {
        ax[i] = 1.0;
        bx[i] = 4.0;
        cx[i] = 1.0;
        dx[i] = 4.0 * knots[i].X + 2.0 * knots[i + 1].X;

        ay[i] = 1.0;
        by[i] = 4.0;
        cy[i] = 1.0;
        dy[i] = 4.0 * knots[i].Y + 2.0 * knots[i + 1].Y;
    }

    const U32 r = n - 1;

    ax[r] = 2.0;
    bx[r] = 7.0;
    cx[r] = 0.0;
    dx[r] = 8.0 * knots[r].X + knots[n].X;

    ay[r] = 2.0;
    by[r] = 7.0;
    cy[r] = 0.0;
    dy[r] = 8.0 * knots[r].Y + knots[n].Y;

    SolveTridiagonal(n, ax, bx, cx, dx, dx);
    SolveTridiagonal(n, ay, by, cy, dy, dy);

    curves.Restart(knots[0]);

    for (U32 i = 0; i < r; i++) {
        const FloatPoint cp0(dx[i], dy[i]);
        const FloatPoint cp1(2.0 * knots[i + 1].X - dx[i + 1],
            2.0 * knots[i + 1].Y - dy[i + 1]);

        curves.CurveTo(cp0, cp1, knots[i + 1]);
    }

    const FloatPoint cp0(dx[r], dy[r]);
    const FloatPoint cp1((knots[n].X + dx[r]) / 2.0,
        (knots[n].Y + dy[r]) / 2.0);

    curves.CurveTo(cp0, cp1, knots[n]);

    Dealloc(dy);
    Dealloc(cy);
    Dealloc(by);
    Dealloc(ay);
    Dealloc(dx);
    Dealloc(cx);
    Dealloc(bx);
    Dealloc(ax);
}


void SolveCubicC2Cyclic(const U32 count, const FloatPoint *knots,
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
            DoSolveCubicC2Cyclic(count, knots, curves);
        }
    }
}


void SolveCubicC2Acyclic(const U32 count, const FloatPoint *knots,
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
            DoSolveCubicC2Acyclic(count, knots, curves);
        }
    }
}
