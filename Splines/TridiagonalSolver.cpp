
#include "TridiagonalSolver.h"


void SolveTridiagonal(const U32 count, const double *const a,
    const double *const b, const double *const c, const double *const d,
    double *const x) noexcept
{
    ASSERT(count >= 2U);
    ASSERT(a != nullptr);
    ASSERT(b != nullptr);
    ASSERT(c != nullptr);
    ASSERT(d != nullptr);
    ASSERT(x != nullptr);

    // Forward sweep.
    //
    // It writes the superdiagonal and the right side of the eliminated system
    // to their own arrays, thus x can share memory with the inputs - x is
    // written only after the sweep.
    double *const superdiagonal = AllocElements<double>(count);
    double *const rightSide = AllocElements<double>(count);

    double previousSuperdiagonal = c[0] / b[0];
    double previousRightSide = d[0] / b[0];

    superdiagonal[0] = previousSuperdiagonal;
    rightSide[0] = previousRightSide;

    for (U32 index = 1U; index < count; index++) {
        const double denominator = b[index] - a[index] * previousSuperdiagonal;

        previousSuperdiagonal = c[index] / denominator;
        previousRightSide = (d[index] - a[index] * previousRightSide) /
            denominator;

        superdiagonal[index] = previousSuperdiagonal;
        rightSide[index] = previousRightSide;
    }

    // Back substitution from the last equation.
    double previousSolution = previousRightSide;

    x[count - 1U] = previousSolution;

    for (U32 index = count - 1U; index > 0U; index--) {
        const U32 previous = index - 1U;

        previousSolution = rightSide[previous] -
            superdiagonal[previous] * previousSolution;

        x[previous] = previousSolution;
    }

    Dealloc(rightSide);
    Dealloc(superdiagonal);
}


void SolveTridiagonalCyclic(const U32 count, const double *const a,
    const double *const b, const double *const c, const double *const d,
    double *const x) noexcept
{
    ASSERT(count >= 2U);
    ASSERT(a != nullptr);
    ASSERT(b != nullptr);
    ASSERT(c != nullptr);
    ASSERT(d != nullptr);
    ASSERT(x != nullptr);

    const U32 last = count - 1U;

    // With two equations each cyclic link joins the same two unknowns as the
    // diagonal next to it, thus the links fold into the subdiagonal and the
    // superdiagonal.
    if (count == 2U) {
        const double subdiagonal[2] = {
            0.0, a[1] + c[1]
        };

        const double superdiagonal[2] = {
            a[0] + c[0], 0.0
        };

        SolveTridiagonal(2U, subdiagonal, b, superdiagonal, d, x);

        return;
    }

    const double firstSubdiagonal = a[0];
    const double lastSuperdiagonal = c[last];

    // Without cyclic links the system is plain tridiagonal.
    if (FuzzyIsZero(firstSubdiagonal) and FuzzyIsZero(lastSuperdiagonal)) {
        SolveTridiagonal(count, a, b, c, d, x);

        return;
    }

    // Sherman-Morrison. The cyclic links are a rank one update of a
    // tridiagonal matrix with a modified diagonal. The update vector is
    // solved in place of a right side, and the two solutions combine.
    double *const update = AllocElements<double>(count);
    double *const diagonal = AllocElements<double>(count);

    memcpy(diagonal, b, SizeOf<double>() * count);

    diagonal[0] -= firstSubdiagonal;
    diagonal[last] -= lastSuperdiagonal;

    memset(update, 0, SizeOf<double>() * count);

    update[0] = firstSubdiagonal;
    update[last] = lastSuperdiagonal;

    SolveTridiagonal(count, a, diagonal, c, update, update);
    SolveTridiagonal(count, a, diagonal, c, d, x);

    const double factor = (x[0] + x[last]) / (1.0 + update[0] + update[last]);

    for (U32 index = 0U; index < count; index++) {
        x[index] -= factor * update[index];
    }

    Dealloc(diagonal);
    Dealloc(update);
}
