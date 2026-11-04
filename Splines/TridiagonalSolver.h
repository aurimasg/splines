
#pragma once


#include "Base.h"


/**
 * Solves a tridiagonal system of linear equations with the Thomas
 * algorithm.
 *
 * See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 *
 * @param count Number of equations. It must be at least 2.
 *
 * @param a Subdiagonal, at least count values. The first value is not read.
 * It must not be null.
 *
 * @param b Main diagonal, at least count values. It must not be null.
 *
 * @param c Superdiagonal, at least count values. The last value is not
 * read. It must not be null.
 *
 * @param d Right side, at least count values. It must not be null.
 *
 * @param x Receives the count values of the solution. It can be the same
 * array as any of the inputs. It must not be null.
 */
extern void SolveTridiagonal(const U32 count, const double *const a,
    const double *const b, const double *const c, const double *const d,
    double *const x) noexcept;


/**
 * Solves a cyclic tridiagonal system of linear equations with the Thomas
 * algorithm and the Sherman-Morrison formula. Unlike SolveTridiagonal, it
 * reads the first value of a and the last value of c: they link the first
 * and the last equations.
 *
 * See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 * See https://en.wikipedia.org/wiki/Sherman–Morrison_formula
 *
 * @param count Number of equations. It must be at least 2.
 *
 * @param a Subdiagonal, at least count values. It must not be null.
 *
 * @param b Main diagonal, at least count values. It must not be null.
 *
 * @param c Superdiagonal, at least count values. It must not be null.
 *
 * @param d Right side, at least count values. It must not be null.
 *
 * @param x Receives the count values of the solution. It can be the same
 * array as any of the inputs. It must not be null.
 */
extern void SolveTridiagonalCyclic(const U32 count, const double *const a,
    const double *const b, const double *const c, const double *const d,
    double *const x) noexcept;
