
/**
 * Solves a spline through a fixed set of knots and prints it as an SVG
 * document to standard output.
 *
 * Usage
 *
 *     splines [-s TYPE] [-c]
 *
 * -s TYPE selects the spline type by its short name. See SplineTypeNames
 * for the list. The default is s9. -c closes the spline so that the final
 * curve returns to the first knot.
 */

#include "../Splines/SplineCubicC2.h"
#include "../Splines/SplineCubicG2.h"
#include "../Splines/SplineHobby.h"
#include "../Splines/SplineKappa.h"
#include "../Splines/SplineS9.h"

#include <cstdio>
#include <cstring>


/**
 * Knot positions in SVG user units, with the Y axis pointing down. They
 * trace three loops along an arc, like a cursive stroke: a lead-in, then
 * for each loop its base, its top and its left side, and a final base.
 * Consecutive knots are distinct, as every solver requires, and the final
 * knot is distinct from the first so that closed splines are valid too.
 */
static constexpr FloatPoint Knots[] = {
    FloatPoint( 40.0, 229.0),
    FloatPoint( 98.0, 208.0),
    FloatPoint( 64.0,  76.0),
    FloatPoint( 52.0, 163.0),
    FloatPoint(184.0, 171.0),
    FloatPoint(225.0,  40.0),
    FloatPoint(169.0, 107.0),
    FloatPoint(277.0, 185.0),
    FloatPoint(381.0,  95.0),
    FloatPoint(298.0, 123.0),
    FloatPoint(348.0, 245.0)
};


static constexpr U32 KnotCount = U32(sizeof(Knots) / sizeof(Knots[0]));


// The open Kappa spline needs the most knots of all solvers.
STATIC_ASSERT(KnotCount >= 4U);


static constexpr double HobbyOmega = 1.0;


enum class SplineType {
    S9,
    CubicC2,
    CubicG2,
    Hobby,
    Kappa
};


struct SplineTypeName final {
    const char *Name;
    SplineType Type;
    const char *Description;
};


static constexpr SplineTypeName SplineTypeNames[] = {
    { "s9",    SplineType::S9,      "S9 spline, the default" },
    { "c2",    SplineType::CubicC2, "Cubic spline with C2 continuity" },
    { "g2",    SplineType::CubicG2, "Cubic spline with G2 continuity" },
    { "hobby", SplineType::Hobby,   "Hobby spline" },
    { "kappa", SplineType::Kappa,   "Kappa spline" }
};


static constexpr U32 SplineTypeNameCount =
    U32(sizeof(SplineTypeNames) / sizeof(SplineTypeNames[0]));


static bool TryFindSplineType(const char *const name, SplineType &type)
{
    for (U32 i = 0U; i < SplineTypeNameCount; i++) {
        if (strcmp(name, SplineTypeNames[i].Name) == 0) {
            type = SplineTypeNames[i].Type;

            return true;
        }
    }

    return false;
}


static void PrintUsage(FILE *const file, const char *const program)
{
    fprintf(file, "Usage: %s [-s TYPE] [-c]\n", program);
    fprintf(file, "\n");
    fprintf(file, "Solves a spline through built-in knots and prints it as SVG.\n");
    fprintf(file, "\n");
    fprintf(file, "  -s TYPE  Spline type. One of:\n");

    for (U32 i = 0U; i < SplineTypeNameCount; i++) {
        fprintf(file, "           %-6s %s\n", SplineTypeNames[i].Name,
            SplineTypeNames[i].Description);
    }

    fprintf(file, "  -c       Close the spline.\n");
    fprintf(file, "  -h       Print this help.\n");
}


template <typename T>
static FloatPoint PointOnCurve(const CurveArray<T> &curves, const U32 index,
    const double t)
{
    static constexpr U32 Degree = CurveArray<T>::Degree;

    FloatPoint p[Degree + 1U];

    for (U32 i = 0U; i <= Degree; i++) {
        p[i] = curves.GetPointAtIndex(index * Degree + i);
    }

    // De Casteljau.
    for (U32 r = 1U; r <= Degree; r++) {
        for (U32 i = 0U; i <= Degree - r; i++) {
            p[i] = InterpolateLinear(p[i], p[i + 1U], t);
        }
    }

    return p[0];
}


template <typename T>
static void PrintSvg(const CurveArray<T> &curves, const char *const typeName,
    const bool isClosed)
{
    static constexpr U32 Degree = CurveArray<T>::Degree;
    static constexpr char Command = Degree == 3U ? 'C' : 'Q';
    static constexpr double StrokeWidth = 2.0;

    static constexpr U32 SampleCount = 32U;

    const U32 curveCount = curves.GetCurveCount();

    ASSERT(curveCount > 0U);

    // Bounds of the curves themselves, sampled, since control points can
    // lie far outside them.
    FloatPoint min = curves.GetPointAtIndex(0U);
    FloatPoint max = min;

    for (U32 i = 0U; i < curveCount; i++) {
        for (U32 k = 1U; k <= SampleCount; k++) {
            const FloatPoint p = PointOnCurve(curves, i,
                double(k) / double(SampleCount));

            min.X = Min(min.X, p.X);
            min.Y = Min(min.Y, p.Y);
            max.X = Max(max.X, p.X);
            max.Y = Max(max.Y, p.Y);
        }
    }

    const double centerX = (min.X + max.X) / 2.0;
    const double centerY = (min.Y + max.Y) / 2.0;

    printf("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    printf("<!-- %s spline, %s, %u knots -->\n", typeName,
        isClosed ? "closed" : "open", unsigned(KnotCount));

    printf("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"100%%\" "
        "height=\"100%%\">\n");
    printf("  <svg x=\"50%%\" y=\"50%%\" overflow=\"visible\">\n");

    const FloatPoint start = curves.GetPointAtIndex(0U);

    printf("    <path transform=\"translate(%.3f %.3f)\" d=\"M %.3f %.3f\n",
        -centerX, -centerY, start.X, start.Y);

    for (U32 i = 0U; i < curveCount; i++) {
        printf("      %c", Command);

        for (U32 j = 1U; j <= Degree; j++) {
            const FloatPoint p = curves.GetPointAtIndex(i * Degree + j);

            printf(" %.3f %.3f", p.X, p.Y);
        }

        printf("\n");
    }

    if (isClosed) {
        printf("      Z\n");
    }

    printf("    \" fill=\"none\" stroke=\"black\" stroke-width=\"%.1f\" "
        "stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n",
        StrokeWidth);
    printf("  </svg>\n");
    printf("</svg>\n");
}


template <typename T>
static int PrintSolvedSpline(const CurveArray<T> &curves,
    const char *const typeName, const bool isClosed)
{
    if (curves.GetCurveCount() == 0U) {
        fprintf(stderr, "Failed to solve %s %s spline through %u knots.\n",
            isClosed ? "closed" : "open", typeName, unsigned(KnotCount));

        return 1;
    }

    PrintSvg(curves, typeName, isClosed);

    return 0;
}


int main(int argc, char **argv)
{
    const char *const program = argc > 0 ? argv[0] : "splines";
    const char *typeName = SplineTypeNames[0].Name;
    SplineType type = SplineTypeNames[0].Type;
    bool isClosed = false;

    for (int i = 1; i < argc; i++) {
        const char *const arg = argv[i];

        if (strcmp(arg, "-c") == 0) {
            isClosed = true;
        } else if (strcmp(arg, "-s") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "%s: option -s needs a spline type.\n",
                    program);
                PrintUsage(stderr, program);

                return 1;
            }

            typeName = argv[++i];

            if (!TryFindSplineType(typeName, type)) {
                fprintf(stderr, "%s: unknown spline type '%s'.\n", program,
                    typeName);
                PrintUsage(stderr, program);

                return 1;
            }
        } else if (strcmp(arg, "-h") == 0) {
            PrintUsage(stdout, program);

            return 0;
        } else {
            fprintf(stderr, "%s: unknown option '%s'.\n", program, arg);
            PrintUsage(stderr, program);

            return 1;
        }
    }

    CurveArray<CubicCurve> cubics;
    CurveArray<QuadraticCurve> quadratics;

    switch (type) {
        case SplineType::S9:
            if (isClosed) {
                S9SolveCyclic(KnotCount, Knots, cubics);
            } else {
                S9SolveAcyclic(KnotCount, Knots, cubics);
            }

            break;
        case SplineType::CubicC2:
            if (isClosed) {
                SolveCubicC2Cyclic(KnotCount, Knots, cubics);
            } else {
                SolveCubicC2Acyclic(KnotCount, Knots, cubics);
            }

            break;
        case SplineType::CubicG2:
            if (isClosed) {
                SolveCubicG2Cyclic(KnotCount, Knots, cubics);
            } else {
                SolveCubicG2Acyclic(KnotCount, Knots, cubics);
            }

            break;
        case SplineType::Hobby:
            if (isClosed) {
                SolveHobbyCyclic(KnotCount, Knots, cubics);
            } else {
                SolveHobbyAcyclic(KnotCount, Knots, cubics, HobbyOmega);
            }

            break;
        case SplineType::Kappa:
            if (isClosed) {
                SolveKappaCyclic(Knots, KnotCount, quadratics);
            } else {
                SolveKappaAcyclic(Knots, KnotCount, quadratics);
            }

            break;
    }

    if (type == SplineType::Kappa) {
        return PrintSolvedSpline(quadratics, typeName, isClosed);
    }

    return PrintSolvedSpline(cubics, typeName, isClosed);
}
