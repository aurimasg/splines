
#include "DeduplicatePoints.h"


U32 DeduplicatePoints(const U32 count, const FloatPoint *s,
    FloatPoint *d, const double maximumSquaredDistance)
{
    if (count < 1) {
        return 0;
    }

    ASSERT(s != nullptr);
    ASSERT(d != nullptr);

    const double epsilon = Clamp(maximumSquaredDistance, 1e-15, 1e+10);

    FloatPoint prev = s[0];
    U32 index = 1;

    d[0] = prev;

    for (U32 i = 1; i < count; i++) {
        const FloatPoint p = s[i];
        const double distance = p.DistanceToSquared(prev);

        if (distance > epsilon) {
            d[index++] = p;

            prev = p;
        }
    }

    return index;
}
