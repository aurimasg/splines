
#include "FloatPoint.h"


FloatPoint FloatPoint::UnitVector() const
{
    const double mag2 = LengthSquared();

    if (mag2 != 0.0 and mag2 != 1.0) {
        const double length = Sqrt(mag2);

        return FloatPoint(X / length, Y / length);
    }

    return *this;
}
