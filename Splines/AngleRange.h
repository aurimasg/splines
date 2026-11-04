
#pragma once


#include "Base.h"


/**
 * Minimum-maximum angle pair representing a range of allowed angles for
 * solver.
 */
struct AngleRange final {
    AngleRange(const double min, const double max);

    double MinAngle = 0;
    double MaxAngle = 0;

private:

    AngleRange() = delete;
};


inline AngleRange::AngleRange(const double min, const double max)
:   MinAngle(min),
    MaxAngle(max)
{
    ASSERT(min <= max);
}
