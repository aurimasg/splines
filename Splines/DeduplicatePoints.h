
#pragma once


#include "Base.h"
#include "FloatPoint.h"


static constexpr double DuplicatePointMaximumSquaredDistance = 1e-7;


extern U32 DeduplicatePoints(const U32 count, const FloatPoint *s,
    FloatPoint *d,
    const double maximumSquaredDistance = DuplicatePointMaximumSquaredDistance);
