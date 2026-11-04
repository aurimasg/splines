
#pragma once


#include "Base.h"
#include "CubicCurve.h"
#include "FloatPoint.h"
#include "QuadraticCurve.h"


/**
 * Owns a connected sequence of curves in one point array.
 */
template <typename T>
struct CurveArray final {
private:

    /**
     * Returns degree of a supported curve type.
     */
    static constexpr U32 GetDegree() noexcept;

public:

    /**
     * Number of points appended for each curve, excluding its start point.
     */
    static constexpr U32 Degree = GetDegree();

    STATIC_ASSERT(__is_same(T, QuadraticCurve) or
        __is_same(T, CubicCurve));


    /**
     * Constructs an empty array without an allocation.
     */
    constexpr CurveArray() noexcept;


    /**
     * Releases all storage owned by this array.
     */
   ~CurveArray() noexcept;


    /**
     * Returns the number of complete curves. An empty array and an array
     * with only a start point both have zero curves.
     */
    constexpr U32 GetCurveCount() const noexcept;


    /**
     * Returns the number of stored points, including the start point.
     */
    constexpr U32 GetPointCount() const noexcept;


    /**
     * Returns a copy of one stored point.
     *
     * @param index Point index.
     */
    FloatPoint GetPointAtIndex(const U32 index) const noexcept;


    /**
     * Replaces one stored point without changing the number of curves.
     * A shared endpoint changes both adjacent curves.
     *
     * @param index Point index.
     *
     * @param point Point to copy.
     */
    void SetPointAtIndex(const U32 index, const FloatPoint point) noexcept;


    /**
     * Returns a curve value with copies of its stored points.
     *
     * @param index Curve index.
     */
    T GetCurveAtIndex(const U32 index) const noexcept;


    /**
     * Removes all points and keeps allocated storage. Call Restart before
     * you append another curve.
     */
    void Erase() noexcept;


    /**
     * Replaces all points with one start point. Keeps allocated storage
     * for reuse. The resulting curve count is zero.
     *
     * @param point Start point to copy. Both coordinates must be finite.
     */
    void Restart(const FloatPoint point) noexcept;


    /**
     * Appends one curve from the current endpoint.
     *
     * @param points Exactly Degree FloatPoint values.
     */
    template <typename... PointTypes>
    void CurveTo(const PointTypes... points) noexcept
        requires (sizeof...(PointTypes) == Degree and
            (__is_same(PointTypes, FloatPoint) and ...));

private:

    void ReserveForPointCount(const U32 pointCount) noexcept;

private:

    FloatPoint *mPoints = nullptr;
    U32 mPointCount = 0U;
    U32 mCapacity = 0U;

private:

    DISABLE_COPY_AND_ASSIGN(CurveArray);
};


template <>
inline constexpr U32 CurveArray<QuadraticCurve>::GetDegree() noexcept {
    return 2U;
}


template <>
inline constexpr U32 CurveArray<CubicCurve>::GetDegree() noexcept {
    return 3U;
}


template <typename T>
inline constexpr CurveArray<T>::CurveArray() noexcept = default;


template <typename T>
inline CurveArray<T>::~CurveArray() noexcept {
    delete [] mPoints;
}


template <typename T>
inline constexpr U32 CurveArray<T>::GetCurveCount() const noexcept {
    if (mPointCount == 0U) {
        return 0U;
    }

    return (mPointCount - 1U) / Degree;
}


template <typename T>
inline constexpr U32 CurveArray<T>::GetPointCount() const noexcept {
    return mPointCount;
}


template <typename T>
inline FloatPoint CurveArray<T>::GetPointAtIndex(const U32 index) const noexcept {
    ASSERT(index < mPointCount);

    return mPoints[index];
}


template <typename T>
inline void CurveArray<T>::SetPointAtIndex(const U32 index, const FloatPoint point) noexcept {
    ASSERT(index < mPointCount);
    ASSERT(NumberIsFinite(point.X));
    ASSERT(NumberIsFinite(point.Y));

    mPoints[index] = point;
}


template <>
inline CubicCurve CurveArray<CubicCurve>::GetCurveAtIndex(const U32 index) const noexcept {
    ASSERT(index < GetCurveCount());

    const FloatPoint *const points = mPoints + index * Degree;

    return CubicCurve(points[0], points[1], points[2], points[3]);
}


template <typename T>
inline void CurveArray<T>::Erase() noexcept {
    mPointCount = 0U;
}


template <typename T>
inline void CurveArray<T>::Restart(const FloatPoint point) noexcept {
    ASSERT(NumberIsFinite(point.X));
    ASSERT(NumberIsFinite(point.Y));

    ReserveForPointCount(1U);

    mPoints[0] = point;
    mPointCount = 1U;
}


template <typename T>
template <typename... PointTypes>
inline void CurveArray<T>::CurveTo(const PointTypes... points) noexcept
    requires (sizeof...(PointTypes) == Degree and
        (__is_same(PointTypes, FloatPoint) and ...))
{
    ASSERT(mPointCount > 0U);
    ASSERT(mPointCount <= UINT32_MAX - Degree);

    const FloatPoint appendedPoints[Degree] = { points... };
    const U32 pointCount = mPointCount + Degree;

    ReserveForPointCount(pointCount);

    for (U32 index = 0U; index < Degree; index++) {
        const FloatPoint point = appendedPoints[index];

        ASSERT(NumberIsFinite(point.X));
        ASSERT(NumberIsFinite(point.Y));

        mPoints[mPointCount + index] = point;
    }

    mPointCount = pointCount;
}


template <typename T>
inline void CurveArray<T>::ReserveForPointCount(const U32 pointCount) noexcept {
    if (pointCount <= mCapacity) {
        return;
    }

    U32 capacity = UINT32_MAX;

    if (mCapacity <= UINT32_MAX / 2U) {
        capacity = Max(mCapacity * 2U, 32U);
    }

    capacity = Max(capacity, pointCount);

    FloatPoint *const points = new FloatPoint[capacity];

    if (mPointCount > 0U) {
        STATIC_ASSERT(IsTriviallyCopyable<FloatPoint>());

        const USize byteCount = USize(mPointCount) * sizeof(FloatPoint);

        memcpy(points, mPoints, byteCount);
    }

    delete [] mPoints;

    mPoints = points;
    mCapacity = capacity;
}
