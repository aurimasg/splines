
#pragma once


#include <cstdint>
#include <cstring>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <new>


#ifdef DEBUG
#define ASSERT(p) assert(p)
#else
#define ASSERT(p)
#endif


#define STATIC_ASSERT(p) static_assert(p, "")


#define DISABLE_COPY_AND_ASSIGN(c) \
    c(const c &a) = delete; \
    void operator=(const c &a);


using I8 = int8_t;
using U8 = uint8_t;
using U32 = uint32_t;
using USize = size_t;
using Size_t = size_t;


template <typename T>
static constexpr U32 SizeOf() {
    return static_cast<U32>(sizeof(T));
}


STATIC_ASSERT(SizeOf<I8>() == 1);
STATIC_ASSERT(SizeOf<U8>() == 1);
STATIC_ASSERT(SizeOf<U32>() == 4);


/**
 * Returns absolute of a given value.
 */
template <typename T>
static constexpr T Abs(const T t) {
    return t >= 0 ? t : -t;
}


/**
 * Returns greater of two given values.
 */
template <typename T>
static constexpr T Max(const T a, const T b) {
    return a < b ? b : a;
}


/**
 * Returns smaller of two given values.
 */
template <typename T>
static constexpr T Min(const T a, const T b) {
    return a < b ? a : b;
}


/**
 * Returns value clamped to range between minimum and maximum values.
 */
template <typename T>
static constexpr T Clamp(const T val, const T min, const T max) {
    return val > max ? max : val < min ? min : val;
}


/**
 * Returns true if two given numbers are considered equal.
 */
static inline bool FuzzyIsEqual(const double a, const double b) {
    return (fabs(a - b) < DBL_EPSILON);
}


/**
 * Returns true if a number can be considered being equal to zero.
 */
static inline bool FuzzyIsZero(const double d) {
    return fabs(d) < DBL_EPSILON;
}


/**
 * Linearly interpolate between A and B.
 * If t is 0, returns A.
 * If t is 1, returns B.
 * If t is something else, returns value linearly interpolated between A and B.
 */
template <typename T, typename V>
static inline T InterpolateLinear(const T A, const T B, const V t) {
    return A + ((B - A) * t);
}


/**
 * Returns square of a given value.
 */
template <typename T>
static constexpr double Square(const T v) {
    return v * v;
}


/**
 * Returns square root of a given number.
 */
static inline double Sqrt(const double v) {
    return sqrt(v);
}


/**
 * Calculates sine of a given number.
 */
static inline double Sin(const double v) {
    return sin(v);
}


/**
 * Calculate cosine of a given number.
 */
static inline double Cos(const double v) {
    return cos(v);
}


/**
 * Returns the angle whose cosine is the specified number. Returned angle is
 * in radians.
 *
 * @param v A number representing a cosine which must be greater than or
 * equal to -1 and less than or equal to 1.
 */
static inline double Acos(const double v) {
    ASSERT(v <= 1.0);
    ASSERT(v >= -1.0);

    return acos(v);
}


/**
 * Returns the angle expressed in radians between the positive X-axis and the
 * ray cast from 0, 0 to the point x, y.
 */
static inline double Atan2(const double y, const double x) {
    return atan2(y, x);
}


static inline double Hypot(const double x, const double y) {
    return hypot(x, y);
}


/**
 * Returns true if a given double precision floating point value is not a
 * number.
 */
static constexpr bool IsNaN(const double x) {
    return x != x;
}


/**
 * Returns true if a given double precision floating point number is finite.
 */
static inline bool NumberIsFinite(const double x) {
    // 0 × finite → 0
    // 0 × infinity → NaN
    // 0 × NaN → NaN
    const double p = x * 0;

    return !IsNaN(p);
}


/**
 * Ratio of the circumference of a circle to its diameter.
 */
static constexpr double Pi = 3.141592653589793238462643383279502884;


static inline void *Alloc(const Size_t size) {
    return malloc(size);
}


static inline void Dealloc(void *p) {
    free(p);
}


template <typename T>
static inline T *AllocElements(const Size_t count) {
    ASSERT(count > 0);

    return static_cast<T *>(Alloc(SizeOf<T>() * count));
}


template <typename T>
static inline T *AllocElementsInit(const Size_t count) {
    ASSERT(count > 0);

    T *p = static_cast<T *>(Alloc(SizeOf<T>() * count));

    for (Size_t i = 0; i < count; i++) {
        (void) new (p + i) T();
    }

    return p;
}


template <typename T>
constexpr inline bool IsTriviallyCopyable() noexcept {
    return __is_trivially_copyable(T);
}
