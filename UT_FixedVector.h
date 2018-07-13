/*
 * Copyright (c) 2018 Side Effects Software Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * COMMENTS:
 *      A vector class templated on its size and data type.
 */

#pragma once

#ifndef __UT_FixedVector__
#define __UT_FixedVector__

#include "SYS_Math.h"
#include "SYS_Types.h"

template<typename T,exint SIZE,bool INSTANTIATED=false>
class UT_FixedVector
{
public:
    typedef UT_FixedVector<T,SIZE,INSTANTIATED> ThisType;
    typedef T value_type;
    typedef T theType;
    static const exint theSize = SIZE;

    T vec[SIZE];

    SYS_FORCE_INLINE UT_FixedVector() = default;

    /// Initializes every component to the same value
    SYS_FORCE_INLINE explicit UT_FixedVector(T that) noexcept
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] = that;
    }

    SYS_FORCE_INLINE UT_FixedVector(const ThisType &that) = default;
    SYS_FORCE_INLINE UT_FixedVector(ThisType &&that) = default;

    /// Converts vector of S into vector of T,
    /// or just copies if same type.
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE UT_FixedVector(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) noexcept
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] = that[i];
    }

    template<typename S>
    SYS_FORCE_INLINE UT_FixedVector(const S that[SIZE]) noexcept
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] = that[i];
    }

    SYS_FORCE_INLINE const T &operator[](exint i) const noexcept
    {
        UT_ASSERT_P(i >= 0 && i < SIZE);
        return vec[i];
    }
    SYS_FORCE_INLINE T &operator[](exint i) noexcept
    {
        UT_ASSERT_P(i >= 0 && i < SIZE);
        return vec[i];
    }

    SYS_FORCE_INLINE constexpr const T *data() const noexcept
    {
        return vec;
    }
    SYS_FORCE_INLINE T *data() noexcept
    {
        return vec;
    }

    SYS_FORCE_INLINE ThisType &operator=(const ThisType &that) = default;
    SYS_FORCE_INLINE ThisType &operator=(ThisType &&that) = default;

    template <typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE ThisType &operator=(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) noexcept
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] = that[i];
        return *this;
    }
    SYS_FORCE_INLINE const ThisType &operator=(T that) noexcept
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] = that;
        return *this;
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE void operator+=(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that)
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] += that[i];
    }
    SYS_FORCE_INLINE void operator+=(T that)
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] += that;
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE auto operator+(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) const -> UT_FixedVector<decltype(vec[0]+that[0]),SIZE>
    {
        using Type = decltype(vec[0]+that[0]);
        UT_FixedVector<Type,SIZE> result;
        for (exint i = 0; i < SIZE; ++i)
            result[i] = vec[i] + that[i];
        return result;
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE void operator-=(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that)
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] -= that[i];
    }
    SYS_FORCE_INLINE void operator-=(T that)
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] -= that;
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE auto operator-(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) const -> UT_FixedVector<decltype(vec[0]-that[0]),SIZE>
    {
        using Type = decltype(vec[0]-that[0]);
        UT_FixedVector<Type,SIZE> result;
        for (exint i = 0; i < SIZE; ++i)
            result[i] = vec[i] - that[i];
        return result;
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE void operator*=(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that)
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] *= that[i];
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE auto operator*(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) const -> UT_FixedVector<decltype(vec[0]*that[0]),SIZE>
    {
        using Type = decltype(vec[0]*that[0]);
        UT_FixedVector<Type,SIZE> result;
        for (exint i = 0; i < SIZE; ++i)
            result[i] = vec[i] * that[i];
        return result;
    }
    SYS_FORCE_INLINE void operator*=(T that)
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] *= that;
    }
    SYS_FORCE_INLINE UT_FixedVector<T,SIZE> operator*(T that) const
    {
        UT_FixedVector<T,SIZE> result;
        for (exint i = 0; i < SIZE; ++i)
            result[i] = vec[i] * that;
        return result;
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE void operator/=(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that)
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] /= that[i];
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE auto operator/(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) const -> UT_FixedVector<decltype(vec[0]/that[0]),SIZE>
    {
        using Type = decltype(vec[0]/that[0]);
        UT_FixedVector<Type,SIZE> result;
        for (exint i = 0; i < SIZE; ++i)
            result[i] = vec[i] / that[i];
        return result;
    }

    SYS_FORCE_INLINE void operator/=(T that)
    {
        if (std::is_integral<T>::value)
        {
            for (exint i = 0; i < SIZE; ++i)
                vec[i] /= that;
        }
        else
        {
            that = 1/that;
            for (exint i = 0; i < SIZE; ++i)
                vec[i] *= that;
        }
    }
    SYS_FORCE_INLINE UT_FixedVector<T,SIZE> operator/(T that) const
    {
        UT_FixedVector<T,SIZE> result;
        if (std::is_integral<T>::value)
        {
            for (exint i = 0; i < SIZE; ++i)
                result[i] = vec[i] / that;
        }
        else
        {
            that = 1/that;
            for (exint i = 0; i < SIZE; ++i)
                result[i] = vec[i] * that;
        }
        return result;
    }
    SYS_FORCE_INLINE void negate()
    {
        for (exint i = 0; i < SIZE; ++i)
            vec[i] = -vec[i];
    }

    SYS_FORCE_INLINE UT_FixedVector<T,SIZE> operator-() const
    {
        UT_FixedVector<T,SIZE> result;
        for (exint i = 0; i < SIZE; ++i)
            result[i] = -vec[i];
        return result;
    }

    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE bool operator==(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) const noexcept
    {
        for (exint i = 0; i < SIZE; ++i)
        {
            if (vec[i] != T(that[i]))
                return false;
        }
        return true;
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE bool operator!=(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) const noexcept
    {
        return !(*this==that);
    }
    SYS_FORCE_INLINE bool isZero() const noexcept
    {
        for (exint i = 0; i < SIZE; ++i)
        {
            if (vec[i] != T(0))
                return false;
        }
        return true;
    }
    SYS_FORCE_INLINE T maxComponent() const
    {
        T v = vec[0];
        for (exint i = 1; i < SIZE; ++i)
            v = (vec[i] > v) ? vec[i] : v;
        return v;
    }
    SYS_FORCE_INLINE T minComponent() const
    {
        T v = vec[0];
        for (exint i = 1; i < SIZE; ++i)
            v = (vec[i] < v) ? vec[i] : v;
        return v;
    }
    SYS_FORCE_INLINE T avgComponent() const
    {
        T v = vec[0];
        for (exint i = 1; i < SIZE; ++i)
            v += vec[i];
        return v / SIZE;
    }

    SYS_FORCE_INLINE T length2() const noexcept
    {
        T a0(vec[0]);
        T result(a0*a0);
        for (exint i = 1; i < SIZE; ++i)
        {
            T ai(vec[i]);
            result += ai*ai;
        }
        return result;
    }
    SYS_FORCE_INLINE T length() const
    {
        T len2 = length2();
        return SYSsqrt(len2);
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE auto dot(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) const -> decltype(vec[0]*that[0])
    {
        using TheType = decltype(vec[0]*that.vec[0]);
        TheType result(vec[0]*that[0]);
        for (exint i = 1; i < SIZE; ++i)
            result += vec[i]*that[i];
        return result;
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE auto distance2(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) const -> decltype(vec[0]-that[0])
    {
        using TheType = decltype(vec[0]-that[0]);
        TheType v(vec[0] - that[0]);
        TheType result(v*v);
        for (exint i = 1; i < SIZE; ++i)
        {
            v = vec[i] - that[i];
            result += v*v;
        }
        return result;
    }
    template<typename S,bool S_INSTANTIATED>
    SYS_FORCE_INLINE auto distance(const UT_FixedVector<S,SIZE,S_INSTANTIATED> &that) const -> decltype(vec[0]-that[0])
    {
        auto dist2 = distance2(that);
        return SYSsqrt(dist2);
    }

    SYS_FORCE_INLINE T normalize()
    {
        T len2 = length2();
        if (len2 == T(0))
            return T(0);
        if (len2 == T(1))
            return T(1);
        T len = SYSsqrt(len2);
        // Check if the square root is equal 1.  sqrt(1+dx) ~ 1+dx/2,
        // so it may get rounded to 1 when it wasn't 1 before.
        if (len != T(1))
            (*this) /= len;
        return len;
    }
};

/// NOTE: Strictly speaking, this should use decltype(that*a[0]),
///       but in the interests of avoiding accidental precision escalation,
///       it uses T.
template<typename T,exint SIZE,bool INSTANTIATED,typename S>
SYS_FORCE_INLINE UT_FixedVector<T,SIZE> operator*(const S &that,const UT_FixedVector<T,SIZE,INSTANTIATED> &a)
{
    T t(that);
    UT_FixedVector<T,SIZE> result;
    for (exint i = 0; i < SIZE; ++i)
        result[i] = t * a[i];
    return result;
}

template<typename T, exint SIZE, bool INSTANTIATED, typename S, bool S_INSTANTIATED>
SYS_FORCE_INLINE auto
dot(const UT_FixedVector<T,SIZE,INSTANTIATED> &a, const UT_FixedVector<S,SIZE,S_INSTANTIATED> &b) -> decltype(a[0]*b[0])
{
    return a.dot(b);
}

template<typename T, exint SIZE, bool INSTANTIATED, typename S, bool S_INSTANTIATED>
SYS_FORCE_INLINE auto
SYSmin(const UT_FixedVector<T,SIZE,INSTANTIATED> &a, const UT_FixedVector<S,SIZE,S_INSTANTIATED> &b) -> UT_FixedVector<decltype(a[0]+b[1]), SIZE>
{
    using Type = decltype(a[0]+b[1]);
    UT_FixedVector<Type, SIZE> result;
    for (exint i = 0; i < SIZE; ++i)
        result[i] = SYSmin(Type(a[i]), Type(b[i]));
    return result;
}

template<typename T, exint SIZE, bool INSTANTIATED, typename S, bool S_INSTANTIATED>
SYS_FORCE_INLINE auto
SYSmax(const UT_FixedVector<T,SIZE,INSTANTIATED> &a, const UT_FixedVector<S,SIZE,S_INSTANTIATED> &b) -> UT_FixedVector<decltype(a[0]+b[1]), SIZE>
{
    using Type = decltype(a[0]+b[1]);
    UT_FixedVector<Type, SIZE> result;
    for (exint i = 0; i < SIZE; ++i)
        result[i] = SYSmax(Type(a[i]), Type(b[i]));
    return result;
}

template<typename T>
struct UT_FixedVectorTraits
{
    typedef UT_FixedVector<T,1> FixedVectorType;
    typedef T DataType;
    static const exint TupleSize = 1;
    static const bool isVectorType = false;
};

template<typename T,exint SIZE,bool INSTANTIATED>
struct UT_FixedVectorTraits<UT_FixedVector<T,SIZE,INSTANTIATED> >
{
    typedef UT_FixedVector<T,SIZE,INSTANTIATED> FixedVectorType;
    typedef T DataType;
    static const exint TupleSize = SIZE;
    static const bool isVectorType = true;
};

#endif
