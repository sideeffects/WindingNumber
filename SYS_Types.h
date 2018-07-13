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
 *      Common type definitions.
 */

#pragma once

#ifndef __SYS_Types__
#define __SYS_Types__

/* Include system types */
#include <limits>
#include <type_traits>
#include <sys/types.h>

/*
 * Integer types
 */
typedef signed char	int8;
typedef	unsigned char	uint8;
typedef short		int16;
typedef unsigned short	uint16;
typedef	int		int32;
typedef unsigned int	uint32;

#ifndef MBSD
typedef unsigned int	uint;
#endif

/*
 * Avoid using uint64.
 * The extra bit of precision is NOT worth the cost in pain and suffering
 * induced by use of unsigned.
 */
#if defined(_WIN32)
    typedef __int64		int64;
    typedef unsigned __int64	uint64;
#elif defined(MBSD)
    // On MBSD, int64/uint64 are also defined in the system headers so we must
    // declare these in the same way or else we get conflicts.
    #include <stdint.h>
    typedef int64_t		int64;
    typedef uint64_t		uint64;
#elif defined(AMD64)
    typedef long		int64;
    typedef unsigned long	uint64;
#else
    typedef long long		int64;
    typedef unsigned long long	uint64;
#endif

/// The problem with int64 is that it implies that it is a fixed 64-bit quantity
/// that is saved to disk. Therefore, we need another integral type for
/// indexing our arrays.
typedef int64 exint;

/// Mark function to be inlined. If this is done, taking the address of such
/// a function is not allowed.
#if defined(__GNUC__) || defined(__clang__)
#define SYS_FORCE_INLINE	__attribute__ ((always_inline)) inline
#elif defined(_MSC_VER)
#define SYS_FORCE_INLINE	__forceinline
#else
#define SYS_FORCE_INLINE	inline
#endif

/// Floating Point Types
typedef float   fpreal32;
typedef double  fpreal64;

/// SYS_FPRealUnionT for type-safe casting with integral types
template <typename T>
union SYS_FPRealUnionT;

template <>
union SYS_FPRealUnionT<fpreal32>
{
    typedef int32	int_type;
    typedef uint32	uint_type;
    typedef fpreal32	fpreal_type;

    enum {
	EXPONENT_BITS = 8,
	MANTISSA_BITS = 23,
	EXPONENT_BIAS = 127 };

    int_type		ival;
    uint_type		uval;
    fpreal_type		fval;
    
    struct
    {
	uint_type mantissa_val: 23;
	uint_type exponent_val: 8;
	uint_type sign_val: 1;
    };
};

template <>
union SYS_FPRealUnionT<fpreal64>
{
    typedef int64	int_type;
    typedef uint64	uint_type;
    typedef fpreal64	fpreal_type;

    enum {
	EXPONENT_BITS = 11,
	MANTISSA_BITS = 52,
	EXPONENT_BIAS = 1023 };

    int_type		ival;
    uint_type		uval;
    fpreal_type		fval;
    
    struct
    {
	uint_type mantissa_val: 52;
	uint_type exponent_val: 11;
	uint_type sign_val: 1;
    };
};

typedef union SYS_FPRealUnionT<fpreal32>    SYS_FPRealUnionF;
typedef union SYS_FPRealUnionT<fpreal64>    SYS_FPRealUnionD;

/// Asserts are disabled
/// @{
#define UT_ASSERT_P(ZZ)         ((void)0)
#define UT_ASSERT(ZZ)           ((void)0)
#define UT_ASSERT_MSG_P(ZZ, MM) ((void)0)
#define UT_ASSERT_MSG(ZZ, MM)   ((void)0)
/// @}

#endif
