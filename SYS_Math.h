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
 *      Miscellaneous math functions. 
 */

#pragma once

#ifndef __SYS_Math__
#define __SYS_Math__

#include "SYS_Types.h"

#include <float.h>
#include <limits>
#include <math.h>

// NOTE:
// These have been carefully written so that in the case of equality
// we always return the first parameter.  This is so that NANs in
// in the second parameter are suppressed.
#define h_min(a, b)	(((a) > (b)) ? (b) : (a))
#define h_max(a, b)	(((a) < (b)) ? (b) : (a))
// DO NOT CHANGE THE ABOVE WITHOUT READING THE COMMENT
#define h_abs(a)	(((a) > 0) ? (a) : -(a))

static constexpr inline  int16 SYSmin(int16 a, int16 b)		{ return h_min(a,b); }
static constexpr inline  int16 SYSmax(int16 a, int16 b)		{ return h_max(a,b); }
static constexpr inline  int16 SYSabs(int16 a)			{ return h_abs(a); }
static constexpr inline  int32 SYSmin(int32 a, int32 b)		{ return h_min(a,b); }
static constexpr inline  int32 SYSmax(int32 a, int32 b)		{ return h_max(a,b); }
static constexpr inline  int32 SYSabs(int32 a)			{ return h_abs(a); }
static constexpr inline  int64 SYSmin(int64 a, int64 b)		{ return h_min(a,b); }
static constexpr inline  int64 SYSmax(int64 a, int64 b)		{ return h_max(a,b); }
static constexpr inline  int64 SYSmin(int32 a, int64 b)		{ return h_min(a,b); }
static constexpr inline  int64 SYSmax(int32 a, int64 b)		{ return h_max(a,b); }
static constexpr inline  int64 SYSmin(int64 a, int32 b)		{ return h_min(a,b); }
static constexpr inline  int64 SYSmax(int64 a, int32 b)		{ return h_max(a,b); }
static constexpr inline  int64 SYSabs(int64 a)			{ return h_abs(a); }
static constexpr inline uint16 SYSmin(uint16 a, uint16 b)		{ return h_min(a,b); }
static constexpr inline uint16 SYSmax(uint16 a, uint16 b)		{ return h_max(a,b); }
static constexpr inline uint32 SYSmin(uint32 a, uint32 b)		{ return h_min(a,b); }
static constexpr inline uint32 SYSmax(uint32 a, uint32 b)		{ return h_max(a,b); }
static constexpr inline uint64 SYSmin(uint64 a, uint64 b)		{ return h_min(a,b); }
static constexpr inline uint64 SYSmax(uint64 a, uint64 b)		{ return h_max(a,b); }
static constexpr inline fpreal32 SYSmin(fpreal32 a, fpreal32 b)	{ return h_min(a,b); }
static constexpr inline fpreal32 SYSmax(fpreal32 a, fpreal32 b)	{ return h_max(a,b); }
static constexpr inline fpreal64 SYSmin(fpreal64 a, fpreal64 b)	{ return h_min(a,b); }
static constexpr inline fpreal64 SYSmax(fpreal64 a, fpreal64 b)	{ return h_max(a,b); }

// Some systems have size_t as a seperate type from uint.  Some don't.
#if (defined(LINUX) && defined(IA64)) || defined(MBSD)
static constexpr inline size_t SYSmin(size_t a, size_t b)		{ return h_min(a,b); }
static constexpr inline size_t SYSmax(size_t a, size_t b)		{ return h_max(a,b); }
#endif

#undef h_min
#undef h_max
#undef h_abs

#define h_clamp(val, min, max, tol)	\
	    ((val <= min+tol) ? min : ((val >= max-tol) ? max : val))

    static constexpr inline int
    SYSclamp(int v, int min, int max)
	{ return h_clamp(v, min, max, 0); }

    static constexpr inline uint
    SYSclamp(uint v, uint min, uint max)
	{ return h_clamp(v, min, max, 0); }

    static constexpr inline int64
    SYSclamp(int64 v, int64 min, int64 max)
	{ return h_clamp(v, min, max, int64(0)); }

    static constexpr inline uint64
    SYSclamp(uint64 v, uint64 min, uint64 max)
	{ return h_clamp(v, min, max, uint64(0)); }

    static constexpr inline fpreal32
    SYSclamp(fpreal32 v, fpreal32 min, fpreal32 max, fpreal32 tol=(fpreal32)0)
	{ return h_clamp(v, min, max, tol); }

    static constexpr inline fpreal64
    SYSclamp(fpreal64 v, fpreal64 min, fpreal64 max, fpreal64 tol=(fpreal64)0)
	{ return h_clamp(v, min, max, tol); }

#undef h_clamp

static inline fpreal64 SYSsqrt(fpreal64 arg)
{ return ::sqrt(arg); }
static inline fpreal32 SYSsqrt(fpreal32 arg)
{ return ::sqrtf(arg); }
static inline fpreal64 SYSatan2(fpreal64 a, fpreal64 b)
{ return ::atan2(a, b); }
static inline fpreal32 SYSatan2(fpreal32 a, fpreal32 b)
{ return ::atan2(a, b); }

static inline fpreal32 SYSabs(fpreal32 a) { return ::fabsf(a); }
static inline fpreal64 SYSabs(fpreal64 a) { return ::fabs(a); }

#endif
