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
 *      SIMD wrapper functions for SSE instructions
 */

#pragma once

#ifndef __VM_SSEFunc__
#define __VM_SSEFunc__

#include "SYS_Types.h"

#if defined(_MSC_VER)
    #pragma warning(push)
    #pragma warning(disable:4799)
#endif

#define CPU_HAS_SIMD_INSTR	1
#define VM_SSE_STYLE		1

#include <emmintrin.h>
typedef __m128	v4sf;
typedef __m128i	v4si;

#if defined(__SSE4_1__)
#define VM_SSE41_STYLE		1
#include <smmintrin.h>
#endif

#if defined(_MSC_VER)
    #pragma warning(pop)
#endif

// Plain casting (no conversion)
// MSVC has problems casting between __m128 and __m128i, so we implement a
// custom casting routine specifically for windows.

#if defined(_MSC_VER)

static SYS_FORCE_INLINE v4sf
vm_v4sf(const v4si &a)
{
    union {
	v4si ival;
	v4sf fval;
    };
    ival = a;
    return fval;
}

static SYS_FORCE_INLINE v4si
vm_v4si(const v4sf &a)
{
    union {
	v4si ival;
	v4sf fval;
    };
    fval = a;
    return ival;
}

#define V4SF(A)		vm_v4sf(A)
#define V4SI(A)		vm_v4si(A)

#else

#define V4SF(A)		(v4sf)A
#define V4SI(A)		(v4si)A

#endif

#define VM_SHUFFLE_MASK(a0,a1, b0,b1)	((b1)<<6|(b0)<<4 | (a1)<<2|(a0))

template <int mask>
static SYS_FORCE_INLINE v4sf
vm_shuffle(const v4sf &a, const v4sf &b)
{
    return _mm_shuffle_ps(a, b, mask);
}

template <int mask>
static SYS_FORCE_INLINE v4si
vm_shuffle(const v4si &a, const v4si &b)
{
    return V4SI(_mm_shuffle_ps(V4SF(a), V4SF(b), mask));
}

template <int A, int B, int C, int D, typename T>
static SYS_FORCE_INLINE T
vm_shuffle(const T &a, const T &b)
{
    return vm_shuffle<VM_SHUFFLE_MASK(A,B,C,D)>(a, b);
}

template <int mask, typename T>
static SYS_FORCE_INLINE T
vm_shuffle(const T &a)
{
    return vm_shuffle<mask>(a, a);
}

template <int A, int B, int C, int D, typename T>
static SYS_FORCE_INLINE T
vm_shuffle(const T &a)
{
    return vm_shuffle<A,B,C,D>(a, a);
}

#if defined(VM_SSE41_STYLE)

static SYS_FORCE_INLINE v4si
vm_insert(const v4si v, int32 a, int n)
{
    switch (n)
    {
    case 0: return _mm_insert_epi32(v, a, 0);
    case 1: return _mm_insert_epi32(v, a, 1);
    case 2: return _mm_insert_epi32(v, a, 2);
    case 3: return _mm_insert_epi32(v, a, 3);
    }
    return v;
}

static SYS_FORCE_INLINE v4sf
vm_insert(const v4sf v, float a, int n)
{
    switch (n)
    {
    case 0: return _mm_insert_ps(v, _mm_set_ss(a), _MM_MK_INSERTPS_NDX(0,0,0));
    case 1: return _mm_insert_ps(v, _mm_set_ss(a), _MM_MK_INSERTPS_NDX(0,1,0));
    case 2: return _mm_insert_ps(v, _mm_set_ss(a), _MM_MK_INSERTPS_NDX(0,2,0));
    case 3: return _mm_insert_ps(v, _mm_set_ss(a), _MM_MK_INSERTPS_NDX(0,3,0));
    }
    return v;
}

static SYS_FORCE_INLINE int
vm_extract(const v4si v, int n)
{
    switch (n)
    {
    case 0: return _mm_extract_epi32(v, 0);
    case 1: return _mm_extract_epi32(v, 1);
    case 2: return _mm_extract_epi32(v, 2);
    case 3: return _mm_extract_epi32(v, 3);
    }
    return 0;
}

static SYS_FORCE_INLINE float
vm_extract(const v4sf v, int n)
{
    SYS_FPRealUnionF	tmp;
    switch (n)
    {
    case 0: tmp.ival = _mm_extract_ps(v, 0); break;
    case 1: tmp.ival = _mm_extract_ps(v, 1); break;
    case 2: tmp.ival = _mm_extract_ps(v, 2); break;
    case 3: tmp.ival = _mm_extract_ps(v, 3); break;
    }
    return tmp.fval;
}

#else

static SYS_FORCE_INLINE v4si
vm_insert(const v4si v, int32 a, int n)
{
    union { v4si vector; int32 comp[4]; };
    vector = v;
    comp[n] = a;
    return vector;
}

static SYS_FORCE_INLINE v4sf
vm_insert(const v4sf v, float a, int n)
{
    union { v4sf vector; float comp[4]; };
    vector = v;
    comp[n] = a;
    return vector;
}

static SYS_FORCE_INLINE int
vm_extract(const v4si v, int n)
{
    union { v4si vector; int32 comp[4]; };
    vector = v;
    return comp[n];
}

static SYS_FORCE_INLINE float
vm_extract(const v4sf v, int n)
{
    union { v4sf vector; float comp[4]; };
    vector = v;
    return comp[n];
}

#endif

static SYS_FORCE_INLINE v4sf
vm_splats(float a)
{
    return _mm_set1_ps(a);
}

static SYS_FORCE_INLINE v4si
vm_splats(uint32 a)
{
    SYS_FPRealUnionF	tmp;
    tmp.uval = a;
    return V4SI(vm_splats(tmp.fval));
}

static SYS_FORCE_INLINE v4si
vm_splats(int32 a)
{
    SYS_FPRealUnionF	tmp;
    tmp.ival = a;
    return V4SI(vm_splats(tmp.fval));
}

static SYS_FORCE_INLINE v4sf
vm_splats(float a, float b, float c, float d)
{
    return vm_shuffle<0,2,0,2>(
	    vm_shuffle<0>(_mm_set_ss(a), _mm_set_ss(b)),
	    vm_shuffle<0>(_mm_set_ss(c), _mm_set_ss(d)));
}

static SYS_FORCE_INLINE v4si
vm_splats(uint32 a, uint32 b, uint32 c, uint32 d)
{
    SYS_FPRealUnionF	af, bf, cf, df;
    af.uval = a;
    bf.uval = b;
    cf.uval = c;
    df.uval = d;
    return V4SI(vm_splats(af.fval, bf.fval, cf.fval, df.fval));
}

static SYS_FORCE_INLINE v4si
vm_splats(int32 a, int32 b, int32 c, int32 d)
{
    SYS_FPRealUnionF	af, bf, cf, df;
    af.ival = a;
    bf.ival = b;
    cf.ival = c;
    df.ival = d;
    return V4SI(vm_splats(af.fval, bf.fval, cf.fval, df.fval));
}

static SYS_FORCE_INLINE v4si
vm_load(const int32 v[4])
{
    return V4SI(_mm_loadu_ps((const float *)v));
}

static SYS_FORCE_INLINE v4sf
vm_load(const float v[4])
{
    return _mm_loadu_ps(v);
}

static SYS_FORCE_INLINE void
vm_store(float dst[4], v4sf value)
{
    _mm_storeu_ps(dst, value);
}

static SYS_FORCE_INLINE v4sf
vm_negate(v4sf a)
{
    return _mm_sub_ps(_mm_setzero_ps(), a);
}

static SYS_FORCE_INLINE v4sf
vm_abs(v4sf a)
{
    return _mm_max_ps(a, vm_negate(a));
}

static SYS_FORCE_INLINE v4sf
vm_fdiv(v4sf a, v4sf b)
{
    return _mm_mul_ps(a, _mm_rcp_ps(b));
}

static SYS_FORCE_INLINE v4sf
vm_fsqrt(v4sf a)
{
    return _mm_rcp_ps(_mm_rsqrt_ps(a));
}

static SYS_FORCE_INLINE v4sf
vm_madd(v4sf a, v4sf b, v4sf c)
{
    return _mm_add_ps(_mm_mul_ps(a, b), c);
}

static const v4si	theSSETrue = vm_splats(0xFFFFFFFF);

static SYS_FORCE_INLINE bool
vm_allbits(const v4si &a)
{
    return _mm_movemask_ps(V4SF(_mm_cmpeq_epi32(a, theSSETrue))) == 0xF;
}


#define VM_EXTRACT	vm_extract
#define VM_INSERT	vm_insert
#define VM_SPLATS	vm_splats
#define VM_LOAD		vm_load
#define VM_STORE	vm_store

#define VM_CMPLT(A,B)	V4SI(_mm_cmplt_ps(A,B))
#define VM_CMPLE(A,B)	V4SI(_mm_cmple_ps(A,B))
#define VM_CMPGT(A,B)	V4SI(_mm_cmpgt_ps(A,B))
#define VM_CMPGE(A,B)	V4SI(_mm_cmpge_ps(A,B))
#define VM_CMPEQ(A,B)	V4SI(_mm_cmpeq_ps(A,B))
#define VM_CMPNE(A,B)	V4SI(_mm_cmpneq_ps(A,B))

#define VM_ICMPLT	_mm_cmplt_epi32
#define VM_ICMPGT	_mm_cmpgt_epi32
#define VM_ICMPEQ	_mm_cmpeq_epi32

#define VM_IADD		_mm_add_epi32
#define VM_ISUB		_mm_sub_epi32

#define VM_ADD		_mm_add_ps
#define VM_SUB		_mm_sub_ps
#define VM_MUL		_mm_mul_ps
#define VM_DIV		_mm_div_ps
#define VM_SQRT		_mm_sqrt_ps
#define VM_ISQRT	_mm_rsqrt_ps
#define VM_INVERT	_mm_rcp_ps
#define VM_ABS		vm_abs

#define VM_FDIV		vm_fdiv
#define VM_NEG		vm_negate
#define VM_FSQRT	vm_fsqrt
#define VM_MADD		vm_madd

#define VM_MIN		_mm_min_ps
#define VM_MAX		_mm_max_ps

#define VM_AND		_mm_and_si128
#define VM_ANDNOT	_mm_andnot_si128
#define VM_OR		_mm_or_si128
#define VM_XOR		_mm_xor_si128

#define VM_ALLBITS	vm_allbits

#define VM_SHUFFLE	vm_shuffle

// Integer to float conversions
#define VM_SSE_ROUND_MASK	0x6000
#define VM_SSE_ROUND_ZERO	0x6000
#define VM_SSE_ROUND_UP		0x4000
#define VM_SSE_ROUND_DOWN	0x2000
#define VM_SSE_ROUND_NEAR	0x0000

#define GETROUND()	(_mm_getcsr()&VM_SSE_ROUND_MASK)
#define SETROUND(x)	(_mm_setcsr(x|(_mm_getcsr()&~VM_SSE_ROUND_MASK)))

// The P functions must be invoked before FLOOR, the E functions invoked
// afterwards to reset the state.

#define VM_P_FLOOR()	uint rounding = GETROUND(); \
			    SETROUND(VM_SSE_ROUND_DOWN);
#define VM_FLOOR	_mm_cvtps_epi32
#define VM_INT		_mm_cvttps_epi32
#define VM_E_FLOOR()	SETROUND(rounding);

// Float to integer conversion
#define VM_IFLOAT	_mm_cvtepi32_ps

#endif
