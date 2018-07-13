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
 *      SIMD wrapper classes for 4 floats or 4 ints
 */

#pragma once

#ifndef __HDK_VM_SIMD__
#define __HDK_VM_SIMD__

#include "SYS_Math.h"
#include <cstdint>

//#define FORCE_NON_SIMD

#include "VM_SSEFunc.h"

class v4uf;

class v4uu {
public:
    SYS_FORCE_INLINE v4uu() {}
    SYS_FORCE_INLINE v4uu(const v4si &v) : vector(v) {}
    SYS_FORCE_INLINE v4uu(const v4uu &v) : vector(v.vector) {}
    explicit SYS_FORCE_INLINE v4uu(int32 v) { vector = VM_SPLATS(v); }
    explicit SYS_FORCE_INLINE v4uu(const int32 v[4])
    { vector = VM_LOAD(v); }
    SYS_FORCE_INLINE v4uu(int32 a, int32 b, int32 c, int32 d)
    { vector = VM_SPLATS(a, b, c, d); }

    // Assignment
    SYS_FORCE_INLINE v4uu operator=(int32 v)
    { vector = v4uu(v).vector; return *this; }
    SYS_FORCE_INLINE v4uu operator=(v4si v)
    { vector = v; return *this; }
    SYS_FORCE_INLINE v4uu operator=(const v4uu &v)
    { vector = v.vector; return *this; }

    SYS_FORCE_INLINE void condAssign(const v4uu &val, const v4uu &c)
    { *this = (c & val) | ((!c) & *this); }

    // Comparison
    SYS_FORCE_INLINE v4uu operator == (const v4uu &v) const
    { return v4uu(VM_ICMPEQ(vector, v.vector)); }
    SYS_FORCE_INLINE v4uu operator != (const v4uu &v) const
    { return ~(*this == v); }
    SYS_FORCE_INLINE v4uu operator >  (const v4uu &v) const
    { return v4uu(VM_ICMPGT(vector, v.vector)); }
    SYS_FORCE_INLINE v4uu operator <  (const v4uu &v) const
    { return v4uu(VM_ICMPLT(vector, v.vector)); }
    SYS_FORCE_INLINE v4uu operator >= (const v4uu &v) const
    { return ~(*this < v); }
    SYS_FORCE_INLINE v4uu operator <= (const v4uu &v) const
    { return ~(*this > v); }

    SYS_FORCE_INLINE v4uu operator == (int32 v) const { return *this == v4uu(v); }
    SYS_FORCE_INLINE v4uu operator != (int32 v) const { return *this != v4uu(v); }
    SYS_FORCE_INLINE v4uu operator >  (int32 v) const { return *this > v4uu(v); }
    SYS_FORCE_INLINE v4uu operator <  (int32 v) const { return *this < v4uu(v); }
    SYS_FORCE_INLINE v4uu operator >= (int32 v) const { return *this >= v4uu(v); }
    SYS_FORCE_INLINE v4uu operator <= (int32 v) const { return *this <= v4uu(v); }

    // Basic math
    SYS_FORCE_INLINE v4uu operator+(const v4uu &r) const
    { return v4uu(VM_IADD(vector, r.vector)); }
    SYS_FORCE_INLINE v4uu operator-(const v4uu &r) const
    { return v4uu(VM_ISUB(vector, r.vector)); }
    SYS_FORCE_INLINE v4uu operator+=(const v4uu &r) { return (*this = *this + r); }
    SYS_FORCE_INLINE v4uu operator-=(const v4uu &r) { return (*this = *this - r); }
    SYS_FORCE_INLINE v4uu operator+(int32 r) const { return *this + v4uu(r); }
    SYS_FORCE_INLINE v4uu operator-(int32 r) const { return *this - v4uu(r); }
    SYS_FORCE_INLINE v4uu operator+=(int32 r) { return (*this = *this + r); }
    SYS_FORCE_INLINE v4uu operator-=(int32 r) { return (*this = *this - r); }

    // logical/bitwise

    SYS_FORCE_INLINE v4uu operator||(const v4uu &r) const
    { return v4uu(VM_OR(vector, r.vector)); }
    SYS_FORCE_INLINE v4uu operator&&(const v4uu &r) const
    { return v4uu(VM_AND(vector, r.vector)); }
    SYS_FORCE_INLINE v4uu operator^(const v4uu &r) const
    { return v4uu(VM_XOR(vector, r.vector)); }
    SYS_FORCE_INLINE v4uu operator!() const
    { return *this == v4uu(0); }

    SYS_FORCE_INLINE v4uu operator|(const v4uu &r) const { return *this || r; }
    SYS_FORCE_INLINE v4uu operator&(const v4uu &r) const { return *this && r; }
    SYS_FORCE_INLINE v4uu operator~() const
    { return *this ^ v4uu(0xFFFFFFFF); }

    // component
    SYS_FORCE_INLINE int32 operator[](int idx) const { return VM_EXTRACT(vector, idx); }
    SYS_FORCE_INLINE void setComp(int idx, int32 v) { vector = VM_INSERT(vector, v, idx); }

    v4uf toFloat() const;

public:
    v4si vector;
};

class v4uf {
public:
    SYS_FORCE_INLINE v4uf() {}
    SYS_FORCE_INLINE v4uf(const v4sf &v) : vector(v) {}
    SYS_FORCE_INLINE v4uf(const v4uf &v) : vector(v.vector) {}
    explicit SYS_FORCE_INLINE v4uf(float v) { vector = VM_SPLATS(v); }
    explicit SYS_FORCE_INLINE v4uf(const float v[4])
    { vector = VM_LOAD(v); }
    SYS_FORCE_INLINE v4uf(float a, float b, float c, float d)
    { vector = VM_SPLATS(a, b, c, d); }

    // Assignment
    SYS_FORCE_INLINE v4uf operator=(float v)
    { vector = v4uf(v).vector; return *this; }
    SYS_FORCE_INLINE v4uf operator=(v4sf v)
    { vector = v; return *this; }
    SYS_FORCE_INLINE v4uf operator=(const v4uf &v)
    { vector = v.vector; return *this; }

    SYS_FORCE_INLINE void condAssign(const v4uf &val, const v4uu &c)
    { *this = (val & c) | (*this & ~c); }

    // Comparison
    SYS_FORCE_INLINE v4uu operator == (const v4uf &v) const
    { return v4uu(VM_CMPEQ(vector, v.vector)); }
    SYS_FORCE_INLINE v4uu operator != (const v4uf &v) const
    { return v4uu(VM_CMPNE(vector, v.vector)); }
    SYS_FORCE_INLINE v4uu operator >  (const v4uf &v) const
    { return v4uu(VM_CMPGT(vector, v.vector)); }
    SYS_FORCE_INLINE v4uu operator <  (const v4uf &v) const
    { return v4uu(VM_CMPLT(vector, v.vector)); }
    SYS_FORCE_INLINE v4uu operator >= (const v4uf &v) const
    { return v4uu(VM_CMPGE(vector, v.vector)); }
    SYS_FORCE_INLINE v4uu operator <= (const v4uf &v) const
    { return v4uu(VM_CMPLE(vector, v.vector)); }

    SYS_FORCE_INLINE v4uu operator == (float v) const { return *this == v4uf(v); }
    SYS_FORCE_INLINE v4uu operator != (float v) const { return *this != v4uf(v); }
    SYS_FORCE_INLINE v4uu operator >  (float v) const { return *this > v4uf(v); }
    SYS_FORCE_INLINE v4uu operator <  (float v) const { return *this < v4uf(v); }
    SYS_FORCE_INLINE v4uu operator >= (float v) const { return *this >= v4uf(v); }
    SYS_FORCE_INLINE v4uu operator <= (float v) const { return *this <= v4uf(v); }


    // Basic math
    SYS_FORCE_INLINE v4uf operator+(const v4uf &r) const
    { return v4uf(VM_ADD(vector, r.vector)); }
    SYS_FORCE_INLINE v4uf operator-(const v4uf &r) const
    { return v4uf(VM_SUB(vector, r.vector)); }
    SYS_FORCE_INLINE v4uf operator-() const
    { return v4uf(VM_NEG(vector)); }
    SYS_FORCE_INLINE v4uf operator*(const v4uf &r) const
    { return v4uf(VM_MUL(vector, r.vector)); }
    SYS_FORCE_INLINE v4uf operator/(const v4uf &r) const
    { return v4uf(VM_DIV(vector, r.vector)); }

    SYS_FORCE_INLINE v4uf operator+=(const v4uf &r) { return (*this = *this + r); }
    SYS_FORCE_INLINE v4uf operator-=(const v4uf &r) { return (*this = *this - r); }
    SYS_FORCE_INLINE v4uf operator*=(const v4uf &r) { return (*this = *this * r); }
    SYS_FORCE_INLINE v4uf operator/=(const v4uf &r) { return (*this = *this / r); }

    SYS_FORCE_INLINE v4uf operator+(float r) const { return *this + v4uf(r); }
    SYS_FORCE_INLINE v4uf operator-(float r) const { return *this - v4uf(r); }
    SYS_FORCE_INLINE v4uf operator*(float r) const { return *this * v4uf(r); }
    SYS_FORCE_INLINE v4uf operator/(float r) const { return *this / v4uf(r); }
    SYS_FORCE_INLINE v4uf operator+=(float r) { return (*this = *this + r); }
    SYS_FORCE_INLINE v4uf operator-=(float r) { return (*this = *this - r); }
    SYS_FORCE_INLINE v4uf operator*=(float r) { return (*this = *this * r); }
    SYS_FORCE_INLINE v4uf operator/=(float r) { return (*this = *this / r); }

    // logical/bitwise

    SYS_FORCE_INLINE v4uf operator||(const v4uu &r) const
    { return v4uf(V4SF(VM_OR(V4SI(vector), r.vector))); }
    SYS_FORCE_INLINE v4uf operator&&(const v4uu &r) const
    { return v4uf(V4SF(VM_AND(V4SI(vector), r.vector))); }
    SYS_FORCE_INLINE v4uf operator^(const v4uu &r) const
    { return v4uf(V4SF(VM_XOR(V4SI(vector), r.vector))); }
    SYS_FORCE_INLINE v4uf operator!() const
    { return v4uf(V4SF((*this == v4uf(0.0F)).vector)); }

    SYS_FORCE_INLINE v4uf operator||(const v4uf &r) const
    { return v4uf(V4SF(VM_OR(V4SI(vector), V4SI(r.vector)))); }
    SYS_FORCE_INLINE v4uf operator&&(const v4uf &r) const
    { return v4uf(V4SF(VM_AND(V4SI(vector), V4SI(r.vector)))); }
    SYS_FORCE_INLINE v4uf operator^(const v4uf &r) const
    { return v4uf(V4SF(VM_XOR(V4SI(vector), V4SI(r.vector)))); }

    SYS_FORCE_INLINE v4uf operator|(const v4uu &r) const { return *this || r; }
    SYS_FORCE_INLINE v4uf operator&(const v4uu &r) const { return *this && r; }
    SYS_FORCE_INLINE v4uf operator~() const
    { return *this ^ v4uu(0xFFFFFFFF); }

    SYS_FORCE_INLINE v4uf operator|(const v4uf &r) const { return *this || r; }
    SYS_FORCE_INLINE v4uf operator&(const v4uf &r) const { return *this && r; }

    // component
    SYS_FORCE_INLINE float operator[](int idx) const { return VM_EXTRACT(vector, idx); }
    SYS_FORCE_INLINE void setComp(int idx, float v) { vector = VM_INSERT(vector, v, idx); }

    // more math
    SYS_FORCE_INLINE v4uf abs() const { return v4uf(VM_ABS(vector)); }
    SYS_FORCE_INLINE v4uf clamp(const v4uf &low, const v4uf &high) const
    { return v4uf(
        VM_MIN(VM_MAX(vector, low.vector), high.vector)); }
    SYS_FORCE_INLINE v4uf clamp(float low, float high) const
    { return v4uf(VM_MIN(VM_MAX(vector,
        v4uf(low).vector), v4uf(high).vector)); }
    SYS_FORCE_INLINE v4uf recip() const { return v4uf(VM_INVERT(vector)); }

    /// This is a lie, it is a signed int.
    SYS_FORCE_INLINE v4uu toUnsignedInt() const { return VM_INT(vector); }
    SYS_FORCE_INLINE v4uu toSignedInt() const { return VM_INT(vector); }

    v4uu floor() const
    {
        VM_P_FLOOR();
        v4uu result = VM_FLOOR(vector);
        VM_E_FLOOR();
        return result;
    }

    /// Returns the integer part of this float, this becomes the
    /// 0..1 fractional component.
    v4uu splitFloat()
    {
        v4uu base = toSignedInt();
        *this -= base.toFloat();
        return base;
    }

    template <int A, int B, int C, int D>
    SYS_FORCE_INLINE v4uf swizzle() const
    { 
        return VM_SHUFFLE<A,B,C,D>(vector);
    }

    SYS_FORCE_INLINE v4uu isFinite() const
    {
        // If the exponent is the maximum value, it's either infinite or NaN.
        const v4si mask = VM_SPLATS(0x7F800000);
        return ~v4uu(VM_ICMPEQ(VM_AND(V4SI(vector), mask), mask));
    }

public:
    v4sf vector;
};

SYS_FORCE_INLINE v4uf
v4uu::toFloat() const
{
    return v4uf(VM_IFLOAT(vector));
}

//
// Custom vector operations
//

static SYS_FORCE_INLINE v4uf
sqrt(const v4uf &a)
{
    return v4uf(VM_SQRT(a.vector));
}

static SYS_FORCE_INLINE v4uf
fabs(const v4uf &a)
{
    return a.abs();
}

// Use this operation to mask disabled values to 0
// rval = !a ? b : 0;

static SYS_FORCE_INLINE v4uf
andn(const v4uu &a, const v4uf &b)
{
    return v4uf(V4SF(VM_ANDNOT(a.vector, V4SI(b.vector))));
}

static SYS_FORCE_INLINE v4uu
andn(const v4uu &a, const v4uu &b)
{
    return v4uu(VM_ANDNOT(a.vector, b.vector));
}

// rval = a ? b : c;
static SYS_FORCE_INLINE v4uf
ternary(const v4uu &a, const v4uf &b, const v4uf &c)
{
    return (b & a) | andn(a, c);
}

static SYS_FORCE_INLINE v4uu
ternary(const v4uu &a, const v4uu &b, const v4uu &c)
{
    return (b & a) | andn(a, c);
}

// rval = !(a && b)
static SYS_FORCE_INLINE v4uu
nand(const v4uu &a, const v4uu &b)
{
    return !v4uu(VM_AND(a.vector, b.vector));
}

static SYS_FORCE_INLINE v4uf
vmin(const v4uf &a, const v4uf &b)
{
    return v4uf(VM_MIN(a.vector, b.vector));
}

static SYS_FORCE_INLINE v4uf
vmax(const v4uf &a, const v4uf &b)
{
    return v4uf(VM_MAX(a.vector, b.vector));
}

static SYS_FORCE_INLINE v4uf
clamp(const v4uf &a, const v4uf &b, const v4uf &c)
{
    return vmax(vmin(a, c), b);
}

static SYS_FORCE_INLINE v4uf
clamp(const v4uf &a, float b, float c)
{
    return vmax(vmin(a, v4uf(c)), v4uf(b));
}

static SYS_FORCE_INLINE bool
allbits(const v4uu &a)
{
    return vm_allbits(a.vector);
}

static SYS_FORCE_INLINE bool
anybits(const v4uu &a)
{
    return !allbits(~a);
}

static SYS_FORCE_INLINE v4uf
madd(const v4uf &v, const v4uf &f, const v4uf &a)
{
    return v4uf(VM_MADD(v.vector, f.vector, a.vector));
}

static SYS_FORCE_INLINE v4uf
madd(const v4uf &v, float f, float a)
{
    return v4uf(VM_MADD(v.vector, v4uf(f).vector, v4uf(a).vector));
}

static SYS_FORCE_INLINE v4uf
madd(const v4uf &v, float f, const v4uf &a)
{
    return v4uf(VM_MADD(v.vector, v4uf(f).vector, a.vector));
}

static SYS_FORCE_INLINE v4uf
msub(const v4uf &v, const v4uf &f, const v4uf &s)
{
    return madd(v, f, -s);
}

static SYS_FORCE_INLINE v4uf
msub(const v4uf &v, float f, float s)
{
    return madd(v, f, -s);
}

static SYS_FORCE_INLINE v4uf
lerp(const v4uf &a, const v4uf &b, const v4uf &w)
{
    v4uf w1 = v4uf(1.0F) - w;
    return madd(a, w1, b*w);
}

static SYS_FORCE_INLINE v4uf
luminance(const v4uf &r, const v4uf &g, const v4uf &b,
    float rw, float gw, float bw)
{
    return v4uf(madd(r, v4uf(rw), madd(g, v4uf(gw), b * bw)));
}

static SYS_FORCE_INLINE float
dot3(const v4uf &a, const v4uf &b)
{
    v4uf res = a*b;
    return res[0] + res[1] + res[2];
}

static SYS_FORCE_INLINE float
dot4(const v4uf &a, const v4uf &b)
{
    v4uf res = a*b;
    return res[0] + res[1] + res[2] + res[3];
}

static SYS_FORCE_INLINE float
length(const v4uf &a)
{
    return SYSsqrt(dot3(a, a));
}

static SYS_FORCE_INLINE v4uf
normalize(const v4uf &a)
{
    return a / length(a);
}

static SYS_FORCE_INLINE v4uf
cross(const v4uf &a, const v4uf &b)
{
    return v4uf(a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0], 0);
}

// Currently there is no specific support for signed integers
typedef v4uu v4ui;

// Assuming that ptr is an array of elements of type STYPE, this operation
// will return the index of the first element that is aligned to (1<<ASIZE)
// bytes.
#define VM_ALIGN(ptr, ASIZE, STYPE)	\
		((((1<<ASIZE)-(intptr_t)ptr)&((1<<ASIZE)-1))/sizeof(STYPE))

#endif
