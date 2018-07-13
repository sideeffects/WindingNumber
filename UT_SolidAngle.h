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
 *      Functions and structures for computing solid angles.
 */

#pragma once

#ifndef __HDK_UT_SolidAngle_h__
#define __HDK_UT_SolidAngle_h__

#include "UT_BVH.h"

#include "UT_FixedVector.h"
#include "SYS_Math.h"
#include <memory>

namespace HDK_Sample {

template<typename T>
using UT_Vector2T = UT_FixedVector<T,2>;
template<typename T>
using UT_Vector3T = UT_FixedVector<T,3>;

template <typename T>
SYS_FORCE_INLINE T cross(const UT_Vector2T<T> &v1, const UT_Vector2T<T> &v2)
{
    return v1[0]*v2[1] - v1[1]*v2[0];
}

template <typename T>
SYS_FORCE_INLINE
UT_Vector3T<T> cross(const UT_Vector3T<T> &v1, const UT_Vector3T<T> &v2)
{
    UT_Vector3T<T> result;
    // compute the cross product:
    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return result;
}

/// Returns the signed solid angle subtended by triangle abc
/// from query point.
///
/// WARNING: This uses the right-handed normal convention, whereas most of
///          Houdini uses the left-handed normal convention, so either
///          negate the output, or swap b and c if you want it to be
///          positive inside and negative outside.
template<typename T>
T UTsignedSolidAngleTri(
    const UT_Vector3T<T> &a,
    const UT_Vector3T<T> &b,
    const UT_Vector3T<T> &c,
    const UT_Vector3T<T> &query)
{
    // Make a, b, and c relative to query
    UT_Vector3T<T> qa = a-query;
    UT_Vector3T<T> qb = b-query;
    UT_Vector3T<T> qc = c-query;

    const T alength = qa.length();
    const T blength = qb.length();
    const T clength = qc.length();

    // If any triangle vertices are coincident with query,
    // query is on the surface, which we treat as no solid angle.
    if (alength == 0 || blength == 0 || clength == 0)
        return T(0);

    // Normalize the vectors
    qa /= alength;
    qb /= blength;
    qc /= clength;

    // The formula on Wikipedia has roughly dot(qa,cross(qb,qc)),
    // but that's unstable when qa, qb, and qc are very close,
    // (e.g. if the input triangle was very far away).
    // This should be equivalent, but more stable.
    const T numerator = dot(qa, cross(qb-qa, qc-qa));

    // If numerator is 0, regardless of denominator, query is on the
    // surface, which we treat as no solid angle.
    if (numerator == 0)
        return T(0);

    const T denominator = T(1) + dot(qa,qb) + dot(qa,qc) + dot(qb,qc);

    return T(2)*SYSatan2(numerator, denominator);
}

template<typename T>
T UTsignedSolidAngleQuad(
    const UT_Vector3T<T> &a,
    const UT_Vector3T<T> &b,
    const UT_Vector3T<T> &c,
    const UT_Vector3T<T> &d,
    const UT_Vector3T<T> &query)
{
    // Make a, b, c, and d relative to query
    UT_Vector3T<T> v[4] = {
        a-query,
        b-query,
        c-query,
        d-query
    };

    const T lengths[4] = {
        v[0].length(),
        v[1].length(),
        v[2].length(),
        v[3].length()
    };

    // If any quad vertices are coincident with query,
    // query is on the surface, which we treat as no solid angle.
    // We could add the contribution from the non-planar part,
    // but in the context of a mesh, we'd still miss some, like
    // we do in the triangle case.
    if (lengths[0] == T(0) || lengths[1] == T(0) || lengths[2] == T(0) || lengths[3] == T(0))
        return T(0);

    // Normalize the vectors
    v[0] /= lengths[0];
    v[1] /= lengths[1];
    v[2] /= lengths[2];
    v[3] /= lengths[3];

    // Compute (unnormalized, but consistently-scaled) barycentric coordinates
    // for the query point inside the tetrahedron of points.
    // If 0 or 4 of the coordinates are positive, (or slightly negative), the
    // query is (approximately) inside, so the choice of triangulation matters.
    // Otherwise, the triangulation doesn't matter.

    const UT_Vector3T<T> diag02 = v[2]-v[0];
    const UT_Vector3T<T> diag13 = v[3]-v[1];
    const UT_Vector3T<T> v01 = v[1]-v[0];
    const UT_Vector3T<T> v23 = v[3]-v[2];

    T bary[4];
    bary[0] = dot(v[3],cross(v23,diag13));
    bary[1] = -dot(v[2],cross(v23,diag02));
    bary[2] = -dot(v[1],cross(v01,diag13));
    bary[3] = dot(v[0],cross(v01,diag02));

    const T dot01 = dot(v[0],v[1]);
    const T dot12 = dot(v[1],v[2]);
    const T dot23 = dot(v[2],v[3]);
    const T dot30 = dot(v[3],v[0]);

    T omega = T(0);

    // Equation of a bilinear patch in barycentric coordinates of its
    // tetrahedron is x0*x2 = x1*x3.  Less is one side; greater is other.
    if (bary[0]*bary[2] < bary[1]*bary[3])
    {
        // Split 0-2: triangles 0,1,2 and 0,2,3
        const T numerator012 = bary[3];
        const T numerator023 = bary[1];
        const T dot02 = dot(v[0],v[2]);

        // If numerator is 0, regardless of denominator, query is on the
        // surface, which we treat as no solid angle.
        if (numerator012 != T(0))
        {
            const T denominator012 = T(1) + dot01 + dot12 + dot02;
            omega = SYSatan2(numerator012, denominator012);
        }
        if (numerator023 != T(0))
        {
            const T denominator023 = T(1) + dot02 + dot23 + dot30;
            omega += SYSatan2(numerator023, denominator023);
        }
    }
    else
    {
        // Split 1-3: triangles 0,1,3 and 1,2,3
        const T numerator013 = -bary[2];
        const T numerator123 = -bary[0];
        const T dot13 = dot(v[1],v[3]);

        // If numerator is 0, regardless of denominator, query is on the
        // surface, which we treat as no solid angle.
        if (numerator013 != T(0))
        {
            const T denominator013 = T(1) + dot01 + dot13 + dot30;
            omega = SYSatan2(numerator013, denominator013);
        }
        if (numerator123 != T(0))
        {
            const T denominator123 = T(1) + dot12 + dot23 + dot13;
            omega += SYSatan2(numerator123, denominator123);
        }
    }
    return T(2)*omega;
}

/// Class for quickly approximating signed solid angle of a large mesh
/// from many query points.  This is useful for computing the
/// generalized winding number at many points.
///
/// NOTE: This is currently only instantiated for <float,float>.
template<typename T,typename S>
class UT_SolidAngle
{
public:
    /// This is outlined so that we don't need to include UT_BVHImpl.h
    UT_SolidAngle();
    /// This is outlined so that we don't need to include UT_BVHImpl.h
    ~UT_SolidAngle();

    /// NOTE: This does not take ownership over triangle_points or positions,
    ///       but does keep pointers to them, so the caller must keep them in
    ///       scope for the lifetime of this structure.
    UT_SolidAngle(
        const int ntriangles,
        const int *const triangle_points,
        const int npoints,
        const UT_Vector3T<S> *const positions,
        const int order = 2)
        : UT_SolidAngle()
    { init(ntriangles, triangle_points, npoints, positions, order); }

    /// Initialize the tree and data.
    /// NOTE: It is safe to call init on a UT_SolidAngle that has had init
    ///       called on it before, to re-initialize it.
    void init(
        const int ntriangles,
        const int *const triangle_points,
        const int npoints,
        const UT_Vector3T<S> *const positions,
        const int order = 2);

    /// Frees myTree and myData, and clears the rest.
    void clear();

    /// Returns true if this is clear
    bool isClear() const
    { return myNTriangles == 0; }

    /// Returns an approximation of the signed solid angle of the mesh from the specified query_point
    /// accuracy_scale is the value of (maxP/q) beyond which the approximation of the box will be used.
    T computeSolidAngle(const UT_Vector3T<T> &query_point, const T accuracy_scale = T(2.0)) const;

private:
    struct BoxData;

    static constexpr uint BVH_N = 4;
    UT_BVH<BVH_N> myTree;
    int myNBoxes;
    int myOrder;
    std::unique_ptr<BoxData[]> myData;
    int myNTriangles;
    const int *myTrianglePoints;
    int myNPoints;
    const UT_Vector3T<S> *myPositions;
};

template<typename T>
T UTsignedAngleSegment(
    const UT_Vector2T<T> &a,
    const UT_Vector2T<T> &b,
    const UT_Vector2T<T> &query)
{
    // Make a and b relative to query
    UT_Vector2T<T> qa = a-query;
    UT_Vector2T<T> qb = b-query;

    // If any segment vertices are coincident with query,
    // query is on the segment, which we treat as no angle.
    if (qa.isZero() || qb.isZero())
        return T(0);

    // numerator = |qa||qb|sin(theta)
    const T numerator = cross(qa, qb);

    // If numerator is 0, regardless of denominator, query is on the
    // surface, which we treat as no solid angle.
    if (numerator == 0)
        return T(0);

    // denominator = |qa||qb|cos(theta)
    const T denominator = dot(qa,qb);

    // numerator/denominator = tan(theta)
    return SYSatan2(numerator, denominator);
}

/// Class for quickly approximating signed subtended angle of a large curve
/// from many query points.  This is useful for computing the
/// generalized winding number at many points.
///
/// NOTE: This is currently only instantiated for <float,float>.
template<typename T,typename S>
class UT_SubtendedAngle
{
public:
    /// This is outlined so that we don't need to include UT_BVHImpl.h
    UT_SubtendedAngle();
    /// This is outlined so that we don't need to include UT_BVHImpl.h
    ~UT_SubtendedAngle();

    /// NOTE: This does not take ownership over segment_points or positions,
    ///       but does keep pointers to them, so the caller must keep them in
    ///       scope for the lifetime of this structure.
    UT_SubtendedAngle(
        const int nsegments,
        const int *const segment_points,
        const int npoints,
        const UT_Vector2T<S> *const positions,
        const int order = 2)
        : UT_SubtendedAngle()
    { init(nsegments, segment_points, npoints, positions, order); }

    /// Initialize the tree and data.
    /// NOTE: It is safe to call init on a UT_SolidAngle that has had init
    ///       called on it before, to re-initialize it.
    void init(
        const int nsegments,
        const int *const segment_points,
        const int npoints,
        const UT_Vector2T<S> *const positions,
        const int order = 2);

    /// Frees myTree and myData, and clears the rest.
    void clear();

    /// Returns true if this is clear
    bool isClear() const
    { return myNSegments == 0; }

    /// Returns an approximation of the signed solid angle of the mesh from the specified query_point
    /// accuracy_scale is the value of (maxP/q) beyond which the approximation of the box will be used.
    T computeAngle(const UT_Vector2T<T> &query_point, const T accuracy_scale = T(2.0)) const;

private:
    struct BoxData;

    static constexpr uint BVH_N = 4;
    UT_BVH<BVH_N> myTree;
    int myNBoxes;
    int myOrder;
    std::unique_ptr<BoxData[]> myData;
    int myNSegments;
    const int *mySegmentPoints;
    int myNPoints;
    const UT_Vector2T<S> *myPositions;
};

} // End HDK_Sample namespace
#endif
