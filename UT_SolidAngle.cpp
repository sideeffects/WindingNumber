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

#include "UT_SolidAngle.h"
#include "UT_BVHImpl.h"

#include "UT_SmallArray.h"
#include "UT_FixedVector.h"
#include "VM_SIMD.h"
#include "SYS_Types.h"
#include <type_traits>
#include <utility>

#define SOLID_ANGLE_TIME_PRECOMPUTE 0

#if SOLID_ANGLE_TIME_PRECOMPUTE
#include <UT/UT_StopWatch.h>
#endif

#define SOLID_ANGLE_DEBUG 0
#if SOLID_ANGLE_DEBUG
#include <UT/UT_Debug.h>
#endif

#define TAYLOR_SERIES_ORDER 2

namespace HDK_Sample {

template<typename T,typename S>
struct UT_SolidAngle<T,S>::BoxData
{
    void clear()
    {
        // Set everything to zero
        memset(this,0,sizeof(*this));
    }

    using Type  = typename std::conditional<BVH_N==4 && std::is_same<T,float>::value, v4uf, UT_FixedVector<T,BVH_N>>::type;
    using SType = typename std::conditional<BVH_N==4 && std::is_same<S,float>::value, v4uf, UT_FixedVector<S,BVH_N>>::type;

    /// An upper bound on the squared distance from myAverageP to the farthest point in the box.
    SType myMaxPDist2;

    /// Centre of mass of the mesh surface in this box
    UT_FixedVector<Type,3> myAverageP;

    /// Unnormalized, area-weighted normal of the mesh in this box
    UT_FixedVector<Type,3> myN;

#if TAYLOR_SERIES_ORDER >= 1
    /// Values for Omega_1
    /// @{
    UT_FixedVector<Type,3> myNijDiag;  // Nxx, Nyy, Nzz
    Type myNxy_Nyx;               // Nxy+Nyx
    Type myNyz_Nzy;               // Nyz+Nzy
    Type myNzx_Nxz;               // Nzx+Nxz
    /// @}
#endif

#if TAYLOR_SERIES_ORDER >= 2
    /// Values for Omega_2
    /// @{
    UT_FixedVector<Type,3> myNijkDiag; // Nxxx, Nyyy, Nzzz
    Type mySumPermuteNxyz;        // (Nxyz+Nxzy+Nyzx+Nyxz+Nzxy+Nzyx) = 2*(Nxyz+Nyzx+Nzxy)
    Type my2Nxxy_Nyxx; // Nxxy+Nxyx+Nyxx = 2Nxxy+Nyxx
    Type my2Nxxz_Nzxx; // Nxxz+Nxzx+Nzxx = 2Nxxz+Nzxx
    Type my2Nyyz_Nzyy; // Nyyz+Nyzy+Nzyy = 2Nyyz+Nzyy
    Type my2Nyyx_Nxyy; // Nyyx+Nyxy+Nxyy = 2Nyyx+Nxyy
    Type my2Nzzx_Nxzz; // Nzzx+Nzxz+Nxzz = 2Nzzx+Nxzz
    Type my2Nzzy_Nyzz; // Nzzy+Nzyz+Nyzz = 2Nzzy+Nyzz
    /// @}
#endif
};

template<typename T,typename S>
UT_SolidAngle<T,S>::UT_SolidAngle()
    : myTree()
    , myNBoxes(0)
    , myOrder(2)
    , myData(nullptr)
    , myNTriangles(0)
    , myTrianglePoints(nullptr)
    , myNPoints(0)
    , myPositions(nullptr)
{}

template<typename T,typename S>
UT_SolidAngle<T,S>::~UT_SolidAngle()
{
    // Default destruction works, but this needs to be outlined
    // to avoid having to include UT_BVHImpl.h in the header,
    // (for the UT_UniquePtr destructor.)
}

template<typename T,typename S>
void UT_SolidAngle<T,S>::init(
    const int ntriangles,
    const int *const triangle_points,
    const int npoints,
    const UT_Vector3T<S> *const positions,
    const int order)
{
#if SOLID_ANGLE_DEBUG
    UTdebugFormat("");
    UTdebugFormat("");
    UTdebugFormat("Building BVH for {} ntriangles on {} points:", ntriangles, npoints);
#endif
    myOrder = order;
    myNTriangles = ntriangles;
    myTrianglePoints = triangle_points;
    myNPoints = npoints;
    myPositions = positions;

#if SOLID_ANGLE_TIME_PRECOMPUTE
    UT_StopWatch timer;
    timer.start();
#endif
    UT_SmallArray<UT::Box<S,3>> triangle_boxes;
    triangle_boxes.setSizeNoInit(ntriangles);
    if (ntriangles < 16*1024)
    {
        const int *cur_triangle_points = triangle_points;
        for (int i = 0; i < ntriangles; ++i, cur_triangle_points += 3)
        {
            UT::Box<S,3> &box = triangle_boxes[i];
            box.initBounds(positions[cur_triangle_points[0]]);
            box.enlargeBounds(positions[cur_triangle_points[1]]);
            box.enlargeBounds(positions[cur_triangle_points[2]]);
        }
    }
    else
    {
        UTparallelFor(UT_BlockedRange<int>(0,ntriangles), [triangle_points,&triangle_boxes,positions](const UT_BlockedRange<int> &r)
        {
            const int *cur_triangle_points = triangle_points + exint(r.begin())*3;
            for (int i = r.begin(), end = r.end(); i < end; ++i, cur_triangle_points += 3)
            {
                UT::Box<S,3> &box = triangle_boxes[i];
                box.initBounds(positions[cur_triangle_points[0]]);
                box.enlargeBounds(positions[cur_triangle_points[1]]);
                box.enlargeBounds(positions[cur_triangle_points[2]]);
            }
        });
    }
#if SOLID_ANGLE_TIME_PRECOMPUTE
    double time = timer.stop();
    UTdebugFormat("{} s to create bounding boxes.", time);
    timer.start();
#endif
    myTree.template init<UT::BVH_Heuristic::BOX_AREA,S,3>(triangle_boxes.array(), ntriangles);
#if SOLID_ANGLE_TIME_PRECOMPUTE
    time = timer.stop();
    UTdebugFormat("{} s to initialize UT_BVH structure.  {} nodes", time, myTree.getNumNodes());
#endif

    //myTree.debugDump();

    const int nnodes = myTree.getNumNodes();

    myNBoxes = nnodes;
    BoxData *box_data = new BoxData[nnodes];
    myData.reset(box_data);

    // Some data are only needed during initialization.
    struct LocalData
    {
        // Bounding box
        UT::Box<S,3> myBox;

        // P and N are needed from each child for computing Nij.
        UT_Vector3T<T> myAverageP;
        UT_Vector3T<T> myAreaP;
        UT_Vector3T<T> myN;

        // Unsigned area is needed for computing the average position.
        T myArea;

#if TAYLOR_SERIES_ORDER >= 1
        // These are needed for computing Nijk.
        UT_Vector3T<T> myNijDiag;
        T myNxy; T myNyx;
        T myNyz; T myNzy;
        T myNzx; T myNxz;
#endif

#if TAYLOR_SERIES_ORDER >= 2
        UT_Vector3T<T> myNijkDiag; // Nxxx, Nyyy, Nzzz
        T mySumPermuteNxyz; // (Nxyz+Nxzy+Nyzx+Nyxz+Nzxy+Nzyx) = 2*(Nxyz+Nyzx+Nzxy)
        T my2Nxxy_Nyxx;     // Nxxy+Nxyx+Nyxx = 2Nxxy+Nyxx
        T my2Nxxz_Nzxx;     // Nxxz+Nxzx+Nzxx = 2Nxxz+Nzxx
        T my2Nyyz_Nzyy;     // Nyyz+Nyzy+Nzyy = 2Nyyz+Nzyy
        T my2Nyyx_Nxyy;     // Nyyx+Nyxy+Nxyy = 2Nyyx+Nxyy
        T my2Nzzx_Nxzz;     // Nzzx+Nzxz+Nxzz = 2Nzzx+Nxzz
        T my2Nzzy_Nyzz;     // Nzzy+Nzyz+Nyzz = 2Nzzy+Nyzz
#endif
    };

    struct PrecomputeFunctors
    {
        BoxData *const myBoxData;
        const UT::Box<S,3> *const myTriangleBoxes;
        const int *const myTrianglePoints;
        const UT_Vector3T<S> *const myPositions;
        const int myOrder;

        PrecomputeFunctors(
            BoxData *box_data,
            const UT::Box<S,3> *triangle_boxes,
            const int *triangle_points,
            const UT_Vector3T<S> *positions,
            const int order)
            : myBoxData(box_data)
            , myTriangleBoxes(triangle_boxes)
            , myTrianglePoints(triangle_points)
            , myPositions(positions)
            , myOrder(order)
        {}
        constexpr SYS_FORCE_INLINE bool pre(const int nodei, LocalData *data_for_parent) const
        {
            return true;
        }
        void item(const int itemi, const int parent_nodei, LocalData &data_for_parent) const
        {
            const UT_Vector3T<S> *const positions = myPositions;
            const int *const cur_triangle_points = myTrianglePoints + 3*itemi;
            const UT_Vector3T<T> a = positions[cur_triangle_points[0]];
            const UT_Vector3T<T> b = positions[cur_triangle_points[1]];
            const UT_Vector3T<T> c = positions[cur_triangle_points[2]];
            const UT_Vector3T<T> ab = b-a;
            const UT_Vector3T<T> ac = c-a;

            const UT::Box<S,3> &triangle_box = myTriangleBoxes[itemi];
            data_for_parent.myBox.initBounds(triangle_box.getMin(), triangle_box.getMax());

            // Area-weighted normal (unnormalized)
            const UT_Vector3T<T> N = T(0.5)*cross(ab,ac);
            const T area2 = N.length2();
            const T area = SYSsqrt(area2);
            const UT_Vector3T<T> P = (a+b+c)/3;
            data_for_parent.myAverageP = P;
            data_for_parent.myAreaP = P*area;
            data_for_parent.myN = N;
#if SOLID_ANGLE_DEBUG
            UTdebugFormat("");
            UTdebugFormat("Triangle {}: P = {}; N = {}; area = {}", itemi, P, N, area);
            UTdebugFormat("             box = {}", data_for_parent.myBox);
#endif

            data_for_parent.myArea = area;
#if TAYLOR_SERIES_ORDER >= 1
            const int order = myOrder;
            if (order < 1)
                return;

            // NOTE: Due to P being at the centroid, triangles have Nij = 0
            //       contributions to Nij.
            data_for_parent.myNijDiag = T(0);
            data_for_parent.myNxy = 0; data_for_parent.myNyx = 0;
            data_for_parent.myNyz = 0; data_for_parent.myNzy = 0;
            data_for_parent.myNzx = 0; data_for_parent.myNxz = 0;
#endif

#if TAYLOR_SERIES_ORDER >= 2
            if (order < 2)
                return;

            // If it's zero-length, the results are zero, so we can skip.
            if (area == 0)
            {
                data_for_parent.myNijkDiag = T(0);
                data_for_parent.mySumPermuteNxyz = 0;
                data_for_parent.my2Nxxy_Nyxx = 0;
                data_for_parent.my2Nxxz_Nzxx = 0;
                data_for_parent.my2Nyyz_Nzyy = 0;
                data_for_parent.my2Nyyx_Nxyy = 0;
                data_for_parent.my2Nzzx_Nxzz = 0;
                data_for_parent.my2Nzzy_Nyzz = 0;
                return;
            }

            // We need to use the NORMALIZED normal to multiply the integrals by.
            UT_Vector3T<T> n = N/area;

            // Figure out the order of a, b, and c in x, y, and z
            // for use in computing the integrals for Nijk.
            UT_Vector3T<T> values[3] = {a, b, c};

            int order_x[3] = {0,1,2};
            if (a[0] > b[0])
                std::swap(order_x[0],order_x[1]);
            if (values[order_x[0]][0] > c[0])
                std::swap(order_x[0],order_x[2]);
            if (values[order_x[1]][0] > values[order_x[2]][0])
                std::swap(order_x[1],order_x[2]);
            T dx = values[order_x[2]][0] - values[order_x[0]][0];

            int order_y[3] = {0,1,2};
            if (a[1] > b[1])
                std::swap(order_y[0],order_y[1]);
            if (values[order_y[0]][1] > c[1])
                std::swap(order_y[0],order_y[2]);
            if (values[order_y[1]][1] > values[order_y[2]][1])
                std::swap(order_y[1],order_y[2]);
            T dy = values[order_y[2]][1] - values[order_y[0]][1];

            int order_z[3] = {0,1,2};
            if (a[2] > b[2])
                std::swap(order_z[0],order_z[1]);
            if (values[order_z[0]][2] > c[2])
                std::swap(order_z[0],order_z[2]);
            if (values[order_z[1]][2] > values[order_z[2]][2])
                std::swap(order_z[1],order_z[2]);
            T dz = values[order_z[2]][2] - values[order_z[0]][2];

            auto &&compute_integrals = [](
                const UT_Vector3T<T> &a,
                const UT_Vector3T<T> &b,
                const UT_Vector3T<T> &c,
                const UT_Vector3T<T> &P,
                T *integral_ii,
                T *integral_ij,
                T *integral_ik,
                const int i)
            {
#if SOLID_ANGLE_DEBUG
                UTdebugFormat("             Splitting on {}; a = {}; b = {}; c = {}", char('x'+i), a, b, c);
#endif
                // NOTE: a, b, and c must be in order of the i axis.
                // We're splitting the triangle at the middle i coordinate.
                const UT_Vector3T<T> oab = b - a;
                const UT_Vector3T<T> oac = c - a;
                const UT_Vector3T<T> ocb = b - c;
                UT_ASSERT_MSG_P(oac[i] > 0, "This should have been checked by the caller.");
                const T t = oab[i]/oac[i];
                UT_ASSERT_MSG_P(t >= 0 && t <= 1, "Either sorting must have gone wrong, or there are input NaNs.");

                const int j = (i==2) ? 0 : (i+1);
                const int k = (j==2) ? 0 : (j+1);
                const T jdiff = t*oac[j] - oab[j];
                const T kdiff = t*oac[k] - oab[k];
                UT_Vector3T<T> cross_a;
                cross_a[0] = (jdiff*oab[k] - kdiff*oab[j]);
                cross_a[1] = kdiff*oab[i];
                cross_a[2] = jdiff*oab[i];
                UT_Vector3T<T> cross_c;
                cross_c[0] = (jdiff*ocb[k] - kdiff*ocb[j]);
                cross_c[1] = kdiff*ocb[i];
                cross_c[2] = jdiff*ocb[i];
                const T area_scale_a = cross_a.length();
                const T area_scale_c = cross_c.length();
                const T Pai = a[i] - P[i];
                const T Pci = c[i] - P[i];

                // Integral over the area of the triangle of (pi^2)dA,
                // by splitting the triangle into two at b, the a side
                // and the c side.
                const T int_ii_a = area_scale_a*(T(0.5)*Pai*Pai + T(2.0/3.0)*Pai*oab[i] + T(0.25)*oab[i]*oab[i]);
                const T int_ii_c = area_scale_c*(T(0.5)*Pci*Pci + T(2.0/3.0)*Pci*ocb[i] + T(0.25)*ocb[i]*ocb[i]);
                *integral_ii = int_ii_a + int_ii_c;
#if SOLID_ANGLE_DEBUG
                UTdebugFormat("             integral_{}{}_a = {}; integral_{}{}_c = {}", char('x'+i), char('x'+i), int_ii_a, char('x'+i), char('x'+i), int_ii_c);
#endif

                int jk = j;
                T *integral = integral_ij;
                T diff = jdiff;
                while (true) // This only does 2 iterations, one for j and one for k
                {
                    if (integral)
                    {
                        T obmidj = b[jk] + T(0.5)*diff;
                        T oabmidj = obmidj - a[jk];
                        T ocbmidj = obmidj - c[jk];
                        T Paj = a[jk] - P[jk];
                        T Pcj = c[jk] - P[jk];
                        // Integral over the area of the triangle of (pi*pj)dA
                        const T int_ij_a = area_scale_a*(T(0.5)*Pai*Paj + T(1.0/3.0)*Pai*oabmidj + T(1.0/3.0)*Paj*oab[i] + T(0.25)*oab[i]*oabmidj);
                        const T int_ij_c = area_scale_c*(T(0.5)*Pci*Pcj + T(1.0/3.0)*Pci*ocbmidj + T(1.0/3.0)*Pcj*ocb[i] + T(0.25)*ocb[i]*ocbmidj);
                        *integral = int_ij_a + int_ij_c;
#if SOLID_ANGLE_DEBUG
                        UTdebugFormat("             integral_{}{}_a = {}; integral_{}{}_c = {}", char('x'+i), char('x'+jk), int_ij_a, char('x'+i), char('x'+jk), int_ij_c);
#endif
                    }
                    if (jk == k)
                        break;
                    jk = k;
                    integral = integral_ik;
                    diff = kdiff;
                }
            };

            T integral_xx = 0;
            T integral_xy = 0;
            T integral_yy = 0;
            T integral_yz = 0;
            T integral_zz = 0;
            T integral_zx = 0;
            // Note that if the span of any axis is zero, the integral must be zero,
            // since there's a factor of (p_i-P_i), i.e. value minus average,
            // and every value must be equal to the average, giving zero.
            if (dx > 0)
            {
                compute_integrals(
                    values[order_x[0]], values[order_x[1]], values[order_x[2]], P,
                    &integral_xx, ((dx >= dy && dy > 0) ? &integral_xy : nullptr), ((dx >= dz && dz > 0) ? &integral_zx : nullptr), 0);
            }
            if (dy > 0)
            {
                compute_integrals(
                    values[order_y[0]], values[order_y[1]], values[order_y[2]], P,
                    &integral_yy, ((dy >= dz && dz > 0) ? &integral_yz : nullptr), ((dx < dy && dx > 0) ? &integral_xy : nullptr), 1);
            }
            if (dz > 0)
            {
                compute_integrals(
                    values[order_z[0]], values[order_z[1]], values[order_z[2]], P,
                    &integral_zz, ((dx < dz && dx > 0) ? &integral_zx : nullptr), ((dy < dz && dy > 0) ? &integral_yz : nullptr), 2);
            }

            UT_Vector3T<T> Niii;
            Niii[0] = integral_xx;
            Niii[1] = integral_yy;
            Niii[2] = integral_zz;
            Niii *= n;
            data_for_parent.myNijkDiag = Niii;
            data_for_parent.mySumPermuteNxyz = 2*(n[0]*integral_yz + n[1]*integral_zx + n[2]*integral_xy);
            T Nxxy = n[0]*integral_xy;
            T Nxxz = n[0]*integral_zx;
            T Nyyz = n[1]*integral_yz;
            T Nyyx = n[1]*integral_xy;
            T Nzzx = n[2]*integral_zx;
            T Nzzy = n[2]*integral_yz;
            data_for_parent.my2Nxxy_Nyxx = 2*Nxxy + n[1]*integral_xx;
            data_for_parent.my2Nxxz_Nzxx = 2*Nxxz + n[2]*integral_xx;
            data_for_parent.my2Nyyz_Nzyy = 2*Nyyz + n[2]*integral_yy;
            data_for_parent.my2Nyyx_Nxyy = 2*Nyyx + n[0]*integral_yy;
            data_for_parent.my2Nzzx_Nxzz = 2*Nzzx + n[0]*integral_zz;
            data_for_parent.my2Nzzy_Nyzz = 2*Nzzy + n[1]*integral_zz;
#if SOLID_ANGLE_DEBUG
            UTdebugFormat("             integral_xx = {}; yy = {}; zz = {}", integral_xx, integral_yy, integral_zz);
            UTdebugFormat("             integral_xy = {}; yz = {}; zx = {}", integral_xy, integral_yz, integral_zx);
#endif
#endif
        }

        void post(const int nodei, const int parent_nodei, LocalData *data_for_parent, const int nchildren, const LocalData *child_data_array) const
        {
            // NOTE: Although in the general case, data_for_parent may be null for the root call,
            //       this functor assumes that it's non-null, so the call below must pass a non-null pointer.

            BoxData &current_box_data = myBoxData[nodei];

            UT_Vector3T<T> N = child_data_array[0].myN;
            ((T*)&current_box_data.myN[0])[0] = N[0];
            ((T*)&current_box_data.myN[1])[0] = N[1];
            ((T*)&current_box_data.myN[2])[0] = N[2];
            UT_Vector3T<T> areaP = child_data_array[0].myAreaP;
            T area = child_data_array[0].myArea;
            UT_Vector3T<T> local_P = child_data_array[0].myAverageP;
            ((T*)&current_box_data.myAverageP[0])[0] = local_P[0];
            ((T*)&current_box_data.myAverageP[1])[0] = local_P[1];
            ((T*)&current_box_data.myAverageP[2])[0] = local_P[2];
            for (int i = 1; i < nchildren; ++i)
            {
                const UT_Vector3T<T> local_N = child_data_array[i].myN;
                N += local_N;
                ((T*)&current_box_data.myN[0])[i] = local_N[0];
                ((T*)&current_box_data.myN[1])[i] = local_N[1];
                ((T*)&current_box_data.myN[2])[i] = local_N[2];
                areaP += child_data_array[i].myAreaP;
                area += child_data_array[i].myArea;
                const UT_Vector3T<T> local_P = child_data_array[i].myAverageP;
                ((T*)&current_box_data.myAverageP[0])[i] = local_P[0];
                ((T*)&current_box_data.myAverageP[1])[i] = local_P[1];
                ((T*)&current_box_data.myAverageP[2])[i] = local_P[2];
            }
            for (int i = nchildren; i < BVH_N; ++i)
            {
                // Set to zero, just to avoid false positives for uses of uninitialized memory.
                ((T*)&current_box_data.myN[0])[i] = 0;
                ((T*)&current_box_data.myN[1])[i] = 0;
                ((T*)&current_box_data.myN[2])[i] = 0;
                ((T*)&current_box_data.myAverageP[0])[i] = 0;
                ((T*)&current_box_data.myAverageP[1])[i] = 0;
                ((T*)&current_box_data.myAverageP[2])[i] = 0;
            }
            data_for_parent->myN = N;
            data_for_parent->myAreaP = areaP;
            data_for_parent->myArea = area;

            UT::Box<S,3> box(child_data_array[0].myBox);
            for (int i = 1; i < nchildren; ++i)
                box.enlargeBounds(child_data_array[i].myBox);

            // Normalize P
            UT_Vector3T<T> averageP;
            if (area > 0)
                averageP = areaP/area;
            else
                averageP = T(0.5)*(box.getMin() + box.getMax());
            data_for_parent->myAverageP = averageP;

            data_for_parent->myBox = box;

            for (int i = 0; i < nchildren; ++i)
            {
                const UT::Box<S,3> &local_box(child_data_array[i].myBox);
                const UT_Vector3T<T> &local_P = child_data_array[i].myAverageP;
                const UT_Vector3T<T> maxPDiff = SYSmax(local_P-UT_Vector3T<T>(local_box.getMin()), UT_Vector3T<T>(local_box.getMax())-local_P);
                ((T*)&current_box_data.myMaxPDist2)[i] = maxPDiff.length2();
            }
            for (int i = nchildren; i < BVH_N; ++i)
            {
                // This child is non-existent.  If we set myMaxPDist2 to infinity, it will never
                // use the approximation, and the traverseVector function can check for EMPTY.
                ((T*)&current_box_data.myMaxPDist2)[i] = std::numeric_limits<T>::infinity();
            }

#if TAYLOR_SERIES_ORDER >= 1
            const int order = myOrder;
            if (order >= 1)
            {
                // We now have the current box's P, so we can adjust Nij and Nijk
                data_for_parent->myNijDiag = child_data_array[0].myNijDiag;
                data_for_parent->myNxy = 0;
                data_for_parent->myNyx = 0;
                data_for_parent->myNyz = 0;
                data_for_parent->myNzy = 0;
                data_for_parent->myNzx = 0;
                data_for_parent->myNxz = 0;
#if TAYLOR_SERIES_ORDER >= 2
                data_for_parent->myNijkDiag = child_data_array[0].myNijkDiag;
                data_for_parent->mySumPermuteNxyz = child_data_array[0].mySumPermuteNxyz;
                data_for_parent->my2Nxxy_Nyxx = child_data_array[0].my2Nxxy_Nyxx;
                data_for_parent->my2Nxxz_Nzxx = child_data_array[0].my2Nxxz_Nzxx;
                data_for_parent->my2Nyyz_Nzyy = child_data_array[0].my2Nyyz_Nzyy;
                data_for_parent->my2Nyyx_Nxyy = child_data_array[0].my2Nyyx_Nxyy;
                data_for_parent->my2Nzzx_Nxzz = child_data_array[0].my2Nzzx_Nxzz;
                data_for_parent->my2Nzzy_Nyzz = child_data_array[0].my2Nzzy_Nyzz;
#endif

                for (int i = 1; i < nchildren; ++i)
                {
                    data_for_parent->myNijDiag += child_data_array[i].myNijDiag;
#if TAYLOR_SERIES_ORDER >= 2
                    data_for_parent->myNijkDiag += child_data_array[i].myNijkDiag;
                    data_for_parent->mySumPermuteNxyz += child_data_array[i].mySumPermuteNxyz;
                    data_for_parent->my2Nxxy_Nyxx += child_data_array[i].my2Nxxy_Nyxx;
                    data_for_parent->my2Nxxz_Nzxx += child_data_array[i].my2Nxxz_Nzxx;
                    data_for_parent->my2Nyyz_Nzyy += child_data_array[i].my2Nyyz_Nzyy;
                    data_for_parent->my2Nyyx_Nxyy += child_data_array[i].my2Nyyx_Nxyy;
                    data_for_parent->my2Nzzx_Nxzz += child_data_array[i].my2Nzzx_Nxzz;
                    data_for_parent->my2Nzzy_Nyzz += child_data_array[i].my2Nzzy_Nyzz;
#endif
                }
                for (int j = 0; j < 3; ++j)
                    ((T*)&current_box_data.myNijDiag[j])[0] = child_data_array[0].myNijDiag[j];
                ((T*)&current_box_data.myNxy_Nyx)[0] = child_data_array[0].myNxy + child_data_array[0].myNyx;
                ((T*)&current_box_data.myNyz_Nzy)[0] = child_data_array[0].myNyz + child_data_array[0].myNzy;
                ((T*)&current_box_data.myNzx_Nxz)[0] = child_data_array[0].myNzx + child_data_array[0].myNxz;
                for (int j = 0; j < 3; ++j)
                    ((T*)&current_box_data.myNijkDiag[j])[0] = child_data_array[0].myNijkDiag[j];
                ((T*)&current_box_data.mySumPermuteNxyz)[0] = child_data_array[0].mySumPermuteNxyz;
                ((T*)&current_box_data.my2Nxxy_Nyxx)[0] = child_data_array[0].my2Nxxy_Nyxx;
                ((T*)&current_box_data.my2Nxxz_Nzxx)[0] = child_data_array[0].my2Nxxz_Nzxx;
                ((T*)&current_box_data.my2Nyyz_Nzyy)[0] = child_data_array[0].my2Nyyz_Nzyy;
                ((T*)&current_box_data.my2Nyyx_Nxyy)[0] = child_data_array[0].my2Nyyx_Nxyy;
                ((T*)&current_box_data.my2Nzzx_Nxzz)[0] = child_data_array[0].my2Nzzx_Nxzz;
                ((T*)&current_box_data.my2Nzzy_Nyzz)[0] = child_data_array[0].my2Nzzy_Nyzz;
                for (int i = 1; i < nchildren; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                        ((T*)&current_box_data.myNijDiag[j])[i] = child_data_array[i].myNijDiag[j];
                    ((T*)&current_box_data.myNxy_Nyx)[i] = child_data_array[i].myNxy + child_data_array[i].myNyx;
                    ((T*)&current_box_data.myNyz_Nzy)[i] = child_data_array[i].myNyz + child_data_array[i].myNzy;
                    ((T*)&current_box_data.myNzx_Nxz)[i] = child_data_array[i].myNzx + child_data_array[i].myNxz;
                    for (int j = 0; j < 3; ++j)
                        ((T*)&current_box_data.myNijkDiag[j])[i] = child_data_array[i].myNijkDiag[j];
                    ((T*)&current_box_data.mySumPermuteNxyz)[i] = child_data_array[i].mySumPermuteNxyz;
                    ((T*)&current_box_data.my2Nxxy_Nyxx)[i] = child_data_array[i].my2Nxxy_Nyxx;
                    ((T*)&current_box_data.my2Nxxz_Nzxx)[i] = child_data_array[i].my2Nxxz_Nzxx;
                    ((T*)&current_box_data.my2Nyyz_Nzyy)[i] = child_data_array[i].my2Nyyz_Nzyy;
                    ((T*)&current_box_data.my2Nyyx_Nxyy)[i] = child_data_array[i].my2Nyyx_Nxyy;
                    ((T*)&current_box_data.my2Nzzx_Nxzz)[i] = child_data_array[i].my2Nzzx_Nxzz;
                    ((T*)&current_box_data.my2Nzzy_Nyzz)[i] = child_data_array[i].my2Nzzy_Nyzz;
                }
                for (int i = nchildren; i < BVH_N; ++i)
                {
                    // Set to zero, just to avoid false positives for uses of uninitialized memory.
                    for (int j = 0; j < 3; ++j)
                        ((T*)&current_box_data.myNijDiag[j])[i] = 0;
                    ((T*)&current_box_data.myNxy_Nyx)[i] = 0;
                    ((T*)&current_box_data.myNyz_Nzy)[i] = 0;
                    ((T*)&current_box_data.myNzx_Nxz)[i] = 0;
                    for (int j = 0; j < 3; ++j)
                        ((T*)&current_box_data.myNijkDiag[j])[i] = 0;
                    ((T*)&current_box_data.mySumPermuteNxyz)[i] = 0;
                    ((T*)&current_box_data.my2Nxxy_Nyxx)[i] = 0;
                    ((T*)&current_box_data.my2Nxxz_Nzxx)[i] = 0;
                    ((T*)&current_box_data.my2Nyyz_Nzyy)[i] = 0;
                    ((T*)&current_box_data.my2Nyyx_Nxyy)[i] = 0;
                    ((T*)&current_box_data.my2Nzzx_Nxzz)[i] = 0;
                    ((T*)&current_box_data.my2Nzzy_Nyzz)[i] = 0;
                }

                for (int i = 0; i < nchildren; ++i)
                {
                    const LocalData &child_data = child_data_array[i];
                    UT_Vector3T<T> displacement = child_data.myAverageP - UT_Vector3T<T>(data_for_parent->myAverageP);
                    UT_Vector3T<T> N = child_data.myN;

                    // Adjust Nij for the change in centre P
                    data_for_parent->myNijDiag += N*displacement;
                    T Nxy = child_data.myNxy + N[0]*displacement[1];
                    T Nyx = child_data.myNyx + N[1]*displacement[0];
                    T Nyz = child_data.myNyz + N[1]*displacement[2];
                    T Nzy = child_data.myNzy + N[2]*displacement[1];
                    T Nzx = child_data.myNzx + N[2]*displacement[0];
                    T Nxz = child_data.myNxz + N[0]*displacement[2];

                    data_for_parent->myNxy += Nxy;
                    data_for_parent->myNyx += Nyx;
                    data_for_parent->myNyz += Nyz;
                    data_for_parent->myNzy += Nzy;
                    data_for_parent->myNzx += Nzx;
                    data_for_parent->myNxz += Nxz;

#if TAYLOR_SERIES_ORDER >= 2
                    if (order >= 2)
                    {
                        // Adjust Nijk for the change in centre P
                        data_for_parent->myNijkDiag += T(2)*displacement*child_data.myNijDiag + displacement*displacement*child_data.myN;
                        data_for_parent->mySumPermuteNxyz += (displacement[0]*(Nyz+Nzy) + displacement[1]*(Nzx+Nxz) + displacement[2]*(Nxy+Nyx));
                        data_for_parent->my2Nxxy_Nyxx +=
                            2*(displacement[1]*child_data.myNijDiag[0] + displacement[0]*child_data.myNxy + N[0]*displacement[0]*displacement[1])
                            + 2*child_data.myNyx*displacement[0] + N[1]*displacement[0]*displacement[0];
                        data_for_parent->my2Nxxz_Nzxx +=
                            2*(displacement[2]*child_data.myNijDiag[0] + displacement[0]*child_data.myNxz + N[0]*displacement[0]*displacement[2])
                            + 2*child_data.myNzx*displacement[0] + N[2]*displacement[0]*displacement[0];
                        data_for_parent->my2Nyyz_Nzyy +=
                            2*(displacement[2]*child_data.myNijDiag[1] + displacement[1]*child_data.myNyz + N[1]*displacement[1]*displacement[2])
                            + 2*child_data.myNzy*displacement[1] + N[2]*displacement[1]*displacement[1];
                        data_for_parent->my2Nyyx_Nxyy +=
                            2*(displacement[0]*child_data.myNijDiag[1] + displacement[1]*child_data.myNyx + N[1]*displacement[1]*displacement[0])
                            + 2*child_data.myNxy*displacement[1] + N[0]*displacement[1]*displacement[1];
                        data_for_parent->my2Nzzx_Nxzz +=
                            2*(displacement[0]*child_data.myNijDiag[2] + displacement[2]*child_data.myNzx + N[2]*displacement[2]*displacement[0])
                            + 2*child_data.myNxz*displacement[2] + N[0]*displacement[2]*displacement[2];
                        data_for_parent->my2Nzzy_Nyzz +=
                            2*(displacement[1]*child_data.myNijDiag[2] + displacement[2]*child_data.myNzy + N[2]*displacement[2]*displacement[1])
                            + 2*child_data.myNyz*displacement[2] + N[1]*displacement[2]*displacement[2];
                    }
#endif
                }
            }
#endif
#if SOLID_ANGLE_DEBUG
            UTdebugFormat("");
            UTdebugFormat("Node {}: nchildren = {}; maxP = {}", nodei, nchildren, SYSsqrt(current_box_data.myMaxPDist2));
            UTdebugFormat("         P = {}; N = {}", current_box_data.myAverageP, current_box_data.myN);
#if TAYLOR_SERIES_ORDER >= 1
            UTdebugFormat("         Nii = {}", current_box_data.myNijDiag);
            UTdebugFormat("         Nxy+Nyx = {}; Nyz+Nzy = {}; Nyz+Nzy = {}", current_box_data.myNxy_Nyx, current_box_data.myNyz_Nzy, current_box_data.myNzx_Nxz);
#if TAYLOR_SERIES_ORDER >= 2
            UTdebugFormat("         Niii = {}; 2(Nxyz+Nyzx+Nzxy) = {}", current_box_data.myNijkDiag, current_box_data.mySumPermuteNxyz);
            UTdebugFormat("         2Nxxy+Nyxx = {}; 2Nxxz+Nzxx = {}", current_box_data.my2Nxxy_Nyxx, current_box_data.my2Nxxz_Nzxx);
            UTdebugFormat("         2Nyyz+Nzyy = {}; 2Nyyx+Nxyy = {}", current_box_data.my2Nyyz_Nzyy, current_box_data.my2Nyyx_Nxyy);
            UTdebugFormat("         2Nzzx+Nxzz = {}; 2Nzzy+Nyzz = {}", current_box_data.my2Nzzx_Nxzz, current_box_data.my2Nzzy_Nyzz);
#endif
#endif
#endif
        }
    };

#if SOLID_ANGLE_TIME_PRECOMPUTE
    timer.start();
#endif
    const PrecomputeFunctors functors(box_data, triangle_boxes.array(), triangle_points, positions, order);
    // NOTE: post-functor relies on non-null data_for_parent, so we have to pass one.
    LocalData local_data;
    myTree.template traverseParallel<LocalData>(4096, functors, &local_data);
    //myTree.template traverse<LocalData>(functors);
#if SOLID_ANGLE_TIME_PRECOMPUTE
    time = timer.stop();
    UTdebugFormat("{} s to precompute coefficients.", time);
#endif
}

template<typename T,typename S>
void UT_SolidAngle<T, S>::clear()
{
    myTree.clear();
    myNBoxes = 0;
    myOrder = 2;
    myData.reset();
    myNTriangles = 0;
    myTrianglePoints = nullptr;
    myNPoints = 0;
    myPositions = nullptr;
}

template<typename T,typename S>
T UT_SolidAngle<T, S>::computeSolidAngle(const UT_Vector3T<T> &query_point, const T accuracy_scale) const
{
    const T accuracy_scale2 = accuracy_scale*accuracy_scale;

    struct SolidAngleFunctors
    {
        const BoxData *const myBoxData;
        const UT_Vector3T<T> myQueryPoint;
        const T myAccuracyScale2;
        const UT_Vector3T<S> *const myPositions;
        const int *const myTrianglePoints;
        const int myOrder;

        SolidAngleFunctors(
            const BoxData *const box_data,
            const UT_Vector3T<T> &query_point,
            const T accuracy_scale2,
            const int order,
            const UT_Vector3T<S> *const positions,
            const int *const triangle_points)
            : myBoxData(box_data)
            , myQueryPoint(query_point)
            , myAccuracyScale2(accuracy_scale2)
            , myOrder(order)
            , myPositions(positions)
            , myTrianglePoints(triangle_points)
        {}
        uint pre(const int nodei, T *data_for_parent) const
        {
            const BoxData &data = myBoxData[nodei];
            const typename BoxData::Type maxP2 = data.myMaxPDist2;
            UT_FixedVector<typename BoxData::Type,3> q;
            q[0] = typename BoxData::Type(myQueryPoint[0]);
            q[1] = typename BoxData::Type(myQueryPoint[1]);
            q[2] = typename BoxData::Type(myQueryPoint[2]);
            q -= data.myAverageP;
            const typename BoxData::Type qlength2 = q[0]*q[0] + q[1]*q[1] + q[2]*q[2];

            // If the query point is within a factor of accuracy_scale of the box radius,
            // it's assumed to be not a good enough approximation, so it needs to descend.
            // TODO: Is there a way to estimate the error?
            static_assert((std::is_same<typename BoxData::Type,v4uf>::value), "FIXME: Implement support for other tuple types!");
            v4uu descend_mask = (qlength2 <= maxP2*myAccuracyScale2);
            uint descend_bitmask = _mm_movemask_ps(V4SF(descend_mask.vector));
            constexpr uint allchildbits = ((uint(1)<<BVH_N)-1);
            if (descend_bitmask == allchildbits)
            {
                *data_for_parent = 0;
                return allchildbits;
            }

            // qlength2 must be non-zero, since it's strictly greater than something.
            // We still need to be careful for NaNs, though, because the 4th power might cause problems.
            const typename BoxData::Type qlength_m2 = typename BoxData::Type(1.0)/qlength2;
            const typename BoxData::Type qlength_m1 = sqrt(qlength_m2);

            // Normalize q to reduce issues with overflow/underflow, since we'd need the 7th power
            // if we didn't normalize, and (1e-6)^-7 = 1e42, which overflows single-precision.
            q *= qlength_m1;

            typename BoxData::Type Omega_approx = -qlength_m2*dot(q,data.myN);
#if TAYLOR_SERIES_ORDER >= 1
            const int order = myOrder;
            if (order >= 1)
            {
                const UT_FixedVector<typename BoxData::Type,3> q2 = q*q;
                const typename BoxData::Type qlength_m3 = qlength_m2*qlength_m1;
                const typename BoxData::Type Omega_1 =
                    qlength_m3*(data.myNijDiag[0] + data.myNijDiag[1] + data.myNijDiag[2]
                        -typename BoxData::Type(3.0)*(dot(q2,data.myNijDiag) +
                            q[0]*q[1]*data.myNxy_Nyx +
                            q[0]*q[2]*data.myNzx_Nxz +
                            q[1]*q[2]*data.myNyz_Nzy));
                Omega_approx += Omega_1;
#if TAYLOR_SERIES_ORDER >= 2
                if (order >= 2)
                {
                    const UT_FixedVector<typename BoxData::Type,3> q3 = q2*q;
                    const typename BoxData::Type qlength_m4 = qlength_m2*qlength_m2;
                    typename BoxData::Type temp0[3] = {
                        data.my2Nyyx_Nxyy+data.my2Nzzx_Nxzz,
                        data.my2Nzzy_Nyzz+data.my2Nxxy_Nyxx,
                        data.my2Nxxz_Nzxx+data.my2Nyyz_Nzyy
                    };
                    typename BoxData::Type temp1[3] = {
                        q[1]*data.my2Nxxy_Nyxx + q[2]*data.my2Nxxz_Nzxx,
                        q[2]*data.my2Nyyz_Nzyy + q[0]*data.my2Nyyx_Nxyy,
                        q[0]*data.my2Nzzx_Nxzz + q[1]*data.my2Nzzy_Nyzz
                    };
                    const typename BoxData::Type Omega_2 =
                        qlength_m4*(typename BoxData::Type(1.5)*dot(q, typename BoxData::Type(3)*data.myNijkDiag + UT_FixedVector<typename BoxData::Type,3>(temp0))
                            -typename BoxData::Type(7.5)*(dot(q3,data.myNijkDiag) + q[0]*q[1]*q[2]*data.mySumPermuteNxyz + dot(q2, UT_FixedVector<typename BoxData::Type,3>(temp1))));
                    Omega_approx += Omega_2;
                }
#endif
            }
#endif

            // If q is so small that we got NaNs and we just have a
            // small bounding box, it needs to descend.
            const v4uu mask = Omega_approx.isFinite() & ~descend_mask;
            Omega_approx = Omega_approx & mask;
            descend_bitmask = (~_mm_movemask_ps(V4SF(mask.vector))) & allchildbits;

            T sum = Omega_approx[0];
            for (int i = 1; i < BVH_N; ++i)
                sum += Omega_approx[i];
            *data_for_parent = sum;

            return descend_bitmask;
        }
        void item(const int itemi, const int parent_nodei, T &data_for_parent) const
        {
            const UT_Vector3T<S> *const positions = myPositions;
            const int *const cur_triangle_points = myTrianglePoints + 3*itemi;
            const UT_Vector3T<T> a = positions[cur_triangle_points[0]];
            const UT_Vector3T<T> b = positions[cur_triangle_points[1]];
            const UT_Vector3T<T> c = positions[cur_triangle_points[2]];

            data_for_parent = UTsignedSolidAngleTri(a, b, c, myQueryPoint);
        }
        SYS_FORCE_INLINE void post(const int nodei, const int parent_nodei, T *data_for_parent, const int nchildren, const T *child_data_array, const uint descend_bits) const
        {
            T sum = (descend_bits&1) ? child_data_array[0] : 0;
            for (int i = 1; i < nchildren; ++i)
                sum += ((descend_bits>>i)&1) ? child_data_array[i] : 0;

            *data_for_parent += sum;
        }
    };
    const SolidAngleFunctors functors(myData.get(), query_point, accuracy_scale2, myOrder, myPositions, myTrianglePoints);

    T sum;
    myTree.traverseVector(functors, &sum);
    return sum;
}

template<typename T,typename S>
struct UT_SubtendedAngle<T,S>::BoxData
{
    void clear()
    {
        // Set everything to zero
        memset(this,0,sizeof(*this));
    }

    using Type  = typename std::conditional<BVH_N==4 && std::is_same<T,float>::value, v4uf, UT_FixedVector<T,BVH_N>>::type;
    using SType = typename std::conditional<BVH_N==4 && std::is_same<S,float>::value, v4uf, UT_FixedVector<S,BVH_N>>::type;

    /// An upper bound on the squared distance from myAverageP to the farthest point in the box.
    SType myMaxPDist2;

    /// Centre of mass of the mesh surface in this box
    UT_FixedVector<Type,2> myAverageP;

    /// Unnormalized, area-weighted normal of the mesh in this box
    UT_FixedVector<Type,2> myN;

    /// Values for Omega_1
    /// @{
    UT_FixedVector<Type,2> myNijDiag;  // Nxx, Nyy
    Type myNxy_Nyx;               // Nxy+Nyx
    /// @}

    /// Values for Omega_2
    /// @{
    UT_FixedVector<Type,2> myNijkDiag; // Nxxx, Nyyy
    Type my2Nxxy_Nyxx; // Nxxy+Nxyx+Nyxx = 2Nxxy+Nyxx
    Type my2Nyyx_Nxyy; // Nyyx+Nyxy+Nxyy = 2Nyyx+Nxyy
    /// @}
};

template<typename T,typename S>
UT_SubtendedAngle<T,S>::UT_SubtendedAngle()
    : myTree()
    , myNBoxes(0)
    , myOrder(2)
    , myData(nullptr)
    , myNSegments(0)
    , mySegmentPoints(nullptr)
    , myNPoints(0)
    , myPositions(nullptr)
{}

template<typename T,typename S>
UT_SubtendedAngle<T,S>::~UT_SubtendedAngle()
{
    // Default destruction works, but this needs to be outlined
    // to avoid having to include UT_BVHImpl.h in the header,
    // (for the UT_UniquePtr destructor.)
}

template<typename T,typename S>
void UT_SubtendedAngle<T,S>::init(
    const int nsegments,
    const int *const segment_points,
    const int npoints,
    const UT_Vector2T<S> *const positions,
    const int order)
{
#if SOLID_ANGLE_DEBUG
    UTdebugFormat("");
    UTdebugFormat("");
    UTdebugFormat("Building BVH for {} segments on {} points:", nsegments, npoints);
#endif
    myOrder = order;
    myNSegments = nsegments;
    mySegmentPoints = segment_points;
    myNPoints = npoints;
    myPositions = positions;

#if SOLID_ANGLE_TIME_PRECOMPUTE
    UT_StopWatch timer;
    timer.start();
#endif
    UT_SmallArray<UT::Box<S,2>> segment_boxes;
    segment_boxes.setSizeNoInit(nsegments);
    if (nsegments < 16*1024)
    {
        const int *cur_segment_points = segment_points;
        for (int i = 0; i < nsegments; ++i, cur_segment_points += 2)
        {
            UT::Box<S,2> &box = segment_boxes[i];
            box.initBounds(positions[cur_segment_points[0]]);
            box.enlargeBounds(positions[cur_segment_points[1]]);
        }
    }
    else
    {
        UTparallelFor(UT_BlockedRange<int>(0,nsegments), [segment_points,&segment_boxes,positions](const UT_BlockedRange<int> &r)
        {
            const int *cur_segment_points = segment_points + exint(r.begin())*2;
            for (int i = r.begin(), end = r.end(); i < end; ++i, cur_segment_points += 2)
            {
                UT::Box<S,2> &box = segment_boxes[i];
                box.initBounds(positions[cur_segment_points[0]]);
                box.enlargeBounds(positions[cur_segment_points[1]]);
            }
        });
    }
#if SOLID_ANGLE_TIME_PRECOMPUTE
    double time = timer.stop();
    UTdebugFormat("{} s to create bounding boxes.", time);
    timer.start();
#endif
    myTree.template init<UT::BVH_Heuristic::BOX_AREA,S,2>(segment_boxes.array(), nsegments);
#if SOLID_ANGLE_TIME_PRECOMPUTE
    time = timer.stop();
    UTdebugFormat("{} s to initialize UT_BVH structure.  {} nodes", time, myTree.getNumNodes());
#endif

    //myTree.debugDump();

    const int nnodes = myTree.getNumNodes();

    myNBoxes = nnodes;
    BoxData *box_data = new BoxData[nnodes];
    myData.reset(box_data);

    // Some data are only needed during initialization.
    struct LocalData
    {
        // Bounding box
        UT::Box<S,2> myBox;

        // P and N are needed from each child for computing Nij.
        UT_Vector2T<T> myAverageP;
        UT_Vector2T<T> myLengthP;
        UT_Vector2T<T> myN;

        // Unsigned length is needed for computing the average position.
        T myLength;

        // These are needed for computing Nijk.
        UT_Vector2T<T> myNijDiag;
        T myNxy; T myNyx;

        UT_Vector2T<T> myNijkDiag; // Nxxx, Nyyy
        T my2Nxxy_Nyxx;     // Nxxy+Nxyx+Nyxx = 2Nxxy+Nyxx
        T my2Nyyx_Nxyy;     // Nyyx+Nyxy+Nxyy = 2Nyyx+Nxyy
    };

    struct PrecomputeFunctors
    {
        BoxData *const myBoxData;
        const UT::Box<S,2> *const mySegmentBoxes;
        const int *const mySegmentPoints;
        const UT_Vector2T<S> *const myPositions;
        const int myOrder;

        PrecomputeFunctors(
            BoxData *box_data,
            const UT::Box<S,2> *segment_boxes,
            const int *segment_points,
            const UT_Vector2T<S> *positions,
            const int order)
            : myBoxData(box_data)
            , mySegmentBoxes(segment_boxes)
            , mySegmentPoints(segment_points)
            , myPositions(positions)
            , myOrder(order)
        {}
        constexpr SYS_FORCE_INLINE bool pre(const int nodei, LocalData *data_for_parent) const
        {
            return true;
        }
        void item(const int itemi, const int parent_nodei, LocalData &data_for_parent) const
        {
            const UT_Vector2T<S> *const positions = myPositions;
            const int *const cur_segment_points = mySegmentPoints + 2*itemi;
            const UT_Vector2T<T> a = positions[cur_segment_points[0]];
            const UT_Vector2T<T> b = positions[cur_segment_points[1]];
            const UT_Vector2T<T> ab = b-a;

            const UT::Box<S,2> &segment_box = mySegmentBoxes[itemi];
            data_for_parent.myBox = segment_box;

            // Length-weighted normal (unnormalized)
            UT_Vector2T<T> N;
            N[0] = ab[1];
            N[1] = -ab[0];
            const T length2 = ab.length2();
            const T length = SYSsqrt(length2);
            const UT_Vector2T<T> P = T(0.5)*(a+b);
            data_for_parent.myAverageP = P;
            data_for_parent.myLengthP = P*length;
            data_for_parent.myN = N;
#if SOLID_ANGLE_DEBUG
            UTdebugFormat("");
            UTdebugFormat("Triangle {}: P = {}; N = {}; length = {}", itemi, P, N, length);
            UTdebugFormat("             box = {}", data_for_parent.myBox);
#endif

            data_for_parent.myLength = length;
            const int order = myOrder;
            if (order < 1)
                return;

            // NOTE: Due to P being at the centroid, segments have Nij = 0
            //       contributions to Nij.
            data_for_parent.myNijDiag = T(0);
            data_for_parent.myNxy = 0; data_for_parent.myNyx = 0;

            if (order < 2)
                return;

            // If it's zero-length, the results are zero, so we can skip.
            if (length == 0)
            {
                data_for_parent.myNijkDiag = T(0);
                data_for_parent.my2Nxxy_Nyxx = 0;
                data_for_parent.my2Nyyx_Nxyy = 0;
                return;
            }

            T integral_xx = ab[0]*ab[0]/T(12);
            T integral_xy = ab[0]*ab[1]/T(12);
            T integral_yy = ab[1]*ab[1]/T(12);
            data_for_parent.myNijkDiag[0] = integral_xx*N[0];
            data_for_parent.myNijkDiag[1] = integral_yy*N[1];
            T Nxxy = N[0]*integral_xy;
            T Nyxx = N[1]*integral_xx;
            T Nyyx = N[1]*integral_xy;
            T Nxyy = N[0]*integral_yy;
            data_for_parent.my2Nxxy_Nyxx = 2*Nxxy + Nyxx;
            data_for_parent.my2Nyyx_Nxyy = 2*Nyyx + Nxyy;
#if SOLID_ANGLE_DEBUG
            UTdebugFormat("             integral_xx = {}; yy = {}", integral_xx, integral_yy);
            UTdebugFormat("             integral_xy = {}", integral_xy);
#endif
        }

        void post(const int nodei, const int parent_nodei, LocalData *data_for_parent, const int nchildren, const LocalData *child_data_array) const
        {
            // NOTE: Although in the general case, data_for_parent may be null for the root call,
            //       this functor assumes that it's non-null, so the call below must pass a non-null pointer.

            BoxData &current_box_data = myBoxData[nodei];

            UT_Vector2T<T> N = child_data_array[0].myN;
            ((T*)&current_box_data.myN[0])[0] = N[0];
            ((T*)&current_box_data.myN[1])[0] = N[1];
            UT_Vector2T<T> lengthP = child_data_array[0].myLengthP;
            T length = child_data_array[0].myLength;
            const UT_Vector2T<T> local_P = child_data_array[0].myAverageP;
            ((T*)&current_box_data.myAverageP[0])[0] = local_P[0];
            ((T*)&current_box_data.myAverageP[1])[0] = local_P[1];
            for (int i = 1; i < nchildren; ++i)
            {
                const UT_Vector2T<T> local_N = child_data_array[i].myN;
                N += local_N;
                ((T*)&current_box_data.myN[0])[i] = local_N[0];
                ((T*)&current_box_data.myN[1])[i] = local_N[1];
                lengthP += child_data_array[i].myLengthP;
                length += child_data_array[i].myLength;
                const UT_Vector2T<T> local_P = child_data_array[i].myAverageP;
                ((T*)&current_box_data.myAverageP[0])[i] = local_P[0];
                ((T*)&current_box_data.myAverageP[1])[i] = local_P[1];
            }
            for (int i = nchildren; i < BVH_N; ++i)
            {
                // Set to zero, just to avoid false positives for uses of uninitialized memory.
                ((T*)&current_box_data.myN[0])[i] = 0;
                ((T*)&current_box_data.myN[1])[i] = 0;
                ((T*)&current_box_data.myAverageP[0])[i] = 0;
                ((T*)&current_box_data.myAverageP[1])[i] = 0;
            }
            data_for_parent->myN = N;
            data_for_parent->myLengthP = lengthP;
            data_for_parent->myLength = length;

            UT::Box<S,2> box(child_data_array[0].myBox);
            for (int i = 1; i < nchildren; ++i)
                box.combine(child_data_array[i].myBox);

            // Normalize P
            UT_Vector2T<T> averageP;
            if (length > 0)
                averageP = lengthP/length;
            else
                averageP = T(0.5)*(box.getMin() + box.getMax());
            data_for_parent->myAverageP = averageP;

            data_for_parent->myBox = box;

            for (int i = 0; i < nchildren; ++i)
            {
                const UT::Box<S,2> &local_box(child_data_array[i].myBox);
                const UT_Vector2T<T> &local_P = child_data_array[i].myAverageP;
                const UT_Vector2T<T> maxPDiff = SYSmax(local_P-UT_Vector2T<T>(local_box.getMin()), UT_Vector2T<T>(local_box.getMax())-local_P);
                ((T*)&current_box_data.myMaxPDist2)[i] = maxPDiff.length2();
            }
            for (int i = nchildren; i < BVH_N; ++i)
            {
                // This child is non-existent.  If we set myMaxPDist2 to infinity, it will never
                // use the approximation, and the traverseVector function can check for EMPTY.
                ((T*)&current_box_data.myMaxPDist2)[i] = std::numeric_limits<T>::infinity();
            }

            const int order = myOrder;
            if (order >= 1)
            {
                // We now have the current box's P, so we can adjust Nij and Nijk
                data_for_parent->myNijDiag = child_data_array[0].myNijDiag;
                data_for_parent->myNxy = 0;
                data_for_parent->myNyx = 0;
                data_for_parent->myNijkDiag = child_data_array[0].myNijkDiag;
                data_for_parent->my2Nxxy_Nyxx = child_data_array[0].my2Nxxy_Nyxx;
                data_for_parent->my2Nyyx_Nxyy = child_data_array[0].my2Nyyx_Nxyy;

                for (int i = 1; i < nchildren; ++i)
                {
                    data_for_parent->myNijDiag += child_data_array[i].myNijDiag;
                    data_for_parent->myNijkDiag += child_data_array[i].myNijkDiag;
                    data_for_parent->my2Nxxy_Nyxx += child_data_array[i].my2Nxxy_Nyxx;
                    data_for_parent->my2Nyyx_Nxyy += child_data_array[i].my2Nyyx_Nxyy;
                }
                for (int j = 0; j < 2; ++j)
                    ((T*)&current_box_data.myNijDiag[j])[0] = child_data_array[0].myNijDiag[j];
                ((T*)&current_box_data.myNxy_Nyx)[0] = child_data_array[0].myNxy + child_data_array[0].myNyx;
                for (int j = 0; j < 2; ++j)
                    ((T*)&current_box_data.myNijkDiag[j])[0] = child_data_array[0].myNijkDiag[j];
                ((T*)&current_box_data.my2Nxxy_Nyxx)[0] = child_data_array[0].my2Nxxy_Nyxx;
                ((T*)&current_box_data.my2Nyyx_Nxyy)[0] = child_data_array[0].my2Nyyx_Nxyy;
                for (int i = 1; i < nchildren; ++i)
                {
                    for (int j = 0; j < 2; ++j)
                        ((T*)&current_box_data.myNijDiag[j])[i] = child_data_array[i].myNijDiag[j];
                    ((T*)&current_box_data.myNxy_Nyx)[i] = child_data_array[i].myNxy + child_data_array[i].myNyx;
                    for (int j = 0; j < 2; ++j)
                        ((T*)&current_box_data.myNijkDiag[j])[i] = child_data_array[i].myNijkDiag[j];
                    ((T*)&current_box_data.my2Nxxy_Nyxx)[i] = child_data_array[i].my2Nxxy_Nyxx;
                    ((T*)&current_box_data.my2Nyyx_Nxyy)[i] = child_data_array[i].my2Nyyx_Nxyy;
                }
                for (int i = nchildren; i < BVH_N; ++i)
                {
                    // Set to zero, just to avoid false positives for uses of uninitialized memory.
                    for (int j = 0; j < 2; ++j)
                        ((T*)&current_box_data.myNijDiag[j])[i] = 0;
                    ((T*)&current_box_data.myNxy_Nyx)[i] = 0;
                    for (int j = 0; j < 2; ++j)
                        ((T*)&current_box_data.myNijkDiag[j])[i] = 0;
                    ((T*)&current_box_data.my2Nxxy_Nyxx)[i] = 0;
                    ((T*)&current_box_data.my2Nyyx_Nxyy)[i] = 0;
                }

                for (int i = 0; i < nchildren; ++i)
                {
                    const LocalData &child_data = child_data_array[i];
                    UT_Vector2T<T> displacement = child_data.myAverageP - UT_Vector2T<T>(data_for_parent->myAverageP);
                    UT_Vector2T<T> N = child_data.myN;

                    // Adjust Nij for the change in centre P
                    data_for_parent->myNijDiag += N*displacement;
                    T Nxy = child_data.myNxy + N[0]*displacement[1];
                    T Nyx = child_data.myNyx + N[1]*displacement[0];

                    data_for_parent->myNxy += Nxy;
                    data_for_parent->myNyx += Nyx;

                    if (order >= 2)
                    {
                        // Adjust Nijk for the change in centre P
                        data_for_parent->myNijkDiag += T(2)*displacement*child_data.myNijDiag + displacement*displacement*child_data.myN;
                        data_for_parent->my2Nxxy_Nyxx +=
                            2*(displacement[1]*child_data.myNijDiag[0] + displacement[0]*child_data.myNxy + N[0]*displacement[0]*displacement[1])
                            + 2*child_data.myNyx*displacement[0] + N[1]*displacement[0]*displacement[0];
                        data_for_parent->my2Nyyx_Nxyy +=
                            2*(displacement[0]*child_data.myNijDiag[1] + displacement[1]*child_data.myNyx + N[1]*displacement[1]*displacement[0])
                            + 2*child_data.myNxy*displacement[1] + N[0]*displacement[1]*displacement[1];
                    }
                }
            }
#if SOLID_ANGLE_DEBUG
            UTdebugFormat("");
            UTdebugFormat("Node {}: nchildren = {}; maxP = {}", nodei, nchildren, SYSsqrt(current_box_data.myMaxPDist2));
            UTdebugFormat("         P = {}; N = {}", current_box_data.myAverageP, current_box_data.myN);
            UTdebugFormat("         Nii = {}", current_box_data.myNijDiag);
            UTdebugFormat("         Nxy+Nyx = {}", current_box_data.myNxy_Nyx);
            UTdebugFormat("         Niii = {}", current_box_data.myNijkDiag);
            UTdebugFormat("         2Nxxy+Nyxx = {}; 2Nyyx+Nxyy = {}", current_box_data.my2Nxxy_Nyxx, current_box_data.my2Nyyx_Nxyy);
#endif
        }
    };

#if SOLID_ANGLE_TIME_PRECOMPUTE
    timer.start();
#endif
    const PrecomputeFunctors functors(box_data, segment_boxes.array(), segment_points, positions, order);
    // NOTE: post-functor relies on non-null data_for_parent, so we have to pass one.
    LocalData local_data;
    myTree.template traverseParallel<LocalData>(4096, functors, &local_data);
    //myTree.template traverse<LocalData>(functors);
#if SOLID_ANGLE_TIME_PRECOMPUTE
    time = timer.stop();
    UTdebugFormat("{} s to precompute coefficients.", time);
#endif
}

template<typename T,typename S>
void UT_SubtendedAngle<T, S>::clear()
{
    myTree.clear();
    myNBoxes = 0;
    myOrder = 2;
    myData.reset();
    myNSegments = 0;
    mySegmentPoints = nullptr;
    myNPoints = 0;
    myPositions = nullptr;
}

template<typename T,typename S>
T UT_SubtendedAngle<T, S>::computeAngle(const UT_Vector2T<T> &query_point, const T accuracy_scale) const
{
    const T accuracy_scale2 = accuracy_scale*accuracy_scale;

    struct AngleFunctors
    {
        const BoxData *const myBoxData;
        const UT_Vector2T<T> myQueryPoint;
        const T myAccuracyScale2;
        const UT_Vector2T<S> *const myPositions;
        const int *const mySegmentPoints;
        const int myOrder;

        AngleFunctors(
            const BoxData *const box_data,
            const UT_Vector2T<T> &query_point,
            const T accuracy_scale2,
            const int order,
            const UT_Vector2T<S> *const positions,
            const int *const segment_points)
            : myBoxData(box_data)
            , myQueryPoint(query_point)
            , myAccuracyScale2(accuracy_scale2)
            , myOrder(order)
            , myPositions(positions)
            , mySegmentPoints(segment_points)
        {}
        uint pre(const int nodei, T *data_for_parent) const
        {
            const BoxData &data = myBoxData[nodei];
            const typename BoxData::Type maxP2 = data.myMaxPDist2;
            UT_FixedVector<typename BoxData::Type,2> q;
            q[0] = typename BoxData::Type(myQueryPoint[0]);
            q[1] = typename BoxData::Type(myQueryPoint[1]);
            q -= data.myAverageP;
            const typename BoxData::Type qlength2 = q[0]*q[0] + q[1]*q[1];

            // If the query point is within a factor of accuracy_scale of the box radius,
            // it's assumed to be not a good enough approximation, so it needs to descend.
            // TODO: Is there a way to estimate the error?
            static_assert((std::is_same<typename BoxData::Type,v4uf>::value), "FIXME: Implement support for other tuple types!");
            v4uu descend_mask = (qlength2 <= maxP2*myAccuracyScale2);
            uint descend_bitmask = _mm_movemask_ps(V4SF(descend_mask.vector));
            constexpr uint allchildbits = ((uint(1)<<BVH_N)-1);
            if (descend_bitmask == allchildbits)
            {
                *data_for_parent = 0;
                return allchildbits;
            }

            // qlength2 must be non-zero, since it's strictly greater than something.
            // We still need to be careful for NaNs, though, because the 4th power might cause problems.
            const typename BoxData::Type qlength_m2 = typename BoxData::Type(1.0)/qlength2;
            const typename BoxData::Type qlength_m1 = sqrt(qlength_m2);

            // Normalize q to reduce issues with overflow/underflow, since we'd need the 6th power
            // if we didn't normalize, and (1e-7)^-6 = 1e42, which overflows single-precision.
            q *= qlength_m1;

            typename BoxData::Type Omega_approx = -qlength_m1*dot(q,data.myN);
            const int order = myOrder;
            if (order >= 1)
            {
                const UT_FixedVector<typename BoxData::Type,2> q2 = q*q;
                const typename BoxData::Type Omega_1 =
                    qlength_m2*(data.myNijDiag[0] + data.myNijDiag[1]
                        -typename BoxData::Type(2.0)*(dot(q2,data.myNijDiag) +
                            q[0]*q[1]*data.myNxy_Nyx));
                Omega_approx += Omega_1;
                if (order >= 2)
                {
                    const UT_FixedVector<typename BoxData::Type,2> q3 = q2*q;
                    const typename BoxData::Type qlength_m3 = qlength_m2*qlength_m1;
                    typename BoxData::Type temp0[2] = {
                        data.my2Nyyx_Nxyy,
                        data.my2Nxxy_Nyxx
                    };
                    typename BoxData::Type temp1[2] = {
                        q[1]*data.my2Nxxy_Nyxx,
                        q[0]*data.my2Nyyx_Nxyy
                    };
                    const typename BoxData::Type Omega_2 =
                        qlength_m3*(dot(q, typename BoxData::Type(3)*data.myNijkDiag + UT_FixedVector<typename BoxData::Type,2>(temp0))
                            -typename BoxData::Type(4.0)*(dot(q3,data.myNijkDiag) + dot(q2, UT_FixedVector<typename BoxData::Type,2>(temp1))));
                    Omega_approx += Omega_2;
                }
            }

            // If q is so small that we got NaNs and we just have a
            // small bounding box, it needs to descend.
            const v4uu mask = Omega_approx.isFinite() & ~descend_mask;
            Omega_approx = Omega_approx & mask;
            descend_bitmask = (~_mm_movemask_ps(V4SF(mask.vector))) & allchildbits;

            T sum = Omega_approx[0];
            for (int i = 1; i < BVH_N; ++i)
                sum += Omega_approx[i];
            *data_for_parent = sum;

            return descend_bitmask;
        }
        void item(const int itemi, const int parent_nodei, T &data_for_parent) const
        {
            const UT_Vector2T<S> *const positions = myPositions;
            const int *const cur_segment_points = mySegmentPoints + 2*itemi;
            const UT_Vector2T<T> a = positions[cur_segment_points[0]];
            const UT_Vector2T<T> b = positions[cur_segment_points[1]];

            data_for_parent = UTsignedAngleSegment(a, b, myQueryPoint);
        }
        SYS_FORCE_INLINE void post(const int nodei, const int parent_nodei, T *data_for_parent, const int nchildren, const T *child_data_array, const uint descend_bits) const
        {
            T sum = (descend_bits&1) ? child_data_array[0] : 0;
            for (int i = 1; i < nchildren; ++i)
                sum += ((descend_bits>>i)&1) ? child_data_array[i] : 0;

            *data_for_parent += sum;
        }
    };
    const AngleFunctors functors(myData.get(), query_point, accuracy_scale2, myOrder, myPositions, mySegmentPoints);

    T sum;
    myTree.traverseVector(functors, &sum);
    return sum;
}

// Instantiate our templates.
template class UT_SolidAngle<fpreal32,fpreal32>;
// FIXME: The SIMD parts will need to be handled differently in order to support fpreal64.
//template class UT_SolidAngle<fpreal64,fpreal32>;
//template class UT_SolidAngle<fpreal64,fpreal64>;
template class UT_SubtendedAngle<fpreal32,fpreal32>;
//template class UT_SubtendedAngle<fpreal64,fpreal32>;
//template class UT_SubtendedAngle<fpreal64,fpreal64>;

} // End HDK_Sample namespace
