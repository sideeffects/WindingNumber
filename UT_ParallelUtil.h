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
 *      Simple wrappers on tbb interface
 */

#ifndef __UT_ParallelUtil__
#define __UT_ParallelUtil__

#include "UT_Array.h"

#include <thread> // This is just included for std::thread::hardware_concurrency()
namespace UT_Thread { int getNumProcessors() {
    return std::thread::hardware_concurrency();
}}

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
//namespace tbb { class split; }

/// Declare prior to use.
template <typename T> 
using UT_BlockedRange = tbb::blocked_range<T>;

// Default implementation that calls range.size()
template< typename RANGE >
struct UT_EstimatorNumItems
{
    UT_EstimatorNumItems() {}

    size_t operator()(const RANGE& range) const
    {
	return range.size();
    }
};

/// This is needed by UT_CoarsenedRange
template <typename RANGE>
inline size_t UTestimatedNumItems(const RANGE& range)
{
    return UT_EstimatorNumItems<RANGE>()(range);
}

/// UT_CoarsenedRange: This should be used only inside 
/// UT_ParallelFor and UT_ParallelReduce
/// This class wraps an existing range with a new range.
/// This allows us to use simple_partitioner, rather than
/// auto_partitioner, which has disastrous performance with
/// the default grain size in ttb 4.
template< typename RANGE >
class UT_CoarsenedRange : public RANGE
{
public:
    // Compiler-generated versions are fine:
    // ~UT_CoarsenedRange();
    // UT_CoarsenedRange(const UT_CoarsenedRange&);

    // Split into two sub-ranges:
    UT_CoarsenedRange(UT_CoarsenedRange& range, tbb::split spl) :
        RANGE(range, spl),
        myGrainSize(range.myGrainSize)
    {        
    }

    // Inherited: bool empty() const

    bool is_divisible() const
    {
        return 
            RANGE::is_divisible() &&
            (UTestimatedNumItems(static_cast<const RANGE&>(*this)) > myGrainSize);
    }

private:
    size_t myGrainSize;

    UT_CoarsenedRange(const RANGE& base_range, const size_t grain_size) :
        RANGE(base_range),
        myGrainSize(grain_size)
    {        
    }

    template <typename Range, typename Body>
    friend void UTparallelFor(
        const Range &range, const Body &body,
        const int subscribe_ratio, const int min_grain_size
    );
};

/// Run the @c body function over a range in parallel.
/// UTparallelFor attempts to spread the range out over at most 
/// subscribe_ratio * num_processor tasks.
/// The factor subscribe_ratio can be used to help balance the load.
/// UTparallelFor() uses tbb for its implementation.
/// The used grain size is the maximum of min_grain_size and
/// if UTestimatedNumItems(range) / (subscribe_ratio * num_processor).
/// If subscribe_ratio == 0, then a grain size of min_grain_size will be used.
/// A range can be split only when UTestimatedNumItems(range) exceeds the
/// grain size the range is divisible. 

///
/// Requirements for the Range functor are:
///   - the requirements of the tbb Range Concept
///   - UT_estimatorNumItems<Range> must return the the estimated number of work items
///     for the range. When Range::size() is not the correct estimate, then a 
///     (partial) specialization of UT_estimatorNumItemsimatorRange must be provided
///     for the type Range.
///
/// Requirements for the Body function are:
///  - @code Body(const Body &); @endcode @n
///	Copy Constructor
///  - @code Body()::~Body(); @endcode @n
///	Destructor
///  - @code void Body::operator()(const Range &range) const; @endcode
///	Function call to perform operation on the range.  Note the operator is
///	@b const.
///
/// The requirements for a Range object are:
///  - @code Range::Range(const Range&); @endcode @n
///	Copy constructor
///  - @code Range::~Range(); @endcode @n
///	Destructor
///  - @code bool Range::is_divisible() const; @endcode @n
///	True if the range can be partitioned into two sub-ranges
///  - @code bool Range::empty() const; @endcode @n
///	True if the range is empty
///  - @code Range::Range(Range &r, UT_Split) const; @endcode @n
///	Split the range @c r into two sub-ranges (i.e. modify @c r and *this)
///
/// Example: @code
///     class Square {
///     public:
///         Square(double *data) : myData(data) {}
///         ~Square();
///         void operator()(const UT_BlockedRange<int64> &range) const
///         {
///             for (int64 i = range.begin(); i != range.end(); ++i)
///                 myData[i] *= myData[i];
///         }
///         double *myData;
///     };
///     ...
///
///     void
///     parallel_square(double *array, int64 length)
///     {
///         UTparallelFor(UT_BlockedRange<int64>(0, length), Square(array));
///     }
/// @endcode
///	
/// @see UTparallelReduce(), UT_BlockedRange()

template <typename Range, typename Body>
void UTparallelFor(
    const Range &range, const Body &body,
    const int subscribe_ratio = 2,
    const int min_grain_size = 1
)
{
    const size_t num_processors( UT_Thread::getNumProcessors() );

    UT_ASSERT( num_processors >= 1 );
    UT_ASSERT( min_grain_size >= 1 );
    UT_ASSERT( subscribe_ratio >= 0 );

    const size_t est_range_size( UTestimatedNumItems(range) );

    // Don't run on an empty range!
    if (est_range_size == 0)
        return;

    // Avoid tbb overhead if entire range needs to be single threaded
    if (num_processors == 1 || est_range_size <= min_grain_size)
    {
        body(range);
        return;
    }

    size_t grain_size(min_grain_size);
    if( subscribe_ratio > 0 )
        grain_size = std::max(
                         grain_size, 
                         est_range_size / (subscribe_ratio * num_processors)
                     );

    UT_CoarsenedRange< Range > coarsened_range(range, grain_size);

    tbb::parallel_for(coarsened_range, body, tbb::simple_partitioner());
}

/// Version of UTparallelFor that is tuned for the case where the range
/// consists of lightweight items, for example,
/// float additions or matrix-vector multiplications.
template <typename Range, typename Body>
void
UTparallelForLightItems(const Range &range, const Body &body)
{
    UTparallelFor(range, body, 2, 1024);
}

/// UTserialFor can be used as a debugging tool to quickly replace a parallel
/// for with a serial for.
template <typename Range, typename Body>
void UTserialFor(const Range &range, const Body &body)
	{ body(range); }

#endif
