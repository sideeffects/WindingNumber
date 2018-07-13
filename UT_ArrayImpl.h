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
 *      This is meant to be included by UT_Array.h and includes
 *      the template implementations needed by external code.
 */

#pragma once

#ifndef __UT_ARRAYIMPL_H_INCLUDED__
#define __UT_ARRAYIMPL_H_INCLUDED__

#include "SYS_Math.h"
#include "SYS_Types.h"

#include <algorithm>
#include <utility>
#include <stdlib.h>
#include <string.h>


// Implemented in UT_Array.C
extern void ut_ArrayImplFree(void *p);


template <typename T>
UT_Array<T>::UT_Array(const UT_Array<T> &a)
    : myCapacity(a.size()), mySize(a.size())
{
    if (myCapacity)
    {
	myData = allocateCapacity(myCapacity);
	copyConstructRange(myData, a.array(), mySize);
    }
    else
    {
	myData = nullptr;
    }
}

template <typename T>
UT_Array<T>::UT_Array(std::initializer_list<T> init)
    : myCapacity(init.size()), mySize(init.size())
{
    if (myCapacity)
    {
	myData = allocateCapacity(myCapacity);
	copyConstructRange(myData, init.begin(), mySize);
    }
    else
    {
	myData = nullptr;
    }
}

template <typename T>
UT_Array<T>::UT_Array(UT_Array<T> &&a) noexcept
{
    if (!a.isHeapBuffer())
    {
	myData = nullptr;
	myCapacity = 0;
	mySize = 0;
	operator=(std::move(a));
	return;
    }

    myCapacity = a.myCapacity;
    mySize = a.mySize;
    myData = a.myData;
    a.myCapacity = a.mySize = 0;
    a.myData = nullptr;
}


template <typename T>
UT_Array<T>::~UT_Array()
{
    // NOTE: We call setCapacity to ensure that we call trivialDestructRange,
    //       then call free on myData.
    setCapacity(0);
}

template <typename T>
inline T *
UT_Array<T>::allocateCapacity(exint capacity)
{
    T *data = (T *)malloc(capacity * sizeof(T));
    // Avoid degenerate case if we happen to be aliased the wrong way
    if (!isHeapBuffer(data))
    {
	T *prev = data;
	data = (T *)malloc(capacity * sizeof(T));
	ut_ArrayImplFree(prev);
    }
    return data;
}

template <typename T>
void
UT_Array<T>::swap( UT_Array<T> &other )
{
    std::swap( myData, other.myData );
    std::swap( myCapacity, other.myCapacity );
    std::swap( mySize, other.mySize );
}


template <typename T>
exint	
UT_Array<T>::insert(exint index)
{
    if (index >= mySize)
    {
	bumpCapacity(index + 1);

	trivialConstructRange(myData + mySize, index - mySize + 1);

	mySize = index+1;
	return index;
    }
    bumpCapacity(mySize + 1);

    UT_ASSERT_P(index >= 0);
    ::memmove((void *)&myData[index+1], (void *)&myData[index],
              ((mySize-index)*sizeof(T)));

    trivialConstruct(myData[index]);

    mySize++;
    return index;
}

template <typename T>
template <typename S>
exint
UT_Array<T>::appendImpl(S &&s)
{
    if (mySize == myCapacity)
    {
	exint idx = safeIndex(s);

        // NOTE: UTbumpAlloc always returns a strictly larger value.
	setCapacity(UTbumpAlloc(myCapacity));
        if (idx >= 0)
            construct(myData[mySize], std::forward<S>(myData[idx]));
        else
            construct(myData[mySize], std::forward<S>(s));
    }
    else
    {
	construct(myData[mySize], std::forward<S>(s));
    }
    return mySize++;
}

template <typename T>
template <typename... S>
exint
UT_Array<T>::emplace_back(S&&... s)
{
    if (mySize == myCapacity)
	setCapacity(UTbumpAlloc(myCapacity));

    construct(myData[mySize], std::forward<S>(s)...);
    return mySize++;
}

template <typename T>
void
UT_Array<T>::append(const T *pt, exint count)
{
    bumpCapacity(mySize + count);
    copyConstructRange(myData + mySize, pt, count);
    mySize += count;
}

template <typename T>
void
UT_Array<T>::appendMultiple(const T &t, exint count)
{
    UT_ASSERT_P(count >= 0);
    if (count <= 0)
	return;
    if (mySize + count >= myCapacity)
    {
	exint tidx = safeIndex(t);

        bumpCapacity(mySize + count);

	for (exint i = 0; i < count; i++)
	    copyConstruct(myData[mySize+i], tidx >= 0 ? myData[tidx] : t);
    }
    else
    {
	for (exint i = 0; i < count; i++)
	    copyConstruct(myData[mySize+i], t);
    }
    mySize += count;
}

template <typename T>
exint
UT_Array<T>::concat(const UT_Array<T> &a)
{
    bumpCapacity(mySize + a.mySize);
    copyConstructRange(myData + mySize, a.myData, a.mySize);
    mySize += a.mySize;

    return mySize;
}

template <typename T>
exint
UT_Array<T>::multipleInsert(exint beg_index, exint count)
{
    exint end_index = beg_index + count;

    if (beg_index >= mySize)
    {
	bumpCapacity(end_index);

	trivialConstructRange(myData + mySize, end_index - mySize);

	mySize = end_index;
	return beg_index;
    }
    bumpCapacity(mySize+count);

    ::memmove((void *)&myData[end_index], (void *)&myData[beg_index],
              ((mySize-beg_index)*sizeof(T)));
    mySize += count;

    trivialConstructRange(myData + beg_index, count);

    return beg_index;
}

template <typename T>
template <typename S>
exint
UT_Array<T>::insertImpl(S &&s, exint index)
{
    if (index == mySize)
    {
        // This case avoids an extraneous call to trivialConstructRange()
        // which the compiler may not optimize out.
        (void) appendImpl(std::forward<S>(s));
    }
    else if (index > mySize)
    {
	exint src_i = safeIndex(s);

	bumpCapacity(index + 1);

	trivialConstructRange(myData + mySize, index - mySize);

        if (src_i >= 0)
	    construct(myData[index], std::forward<S>(myData[src_i]));
        else
	    construct(myData[index], std::forward<S>(s));

	mySize = index + 1;
    }
    else // (index < mySize)
    {
	exint src_i = safeIndex(s);

        bumpCapacity(mySize + 1);

        ::memmove((void *)&myData[index+1], (void *)&myData[index],
                  ((mySize-index)*sizeof(T)));

        if (src_i >= index)
            ++src_i;

        if (src_i >= 0)
	    construct(myData[index], std::forward<S>(myData[src_i]));
        else
	    construct(myData[index], std::forward<S>(s));

        ++mySize;
    }

    return index;
}

template <typename T>
exint
UT_Array<T>::removeAt(exint idx)
{
    trivialDestruct(myData[idx]);
    if (idx != --mySize)
    {
	::memmove((void *)&myData[idx], (void *)&myData[idx+1],
	          ((mySize-idx)*sizeof(T)));
    }

    return idx;
}

template <typename T>
void
UT_Array<T>::removeRange(exint begin_i, exint end_i)
{
    UT_ASSERT(begin_i <= end_i);
    UT_ASSERT(end_i <= size());
    if (end_i < size())
    {
	trivialDestructRange(myData + begin_i, end_i - begin_i);
	::memmove((void *)&myData[begin_i], (void *)&myData[end_i],
	          (mySize - end_i)*sizeof(T));
    }
    setSize(mySize - (end_i - begin_i));
}

template <typename T>
void
UT_Array<T>::extractRange(exint begin_i, exint end_i, UT_Array<T>& dest)
{
    UT_ASSERT_P(begin_i >= 0);
    UT_ASSERT_P(begin_i <= end_i);
    UT_ASSERT_P(end_i <= size());
    UT_ASSERT(this != &dest);

    exint nelements = end_i - begin_i;

    // grow the raw array if necessary.
    dest.setCapacityIfNeeded(nelements);

    ::memmove((void*)dest.myData, (void*)&myData[begin_i],
              nelements * sizeof(T));
    dest.mySize = nelements;

    // we just asserted this was true, but just in case
    if (this != &dest)
    {
        if (end_i < size())
        {
            ::memmove((void*)&myData[begin_i], (void*)&myData[end_i],
                      (mySize - end_i) * sizeof(T));
        }
        setSize(mySize - nelements);
    }
}

template <typename T>
void
UT_Array<T>::move(exint srcIdx, exint destIdx, exint howMany)
{
    // Make sure all the parameters are valid.
    if( srcIdx < 0 )
	srcIdx = 0;
    if( destIdx < 0 )
	destIdx = 0;
    // If we are told to move a set of elements that would extend beyond the
    // end of the current array, trim the group.
    if( srcIdx + howMany > size() )
	howMany = size() - srcIdx;
    // If the destIdx would have us move the source beyond the end of the
    // current array, move the destIdx back.
    if( destIdx + howMany > size() )
	destIdx = size() - howMany;
    if( srcIdx != destIdx && howMany > 0 )
    {
	void		**tmp = 0;
	exint	  	savelen;

	savelen = SYSabs(srcIdx - destIdx);
	tmp = (void **)::malloc(savelen*sizeof(T));
	if( srcIdx > destIdx && howMany > 0 )
	{
	    // We're moving the group backwards. Save all the stuff that
	    // we would overwrite, plus everything beyond that to the
	    // start of the source group. Then move the source group, then
	    // tack the saved data onto the end of the moved group.
	    ::memcpy(tmp, (void *)&myData[destIdx],  (savelen*sizeof(T)));
	    ::memmove((void *)&myData[destIdx], (void *)&myData[srcIdx],
	              (howMany*sizeof(T)));
	    ::memcpy((void *)&myData[destIdx+howMany], tmp, (savelen*sizeof(T)));
	}
	if( srcIdx < destIdx && howMany > 0 )
	{
	    // We're moving the group forwards. Save from the end of the
	    // group being moved to the end of the where the destination
	    // group will end up. Then copy the source to the destination.
	    // Then move back up to the original source location and drop
	    // in our saved data.
	    ::memcpy(tmp, (void *)&myData[srcIdx+howMany],  (savelen*sizeof(T)));
	    ::memmove((void *)&myData[destIdx], (void *)&myData[srcIdx],
	              (howMany*sizeof(T)));
	    ::memcpy((void *)&myData[srcIdx], tmp, (savelen*sizeof(T)));
	}
	::free(tmp);
    }
}

template <typename T>
template <typename IsEqual>
exint
UT_Array<T>::removeIf(IsEqual is_equal)
{
    // Move dst to the first element to remove.
    exint dst;
    for (dst = 0; dst < mySize; dst++)
    {
	if (is_equal(myData[dst]))
	    break;
    }
    // Now start looking at all the elements past the first one to remove.
    for (exint idx = dst+1; idx < mySize; idx++)
    {
	if (!is_equal(myData[idx]))
	{
	    UT_ASSERT(idx != dst);
	    myData[dst] = myData[idx];
	    dst++;
	}
	// On match, ignore.
    }
    // New size
    mySize = dst;
    return mySize;
}

template <typename T>
void
UT_Array<T>::cycle(exint howMany)
{
    char	*tempPtr;
    exint	 numShift;	//  The number of items we shift
    exint   	 remaining;	//  mySize - numShift

    if (howMany == 0 || mySize < 1) return;

    numShift = howMany % (exint)mySize;
    if (numShift < 0) numShift += mySize;
    remaining = mySize - numShift;
    tempPtr = new char[numShift*sizeof(T)];

    ::memmove(tempPtr, (void *)&myData[remaining], (numShift * sizeof(T)));
    ::memmove((void *)&myData[numShift], (void *)&myData[0], (remaining * sizeof(T)));
    ::memmove((void *)&myData[0], tempPtr, (numShift * sizeof(T)));

    delete [] tempPtr;
}

template <typename T>
void
UT_Array<T>::constant(const T &value)
{
    for (exint i = 0; i < mySize; i++)
    {
	myData[i] = value;
    }
}

template <typename T>
void
UT_Array<T>::zero()
{
    if (isPOD())
	::memset((void *)myData, 0, mySize*sizeof(T));
    else
	trivialConstructRange(myData, mySize);
}

template <typename T>
void		
UT_Array<T>::setCapacity(exint capacity)
{
    // Do nothing when new capacity is the same as the current
    if (capacity == myCapacity)
	return;

    // Special case for non-heap buffers
    if (!isHeapBuffer())
    {
	if (capacity < mySize)
	{
	    // Destroy the extra elements without changing myCapacity
	    trivialDestructRange(myData + capacity, mySize - capacity);
	    mySize = capacity;
	}
	else if (capacity > myCapacity)
	{
	    T *prev = myData;
	    myData = (T *)malloc(sizeof(T) * capacity);
	    // myData is safe because we're already a stack buffer
	    UT_ASSERT_P(isHeapBuffer());
	    if (mySize > 0)
		memcpy((void *)myData, (void *)prev, sizeof(T) * mySize);
	    myCapacity = capacity;
	}
	else 
	{
	    // Keep myCapacity unchanged in this case
	    UT_ASSERT_P(capacity >= mySize && capacity <= myCapacity);
	}
	return;
    }

    if (capacity == 0)
    {
	if (myData)
	{
	    trivialDestructRange(myData, mySize);
	    free(myData);
	}
	myData     = 0;
	myCapacity = 0;
        mySize = 0;
	return;
    }

    if (capacity < mySize)
    {
	trivialDestructRange(myData + capacity, mySize - capacity);
	mySize = capacity;
    }

    if (myData)
	myData = (T *)realloc(myData, capacity*sizeof(T));
    else
	myData = (T *)malloc(sizeof(T) * capacity);

    // Avoid degenerate case if we happen to be aliased the wrong way
    if (!isHeapBuffer())
    {
	T *prev = myData;
	myData = (T *)malloc(sizeof(T) * capacity);
	if (mySize > 0)
	    memcpy((void *)myData, (void *)prev, sizeof(T) * mySize);
	ut_ArrayImplFree(prev);
    }

    myCapacity = capacity;
    UT_ASSERT(myData);
}

template <typename T>
UT_Array<T> &
UT_Array<T>::operator=(const UT_Array<T> &a)
{
    if (this == &a)
	return *this;

    // Grow the raw array if necessary.
    setCapacityIfNeeded(a.size());

    // Make sure destructors and constructors are called on all elements
    // being removed/added.
    trivialDestructRange(myData, mySize);
    copyConstructRange(myData, a.myData, a.size());

    mySize = a.size();

    return *this;
}

template <typename T>
UT_Array<T> &
UT_Array<T>::operator=(std::initializer_list<T> a)
{
    const exint new_size = a.size();

    // Grow the raw array if necessary.
    setCapacityIfNeeded(new_size);

    // Make sure destructors and constructors are called on all elements
    // being removed/added.
    trivialDestructRange(myData, mySize);

    copyConstructRange(myData, a.begin(), new_size);

    mySize = new_size;

    return *this;
}

template <typename T>
UT_Array<T> &
UT_Array<T>::operator=(UT_Array<T> &&a)
{
    if (!a.isHeapBuffer())
    {
	// Cannot steal from non-heap buffers
	clear();
	const exint n = a.size();
	setCapacityIfNeeded(n);
	if (isPOD())
	{
	    if (n > 0)
		memcpy(myData, a.myData, n * sizeof(T));
	}
	else
	{
	    for (exint i = 0; i < n; ++i)
		new (&myData[i]) T(std::move(a.myData[i]));
	}
	mySize = a.mySize;
	a.mySize = 0;
	return *this;
    }
    // else, just steal even if we're a small buffer

    // Destroy all the elements we're currently holding.
    if (myData)
    {
	trivialDestructRange(myData, mySize);
	if (isHeapBuffer())
	    ::free(myData);
    }
    
    // Move the contents of the other array to us and empty the other container
    // so that it destructs cleanly.
    myCapacity = a.myCapacity;
    mySize = a.mySize;
    myData = a.myData;
    a.myCapacity = a.mySize = 0;
    a.myData = nullptr;

    return *this;
}


template <typename T>
bool
UT_Array<T>::operator==(const UT_Array<T> &a) const
{
    if (this == &a) return true;
    if (mySize != a.size()) return false;
    for (exint i = 0; i < mySize; i++)
	if (!(myData[i] == a(i))) return false;
    return true;
}

template <typename T>
bool
UT_Array<T>::operator!=(const UT_Array<T> &a) const
{
    return (!operator==(a));
}

#endif // __UT_ARRAYIMPL_H_INCLUDED__
