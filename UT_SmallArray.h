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
 *      Special case for arrays that are usually small,
 *      to avoid a heap allocation when the array really is small.
 */

#pragma once

#ifndef __UT_SMALLARRAY_H_INCLUDED__
#define __UT_SMALLARRAY_H_INCLUDED__

#include "UT_Array.h"
#include "SYS_Types.h"
#include <utility>
#include <stddef.h>

/// An array class with the small buffer optimization, making it ideal for
/// cases when you know it will only contain a few elements at the expense of
/// increasing the object size by MAX_BYTES (subject to alignment).
template <typename T, size_t MAX_BYTES = 64>
class UT_SmallArray : public UT_Array<T>
{
    // As many elements that fit into MAX_BYTES with 1 item minimum
    enum { MAX_ELEMS = MAX_BYTES/sizeof(T) < 1 ? 1 : MAX_BYTES/sizeof(T) };

public:

// gcc falsely warns about our use of offsetof() on non-POD types. We can't
// easily suppress this because it has to be done in the caller at
// instantiation time. Instead, punt to a runtime check instead.
#if defined(__clang__) || defined(_MSC_VER)
    #define UT_SMALL_ARRAY_SIZE_ASSERT()    \
        using ThisT = UT_SmallArray<T,MAX_BYTES>; \
	static_assert(offsetof(ThisT, myBuffer) == sizeof(UT_Array<T>), \
            "In order for UT_Array's checks for whether it needs to free the buffer to work, " \
            "the buffer must be exactly following the base class memory.")
#else
    #define UT_SMALL_ARRAY_SIZE_ASSERT()    \
	UT_ASSERT_P(!UT_Array<T>::isHeapBuffer());
#endif

    /// Default construction
    UT_SmallArray()
	: UT_Array<T>(/*capacity*/0)
    {
	UT_Array<T>::unsafeShareData((T*)myBuffer, 0, MAX_ELEMS);
	UT_SMALL_ARRAY_SIZE_ASSERT();
    }
    
    /// Copy constructor
    /// @{
    explicit UT_SmallArray(const UT_Array<T> &copy)
	: UT_Array<T>(/*capacity*/0)
    {
	UT_Array<T>::unsafeShareData((T*)myBuffer, 0, MAX_ELEMS);
	UT_SMALL_ARRAY_SIZE_ASSERT();
	UT_Array<T>::operator=(copy);
    }
    explicit UT_SmallArray(const UT_SmallArray<T,MAX_BYTES> &copy)
	: UT_Array<T>(/*capacity*/0)
    {
	UT_Array<T>::unsafeShareData((T*)myBuffer, 0, MAX_ELEMS);
	UT_SMALL_ARRAY_SIZE_ASSERT();
	UT_Array<T>::operator=(copy);
    }
    /// @}

    /// Move constructor
    /// @{
    UT_SmallArray(UT_Array<T> &&movable) noexcept
    {
	UT_Array<T>::unsafeShareData((T*)myBuffer, 0, MAX_ELEMS);
	UT_SMALL_ARRAY_SIZE_ASSERT();
	UT_Array<T>::operator=(std::move(movable));
    }
    UT_SmallArray(UT_SmallArray<T,MAX_BYTES> &&movable) noexcept
    {
	UT_Array<T>::unsafeShareData((T*)myBuffer, 0, MAX_ELEMS);
	UT_SMALL_ARRAY_SIZE_ASSERT();
	UT_Array<T>::operator=(std::move(movable));
    }
    /// @}

    /// Initializer list constructor
    explicit UT_SmallArray(std::initializer_list<T> init)
    {
        UT_Array<T>::unsafeShareData((T*)myBuffer, 0, MAX_ELEMS);
        UT_SMALL_ARRAY_SIZE_ASSERT();
        UT_Array<T>::operator=(init);
    }

#undef UT_SMALL_ARRAY_SIZE_ASSERT

    /// Assignment operator
    /// @{
    UT_SmallArray<T,MAX_BYTES> &
    operator=(const UT_SmallArray<T,MAX_BYTES> &copy)
    {
	UT_Array<T>::operator=(copy);
	return *this;
    }
    UT_SmallArray<T,MAX_BYTES> &
    operator=(const UT_Array<T> &copy)
    {
	UT_Array<T>::operator=(copy);
	return *this;
    }
    /// @}

    /// Move operator
    /// @{
    UT_SmallArray<T,MAX_BYTES> &
    operator=(UT_SmallArray<T,MAX_BYTES> &&movable)
    {
	UT_Array<T>::operator=(std::move(movable));
	return *this;
    }
    UT_SmallArray<T,MAX_BYTES> &
    operator=(UT_Array<T> &&movable)
    {
        UT_Array<T>::operator=(std::move(movable));
        return *this;
    }
    /// @}

    UT_SmallArray<T,MAX_BYTES> &
    operator=(std::initializer_list<T> src)
    {
        UT_Array<T>::operator=(src);
        return *this;
    }
private:
    alignas(T) char myBuffer[MAX_ELEMS*sizeof(T)];
};

#endif // __UT_SMALLARRAY_H_INCLUDED__
