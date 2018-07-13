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
 *      This is the array class implementation used by almost everything here.
 */

#pragma once

#ifndef __UT_ARRAY_H_INCLUDED__
#define __UT_ARRAY_H_INCLUDED__

#include "SYS_Types.h"

#include <algorithm>
#include <functional>
#include <type_traits>
#include <string.h>

 /// This routine describes how to change the size of an array.
 /// It must increase the current_size by at least one!
 ///
 /// Current expected sequence of small sizes:
 ///    4,    8,   16,   32,   48,   64,   80,   96,  112,
 ///  128,  256,  384,  512,  640,  768,  896, 1024,
 /// (increases by approx factor of 1.125 each time after this)
template <typename T>
static inline T
UTbumpAlloc(T current_size)
{
    // NOTE: These must be powers of two.  See below.
    constexpr T SMALL_ALLOC(16);
    constexpr T BIG_ALLOC(128);

    // For small values, we increment by fixed amounts.  For
    // large values, we increment by one eighth of the current size.
    // This prevents n^2 behaviour with allocation one element at a time.
    // A factor of 1/8 will waste 1/16 the memory on average, and will
    // double the size of the array in approximately 6 reallocations.
    if (current_size < T(8))
    {
        return (current_size < T(4)) ? T(4) : T(8);
    }
    if (current_size < T(BIG_ALLOC))
    {
        // Snap up to next multiple of SMALL_ALLOC (must be power of 2)
        return (current_size + T(SMALL_ALLOC)) & ~T(SMALL_ALLOC-1);
    }
    if (current_size < T(BIG_ALLOC * 8))
    {
        // Snap up to next multiple of BIG_ALLOC (must be power of 2)
        return (current_size + T(BIG_ALLOC)) & ~T(BIG_ALLOC-1);
    }

    T bump = current_size >> 3; // Divided by 8.
    current_size += bump;
    return current_size;
}

template <typename T>
class UT_Array
{
public:
    typedef T value_type;

    typedef int (*Comparator)(const T *, const T *);

    /// Copy constructor. It duplicates the data.
    /// It's marked explicit so that it's not accidentally passed by value.
    /// You can always pass by reference and then copy it, if needed.
    /// If you have a line like:
    /// UT_Array<int> a = otherarray;
    /// and it really does need to copy instead of referencing,
    /// you can rewrite it as:
    /// UT_Array<int> a(otherarray);
    explicit UT_Array(const UT_Array<T> &a);
    
    /// Move constructor. Steals the working data from the original.
    UT_Array(UT_Array<T> &&a) noexcept;
    
    /// Construct based on given capacity and size
    UT_Array(exint capacity, exint size)
    {
        myData = capacity ? allocateCapacity(capacity) : NULL;
        if (capacity < size)
            size = capacity;
        mySize = size;
        myCapacity = capacity;
        trivialConstructRange(myData, mySize);
    }

    /// Construct based on given capacity with a size of 0
    explicit UT_Array(exint capacity = 0) : myCapacity(capacity), mySize(0)
    {
        myData = capacity ? allocateCapacity(capacity) : NULL;
    }

    /// Construct with the contents of an initializer list
    explicit UT_Array(std::initializer_list<T> init);

    ~UT_Array();

    void	    swap(UT_Array<T> &other);

    /// Append an element to the current elements and return its index in the
    /// array, or insert the element at a specified position; if necessary,
    /// insert() grows the array to accommodate the element. The insert
    /// methods use the assignment operator '=' to place the element into the
    /// right spot; be aware that '=' works differently on objects and pointers.
    /// The test for duplicates uses the logical equal operator '=='; as with
    /// '=', the behaviour of the equality operator on pointers versus objects
    /// is not the same.
    /// Use the subscript operators instead of insert() if you are appending
    /// to the array, or if  you don't mind overwriting the element already 
    /// inserted at the given index.
    exint           append(void) { return insert(mySize); }
    exint           append(const T &t) { return appendImpl(t); }
    exint           append(T &&t) { return appendImpl(std::move(t)); }
    void            append(const T *pt, exint count);
    void	    appendMultiple(const T &t, exint count);
    exint	    insert(exint index);
    exint	    insert(const T &t, exint i)
                        { return insertImpl(t, i); }
    exint	    insert(T &&t, exint i)
                        { return insertImpl(std::move(t), i); }

    /// Adds a new element to the array (resizing if necessary) and forwards
    /// the given arguments to T's constructor.
    /// NOTE: Unlike append(), the arguments cannot reference any existing
    /// elements in the array. Checking for and handling such cases would
    /// remove most of the performance gain versus append(T(...)). Debug builds
    /// will assert that the arguments are valid.
    template <typename... S>
    exint	    emplace_back(S&&... s);

    /// Takes another T array and concatenate it onto my end
    exint	    concat(const UT_Array<T> &a);

    /// Insert an element "count" times at the given index. Return the index.
    exint	    multipleInsert(exint index, exint count);

    /// An alias for unique element insertion at a certain index. Also used by
    /// the other insertion methods.
    exint	    insertAt(const T &t, exint index)
                        { return insertImpl(t, index); }

    /// Return true if given index is valid.
    bool	    isValidIndex(exint index) const
			{ return (index >= 0 && index < mySize); }

    /// Remove one element from the array given its
    /// position in the list, and fill the gap by shifting the elements down
    /// by one position.  Return the index of the element removed or -1 if
    /// the index was out of bounds.
    exint	    removeIndex(exint index)
    {
        return isValidIndex(index) ? removeAt(index) : -1;
    }
    void	    removeLast()
    {
        if (mySize) removeAt(mySize-1);
    }

    /// Remove the range [begin_i,end_i) of elements from the array.
    void	    removeRange(exint begin_i, exint end_i);

    /// Remove the range [begin_i, end_i) of elements from this array and place
    /// them in the dest array, shrinking/growing the dest array as necessary.
    void            extractRange(exint begin_i, exint end_i,
                                 UT_Array<T>& dest);

    /// Removes all matching elements from the list, shuffling down and changing
    /// the size appropriately.
    /// Returns the number of elements left.
    template <typename IsEqual>
    exint	    removeIf(IsEqual is_equal);

    /// Remove all matching elements. Also sets the capacity of the array.
    template <typename IsEqual>
    void	    collapseIf(IsEqual is_equal)
    {
        removeIf(is_equal);
        setCapacity(size());
    }

    /// Move howMany objects starting at index srcIndex to destIndex;
    /// This method will remove the elements at [srcIdx, srcIdx+howMany) and
    /// then insert them at destIdx.  This method can be used in place of
    /// the old shift() operation.
    void	    move(exint srcIdx, exint destIdx, exint howMany);

    /// Cyclically shifts the entire array by howMany
    void	    cycle(exint howMany);

    /// Quickly set the array to a single value.
    void	    constant(const T &v);
    /// Zeros the array if a POD type, else trivial constructs if a class type.
    void	    zero();

    /// The fastest search possible, which does pointer arithmetic to find the
    /// index of the element. WARNING: index() does no out-of-bounds checking.
    exint	    index(const T &t) const { return &t - myData; }
    exint	    safeIndex(const T &t) const
    {
        return (&t >= myData && &t < (myData + mySize))
            ? &t - myData : -1;
    }

    /// Set the capacity of the array, i.e. grow it or shrink it. The
    /// function copies the data after reallocating space for the array.
    void            setCapacity(exint newcapacity);
    void            setCapacityIfNeeded(exint mincapacity)
    {
        if (capacity() < mincapacity)
            setCapacity(mincapacity);
    }
    /// If the capacity is smaller than mincapacity, expand the array
    /// to at least mincapacity and to at least a constant factor of the
    /// array's previous capacity, to avoid having a linear number of
    /// reallocations in a linear number of calls to bumpCapacity.
    void            bumpCapacity(exint mincapacity)
    {
        if (capacity() >= mincapacity)
            return;
        // The following 4 lines are just
        // SYSmax(mincapacity, UTbumpAlloc(capacity())), avoiding SYSmax
        exint bumped = UTbumpAlloc(capacity());
        exint newcapacity = mincapacity;
        if (bumped > mincapacity)
            newcapacity = bumped;
        setCapacity(newcapacity);
    }

    /// First bumpCapacity to ensure that there's space for newsize,
    /// expanding either not at all or by at least a constant factor
    /// of the array's previous capacity,
    /// then set the size to newsize.
    void            bumpSize(exint newsize)
    {
        bumpCapacity(newsize);
        setSize(newsize);
    }
    /// NOTE: bumpEntries() will be deprecated in favour of bumpSize() in a
    ///       future version.
    void            bumpEntries(exint newsize)
    {
        bumpSize(newsize);
    }

    /// Query the capacity, i.e. the allocated length of the array.
    /// NOTE: capacity() >= size().
    exint           capacity() const { return myCapacity; }
    /// Query the size, i.e. the number of occupied elements in the array.
    /// NOTE: capacity() >= size().
    exint           size() const     { return mySize; }
    /// Alias of size().  size() is preferred.
    exint           entries() const  { return mySize; }
    /// Returns true iff there are no occupied elements in the array.
    bool            isEmpty() const  { return mySize==0; }

    /// Set the size, the number of occupied elements in the array.
    /// NOTE: This will not do bumpCapacity, so if you call this
    ///       n times to increase the size, it may take
    ///       n^2 time.
    void            setSize(exint newsize)
    {
        if (newsize < 0)
            newsize = 0;
        if (newsize == mySize)
            return;
        setCapacityIfNeeded(newsize);
        if (mySize > newsize)
            trivialDestructRange(myData + newsize, mySize - newsize);
        else // newsize > mySize
            trivialConstructRange(myData + mySize, newsize - mySize);
        mySize = newsize;
    }
    /// Alias of setSize().  setSize() is preferred.
    void            entries(exint newsize)
    {
        setSize(newsize);
    }
    /// Set the size, but unlike setSize(newsize), this function
    /// will not initialize new POD elements to zero. Non-POD data types
    /// will still have their constructors called.
    /// This function is faster than setSize(ne) if you intend to fill in
    /// data for all elements.
    void            setSizeNoInit(exint newsize)
    {
        if (newsize < 0)
            newsize = 0;
        if (newsize == mySize)
            return;
        setCapacityIfNeeded(newsize);
        if (mySize > newsize)
            trivialDestructRange(myData + newsize, mySize - newsize);
        else if (!isPOD()) // newsize > mySize
            trivialConstructRange(myData + mySize, newsize - mySize);
        mySize = newsize;
    }

    /// Decreases, but never expands, to the given maxsize.
    void            truncate(exint maxsize)
    {
        if (maxsize >= 0 && size() > maxsize)
            setSize(maxsize);
    }
    /// Resets list to an empty list.
    void            clear() {
        // Don't call setSize(0) since that would require a valid default
        // constructor.
        trivialDestructRange(myData, mySize);
        mySize = 0;
    }

    /// Assign array a to this array by copying each of a's elements with
    /// memcpy for POD types, and with copy construction for class types.
    UT_Array<T> &   operator=(const UT_Array<T> &a);

    /// Replace the contents with those from the initializer_list ilist
    UT_Array<T> &   operator=(std::initializer_list<T> ilist);

    /// Move the contents of array a to this array. 
    UT_Array<T> &   operator=(UT_Array<T> &&a);

    /// Compare two array and return true if they are equal and false otherwise.
    /// Two elements are checked against each other using operator '==' or
    /// compare() respectively.
    /// NOTE: The capacities of the arrays are not checked when
    ///       determining whether they are equal.
    bool            operator==(const UT_Array<T> &a) const;
    bool            operator!=(const UT_Array<T> &a) const;
     
    /// Subscript operator
    /// NOTE: This does NOT do any bounds checking unless paranoid
    ///       asserts are enabled.
    T &		    operator()(exint i)
    {
        UT_ASSERT_P(i >= 0 && i < mySize);
        return myData[i];
    }
    /// Const subscript operator
    /// NOTE: This does NOT do any bounds checking unless paranoid
    ///       asserts are enabled.
    const T &	    operator()(exint i) const
    {
	UT_ASSERT_P(i >= 0 && i < mySize);
	return myData[i];
    }

    /// Subscript operator
    /// NOTE: This does NOT do any bounds checking unless paranoid
    ///       asserts are enabled.
    T &		    operator[](exint i)
    {
        UT_ASSERT_P(i >= 0 && i < mySize);
        return myData[i];
    }
    /// Const subscript operator
    /// NOTE: This does NOT do any bounds checking unless paranoid
    ///       asserts are enabled.
    const T &	    operator[](exint i) const
    {
	UT_ASSERT_P(i >= 0 && i < mySize);
	return myData[i];
    }
    
    /// forcedRef(exint) will grow the array if necessary, initializing any
    /// new elements to zero for POD types and default constructing for
    /// class types.
    T &             forcedRef(exint i)
    {
        UT_ASSERT_P(i >= 0);
        if (i >= mySize)
            bumpSize(i+1);
        return myData[i];
    }

    /// forcedGet(exint) does NOT grow the array, and will return default
    /// objects for out of bound array indices.
    T               forcedGet(exint i) const
    {
        return (i >= 0 && i < mySize) ? myData[i] : T();
    }

    T &		    last()
    {
        UT_ASSERT_P(mySize);
        return myData[mySize-1];
    }
    const T &	    last() const
    {
        UT_ASSERT_P(mySize);
        return myData[mySize-1];
    }

    T *		    getArray() const		    { return myData; }
    const T *	    getRawArray() const		    { return myData; }

    T *		    array()			    { return myData; }
    const T *	    array() const		    { return myData; }

    T *		    data()			    { return myData; }
    const T *	    data() const		    { return myData; }

    /// This method allows you to swap in a new raw T array, which must be
    /// the same size as myCapacity. Use caution with this method.
    T *		    aliasArray(T *newdata)
    { T *data = myData; myData = newdata; return data; }

    template <typename IT, bool FORWARD>
    class base_iterator : 
	public std::iterator<std::random_access_iterator_tag, T, exint> 
    {
        public:
	    typedef IT&		reference;
	    typedef IT*		pointer;
	
	    // Note: When we drop gcc 4.4 support and allow range-based for
	    // loops, we should also drop atEnd(), which means we can drop
	    // myEnd here.
	    base_iterator() : myCurrent(NULL), myEnd(NULL) {}
	    
	      // Allow iterator to const_iterator conversion
	    template<typename EIT>
	    base_iterator(const base_iterator<EIT, FORWARD> &src)
		: myCurrent(src.myCurrent), myEnd(src.myEnd) {}
	    
	    pointer	operator->() const 
			{ return FORWARD ? myCurrent : myCurrent - 1; }
	    
	    reference	operator*() const
			{ return FORWARD ? *myCurrent : myCurrent[-1]; }

	    reference	item() const
			{ return FORWARD ? *myCurrent : myCurrent[-1]; }
	    
	    reference	operator[](exint n) const
			{ return FORWARD ? myCurrent[n] : myCurrent[-n - 1]; } 

	    /// Pre-increment operator
            base_iterator &operator++()
			{
        		    if (FORWARD) ++myCurrent; else --myCurrent;
			    return *this;
			}
	    /// Post-increment operator
	    base_iterator operator++(int)
			{
			    base_iterator tmp = *this;
        		    if (FORWARD) ++myCurrent; else --myCurrent;
        		    return tmp;
			}
	    /// Pre-decrement operator
	    base_iterator &operator--()
			{
			    if (FORWARD) --myCurrent; else ++myCurrent;
			    return *this;
			}
	    /// Post-decrement operator
	    base_iterator operator--(int)
			{
			    base_iterator tmp = *this;
			    if (FORWARD) --myCurrent; else ++myCurrent;
			    return tmp;
			}

	    base_iterator &operator+=(exint n)   
			{
			    if (FORWARD)
				myCurrent += n;
			    else
				myCurrent -= n;
			    return *this;
			}
            base_iterator operator+(exint n) const
			{
			    if (FORWARD)
				return base_iterator(myCurrent + n, myEnd);
			    else
				return base_iterator(myCurrent - n, myEnd);
			}
	    
            base_iterator &operator-=(exint n)
        		{ return (*this) += (-n); }
            base_iterator operator-(exint n) const
			{ return (*this) + (-n); }
            
	    bool	 atEnd() const		{ return myCurrent == myEnd; }
	    void	 advance()		{ this->operator++(); }
	    
	    // Comparators
	    template<typename ITR, bool FR>
	    bool 	 operator==(const base_iterator<ITR, FR> &r) const
			 { return myCurrent == r.myCurrent; }
	    
	    template<typename ITR, bool FR>
	    bool 	 operator!=(const base_iterator<ITR, FR> &r) const
			 { return myCurrent != r.myCurrent; }
	    
	    template<typename ITR>
	    bool	 operator<(const base_iterator<ITR, FORWARD> &r) const
	    {
		if (FORWARD) 
		    return myCurrent < r.myCurrent;
		else
		    return r.myCurrent < myCurrent;
	    }
	    
	    template<typename ITR>
	    bool	 operator>(const base_iterator<ITR, FORWARD> &r) const
	    {
		if (FORWARD) 
		    return myCurrent > r.myCurrent;
		else
		    return r.myCurrent > myCurrent;
	    }

	    template<typename ITR>
	    bool	 operator<=(const base_iterator<ITR, FORWARD> &r) const
	    {
		if (FORWARD) 
		    return myCurrent <= r.myCurrent;
		else
		    return r.myCurrent <= myCurrent;
	    }

	    template<typename ITR>
	    bool	 operator>=(const base_iterator<ITR, FORWARD> &r) const
	    {
		if (FORWARD) 
		    return myCurrent >= r.myCurrent;
		else
		    return r.myCurrent >= myCurrent;
	    }
	    
	    // Difference operator for std::distance
	    template<typename ITR>
	    exint	 operator-(const base_iterator<ITR, FORWARD> &r) const
	    {
		if (FORWARD) 
		    return exint(myCurrent - r.myCurrent);
		else
		    return exint(r.myCurrent - myCurrent);
	    }
	    
	    
        protected:
	    friend class UT_Array<T>;
	    base_iterator(IT *c, IT *e) : myCurrent(c), myEnd(e) {}
	private:

	    IT			*myCurrent;
	    IT			*myEnd;
    };
    
    typedef base_iterator<T, true>		iterator;
    typedef base_iterator<const T, true>	const_iterator;
    typedef base_iterator<T, false>		reverse_iterator;
    typedef base_iterator<const T, false>	const_reverse_iterator;
    typedef const_iterator	traverser; // For backward compatibility

    /// Begin iterating over the array.  The contents of the array may be 
    /// modified during the traversal.
    iterator		begin()
			{
			    return iterator(myData, myData + mySize);
			}
    /// End iterator.
    iterator		end()
			{
			    return iterator(myData + mySize,
					    myData + mySize);
			}

    /// Begin iterating over the array.  The array may not be modified during
    /// the traversal.
    const_iterator	begin() const
			{
			    return const_iterator(myData, myData + mySize);
			}
    /// End const iterator.  Consider using it.atEnd() instead.
    const_iterator	end() const
			{
			    return const_iterator(myData + mySize,
						  myData + mySize);
			}
    
    /// Begin iterating over the array in reverse. 
    reverse_iterator	rbegin()
			{
			    return reverse_iterator(myData + mySize,
			                            myData);
			}
    /// End reverse iterator.
    reverse_iterator	rend()
			{
			    return reverse_iterator(myData, myData);
			}
    /// Begin iterating over the array in reverse. 
    const_reverse_iterator rbegin() const
			{
			    return const_reverse_iterator(myData + mySize,
							  myData);
			}
    /// End reverse iterator.  Consider using it.atEnd() instead.
    const_reverse_iterator rend() const
			{
			    return const_reverse_iterator(myData, myData);
			}
    
    /// Remove item specified by the reverse_iterator.
    void		removeItem(const reverse_iterator &it)
			{
			    removeAt(&it.item() - myData);
			}


    /// Very dangerous methods to share arrays.
    /// The array is not aware of the sharing, so ensure you clear
    /// out the array prior a destructor or setCapacity operation.
    void	    unsafeShareData(UT_Array<T> &src)
		    {
			myData = src.myData;
			myCapacity = src.myCapacity;
			mySize = src.mySize;
		    }
    void	    unsafeShareData(T *src, exint srcsize)
		    {
			myData = src;
			myCapacity = srcsize;
			mySize = srcsize;
		    }
    void	    unsafeShareData(T *src, exint size, exint capacity)
		    {
			myData = src;
			mySize = size;
			myCapacity = capacity;
		    }
    void	    unsafeClearData()
		    {
			myData = NULL;
			myCapacity = 0;
			mySize = 0;
		    }

    /// Returns true if the data used by the array was allocated on the heap.
    inline bool	    isHeapBuffer() const
		    {
			return (myData != (T *)(((char*)this) + sizeof(*this)));
		    }
    inline bool	    isHeapBuffer(T* data) const
		    {
			return (data != (T *)(((char*)this) + sizeof(*this)));
		    }

protected:
    // Check whether T may have a constructor, destructor, or copy
    // constructor.  This test is conservative in that some POD types will
    // not be recognized as POD by this function. To mark your type as POD,
    // use the SYS_DECLARE_IS_POD() macro in SYS_TypeDecorate.h.
    static constexpr SYS_FORCE_INLINE bool isPOD()
    {
        return std::is_pod<T>::value;
    }

    /// Implements both append(const T &) and append(T &&) via perfect
    /// forwarding. Unlike the variadic emplace_back(), its argument may be a
    /// reference to another element in the array.
    template <typename S>
    exint           appendImpl(S &&s);

    /// Similar to appendImpl() but for insertion.
    template <typename S>
    exint           insertImpl(S &&s, exint index);

    // Construct the given type
    template <typename... S>
    static void	    construct(T &dst, S&&... s)
		    {
                        new (&dst) T(std::forward<S>(s)...);
		    }

    // Copy construct the given type
    static void	    copyConstruct(T &dst, const T &src)
		    {
			if (isPOD())
			    dst = src;
			else
			    new (&dst) T(src);
		    }
    static void	    copyConstructRange(T *dst, const T *src, exint n)
		    {
			if (isPOD())
                        {
                            if (n > 0)
                            {
                                ::memcpy((void *)dst, (const void *)src,
                                         n * sizeof(T));
                            }
                        }
			else
			{
			    for (exint i = 0; i < n; i++)
				new (&dst[i]) T(src[i]);
			}
		    }

    /// Element Constructor
    static void	    trivialConstruct(T &dst)
		    {
			if (!isPOD())
			    new (&dst) T();
			else
			    memset((void *)&dst, 0, sizeof(T));
		    }
    static void	    trivialConstructRange(T *dst, exint n)
		    {
			if (!isPOD())
			{
			    for (exint i = 0; i < n; i++)
				new (&dst[i]) T();
			}
			else if (n == 1)
			{
			    // Special case for n == 1. If the size parameter
			    // passed to memset is known at compile time, this
			    // function call will be inlined. This results in
			    // much faster performance than a real memset
			    // function call which is required in the case
			    // below, where n is not known until runtime.
			    // This makes calls to append() much faster.
			    memset((void *)dst, 0, sizeof(T));
			}
			else
			    memset((void *)dst, 0, sizeof(T) * n);
		    }

    /// Element Destructor
    static void	    trivialDestruct(T &dst)
		    {
			if (!isPOD())
			    dst.~T();
		    }
    static void	    trivialDestructRange(T *dst, exint n)
		    {
			if (!isPOD())
			{
			    for (exint i = 0; i < n; i++)
				dst[i].~T();
			}
		    }

private:
    /// Pointer to the array of elements of type T
    T *myData;

    /// The number of elements for which we have allocated memory
    exint myCapacity;

    /// The actual number of valid elements in the array
    exint mySize;

    // The guts of the remove() methods.
    exint	    removeAt(exint index);

    T *		    allocateCapacity(exint num_items);
};

#include "UT_ArrayImpl.h"

#endif // __UT_ARRAY_H_INCLUDED__
