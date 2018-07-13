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
 *      A wrapper function for the "free" function, used by UT_(Small)Array
 */

#include "UT_Array.h"
#include "UT_SmallArray.h"
#include <stdlib.h>

// This needs to be here or else the warning suppression doesn't work because
// the templated calling code won't otherwise be compiled until after we've
// already popped the warning.state. So we just always disable this at file
// scope here.
#if defined(__GNUC__) && !defined(__clang__)
    _Pragma("GCC diagnostic push")
    _Pragma("GCC diagnostic ignored \"-Wfree-nonheap-object\"")
#endif
void ut_ArrayImplFree(void *p)
{
    free(p);
}
#if defined(__GNUC__) && !defined(__clang__)
    _Pragma("GCC diagnostic pop")
#endif
