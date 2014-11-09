/*
 * karect version 1.0
 * Copyright (c) 2014 KAUST All Rights Reserved.
 * Author: Amin Allam
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __MEMORY_UTILITIES
#define __MEMORY_UTILITIES

inline void* AllocateHeap(long long n)
{
	return malloc(n);
}

inline void FreeHeap(void* a, long long n)
{
	free(a);
}

void PrintMemory(long long v)
{
	printf("%lld[", v);
	int GBs=v/(1024*1024*1024); v%=(1024*1024*1024);
	int MBs=v/(1024*1024); v%=(1024*1024);
	int KBs=v/1024; v%=1024;
	int Bs=v;
	int first=1;
	if(GBs) {if(!first) printf(":"); printf("%dg", GBs); first=0;}
	if(MBs) {if(!first) printf(":"); printf("%dm", MBs); first=0;}
	if(KBs) {if(!first) printf(":"); printf("%dk", KBs); first=0;}
	if(Bs) {if(!first) printf(":"); printf("%db", Bs); first=0;}
	printf("]");
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
