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

enum{MATCH_EDIT_DIST, MATCH_SUBSTITUTE_ONLY, MATCH_INSERT_DELETE_ONLY};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum{ED_NONE=0, ED_SUBSTITUTE=1, ED_DELETE=2, ED_INSERT=4};
enum{ED_MATCH=ED_SUBSTITUTE};

struct ShortNode
{
	short	ia;
	char	type;
	char	ch;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define ED_INF 0x7FFFFFFE // ED_INF+1 is still signed integer

// T is a working array of size (nb+1)
int EditDistance(int* T, const char* a, int na, const char* b, int nb, int match_type)
{
	const int nosub=(match_type==MATCH_INSERT_DELETE_ONLY);

	if(match_type==MATCH_SUBSTITUTE_ONLY)
	{
		if(nb<na)
		{
			int bh=nb;
			int j;
			for(j=0;j<na-nb+1;j++)
			{
				int i,h=0; for(i=0;i<nb;i++) h+=(a[i+j]!=b[i]);
				if(h<bh) bh=h;
			}
			return bh+(na-nb);
		}
		if(na!=nb) return na;
		int i,h=0; for(i=0;i<na;i++) h+=(a[i]!=b[i]);
		return h;
	}

	int ia, ib;

	for(ib=0;ib<=nb;ib++)
		T[ib]=ib;

	for(ia=1;ia<=na;ia++)
	{
		int pT=ia-1, nT=ia;
		for(ib=1;ib<=nb;ib++)
		{
			if(b[ib-1]==a[ia-1]) nT=pT;
			else
			{
				if(pT+nosub<nT) nT=pT+nosub;
				if(T[ib]<nT) nT=T[ib];
				nT++;
			}
			pT=T[ib];
			T[ib]=nT;
		}
	}

	return T[nb];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ApproxEditDist(int* T, int* H, int* D, const char* a, int na, const char* b, int nb, int threshold, int match_type)
{
	const int nosub=(match_type==MATCH_INSERT_DELETE_ONLY);
	const int sub_cost=1+nosub;

	int oo=threshold+1;
	int ia, ib;

	int mid_a=na/2;

	int ia_st=0;
	int ia_en=na;

	ia=ia_st;

	int ib_st=ia-threshold; if(ib_st<0) ib_st=0;
	int ib_en=ia+threshold; if(ib_en>nb) ib_en=nb;

	for(ib=ib_st;ib<=ib_en;ib++)
		T[ib]=ib;

	ia++;

	for(;ia<=mid_a;ia++)
	{
		int prev_ib_st=ib_st;
		int prev_ib_en=ib_en;

		ib_st=ia-threshold; if(ib_st<0) ib_st=0;
		ib_en=ia+threshold; if(ib_en>nb) ib_en=nb;

		ib=ib_st;

		int prev_T_ibm1=-1; // T[ca-1][ib-1]

		if(ib_st==0)
		{
			prev_T_ibm1=T[ib];

			T[ib]=ia;
			ib++;
		}
		else
		{
			prev_T_ibm1=T[ib-1];

			int d1=prev_T_ibm1; if(a[ia-1]!=b[ib-1]) d1+=sub_cost;
			int d2=T[ib]+1; if(ib>prev_ib_en) d2=oo;

			int ret=d1; if(d2<ret) ret=d2;
			prev_T_ibm1=T[ib]; T[ib]=ret;

			ib++;
		}

		if(ib-1<prev_ib_st)
			prev_T_ibm1=oo;

		for(;ib<ib_en;ib++)
		{
			int d1=prev_T_ibm1; if(a[ia-1]!=b[ib-1]) d1+=sub_cost;

			int ret=T[ib-1]; if(T[ib]<ret) ret=T[ib]; ++ret; if(d1<ret) ret=d1;

			prev_T_ibm1=T[ib];
			T[ib]=ret;
		}

		ib=ib_en;

		if(ib>ib_st)
		{
			int d1=prev_T_ibm1; if(a[ia-1]!=b[ib-1]) d1+=sub_cost;
			int d2=T[ib]+1; if(ib>prev_ib_en) d2=oo;
			int d3=T[ib-1]+1;

			int ret=d1; if(d2<ret) ret=d2; if(d3<ret) ret=d3;

			prev_T_ibm1=T[ib]; T[ib]=ret;
		}
	}

	for(ib=ib_st;ib<=ib_en;ib++)
	{
		H[ib]=ib;
		D[ib]=T[ib];
	}

	for(;ia<=ia_en;ia++)
	{
		int prev_ib_st=ib_st;
		int prev_ib_en=ib_en;

		ib_st=ia-threshold; if(ib_st<0) ib_st=0;
		ib_en=ia+threshold; if(ib_en>nb) ib_en=nb;

		ib=ib_st;

		int prev_T_ibm1=-1; // T[ca-1][ib-1]
		int prev_H_ibm1=-1;
		int prev_D_ibm1=-1;

		if(ib_st==0)
		{
			prev_T_ibm1=T[ib];
			prev_H_ibm1=H[ib]; prev_D_ibm1=D[ib];

			T[ib]=ia;
			ib++;
		}
		else
		{
			prev_T_ibm1=T[ib-1];
			prev_H_ibm1=H[ib-1]; prev_D_ibm1=D[ib-1];

			int d1=prev_T_ibm1; if(a[ia-1]!=b[ib-1]) d1+=sub_cost;
			int d2=T[ib]+1; if(ib>prev_ib_en) d2=oo;

			int ret=d1, nH=prev_H_ibm1, nD=prev_D_ibm1;
			if(d2<ret) {ret=d2; nH=H[ib]; nD=D[ib];}

			prev_T_ibm1=T[ib]; T[ib]=ret;
			prev_H_ibm1=H[ib]; H[ib]=nH; prev_D_ibm1=D[ib]; D[ib]=nD;

			ib++;
		}

		if(ib-1<prev_ib_st)
			prev_T_ibm1=oo;

		for(;ib<ib_en;ib++)
		{
			int d1=prev_T_ibm1; if(a[ia-1]!=b[ib-1]) d1+=sub_cost;
			int d2=T[ib]+1;
			int d3=T[ib-1]+1;

			int ret=d1, nH=prev_H_ibm1, nD=prev_D_ibm1;
			if(d2<ret) {ret=d2; nH=H[ib]; nD=D[ib];}
			if(d3<ret) {ret=d3; nH=H[ib-1]; nD=D[ib-1];}

			prev_T_ibm1=T[ib];
			T[ib]=ret;

			prev_H_ibm1=H[ib]; H[ib]=nH; prev_D_ibm1=D[ib]; D[ib]=nD;
		}

		ib=ib_en;

		if(ib>ib_st)
		{
			int d1=prev_T_ibm1; if(a[ia-1]!=b[ib-1]) d1+=sub_cost;
			int d2=T[ib]+1; if(ib>prev_ib_en) d2=oo;
			int d3=T[ib-1]+1;

			int ret=d1, nH=prev_H_ibm1, nD=prev_D_ibm1;
			if(d2<ret) {ret=d2; nH=H[ib]; nD=D[ib];}
			if(d3<ret) {ret=d3; nH=H[ib-1]; nD=D[ib-1];}

			prev_T_ibm1=T[ib]; T[ib]=ret;
			prev_H_ibm1=H[ib]; H[ib]=nH; prev_D_ibm1=D[ib]; D[ib]=nD;
		}
	}

	int final_ret=oo;

	if(nb>=ib_st && nb<=ib_en)
		final_ret=T[nb];

	return final_ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ApproxOptAlignNode
{
	int st_a;
	const char* a; int na;
	const char* b; int nb;
	int threshold;
	ApproxOptAlignNode* prev_node; // used for stack
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ApproxGetOptimalAlignment(int* T, int st_a, const char* a, int na, const char* b, int nb,
								int threshold, int match_type, vector<ShortNode>& v)
{
	int i;
	ShortNode sn;

	Stack<ApproxOptAlignNode> stack;
	InitStack(&stack);

	ApproxOptAlignNode* cur_node=(ApproxOptAlignNode*)AllocateHeap(sizeof(ApproxOptAlignNode));

	cur_node->st_a=st_a;
	cur_node->a=a; cur_node->na=na;
	cur_node->b=b; cur_node->nb=nb;
	cur_node->threshold=threshold;
	Push(&stack, cur_node);

	while(1)
	{
		cur_node=Pop(&stack);
		if(!cur_node) break;

		int st_a=cur_node->st_a;
		const char* a=cur_node->a; int na=cur_node->na;
		const char* b=cur_node->b; int nb=cur_node->nb;
		int threshold=cur_node->threshold;

		FreeHeap(cur_node, sizeof(ApproxOptAlignNode));

		if(nb==0) {for(i=0;i<na;i++) {sn.ia=st_a+i; sn.type=ED_DELETE; sn.ch='-'; v.push_back(sn);} continue;}

		if(na==1)
		{
			for(i=0;i<nb;i++) if(a[0]==b[i]) break;
			if(i<nb)
			{
				int m=i;
				for(i=0;i<m;i++) {sn.ia=st_a; sn.type=ED_INSERT; sn.ch=b[i]; v.push_back(sn);}
				for(i=m+1;i<nb;i++) {sn.ia=st_a+1; sn.type=ED_INSERT; sn.ch=b[i]; v.push_back(sn);}
			}
			else
			{
				sn.ia=st_a; sn.type=ED_SUBSTITUTE; sn.ch=b[0]; v.push_back(sn);
				for(i=1;i<nb;i++) {sn.ia=st_a+1; sn.type=ED_INSERT; sn.ch=b[i]; v.push_back(sn);}
			}
			continue;
		}

		int* H=T+(nb+1);
		int* D=H+(nb+1);
		int ret=ApproxEditDist(T, H, D, a, na, b, nb, threshold, match_type);

		int mid_a=na/2;

		if(ret<0)
		{
			printf("Error Approx Ret\n");
			fflush(NULL);
		}

		int mid_b=H[nb];

		if(mid_b<0 || mid_b>nb)
		{
			printf("Error H\n");
			fflush(NULL);
		}

		ApproxOptAlignNode* left_node=(ApproxOptAlignNode*)AllocateHeap(sizeof(ApproxOptAlignNode));
		ApproxOptAlignNode* right_node=(ApproxOptAlignNode*)AllocateHeap(sizeof(ApproxOptAlignNode));

		left_node->st_a=st_a;
		left_node->a=a; left_node->na=mid_a;
		left_node->b=b; left_node->nb=mid_b;
		left_node->threshold=D[nb];

		right_node->st_a=st_a+mid_a;
		right_node->a=a+mid_a; right_node->na=na-mid_a;
		right_node->b=b+mid_b; right_node->nb=nb-mid_b;
		right_node->threshold=T[nb]-D[nb];

		Push(&stack, right_node);
		Push(&stack, left_node);
	}

	DestroyStack(&stack);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// T is a working array of size (nb+1)
int ApproxEditDistanceSearchBInsideA(int* T, const char* a, int na, const char* b, int nb, int k, int match_type, vector<int>* pend_locs) // k = Max allowed ED
{
	const int nosub=(match_type==MATCH_INSERT_DELETE_ONLY);

	int final_ret=k+1;

	int ia, ib;
	int last_ib=k+1;
	if(last_ib>nb) last_ib=nb;

	for(ib=0;ib<=nb;ib++)
		T[ib]=ib;

	for(ia=1;ia<=na;ia++)
	{
		int pT=0, nT=0;
		for(ib=1;ib<=last_ib;ib++)
		{
			if(b[ib-1]==a[ia-1]) nT=pT;
			else
			{
				if(pT+nosub<nT) nT=pT+nosub;
				if(T[ib]<nT) nT=T[ib];
				nT++;
			}
			pT=T[ib];
			T[ib]=nT;
		}
		while(T[last_ib]>k) last_ib--;
		if(last_ib==nb)
		{
			if(T[last_ib]<final_ret) {final_ret=T[last_ib]; if(pend_locs) pend_locs->clear();}
			if(T[last_ib]<=final_ret && pend_locs) pend_locs->push_back(ia);
		}
		else last_ib++;
	}

	if(final_ret>k) return -1;
	return final_ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// T is a working array of size (nb+1)
int ApproxEditDistanceFlexibleEndA(int* T, const char* a, int na, const char* b, int nb, int k, int match_type, vector<int>* pend_locs) // k = Max allowed ED
{
	const int nosub=(match_type==MATCH_INSERT_DELETE_ONLY);

	if(match_type==MATCH_SUBSTITUTE_ONLY)
	{
		int i, ret=0;
		for(i=0;i<na && i<nb;i++)
		{
			ret+=(a[i]!=b[i]);
			if(ret>k) return -1;
		}
		if(pend_locs) {pend_locs->clear(); pend_locs->push_back(na);} // line actually not used
		return ret;
	}

	int final_ret=k+1;

	int ia, ib;
	int last_ib=k+1;
	if(last_ib>nb) last_ib=nb;

	for(ib=0;ib<=nb;ib++)
		T[ib]=ib;

	for(ia=1;ia<=na;ia++)
	{
		int pT=ia-1, nT=ia;
		for(ib=1;ib<=last_ib;ib++)
		{
			if(b[ib-1]==a[ia-1]) nT=pT;
			else
			{
				if(pT+nosub<nT) nT=pT+nosub;
				if(T[ib]<nT) nT=T[ib];
				nT++;
			}
			pT=T[ib];
			T[ib]=nT;
		}
		while(T[last_ib]>k) last_ib--;
		if(last_ib==nb)
		{
			if(T[last_ib]<final_ret) {final_ret=T[last_ib]; if(pend_locs) pend_locs->clear();}
			if(T[last_ib]<=final_ret && pend_locs) pend_locs->push_back(ia);
		}
		else last_ib++;
	}

	if(final_ret>k) return -1;
	return final_ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ApproxEditDist(int* T, const char* a, int na, const char* b, int nb, int shift1, int shift2, int threshold, int match_type,
					bool flex_left, bool flex_right, int& last_a, int& last_b)
{
	const int nosub=(match_type==MATCH_INSERT_DELETE_ONLY);
	const int sub_cost=1+nosub;

	int oo=threshold+1;
	int ia, ib;

	int ia_st=0+shift1-threshold; if(ia_st<0) ia_st=0;
	int ia_en=nb+shift2+threshold; if(ia_en>na) ia_en=na;

	ia=ia_st;

	int ib_st=ia-shift2-threshold; if(ib_st<0) ib_st=0;
	int ib_en=ia-shift1+threshold; if(ib_en>nb) ib_en=nb;

	int final_ret=oo;

	for(ib=ib_st;ib<=ib_en;ib++)
	{
		T[ib]=ib;
		if(flex_left && ia==0) T[ib]=0;
	}

	ib=ib_en;
	if(ib==nb && flex_right && T[ib]<=final_ret)
	{
		last_a=ia; last_b=ib;
		final_ret=T[ib];
	}

	ia++;

	for(;ia<=ia_en;ia++)
	{
		int prev_ib_st=ib_st;
		int prev_ib_en=ib_en;

		ib_st=ia-shift2-threshold; if(ib_st<0) ib_st=0;
		ib_en=ia-shift1+threshold; if(ib_en>nb) ib_en=nb;

		int m=oo; // make it = final_ret
		ib=ib_st;

		int prev_T_ibm1=-1; // T[ca-1][ib-1]

		if(ib_st==0)
		{
			prev_T_ibm1=T[ib];

			T[ib]=ia;
			if(flex_left) T[ib]=0;

			if(T[ib]<m) m=T[ib];

			ib++;
		}
		else
		{
			prev_T_ibm1=T[ib-1];

			int d1=prev_T_ibm1; if(a[ia-1]!=b[ib-1]) d1+=sub_cost;
			int d2=T[ib]+1; if(ib>prev_ib_en) d2=oo;
			int ret=d1; if(d2<ret) ret=d2;
			prev_T_ibm1=T[ib]; T[ib]=ret; if(ret<m) m=ret;

			ib++;
		}

		if(ib-1<prev_ib_st)// || ib-1>prev_ib_en)
			prev_T_ibm1=oo;

		for(;ib<ib_en;ib++)
		{
			int d1=prev_T_ibm1; if(a[ia-1]!=b[ib-1]) d1+=sub_cost;

			int ret=T[ib-1]; if(T[ib]<ret) ret=T[ib]; ret++; if(d1<ret) ret=d1;

			prev_T_ibm1=T[ib];
			T[ib]=ret;

			if(ret<m) m=ret;
		}

		ib=ib_en;

		if(ib>ib_st)
		{
			int d1=prev_T_ibm1; if(a[ia-1]!=b[ib-1]) d1+=sub_cost;//+(1+nosub)*(a[ia-1]!=b[ib-1]);
			int d2=T[ib]+1; if(ib>prev_ib_en) d2=oo;
			int d3=T[ib-1]+1;
			int ret=d1; if(d2<ret) ret=d2; if(d3<ret) ret=d3;
			prev_T_ibm1=T[ib]; T[ib]=ret;
			if(ret<m) m=ret;
		}

		if(ib==nb && flex_right && T[ib]<=final_ret)
		{
			last_a=ia; last_b=ib;
			final_ret=T[ib];
		}

		if(ia==na)
		{
			if(flex_right)
			{
				for(ib=ib_st;ib<=ib_en;ib++)
				{
					if(T[ib]<=final_ret)
					{
						last_a=ia; last_b=ib;
						final_ret=T[ib];
					}
				}
			}
			else if(nb>=ib_st && nb<=ib_en)
			{
				last_a=ia; last_b=nb;
				final_ret=T[nb];
			}
		}

		if(m>threshold && final_ret>threshold)
			return -1;
	}

	if(final_ret>threshold) return -1;
	return final_ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ApproxEditDistBackward(int* T, const char* a, int na, const char* b, int nb, int shift1, int shift2, int threshold, int match_type,
					bool flex_left, bool flex_right, int& first_a, int& first_b)
{
	const int nosub=(match_type==MATCH_INSERT_DELETE_ONLY);
	const int sub_cost=1+nosub;

	int oo=threshold+1;
	int ia, ib;

	int ia_st=0+shift1-threshold; if(ia_st<0) ia_st=0;
	int ia_en=nb+shift2+threshold; if(ia_en>na) ia_en=na;

	ia=ia_en;

	int ib_st=ia-shift2-threshold; if(ib_st<0) ib_st=0;
	int ib_en=ia-shift1+threshold; if(ib_en>nb) ib_en=nb;

	int final_ret=oo;

	for(ib=ib_en;ib>=ib_st;ib--)
	{
		T[ib]=nb-ib;
		if(flex_left && ia==na) T[ib]=0;
	}

	ib=ib_st;

	if(flex_right && ib==0 && T[ib]<=final_ret)
	{
		first_a=ia; first_b=ib;
		final_ret=T[ib];
	}

	ia--;

	for(;ia>=ia_st;ia--)
	{
		int prev_ib_st=ib_st;
		int prev_ib_en=ib_en;

		ib_st=ia-shift2-threshold; if(ib_st<0) ib_st=0;
		ib_en=ia-shift1+threshold; if(ib_en>nb) ib_en=nb;

		int m=oo; // make it = final_ret
		ib=ib_en;

		int prev_T_ibm1=-1; // T[ca-1][ib+1]

		if(ib_en==nb)
		{
			prev_T_ibm1=T[ib];

			T[ib]=na-ia;
			if(flex_left) T[ib]=0;

			// moved out
			if(T[ib]<m) m=T[ib];

			ib--;
		}
		else
		{
			prev_T_ibm1=T[ib+1];

			int d1=prev_T_ibm1; if(a[ia]!=b[ib]) d1+=sub_cost;
			int d2=T[ib]+1; if(ib<prev_ib_st) d2=oo;
			int ret=d1; if(d2<ret) ret=d2;
			prev_T_ibm1=T[ib]; T[ib]=ret; if(ret<m) m=ret;

			ib--;
		}

		if(ib+1>prev_ib_en)
			prev_T_ibm1=oo;

		for(;ib>ib_st;ib--)
		{
			int d1=prev_T_ibm1; if(a[ia]!=b[ib]) d1+=sub_cost;
			int ret=T[ib+1]; if(T[ib]<ret) ret=T[ib]; ++ret; if(d1<ret) ret=d1;
			prev_T_ibm1=T[ib]; T[ib]=ret; if(ret<m) m=ret;
		}

		ib=ib_st;

		if(ib<ib_en)
		{
			int d1=prev_T_ibm1; if(a[ia]!=b[ib]) d1+=sub_cost;
			int d2=T[ib]+1; if(ib<prev_ib_st) d2=oo;
			int d3=T[ib+1]+1;
			int ret=d1; if(d2<ret) ret=d2; if(d3<ret) ret=d3;
			prev_T_ibm1=T[ib]; T[ib]=ret; if(ret<m) m=ret;
		}

		if(flex_right && ib==0 && T[ib]<=final_ret)
		{
			first_a=ia; first_b=ib;
			final_ret=T[ib];
		}

		if(ia==0)
		{
			if(flex_right)
			{
				for(ib=ib_en;ib>=ib_st;ib--)
				{
					if(T[ib]<=final_ret)
					{
						first_a=ia; first_b=ib;
						final_ret=T[ib];
					}
				}
			}
			else if(!flex_right && ia==0 && 0>=ib_st && 0<=ib_en)
			{
				first_a=ia; first_b=0;
				final_ret=T[0];
			}
		}

		if(m>threshold && final_ret>threshold)
			return -1;
	}

	if(final_ret>threshold) return -1;
	return final_ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
