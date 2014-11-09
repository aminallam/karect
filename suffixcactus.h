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

#include "stack.h"
#include "suffixarray_induced.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ResultLocation
{
	int shift;
	int start_suff_ind;
	int end_suff_ind;
};

struct SuffixFilterInfo
{
	int num_needed_parts;
	int* part_end;

	vector<ResultLocation>* pres;

	int pat_start;
	int start_part;
};

struct EditDistInfo
{
	int max_edit_dist;
	vector<int>* res;
	int max_part_res;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SuffixCactus
{
	const char*	seq;
	int			seq_len;

	int*		suff;
	int*		lcp;
	int*		sibling;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Print(SuffixCactus* suff_cactus)
{
	int i;
	for(i=0;i<suff_cactus->seq_len;i++)
	{
		printf("%d %d %d\n", suff_cactus->suff[i], suff_cactus->lcp[i], suff_cactus->sibling[i]);
	}
	printf("--------------------------------------------------------\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ConstructSibling(int* lcp, int* sibling, int len)
{
	sibling[0]=0;
	int s, skm1=-1, si, sip1;

	for(s=0;s<len-1;s++)
	{
		if(lcp[s]<=lcp[s+1])
		{
			sibling[s+1]=skm1;
			skm1=s;
		}
		else
		{
			sip1=s;
			si=skm1;

			while(lcp[si]>lcp[s+1])
			{
				int r=sibling[si+1];
				sibling[si+1]=sip1;
				sip1=si;
				si=r;
			}

			sibling[s+1]=sip1;
			skm1=si;
		}
	}

	sip1=len-1;
	si=skm1;

	while(si>=0)
	{
		int r=sibling[si+1];
		sibling[si+1]=sip1;
		sip1=si;
		si=r;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Build(SuffixCactus* suff_cactus, const char* _seq, int _seq_len)
{
	suff_cactus->seq=_seq;
	suff_cactus->seq_len=_seq_len;

	suff_cactus->suff=(int*)AllocateHeap(sizeof(int)*(suff_cactus->seq_len));

	int* lcp=(int*)AllocateHeap(sizeof(int)*(suff_cactus->seq_len));

	sais(suff_cactus->seq, suff_cactus->suff, lcp, suff_cactus->seq_len); // (With Lcp)

	suff_cactus->lcp=lcp;
	suff_cactus->sibling=(int*)AllocateHeap(sizeof(int)*(suff_cactus->seq_len));
	ConstructSibling(suff_cactus->lcp, suff_cactus->sibling, suff_cactus->seq_len);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Destroy(SuffixCactus* suff_cactus)
{
	FreeHeap(suff_cactus->suff, sizeof(int)*(suff_cactus->seq_len));
	FreeHeap(suff_cactus->lcp, sizeof(int)*(suff_cactus->seq_len));
	FreeHeap(suff_cactus->sibling, sizeof(int)*(suff_cactus->seq_len));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int first_child(int s, int seq_len, int* sibling)
{
	if(s+1>=seq_len || sibling[s+1]<s+1) return -1;
	return sibling[s+1];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int next_sibling(int s, int* sibling)
{
	if(sibling[s]>=s) return -1;
	return sibling[s];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Return After Last Position of the Matching Interval (not used currently)
int GetLastPosition(SuffixCactus* suff_cactus, int suff_ind, int suff_pos, int allow_equal)
{
	int parent=suff_ind;
	while(parent>=0)
	{
		int cur_suff_ind=first_child(parent, suff_cactus->seq_len, suff_cactus->sibling);
		int v=-1;
		while(cur_suff_ind>=0)
		{
			if(((!allow_equal && suff_cactus->lcp[cur_suff_ind]<=suff_pos) || (allow_equal && suff_cactus->lcp[cur_suff_ind]<suff_pos)) && cur_suff_ind>suff_ind)
				v=cur_suff_ind;
			cur_suff_ind=next_sibling(cur_suff_ind, suff_cactus->sibling);
		}
		if(v>=0) return v;
		// get parent
		int last_sibling=parent;
		cur_suff_ind=parent;
		while(cur_suff_ind>=0)
		{
			cur_suff_ind=next_sibling(cur_suff_ind, suff_cactus->sibling);
			if(cur_suff_ind>=0) last_sibling=cur_suff_ind;
		}
		parent=last_sibling-1;
	}
	return suff_cactus->seq_len;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int Match(SuffixCactus* suff_cactus, const char* pat, int pat_len, vector<int>& res, int get_size_only)
{
	res.clear();

	const char* seq=suff_cactus->seq;
	int seq_len=suff_cactus->seq_len;
	int* suff=suff_cactus->suff;
	int* lcp=suff_cactus->lcp;
	int* sibling=suff_cactus->sibling;

	int cur_suff_ind=0;
	int cur_pat=0;

	while(cur_suff_ind>=0)
	{
		while(cur_pat<pat_len && suff[cur_suff_ind]+cur_pat<seq_len && seq[suff[cur_suff_ind]+cur_pat]==pat[cur_pat]) cur_pat++;
		if(cur_pat==pat_len) break;

		cur_suff_ind=first_child(cur_suff_ind, seq_len, sibling);

		while(cur_suff_ind>=0)
		{
			if(lcp[cur_suff_ind]>=cur_pat)
				break;
			cur_suff_ind=next_sibling(cur_suff_ind, sibling);
		}

		while(cur_suff_ind>=0 && lcp[cur_suff_ind]==cur_pat && pat[cur_pat]!=seq[suff[cur_suff_ind]+cur_pat])
			cur_suff_ind=first_child(cur_suff_ind, seq_len, sibling);
	}

	int num_matches=0;

	if(cur_pat==pat_len)
	{
		if(get_size_only)
		{
			num_matches=GetLastPosition(suff_cactus, cur_suff_ind, pat_len, 1);
			num_matches-=cur_suff_ind;
		}
		else
		{
			num_matches=cur_suff_ind;
			do{
				//printf("Match at = %d\n", suff[cur_suff_ind]);
				//res.add(suff[cur_suff_ind]);
				res.push_back(suff[cur_suff_ind]);
				cur_suff_ind++;
			}
			while(cur_suff_ind<seq_len && lcp[cur_suff_ind]>=pat_len);
			num_matches=cur_suff_ind-num_matches;
		}

		//if(cur_suff_ind!=a) printf("%d %d Error\n", cur_suff_ind, a);
	}

	return num_matches;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct CactusNodeApprox
{
	int					cur_suff_ind;
	int					cur_suff_pos;
	int					cur_pat;
	int					last_sibling;
	int					edit_dist;
	int					cur_filter_part; // used for filter only
	char				last_op;
	char				last_op_char;
	CactusNodeApprox*	prev_node;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MatchApprox(SuffixCactus* suff_cactus, const char* pat, int pat_len, int match_type, EditDistInfo* edit_info, SuffixFilterInfo* filter_info)
{
	if(!edit_info && !filter_info) return;

	if(edit_info) edit_info->res->clear();

	const char* seq=suff_cactus->seq;
	int seq_len=suff_cactus->seq_len;
	int* suff=suff_cactus->suff;
	int* lcp=suff_cactus->lcp;
	int* sibling=suff_cactus->sibling;

	Stack<CactusNodeApprox> stack;
	InitStack(&stack);

	CactusNodeApprox* stack_node=(CactusNodeApprox*)AllocateHeap(sizeof(CactusNodeApprox));

	stack_node->cur_pat=0;
	if(filter_info) stack_node->cur_pat=filter_info->pat_start;

	stack_node->cur_suff_ind=0;
	stack_node->cur_suff_pos=0;
	stack_node->last_sibling=-2;

	stack_node->edit_dist=0;

	stack_node->cur_filter_part=0;
	if(filter_info) stack_node->cur_filter_part=filter_info->start_part;

	stack_node->last_op=OP_NONE;
	stack_node->last_op_char=0;

	Push(&stack, stack_node);

	int cur_suff_ind=-1;
	if(stack_node->last_sibling==-2) cur_suff_ind=first_child(stack_node->cur_suff_ind, seq_len, sibling);
	else cur_suff_ind=stack_node->last_sibling;

	while(cur_suff_ind>=0)
	{
		if(lcp[cur_suff_ind]>=stack_node->cur_suff_pos)
			break;
		cur_suff_ind=next_sibling(cur_suff_ind, sibling);
	}

	stack_node->last_sibling=cur_suff_ind;

	while(cur_suff_ind>=0 && lcp[cur_suff_ind]==stack_node->cur_suff_pos)
	{
		CactusNodeApprox* cur_stack_node=(CactusNodeApprox*)AllocateHeap(sizeof(CactusNodeApprox));
		*cur_stack_node=*stack_node;

		cur_stack_node->cur_suff_ind=cur_suff_ind;
		cur_stack_node->last_sibling=-2;

		Push(&stack, cur_stack_node);

		cur_suff_ind=first_child(cur_suff_ind, seq_len, sibling);
	}

	while(1)
	{
		stack_node=Pop(&stack);
		if(!stack_node) break;

		if(edit_info)
		{
			if(stack_node->cur_pat==pat_len)
			{
				if(stack_node->edit_dist==edit_info->max_edit_dist)
				{
					vector<int>& res=*edit_info->res;

					int num_added=0;

					int cur_suff_ind=stack_node->cur_suff_ind;
					do{
						res.push_back(suff[cur_suff_ind]);

						cur_suff_ind++;

						num_added++;
					}
					while(cur_suff_ind<seq_len && lcp[cur_suff_ind]>stack_node->cur_suff_pos && num_added<edit_info->max_part_res);
				}

				FreeHeap(stack_node, sizeof(CactusNodeApprox));
				continue;
			}
		}
		else
		{
			if(stack_node->cur_pat==filter_info->part_end[stack_node->cur_filter_part])
			{
				if(stack_node->cur_filter_part==filter_info->start_part+filter_info->num_needed_parts-1)
				{
					vector<ResultLocation>& res=*filter_info->pres;

					ResultLocation rloc;
					rloc.shift=-filter_info->pat_start;
					rloc.start_suff_ind=stack_node->cur_suff_ind;
					rloc.end_suff_ind=GetLastPosition(suff_cactus, stack_node->cur_suff_ind, stack_node->cur_suff_pos, 0);

					res.push_back(rloc);

					FreeHeap(stack_node, sizeof(CactusNodeApprox));
					continue;
				}
				else
				{
					stack_node->cur_filter_part++;
				}
			}
		}

		int find_match_char=0;

		{
			if(suff[stack_node->cur_suff_ind]+stack_node->cur_suff_pos<seq_len)
			{
				CactusNodeApprox* new_stack_node=(CactusNodeApprox*)AllocateHeap(sizeof(CactusNodeApprox));
				*new_stack_node=*stack_node;

				if(seq[suff[stack_node->cur_suff_ind]+stack_node->cur_suff_pos]!=pat[stack_node->cur_pat])
				{
					new_stack_node->last_op=OP_SUBSTITUTE;
					new_stack_node->edit_dist++;
				}
				else
				{
					new_stack_node->last_op=OP_MATCH;
					find_match_char=1;
				}

				new_stack_node->last_op_char=seq[suff[stack_node->cur_suff_ind]+stack_node->cur_suff_pos];

				new_stack_node->cur_suff_pos++;
				new_stack_node->cur_pat++;

				if(seq[suff[stack_node->cur_suff_ind]+stack_node->cur_suff_pos]=='$' ||
					((filter_info && new_stack_node->edit_dist>stack_node->cur_filter_part-filter_info->start_part) || (edit_info && new_stack_node->edit_dist>edit_info->max_edit_dist)))
					FreeHeap(new_stack_node, sizeof(CactusNodeApprox));
				else
				{
					if(new_stack_node->last_op==OP_SUBSTITUTE && (stack_node->last_op==OP_INSERT || stack_node->last_op==OP_DELETE)) // || new_stack_node->last_op==OP_MATCH && stack_node->last_op==OP_DELETE && stack_node->last_op_char==new_stack_node->last_op_char)
					{
						FreeHeap(new_stack_node, sizeof(CactusNodeApprox));
					}
					else
					{
						Push(&stack, new_stack_node);

						int cur_suff_ind=-1;
						if(stack_node->last_sibling==-2) cur_suff_ind=first_child(stack_node->cur_suff_ind, seq_len, sibling);
						else cur_suff_ind=stack_node->last_sibling;

						while(cur_suff_ind>=0)
						{
							if(lcp[cur_suff_ind]>=new_stack_node->cur_suff_pos)
								break;
							cur_suff_ind=next_sibling(cur_suff_ind, sibling);
						}

						stack_node->last_sibling=cur_suff_ind;

						while(cur_suff_ind>=0 && lcp[cur_suff_ind]==new_stack_node->cur_suff_pos)
						{
							CactusNodeApprox* cur_stack_node=(CactusNodeApprox*)AllocateHeap(sizeof(CactusNodeApprox));
							*cur_stack_node=*new_stack_node;

							cur_stack_node->cur_suff_ind=cur_suff_ind;
							cur_stack_node->last_sibling=-2;

							Push(&stack, cur_stack_node);

							cur_suff_ind=first_child(cur_suff_ind, seq_len, sibling);
						}
					}
				}
			}
		}

		if(match_type==MATCH_SUBSTITUTE_ONLY)
		{
			FreeHeap(stack_node, sizeof(CactusNodeApprox));
			continue;
		}

		// Delete Case

		if(stack_node->last_op!=OP_INSERT)
		{
			//if(stack_node->cur_pat<pat_len && (max_errors && stack_node->edit_dist<max_errors[stack_node->cur_pat] || !max_errors && stack_node->edit_dist<max_edit_dist))
			if(stack_node->cur_pat<pat_len && ((filter_info && stack_node->edit_dist<stack_node->cur_filter_part-filter_info->start_part) || (edit_info && stack_node->edit_dist<edit_info->max_edit_dist)))
			{
				CactusNodeApprox* new_stack_node=(CactusNodeApprox*)AllocateHeap(sizeof(CactusNodeApprox));
				*new_stack_node=*stack_node;

				new_stack_node->cur_pat++;
				new_stack_node->edit_dist++;

				new_stack_node->last_op=OP_DELETE;
				new_stack_node->last_op_char=pat[stack_node->cur_pat];

				Push(&stack, new_stack_node);
			}
		}

		// dont insert a character in your pattern, while the next pattern char is same
		if(find_match_char)
		{
			FreeHeap(stack_node, sizeof(CactusNodeApprox));
			continue;
		}

		// Insert Case

		if(stack_node->last_op!=OP_DELETE)
		{
			if((filter_info && stack_node->edit_dist<stack_node->cur_filter_part-filter_info->start_part) || (edit_info && stack_node->edit_dist<edit_info->max_edit_dist))
			{
				if(suff[stack_node->cur_suff_ind]+stack_node->cur_suff_pos<seq_len && seq[suff[stack_node->cur_suff_ind]+stack_node->cur_suff_pos]!='$')
				{
					CactusNodeApprox* new_stack_node=(CactusNodeApprox*)AllocateHeap(sizeof(CactusNodeApprox));
					*new_stack_node=*stack_node;

					new_stack_node->cur_suff_pos++;
					new_stack_node->edit_dist++;

					new_stack_node->last_op=OP_INSERT;
					new_stack_node->last_op_char=seq[suff[stack_node->cur_suff_ind]+stack_node->cur_suff_pos];

					Push(&stack, new_stack_node);

					int cur_suff_ind=-1;
					if(stack_node->last_sibling==-2) cur_suff_ind=first_child(stack_node->cur_suff_ind, seq_len, sibling);
					else cur_suff_ind=stack_node->last_sibling;

					while(cur_suff_ind>=0)
					{
						if(lcp[cur_suff_ind]>=new_stack_node->cur_suff_pos)
							break;
						cur_suff_ind=next_sibling(cur_suff_ind, sibling);
					}

					stack_node->last_sibling=cur_suff_ind;

					while(cur_suff_ind>=0 && lcp[cur_suff_ind]==new_stack_node->cur_suff_pos)
					{
						CactusNodeApprox* cur_stack_node=(CactusNodeApprox*)AllocateHeap(sizeof(CactusNodeApprox));
						*cur_stack_node=*new_stack_node;

						cur_stack_node->cur_suff_ind=cur_suff_ind;
						cur_stack_node->last_sibling=-2;

						Push(&stack, cur_stack_node);

						cur_suff_ind=first_child(cur_suff_ind, seq_len, sibling);
					}
				}
			}
		}

		FreeHeap(stack_node, sizeof(CactusNodeApprox));
	}

	DestroyStack(&stack);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
