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
// used in suffix tree and suffix cactus
enum{OP_NONE, OP_MATCH, OP_SUBSTITUTE, OP_DELETE, OP_INSERT};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum{SMALL_FIRST, LARGE_FIRST, SMALL_ENDS, LARGE_ENDS};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ConstructPartSizes(int pat_len, int num_parts, int partitioning_type, int* part_size, int* part_end, int* part_end_rev)
{
	int i;

	if(partitioning_type==SMALL_FIRST)
	{
		for(i=0;i<num_parts;i++) part_size[i]=(pat_len+i)/num_parts; // small parts first
	}
	else if(partitioning_type==LARGE_FIRST)
	{
		for(i=0;i<num_parts;i++) part_size[i]=(pat_len+num_parts-1-i)/num_parts; // large parts first
	}
	else if(partitioning_type==SMALL_ENDS)
	{
		int max_part_size=(pat_len+num_parts-1)/num_parts;
		int min_part_size=max_part_size; if(pat_len/max_part_size<num_parts) min_part_size--;

		// assign each part to min_size
		for(i=0;i<num_parts;i++) part_size[i]=min_part_size;

		int rem=pat_len-min_part_size*num_parts;

		//if(min_part_size<max_part_size)
		if(rem)
		{
			// start from mid and go zigzag assign remaining

			int mid=num_parts/2;
			part_size[mid]++; rem--;

			int sh=1;
			while(rem)
			{
				part_size[mid-sh]++; rem--;
				if(rem) {part_size[mid+sh]++; rem--;}
				sh++;
			}
		}
	}
	else if(partitioning_type==LARGE_ENDS)
	{
		int max_part_size=(pat_len+num_parts-1)/num_parts;
		int min_part_size=max_part_size; if(pat_len/max_part_size<num_parts) min_part_size--;

		// assign each part to min_size
		for(i=0;i<num_parts;i++) part_size[i]=max_part_size;

		int rem=max_part_size*num_parts-pat_len;

		if(rem)
		{
			// start from mid and go zigzag assign remaining

			int mid=num_parts/2;
			part_size[mid]--; rem--;

			int sh=1;
			while(rem)
			{
				part_size[mid-sh]--; rem--;
				if(rem) {part_size[mid+sh]--; rem--;}
				sh++;
			}
		}
	}

	for(i=0;i<num_parts;i++) part_end[i]=part_size[i]+(i?part_end[i-1]:0);
	for(i=0;i<num_parts;i++) part_end_rev[i]=part_size[num_parts-1-i]+(i?part_end_rev[i-1]:0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct CandidateShift
{
	long long start_read;
	long long end_read;

	int shift1;
	int shift2;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct PossibleMatch
{
	char* str;
	char* qual;
	int str_len; // if =0 means not computed yet (but this will not happen for now)
	int shift1, shift2;
	int max_overlap; // not computed here

	bool operator < (const PossibleMatch& pm) const
	{
		return max_overlap>pm.max_overlap;
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ComputeMaxOverlap(int na, int nb, int shift1, int shift2)
{
	int cur_max_overlap_size=na; if(nb<na) cur_max_overlap_size=nb;

	if(shift1>0)
	{
		cur_max_overlap_size=na-shift1;
		if(nb<cur_max_overlap_size) cur_max_overlap_size=nb;
	}
	else if(shift2<0)
	{
		cur_max_overlap_size=nb+shift2;
		if(na<cur_max_overlap_size) cur_max_overlap_size=na;
	}

	return cur_max_overlap_size;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int					kmer_entry_size;

unsigned int		set1_mask;
unsigned int		get1_mask;

int					set2_shift;
int					set2_num_lowest_bits;
unsigned int		set2_mask;

int					set3_shift;
int					set3_num_lowest_bits;
unsigned long long	set3_mask;

template <class TypeKmer, class TypeIndex>
void InitializeKmerEntry(int num_bits_kmer, int num_bits_index)
{
	num_bits_kmer=((num_bits_kmer+3)/4)*4;
	num_bits_index=((num_bits_index+7)/8)*8;
	kmer_entry_size=(2*num_bits_kmer+num_bits_index)/8;

	int num_highest_bits=sizeof(TypeKmer)*8-num_bits_kmer;
	set1_mask=((1<<num_highest_bits)-1)<<num_bits_kmer;

	get1_mask=(1<<num_bits_kmer)-1;

	set2_shift=(num_bits_kmer+num_bits_kmer)/8-sizeof(TypeKmer);
	set2_num_lowest_bits=sizeof(TypeKmer)*8-num_bits_kmer;
	set2_mask=(1<<set2_num_lowest_bits)-1;

	set3_shift=(2*num_bits_kmer+num_bits_index)/8-sizeof(TypeIndex);
	set3_num_lowest_bits=sizeof(TypeIndex)*8-num_bits_index;
	set3_mask=(((long long)1)<<set3_num_lowest_bits)-1;
}

inline void SetRight(char* entry, unsigned int a) // num_bits should be multiple of 4
{
	unsigned int* poa=(unsigned int*)entry;
	unsigned int highest_bits=(*poa)&set1_mask;
	*poa=a|highest_bits;
}

inline unsigned int GetRight(char* entry)
{
	unsigned int poa=*(unsigned int*)entry;
	return poa&get1_mask;
}

inline void SetLeft(char* entry, unsigned int a) // num_bits is the same for first element and second element
{
	unsigned int* poa=(unsigned int*)(entry+set2_shift);
	unsigned int lowest_bits=(*poa)&set2_mask;
	*poa=(a<<set2_num_lowest_bits)|lowest_bits;
}

inline unsigned int GetLeft(char* entry)
{
	unsigned int poa=*(unsigned int*)(entry+set2_shift);
	return poa>>set2_num_lowest_bits;
}

inline void SetIndex(char* entry, unsigned long long a)
{
	unsigned long long* poa=(unsigned long long*)(entry+set3_shift); // it starts actually at data+3, we will get this and ignore the lowest byte
	unsigned long long lowest_bits=(*poa)&set3_mask;
	*poa=(a<<set3_num_lowest_bits)|lowest_bits;
}

inline unsigned long long GetIndex(char* entry)
{
	unsigned long long poa=*(unsigned long long*)(entry+set3_shift);
	return poa>>set3_num_lowest_bits;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SuffixFilter
{
	char*			all_reads;
	char*			all_quals;
	long long*		all_reads_starts;

	long long*		kmer_sizes; // Make it long long
	long long*		kmer_starts; // Make it long long
	char*			kmer_adj;

	long long		total_reads_size;
	long long		num_reads;
	long long		num_kmers;
	int				max_read_len;

	int				min_overlap;
	double			min_overlap_percentage;
	double			max_overlap_error_rate;
	double			truncate_factor;
	double			min_straight_edge_val;
	double			min_straight_edge_percentage;
	double			min_read_weight;
	int				max_mul_len_matches;
	int				max_kmer_slots;

	int				max_part_res;
	int				match_type;

	int				kmer_size;
	int				max_kmer_errors;
	int				kmer_with_errors_size;
	int				read_index_bits;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BuildSmallSubTable(SuffixFilter* suff_filter, unsigned int cur_start_kmer, unsigned int cur_end_kmer, long long cur_start_read, long long cur_end_read)
{
	int kmer_size=suff_filter->kmer_size;
	int kmer_with_errors_size=suff_filter->kmer_with_errors_size;

	long long i;
	int j;

	for(i=cur_start_read;i<cur_end_read;i++)
	{
		long long st_read=suff_filter->all_reads_starts[i];
		long long en_read=suff_filter->all_reads_starts[i+1]-1;
		char* read=suff_filter->all_reads+st_read;
		int len=en_read-st_read;
		if(len<kmer_size) continue;

		unsigned int kmer=GetKmer(read, kmer_size, kmer_size);

		unsigned int right_kmer=GetKmer(read+kmer_size, len-kmer_size, kmer_with_errors_size);
		unsigned int left_kmer=0;//GetRevKmer(read-1, 0, KMER_WITH_ERRORS_SIZE);

		for(j=0;j<len-kmer_size+1;j++)
		{
			if(j)
			{
				int new_char=CharToInt[(int)read[j+kmer_size-1]];
				kmer=(kmer>>BITS_PER_SYMBOL)|(new_char<<((kmer_size-1)*BITS_PER_SYMBOL));

				new_char=0; if(j+kmer_size+kmer_with_errors_size-1<len) new_char=CharToInt[(int)read[j+kmer_size+kmer_with_errors_size-1]];
				right_kmer=(right_kmer>>BITS_PER_SYMBOL)|(new_char<<((kmer_with_errors_size-1)*BITS_PER_SYMBOL));

				new_char=CharToInt[(int)read[j-1]];
				left_kmer=((left_kmer<<BITS_PER_SYMBOL)&((1<<(kmer_with_errors_size*BITS_PER_SYMBOL))-1))|new_char;
			}

			if(kmer<cur_start_kmer || kmer>cur_end_kmer || suff_filter->kmer_sizes[kmer]>=suff_filter->max_kmer_slots) continue;

			long long kmer_ind=suff_filter->kmer_starts[kmer]+suff_filter->kmer_sizes[kmer];

			char* entry=suff_filter->kmer_adj+kmer_ind*kmer_entry_size;
			memset(entry, 0, kmer_entry_size);
			SetRight(entry, right_kmer);
			SetLeft(entry, left_kmer);
			SetIndex(entry, st_read+j);

			suff_filter->kmer_sizes[kmer]++;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SolutionBuild
{
	SuffixFilter*	suff_filter;
	unsigned int	start_kmer;
	unsigned int	end_kmer;
	long long 		start_read;
	long long 		end_read;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void* threadfunctionBuild(void* arg)
{
	SolutionBuild* sol=(SolutionBuild*)arg;
	BuildSmallSubTable(sol->suff_filter, sol->start_kmer, sol->end_kmer, sol->start_read, sol->end_read);
	return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BuildSubTable(SuffixFilter* suff_filter, int num_threads, long long cache_block_size,
							long long start_kmer, long long end_kmer, long long num_entries)
{
	long long i;
	long long k;
	int j;

	long long max_num_entries=(num_entries+num_threads-1)/num_threads;

	SolutionBuild** psols=new SolutionBuild*[num_threads];
	pthread_t* thread_ID=new pthread_t[num_threads];
	void** exit_status=new void*[num_threads];
	for(j=0;j<num_threads;j++) psols[j]=new SolutionBuild;

	long long cur_start_read=0;

	while(cur_start_read<suff_filter->num_reads)
	{
		long long cur_read_sizes=0;

		for(i=cur_start_read;cur_read_sizes<=cache_block_size && i<suff_filter->num_reads;i++)
		{
			cur_read_sizes=suff_filter->all_reads_starts[i+1]-suff_filter->all_reads_starts[cur_start_read];
		}

		long long cur_end_read=i;

		int cur_thread=0;
		long long cur_start_kmer=start_kmer;

		while(cur_start_kmer<end_kmer)
		{
			long long cur_num_entries=0;

			for(k=cur_start_kmer;cur_num_entries<=max_num_entries && k<end_kmer;k++)
			{
				cur_num_entries=((k<end_kmer-1)?suff_filter->kmer_starts[k+1]:num_entries)-suff_filter->kmer_starts[cur_start_kmer];
			}
			long long cur_end_kmer=k; // kmer after last

			if(cur_thread==num_threads-1) cur_end_kmer=end_kmer;

			psols[cur_thread]->suff_filter=suff_filter;
			psols[cur_thread]->start_kmer=cur_start_kmer;
			psols[cur_thread]->end_kmer=cur_end_kmer-1;
			psols[cur_thread]->start_read=cur_start_read;
			psols[cur_thread]->end_read=cur_end_read;

			cur_thread++;

			cur_start_kmer=cur_end_kmer;
		}

		for(j=1;j<cur_thread;j++)
		{
			pthread_create(&thread_ID[j], NULL, threadfunctionBuild, psols[j]);
		}

		threadfunctionBuild(psols[0]);

		for(j=1;j<cur_thread;j++)
		{
			pthread_join(thread_ID[j], &exit_status[j]);
		}

		cur_start_read=cur_end_read;
	}

	for(j=0;j<num_threads;j++) delete psols[j];
	delete[] psols;
	delete[] thread_ID;
	delete[] exit_status;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BuildAllSubTables(SuffixFilter* suff_filter, int num_threads, long long cache_block_size, long long total_num_entries)
{
	int kmer_size=suff_filter->kmer_size;
	long long num_kmers=((long long)1)<<(kmer_size*BITS_PER_SYMBOL);

	// Added after using kmer_starts to get sizes
	memset(suff_filter->kmer_sizes, 0, sizeof(long long)*num_kmers);

	long long max_num_entries=total_num_entries;

	long long k;
	long long cur_start_kmer=0;

	while(cur_start_kmer<num_kmers)
	{
		long long cur_num_entries=0;

		for(k=cur_start_kmer;cur_num_entries<=max_num_entries && k<num_kmers;k++)
		{
			cur_num_entries=((k<num_kmers-1)?suff_filter->kmer_starts[k+1]:total_num_entries)-suff_filter->kmer_starts[cur_start_kmer];
		}
		long long cur_end_kmer=k; // kmer after last

		if(cur_num_entries>0) BuildSubTable(suff_filter, num_threads, cache_block_size, cur_start_kmer, cur_end_kmer, cur_num_entries);

		cur_start_kmer=cur_end_kmer;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BuildTable(SuffixFilter* suff_filter, int num_threads, long long cache_block_size)
{
	long long i;
	int j;

	int kmer_size=suff_filter->kmer_size;

	long long num_kmers=((long long)1)<<(kmer_size*BITS_PER_SYMBOL);

	suff_filter->kmer_sizes=(long long*)AllocateHeap(sizeof(long long)*num_kmers);
	suff_filter->kmer_starts=(long long*)AllocateHeap(sizeof(long long)*num_kmers);

	memset(suff_filter->kmer_sizes, 0, sizeof(long long)*num_kmers);

	for(i=0;i<suff_filter->num_reads;i++)
	{
		long long st_read=suff_filter->all_reads_starts[i];
		long long en_read=suff_filter->all_reads_starts[i+1]-1;
		char* read=suff_filter->all_reads+st_read;
		int len=en_read-st_read;
		if(len<kmer_size) continue;

		unsigned int kmer=GetKmer(read, kmer_size, kmer_size);

		for(j=0;j<len-kmer_size+1;j++)
		{
			if(j)
			{
				int new_char=CharToInt[(int)read[j+kmer_size-1]];
				kmer=(kmer>>BITS_PER_SYMBOL)|(new_char<<((kmer_size-1)*BITS_PER_SYMBOL));
			}

			if(suff_filter->kmer_sizes[kmer]<suff_filter->max_kmer_slots) suff_filter->kmer_sizes[kmer]++;
		}
	}

	suff_filter->kmer_starts[0]=0;

	unsigned int k;

	for(k=1;k<num_kmers;k++)
	{
		suff_filter->kmer_starts[k]=suff_filter->kmer_starts[k-1]+suff_filter->kmer_sizes[k-1];
	}

	long long total_num_entries=suff_filter->kmer_starts[k-1]+suff_filter->kmer_sizes[k-1];
	if(total_num_entries==0) return;

	suff_filter->kmer_adj=(char*)AllocateHeap(total_num_entries*kmer_entry_size);

	BuildAllSubTables(suff_filter, num_threads, cache_block_size, total_num_entries);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Result
{
	long long	ref_pos;
	int			cur_pat_start;

	Result() {}
	Result(long long a, int b) {ref_pos=a; cur_pat_start=b;}

	bool operator<(const Result& r)const
	{
		return ref_pos<r.ref_pos;
	}
};

int MatchExact2(SuffixFilter* suff_filter, const char* pat, vector<Result>& res2, int get_size_only)
{
	int kmer_size=suff_filter->kmer_size;

	long long i;
	unsigned int kmer=GetKmer(pat, kmer_size, kmer_size);
	unsigned int right_kmer=GetKmer(pat+kmer_size, kmer_size, kmer_size);
	unsigned int mask=(1<<(kmer_size*BITS_PER_SYMBOL))-1;

	long long st_kmer=suff_filter->kmer_starts[kmer];
	long long en_kmer=st_kmer+suff_filter->kmer_sizes[kmer];

	int num_added=0;

	for(i=st_kmer;i<en_kmer && num_added<suff_filter->max_part_res;i++)
	{
		char* entry=suff_filter->kmer_adj+i*kmer_entry_size;

		if(right_kmer==(GetRight(entry)&mask))
		{
			if(!get_size_only) res2.push_back(Result(GetIndex(entry),0));
			num_added++;
		}
	}

	return num_added;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MatchExact(SuffixFilter* suff_filter, const char* pat, vector<Result>& res2, vector<Result>& res3, int get_size_only)
{
	int kmer_size=suff_filter->kmer_size;

	long long i;
	unsigned int left_kmer=GetRevKmer(pat+kmer_size-1, kmer_size, kmer_size);
	unsigned int kmer=GetKmer(pat+kmer_size, kmer_size, kmer_size);
	unsigned int right_kmer=GetKmer(pat+2*kmer_size, kmer_size, kmer_size);

	unsigned int mask=(1<<(kmer_size*BITS_PER_SYMBOL))-1;

	long long st_kmer=suff_filter->kmer_starts[kmer];
	long long en_kmer=st_kmer+suff_filter->kmer_sizes[kmer];

	int num_added_2[4]={0};
	int num_added_3=0;

	for(i=st_kmer;i<en_kmer && num_added_3<suff_filter->max_part_res;i++)
	{
		char* entry=suff_filter->kmer_adj+i*kmer_entry_size;

		int n=0;
		if(right_kmer==(GetRight(entry)&mask)) n|=1;
		if(left_kmer==(GetLeft(entry)&mask)) n|=2;

		if(n==3 && num_added_3<suff_filter->max_part_res)
		{
			if(!get_size_only) res3.push_back(Result(GetIndex(entry)-kmer_size, 0));
			num_added_3++;
		}
		if(n>=1 && num_added_2[n]<suff_filter->max_part_res)
		{
			if(!get_size_only)
			{
				if(n&1) res2.push_back(Result(GetIndex(entry), kmer_size));
				if(n&2) res2.push_back(Result(GetIndex(entry)-kmer_size, 0)); // put else to avoid repetition (the two cases are very similar, if not the same)
			}
			num_added_2[n]++;
		}
	}
	int num_added=num_added_2[0]+num_added_2[1]+num_added_2[2]+num_added_2[3];
	return num_added;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetCandidates2(SuffixFilter* suff_filter, unsigned short* kmer_flags, const char* pat, vector<Result>& res2, const int max_num_errors, int dir, int shift, bool allow_exact)
{
	long long i;
	int kmer_size=suff_filter->kmer_size;
	int kmer_with_errors_size=suff_filter->kmer_with_errors_size;

	int* T=new int[2*(kmer_size+1)];
	int k;

	unsigned int mask=(1<<(kmer_size*BITS_PER_SYMBOL))-1;

	vector<unsigned int> added_kmers;
	bool exact_done=false;

	if(dir&1)
	{
		exact_done=true;
		unsigned int kmer_a=GetKmer(pat, kmer_size, kmer_size);
		long long st_kmer=suff_filter->kmer_starts[kmer_a];
		long long en_kmer=st_kmer+suff_filter->kmer_sizes[kmer_a];

		for(i=st_kmer;i<en_kmer;i++)
		{
			char* entry=suff_filter->kmer_adj+i*kmer_entry_size;

			unsigned int right_kmer=GetRight(entry);
			unsigned char* pnum_added=(unsigned char*)&kmer_flags[right_kmer&mask];

			if(*pnum_added)
			{
				if((*pnum_added)>1 && (*pnum_added)<=suff_filter->max_part_res)
				{
					res2.push_back(Result(GetIndex(entry), shift));
					(*pnum_added)++;
				}
				continue;
			}

			char* right_kmer_str=new char[kmer_with_errors_size];
			GetStrFromKmer(right_kmer, right_kmer_str, kmer_with_errors_size);

			int ed=ApproxEditDistanceFlexibleEndA(T, right_kmer_str, kmer_with_errors_size, pat+kmer_size, kmer_size, max_num_errors, suff_filter->match_type, 0);

			if(ed>0 || (ed==0 && allow_exact))
			{
				res2.push_back(Result(GetIndex(entry), shift));
				(*pnum_added)++;
			}

			(*pnum_added)++;
			added_kmers.push_back(right_kmer&mask);

			delete[] right_kmer_str;
		}

		for(k=0;k<(int)added_kmers.size();k++) kmer_flags[added_kmers[k]]=0;
		added_kmers.clear();
	}

	if(dir&2)
	{
		unsigned int kmer_b=GetKmer(pat+kmer_size, kmer_size, kmer_size);
		long long st_kmer=suff_filter->kmer_starts[kmer_b];
		long long en_kmer=st_kmer+suff_filter->kmer_sizes[kmer_b];

		for(i=st_kmer;i<en_kmer;i++)
		{
			char* entry=suff_filter->kmer_adj+i*kmer_entry_size;

			unsigned int left_kmer=GetLeft(entry);

			unsigned char* pnum_added=(unsigned char*)&kmer_flags[left_kmer&mask];

			if(*pnum_added)
			{
				if((*pnum_added)>1 && (*pnum_added)<=suff_filter->max_part_res)
				{
					res2.push_back(Result(GetIndex(entry)-kmer_size, shift));
					(*pnum_added)++;
				}
				continue;
			}

			int j;
			char* left_kmer_str=new char[kmer_with_errors_size];
			GetStrFromKmer(left_kmer, left_kmer_str, kmer_with_errors_size);

			char* pat_rev=new char[kmer_size];
			for(j=0;j<kmer_size;j++) pat_rev[j]=pat[kmer_size-1-j];

			int ed=ApproxEditDistanceFlexibleEndA(T, left_kmer_str, kmer_with_errors_size, pat_rev, kmer_size, max_num_errors, suff_filter->match_type, 0);

			if(ed>0 || (ed==0 && !exact_done && allow_exact)) // the ed==0 is avoided here exists because it means that it is added in the previous part
			{
				res2.push_back(Result(GetIndex(entry)-kmer_size, shift));
				(*pnum_added)++;
			}

			(*pnum_added)++;
			added_kmers.push_back(left_kmer&mask);

			delete[] pat_rev;
			delete[] left_kmer_str;
		}

		for(k=0;k<(int)added_kmers.size();k++) kmer_flags[added_kmers[k]]=0;
		added_kmers.clear();
	}

	delete[] T;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetCandidates(SuffixFilter* suff_filter, unsigned short* kmer_flags, const char* pat, vector<Result>& res2, vector<Result>& res3, const int max_num_errors)
{
	long long i;
	int kmer_size=suff_filter->kmer_size;
	int kmer_with_errors_size=suff_filter->kmer_with_errors_size;

	int* T=new int[2*(kmer_size+1)];

	unsigned int pat_left_kmer=GetRevKmer(pat+kmer_size-1, kmer_size, kmer_size);
	unsigned int pat_mid_kmer=GetKmer(pat+kmer_size, kmer_size, kmer_size);
	unsigned int pat_right_kmer=GetKmer(pat+2*kmer_size, kmer_size, kmer_size);

	long long st_kmer=suff_filter->kmer_starts[pat_mid_kmer];
	long long en_kmer=st_kmer+suff_filter->kmer_sizes[pat_mid_kmer];

	unsigned int mask=(1<<(kmer_size*BITS_PER_SYMBOL))-1;

	vector<unsigned int> added_kmers;

	char* right_kmer_str=new char[kmer_with_errors_size];

	for(i=st_kmer;i<en_kmer;i++)
	{
		char* entry=suff_filter->kmer_adj+i*kmer_entry_size;

		int allow_res3=(pat_left_kmer==(GetLeft(entry)&mask));

		unsigned int right_kmer=GetRight(entry);

		unsigned char* num_added=(unsigned char*)&kmer_flags[right_kmer&mask];

		if(num_added[0])
		{
			if(num_added[0]>1)
			{
				if(num_added[0]<=suff_filter->max_part_res || (allow_res3 && num_added[1]<suff_filter->max_part_res))
				{
					res2.push_back(Result(GetIndex(entry), kmer_size));
					num_added[0]++;
				}

				if(allow_res3 && num_added[1]<suff_filter->max_part_res)
				{
					res3.push_back(Result(GetIndex(entry)-kmer_size, 0));
					num_added[1]++;
				}
			}

			continue;
		}

		GetStrFromKmer(right_kmer, right_kmer_str, kmer_with_errors_size);

		int ed=ApproxEditDistanceFlexibleEndA(T, right_kmer_str, kmer_with_errors_size, pat+2*kmer_size, kmer_size, max_num_errors, suff_filter->match_type, 0);

		if(ed>=0)
		{
			res2.push_back(Result(GetIndex(entry), kmer_size));
			num_added[0]++;

			if(allow_res3)
			{
				res3.push_back(Result(GetIndex(entry)-kmer_size, 0));
				num_added[1]++;
			}
		}

		num_added[0]++;
		added_kmers.push_back(right_kmer&mask);
	}

	delete[] right_kmer_str;

	int k;
	for(k=0;k<(int)added_kmers.size();k++) kmer_flags[added_kmers[k]]=0;
	added_kmers.clear();

	char* left_kmer_str=new char[kmer_with_errors_size];

	char* pat_rev=new char[kmer_size];
	for(i=0;i<kmer_size;i++) pat_rev[i]=pat[kmer_size-1-i];

	for(i=st_kmer;i<en_kmer;i++)
	{
		char* entry=suff_filter->kmer_adj+i*kmer_entry_size;

		int allow_res3=(pat_right_kmer==(GetRight(entry)&mask));

		unsigned int left_kmer=GetLeft(entry);

		unsigned char* num_added=(unsigned char*)&kmer_flags[left_kmer&mask];

		if(num_added[0])
		{
			if(num_added[0]>1)
			{
				if(num_added[0]<=suff_filter->max_part_res || (allow_res3 && num_added[1]<suff_filter->max_part_res))
				{
					res2.push_back(Result(GetIndex(entry)-kmer_size, 0));
					num_added[0]++;
				}

				if(allow_res3 && num_added[1]<suff_filter->max_part_res)
				{
					res3.push_back(Result(GetIndex(entry)-kmer_size, 0));
					num_added[1]++;
				}
			}

			continue;
		}

		GetStrFromKmer(left_kmer, left_kmer_str, kmer_with_errors_size);

		int ed=ApproxEditDistanceFlexibleEndA(T, left_kmer_str, kmer_with_errors_size, pat_rev, kmer_size, max_num_errors, suff_filter->match_type, 0);

		if(ed>=0)
		{
			res2.push_back(Result(GetIndex(entry)-kmer_size, 0));
			num_added[0]++;

			if(allow_res3 && ed>0)
			{
				res3.push_back(Result(GetIndex(entry)-kmer_size, 0));
				num_added[1]++;
			}
		}

		num_added[0]++;
		added_kmers.push_back(left_kmer&mask);
	}

	for(k=0;k<(int)added_kmers.size();k++) kmer_flags[added_kmers[k]]=0;
	added_kmers.clear();

	delete[] left_kmer_str;

	delete[] pat_rev;
	delete[] T;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int GetOverheadPerKmer()
{
	return kmer_entry_size;
}

int GetOverheadPerChar(int is_fastq)
{
	return 1+is_fastq;
}

long long GetFixedOverhead(int kmer_size, int num_threads, long long global_max_read_len)
{
	long long num_kmers=((long long)1)<<(kmer_size*BITS_PER_SYMBOL);

	long long overhead=num_kmers*(2*sizeof(long long));

	overhead+=num_threads*num_kmers*sizeof(unsigned short); // kmer_flags
	overhead+=num_threads*3*(global_max_read_len+1)*sizeof(int); // ErrorGraph

	return overhead;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

long long Build(SuffixFilter* suff_filter,
			const char* file_name, int file_block_size, long long cache_block_size, long long start_file, long long max_num_kmers, int insert_complement, //vector<Read>& reads, int total_reads_size,
			int min_overlap, double min_overlap_percentage, double max_overlap_error_rate, double truncate_factor, double min_straight_edge_val, double min_straight_edge_percentage, double min_read_weight, int max_mul_len_matches, int max_kmer_slots,
			int max_part_res, int kmer_size, int max_kmer_errors, int kmer_with_errors_size, int read_index_bits, int match_type, int num_threads, bool is_fastq, int use_qual)
{
	suff_filter->min_overlap=min_overlap;
	suff_filter->min_overlap_percentage=min_overlap_percentage;
	suff_filter->max_overlap_error_rate=max_overlap_error_rate;
	suff_filter->truncate_factor=truncate_factor;
	suff_filter->min_straight_edge_val=min_straight_edge_val;
	suff_filter->min_straight_edge_percentage=min_straight_edge_percentage;
	suff_filter->min_read_weight=min_read_weight;
	suff_filter->max_mul_len_matches=max_mul_len_matches;
	suff_filter->max_kmer_slots=max_kmer_slots;
	suff_filter->match_type=match_type;
	suff_filter->max_part_res=max_part_res;
	suff_filter->kmer_size=kmer_size;
	suff_filter->max_kmer_errors=max_kmer_errors;
	suff_filter->kmer_with_errors_size=kmer_with_errors_size;
	suff_filter->read_index_bits=read_index_bits;

	long long prev_start_file=start_file;
	GetReads(file_name, file_block_size, start_file, insert_complement, kmer_size, max_num_kmers, &suff_filter->all_reads, &suff_filter->all_quals, &suff_filter->all_reads_starts, suff_filter->num_reads, suff_filter->num_kmers, suff_filter->total_reads_size, suff_filter->max_read_len, is_fastq, use_qual);

	printf("Total Reads Size = %lld\n", suff_filter->total_reads_size); fflush(NULL);

	printf("Building Index [%lld - %lld]\n", prev_start_file, start_file); fflush(NULL);

	BuildTable(suff_filter, num_threads, cache_block_size);

	printf("Building Index Finished\n"); fflush(NULL);

	return start_file;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Destroy(SuffixFilter* suff_filter)
{
	long long num_kmers=((long long)1)<<(suff_filter->kmer_size*BITS_PER_SYMBOL);

	FreeHeap(suff_filter->all_reads, suff_filter->total_reads_size+2);
	if(suff_filter->all_quals) FreeHeap(suff_filter->all_quals, suff_filter->total_reads_size+2);
	FreeHeap(suff_filter->all_reads_starts, sizeof(long long)*(suff_filter->num_reads+1));

	FreeHeap(suff_filter->kmer_adj, kmer_entry_size*(suff_filter->kmer_starts[num_kmers-1]+suff_filter->kmer_sizes[num_kmers-1]));

	FreeHeap(suff_filter->kmer_sizes, sizeof(long long)*num_kmers);
	FreeHeap(suff_filter->kmer_starts, sizeof(long long)*num_kmers);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PrintPatternPart(int i, int num_group_parts, int num_parts, const char* pat, int pat_len, int* part_end)
{
	int part_start=0; if(i) part_start=part_end[i-1];
	int j, end=pat_len; if(i+num_group_parts-1<num_parts) end=part_end[i+num_group_parts-1];
	char* part_str=new char[MAX_READ_LEN];
	for(j=part_start;j<end;j++) part_str[j-part_start]=pat[j];
	part_str[end-part_start]=0;
	printf("%d) Part [%s]: ", i, part_str); //filtered_part_res[i].sort(); filtered_part_res[i].print();
	printf("\n");
	delete[] part_str;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<Result>& GetCandidatesAll(SuffixFilter* suff_filter, unsigned short* kmer_flags, const char* pat,
						vector<Result>& res2, vector<Result>& res3, int mode, const int max_num_errors)
{
	// the length of pat is 3*kmer_size

	res2.clear();
	res3.clear();

	int kmer_size=suff_filter->kmer_size;

	if(mode==2)
	{
		if(max_num_errors==0)
		{
			MatchExact2(suff_filter, pat, res2, 0);
		}
		else
		{
			GetCandidates2(suff_filter, kmer_flags, pat, res2, max_num_errors, 3, 0, true);
		}
	}
	else
	{
		if(max_num_errors==0)
		{
			MatchExact(suff_filter, pat, res2, res3, 0);
			if((int)res3.size()>=suff_filter->max_part_res) return res3;
		}
		else
		{
			GetCandidates(suff_filter, kmer_flags, pat, res2, res3, max_num_errors);
			if((int)res3.size()>=suff_filter->max_part_res) return res3;

			GetCandidates2(suff_filter, kmer_flags, pat, res2, max_num_errors, 1, 0, false);
			GetCandidates2(suff_filter, kmer_flags, pat+kmer_size, res2, max_num_errors, 2, kmer_size, false);
		}
	}

	return res2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ComputeFilter(SuffixFilter* suff_filter, unsigned short* kmer_flags,
					const char* pat, int pat_len, vector<PossibleMatch>& possible_matches)//BTreeInMem<BTreeKeyInt<long long>, BTreeValInt<CandidateShift> >& cand_shifts)
{
	int kmer_size=suff_filter->kmer_size;

	int mode=3;
	int target_kmer_size=3*kmer_size;

	if(pat_len<target_kmer_size)
	{
		mode=2;
		target_kmer_size=2*kmer_size;
	}

	if(pat_len<target_kmer_size)
		return;

	int min_kmers_overlap=kmer_size;

	int num_kmers=(pat_len-min_kmers_overlap+target_kmer_size-min_kmers_overlap-1)/(target_kmer_size-min_kmers_overlap);

	int total_kmers_size=num_kmers*target_kmer_size;
	int total_kmer_overlaps_size=total_kmers_size-pat_len;

	int cur_pat_start=0;
	vector<Result> res2, res3;

	vector<Result> all_res;

	int i, j;

	for(i=0;i<num_kmers;i++)
	{
		for(j=0;j<=suff_filter->max_kmer_errors;j++)
		{
			vector<Result>& res=GetCandidatesAll(suff_filter, kmer_flags, pat+cur_pat_start, res2, res3, mode, j);

			if(j<suff_filter->max_kmer_errors && (int)res.size()<suff_filter->max_part_res) continue;

			all_res.reserve(all_res.size()+res.size());

			int k;
			for(k=0;k<(int)res.size();k++)
			{
				Result r;
				//r.ref_pos=res[k];
				//r.cur_pat_start=cur_pat_start;
				r.ref_pos=res[k].ref_pos;
				r.cur_pat_start=cur_pat_start+res[k].cur_pat_start;
				all_res.push_back(r);
			}

			break;
		}

		int cur_overlap=0; if(num_kmers>1) cur_overlap=(total_kmer_overlaps_size+i)/(num_kmers-1);
		cur_pat_start+=target_kmer_size-cur_overlap;
	}

	sort(all_res.begin(), all_res.end());
	int num_res=all_res.size();

	int cur_min_overlap=suff_filter->min_overlap_percentage*pat_len;
	if(cur_min_overlap<suff_filter->min_overlap) cur_min_overlap=suff_filter->min_overlap;

	int min_overlap_size=cur_min_overlap-((suff_filter->match_type==MATCH_SUBSTITUTE_ONLY)?0:suff_filter->max_kmer_errors);

	for(i=0;i<num_res;i++)
	{
		long long ref_pos=all_res[i].ref_pos;
		int cur_pat_start=all_res[i].cur_pat_start;

		if(ref_pos>=0 && ref_pos<suff_filter->total_reads_size && suff_filter->all_reads[ref_pos]!='$' && suff_filter->all_reads[ref_pos]!='#')
		{
			long long end_ref_read=ref_pos;
			long long start_ref_read=ref_pos;

			while(suff_filter->all_reads[start_ref_read]!='$') start_ref_read--; start_ref_read++;
			while(suff_filter->all_reads[end_ref_read]!='$') end_ref_read++;

			bool found_shifts=false;
			int cur_shift1=0;
			int cur_shift2=0;

			if(end_ref_read-ref_pos>=min_overlap_size)
			{
				cur_shift1=(ref_pos-start_ref_read)-cur_pat_start;
				cur_shift2=cur_shift1;
				found_shifts=true;
			}

			for(j=i+1;j<num_res;j++)
			{
				ref_pos=all_res[j].ref_pos;
				cur_pat_start=all_res[j].cur_pat_start;

				if(ref_pos<end_ref_read)
				{
					if(end_ref_read-ref_pos>=min_overlap_size)
					{
						int shift=(ref_pos-start_ref_read)-cur_pat_start;
						if(shift<cur_shift1) cur_shift1=shift;
						if(shift>cur_shift2) cur_shift2=shift;
					}
				}
				else break;
			}

			if(found_shifts)
			{
				PossibleMatch pm;
				pm.str=(char*)(suff_filter->all_reads+start_ref_read);
				pm.qual=0; if(suff_filter->all_quals) pm.qual=(char*)(suff_filter->all_quals+start_ref_read);
				pm.str_len=end_ref_read-start_ref_read;

				pm.shift1=-cur_shift1;
				pm.shift2=-cur_shift2;
				if(pm.shift2<pm.shift1) {int u=pm.shift1;pm.shift1=pm.shift2;pm.shift2=u;}

				pm.max_overlap=ComputeMaxOverlap(pat_len, pm.str_len, pm.shift1, pm.shift2);

				possible_matches.push_back(pm);
			}

			i=j-1;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
