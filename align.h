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

#include "suffixcactus.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class SuffixStruct>
void MatchApproxFast(SuffixStruct* suff_tree, SuffixStruct* suff_tree_rev, char* flag, const char* pat, int pat_len, int match_type, int cur_max_errors, vector<int>& res)
{
	res.clear();
	int num_group_parts=cur_max_errors+1;

	int num_satisfying_parts=4;

	int total_size=suff_tree->seq_len-1;

	int i,j;
	int* part_size=new int[num_group_parts];
	int* part_end=new int[num_group_parts];
	int* part_end_rev=new int[num_group_parts];

	ConstructPartSizes(pat_len, num_group_parts, SMALL_ENDS, part_size, part_end, part_end_rev);

	char* pat_rev=new char[pat_len+1];
	memcpy(pat_rev, pat, pat_len+1);
	Reverse(pat_rev, pat_len);

	SuffixFilterInfo filter_info;

	SuffixFilterInfo filter_info_rev=filter_info;

	filter_info.part_end=part_end;
	filter_info_rev.part_end=part_end_rev;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	vector<ResultLocation> part_res;
	vector<ResultLocation> part_res_rev;

	filter_info.pres=&part_res;
	filter_info_rev.pres=&part_res_rev;

	for(i=0;i<num_group_parts;i++)
	{
		filter_info.num_needed_parts=num_group_parts-i;
		if(filter_info.num_needed_parts>num_satisfying_parts) filter_info.num_needed_parts=num_satisfying_parts;
		filter_info.start_part=i;
		int pat_start=0; if(i) pat_start=part_end[i-1];
		filter_info.pat_start=pat_start;
		MatchApprox(suff_tree, pat, pat_len, match_type, 0, &filter_info);
	}

	for(i=0;i<num_group_parts;i++) // match with filter of 000111222...etc starting from this part
	{
		filter_info_rev.num_needed_parts=num_group_parts-i;
		if(filter_info_rev.num_needed_parts>num_satisfying_parts) filter_info_rev.num_needed_parts=num_satisfying_parts;
		filter_info_rev.start_part=i;
		int pat_start=0; if(i) pat_start=part_end_rev[i-1];
		filter_info_rev.pat_start=pat_start;
		MatchApprox(suff_tree_rev, pat_rev, pat_len, match_type, 0, &filter_info_rev);//, Set<int>& res)
	}

	for(j=0;j<(int)part_res_rev.size();j++)
	{
		int cur_suff_ind=part_res_rev[j].start_suff_ind;
		do{
			int cur_val=suff_tree_rev->suff[cur_suff_ind]+part_res_rev[j].shift;
			if(cur_val<0) cur_val=0; if(cur_val>total_size) cur_val=total_size;

			flag[cur_val]=1;
			cur_suff_ind++;
		}
		while(cur_suff_ind<part_res_rev[j].end_suff_ind);// && cur_suff_ind<part_res_rev[rev_part][j].start_suff_ind+max_part_res);
	}

	for(j=0;j<(int)part_res.size();j++)
	{
		int cur_suff_ind=part_res[j].start_suff_ind;
		do{
			int cur_val=suff_tree->suff[cur_suff_ind]+part_res[j].shift;
			if(cur_val<0) cur_val=0; if(cur_val>total_size) cur_val=total_size;

			if(!(flag[cur_val]&2))
			{
				flag[cur_val]|=2;

				int min_val=-cur_max_errors-pat_len-cur_val+total_size; if(min_val<0) min_val=0;
				int max_val=cur_max_errors-pat_len-cur_val+total_size; if(max_val>total_size) max_val=total_size;

				int k;

				for(k=min_val;k<=max_val;k++) if(flag[k]&1)
				{
					//num_added+=res.add(cur_val);
					res.push_back(cur_val);
					break;
				}
			}

			cur_suff_ind++;
		}
		while(cur_suff_ind<part_res[j].end_suff_ind);
	}

	for(j=0;j<(int)part_res.size();j++)
	{
		int cur_suff_ind=part_res[j].start_suff_ind;
		do{
			int cur_val=suff_tree->suff[cur_suff_ind]+part_res[j].shift;
			if(cur_val<0) cur_val=0; if(cur_val>total_size) cur_val=total_size;

			flag[cur_val]=0;
			cur_suff_ind++;
		}
		while(cur_suff_ind<part_res[j].end_suff_ind);
	}

	for(j=0;j<(int)part_res_rev.size();j++)
	{
		int cur_suff_ind=part_res_rev[j].start_suff_ind;
		do{
			int cur_val=suff_tree_rev->suff[cur_suff_ind]+part_res_rev[j].shift;
			if(cur_val<0) cur_val=0; if(cur_val>total_size) cur_val=total_size;

			flag[cur_val]=0;
			cur_suff_ind++;
		}
		while(cur_suff_ind<part_res_rev[j].end_suff_ind);
	}

	delete[] pat_rev;

	delete[] part_size;
	delete[] part_end;
	delete[] part_end_rev;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ResMap
{
	int start;
	int end;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int GetMappings(int* T, const char* cur_seq, int cur_seq_len, const char* pat, int pat_len, int req_edit, int st, vector<ResMap>& res, int match_type)
{
	int k;

	if(match_type==MATCH_SUBSTITUTE_ONLY)
	{
		int i, edit=req_edit+1;

		for(k=0;k<cur_seq_len-pat_len+1;k++)
		{
			int h=0; for(i=0;i<pat_len;i++) h+=(pat[i]!=cur_seq[k+i]);
			if(h<edit) {edit=h; res.clear();}
			if(h<=edit) {ResMap r; r.start=st+k; r.end=pat_len; res.push_back(r);}
		}

		if(edit>req_edit) return -1;
		return edit;
	}

	vector<int> end_locs;
	int edit=ApproxEditDistanceSearchBInsideA(T, cur_seq, cur_seq_len, pat, pat_len, req_edit, match_type, &end_locs);

	if(edit>=0 || edit<=req_edit)
	{
		char* rev_pat=(char*)AllocateHeap(pat_len);
		char* rev_seq=(char*)AllocateHeap(cur_seq_len);

		memcpy(rev_pat, pat, pat_len);
		Reverse(rev_pat, pat_len);

		memcpy(rev_seq, cur_seq, cur_seq_len);
		Reverse(rev_seq, cur_seq_len);

		for(k=0;k<(int)end_locs.size();k++)
		{
			int cur_end_loc=cur_seq_len-1-(end_locs[k]-1);

			vector<int> start_locs;
			edit=ApproxEditDistanceFlexibleEndA(T, rev_seq+cur_end_loc, cur_seq_len-cur_end_loc, rev_pat, pat_len, req_edit, match_type, &start_locs);

			//if(edit!=i) printf("ErrorA\n"); fflush(NULL);

			int p,q;

			for(p=0;p<(int)start_locs.size();p++)
			{
				int start_loc=cur_seq_len-1-(start_locs[p]-1+cur_end_loc);

				for(q=0;q<(int)res.size();q++)
				{
					if(res[q].start==st+start_loc && res[q].end==end_locs[k]-start_loc)
						break;
				}

				if(q==(int)res.size())
				{
					ResMap r;
					r.start=st+start_loc;
					r.end=end_locs[k]-start_loc;
					res.push_back(r);

					//int u=EditDistance(T, suff_tree->seq+res_starts.back(), res_ends.back(), pat, pat_len);
					//if(u!=i) printf("ErrorB %d/%d\n", u, i); fflush(NULL);
				}
			}
		}

		FreeHeap(rev_pat, pat_len);
		FreeHeap(rev_seq, cur_seq_len);

		//res.push_back(cur_res[j]);
	}

	return edit;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class SuffixStruct>
int GetEditDist(SuffixStruct* suff_tree, SuffixStruct* suff_tree_rev, char* flag, int* T, const char* pat, int pat_len,
					int max_edit_dist, int min_part_size, int match_type, vector<ResMap>& res)
{
	res.clear();

	int i, j, k;

	int mi=max_edit_dist;
	if(min_part_size>0) mi=pat_len/min_part_size-1;
	if(mi<0) mi=0;
	if(mi>max_edit_dist) mi=max_edit_dist;

	for(i=0;i<=mi;i++)
	{
		if(i==0)
		{
			vector<int> r;
			if(Match(suff_tree, pat, pat_len, r, 1)) // exact results are too much, we dont need them actually
				return i;
		}
		else
		{
			vector<int> cur_res;
			MatchApproxFast(suff_tree, suff_tree_rev, flag, pat, pat_len, match_type, i, cur_res);

			bool find_match=false;

			for(j=0;j<(int)cur_res.size();j++)
			{
				int edit=-1;

				if(match_type==MATCH_SUBSTITUTE_ONLY)
				{
					edit=0;
					for(k=0;k<pat_len;k++)
					{
						int pos=cur_res[j]+k;
						if(pos<1 || pos>=suff_tree->seq_len-2 || suff_tree->seq[pos]=='$') {edit=-1; break;}
						edit+=(pat[k]!=suff_tree->seq[pos]);
					}

					if(edit>=0 && edit<i)
					{
						printf("Error Smaller Edit %d/%d\n%s\n", edit, i, pat);
						for(k=cur_res[j];k<cur_res[j]+pat_len;k++) printf("%c", suff_tree->seq[k]); printf("\n");
					}

					if(edit==i) {find_match=true; ResMap r; r.start=cur_res[j]; r.end=pat_len; res.push_back(r);}
				}
				else
				{
					int st=cur_res[j]-i; if(st<1) st=1;
					int en=cur_res[j]+pat_len+i; if(en>suff_tree->seq_len-2) en=suff_tree->seq_len-2;

					//for(k=st;k<cur_res[j];k++) if(suff_tree->seq[k]=='$') {st=k+1; break;}
					//for(k=st;k<en;k++) if(suff_tree->seq[k]=='$') {en=k; break;}

					for(k=st;k<en;k++) if(suff_tree->seq[k]=='$') break;
					if(k<en) {if(en-k<k-st) en=k; else st=k+1;}

					const char* cur_seq=suff_tree->seq+st;
					int cur_seq_len=en-st;

					int edit=GetMappings(T, cur_seq, cur_seq_len, pat, pat_len, i, st, res, match_type);

					if(edit==i) find_match=true;

					if(edit>=0 && edit<i)
					{
						printf("Error Smaller Edit %d/%d\n%s\n", edit, i, pat);
						for(k=st;k<en;k++) printf("%c", suff_tree->seq[k]); printf("\n");
					}
				}
			}

			if(find_match) return i;
		}
	}

	return -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ProcessAlignRead
{
	char*		read_str;
	int			read_len;
	int			edit;

	vector<ResMap> res;

	int*		T;

	int			done;
	pthread_mutex_t	done_mutex;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SolutionAlign
{
	int				ind;	// thread index

	SuffixCactus*	suff_tree;
	SuffixCactus*	suff_tree_rev;
	char*			flag;

	int				max_edit_dist;
	int				min_part_size;
	int				match_type;

	int				cur_num_reads;
	ProcessAlignRead* pr;

	int				dummy[1000]; // be memory write safe

	void Compute()
	{
		int i;

		for(i=0;i<cur_num_reads;i++)
		{
			int is_done=0;
			pthread_mutex_lock(&pr[i].done_mutex);
			is_done=pr[i].done;
			if(is_done==0) pr[i].done=1; // processing
			pthread_mutex_unlock(&pr[i].done_mutex);

			if(is_done==0)
			{
				pr[i].edit=GetEditDist(suff_tree, suff_tree_rev, flag, pr[i].T, pr[i].read_str, pr[i].read_len, max_edit_dist, min_part_size, match_type, pr[i].res);
			}
		}
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void* threadfunctionAlign(void* arg)
{
	SolutionAlign* sol=(SolutionAlign*)arg;
	sol->Compute();
	return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

long long ComputeRefGenomeSize(const char* file_name_ref)
{
	FILE* file_ref=fopen(file_name_ref, "r"); if(!file_ref) return 0;
	int file_block_size=1024*1024;
	BufferedInFile buf_ref; buf_ref.Initialize(file_ref, file_block_size);
	BlockedString bstr; bstr.Initialize(file_block_size);
	while(1) if(!buf_ref.ReadLine(bstr, '>')) break;
	buf_ref.Destroy(); fclose(file_ref);
	long long file_ref_size=bstr.GetSize();
	return file_ref_size;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class SuffixStruct>
void Align(FastqFiles& file_name_org, const char* file_name_ref, const char* file_name_align, int match_type, int accuracy, int circular, int num_threads, int max_reads_per_step, int file_block_size)
{
	FILE* file_ref=fopen(file_name_ref, "r");
	if(!file_ref) return;
	BufferedInFile buf_ref; buf_ref.Initialize(file_ref, file_block_size);
	BlockedString bstr; bstr.Initialize(file_block_size);
	while(1) if(!buf_ref.ReadLine(bstr, '>')) break;
	buf_ref.Destroy();
	fclose(file_ref);
	long long file_ref_size=bstr.GetSize();
	char* str_ref=(char*)AllocateHeap(2*(file_ref_size+circular+1)+2+1); // note the size includes the \0
	int cur_size=0; str_ref[cur_size++]='$';
	bstr.GetString(str_ref+cur_size, file_ref_size);
	bstr.Destroy();
	cur_size+=file_ref_size;
	memcpy(str_ref+cur_size, str_ref+1, circular);
	cur_size+=circular;
	str_ref[cur_size++]='$';
	memcpy(str_ref+cur_size, str_ref+1, file_ref_size+circular);
	ReverseComplement(str_ref+cur_size, file_ref_size+circular);
	cur_size+=file_ref_size+circular;
	str_ref[cur_size++]='$'; str_ref[cur_size++]='#'; str_ref[cur_size]=0;

	char aa[1<<8]={0};

	int ii;
	for(ii=0;ii<cur_size;ii++)
	{
		if(!aa[(int)str_ref[ii]]) printf("[%c]", str_ref[ii]);
		aa[(int)str_ref[ii]]=1;
	}

	printf(" Size = %lld\n", file_ref_size);
	printf("SuffArray Started\n"); fflush(NULL);

	SuffixStruct suff_tree;
	Build(&suff_tree, str_ref, cur_size);

	printf("SuffArray Finished\n"); fflush(NULL);

	char* str_ref_rev=(char*)AllocateHeap(cur_size+1);
	memcpy(str_ref_rev, str_ref, cur_size+1);
	Reverse(str_ref_rev, cur_size-1);

	SuffixStruct suff_tree_rev;
	Build(&suff_tree_rev, str_ref_rev, cur_size);

	printf("SuffArrayRev Finished\n"); fflush(NULL);

	file_name_org.Reset();

	char* buf=new char[MAX_READ_LEN+1];

	FILE* file_fastq_org=0; if(!file_name_org.GetCurFileName()) return;
	file_fastq_org=fopen(file_name_org.GetCurFileName(), "r");

	char first_char_org=GetFirstChar(file_name_org.GetCurFileName(), file_fastq_org, buf, 0);
	if(first_char_org==0) return;

	long long tt=1;
	int ct=0;
	while(tt<cur_size) {tt*=4; ct++;}

	int max_edit_dist=ED_INF;
	int min_part_size=ct-accuracy; // accuracy=4, make the actual program accuracy=5

	SolutionAlign** psols=new SolutionAlign*[num_threads];
	pthread_t* thread_ID=new pthread_t[num_threads];
	void** exit_status=new void*[num_threads];

	int i;

	long long num_prev_reads=0;
	ProcessAlignRead* pr=new ProcessAlignRead[max_reads_per_step];

	for(i=0;i<num_threads;i++)
	{
		psols[i]=new SolutionAlign;
		psols[i]->ind=i;

		psols[i]->pr=pr;

		psols[i]->suff_tree=&suff_tree;
		psols[i]->suff_tree_rev=&suff_tree_rev;
		psols[i]->flag=(char*)AllocateHeap(cur_size);
		memset(psols[i]->flag, 0, cur_size);

		psols[i]->max_edit_dist=max_edit_dist;
		psols[i]->min_part_size=min_part_size;
		psols[i]->match_type=match_type;
	}

	for(i=0;i<max_reads_per_step;i++)
	{
		pthread_mutex_init(&pr[i].done_mutex, NULL);
	}

	FILE* file_align=fopen(file_name_align, "w");

	char* r1=new char[MAX_READ_LEN+1];

	long long total_num_reads=0;
	long long num_aligned_reads=0;
	long long num_excluded_reads=0;

	long long total_num_bases=0;
	long long num_aligned_bases=0;
	long long num_excluded_bases=0;

	while(1)
	{
		bool end_file=false;

		int num_reads=0;

		while(num_reads<max_reads_per_step)
		{
			r1[0]=0;

			int nr1=GetReadSequence(file_fastq_org, first_char_org, r1, buf, 0, 0, 0, 'N'); //info_line_org[1-cur_info]

			if(nr1==0)
			{
				fclose(file_fastq_org);
				file_name_org.FinishCurFile();

				if(file_name_org.GetCurFileName())
				{
					file_fastq_org=fopen(file_name_org.GetCurFileName(), "r");
					first_char_org=GetFirstChar(file_name_org.GetCurFileName(), file_fastq_org, buf, 0);
					if(first_char_org==0) return;
					nr1=GetReadSequence(file_fastq_org, first_char_org, r1, buf, 0, 0, 0, 'N');
					if(nr1==0) return;
				}
			}

			if(nr1==0)
			{
				end_file=true;
				break;
			}

			pr[num_reads].read_len=nr1;
			pr[num_reads].read_str=(char*)AllocateHeap(nr1+1);
			strcpy(pr[num_reads].read_str, r1);
			pr[num_reads].done=0;
			pr[num_reads].T=(int*)AllocateHeap((nr1+1)*sizeof(int));

			num_reads++;
		}

		for(i=0;i<num_threads;i++)
		{
			psols[i]->cur_num_reads=num_reads;
		}

		for(i=1;i<num_threads;i++)
		{
			pthread_create(&thread_ID[i], NULL, threadfunctionAlign, psols[i]);
		}

		threadfunctionAlign(psols[0]);

		for(i=1;i<num_threads;i++)
		{
			pthread_join(thread_ID[i], &exit_status[i]);
		}

		for(i=0;i<num_reads;i++)
		{
			if(total_num_reads && (total_num_reads%100000==0)) {printf("."); fflush(NULL);}

			int e1=pr[i].edit;

			strcpy(r1, pr[i].read_str);
			int nr1=strlen(r1);

			FreeHeap(pr[i].read_str, nr1+1);
			FreeHeap(pr[i].T, (nr1+1)*sizeof(int));

			if(e1>=0)
			{
				int num_res=pr[i].res.size();
				fprintf(file_align, "%s %d %d", r1, e1, num_res);

				if(e1>0)
				{
					if(num_res==0) printf("UnExpected Error\n");
					int j;
					for(j=0;j<num_res;j++)
					{
						fprintf(file_align, " %d:%d", pr[i].res[j].start, pr[i].res[j].end);
					}
				}

				fprintf(file_align, "\n");

				num_aligned_reads++;
				num_aligned_bases+=nr1;
			}
			else
			{
				num_excluded_reads++;
				num_excluded_bases+=nr1;
			}

			total_num_reads++;
			total_num_bases+=nr1;
		}

		num_prev_reads+=num_reads;
		if(end_file) break;
	}

	printf("\n");
	printf("Reference Genome Size = %lld\n", file_ref_size);
	printf("Average Read Len      = %Lf\n", (long double)total_num_bases/total_num_reads);
	printf("Coverage              = %Lf\n", (long double)total_num_bases/file_ref_size);

	printf("\n");

	printf("Total Num Reads       = %lld\n", total_num_reads);
	printf("Num Aligned Reads     = %lld (%lf %%)\n", num_aligned_reads, (double)100.0*num_aligned_reads/total_num_reads);
	printf("Num Excluded Reads    = %lld (%lf %%)\n", num_excluded_reads, (double)100.0*num_excluded_reads/total_num_reads);

	printf("\n");

	printf("Total Num Bases       = %lld\n", total_num_bases);
	printf("Num Aligned Bases     = %lld (%lf %%)\n", num_aligned_bases, (double)100.0*num_aligned_bases/total_num_bases);
	printf("Num Excluded Bases    = %lld (%lf %%)\n", num_excluded_bases, (double)100.0*num_excluded_bases/total_num_bases);

	fflush(NULL);

	Destroy(&suff_tree);
	Destroy(&suff_tree_rev);

	fclose(file_align);

	FreeHeap(str_ref, 2*(file_ref_size+circular+1)+3);
	FreeHeap(str_ref_rev, 2*(file_ref_size+circular+1)+3);

	for(i=0;i<max_reads_per_step;i++)
	{
		pthread_mutex_destroy(&pr[i].done_mutex);
	}

	for(i=0;i<num_threads;i++)
	{
		FreeHeap(psols[i]->flag, cur_size);
		delete psols[i];
	}

	delete[] psols;
	delete[] thread_ID;
	delete[] exit_status;

	delete[] pr;

	delete[] r1;
	delete[] buf;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AlignMain(int argc, char* argv[])
{
	printf("Alignment Started\n"); fflush(NULL);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	int			default_match_type=MATCH_EDIT_DIST;
	int			default_accuracy=5;
	int			default_circular=0;

	int			default_num_threads=16;
	int			default_max_reads_per_step=1000;
	int			default_file_block_size=1024*1024;

	const char* default_input_dir="./";

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	int			match_type=default_match_type; bool set_match_type; set_match_type=false;
	int			accuracy=default_accuracy; //bool_set_accuracy; set_accuracy=false;
	int			circular=default_circular; //bool_set_circular; set_circular=false;

	int			num_threads=default_num_threads; //bool_set_num_threads; set_num_threads=false;
	int			max_reads_per_step=default_max_reads_per_step; //bool_set_max_reads_per_step; set_max_reads_per_step=false;
	int			file_block_size=default_file_block_size; //bool_set_file_block_size; set_file_block_size=false;

	char 		input_dir[MAX_PATH_LEN+1]; strcpy(input_dir, default_input_dir); //bool_set_input_dir; set_input_dir=false;

	vector<FileName> user_input_file_names; bool set_user_input_file_names; set_user_input_file_names=false;
	char		user_ref_genome_file_name[MAX_PATH_LEN+1]; user_ref_genome_file_name[0]=0; bool set_user_ref_genome_file_name; set_user_ref_genome_file_name=false;
	char		user_align_file_name[MAX_PATH_LEN+1]; user_align_file_name[0]=0; bool set_user_align_file_name; set_user_align_file_name=false;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// GetData from user here, be sure that slash is added to directory names

	int i;

	for(i=2;i<argc;i++)
	{
		bool processed=false;
		char* str=argv[i];

		int j, len=strlen(str);
		for(j=0;j<len;j++) if(str[j]=='=') break;

		if(j<len-1)
		{
			char* opt=(char*)malloc(j+1);
			memcpy(opt, str, j); opt[j]=0;
			char* val=str+j+1;

			if(strcmp(opt,"-matchtype")==0)
			{
				if(strcmp(val,"edit")==0) {match_type=MATCH_EDIT_DIST; processed=true; set_match_type=true;}
				else if(strcmp(val,"hamming")==0) {match_type=MATCH_SUBSTITUTE_ONLY; processed=true; set_match_type=true;}
				else if(strcmp(val,"insdel")==0) {match_type=MATCH_INSERT_DELETE_ONLY; processed=true; set_match_type=true;}
			}
			else if(strcmp(opt,"-threads")==0)
			{
				sscanf(val, "%d", &num_threads);
				processed=true; //set_num_threads=true;
			}
			else if(strcmp(opt,"-readsperstep")==0) // ADVANCED
			{
				sscanf(val, "%d", &max_reads_per_step);
				processed=true; //set_max_reads_per_step=true;
			}
			else if(strcmp(opt,"-accuracy")==0) // ADVANCED
			{
				sscanf(val, "%d", &accuracy);
				processed=true; //set_accuracy=true;
			}
			else if(strcmp(opt,"-circular")==0)
			{
				sscanf(val, "%d", &circular);
				processed=true; //set_circular=true;
			}
			else if(strcmp(opt,"-inputdir")==0)
			{
				strcpy(input_dir, val); AddSlash(input_dir);
				processed=true; //set_input_dir=true;
			}
			else if(strcmp(opt,"-inputfile")==0) // REQUIRED
			{
				FileName f; f.s[0]=0;
				strcpy(f.s, val);
				user_input_file_names.push_back(f);
				processed=true; set_user_input_file_names=true;
			}
			else if(strcmp(opt,"-refgenomefile")==0) // REQUIRED
			{
				strcpy(user_ref_genome_file_name, val);
				processed=true; set_user_ref_genome_file_name=true;
			}
			else if(strcmp(opt,"-alignfile")==0) // REQUIRED
			{
				strcpy(user_align_file_name, val);
				processed=true; set_user_align_file_name=true;
			}

			free(opt);
		}

		if(!processed)
		{
			printf("Unknown option [%s].\n", str);
			fflush(NULL);
			return;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(num_threads<1) num_threads=1;
	if(max_reads_per_step<100) max_reads_per_step=100;
	if(circular<0) circular=0;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*******************************************************************************************************

	FastqFiles fastq_file_names;
	char ref_genome_file_name[MAX_PATH_LEN+1];

	FileName file;

	for(i=0;i<(int)user_input_file_names.size();i++)
	{
		GetFileName(user_input_file_names[i].s, input_dir, file.s);
		fastq_file_names.AddFile(file);
	}

	GetFileName(user_ref_genome_file_name, input_dir, ref_genome_file_name);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*******************************************************************************************************

	bool accepted=true;

	if(!set_user_input_file_names)
	{
		printf("Please specify fasta or fastq input file(s) using \"-inputfile=filename\"\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_match_type)
	{
		printf("You must specify explicitly the matching type (the \"-matchtype\" option).\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_user_ref_genome_file_name)
	{
		printf("Please specify input reference genome file using \"-refgenomefile=filename\"\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_user_align_file_name)
	{
		printf("Please specify output alignment file using \"-alignfile=filename\"\n"); fflush(NULL);
		accepted=false;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(i=0;accepted && i<(int)fastq_file_names.file_names.size();i++)
	{
		FILE* file=fopen(fastq_file_names.file_names[i].s, "r");
		if(!file)
		{
			printf("Cannot open input file [%s] for reading.\n", fastq_file_names.file_names[i].s); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(accepted)
	{
		FILE* file=fopen(ref_genome_file_name, "r");
		if(!file)
		{
			printf("Cannot open reference genome file [%s] for reading.\n", ref_genome_file_name); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(accepted)
	{
		FILE* file=fopen(user_align_file_name, "w");
		if(!file)
		{
			printf("Cannot open alignment file [%s] for writing.\n", user_align_file_name); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(!accepted)
		return;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*******************************************************************************************************

	Align<SuffixCactus>(fastq_file_names, ref_genome_file_name, user_align_file_name, match_type, accuracy, circular, num_threads, max_reads_per_step, file_block_size);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Alignment Finished\n\n"); fflush(NULL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AlignHelp()
{
	printf("-Run \"karect -align [options_list]\" for the alignment tool.\n");
	printf("-Available options: (i=integer, f=file, s=directory)\n");
	printf("-Essential options:\n");
	printf("  \"-matchtype=[edit|hamming|insdel]\": Specify the matching type. \"hamming\" allows substitution errors only.\n"
		   "                        \"edit\" allows insertions, deletion, and substitutions with equal costs.\n"
		   "                        \"insdel\" is the same as \"edit\", but the cost of substitutions is doubled.\n");
	printf("  \"-inputfile=f\":     Specify an input fasta/fastq file. This option can be repeated for multiple files.\n");
	printf("  \"-refgenomefile=f\": Specify the file containing the reference genome (to be aligned with).\n");
	printf("  \"-alignfile=f\":     Specify the output alignment file.\n");
	printf("-Basic options:\n");
	printf("  \"-inputdir=s\":      Specify the input directory. Ignored if input file paths are complete [Default=.].\n");
	printf("  \"-threads=i\":       Specify the number of threads [Default=16].\n");
	printf("-Advanced options:\n");
	printf("  \"-circular=i\":      Specify the sequence size to be appended circularly (for circular genomes) [Default=0].\n");
	printf("  \"-accuracy=i\":      Specify the alignment accuracy [Default=5].\n");
	printf("  \"-readsperstep=i\":  Specify the maximum number of processed reads per step [Default=1000].\n");
	printf("-Example:\n"
		   "        ./karect -align -inputdir=/sra_data -inputfile=SRR001666_1.fasta -inputfile=SRR001666_2.fasta\n"
		   "              -refgenomefile=NC_000913.fna -alignfile=./SRR001666_align.txt -matchtype=hamming -threads=12\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
