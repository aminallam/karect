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

enum{CELL_HAPLOID, CELL_DIPLOID};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CorrectRead(SuffixFilter* suff_filter, ErrorGraph* eg, unsigned short* kmer_flags, char* read_str, int read_len, char* read_quality, char* result_str, int& result_len, char* result_quality, int& num_truncated)
{
	vector<PossibleMatch> possible_matches;

	ComputeFilter(suff_filter, kmer_flags, read_str, read_len, possible_matches);

	long long max_cands=suff_filter->max_mul_len_matches; if(read_len) max_cands/=read_len;
	if(max_cands<3*suff_filter->max_part_res) max_cands=3*suff_filter->max_part_res;

	if((long long)possible_matches.size()>max_cands)
	{
		// The difference between with and without sorting is most probably because of the comparison orders (tie case)
		sort(possible_matches.begin(), possible_matches.end());
		possible_matches.resize(max_cands);
	}

	FixReadUsingGraph(eg, kmer_flags, suff_filter->kmer_size, read_str, read_len, read_quality, possible_matches, result_str, result_len, result_quality, num_truncated);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ProcessRead
{
	char*		read_str;
	int			read_len;
	char*		read_quality;

	char*		result_str;
	int			result_len;
	char*		result_quality;

	int			buf_eg_size;
	char*		buf_eg;

	int			done;
	pthread_mutex_t	done_mutex;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Solution
{
	int					ind;	// thread index

	int					cur_num_reads;
	SuffixFilter*		suff_filter;
	ErrorGraph*			eg;

	unsigned short*		kmer_flags;

	ProcessRead*		pr;

	int					dummy[1000]; // be memory write safe

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
				int num_truncated=0;

				bool confirmed=false;

				if(eg->load_mode)
				{
					LoadErrorGraph(eg, pr[i].buf_eg);
					FreeHeap(pr[i].buf_eg, pr[i].buf_eg_size);
					confirmed=eg->confirmed;
				}

				if(confirmed) pr[i].result_len=0;
				else CorrectRead(suff_filter, eg, kmer_flags, pr[i].read_str, pr[i].read_len, pr[i].read_quality, pr[i].result_str, pr[i].result_len, pr[i].result_quality, num_truncated);

				if(eg->save_mode)
				{
					pr[i].buf_eg_size=GetErrorGraphSize(eg);
					pr[i].buf_eg=(char*)AllocateHeap(pr[i].buf_eg_size);
					SaveErrorGraph(eg, pr[i].buf_eg);
				}
				else
				{
					if(pr[i].result_len==0)
					{
						pr[i].result_len=pr[i].read_len; strcpy(pr[i].result_str, pr[i].read_str);
						if(pr[i].read_quality) strcpy(pr[i].result_quality, pr[i].read_quality);
					}
				}
			}
		}
	}
};

void* threadfunction(void* arg)
{
	Solution* sol=(Solution*)arg;
	sol->Compute();
	return NULL;
}

void CorrectErrors(FastqFiles& fastq_file_name, const char* reads_file_name, const char* graph_a_file_name, const char* graph_b_file_name, FastqFiles& res_file_name,
					SuffixFilter* suff_filter, long long num_reads_without_complements, int global_max_read_len, int num_threads, int file_block_size, int max_reads_per_step, bool is_fastq)//, int inc)//, int start_read, int max_num_reads)
{
	int i;

	long long num_prev_reads=0;
	ProcessRead* pr=new ProcessRead[max_reads_per_step];

	FILE* in_graph=0;
	BufferedInFile buf_in_graph;

	if(graph_a_file_name)
	{
		in_graph=fopen(graph_a_file_name, "rb");
		buf_in_graph.Initialize(in_graph, file_block_size);
	}

	FILE* file_reads=fopen(reads_file_name, "rb");
	BufferedInFile buf_reads;
	buf_reads.Initialize(file_reads, file_block_size);

	FILE* out_file=0;
	BufferedOutFile buf_out_file;

	FILE* out_graph=0;
	BufferedOutFile buf_out_graph;

	int cur_info=0;
	char* info_line[2]; info_line[0]=new char[MAX_READ_LEN+1]; info_line[1]=new char[MAX_READ_LEN+1];
	char* buf=new char[MAX_READ_LEN+1];
	char* quality_info=new char[MAX_READ_LEN+1];

	FILE* fastq_file=0;
	char fastq_first_char=0;

	if(graph_b_file_name)
	{
		out_graph=fopen(graph_b_file_name, "wb");
		buf_out_graph.Initialize(out_graph, file_block_size);
	}
	else
	{
		fastq_file_name.Reset();
		res_file_name.Reset();

		if(!fastq_file_name.GetCurFileName()) return;
		fastq_file=fopen(fastq_file_name.GetCurFileName(), "r");

		fastq_first_char=GetFirstChar(fastq_file_name.GetCurFileName(), fastq_file, buf, info_line[cur_info]);
		if(fastq_first_char==0) return;

		if(!res_file_name.GetCurFileName()) return;

		out_file=fopen(res_file_name.GetCurFileName(), "w");
		buf_out_file.Initialize(out_file, file_block_size);
	}

	Solution** psols=new Solution*[num_threads];
	pthread_t* thread_ID=new pthread_t[num_threads];
	void** exit_status=new void*[num_threads];

	for(i=0;i<num_threads;i++)
	{
		psols[i]=new Solution;
		psols[i]->ind=i;

		psols[i]->suff_filter=suff_filter;

		psols[i]->pr=pr;

		ErrorGraph* eg=new ErrorGraph;
		InitErrorGraph(eg, graph_a_file_name!=0, graph_b_file_name!=0, global_max_read_len, suff_filter->min_overlap, suff_filter->min_overlap_percentage, suff_filter->max_overlap_error_rate, suff_filter->truncate_factor, suff_filter->min_straight_edge_val, suff_filter->min_straight_edge_percentage, suff_filter->min_read_weight, suff_filter->max_mul_len_matches, suff_filter->match_type);//, suff_filter->max_num_traces, suff_filter->match_type);
		psols[i]->eg=eg;

		long long num_kmers=((long long)1)<<(suff_filter->kmer_size*BITS_PER_SYMBOL);
		psols[i]->kmer_flags=(unsigned short*)AllocateHeap(num_kmers*sizeof(unsigned short));
		memset(psols[i]->kmer_flags, 0, num_kmers*sizeof(unsigned short));
	}

	for(i=0;i<max_reads_per_step;i++)
	{
		pr[i].read_str=(char*)AllocateHeap(global_max_read_len+1);
		pr[i].read_quality=0; if(is_fastq) pr[i].read_quality=(char*)AllocateHeap(global_max_read_len+1);
		pr[i].result_str=(char*)AllocateHeap(2*global_max_read_len+1);
		pr[i].result_quality=0; if(is_fastq) pr[i].result_quality=(char*)AllocateHeap(2*global_max_read_len+1);
		pthread_mutex_init(&pr[i].done_mutex, NULL);
	}

	char* seq=new char[MAX_READ_LEN+1];

	long long one_percent=(num_reads_without_complements+99)/100;
	long long last_percent=one_percent;

	printf("Each following dot indicates that 1%% of the required processing at current stage has been finished:\n"); fflush(NULL);
	//	    ....................................................................................................

	while(1)
	{
		while(num_prev_reads>=last_percent)
		{
			last_percent+=one_percent;
			printf(".");
			fflush(NULL);
		}

		bool end_file=false;

		int num_reads=0;

		while(num_reads<max_reads_per_step)
		{
			int seq_len=buf_reads.ReadUntil(seq, '$');
			if(seq_len==0) {end_file=true; break;}

			pr[num_reads].read_len=seq_len;
			strcpy(pr[num_reads].read_str, seq);

			if(is_fastq)
			{
				buf_reads.ReadNumBytes(seq, seq_len); seq[seq_len]=0;
				strcpy(pr[num_reads].read_quality, seq);
			}

			pr[num_reads].result_len=0;
			pr[num_reads].done=0;

			if(in_graph)
			{
				pr[num_reads].buf_eg_size=0;
				buf_in_graph.ReadNumBytes(&pr[num_reads].buf_eg_size, sizeof(pr[num_reads].buf_eg_size));
				pr[num_reads].buf_eg=(char*)AllocateHeap(pr[num_reads].buf_eg_size);
				buf_in_graph.ReadNumBytes(pr[num_reads].buf_eg, pr[num_reads].buf_eg_size);
			}

			num_reads++;
		}

		/////////////////////////////////////////////////////////////////////////////

		for(i=0;i<num_threads;i++)
		{
			psols[i]->cur_num_reads=num_reads;
		}

		for(i=1;i<num_threads;i++)
		{
			pthread_create(&thread_ID[i], NULL, threadfunction, psols[i]);
		}

		threadfunction(psols[0]);

		for(i=1;i<num_threads;i++)
		{
			pthread_join(thread_ID[i], &exit_status[i]);
		}

		/////////////////////////////////////////////////////////////////////////////

		if(graph_b_file_name)
		{
			for(i=0;i<num_reads;i++)
			{
				buf_out_graph.WriteNumBytes(&pr[i].buf_eg_size, sizeof(pr[i].buf_eg_size));
				buf_out_graph.WriteNumBytes(pr[i].buf_eg, pr[i].buf_eg_size);
				FreeHeap(pr[i].buf_eg, pr[i].buf_eg_size);
			}
		}
		else
		{
			for(i=0;i<num_reads;i++)
			{
				int seq_len=GetReadSequence(fastq_file, fastq_first_char, seq, buf, info_line[1-cur_info], quality_info, 0, 'N');

				if(seq_len==0)
				{
					fclose(fastq_file);
					fastq_file=0;
					fastq_file_name.FinishCurFile();

					buf_out_file.Destroy();
					fclose(out_file);
					out_file=0;
					res_file_name.FinishCurFile();

					if(fastq_file_name.GetCurFileName())
					{
						fastq_file=fopen(fastq_file_name.GetCurFileName(), "r");
						cur_info=0;
						fastq_first_char=GetFirstChar(fastq_file_name.GetCurFileName(), fastq_file, buf, info_line[cur_info]);
						if(fastq_first_char==0) return;
						seq_len=GetReadSequence(fastq_file, fastq_first_char, seq, buf, info_line[1-cur_info], quality_info, 0, 'N');
						if(seq_len==0) return;

						if(!res_file_name.GetCurFileName()) return;
						out_file=fopen(res_file_name.GetCurFileName(), "w");
						buf_out_file.Initialize(out_file, file_block_size);
					}
				}
				if(seq_len==0) {printf("UnExpected End of File\n"); fflush(NULL); break;}

				char ch_end_line='\n';

				ModifyInfoLength(info_line[cur_info], strlen(pr[i].result_str));
				buf_out_file.WriteNumBytes(info_line[cur_info], strlen(info_line[cur_info])); buf_out_file.WriteNumBytes(&ch_end_line, 1);
				buf_out_file.WriteNumBytes(pr[i].result_str, strlen(pr[i].result_str)); buf_out_file.WriteNumBytes(&ch_end_line, 1);

				if(is_fastq)
				{
					ModifyInfoLength(quality_info, strlen(pr[i].result_str));
					buf_out_file.WriteNumBytes(quality_info, strlen(quality_info)); buf_out_file.WriteNumBytes(&ch_end_line, 1);
					buf_out_file.WriteNumBytes(pr[i].result_quality, strlen(pr[i].result_quality)); buf_out_file.WriteNumBytes(&ch_end_line, 1);
				}

				cur_info=1-cur_info;
			}
		}

		num_prev_reads+=num_reads;

		if(end_file) break;
	}

	for(i=0;i<max_reads_per_step;i++)
	{
		FreeHeap(pr[i].read_str, global_max_read_len+1);
		if(is_fastq) FreeHeap(pr[i].read_quality, global_max_read_len+1);
		FreeHeap(pr[i].result_str, 2*global_max_read_len+1);
		if(is_fastq) FreeHeap(pr[i].result_quality, 2*global_max_read_len+1);

		pthread_mutex_destroy(&pr[i].done_mutex);
	}

	for(i=0;i<num_threads;i++)
	{
		DestroyErrorGraph(psols[i]->eg);
		delete psols[i]->eg;

		long long num_kmers=((long long)1)<<(suff_filter->kmer_size*BITS_PER_SYMBOL);
		FreeHeap(psols[i]->kmer_flags, num_kmers*sizeof(unsigned short));

		delete psols[i];
	}

	delete[] psols;
	delete[] thread_ID;
	delete[] exit_status;

	buf_reads.Destroy();
	fclose(file_reads);

	if(in_graph)
	{
		buf_in_graph.Destroy();
		fclose(in_graph);
	}

	if(out_graph)
	{
		buf_out_graph.Destroy();
		fclose(out_graph);
	}

	if(fastq_file)
	{
		fclose(fastq_file);
	}

	if(out_file)
	{
		buf_out_file.Destroy();
		fclose(out_file);
	}

	delete[] pr;

	delete[] buf;
	delete[] seq;
	delete[] info_line[0]; delete[] info_line[1];
	delete[] quality_info;

	printf("\n"); fflush(NULL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CorrectErrors(long long global_num_reads, long long global_num_kmers, long long global_total_size, int global_max_read_len, bool is_fastq,
					FastqFiles& fastq_file_name, const char* file_name, FastqFiles& res_file_name, int insert_complement,// int start_read, int max_num_reads,
					const char* graph_a_file_name, const char* graph_b_file_name,
					int min_overlap, double min_overlap_percentage, double max_overlap_error_rate, double truncate_factor, double min_straight_edge_val,
					double min_straight_edge_percentage, double min_read_weight, long long max_mul_len_matches, int max_kmer_slots, int max_part_res, int kmer_size, int max_kmer_errors, int kmer_with_errors_size, int read_index_bits,
					int match_type, long long memory_budget, int num_threads, int file_block_size, long long cache_block_size, int max_reads_per_step, int use_qual)
{
	printf("Global Num Reads       = %lld\n", global_num_reads/(insert_complement+1));
	printf("Global Num Bases       = %lld\n", (global_total_size-1-global_num_reads)/(insert_complement+1));
	printf("Global Num Kmers       = %lld\n", global_num_kmers/(insert_complement+1));
	printf("Global Max Read Len    = %d\n", global_max_read_len);
	printf("Global Avg Read Len    = %Lf\n", ((long double)global_total_size-1-global_num_reads)/global_num_reads);

	printf("Min Overlap            = %d\n", min_overlap);
	printf("Min Overlap Percentage = %lf\n", min_overlap_percentage);
	printf("Max Overlap Error Rate = %lf\n", max_overlap_error_rate);
	printf("Truncate Factor        = %lf\n", truncate_factor);
	printf("Min Straight Edge Val  = %lf\n", min_straight_edge_val);
	printf("Min Straight Edge Per  = %lf\n", min_straight_edge_percentage);
	printf("Min Read Weight        = %lf\n", min_read_weight);
	printf("Max Mul Len Matches    = %lld\n", max_mul_len_matches);
	printf("Max Kmer Slots         = %d\n", max_kmer_slots);
	printf("Use Qual               = %d\n", use_qual);

	printf("Match Type             = %s\n", (match_type==MATCH_EDIT_DIST)?"Edit Distance":((match_type==MATCH_SUBSTITUTE_ONLY)?"Hamming Distance":"Insert/Delete Only"));
	printf("Kmer Size              = %d\n", kmer_size);
	printf("Max Kmer Errors        = %d\n", max_kmer_errors);
	printf("Kmer With Errors Size  = %d\n", kmer_with_errors_size);
	printf("Read Index Bits        = %d\n", read_index_bits);

	printf("Num Threads            = %d\n", num_threads);
	printf("Memory Budget          = "); PrintMemory(memory_budget); printf("\n");
	printf("File Block Size        = "); PrintMemory(file_block_size); printf("\n");
	printf("Cache Block Size       = "); PrintMemory(cache_block_size); printf("\n");
	printf("Max Reads Per Step     = %d\n", max_reads_per_step);
	printf("Max Part Res           = %d\n", max_part_res);

	fflush(NULL);

	long long global_needed_memory=global_num_kmers*GetOverheadPerKmer()+global_total_size*GetOverheadPerChar(is_fastq?1:0);

	long long working_budget=(long long)50*1024*1024+2*file_block_size;
	working_budget+=GetFixedOverhead(kmer_size, num_threads, global_max_read_len);

	long double kmer_overhead=(long double)global_needed_memory/global_num_kmers;

	long long max_num_kmers=(memory_budget-working_budget)/kmer_overhead;
	if(max_num_kmers<200) max_num_kmers=200;

	int num_parts=(global_num_kmers+max_num_kmers-1)/max_num_kmers;

	max_num_kmers=(global_num_kmers+num_parts-1)/num_parts;
	if(max_num_kmers<200) max_num_kmers=200;

	long long end_file=(global_total_size-1)/(insert_complement+1);
	if(is_fastq) end_file+=(global_total_size-1-global_num_reads)/(insert_complement+1);

	long long start_file=0;

	int cur_step=0;
	const char* file_name_graph[2]={graph_a_file_name, graph_b_file_name};

	while(start_file<end_file)
	{
		printf("Main CorrectErrors StartFile=%lld\n", start_file);

		SuffixFilter suff_filter;

		start_file=Build(&suff_filter, file_name, file_block_size, cache_block_size, start_file, max_num_kmers, insert_complement,
										min_overlap, min_overlap_percentage, max_overlap_error_rate, truncate_factor, min_straight_edge_val, min_straight_edge_percentage, min_read_weight, max_mul_len_matches, max_kmer_slots, max_part_res,
										kmer_size, max_kmer_errors, kmer_with_errors_size, read_index_bits, match_type, num_threads, is_fastq, use_qual);

		const char* file_graph_a=0; if(cur_step) file_graph_a=file_name_graph[cur_step%2];
		const char* file_graph_b=0; if(start_file<end_file) file_graph_b=file_name_graph[(cur_step+1)%2];

		CorrectErrors(fastq_file_name, file_name, file_graph_a, file_graph_b, res_file_name, &suff_filter, global_num_reads/(insert_complement+1), global_max_read_len, num_threads, file_block_size, max_reads_per_step, is_fastq);//, 1+insert_complement);//, start_read, max_num_reads);

		Destroy(&suff_filter);

		cur_step++;
	}

	if(cur_step)
	{
		FILE* delete_graph=fopen(file_name_graph[0], "wb"); fclose(delete_graph);
		delete_graph=fopen(file_name_graph[1], "wb"); fclose(delete_graph);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool CorrectErrorsMain(FastqFiles& file_name, const char* input_file_name, FastqFiles& res_file_name, int insert_complement,// int start_read, int max_num_reads,
					const char* graph_a_file_name, const char* graph_b_file_name,
					int min_overlap, double min_overlap_percentage, double max_overlap_error_rate, double truncate_factor, double min_straight_edge_val,
					double min_straight_edge_percentage, double min_read_weight, long long max_mul_len_matches, int max_kmer_slots, int kmer_factor, int max_part_res, int kmer_size, int max_kmer_errors, int match_type, int cell_type,
					long long memory_budget, int num_threads, int file_block_size, long long cache_block_size, int max_reads_per_step, int estimate_coverage, int use_qual, double aggressive)
{
	long long global_num_reads;
	long long global_num_kmers;
	long long global_total_size;
	int global_max_read_len;
	bool is_fastq;

	if(!GetGlobalReadStats(file_name, input_file_name, file_block_size, insert_complement, kmer_size,
							global_num_reads, global_num_kmers, global_total_size, global_max_read_len, is_fastq))
		return false;

	double avg_read_len=((long double)global_total_size-1-global_num_reads)/global_num_reads;
	if(min_overlap>0.7*avg_read_len) min_overlap=0.7*avg_read_len;

	if(estimate_coverage)
	{
		double coverage=EstimateCoverage(input_file_name, file_block_size, insert_complement, is_fastq);
		coverage*=(global_total_size-1-global_num_reads)/global_num_kmers;

		printf("Estimated Coverage     = %lf\n", coverage);

		if(coverage<10) coverage=10;
		coverage*=aggressive;
		if(cell_type==CELL_DIPLOID) coverage*=0.5;
		if(coverage<min_straight_edge_val) min_straight_edge_val=coverage;
	}

	if(max_part_res<1.5*min_straight_edge_val) max_part_res=1.5*min_straight_edge_val;

	int read_index_bits=0;
	long long yaa=1;
	while(yaa<global_total_size) {yaa*=2; read_index_bits++;}

	if(kmer_factor>0)
	{
		while(kmer_size<14)
		{
			long long val=1; val=(val<<(2*kmer_size))*kmer_factor;
			if(val>global_num_kmers) break;
			kmer_size++;
			//4^k > total_num_kmers/max_total_num_part
		}
	}

	int kmer_with_errors_size=kmer_size; if(match_type!=MATCH_SUBSTITUTE_ONLY) kmer_with_errors_size+=max_kmer_errors;
	read_index_bits=((read_index_bits+7)/8)*8;
	if(read_index_bits<32) read_index_bits=32;

	InitializeKmerEntry<unsigned int, unsigned long long>(kmer_with_errors_size*BITS_PER_SYMBOL, read_index_bits);

	max_mul_len_matches/=(global_total_size-1-global_num_reads)/global_num_reads;

	CorrectErrors(global_num_reads, global_num_kmers, global_total_size, global_max_read_len, is_fastq, file_name, input_file_name, res_file_name, insert_complement, graph_a_file_name, graph_b_file_name, min_overlap, min_overlap_percentage, max_overlap_error_rate, truncate_factor, min_straight_edge_val, min_straight_edge_percentage, min_read_weight, max_mul_len_matches, max_kmer_slots, max_part_res, kmer_size, max_kmer_errors, kmer_with_errors_size, read_index_bits, match_type, memory_budget, num_threads, file_block_size, cache_block_size, max_reads_per_step, use_qual);

	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CorrectErrorsMain(int argc, char* argv[])
{
	printf("Error Correction Started\n"); fflush(NULL);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	double		default_aggressive=0.42;
	int			default_min_overlap=35;
	double		default_min_overlap_percentage=0.20;
	double		default_max_overlap_error_rate_first_stage=0.25;
	double		default_max_overlap_error_rate_second_stage=0.25;
	double		default_truncate_factor=2.5;
	double		default_min_straight_edge_val=100.0;
	int			default_estimate_coverage=1;
	int			default_use_qual=1;
	double		default_min_straight_edge_percentage=1.0;
	double		default_min_read_weight=0;
	long long	default_max_mul_len_matches=2*1000; default_max_mul_len_matches*=1000*1000;
	int			default_max_kmer_slots=100000;
	int			default_kmer_factor=1000;

	double		default_high_error_rate=0.50;
	double		default_high_error_read_weight=0;

	int			default_insert_complement=1;
	int			default_match_type=MATCH_EDIT_DIST;
	int			default_cell_type=CELL_HAPLOID;
	int			default_kmer_size=9;
	int			default_max_kmer_errors=2;

	long long	default_memory_budget=10*1024*1024; default_memory_budget*=1024*1024;
	int			default_num_threads=16;
	int			default_file_block_size=10; default_file_block_size*=1024*1024;
	long long	default_cache_block_size=128; default_cache_block_size*=1024*1024;
	int			default_max_reads_per_step=1000;
	int			default_max_part_res=30;

	int			default_allow_trimming=0;
	int			default_num_stages=1;

	const char* default_input_dir="./";
	const char* default_working_dir="./";
	const char* default_result_dir="./";
	const char* default_result_prefix="karect_";

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	double		aggressive=default_aggressive;
	int			min_overlap=default_min_overlap; //bool_set_min_overlap; set_min_overlap=false;
	double		min_overlap_percentage=default_min_overlap_percentage; //bool_set_min_overlap_percentage; set_min_overlap_percentage=false;
	double		max_overlap_error_rate_first_stage=default_max_overlap_error_rate_first_stage; //bool_set_max_overlap_error_rate_first_stage; set_max_overlap_error_rate_first_stage=false;
	double		max_overlap_error_rate_second_stage=default_max_overlap_error_rate_second_stage; bool set_max_overlap_error_rate_second_stage; set_max_overlap_error_rate_second_stage=false;
	double		truncate_factor=default_truncate_factor; //bool_set_truncate_factor; set_truncate_factor=false;
	double		min_straight_edge_val=default_min_straight_edge_val; //bool_set_min_straight_edge_val; set_min_straight_edge_val=false;
	int			estimate_coverage=default_estimate_coverage;
	int			use_qual=default_use_qual;
	double		min_straight_edge_percentage=default_min_straight_edge_percentage; //bool_set_min_straight_edge_percentage; set_min_straight_edge_percentage=false;
	double		min_read_weight=default_min_read_weight; //bool_set_min_read_weight; set_min_read_weight=false;
	long long	max_mul_len_matches=default_max_mul_len_matches;
	int			max_kmer_slots=default_max_kmer_slots; //bool_set_max_kmer_slots; set_max_kmer_slots=false; // NOT USED NOW
	int			kmer_factor=default_kmer_factor;

	int			insert_complement=default_insert_complement; //bool_set_insert_complement; set_insert_complement=false;
	int			match_type=default_match_type; bool set_match_type; set_match_type=false;
	int			cell_type=default_cell_type; bool set_cell_type; set_cell_type=false;
	int			kmer_size=default_kmer_size; //bool_set_kmer_size; set_kmer_size=false;
	int			max_kmer_errors=default_max_kmer_errors; //bool_set_max_kmer_errors; set_max_kmer_errors=false;

	long long	memory_budget=default_memory_budget; //bool_set_memory_budget; set_memory_budget=false;
	int			num_threads=default_num_threads; //bool_set_num_threads; set_num_threads=false;
	int			file_block_size=default_file_block_size; //bool_set_file_block_size; set_file_block_size=false;
	long long	cache_block_size=default_cache_block_size; //bool_set_cache_block_size; set_cache_block_size=false;
	int			max_reads_per_step=default_max_reads_per_step; //bool_set_max_reads_per_step; set_max_reads_per_step=false;
	int			max_part_res=default_max_part_res; //bool_set_max_part_res; set_max_part_res=false;

	int			allow_trimming=default_allow_trimming; //bool set_allow_trimming; set_allow_trimming=false;
	int			num_stages=default_num_stages; //bool_set_num_stages; set_num_stages=false;

	char 		input_dir[MAX_PATH_LEN+1]; strcpy(input_dir, default_input_dir); //bool_set_input_dir; set_input_dir=false;
	char 		working_dir[MAX_PATH_LEN+1]; strcpy(working_dir, default_working_dir); //bool_set_working_dir; set_working_dir=false;
	char 		result_dir[MAX_PATH_LEN+1]; strcpy(result_dir, default_result_dir); //bool_set_result_dir; set_result_dir=false;
	char 		result_prefix[MAX_PATH_LEN+1]; strcpy(result_prefix, default_result_prefix); //bool_set_result_prefix; set_result_prefix=false;

	vector<FileName> user_input_file_names; bool set_user_input_file_names; set_user_input_file_names=false;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*******************************************************************************************************
	// GetData from user here, be sure that slash is added to directory names

	int i;

	for(i=2;i<argc;i++)
	{
		bool processed=false;
		char* str=argv[i];

		if(strcmp(str,"-nocomp")==0) {insert_complement=0; processed=true;}// set_insert_complement=true;} // ADVANCED
		else if(strcmp(str,"-higherror")==0)
		{
			max_overlap_error_rate_first_stage=max_overlap_error_rate_second_stage=default_high_error_rate;
			min_read_weight=default_high_error_read_weight;
			processed=true;
		}
		else
		{
			int j, len=strlen(str);
			for(j=0;j<len;j++) if(str[j]=='=') break;

			if(j<len-1)
			{
				char* opt=(char*)malloc(j+1);
				memcpy(opt, str, j); opt[j]=0;
				char* val=str+j+1;

				if(strcmp(opt,"-trim")==0)
				{
					if(strcmp(val,"yes")==0) {allow_trimming=1; processed=true;}// set_allow_trimming=true;}
					else if(strcmp(val,"no")==0) {allow_trimming=0; processed=true;}// set_allow_trimming=true;}
				}
				else if(strcmp(opt,"-estcov")==0)
				{
					if(strcmp(val,"yes")==0) {estimate_coverage=1; processed=true;}// set_allow_trimming=true;}
					else if(strcmp(val,"no")==0) {estimate_coverage=0; processed=true;}// set_allow_trimming=true;}
				}
				else if(strcmp(opt,"-usequal")==0)
				{
					if(strcmp(val,"yes")==0) {use_qual=1; processed=true;}// set_allow_trimming=true;}
					else if(strcmp(val,"no")==0) {use_qual=0; processed=true;}// set_allow_trimming=true;}
				}
				else if(strcmp(opt,"-aggressive")==0)
				{
					sscanf(val, "%lf", &aggressive);
					processed=true; //set_min_overlap_percentage=true;
				}
				else if(strcmp(opt,"-minoverlap")==0)
				{
					sscanf(val, "%d", &min_overlap);
					processed=true; //set_min_overlap=true;
				}
				else if(strcmp(opt,"-minoverlapper")==0)
				{
					sscanf(val, "%lf", &min_overlap_percentage);
					processed=true; //set_min_overlap_percentage=true;
				}
				else if(strcmp(opt,"-minreadweigth")==0)
				{
					sscanf(val, "%lf", &min_read_weight);
					processed=true; //set_min_read_weight=true;
				}
				else if(strcmp(opt,"-errorrate")==0)
				{
					sscanf(val, "%lf", &max_overlap_error_rate_first_stage);
					processed=true; //set_max_overlap_error_rate_first_stage=true;
				}
				else if(strcmp(opt,"-errorratesec")==0)
				{
					sscanf(val, "%lf", &max_overlap_error_rate_second_stage);
					processed=true; set_max_overlap_error_rate_second_stage=true;
				}
				else if(strcmp(opt,"-trimfact")==0)
				{
					sscanf(val, "%lf", &truncate_factor);
					processed=true; //set_truncate_factor=true;
				}
				else if(strcmp(opt,"-reserveval")==0)
				{
					sscanf(val, "%lf", &min_straight_edge_val);
					processed=true; //set_min_straight_edge_val=true;
				}
				else if(strcmp(opt,"-reserveper")==0)
				{
					sscanf(val, "%lf", &min_straight_edge_percentage);
					processed=true; //set_min_straight_edge_percentage=true;
				}
				else if(strcmp(opt,"-matchtype")==0)
				{
					if(strcmp(val,"edit")==0) {match_type=MATCH_EDIT_DIST; processed=true; set_match_type=true;}
					else if(strcmp(val,"hamming")==0) {match_type=MATCH_SUBSTITUTE_ONLY; processed=true; set_match_type=true;}
					else if(strcmp(val,"insdel")==0) {match_type=MATCH_INSERT_DELETE_ONLY; processed=true; set_match_type=true;}
				}
				else if(strcmp(opt,"-celltype")==0)
				{
					if(strcmp(val,"haploid")==0) {cell_type=CELL_HAPLOID; processed=true; set_cell_type=true;}
					else if(strcmp(val,"diploid")==0) {cell_type=CELL_DIPLOID; processed=true; set_cell_type=true;}
				}
				else if(strcmp(opt,"-kmer")==0)
				{
					sscanf(val, "%d", &kmer_size);
					processed=true; //set_kmer_size=true;
				}
				else if(strcmp(opt,"-kmererrors")==0) // ADVANCED
				{
					sscanf(val, "%d", &max_kmer_errors);
					processed=true; //set_max_kmer_errors=true;
				}
				else if(strcmp(opt,"-memory")==0)
				{
					double memory;
					sscanf(val, "%lf", &memory);
					memory_budget=memory*1024; memory_budget*=1024*1024;
					processed=true; //set_memory_budget=true;
				}
				else if(strcmp(opt,"-threads")==0)
				{
					sscanf(val, "%d", &num_threads);
					processed=true; //set_num_threads=true;
				}
				else if(strcmp(opt,"-fbs")==0) // ADVANCED
				{
					int fbs;
					sscanf(val, "%d", &fbs);
					file_block_size=fbs; file_block_size*=1024*1024;
					processed=true; //set_file_block_size=true;
				}
				else if(strcmp(opt,"-cbs")==0) // ADVANCED
				{
					int cbs;
					sscanf(val, "%d", &cbs);
					cache_block_size=cbs; cache_block_size*=1024*1024;
					processed=true; //set_cache_block_size=true;
				}
				else if(strcmp(opt,"-readsperstep")==0) // ADVANCED
				{
					sscanf(val, "%d", &max_reads_per_step);
					processed=true; //set_max_reads_per_step=true;
				}
				else if(strcmp(opt,"-maxlenmatches")==0) // ADVANCED
				{
					int v;
					sscanf(val, "%d", &v);
					max_mul_len_matches=v;
					max_mul_len_matches*=1000*1000;
					processed=true; //set_max_kmer_slots=true;
				}
				else if(strcmp(opt,"-maxkmerslots")==0) // ADVANCED
				{
					sscanf(val, "%d", &max_kmer_slots);
					processed=true; //set_max_kmer_slots=true;
				}
				else if(strcmp(opt,"-kmerfactor")==0) // ADVANCED
				{
					sscanf(val, "%d", &kmer_factor);
					processed=true; //set_kmer_factor=true;
				}
				else if(strcmp(opt,"-maxkmerres")==0) // ADVANCED
				{
					sscanf(val, "%d", &max_part_res);
					processed=true; //set_max_part_res=true;
				}
				else if(strcmp(opt,"-numstages")==0)
				{
					sscanf(val, "%d", &num_stages);
					processed=true; //set_num_stages=true;
				}
				else if(strcmp(opt,"-inputdir")==0)
				{
					strcpy(input_dir, val); AddSlash(input_dir);
					processed=true; //set_input_dir=true;
				}
				else if(strcmp(opt,"-tempdir")==0)
				{
					strcpy(working_dir, val); AddSlash(working_dir);
					processed=true; //set_working_dir=true;
				}
				else if(strcmp(opt,"-resultdir")==0)
				{
					strcpy(result_dir, val); AddSlash(result_dir);
					processed=true; //set_result_dir=true;
				}
				else if(strcmp(opt,"-resultprefix")==0)
				{
					strcpy(result_prefix, val);
					processed=true; //set_result_prefix=true;
				}
				else if(strcmp(opt,"-inputfile")==0) // REQUIRED
				{
					FileName f; f.s[0]=0;
					strcpy(f.s, val);
					user_input_file_names.push_back(f);
					processed=true; set_user_input_file_names=true;
				}

				free(opt);
			}
		}

		if(!processed)
		{
			printf("Unknown option [%s].\n", str);
			fflush(NULL);
			return;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(!allow_trimming)
		truncate_factor=0;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(truncate_factor<0) truncate_factor=0;
	if(min_straight_edge_val<0) min_straight_edge_val=0;
	if(min_straight_edge_percentage<0) min_straight_edge_percentage=0;
	if(memory_budget<100*1024*1024) memory_budget=100*1024*1024;
	if(num_threads<1) num_threads=1;
	if(max_part_res<10) max_part_res=10;
	if(file_block_size<1024*1024) file_block_size=1024*1024;
	if(cache_block_size<1024*1024) cache_block_size=1024*1024;
	if(max_reads_per_step<100) max_reads_per_step=100;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*******************************************************************************************************

	FastqFiles	fastq_file_names;
	FastqFiles	res_temp_file_names;
	FastqFiles	res_file_names;

	char graph_a_file_name[MAX_PATH_LEN+1];
	char graph_b_file_name[MAX_PATH_LEN+1];

	char input_file_name[MAX_PATH_LEN+1];

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	FileName file;

	for(i=0;i<(int)user_input_file_names.size();i++)
	{
		GetFileName(user_input_file_names[i].s, input_dir, file.s);
		fastq_file_names.AddFile(file);

		const char* user_file=fastq_file_names.file_names[i].s;

		int last_slash=strlen(user_file)-1;
		while(last_slash>=0 && user_file[last_slash]!='\\' && user_file[last_slash]!='/') last_slash--;

		int file_name_start=last_slash+1;
		if(file_name_start>=(int)strlen(user_file)) continue;

		file.s[0]=0; strcpy(file.s, working_dir);
		strcpy(file.s+strlen(file.s), "temp_res_");
		strcpy(file.s+strlen(file.s), user_file+file_name_start);
		res_temp_file_names.AddFile(file);

		file.s[0]=0; strcpy(file.s, result_dir);
		strcpy(file.s+strlen(file.s), result_prefix);
		strcpy(file.s+strlen(file.s), user_file+file_name_start);
		res_file_names.AddFile(file);
	}

	int len_working_dir=strlen(working_dir);

	strcpy(graph_a_file_name, working_dir);
	strcpy(graph_a_file_name+len_working_dir, "res_graph_a.txt");

	strcpy(graph_b_file_name, working_dir);
	strcpy(graph_b_file_name+len_working_dir, "res_graph_b.txt");

	strcpy(input_file_name, working_dir);
	strcpy(input_file_name+len_working_dir, "input_file.txt");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool accepted=true;

	if(num_stages==1 && set_max_overlap_error_rate_second_stage)
	{
		printf("You can not set error rate of second stage when \"numstages\" = 1.\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_match_type)
	{
		printf("You must specify explicitly the matching type (the \"-matchtype\" option).\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_cell_type)
	{
		printf("You must specify explicitly the cell type (the \"-celltype\" option).\n"); fflush(NULL);
		accepted=false;
	}
	else if(min_overlap<10)
	{
		printf("Value of \"minoverlap\" must be >= 10.\n"); fflush(NULL);
		accepted=false;
	}
	else if(max_overlap_error_rate_first_stage<0.01 || max_overlap_error_rate_second_stage<0.01 ||
			max_overlap_error_rate_first_stage>0.75 || max_overlap_error_rate_second_stage>0.75)
	{
		printf("Value of \"errorrate\" must be in [0.01, 0.75].\n"); fflush(NULL);
		accepted=false;
	}
	else if(kmer_size<8 || kmer_size>14)
	{
		printf("Value of \"kmer\" must be in [8, 14].\n"); fflush(NULL);
		accepted=false;
	}
	else if(max_kmer_errors<0 || max_kmer_errors>2)
	{
		printf("Value of \"kmererrors\" must be in [0, 2].\n"); fflush(NULL);
		accepted=false;
	}
	else if(num_stages<1 || num_stages>2)
	{
		printf("Value of \"numstages\" must be in [1, 2].\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_user_input_file_names)
	{
		printf("Please specify fasta or fastq input file(s) using \"-inputfile=filename\".\n"); fflush(NULL);
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

	for(i=0;accepted && i<(int)res_temp_file_names.file_names.size();i++)
	{
		FILE* file=fopen(res_temp_file_names.file_names[i].s, "w");
		if(!file)
		{
			printf("Cannot open temporary result file [%s] for writing.\n", res_temp_file_names.file_names[i].s); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	for(i=0;accepted && i<(int)res_file_names.file_names.size();i++)
	{
		FILE* file=fopen(res_file_names.file_names[i].s, "w");
		if(!file)
		{
			printf("Cannot open result file [%s] for writing.\n", res_file_names.file_names[i].s); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(accepted)
	{
		FILE* file=fopen(input_file_name, "wb");
		if(!file)
		{
			printf("Cannot open temporary file [%s] for writing.\n", input_file_name); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(accepted)
	{
		FILE* file=fopen(graph_a_file_name, "wb");
		if(!file)
		{
			printf("Cannot open temporary file [%s] for writing.\n", graph_a_file_name); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(accepted)
	{
		FILE* file=fopen(graph_b_file_name, "wb");
		if(!file)
		{
			printf("Cannot open temporary file [%s] for writing.\n", graph_b_file_name); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(!accepted)
		return;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(num_stages==1)
	{
		if(!CorrectErrorsMain(fastq_file_names, input_file_name, res_file_names, insert_complement, graph_a_file_name, graph_b_file_name, min_overlap, min_overlap_percentage, max_overlap_error_rate_second_stage, truncate_factor, min_straight_edge_val, min_straight_edge_percentage, min_read_weight, max_mul_len_matches, max_kmer_slots, kmer_factor, max_part_res, kmer_size, max_kmer_errors, match_type, cell_type, memory_budget, num_threads, file_block_size, cache_block_size, max_reads_per_step, estimate_coverage, use_qual, aggressive)) return;
	}
	else
	{
		if(!CorrectErrorsMain(fastq_file_names, input_file_name, res_temp_file_names, insert_complement, graph_a_file_name, graph_b_file_name, min_overlap, min_overlap_percentage, max_overlap_error_rate_first_stage, truncate_factor, min_straight_edge_val, min_straight_edge_percentage, min_read_weight, max_mul_len_matches, max_kmer_slots, kmer_factor, max_part_res, kmer_size, max_kmer_errors, match_type, cell_type, memory_budget, num_threads, file_block_size, cache_block_size, max_reads_per_step, estimate_coverage, use_qual, aggressive)) return;

		if(!CorrectErrorsMain(res_temp_file_names, input_file_name, res_file_names, insert_complement, graph_a_file_name, graph_b_file_name, min_overlap, min_overlap_percentage, max_overlap_error_rate_second_stage, truncate_factor, min_straight_edge_val, min_straight_edge_percentage, min_read_weight, max_mul_len_matches, max_kmer_slots, kmer_factor, max_part_res, kmer_size, max_kmer_errors, match_type, cell_type, memory_budget, num_threads, file_block_size, cache_block_size, max_reads_per_step, estimate_coverage, use_qual, aggressive)) return;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Error Correction Finished\n\n"); fflush(NULL);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CorrectErrorsHelp()
{
	printf("-Run \"karect -correct [options_list]\" for the error correction tool.\n");
	printf("-Available options: (i=integer, d=decimal, f=file, s=directory)\n");
	printf("-Essential options:\n");
	printf("  \"-inputfile=f\":      Specify an input fasta/fastq file. This option can be repeated for multiple files\n");
	printf("  \"-celltype=[haploid|diploid]\": Specify the cell type. Use \"haploid\" for bacteria and viruses.\n");
	printf("  \"-matchtype=[edit|hamming|insdel]\": Specify the matching type. \"hamming\" allows substitution errors only.\n"
		   "                         \"edit\" allows insertions, deletion, and substitutions with equal costs.\n"
		   "                         \"insdel\" is the same as \"edit\", but the cost of substitutions is doubled.\n"
		   "                         Use \"hamming\" for Illumina datasets, and \"edit\" for 454 datasets.\n");
	printf("-Basic options:\n");
	printf("  \"-inputdir=s\":       Specify the input directory. Ignored if input file paths are complete [Default=.].\n");
	printf("  \"-resultdir=s\":      Specify a directory to save result file(s) [Default=.].\n");
	printf("  \"-resultprefix=s\":   Specify a prefix string to the result file(s) [Default=karect_].\n");
	printf("  \"-tempdir=s\":        Specify a directory to save temporary output files [Default=.].\n");
	printf("  \"-threads=i\":        Specify the number of threads [Default=16].\n");
	printf("  \"-memory=d\":         Specify an upper bound on the memory that can be used in gigabytes [Default=10240.0].\n");
	printf("-Advanced options:\n");
	printf("  \"-aggressive=d\":     Specify the aggressiveness towards error correction [Default=0.42].\n");
	printf("  \"-numstages=i\":      Specify the number of stages (1 or 2) [Default=1].\n");
	printf("  \"-minoverlap=i\":     Specify the minimum overlap size [Default=35].\n");
	printf("  \"-minoverlapper=d\":  Specify the minimum overlap percentage [Default=0.20].\n");
	printf("  \"-minreadweigth=d\":  Specify the minimum read weight [Default=0].\n");
	printf("  \"-errorrate=d\":      Specify the first stage maximum allowed error rate [Default=0.25].\n");
	printf("  \"-errorratesec=d\":   Specify the second stage maximum allowed error rate [Default=0.25].\n");
	printf("  \"-reserveval=d\":     Specify the minimum reservation value [Default=100.0].\n");
	printf("  \"-estcov=[yes|no]\":  Estimate coverage and use it to adjust the minimum reservation value [Default=yes].\n");
	printf("  \"-usequal=[yes|no]\": Use quality values of candidate reads [Default=yes].\n");
	printf("  \"-higherror\":        Work in high error rate mode (error rate = 0.50).\n");
	printf("  \"-trimfact=d\":       Specify the trimming factor [Default=2.5].\n");
	printf("  \"-reserveper=d\":     Specify the minimum reservation percentage [Default=1.0].\n");
	printf("  \"-kmer=i\":           Specify the minimum kmer size (will increase according to \'-kmerfactor\') [Default=9].\n");
	printf("  \"-trim=[yes|no]\":    Allow/Disallow trimming. Do not allow trimming if evaluating results afterwards, or\n"
		   "                         if you will pass results to an assembler which expects fixed read sizes [Default=no].\n");
	printf("-More advanced options:\n");
	printf("  \"-maxlenmatches=i\":  Specify the maximum number of expected alignment computations (millions) [Default=2000].\n");
	printf("  \"-maxkmerslots=i\":   Specify the maximum number kmer slots to be used [Default=100,000].\n");
	printf("  \"-kmerfactor=i\":     Specify the factor f such that 4^kmersize > total_num_kmers/f [Default=1000].\n");
	printf("  \"-maxkmerres=i\":     Specify the maximum number kmer results to be used [Default=30].\n");
	printf("  \"-kmererrors=i\":     Specify the maximum allowed kmer errors (0,1,2) [Default=2].\n");
	printf("  \"-readsperstep=i\":   Specify the maximum number of processed reads per step [Default=1000].\n");
	printf("  \"-fbs=i\":            Specify file block size in megabytes [Default=10].\n");
	printf("  \"-cbs=i\":            Specify cache block size in megabytes [Default=128].\n");
	printf("-Example:\n"
		   "        ./karect -correct -inputdir=/sra_data -inputfile=SRR001666_1.fasta -inputfile=SRR001666_2.fasta\n"
		   "              -resultprefix=karect_r1_ -celltype=haploid -matchtype=hamming -errorrate=0.25 -threads=12\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
