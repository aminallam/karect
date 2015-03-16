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

struct ProcessEvalRead
{
	char*		org_str;
	int			org_len;
	int			org_edit;
	vector<ResMap> org_res;

	char*		corrected_str;
	int			corrected_len;

	int			corrected_edit; // min edit distance to any of the org_res locations
	int			num_edit_opers; // number of edit operations, that is, edit distance between corrected and org
	int			num_trimmed;
	int			org_edit_after_trim;

	int*		T;

	int			done;
	pthread_mutex_t	done_mutex;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SolutionEval
{
	int				ind;	// thread index

	char*			str_ref;
	int				str_ref_size;

	int				allow_trimming;
	int				match_type;

	int				cur_num_reads;
	ProcessEvalRead* pr;

	int				dummy[1000]; // be memory write safe

	void Compute()
	{
		int i,j,k;

		for(i=0;i<cur_num_reads;i++)
		{
			int is_done=0;
			pthread_mutex_lock(&pr[i].done_mutex);
			is_done=pr[i].done;
			if(is_done==0) pr[i].done=1; // processing
			pthread_mutex_unlock(&pr[i].done_mutex);

			if(is_done==0)
			{
				pr[i].num_trimmed=0;
				pr[i].org_edit_after_trim=pr[i].org_edit;

				if(!allow_trimming)
				{
					if(match_type==MATCH_SUBSTITUTE_ONLY && pr[i].org_len>pr[i].corrected_len)
					{
						pr[i].num_edit_opers=EditDistance(pr[i].T, pr[i].org_str, pr[i].org_len, pr[i].corrected_str, pr[i].corrected_len, match_type);

						pr[i].corrected_edit=ED_INF;
						if(pr[i].org_edit==0) pr[i].corrected_edit=pr[i].num_edit_opers;

						for(j=0;j<(int)pr[i].org_res.size();j++)
						{
							int edit=EditDistance(pr[i].T, str_ref+pr[i].org_res[j].start, pr[i].org_res[j].end, pr[i].corrected_str, pr[i].corrected_len, match_type);
							if(edit<pr[i].corrected_edit) pr[i].corrected_edit=edit;
						}
					}
					else if(match_type==MATCH_SUBSTITUTE_ONLY && pr[i].org_len!=pr[i].corrected_len)
					{
						pr[i].num_edit_opers=pr[i].corrected_edit=pr[i].org_len;
					}
					else
					{
						pr[i].num_edit_opers=EditDistance(pr[i].T, pr[i].org_str, pr[i].org_len, pr[i].corrected_str, pr[i].corrected_len, match_type);

						pr[i].corrected_edit=ED_INF;
						if(pr[i].org_edit==0) pr[i].corrected_edit=pr[i].num_edit_opers;

						for(j=0;j<(int)pr[i].org_res.size();j++)
						{
							int edit=EditDistance(pr[i].T, str_ref+pr[i].org_res[j].start, pr[i].org_res[j].end, pr[i].corrected_str, pr[i].corrected_len, match_type);
							if(edit<pr[i].corrected_edit) pr[i].corrected_edit=edit;
						}
					}
				}
				else
				{
					if(pr[i].org_edit==0)
					{
						vector<ResMap> res;
						int edit=GetMappings(pr[i].T, pr[i].org_str, pr[i].org_len, pr[i].corrected_str, pr[i].corrected_len,
												ED_INF, 0, res, match_type);

						int best_num_trimmed=ED_INF;

						for(k=0;k<(int)res.size();k++)
						{
							int num_trimmed=pr[i].org_len-res[k].end;
							if(num_trimmed<best_num_trimmed) best_num_trimmed=num_trimmed;
						}

						pr[i].num_trimmed=best_num_trimmed;
						pr[i].num_edit_opers=edit;
						pr[i].corrected_edit=edit;
					}
					else
					{
						int best_edit_corrected=ED_INF;
						vector<ResMap> best_res_corrected;

						for(j=0;j<(int)pr[i].org_res.size();j++)
						{
							vector<ResMap> res;
							int edit=GetMappings(pr[i].T, str_ref+pr[i].org_res[j].start, pr[i].org_res[j].end, pr[i].corrected_str, pr[i].corrected_len,
													best_edit_corrected, pr[i].org_res[j].start, res, match_type);
							if(edit>=0 && edit<=best_edit_corrected)
							{
								if(edit<best_edit_corrected) best_res_corrected.clear();
								for(k=0;k<(int)res.size();k++)
								{
									best_edit_corrected=edit; //best_trimming_corrected=num_trimmed;
									best_res_corrected.push_back(res[k]);
								}
							}
						}

						int best_edit_org=ED_INF;
						vector<ResMap> best_res_org;

						for(j=0;j<(int)best_res_corrected.size();j++)
						{
							vector<ResMap> res;
							int edit=GetMappings(pr[i].T, pr[i].org_str, pr[i].org_len, str_ref+best_res_corrected[j].start, best_res_corrected[j].end,
													ED_INF, 0, res, match_type);
							if(edit>=0 && edit<=best_edit_org)
							{
								if(edit<best_edit_org) best_res_org.clear(); //best_trimming_org=ED_INF;}
								for(k=0;k<(int)res.size();k++)
								{
									best_edit_org=edit; //best_trimming_org=num_trimmed;
									best_res_org.push_back(res[k]);
								}
							}
						}

						pr[i].num_trimmed=ED_INF;
						pr[i].corrected_edit=best_edit_corrected;
						pr[i].num_edit_opers=ED_INF;
						pr[i].org_edit_after_trim=best_edit_org;

						for(j=0;j<(int)best_res_org.size();j++)
						{
							int num_edit_opers=EditDistance(pr[i].T, pr[i].org_str+best_res_org[j].start, best_res_org[j].end, pr[i].corrected_str, pr[i].corrected_len, match_type);

							if(num_edit_opers<=pr[i].num_edit_opers)
							{
								int num_trimmed=pr[i].org_len-best_res_org[j].end;
								if(num_edit_opers<pr[i].num_edit_opers || num_trimmed<pr[i].num_trimmed) pr[i].num_trimmed=num_trimmed;
								if(num_edit_opers<pr[i].num_edit_opers || num_edit_opers<pr[i].num_edit_opers) pr[i].num_edit_opers=num_edit_opers;
							}
						}
					}
				}
			}
		}
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void* threadfunctionEval(void* arg)
{
	SolutionEval* sol=(SolutionEval*)arg;
	sol->Compute();
	return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//"TN) num_correct_stayed_correct"
//"TP) num_wrong_correctly_fixed"
//"FP) num_correct_wrongly_fixed"
//"FN) num_wrong_stayed_wrong"

void GetMeasures(double TN, double TP, double FP, double FN,
				double& Recall_Sensitivity, double& Specificity, double& Precision, double& FScore, double& Gain)
{
	Recall_Sensitivity=TP/(TP+FN); // num_wrong_correctly_fixed/num_wrong
	Specificity=TN/(TN+FP); // num_correct_stayed_correct/num_correct
	Precision=TP/(TP+FP); // num_wrong_correctly_fixed/num_changed(num_changed=num_fixed_correctly_or_incorrectly)
	FScore=2*Precision*Recall_Sensitivity/(Precision+Recall_Sensitivity); // combine recall and precision
	Gain=(TP-FP)/(TP+FN); //  combine sensitivity and specificity
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class SuffixStruct>
void Evaluate(FastqFiles& file_name_org, FastqFiles& file_name_fixed, const char* file_name_ref, const char* file_name_align, const char* file_name_eval, int match_type, int allow_trimming, int circular, int num_threads, int max_reads_per_step, int file_block_size)
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

	printf("Reference Genome Size = %lld\n", file_ref_size);

	file_name_org.Reset();
	file_name_fixed.Reset();

	FILE* file_fastq_org=0; if(!file_name_org.GetCurFileName()) return;
	file_fastq_org=fopen(file_name_org.GetCurFileName(), "r");

	FILE* file_fastq_fixed=0; if(!file_name_fixed.GetCurFileName()) return;
	file_fastq_fixed=fopen(file_name_fixed.GetCurFileName(), "r");

	char* buf=new char[MAX_READ_LEN+1];

	char first_char_org=GetFirstChar(file_name_org.GetCurFileName(), file_fastq_org, buf, 0);
	if(first_char_org==0) return;

	char first_char_fixed=GetFirstChar(file_name_fixed.GetCurFileName(), file_fastq_fixed, buf, 0);
	if(first_char_fixed==0) return;

	SolutionEval** psols=new SolutionEval*[num_threads];
	pthread_t* thread_ID=new pthread_t[num_threads];
	void** exit_status=new void*[num_threads];

	int i;

	long long num_prev_reads=0;
	ProcessEvalRead* pr=new ProcessEvalRead[max_reads_per_step];

	for(i=0;i<num_threads;i++)
	{
		psols[i]=new SolutionEval;
		psols[i]->ind=i;

		psols[i]->pr=pr;

		psols[i]->str_ref=str_ref;
		psols[i]->str_ref_size=cur_size;

		psols[i]->allow_trimming=allow_trimming;
		psols[i]->match_type=match_type;
	}

	for(i=0;i<max_reads_per_step;i++)
	{
		pthread_mutex_init(&pr[i].done_mutex, NULL);
	}

	FILE* file_eval=fopen(file_name_eval, "w");

	char* r1=new char[MAX_READ_LEN+1];
	char* r2=new char[MAX_READ_LEN+1];

	long long total_num_reads=0;

	long long total_num_bases=0;
	long long total_num_bases_from_corrected_reads=0;

	long long num_trimmed_bases=0;
	long long num_trimmed_bases_from_corrected_reads=0;

	long long num_correct_stayed_correct=0;
	long long num_wrong_correctly_fixed=0;
	long long num_correct_wrongly_fixed=0;
	long long num_wrong_stayed_wrong=0;

	long long num_correct_stayed_correct_bases=0;
	long long num_wrong_correctly_fixed_bases=0;
	long long num_correct_wrongly_fixed_bases=0;
	long long num_wrong_stayed_wrong_bases=0;

	char* aligned=new char[MAX_READ_LEN+1];
	int aligned_len=-1;
	int aligned_edit=-1;
	vector<ResMap> aligned_res;
	int aligned_num_res=0;

	FILE* file_align=fopen(file_name_align, "r");
	if(!file_align) return;

	if(EOF==fscanf(file_align, "%s %d %d ", aligned, &aligned_edit, &aligned_num_res)) aligned_len=-1;
	else{aligned_len=strlen(aligned); aligned_res.clear(); for(i=0;i<aligned_num_res;i++) {ResMap r; if(EOF!=fscanf(file_align, "%d:%d ", &r.start, &r.end)) aligned_res.push_back(r);}}

	while(aligned_len>0)
	{
		bool end_file=false;

		int num_reads=0;

		while(num_reads<max_reads_per_step)
		{
			r1[0]=0; r2[0]=0;

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

			int nr2=GetReadSequence(file_fastq_fixed, first_char_fixed, r2, buf, 0, 0, 0, 'N');

			if(nr2==0)
			{
				fclose(file_fastq_fixed);
				file_name_fixed.FinishCurFile();

				if(file_name_fixed.GetCurFileName())
				{
					file_fastq_fixed=fopen(file_name_fixed.GetCurFileName(), "r");
					first_char_fixed=GetFirstChar(file_name_fixed.GetCurFileName(), file_fastq_fixed, buf, 0);
					if(first_char_fixed==0) return;
					nr2=GetReadSequence(file_fastq_fixed, first_char_fixed, r2, buf, 0, 0, 0, 'N');
					if(nr2==0) return;
				}
			}

			if(nr1==0 || nr2==0)
			{
				end_file=true;
				break;
			}

			if(strcmp(aligned, r1)==0)
			{
				int maxr=nr1; if(nr2>maxr) maxr=nr2;
				for(i=0;i<aligned_num_res;i++) if(aligned_res[i].end>maxr) maxr=aligned_res[i].end;

				pr[num_reads].org_len=nr1;
				pr[num_reads].org_str=(char*)AllocateHeap(nr1+1);
				strcpy(pr[num_reads].org_str, r1);

				pr[num_reads].corrected_len=nr2;
				pr[num_reads].corrected_str=(char*)AllocateHeap(nr2+1);
				strcpy(pr[num_reads].corrected_str, r2);

				pr[num_reads].done=0;
				pr[num_reads].T=(int*)AllocateHeap((maxr+1)*sizeof(int));

				pr[num_reads].org_edit=aligned_edit;
				pr[num_reads].org_res=aligned_res;

				num_reads++;

				if(EOF==fscanf(file_align, "%s %d %d ", aligned, &aligned_edit, &aligned_num_res)) {aligned_len=-1; end_file=true; break;}
				else{aligned_len=strlen(aligned); aligned_res.clear(); for(i=0;i<aligned_num_res;i++) {ResMap r; if(EOF!=fscanf(file_align, "%d:%d ", &r.start, &r.end)) aligned_res.push_back(r);}}
			}
		}

		for(i=0;i<num_threads;i++)
		{
			psols[i]->cur_num_reads=num_reads;
		}

		for(i=1;i<num_threads;i++)
		{
			pthread_create(&thread_ID[i], NULL, threadfunctionEval, psols[i]);
		}

		threadfunctionEval(psols[0]);

		for(i=1;i<num_threads;i++)
		{
			pthread_join(thread_ID[i], &exit_status[i]);
		}

		for(i=0;i<num_reads;i++)
		{
			int e1=pr[i].org_edit;
			int e2=pr[i].corrected_edit;
			int trimmed=pr[i].num_trimmed;
			int num_opers=pr[i].num_edit_opers;
			strcpy(r1, pr[i].org_str); int nr1=strlen(r1);
			strcpy(r2, pr[i].corrected_str); int nr2=strlen(r2);

			int tn=0, tp=0, fp=0, fn=0;

			int x=2*pr[i].org_edit_after_trim, y=2*e2, z=2*num_opers;

			//x=tp+fn
			//y=fp+fn
			//z=tp+fp

			//x-z=fn-fp
			//x-z+y=2*fn

			fn=(x-z+y)/2;
			tp=x-fn;
			fp=z-tp;

			num_wrong_correctly_fixed_bases+=tp;
			num_correct_wrongly_fixed_bases+=fp;
			num_wrong_stayed_wrong_bases+=fn;
			num_correct_stayed_correct_bases+=tn;
			// TN is somehow undefined in edit distance, and hence specificity is also undefined

			int j, maxr=nr1; if(nr2>maxr) maxr=nr2;
			for(j=0;j<(int)pr[i].org_res.size();j++) if(pr[i].org_res[j].end>maxr) maxr=pr[i].org_res[j].end;

			FreeHeap(pr[i].org_str, nr1+1);
			FreeHeap(pr[i].corrected_str, nr2+1);
			FreeHeap(pr[i].T, (maxr+1)*sizeof(int));

			int type=-1;

			if(total_num_reads && (total_num_reads%100000==0)) {printf("."); fflush(NULL);}

			if(e1==0 && e2==0) {type=0; num_correct_stayed_correct++;}
			if(e1>0 && e2==0) {type=1; num_wrong_correctly_fixed++;}
			if(e1==0 && e2>0) {type=2; num_correct_wrongly_fixed++;}
			if(e1>0 && e2>0) {type=3; num_wrong_stayed_wrong++;}

			num_trimmed_bases+=trimmed;
			if(e2==0) num_trimmed_bases_from_corrected_reads+=trimmed;

			char cha[]="TTFF";
			char chb[]="NPPN";

			fprintf(file_eval, "%lld) %s[%d] %s[%d] {%c%c,%d,%d}\n", total_num_reads, r1, e1, r2, e2, cha[type], chb[type], num_opers, trimmed);

			total_num_reads++;

			total_num_bases+=nr1;
			if(e2==0) total_num_bases_from_corrected_reads+=nr1;
		}

		num_prev_reads+=num_reads;
		if(end_file) break;
	}

	printf("\n");

	num_wrong_correctly_fixed_bases/=2;
	num_correct_wrongly_fixed_bases/=2;
	num_wrong_stayed_wrong_bases/=2;

	num_correct_stayed_correct_bases=total_num_bases-num_trimmed_bases-
			num_wrong_correctly_fixed_bases-num_correct_wrongly_fixed_bases-num_wrong_stayed_wrong_bases;

	if(!allow_trimming)
	{
		if(num_trimmed_bases || num_trimmed_bases_from_corrected_reads)
			printf("Unexpected Error Eval.\n");
	}

	double Recall_Sensitivity, Specificity, Precision, FScore, Gain;

	GetMeasures(num_correct_stayed_correct, num_wrong_correctly_fixed, num_correct_wrongly_fixed, num_wrong_stayed_wrong,
					Recall_Sensitivity, Specificity, Precision, FScore, Gain);

	double Recall_Sensitivity_bases, Specificity_bases, Precision_bases, FScore_bases, Gain_bases;

	GetMeasures(num_correct_stayed_correct_bases, num_wrong_correctly_fixed_bases, num_correct_wrongly_fixed_bases, num_wrong_stayed_wrong_bases,
					Recall_Sensitivity_bases, Specificity_bases, Precision_bases, FScore_bases, Gain_bases);

	if(allow_trimming)
	{
		fprintf(file_eval, "-------------------------------------------------------------------------------\n");
		fprintf(file_eval, "*) num_trimmed_bases                = %lld [of %lld] (%lf %%)\n", num_trimmed_bases, total_num_bases, 100.0*num_trimmed_bases/total_num_bases);
		fprintf(file_eval, "*) num_trimmed_bases_from_corrected = %lld [of %lld] (%lf %%)\n", num_trimmed_bases_from_corrected_reads, total_num_bases_from_corrected_reads, 100.0*num_trimmed_bases_from_corrected_reads/total_num_bases_from_corrected_reads);
	}

	fprintf(file_eval, "-------------------------------------------------------------------------------\n");
	fprintf(file_eval, "Whole Read Statistics:\n");
	fprintf(file_eval, "TN) num_correct_stayed_correct      = %lld [of %lld] (%lf %%)\n", num_correct_stayed_correct, total_num_reads, 100.0*num_correct_stayed_correct/total_num_reads);
	fprintf(file_eval, "TP) num_wrong_correctly_fixed       = %lld [of %lld] (%lf %%)\n", num_wrong_correctly_fixed, total_num_reads, 100.0*num_wrong_correctly_fixed/total_num_reads);
	fprintf(file_eval, "FP) num_correct_wrongly_fixed       = %lld [of %lld] (%lf %%)\n", num_correct_wrongly_fixed, total_num_reads, 100.0*num_correct_wrongly_fixed/total_num_reads);
	fprintf(file_eval, "FN) num_wrong_stayed_wrong          = %lld [of %lld] (%lf %%)\n", num_wrong_stayed_wrong, total_num_reads, 100.0*num_wrong_stayed_wrong/total_num_reads);
	fprintf(file_eval, "*) Recall (Sensitivity)             = %lf %%\n", 100*Recall_Sensitivity);
	fprintf(file_eval, "*) Precision                        = %lf %%\n", 100*Precision);
	fprintf(file_eval, "*) FScore                           = %lf %%\n", 100*FScore);
	fprintf(file_eval, "*) Gain                             = %lf %%\n", 100*Gain);
	fprintf(file_eval, "-------------------------------------------------------------------------------\n");

	if(allow_trimming)
	{
		printf("-------------------------------------------------------------------------------\n");
		printf("*) num_trimmed_bases                = %lld [of %lld] (%lf %%)\n", num_trimmed_bases, total_num_bases, 100.0*num_trimmed_bases/total_num_bases);
		printf("*) num_trimmed_bases_from_corrected = %lld [of %lld] (%lf %%)\n", num_trimmed_bases_from_corrected_reads, total_num_bases_from_corrected_reads, 100.0*num_trimmed_bases_from_corrected_reads/total_num_bases_from_corrected_reads);
	}

	printf("-------------------------------------------------------------------------------\n");
	printf("Whole Read Statistics:\n");
	printf("TN) num_correct_stayed_correct      = %lld [of %lld] (%lf %%)\n", num_correct_stayed_correct, total_num_reads, 100.0*num_correct_stayed_correct/total_num_reads);
	printf("TP) num_wrong_correctly_fixed       = %lld [of %lld] (%lf %%)\n", num_wrong_correctly_fixed, total_num_reads, 100.0*num_wrong_correctly_fixed/total_num_reads);
	printf("FP) num_correct_wrongly_fixed       = %lld [of %lld] (%lf %%)\n", num_correct_wrongly_fixed, total_num_reads, 100.0*num_correct_wrongly_fixed/total_num_reads);
	printf("FN) num_wrong_stayed_wrong          = %lld [of %lld] (%lf %%)\n", num_wrong_stayed_wrong, total_num_reads, 100.0*num_wrong_stayed_wrong/total_num_reads);
	printf("*) Recall (Sensitivity)             = %lf %%\n", 100*Recall_Sensitivity);
	printf("*) Precision                        = %lf %%\n", 100*Precision);
	printf("*) FScore                           = %lf %%\n", 100*FScore);
	printf("*) Gain                             = %lf %%\n", 100*Gain);
	printf("-------------------------------------------------------------------------------\n");

	fprintf(file_eval, "Bases Statistics:\n");
	fprintf(file_eval, "TN) num_correct_stayed_correct      = %lld [of %lld] (%lf %%)\n", num_correct_stayed_correct_bases, total_num_bases-num_trimmed_bases, 100.0*num_correct_stayed_correct_bases/total_num_bases);
	fprintf(file_eval, "TP) num_wrong_correctly_fixed       = %lld [of %lld] (%lf %%)\n", num_wrong_correctly_fixed_bases, total_num_bases-num_trimmed_bases, 100.0*num_wrong_correctly_fixed_bases/total_num_bases);
	fprintf(file_eval, "FP) num_correct_wrongly_fixed       = %lld [of %lld] (%lf %%)\n", num_correct_wrongly_fixed_bases, total_num_bases-num_trimmed_bases, 100.0*num_correct_wrongly_fixed_bases/total_num_bases);
	fprintf(file_eval, "FN) num_wrong_stayed_wrong          = %lld [of %lld] (%lf %%)\n", num_wrong_stayed_wrong_bases, total_num_bases-num_trimmed_bases, 100.0*num_wrong_stayed_wrong_bases/total_num_bases);
	fprintf(file_eval, "*) Recall (Sensitivity)             = %lf %%\n", 100*Recall_Sensitivity_bases);
	fprintf(file_eval, "*) Precision                        = %lf %%\n", 100*Precision_bases);
	fprintf(file_eval, "*) FScore                           = %lf %%\n", 100*FScore_bases);
	fprintf(file_eval, "*) Gain                             = %lf %%\n", 100*Gain_bases);
	fprintf(file_eval, "-------------------------------------------------------------------------------\n");

	printf("Bases Statistics:\n");
	printf("TN) num_correct_stayed_correct      = %lld [of %lld] (%lf %%)\n", num_correct_stayed_correct_bases, total_num_bases-num_trimmed_bases, 100.0*num_correct_stayed_correct_bases/total_num_bases);
	printf("TP) num_wrong_correctly_fixed       = %lld [of %lld] (%lf %%)\n", num_wrong_correctly_fixed_bases, total_num_bases-num_trimmed_bases, 100.0*num_wrong_correctly_fixed_bases/total_num_bases);
	printf("FP) num_correct_wrongly_fixed       = %lld [of %lld] (%lf %%)\n", num_correct_wrongly_fixed_bases, total_num_bases-num_trimmed_bases, 100.0*num_correct_wrongly_fixed_bases/total_num_bases);
	printf("FN) num_wrong_stayed_wrong          = %lld [of %lld] (%lf %%)\n", num_wrong_stayed_wrong_bases, total_num_bases-num_trimmed_bases, 100.0*num_wrong_stayed_wrong_bases/total_num_bases);
	printf("*) Recall (Sensitivity)             = %lf %%\n", 100*Recall_Sensitivity_bases);
	printf("*) Precision                        = %lf %%\n", 100*Precision_bases);
	printf("*) FScore                           = %lf %%\n", 100*FScore_bases);
	printf("*) Gain                             = %lf %%\n", 100*Gain_bases);
	printf("-------------------------------------------------------------------------------\n");

	fclose(file_align);
	fclose(file_eval);

	FreeHeap(str_ref, 2*(file_ref_size+circular+1)+2+1);

	for(i=0;i<max_reads_per_step;i++)
	{
		pthread_mutex_destroy(&pr[i].done_mutex);
	}

	for(i=0;i<num_threads;i++)
	{
		delete psols[i];
	}

	delete[] psols;
	delete[] thread_ID;
	delete[] exit_status;

	delete[] pr;

	delete[] buf;
	delete[] r1;
	delete[] r2;
	delete[] aligned;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EvaluateMain(int argc, char* argv[])
{
	printf("Evaluation Started\n"); fflush(NULL);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	int			default_match_type=MATCH_EDIT_DIST;
	int			default_allow_trimming=0;
	int			default_circular=0;

	int			default_num_threads=16;
	int			default_max_reads_per_step=1000;
	int			default_file_block_size=1024*1024;

	const char* default_input_dir="./";
	const char* default_result_dir="./";

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	int			match_type=default_match_type; bool set_match_type; set_match_type=false;
	int			allow_trimming=default_allow_trimming; //bool_set_allow_trimming; set_allow_trimming=false;
	int			circular=default_circular; //bool_set_circular; set_circular=false;

	int			num_threads=default_num_threads; //bool_set_num_threads; set_num_threads=false;
	int			max_reads_per_step=default_max_reads_per_step; //bool_set_max_reads_per_step; set_max_reads_per_step=false;
	int			file_block_size=default_file_block_size; //bool_set_file_block_size; set_file_block_size=false;

	char 		input_dir[MAX_PATH_LEN+1]; strcpy(input_dir, default_input_dir); //bool_set_input_dir; set_input_dir=false;
	char 		result_dir[MAX_PATH_LEN+1]; strcpy(result_dir, default_result_dir); //bool_set_result_dir; set_result_dir=false;

	vector<FileName> user_input_file_names; bool set_user_input_file_names; set_user_input_file_names=false;
	vector<FileName> user_result_file_names; bool set_user_result_file_names; set_user_result_file_names=false;
	char		user_ref_genome_file_name[MAX_PATH_LEN+1]; user_ref_genome_file_name[0]=0; bool set_user_ref_genome_file_name; set_user_ref_genome_file_name=false;
	char		user_align_file_name[MAX_PATH_LEN+1]; user_align_file_name[0]=0; bool set_user_align_file_name; set_user_align_file_name=false;
	char		user_eval_file_name[MAX_PATH_LEN+1]; user_eval_file_name[0]=0; bool set_user_eval_file_name; set_user_eval_file_name=false;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*******************************************************************************************************
	// GetData from user here, be sure that slash is added to directory names

	int i;

	for(i=2;i<argc;i++)
	{
		bool processed=false;
		char* str=argv[i];

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
				else if(strcmp(opt,"-matchtype")==0)
				{
					if(strcmp(val,"edit")==0) {match_type=MATCH_EDIT_DIST; processed=true; set_match_type=true;}
					else if(strcmp(val,"hamming")==0) {match_type=MATCH_SUBSTITUTE_ONLY; processed=true; set_match_type=true;}
					else if(strcmp(val,"insdel")==0) {match_type=MATCH_INSERT_DELETE_ONLY; processed=true; set_match_type=true;}
				}
				else if(strcmp(opt,"-circular")==0)
				{
					sscanf(val, "%d", &circular);
					processed=true; //set_circular=true;
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
				else if(strcmp(opt,"-inputdir")==0)
				{
					strcpy(input_dir, val); AddSlash(input_dir);
					processed=true; //set_input_dir=true;
				}
				else if(strcmp(opt,"-resultdir")==0)
				{
					strcpy(result_dir, val); AddSlash(result_dir);
					processed=true; //set_result_dir=true;
				}
				else if(strcmp(opt,"-inputfile")==0) // REQUIRED
				{
					FileName f; f.s[0]=0;
					strcpy(f.s, val);
					user_input_file_names.push_back(f);
					processed=true; set_user_input_file_names=true;
				}
				else if(strcmp(opt,"-resultfile")==0) // REQUIRED
				{
					FileName f; f.s[0]=0;
					strcpy(f.s, val);
					user_result_file_names.push_back(f);
					processed=true; set_user_result_file_names=true;
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
				else if(strcmp(opt,"-evalfile")==0) // REQUIRED
				{
					strcpy(user_eval_file_name, val);
					processed=true; set_user_eval_file_name=true;
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

	if(num_threads<1) num_threads=1;
	if(max_reads_per_step<100) max_reads_per_step=100;
	if(circular<0) circular=0;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	FastqFiles fastq_file_names;
	FastqFiles res_file_names;

	char ref_genome_file_name[MAX_PATH_LEN+1];

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	FileName file;

	for(i=0;i<(int)user_input_file_names.size();i++)
	{
		GetFileName(user_input_file_names[i].s, input_dir, file.s);
		fastq_file_names.AddFile(file);
	}

	for(i=0;i<(int)user_result_file_names.size();i++)
	{
		GetFileName(user_result_file_names[i].s, result_dir, file.s);
		res_file_names.AddFile(file);
	}

	GetFileName(user_ref_genome_file_name, input_dir, ref_genome_file_name);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

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
	else if(!set_user_result_file_names)
	{
		printf("Please specify fasta or fastq corrected results file(s) using \"-resultfile=filename\"\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_user_ref_genome_file_name)
	{
		printf("Please specify input reference genome file using \"-refgenomefile=filename\"\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_user_align_file_name)
	{
		printf("Please specify input alignment file using \"-alignfile=filename\"\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_user_eval_file_name)
	{
		printf("Please specify output evaluation file using \"-evalfile=filename\"\n"); fflush(NULL);
		accepted=false;
	}
	else if(fastq_file_names.file_names.size()!=res_file_names.file_names.size())
	{
		printf("The number of input file(s) must match the number of corrected results file(s).\n"); fflush(NULL);
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

	for(i=0;accepted && i<(int)res_file_names.file_names.size();i++)
	{
		FILE* file=fopen(res_file_names.file_names[i].s, "r");
		if(!file)
		{
			printf("Cannot open corrected results file [%s] for reading.\n", res_file_names.file_names[i].s); fflush(NULL);
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
		FILE* file=fopen(user_align_file_name, "r");
		if(!file)
		{
			printf("Cannot open alignment file [%s] for reading.\n", user_align_file_name); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(accepted)
	{
		FILE* file=fopen(user_eval_file_name, "w");
		if(!file)
		{
			printf("Cannot open evaluation file [%s] for writing.\n", user_eval_file_name); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(!accepted)
		return;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	Evaluate<SuffixCactus>(fastq_file_names, res_file_names, ref_genome_file_name, user_align_file_name, user_eval_file_name, match_type, allow_trimming, circular, num_threads, max_reads_per_step, file_block_size);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Evaluation Finished\n\n"); fflush(NULL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EvalHelp()
{
	printf("-Run \"karect -eval [options_list]\" for the evaluation tool.\n");
	printf("-Available options: (i=integer, f=file, s=directory)\n");
	printf("-Essential options:\n");
	printf("  \"-matchtype=[edit|hamming|insdel]\": Specify the matching type. \"hamming\" allows substitution errors only.\n"
		   "                        \"edit\" allows insertions, deletion, and substitutions with equal costs.\n"
		   "                        \"insdel\" is the same as \"edit\", but the cost of substitutions is doubled.\n");
	printf("  \"-inputfile=f\":     Specify an input fasta/fastq file. This option can be repeated for multiple files.\n");
	printf("  \"-resultfile=f\":    Specify a result fasta/fastq file (resulting from running \"karect -correct\",\n"
		   "                        or any other correction tool). This option can be repeated for multiple files.\n");
	printf("  \"-refgenomefile=f\": Specify the file containing the reference genome (to be aligned with).\n");
	printf("  \"-alignfile=f\":     Specify the alignment file resulted from running \"karect -align\".\n");
	printf("  \"-evalfile=f\":      Specify the output evaluation file.\n");
	printf("-Basic options:\n");
	printf("  \"-inputdir=s\":      Specify the input files directory. Ignored if file paths are complete [Default=.].\n");
	printf("  \"-resultdir=s\":     Specify the result files directory. Ignored if file paths are complete [Default=.].\n");
	printf("  \"-threads=i\":       Specify the number of threads [Default=16].\n");
	printf("-Advanced options:\n");
	printf("  \"-circular=i\":      Specify the sequence size to be appended circularly (for circular genomes) [Default=0].\n");
	printf("  \"-readsperstep=i\":  Specify the maximum number of processed reads per step [Default=1000].\n");
	printf("-Experimental options:\n");
	printf("  \"-trim=[yes|no]\":   Allow/Disallow trimming. Allow trimming if you used \"-trim=yes\" while correction,\n"
		   "                        however, in this case the evaluation results are experimental [Default=no].\n");
	printf("-Example:\n"
		   "        ./karect -eval -inputdir=/sra_data -inputfile=SRR001666_1.fasta -inputfile=SRR001666_2.fasta\n"
		   "              -resultdir=/results -resultfile=karect_SRR001666_1.fasta -resultfile=karect_SRR001666_2.fasta\n"
		   "              -refgenomefile=NC_000913.fna -alignfile=./SRR001666_align.txt -evalfile=./SRR001666_eval.txt\n"
		   "              -matchtype=hamming -threads=12\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
