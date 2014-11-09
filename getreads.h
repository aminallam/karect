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

#define MAX_READ_LEN	(1<<16)
#define MAX_PATH_LEN	(1<<14)

#define NUM_VALID_CHARS 4
#define BITS_PER_SYMBOL 2

const char QUAL_STR[]="!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
const int QUAL_STR_SIZE=sizeof(QUAL_STR)/sizeof(char);

char			Comp[1<<8];
unsigned int	CharToInt[1<<8];
unsigned int	CharToInt5Val[1<<8];
char			IntToChar[5];
int				CharToQual[1<<8];
double			CharToProb[1<<8];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

char QualToChar(int val)
{
	if(val<0) val=0;
	if(val>=QUAL_STR_SIZE) val=QUAL_STR_SIZE-1;
	return QUAL_STR[val];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ProbToQual(double val)
{
	if(val<=0.2) return 1; // do not use the prob 0 quality score (!)
	if(val>=0.9999) return 40;
	int i=(-10.0*log(1-val)/log(10.0))+0.5;
	return i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AddSlash(char* dir)
{
	int len=strlen(dir);
	if(dir[len-1]!='/' && dir[len-1]!='\\')
	{
		dir[len]='/';
		dir[len+1]=0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct FileName
{
	char s[MAX_PATH_LEN+1];

	FileName(){}

	FileName(const char* p)
	{
		strcpy(s,p);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetFileName(const char* user_file, const char* input_dir, char* file_name)
{
	int user_file_len=strlen(user_file);
	int complete_path=(user_file_len>0 && (user_file[0]=='\\' || user_file[0]=='/')) || (user_file_len>1 && user_file[1]==':');

	file_name[0]=0;
	if(!complete_path) strcpy(file_name, input_dir);
	int len_input_dir=strlen(file_name);
	strcpy(file_name+len_input_dir, user_file);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct FastqFiles
{
	vector<FileName>	file_names;
	int					cur_file;

	void Reset()
	{
		cur_file=0;
	}

	const char* GetCurFileName()
	{
		if(cur_file>=(int)file_names.size()) return 0;
		return file_names[cur_file].s;
	}

	void FinishCurFile()
	{
		cur_file++;
	}

	void AddFile(FileName& file_name)
	{
		file_names.push_back(file_name);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GlobalInitComp()
{
	int i;
	for(i=0;i<(1<<8);i++) Comp[i]=i;
	Comp['A']='T';
	Comp['T']='A';
	Comp['C']='G';
	Comp['G']='C';

	for(i=0;i<(1<<8);i++) CharToInt[i]=0;
	CharToInt['A']=0;
	CharToInt['T']=1;
	CharToInt['C']=2;
	CharToInt['G']=3;

	for(i=0;i<(1<<8);i++) CharToInt5Val[i]=4;
	CharToInt5Val['A']=0;
	CharToInt5Val['T']=1;
	CharToInt5Val['C']=2;
	CharToInt5Val['G']=3;

	IntToChar[0]='A';
	IntToChar[1]='T';
	IntToChar[2]='C';
	IntToChar[3]='G';
	IntToChar[4]='N';

	for(i=0;i<(1<<8);i++) CharToQual[i]=0;
	for(i=0;i<QUAL_STR_SIZE;i++) CharToQual[(unsigned char)QUAL_STR[i]]=i;

	for(i=0;i<(1<<8);i++) CharToProb[i]=0.1;
	for(i=0;i<QUAL_STR_SIZE;i++) CharToProb[(unsigned char)QUAL_STR[i]]=1-pow(10.0, -i/10.0);
	CharToProb[(unsigned char)QUAL_STR[0]]=0.1; // (it should be 0, but we increase it to 0.1)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
void Reverse(char* str, Type nstr)
{
	Type j; for(j=0;j<nstr/2;j++) {char t=str[j];str[j]=str[nstr-1-j];str[nstr-1-j]=t;}
}

template<class Type>
void ReverseComplement(char* str, Type nstr)
{
	Type j; for(j=0;j<nstr/2;j++) {char t=Comp[(unsigned char)str[j]];str[j]=Comp[(unsigned char)str[nstr-1-j]];str[nstr-1-j]=t;}
	if(nstr%2) str[nstr/2]=Comp[(unsigned char)str[nstr/2]];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

char CorrectReadChar(char ch, char unknown_char)
{
	if(ch=='a') return 'A';
	if(ch=='c') return 'C';
	if(ch=='g') return 'G';
	if(ch=='t') return 'T';
	if(ch=='A' || ch=='C' || ch=='G' || ch=='T') return ch;
	return unknown_char;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

char GetFirstChar(const char* file_name, FILE* file_fastq, char* buf, char* info_line)
{
	if(file_fastq==0)
	{
		printf("Cant open file: [%s]\n", file_name);
		fflush(NULL);
		return 0;
	}

	buf[0]=0;

	// detect file type from first character (file can be fasta or fastq)

	if(EOF==fscanf(file_fastq, "%[^\n\r] ", buf))
	{
		printf("Cannot read first line of the file: [%s]\n", file_name);
		fflush(NULL);
		return 0;
	}

	if(buf[0]!='>' && buf[0]!='@')
	{
		printf("Unsupported file format. File starts with unrecognized symbol (%c). Only fasta and fastq formats are supported. [%s]\n", buf[0], file_name);
		fflush(NULL);
		return 0;
	}

	if(info_line) strcpy(info_line, buf);
	return buf[0];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int GetReadSequence(FILE* file_fastq, char first_char, char* seq, char* buf, char* info_line, char* quality_info, char* quality_scores, char unknown_char, bool replace_eol_by_space=false)
{
	int seq_len=0;
	seq[0]=0;
	buf[0]=0;
	if(quality_info) quality_info[0]=0;
	if(quality_scores) quality_scores[0]=0;

	bool file_end=false;

	while(1)
	{
		if(EOF==fscanf(file_fastq, "%[^\n\r] ", seq+seq_len)) {file_end=true; break;}

		if(first_char=='>' && seq[seq_len]=='>')
		{
			if(info_line) strcpy(info_line, seq+seq_len);
			seq[seq_len]=0;
			break;
		}

		if(first_char=='@' && seq[seq_len]=='+')
		{
			if(quality_info) strcpy(quality_info, seq+seq_len);
			int buf_len=0;

			seq[seq_len]=0;
			while(1)
			{
				if(EOF==fscanf(file_fastq, "%[^\n\r] ", buf+buf_len)) {file_end=true; break;}
				if(buf[buf_len]=='@' && buf_len>=seq_len) break; // note that the @ can appear at the beginning of the quality scores line!
				buf_len+=strlen(buf+buf_len);
			}
			if(info_line) strcpy(info_line, buf+buf_len);
			buf[buf_len]=0;
			if(quality_scores) strcpy(quality_scores, buf);
			break;
		}
		seq_len+=strlen(seq+seq_len);
		if(replace_eol_by_space) {seq[seq_len++]=' '; seq[seq_len]=0;}
	}

	if(unknown_char)
	{
		int i;
		for(i=0;i<seq_len;i++)
		{
			seq[i]=CorrectReadChar(seq[i], unknown_char);
		}
	}

	if(!file_end && seq_len==0)
	{
		seq[0]='A'; seq[1]=0;
		if(quality_scores) {quality_scores[0]='I'; quality_scores[1]=0;}
		seq_len=1;
	}
	return seq_len;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int GetGlobalReadStats(FastqFiles& file_name, const char* input_file_name, int file_block_size, int insert_complement, int kmer_size,
						long long& num_reads, long long& num_kmers, long long& total_size, int& max_read_len, bool& is_fastq) // size does not include separators
{
	num_kmers=0;
	max_read_len=0;
	num_reads=0;
	total_size=1;
	is_fastq=false;

	file_name.Reset();

	FILE* file_fastq=0;
	if(!file_name.GetCurFileName()) return 0;
	file_fastq=fopen(file_name.GetCurFileName(), "r");

	FILE* file_input=fopen(input_file_name, "wb");

	BufferedOutFile buf_input;
	buf_input.Initialize(file_input, file_block_size);

	char* buf=new char[MAX_READ_LEN+1];
	char* seq=new char[MAX_READ_LEN+1];
	char* quality=new char[MAX_READ_LEN+1];

	char first_char=GetFirstChar(file_name.GetCurFileName(), file_fastq, buf, 0);
	if(first_char==0) return 0;
	if(first_char=='@') is_fastq=true;

	while(1)
	{
		int seq_len=GetReadSequence(file_fastq, first_char, seq, buf, 0, 0, quality, 'N');
		if(seq_len==0)
		{
			fclose(file_fastq);
			file_name.FinishCurFile();

			if(file_name.GetCurFileName())
			{
				file_fastq=fopen(file_name.GetCurFileName(), "r");
				first_char=GetFirstChar(file_name.GetCurFileName(), file_fastq, buf, 0);
				if(first_char==0) return 0;
				seq_len=GetReadSequence(file_fastq, first_char, seq, buf, 0, 0, quality, 'N');
				if(seq_len==0) return 0;
			}
		}
		if(seq_len==0) break;

		total_size+=(seq_len+1)*(insert_complement+1);
		num_reads+=(insert_complement+1);
		if(seq_len-kmer_size+1>0) num_kmers+=(seq_len-kmer_size+1)*(insert_complement+1);

		if(seq_len>max_read_len) max_read_len=seq_len;

		seq[seq_len]='$';
		buf_input.WriteNumBytes(seq, seq_len+1);

		if(is_fastq && (quality[0]==0 || (int)strlen(quality)!=seq_len))
		{
			printf("Incorrect Fastq Format: Quality string length does not match read length.\n");
			seq[seq_len]=0;
			printf("Read=[%s]\n", seq);
			printf("Qual=[%s]\n", quality);
			fflush(NULL);
			return 0;
		}

		if(is_fastq) buf_input.WriteNumBytes(quality, seq_len);
	}

	buf_input.Destroy();

	fclose(file_input);

	delete[] seq;
	delete[] buf;
	delete[] quality;

	if(num_reads==0) return 0;
	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetReads(const char* file_name, int file_block_size, long long& start_file, int insert_complement, int kmer_size, long long max_num_kmers,
				char** pall_reads, char** pall_quals, long long** pall_reads_starts, long long& num_reads, long long& num_kmers, long long& total_size, int& max_read_len, bool is_fastq, int use_qual)
{
	*pall_quals=0;

	num_kmers=0;
	max_read_len=0;
	num_reads=0;
	total_size=1;

	FILE* file_reads=fopen(file_name, "rb");

	fseek_long(file_reads, start_file);

	BufferedInFile buf_reads;
	buf_reads.Initialize(file_reads, file_block_size);

	char* seq=new char[MAX_READ_LEN+1];

	while(num_kmers<max_num_kmers)
	{
		int seq_len=buf_reads.ReadUntil(seq, '$');
		if(seq_len==0) break;

		total_size+=(seq_len+1)*(insert_complement+1);
		num_reads+=(insert_complement+1);
		if(seq_len-kmer_size+1>0) num_kmers+=(seq_len-kmer_size+1)*(insert_complement+1);

		if(seq_len>max_read_len) max_read_len=seq_len;

		if(is_fastq) buf_reads.ReadNumBytes(seq, seq_len);
	}

	buf_reads.Destroy();

	fseek_long(file_reads, start_file);
	buf_reads.Initialize(file_reads, file_block_size);

	char* all_reads=(char*)AllocateHeap(total_size+2); // 1 for the # and 1 for the \0
	*pall_reads=all_reads;

	if(is_fastq && use_qual) *pall_quals=(char*)AllocateHeap(total_size+2);
	char* all_quals=*pall_quals;

	long long* all_reads_starts=(long long*)AllocateHeap(sizeof(long long)*(num_reads+1));
	*pall_reads_starts=all_reads_starts;

	long long cur_read=0;
	long long cur_str_size=0;

	*(all_reads+cur_str_size)='$';
	if(all_quals) *(all_quals+cur_str_size)='!';
	cur_str_size++;

	while(cur_read<num_reads)
	{
		int read_len=buf_reads.ReadUntil(seq, '$');
		if(read_len==0) break;

		all_reads_starts[cur_read++]=cur_str_size;
		memcpy(all_reads+cur_str_size, seq, read_len);
		cur_str_size+=read_len;
		*(all_reads+cur_str_size)='$';
		cur_str_size++;

		if(insert_complement)
		{
			ReverseComplement(seq, read_len);

			all_reads_starts[cur_read++]=cur_str_size;
			memcpy(all_reads+cur_str_size, seq, read_len);
			cur_str_size+=read_len;
			*(all_reads+cur_str_size)='$';
			cur_str_size++;
		}

		if(is_fastq)
		{
			buf_reads.ReadNumBytes(seq, read_len);

			if(all_quals)
			{
				memcpy(all_quals+cur_str_size-(insert_complement+1)*(read_len+1), seq, read_len);
				*(all_quals+cur_str_size-insert_complement*(read_len+1)-1)='!';

				if(insert_complement)
				{
					Reverse(seq, read_len);
					memcpy(all_quals+cur_str_size-(read_len+1), seq, read_len);
					*(all_quals+cur_str_size-1)='!';
				}
			}
		}
	}

	all_reads_starts[cur_read++]=cur_str_size; // used as sentinel to get the size of the end read fastly

	*(all_reads+cur_str_size)='#';
	if(all_quals) *(all_quals+cur_str_size)='!';
	cur_str_size++;

	*(all_reads+cur_str_size)=0;
	if(all_quals) *(all_quals+cur_str_size)=0;

	// exclude the first $ included in the total size
	start_file+=(total_size-1)/(insert_complement+1);
	if(is_fastq) start_file+=(total_size-1-num_reads)/(insert_complement+1);

	buf_reads.Destroy();

	fclose(file_reads);

	delete[] seq;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int GetKmer(const char* str, int len, int kmer_size)
{
	int j;
	unsigned int cur_kmer=0;

	for(j=0;j<kmer_size && j<len;j++)
	{
		cur_kmer=(CharToInt[(int)str[j]]<<(j*BITS_PER_SYMBOL))|cur_kmer;
	}

	return cur_kmer;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetStrFromKmer(unsigned int kmer, char* str, int kmer_size)
{
	int j;
	unsigned int mask=(1<<BITS_PER_SYMBOL)-1;

	for(j=0;j<kmer_size;j++)
	{
		str[j]=IntToChar[kmer&mask];
		kmer>>=BITS_PER_SYMBOL;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int GetRevKmer(const char* str, int len, int kmer_size)
{
	int j;
	unsigned int cur_kmer=0;

	for(j=0;j<kmer_size && j<len;j++)
	{
		cur_kmer=(CharToInt[(int)str[-j]]<<(j*BITS_PER_SYMBOL))|cur_kmer;
	}

	return cur_kmer;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct KmerCount
{
	unsigned int kmer;
	unsigned int cnt;
	bool operator<(const KmerCount& kc)const {return cnt<kc.cnt;}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double EstimateCoverage(const char* file_name, int file_block_size, int insert_complement, bool is_fastq)
{
	long long i;
	int j;

	const int num_mod_kmers=1000;
	const int basic_kmer_size=10;
	const int extend_kmer_size=7;

	long long num_basic_kmers=((long long)1)<<(basic_kmer_size*BITS_PER_SYMBOL);

	vector<KmerCount> basic_kmers(num_basic_kmers);
	for(i=0;i<num_basic_kmers;i++) {basic_kmers[i].kmer=i; basic_kmers[i].cnt=0;}

	FILE* file_reads=fopen(file_name, "rb");

	BufferedInFile buf_reads;
	buf_reads.Initialize(file_reads, file_block_size);

	char* read=new char[MAX_READ_LEN+1];

	while(1)
	{
		int len=buf_reads.ReadUntil(read, '$');
		if(len==0) break;

		if(len>=basic_kmer_size)
		{
			unsigned int kmer=GetKmer(read, basic_kmer_size, basic_kmer_size);

			for(j=0;j<len-basic_kmer_size+1;j++)
			{
				if(j){
					int new_char=CharToInt[(int)read[j+basic_kmer_size-1]];
					kmer=(kmer>>BITS_PER_SYMBOL)|(new_char<<((basic_kmer_size-1)*BITS_PER_SYMBOL));
				}
				basic_kmers[kmer].cnt++;
			}

			ReverseComplement(read, len);

			kmer=GetKmer(read, basic_kmer_size, basic_kmer_size);

			for(j=0;j<len-basic_kmer_size+1;j++)
			{
				if(j){
					int new_char=CharToInt[(int)read[j+basic_kmer_size-1]];
					kmer=(kmer>>BITS_PER_SYMBOL)|(new_char<<((basic_kmer_size-1)*BITS_PER_SYMBOL));
				}
				basic_kmers[kmer].cnt++;
			}
		}

		if(is_fastq) buf_reads.ReadNumBytes(read, len);
	}

	buf_reads.Destroy();
	fclose(file_reads);

	sort(basic_kmers.begin(), basic_kmers.end());

	unsigned int mod_kmers[num_mod_kmers];

	for(i=0;i<num_mod_kmers;i++)
		mod_kmers[i]=basic_kmers[i+num_basic_kmers/2-num_mod_kmers/2].kmer;

	for(i=0;i<num_basic_kmers;i++) basic_kmers[i].cnt=0;
	for(i=0;i<num_mod_kmers;i++) basic_kmers[mod_kmers[i]].cnt=i+1;

	long long num_extend_kmers=((long long)1)<<(extend_kmer_size*BITS_PER_SYMBOL);

	unsigned int* all_kmer_counts=new unsigned int[num_mod_kmers*num_extend_kmers];
	memset(all_kmer_counts, 0, sizeof(unsigned int)*num_mod_kmers*num_extend_kmers);

	file_reads=fopen(file_name, "rb");

	buf_reads.Initialize(file_reads, file_block_size);

	while(1)
	{
		int len=buf_reads.ReadUntil(read, '$');
		if(len==0) break;

		if(len>=basic_kmer_size+extend_kmer_size)
		{
			unsigned int kmer=GetKmer(read, basic_kmer_size, basic_kmer_size);
			unsigned int extend_kmer=GetKmer(read+basic_kmer_size, extend_kmer_size, extend_kmer_size);

			for(j=0;j<len-basic_kmer_size-extend_kmer_size+1;j++)
			{
				if(j)
				{
					int new_char=CharToInt[(int)read[j+basic_kmer_size-1]];
					kmer=(kmer>>BITS_PER_SYMBOL)|(new_char<<((basic_kmer_size-1)*BITS_PER_SYMBOL));

					new_char=CharToInt[(int)read[j+basic_kmer_size+extend_kmer_size-1]];
					extend_kmer=(extend_kmer>>BITS_PER_SYMBOL)|(new_char<<((extend_kmer_size-1)*BITS_PER_SYMBOL));
				}
				if(basic_kmers[kmer].cnt)
					all_kmer_counts[(basic_kmers[kmer].cnt-1)*num_extend_kmers+extend_kmer]++;
			}

			ReverseComplement(read, len);

			kmer=GetKmer(read, basic_kmer_size, basic_kmer_size);
			extend_kmer=GetKmer(read+basic_kmer_size, extend_kmer_size, extend_kmer_size);

			for(j=0;j<len-basic_kmer_size-extend_kmer_size+1;j++)
			{
				if(j)
				{
					int new_char=CharToInt[(int)read[j+basic_kmer_size-1]];
					kmer=(kmer>>BITS_PER_SYMBOL)|(new_char<<((basic_kmer_size-1)*BITS_PER_SYMBOL));

					new_char=CharToInt[(int)read[j+basic_kmer_size+extend_kmer_size-1]];
					extend_kmer=(extend_kmer>>BITS_PER_SYMBOL)|(new_char<<((extend_kmer_size-1)*BITS_PER_SYMBOL));
				}
				if(basic_kmers[kmer].cnt)
					all_kmer_counts[(basic_kmers[kmer].cnt-1)*num_extend_kmers+extend_kmer]++;
			}
		}

		if(is_fastq) buf_reads.ReadNumBytes(read, len);
	}

	buf_reads.Destroy();
	fclose(file_reads);
	delete[] read;

	sort(all_kmer_counts, all_kmer_counts+num_mod_kmers*num_extend_kmers);

	double coverage=0;

	int st=0;
	while(st<num_mod_kmers*num_extend_kmers && all_kmer_counts[st]<=5) st++;

	if(st<num_mod_kmers*num_extend_kmers)
	{
		for(i=st;i<num_mod_kmers*num_extend_kmers;i++) coverage+=all_kmer_counts[i];
		coverage/=num_mod_kmers*num_extend_kmers-st;
	}

	delete[] all_kmer_counts;

	return coverage;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MergeTwoPairedEndFilesIntoOneFile(const char* file_name_a, const char* file_name_b, const char* file_name_out, int out_fasta, char unknown_char, int adjust_qual, int no_qual_info, int sample)
{
	FILE* file_fastq_a=fopen(file_name_a, "r");
	FILE* file_fastq_b=fopen(file_name_b, "r");
	FILE* file_fasta_out=fopen(file_name_out, "w");

	int cur_info=0;
	char* buf_a=new char[MAX_READ_LEN+1];
	char* buf_b=new char[MAX_READ_LEN+1];
	char* info_line_a[2]; info_line_a[0]=new char[MAX_READ_LEN+1]; info_line_a[1]=new char[MAX_READ_LEN+1];
	char* info_line_b[2]; info_line_b[0]=new char[MAX_READ_LEN+1]; info_line_b[1]=new char[MAX_READ_LEN+1];
	char* quality_info_a=new char[MAX_READ_LEN+1];
	char* quality_info_b=new char[MAX_READ_LEN+1];
	char* quality_scores_a=new char[MAX_READ_LEN+1];
	char* quality_scores_b=new char[MAX_READ_LEN+1];

	char first_char_a=GetFirstChar(file_name_a, file_fastq_a, buf_a, info_line_a[cur_info]);
	if(first_char_a==0) return 0;

	char first_char_b=GetFirstChar(file_name_b, file_fastq_b, buf_b, info_line_b[cur_info]);
	if(first_char_b==0) return 0;

	char* seq_a=new char[MAX_READ_LEN+1];
	char* seq_b=new char[MAX_READ_LEN+1];

	int finished_a=0;
	int finished_b=0;

	long long acc_sample=0;

	while(1)
	{
		int seq_len_a=0; if(!finished_a) seq_len_a=GetReadSequence(file_fastq_a, first_char_a, seq_a, buf_a, info_line_a[1-cur_info], quality_info_a, quality_scores_a, unknown_char);
		int seq_len_b=0; if(!finished_b) seq_len_b=GetReadSequence(file_fastq_b, first_char_b, seq_b, buf_b, info_line_b[1-cur_info], quality_info_b, quality_scores_b, unknown_char);
		if(seq_len_a==0 && seq_len_b==0) break;
		if(seq_len_a==0) finished_a=1; if(seq_len_b==0) finished_b=1;

		if(acc_sample/1000 < (acc_sample+sample)/1000)
		{
			if(out_fasta) info_line_a[cur_info][0]='>';
			if(seq_len_a>0) fprintf(file_fasta_out, "%s\n%s\n", info_line_a[cur_info], seq_a);
			if(!out_fasta && quality_info_a[0]){if(no_qual_info) quality_info_a[1]=0; fprintf(file_fasta_out, "%s\n", quality_info_a);}
			if(adjust_qual && !out_fasta && quality_scores_a[0]=='@') quality_scores_a[0]++;
			if(!out_fasta && quality_scores_a[0]) fprintf(file_fasta_out, "%s\n", quality_scores_a);

			if(out_fasta) info_line_b[cur_info][0]='>';
			if(seq_len_b>0) fprintf(file_fasta_out, "%s\n%s\n", info_line_b[cur_info], seq_b);
			if(!out_fasta && quality_info_b[0]) {if(no_qual_info) quality_info_b[1]=0; fprintf(file_fasta_out, "%s\n", quality_info_b);}
			if(adjust_qual && !out_fasta && quality_scores_b[0]=='@') quality_scores_b[0]++;
			if(!out_fasta && quality_scores_b[0]) fprintf(file_fasta_out, "%s\n", quality_scores_b);
		}

		acc_sample+=sample;

		cur_info=1-cur_info;
	}

	fclose(file_fastq_a);
	fclose(file_fastq_b);
	fclose(file_fasta_out);

	delete[] buf_a;
	delete[] buf_b;
	delete[] info_line_a[0]; delete[] info_line_a[1];
	delete[] info_line_b[0]; delete[] info_line_b[1];
	delete[] seq_a;
	delete[] seq_b;
	delete[] quality_info_a;
	delete[] quality_info_b;
	delete[] quality_scores_a;
	delete[] quality_scores_b;

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ConcatenateFilesIntoOneFile(vector<FileName>& file_names, const char* file_name_out, int out_fasta, int out_qual, int one_seq, char unknown_char, int adjust_qual, int no_qual_info, int out_fastq, int sample)
{
	int max_seq_len=MAX_READ_LEN;
	if(one_seq) max_seq_len=0x7FFFFFFE;

	char* buf=new char[MAX_READ_LEN+1];
	char* info_line[2]; info_line[0]=new char[MAX_READ_LEN+1]; info_line[1]=new char[MAX_READ_LEN+1];
	char* seq=new char[max_seq_len+1];
	char* quality_info=new char[MAX_READ_LEN+1];
	char* quality_scores=new char[MAX_READ_LEN+1];

	bool first_info=true;
	FILE* file_fasta_out=fopen(file_name_out, "w");

	FILE* file_qual_out=0;
	if(out_qual)
	{
		strcpy(buf, file_name_out);
		int len_buf=strlen(buf); while(len_buf>0 && buf[len_buf-1]!='.') len_buf--; len_buf--;
		if(len_buf<=0) len_buf=strlen(buf);
		strcpy(buf+len_buf, ".qual");
		file_qual_out=fopen(buf, "w");
	}

	long long acc_sample=0;

	int i;
	for(i=0;i<(int)file_names.size();i++)
	{
		FILE* file_fastq_in=fopen(file_names[i].s, "r");

		int cur_info=0;

		char first_char=GetFirstChar(file_names[i].s, file_fastq_in, buf, info_line[cur_info]);
		if(first_char==0) return 0;

		while(1)
		{
			int seq_len=GetReadSequence(file_fastq_in, first_char, seq, buf, info_line[1-cur_info], quality_info, quality_scores, unknown_char);
			if(seq_len==0) break;

			if(acc_sample/1000 < (acc_sample+sample)/1000)
			{
				if(out_fasta) info_line[cur_info][0]='>';

				if(out_fastq && (!quality_scores[0] || !quality_info[0]))
				{
					info_line[cur_info][0]='@';
					quality_info[0]='+'; quality_info[1]=0;
					int j; for(j=0;j<seq_len;j++) quality_scores[j]='I';
					quality_scores[seq_len]=0;
				}

				if(!one_seq) fprintf(file_fasta_out, "%s\n%s\n", info_line[cur_info], seq);
				else if(first_info) {fprintf(file_fasta_out, "%s\n%s", info_line[cur_info], seq); first_info=false;}
				else fprintf(file_fasta_out, "%s", seq);

				if(!out_fasta && quality_info[0]) {if(no_qual_info) quality_info[1]=0; fprintf(file_fasta_out, "%s\n", quality_info);}

				if(adjust_qual && !out_fasta && quality_scores[0]=='@') quality_scores[0]++;

				if(!out_fasta && quality_scores[0]) fprintf(file_fasta_out, "%s\n", quality_scores);

				if(out_qual && quality_scores[0])
				{
					if(quality_info[0]==0 || quality_info[1]==0) strcpy(quality_info, info_line[cur_info]);
					quality_info[0]='>';
					fprintf(file_qual_out, "%s\n", quality_info);
					int j;
					for(j=0;quality_scores[j];j++)
					{
						int q=CharToQual[(unsigned char)quality_scores[j]];
						if(j) fprintf(file_qual_out, " %d", q);
						else fprintf(file_qual_out, "%d", q);
					}
					fprintf(file_qual_out, "\n");
				}
			}

			acc_sample+=sample;

			cur_info=1-cur_info;
		}

		fclose(file_fastq_in);
	}

	if(one_seq) fprintf(file_fasta_out, "\n");
	fclose(file_fasta_out);

	if(file_qual_out) fclose(file_qual_out);

	delete[] buf;
	delete[] info_line[0]; delete[] info_line[1];
	delete[] seq;
	delete[] quality_info;
	delete[] quality_scores;

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int GetReadStats(vector<FileName>& file_names)
{
	long long total_num_reads=0;
	long long total_read_sizes=0;
	int max_read_len=0;
	int min_read_len=0x7FFFFFFF;

	int max_seq_len=MAX_READ_LEN;

	char* buf=new char[MAX_READ_LEN+1];
	char* info_line[2]; info_line[0]=new char[MAX_READ_LEN+1]; info_line[1]=new char[MAX_READ_LEN+1];
	char* seq=new char[max_seq_len+1];
	char* quality_info=new char[MAX_READ_LEN+1];
	char* quality_scores=new char[MAX_READ_LEN+1];

	int i;
	for(i=0;i<(int)file_names.size();i++)
	{
		FILE* file_fastq_in=fopen(file_names[i].s, "r");

		int cur_info=0;

		char first_char=GetFirstChar(file_names[i].s, file_fastq_in, buf, info_line[cur_info]);
		if(first_char==0) return 0;

		while(1)
		{
			int seq_len=GetReadSequence(file_fastq_in, first_char, seq, buf, info_line[1-cur_info], quality_info, quality_scores, 0);
			if(seq_len==0) break;

			total_num_reads++;
			total_read_sizes+=seq_len;

			if(seq_len>max_read_len) max_read_len=seq_len;
			if(seq_len<min_read_len) min_read_len=seq_len;

			cur_info=1-cur_info;
		}

		fclose(file_fastq_in);
	}

	delete[] buf;
	delete[] info_line[0]; delete[] info_line[1];
	delete[] seq;
	delete[] quality_info;
	delete[] quality_scores;

	printf("Num Reads        = %lld\n", total_num_reads);
	printf("Size Reads       = %lld\n", total_read_sizes);
	printf("Max Read Len     = %d\n", max_read_len);
	printf("Min Read Len     = %d\n", min_read_len);
	printf("Average Read Len = %Lf\n", (long double)total_read_sizes/total_num_reads);
	fflush(NULL);

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool IsInfoPaired(char* info1, char* info2)
{
	char* a=new char[MAX_READ_LEN+1];
	char* b=new char[MAX_READ_LEN+1];

	int i, k, na=strlen(info1), nb=strlen(info2);
	strcpy(a, info1);
	strcpy(b, info2);

	for(i=0;i<na;i++) if(a[i]==' ' || a[i]=='\t') break; na=i; a[na]=0;
	for(i=0;i<nb;i++) if(b[i]==' ' || b[i]=='\t') break; nb=i; b[nb]=0;

	const char* strip[]={"/1", "/2", "_left", "_right"};
	int num_strips=sizeof(strip)/sizeof(const char*);

	for(k=0;k<num_strips;k++)
	{
		int ns=strlen(strip[k]);
		if(strcmp(a+na-ns, strip[k])==0) {na-=ns; a[na]=0; break;}
	}
	for(k=0;k<num_strips;k++)
	{
		int ns=strlen(strip[k]);
		if(strcmp(b+nb-ns, strip[k])==0) {nb-=ns; b[nb]=0; break;}
	}

	bool paired=(strcmp(a,b)==0);

	delete[] a;
	delete[] b;

	return paired;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int SplitPaired(const char* file_name_in, int out_fasta, char unknown_char, int adjust_qual, int no_qual_info)
{
	char* buf=new char[MAX_READ_LEN+1];
	char* info_line[2]; info_line[0]=new char[MAX_READ_LEN+1]; info_line[1]=new char[MAX_READ_LEN+1];
	char* seq=new char[MAX_READ_LEN+1];
	char* quality_info=new char[MAX_READ_LEN+1];
	char* quality_scores=new char[MAX_READ_LEN+1];

	const char* fragment_str="_frag";
	const char* paired_str="_pair";

	strcpy(buf, file_name_in); int len_buf=strlen(buf); while(len_buf>0 && buf[len_buf-1]!='.') len_buf--; len_buf--;
	if(len_buf<=0) len_buf=strlen(buf); strcpy(buf+len_buf, fragment_str); strcpy(buf+len_buf+strlen(fragment_str), file_name_in+len_buf);
	FILE* file_frag_out=fopen(buf, "w");

	strcpy(buf, file_name_in); len_buf=strlen(buf); while(len_buf>0 && buf[len_buf-1]!='.') len_buf--; len_buf--;
	if(len_buf<=0) len_buf=strlen(buf); strcpy(buf+len_buf, paired_str); strcpy(buf+len_buf+strlen(paired_str), file_name_in+len_buf);
	FILE* file_paired_out=fopen(buf, "w");

	FILE* file_fastq_in=fopen(file_name_in, "r");

	int cur_info=0;

	char first_char=GetFirstChar(file_name_in, file_fastq_in, buf, info_line[cur_info]);
	if(first_char==0) return 0;

	bool last_paired=false;

	while(1)
	{
		int seq_len=GetReadSequence(file_fastq_in, first_char, seq, buf, info_line[1-cur_info], quality_info, quality_scores, unknown_char);
		if(seq_len==0) break;

		FILE* cur_file=file_frag_out;
		bool cur_paired=false;
		if(last_paired || IsInfoPaired(info_line[cur_info], info_line[1-cur_info])) {cur_file=file_paired_out; cur_paired=true;}
		if(!last_paired && cur_paired) last_paired=true; else last_paired=false;

		if(out_fasta) info_line[cur_info][0]='>';
		if(seq_len>0) fprintf(cur_file, "%s\n%s\n", info_line[cur_info], seq);

		if(!out_fasta && quality_info[0]) {if(no_qual_info) quality_info[1]=0; fprintf(cur_file, "%s\n", quality_info);}

		if(adjust_qual && !out_fasta && quality_scores[0]=='@') quality_scores[0]++;

		if(!out_fasta && quality_scores[0]) fprintf(cur_file, "%s\n", quality_scores);

		cur_info=1-cur_info;
	}

	fclose(file_fastq_in);

	fclose(file_frag_out);
	fclose(file_paired_out);

	delete[] buf;
	delete[] info_line[0]; delete[] info_line[1];
	delete[] seq;
	delete[] quality_info;
	delete[] quality_scores;

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int SplitPairsTwoFiles(const char* file_name_in, int out_fasta, char unknown_char, int adjust_qual, int no_qual_info)
{
	char* buf=new char[MAX_READ_LEN+1];
	char* info_line[2]; info_line[0]=new char[MAX_READ_LEN+1]; info_line[1]=new char[MAX_READ_LEN+1];
	char* seq=new char[MAX_READ_LEN+1];
	char* quality_info=new char[MAX_READ_LEN+1];
	char* quality_scores=new char[MAX_READ_LEN+1];

	const char* pair_1_str="_1";
	const char* pair_2_str="_2";

	strcpy(buf, file_name_in); int len_buf=strlen(buf); while(len_buf>0 && buf[len_buf-1]!='.') len_buf--; len_buf--;
	if(len_buf<=0) len_buf=strlen(buf); strcpy(buf+len_buf, pair_1_str); strcpy(buf+len_buf+strlen(pair_1_str), file_name_in+len_buf);
	FILE* file_pair_1_out=fopen(buf, "w");

	strcpy(buf, file_name_in); len_buf=strlen(buf); while(len_buf>0 && buf[len_buf-1]!='.') len_buf--; len_buf--;
	if(len_buf<=0) len_buf=strlen(buf); strcpy(buf+len_buf, pair_2_str); strcpy(buf+len_buf+strlen(pair_2_str), file_name_in+len_buf);
	FILE* file_pair_2_out=fopen(buf, "w");

	FILE* file_fastq_in=fopen(file_name_in, "r");

	int cur_info=0;

	char first_char=GetFirstChar(file_name_in, file_fastq_in, buf, info_line[cur_info]);
	if(first_char==0) return 0;

	while(1)
	{
		int seq_len=GetReadSequence(file_fastq_in, first_char, seq, buf, info_line[1-cur_info], quality_info, quality_scores, unknown_char);
		if(seq_len==0) break;

		FILE* cur_file=file_pair_1_out;
		if(cur_info) cur_file=file_pair_2_out;

		if(out_fasta) info_line[cur_info][0]='>';
		if(seq_len>0) fprintf(cur_file, "%s\n%s\n", info_line[cur_info], seq);

		if(!out_fasta && quality_info[0]) {if(no_qual_info) quality_info[1]=0; fprintf(cur_file, "%s\n", quality_info);}

		if(adjust_qual && !out_fasta && quality_scores[0]=='@') quality_scores[0]++;

		if(!out_fasta && quality_scores[0]) fprintf(cur_file, "%s\n", quality_scores);

		cur_info=1-cur_info;
	}

	fclose(file_fastq_in);

	fclose(file_pair_1_out);
	fclose(file_pair_2_out);

	delete[] buf;
	delete[] info_line[0]; delete[] info_line[1];
	delete[] seq;
	delete[] quality_info;
	delete[] quality_scores;

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int GetPairIDs(const char* file_name_in)
{
	char* buf=new char[MAX_READ_LEN+1];
	char* info_line[2]; info_line[0]=new char[MAX_READ_LEN+1]; info_line[1]=new char[MAX_READ_LEN+1];
	char* seq=new char[MAX_READ_LEN+1];
	char* quality_info=new char[MAX_READ_LEN+1];
	char* quality_scores=new char[MAX_READ_LEN+1];

	const char* pairids_str="_pids";

	strcpy(buf, file_name_in); int len_buf=strlen(buf); while(len_buf>0 && buf[len_buf-1]!='.') len_buf--; len_buf--;
	if(len_buf<=0) len_buf=strlen(buf); strcpy(buf+len_buf, pairids_str); strcpy(buf+len_buf+strlen(pairids_str), file_name_in+len_buf);
	FILE* file_pairids_out=fopen(buf, "w");

	FILE* file_fastq_in=fopen(file_name_in, "r");

	int cur_info=0;

	char first_char=GetFirstChar(file_name_in, file_fastq_in, buf, info_line[cur_info]);
	if(first_char==0) return 0;

	while(1)
	{
		int seq_len=GetReadSequence(file_fastq_in, first_char, seq, buf, info_line[1-cur_info], quality_info, quality_scores, 0);
		if(seq_len==0) break;

		int j, n=strlen(info_line[cur_info]);
		for(j=1;j<n;j++) if(info_line[cur_info][j]==' ' || info_line[cur_info][j]=='\t') {info_line[cur_info][j]=0; break;}

		if(cur_info==0) fprintf(file_pairids_out, "%s ", info_line[cur_info]+1);
		else fprintf(file_pairids_out, "%s\n", info_line[cur_info]+1);

		cur_info=1-cur_info;
	}

	fclose(file_fastq_in);

	fclose(file_pairids_out);

	delete[] buf;
	delete[] info_line[0]; delete[] info_line[1];
	delete[] seq;
	delete[] quality_info;
	delete[] quality_scores;

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int GetHShrecInt(char* s)
{
	int i=0;
	while(i<3) {if(*s==' ') i++; s++;}
	i=0;
	while((*s>='0') && (*s<='9')) {i=i*10+((*s)-'0'); s++;}
	return i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MergeHshrecFiles(const char* file_name_a, const char* file_name_b, const char* file_name_out, char unknown_char)
{
	FILE* file_fastq_a=fopen(file_name_a, "r");
	FILE* file_fastq_b=fopen(file_name_b, "r");
	FILE* file_fasta_out=fopen(file_name_out, "w");

	int cur_info_a=0;
	int cur_info_b=0;
	char* buf_a=new char[MAX_READ_LEN+1];
	char* buf_b=new char[MAX_READ_LEN+1];
	char* info_line_a[2]; info_line_a[0]=new char[MAX_READ_LEN+1]; info_line_a[1]=new char[MAX_READ_LEN+1];
	char* info_line_b[2]; info_line_b[0]=new char[MAX_READ_LEN+1]; info_line_b[1]=new char[MAX_READ_LEN+1];

	char first_char_a=GetFirstChar(file_name_a, file_fastq_a, buf_a, info_line_a[cur_info_a]);
	if(first_char_a==0) return 0;

	char first_char_b=GetFirstChar(file_name_b, file_fastq_b, buf_b, info_line_b[cur_info_b]);
	if(first_char_b==0) return 0;

	char* seq_a=new char[MAX_READ_LEN+1];
	char* seq_b=new char[MAX_READ_LEN+1];

	int finished_a=0;
	int finished_b=0;

	int seq_len_a=0;
	int seq_len_b=0;

	int ia=GetHShrecInt(info_line_a[cur_info_a]);
	int ib=GetHShrecInt(info_line_b[cur_info_b]);

	while(1)
	{
		if(ia<ib)
		{
			if(!finished_a)
			{
				seq_len_a=GetReadSequence(file_fastq_a, first_char_a, seq_a, buf_a, info_line_a[1-cur_info_a], 0, 0, unknown_char);
				if(seq_len_a>0) fprintf(file_fasta_out, "%s\n%s\n", info_line_a[cur_info_a], seq_a);
				if(seq_len_a==0) finished_a=1;
				cur_info_a=1-cur_info_a;
				if(seq_len_a==0) ia=0x7FFFFFFF;
				else ia=GetHShrecInt(info_line_a[cur_info_a]);
			}
		}
		else
		{
			if(!finished_b)
			{
				seq_len_b=GetReadSequence(file_fastq_b, first_char_b, seq_b, buf_b, info_line_b[1-cur_info_b], 0, 0, unknown_char);
				if(seq_len_b>0) fprintf(file_fasta_out, "%s\n%s\n", info_line_b[cur_info_b], seq_b);
				if(seq_len_b==0) finished_b=1;
				cur_info_b=1-cur_info_b;
				if(seq_len_b==0) ib=0x7FFFFFFF;
				else ib=GetHShrecInt(info_line_b[cur_info_b]);
			}
		}

		if(finished_a && finished_b) break;
	}

	fclose(file_fastq_a);
	fclose(file_fastq_b);
	fclose(file_fasta_out);

	delete[] buf_a;
	delete[] buf_b;
	delete[] info_line_a[0]; delete[] info_line_a[1];
	delete[] info_line_b[0]; delete[] info_line_b[1];
	delete[] seq_a;
	delete[] seq_b;

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MergeQuakeFiles(const char* file_name_o, const char* file_name_a, const char* file_name_b, const char* file_name_out, char unknown_char)
{
	FILE* file_fastq_o=fopen(file_name_o, "r");
	FILE* file_fastq_a=fopen(file_name_a, "r");
	FILE* file_fastq_b=fopen(file_name_b, "r");
	FILE* file_fasta_out=fopen(file_name_out, "w");

	int cur_info_o=0;
	int cur_info_a=0;
	int cur_info_b=0;
	char* buf_o=new char[MAX_READ_LEN+1];
	char* buf_a=new char[MAX_READ_LEN+1];
	char* buf_b=new char[MAX_READ_LEN+1];
	char* info_line_o[2]; info_line_o[0]=new char[MAX_READ_LEN+1]; info_line_o[1]=new char[MAX_READ_LEN+1];
	char* info_line_a[2]; info_line_a[0]=new char[MAX_READ_LEN+1]; info_line_a[1]=new char[MAX_READ_LEN+1];
	char* info_line_b[2]; info_line_b[0]=new char[MAX_READ_LEN+1]; info_line_b[1]=new char[MAX_READ_LEN+1];

	char first_char_o=GetFirstChar(file_name_o, file_fastq_o, buf_o, info_line_o[cur_info_o]);
	if(first_char_o==0) return 0;

	char first_char_a=GetFirstChar(file_name_a, file_fastq_a, buf_a, info_line_a[cur_info_a]);
	if(first_char_a==0) return 0;

	char first_char_b=GetFirstChar(file_name_b, file_fastq_b, buf_b, info_line_b[cur_info_b]);
	if(first_char_b==0) return 0;

	char* seq_o=new char[MAX_READ_LEN+1];
	char* seq_a=new char[MAX_READ_LEN+1];
	char* seq_b=new char[MAX_READ_LEN+1];

	char* quality_info_a=new char[MAX_READ_LEN+1];
	char* quality_scores_a=new char[MAX_READ_LEN+1];

	char* quality_info_b=new char[MAX_READ_LEN+1];
	char* quality_scores_b=new char[MAX_READ_LEN+1];

	int finished_a=0;
	int finished_b=0;

	int seq_len_a=0;
	int seq_len_b=0;

	GetReadSequence(file_fastq_o, first_char_o, seq_o, buf_o, info_line_o[1-cur_info_o], 0, 0, unknown_char);

	seq_len_a=GetReadSequence(file_fastq_a, first_char_a, seq_a, buf_a, info_line_a[1-cur_info_a], quality_info_a, quality_scores_a, unknown_char);
	if(seq_len_a==0) finished_a=1;

	seq_len_b=GetReadSequence(file_fastq_b, first_char_b, seq_b, buf_b, info_line_b[1-cur_info_b], quality_info_b, quality_scores_b, unknown_char);
	if(seq_len_b==0) finished_b=1;

	while(!finished_a || !finished_b)
	{
		if(!finished_a && strncmp(info_line_o[cur_info_o], info_line_a[cur_info_a], strlen(info_line_o[cur_info_o]))==0)
		{
			fprintf(file_fasta_out, "%s\n%s\n", info_line_a[cur_info_a], seq_a);

			if(quality_info_a[0])
				fprintf(file_fasta_out, "%s\n%s\n", quality_info_a, quality_scores_a);

			cur_info_a=1-cur_info_a;
			seq_len_a=GetReadSequence(file_fastq_a, first_char_a, seq_a, buf_a, info_line_a[1-cur_info_a], quality_info_a, quality_scores_a, unknown_char);
			if(seq_len_a==0) finished_a=1;

			cur_info_o=1-cur_info_o;
			GetReadSequence(file_fastq_o, first_char_o, seq_o, buf_o, info_line_o[1-cur_info_o], 0, 0, unknown_char);
		}
		else if(!finished_b && strncmp(info_line_o[cur_info_o], info_line_b[cur_info_b], strlen(info_line_o[cur_info_o]))==0)
		{
			fprintf(file_fasta_out, "%s\n%s\n", info_line_b[cur_info_b], seq_b);

			if(quality_info_b[0])
				fprintf(file_fasta_out, "%s\n%s\n", quality_info_b, quality_scores_b);

			cur_info_b=1-cur_info_b;
			seq_len_b=GetReadSequence(file_fastq_b, first_char_b, seq_b, buf_b, info_line_b[1-cur_info_b], quality_info_b, quality_scores_b, unknown_char);
			if(seq_len_b==0) finished_b=1;

			cur_info_o=1-cur_info_o;
			GetReadSequence(file_fastq_o, first_char_o, seq_o, buf_o, info_line_o[1-cur_info_o], 0, 0, unknown_char);
		}
		else
		{
			printf("Error Unknown Header\n");
			printf("o=%s\na=%s\nb=%s\n", info_line_o[cur_info_o], info_line_a[cur_info_a], info_line_b[cur_info_b]);
			fflush(NULL);
			exit(0);
		}
	}

	fclose(file_fastq_o);
	fclose(file_fastq_a);
	fclose(file_fastq_b);
	fclose(file_fasta_out);

	delete[] buf_o;
	delete[] buf_a;
	delete[] buf_b;
	delete[] info_line_o[0]; delete[] info_line_o[1];
	delete[] info_line_a[0]; delete[] info_line_a[1];
	delete[] info_line_b[0]; delete[] info_line_b[1];
	delete[] seq_o;
	delete[] seq_a;
	delete[] seq_b;

	delete[] quality_info_a;
	delete[] quality_scores_a;

	delete[] quality_info_b;
	delete[] quality_scores_b;

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool AdjustQualString(char* qual)
{
	char* nq=qual;
	while(*qual)
	{
		while(*qual && (*qual<'0' || *qual>'9')) qual++;
		if(!(*qual)) break;
		int num=*qual-'0'; qual++;
		while(*qual>='0' && *qual<='9') {num=num*10+*qual-'0'; qual++;}
		if(num<0 || num>=QUAL_STR_SIZE) {printf("Unknown quality value (%d)\n", num); fflush(NULL); return false;}
		*nq++=QUAL_STR[num];
	}
	*nq=0;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MergeFastaQualFiles(const char* file_name_a, const char* file_name_b, const char* file_name_out, char unknown_char, int no_qual_info)
{
	FILE* file_fastq_a=fopen(file_name_a, "r");
	FILE* file_fastq_b=fopen(file_name_b, "r");
	FILE* file_fasta_out=fopen(file_name_out, "w");

	int cur_info=0;
	char* buf_a=new char[MAX_READ_LEN+1];
	char* buf_b=new char[MAX_READ_LEN+1];
	char* info_line_a[2]; info_line_a[0]=new char[MAX_READ_LEN+1]; info_line_a[1]=new char[MAX_READ_LEN+1];
	char* info_line_b[2]; info_line_b[0]=new char[MAX_READ_LEN+1]; info_line_b[1]=new char[MAX_READ_LEN+1];

	char first_char_a=GetFirstChar(file_name_a, file_fastq_a, buf_a, info_line_a[cur_info]);
	if(first_char_a==0) return 0;

	char first_char_b=GetFirstChar(file_name_b, file_fastq_b, buf_b, info_line_b[cur_info]);
	if(first_char_b==0) return 0;

	char* seq_a=new char[MAX_READ_LEN+1];
	char* seq_b=new char[MAX_READ_LEN+1];

	int seq_len_a=0;
	int seq_len_b=0;

	while(1)
	{
		seq_len_a=GetReadSequence(file_fastq_a, first_char_a, seq_a, buf_a, info_line_a[1-cur_info], 0, 0, unknown_char);
		seq_len_b=GetReadSequence(file_fastq_b, first_char_b, seq_b, buf_b, info_line_b[1-cur_info], 0, 0, 0, true);
		if(seq_len_a==0 || seq_len_b==0) break;

		info_line_a[cur_info][0]='@';
		fprintf(file_fasta_out, "%s\n%s\n", info_line_a[cur_info], seq_a);

		AdjustQualString(seq_b);

		info_line_b[cur_info][0]='+';
		if(no_qual_info) info_line_b[cur_info][1]=0;
		fprintf(file_fasta_out, "%s\n%s\n", info_line_b[cur_info], seq_b);

		cur_info=1-cur_info;
	}

	fclose(file_fastq_a);
	fclose(file_fastq_b);
	fclose(file_fasta_out);

	delete[] buf_a;
	delete[] buf_b;
	delete[] info_line_a[0]; delete[] info_line_a[1];
	delete[] info_line_b[0]; delete[] info_line_b[1];
	delete[] seq_a;
	delete[] seq_b;

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int IntToString(int num, char* str)
{
	int str_len=0;
	if(num==0) {str[str_len++]='0'; str[str_len]=0; return str_len;}
	while(num!=0)
	{
		int nlet=num%10;
		str[str_len++]=(nlet+'0');
		num=num/10;
	}
	str[str_len]=0;
	Reverse(str, str_len);
	return str_len;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ModifyInfoLength(char* info, int new_int)
{
	char* new_info=new char[MAX_READ_LEN+1];
	int new_info_len=0;

	int i, j=-1, k, info_len=strlen(info);
	const char* str_length[]={" len=", " length="};
	int num_str_length=sizeof(str_length)/sizeof(const char*);

	bool found=false;

	for(k=0;k<num_str_length && !found;k++)
	{
		for(i=0;i<info_len;i++)
		{
			int str_length_len=strlen(str_length[k]);

			for(j=0;j<str_length_len && i+j<info_len;j++)
				if(info[i+j]!=str_length[k][j])
					break;

			if(j==str_length_len)
			{
				j+=i;
				found=true;
				break;
			}
		}
	}

	if(found && j>=0)
	{
		memcpy(new_info, info, j);
		new_info_len=j;
		new_info_len+=IntToString(new_int, new_info+new_info_len);
		new_info[new_info_len]=0;

		for(;j<info_len;j++) if(info[j]<'0' || info[j]>'9') break;
		if(j<info_len) strcpy(new_info+new_info_len, info+j);

		strcpy(info, new_info);

		//printf("%s\n\n", new_info);
	}

	delete[] new_info;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ReplaceInfoSecondByInfoFirst(const char* file_name_a, const char* file_name_b, const char* file_name_out, int out_fasta, char unknown_char, int adjust_qual, int no_qual_info)
{
	FILE* file_fastq_a=fopen(file_name_a, "r");
	FILE* file_fastq_b=fopen(file_name_b, "r");
	FILE* file_fasta_out=fopen(file_name_out, "w");

	int cur_info=0;
	char* buf_a=new char[MAX_READ_LEN+1];
	char* buf_b=new char[MAX_READ_LEN+1];
	char* info_line_a[2]; info_line_a[0]=new char[MAX_READ_LEN+1]; info_line_a[1]=new char[MAX_READ_LEN+1];
	char* info_line_b[2]; info_line_b[0]=new char[MAX_READ_LEN+1]; info_line_b[1]=new char[MAX_READ_LEN+1];
	char* quality_info_a=new char[MAX_READ_LEN+1];
	char* quality_info_b=new char[MAX_READ_LEN+1];
	char* quality_scores_a=new char[MAX_READ_LEN+1];
	char* quality_scores_b=new char[MAX_READ_LEN+1];

	char first_char_a=GetFirstChar(file_name_a, file_fastq_a, buf_a, info_line_a[cur_info]);
	if(first_char_a==0) return 0;

	char first_char_b=GetFirstChar(file_name_b, file_fastq_b, buf_b, info_line_b[cur_info]);
	if(first_char_b==0) return 0;

	char* seq_a=new char[MAX_READ_LEN+1];
	char* seq_b=new char[MAX_READ_LEN+1];

	int finished_a=0;
	int finished_b=0;

	while(1)
	{
		int seq_len_a=0; if(!finished_a) seq_len_a=GetReadSequence(file_fastq_a, first_char_a, seq_a, buf_a, info_line_a[1-cur_info], quality_info_a, quality_scores_a, unknown_char);
		int seq_len_b=0; if(!finished_b) seq_len_b=GetReadSequence(file_fastq_b, first_char_b, seq_b, buf_b, info_line_b[1-cur_info], quality_info_b, quality_scores_b, unknown_char);
		if(seq_len_a==0 && seq_len_b==0) break;
		if(seq_len_a==0) finished_a=1; if(seq_len_b==0) finished_b=1;

		if(out_fasta || first_char_b=='>') info_line_a[cur_info][0]='>';
		if(seq_len_b>0) {ModifyInfoLength(info_line_a[cur_info], seq_len_b); fprintf(file_fasta_out, "%s\n%s\n", info_line_a[cur_info], seq_b);}

		if(!out_fasta && first_char_b=='@' && quality_scores_b[0])
		{
			if(quality_info_a[0]){if(no_qual_info) quality_info_a[1]=0; ModifyInfoLength(quality_info_a, seq_len_b); fprintf(file_fasta_out, "%s\n", quality_info_a);}

			if(adjust_qual && quality_scores_b[0]=='@') quality_scores_b[0]++;
			if(quality_scores_b[0]) fprintf(file_fasta_out, "%s\n", quality_scores_b);
		}

		cur_info=1-cur_info;
	}

	fclose(file_fastq_a);
	fclose(file_fastq_b);
	fclose(file_fasta_out);

	delete[] buf_a;
	delete[] buf_b;
	delete[] info_line_a[0]; delete[] info_line_a[1];
	delete[] info_line_b[0]; delete[] info_line_b[1];
	delete[] seq_a;
	delete[] seq_b;
	delete[] quality_info_a;
	delete[] quality_info_b;
	delete[] quality_scores_a;
	delete[] quality_scores_b;

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MergeFilesMain(int argc, char* argv[])
{
	printf("Merge Started\n"); fflush(NULL);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	int			default_paired=0;
	int			default_sample=1000;
	int			default_out_fasta=0;
	int			default_out_fastq=0;
	int			default_out_qual=0;
	int			default_one_seq=0;
	int			default_hshrec=0;
	int			default_quake=0;
	int			default_fasta_qual=0;
	int			default_replace_info=0;
	int			default_split_paired=0;
	int			default_split_pairs=0;
	int			default_out_pair_ids=0;
	int			default_adjust_qual=0;
	int			default_no_qual_info=0;
	char		default_unknown_char=0;

	const char* default_input_dir="./";

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	int			paired=default_paired; //bool_set_paired; set_paired=false;
	int			sample=default_sample;
	int			out_fasta=default_out_fasta; //bool_set_out_fasta; set_out_fasta=false;
	int			out_fastq=default_out_fastq; //bool_set_out_fastq; set_out_fastq=false;
	int			out_qual=default_out_qual; //bool_set_out_qual; set_out_qual=false;
	int			one_seq=default_one_seq; //bool_set_one_seq; set_one_seq=false;
	int			hshrec=default_hshrec; //bool_set_hshrec; set_hshrec=false;
	int			quake=default_quake; //bool_set_quake; set_quake=false;
	int			fasta_qual=default_fasta_qual; //bool_set_fasta_qual; set_fasta_qual=false;
	int			replace_info=default_replace_info; //bool_set_replace_info; set_replace_info=false;
	int			split_paired=default_split_paired; //bool_split_paired; set_split_paired=false;
	int			split_pairs=default_split_pairs; //bool_split_pairs; set_split_pairs=false;
	int			out_pair_ids=default_out_pair_ids;
	int			adjust_qual=default_adjust_qual; //bool_set_adjust_qual; set_adjust_qual=false;
	int			no_qual_info=default_no_qual_info; //bool_set_no_qual_info; set_no_qual_info=false;
	char		unknown_char=default_unknown_char; //bool_set_unknown_char; set_unknown_char=false;

	char 		input_dir[MAX_PATH_LEN+1]; strcpy(input_dir, default_input_dir); //bool_set_input_dir; set_input_dir=false;

	vector<FileName> user_input_file_names; bool set_user_input_file_names; set_user_input_file_names=false;
	char		user_merged_file_name[MAX_PATH_LEN+1]; user_merged_file_name[0]=0; bool set_user_merged_file_name; set_user_merged_file_name=false;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*******************************************************************************************************
	// GetData from user here, be sure that slash is added to directory names

	int i;

	for(i=2;i<argc;i++)
	{
		bool processed=false;
		char* str=argv[i];

		if(strcmp(str,"-paired")==0) {paired=1; processed=true;}// set_paired=true;}
		else if(strcmp(str,"-adjustqual")==0) {adjust_qual=1; processed=true;} //set_adjust_qual=true;}
		else if(strcmp(str,"-noqualinfo")==0) {no_qual_info=1; processed=true;} //no_qual_info=true;}
		else if(strcmp(str,"-hshrec")==0) {hshrec=1; processed=true;}// set_hshrec=true;}
		else if(strcmp(str,"-quake")==0) {quake=1; processed=true;}// set_quake=true;}
		else if(strcmp(str,"-fastaqual")==0) {fasta_qual=1; processed=true;}// set_fasta_qual=true;}
		else if(strcmp(str,"-replaceinfo")==0) {replace_info=1; processed=true;}// set_replace_info=true;}
		else if(strcmp(str,"-splitpaired")==0) {split_paired=1; processed=true;}// set_split_paired=true;}
		else if(strcmp(str,"-splitpairs")==0) {split_pairs=1; processed=true;}// set_split_paired=true;}
		else if(strcmp(str,"-outpairids")==0) {out_pair_ids=1; processed=true;}// set_split_paired=true;}
		else if(strcmp(str,"-oneseq")==0) {one_seq=1; processed=true;}// set_one_seq=true;}
		else if(strcmp(str,"-outfasta")==0) {out_fasta=1; processed=true;}// set_out_fasta=true;}
		else if(strcmp(str,"-outfastq")==0) {out_fastq=1; processed=true;}// set_out_fastq=true;}
		else if(strcmp(str,"-outqual")==0) {out_qual=1; processed=true;}// set_out_qual=true;}
		else
		{
			int j, len=strlen(str);
			for(j=0;j<len;j++) if(str[j]=='=') break;

			if(j<len-1)
			{
				char* opt=(char*)malloc(j+1);
				memcpy(opt, str, j); opt[j]=0;
				char* val=str+j+1;

				if(strcmp(opt,"-sample")==0)
				{
					sscanf(val, "%d", &sample);
					processed=true; //set_sample=true;
				}
				else if(strcmp(opt,"-unknownchar")==0)
				{
					sscanf(val, "%c", &unknown_char);
					processed=true; //set_unknown_char=true;
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
				else if(strcmp(opt,"-mergedfile")==0) // REQUIRED
				{
					strcpy(user_merged_file_name, val);
					processed=true; set_user_merged_file_name=true;
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
	//*******************************************************************************************************

	vector<FileName> input_file_names;
	char merged_file_name[MAX_PATH_LEN+1];

	for(i=0;i<(int)user_input_file_names.size();i++)
	{
		FileName file;
		GetFileName(user_input_file_names[i].s, input_dir, file.s);
		input_file_names.push_back(file);
	}

	GetFileName(user_merged_file_name, input_dir, merged_file_name);

	if(sample<=0 || sample>1000)
		sample=1000;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*******************************************************************************************************

	// check erroneous and completeness of values here
	// check opening files here

	bool accepted=true;

	if(!set_user_input_file_names)
	{
		printf("Please specify fasta or fastq input file(s) using \"-inputfile=filename\"\n"); fflush(NULL);
		accepted=false;
	}
	else if(!set_user_merged_file_name && !split_paired && !split_pairs && !out_pair_ids)
	{
		printf("Please specify output merged file using \"-mergedfile=filename\"\n"); fflush(NULL);
		accepted=false;
	}
	else if(set_user_merged_file_name && (split_paired || split_pairs || out_pair_ids))
	{
		printf("Merged file option can not be used with \"-splitpaired\"\n"); fflush(NULL);
		accepted=false;
	}
	else if(paired && input_file_names.size()!=2)
	{
		printf("Number of input files must be 2 for the \"-paired\" option."); fflush(NULL);
		accepted=false;
	}
	else if(paired && one_seq)
	{
		printf("You cannot use \"-oneseq\" and \"-paired\" options at the same time."); fflush(NULL);
		accepted=false;
	}
	else if(hshrec && input_file_names.size()!=2)
	{
		printf("Number of input files must be 2 for the \"-hshrec\" option."); fflush(NULL);
		accepted=false;
	}
	else if(quake && input_file_names.size()!=3)
	{
		printf("Number of input files must be 3 for the \"-quake\" option."); fflush(NULL);
		accepted=false;
	}
	else if(fasta_qual && input_file_names.size()!=2)
	{
		printf("Number of input files must be 2 for the \"-fastaqual\" option."); fflush(NULL);
		accepted=false;
	}
	else if(replace_info && input_file_names.size()!=2)
	{
		printf("Number of input files must be 2 for the \"-replaceinfo\" option."); fflush(NULL);
		accepted=false;
	}
	else if(out_qual && split_paired)
	{
		printf("You can not use \"outqual\" option with \"splitpaired\" option."); fflush(NULL);
		accepted=false;
	}
	else if(split_paired && input_file_names.size()!=1)
	{
		printf("You can use only one file with \"splitpaired\" option."); fflush(NULL);
		accepted=false;
	}
	else if(split_pairs && input_file_names.size()!=1)
	{
		printf("You can use only one file with \"splitpairs\" option."); fflush(NULL);
		accepted=false;
	}
	else if(out_pair_ids && input_file_names.size()!=1)
	{
		printf("You can use only one file with \"outpairids\" option."); fflush(NULL);
		accepted=false;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(i=0;accepted && i<(int)input_file_names.size();i++)
	{
		FILE* file=fopen(input_file_names[i].s, "r");
		if(!file)
		{
			printf("Cannot open input file [%s] for reading.\n", input_file_names[i].s); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(accepted && !split_paired && !split_pairs && !out_pair_ids)
	{
		FILE* file=fopen(merged_file_name, "w");
		if(!file)
		{
			printf("Cannot open merged file [%s] for writing.\n", merged_file_name); fflush(NULL);
			accepted=false;
		}
		else fclose(file);
	}

	if(!accepted)
		return;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*******************************************************************************************************

	if(hshrec)
	{
		MergeHshrecFiles(input_file_names[0].s, input_file_names[1].s, merged_file_name, unknown_char);
	}
	else if(quake)
	{
		MergeQuakeFiles(input_file_names[0].s, input_file_names[1].s, input_file_names[2].s, merged_file_name, unknown_char);
	}
	else if(fasta_qual)
	{
		MergeFastaQualFiles(input_file_names[0].s, input_file_names[1].s, merged_file_name, unknown_char, no_qual_info);
	}
	else if(replace_info)
	{
		ReplaceInfoSecondByInfoFirst(input_file_names[0].s, input_file_names[1].s, merged_file_name, out_fasta, unknown_char, adjust_qual, no_qual_info);
	}
	else if(paired)
	{
		MergeTwoPairedEndFilesIntoOneFile(input_file_names[0].s, input_file_names[1].s, merged_file_name, out_fasta, unknown_char, adjust_qual, no_qual_info, sample);
	}
	else if(split_paired)
	{
		SplitPaired(input_file_names[0].s, out_fasta, unknown_char, adjust_qual, no_qual_info);
	}
	else if(split_pairs)
	{
		SplitPairsTwoFiles(input_file_names[0].s, out_fasta, unknown_char, adjust_qual, no_qual_info);
	}
	else if(out_pair_ids)
	{
		GetPairIDs(input_file_names[0].s);
	}
	else
	{
		ConcatenateFilesIntoOneFile(input_file_names, merged_file_name, out_fasta, out_qual, one_seq, unknown_char, adjust_qual, no_qual_info, out_fastq, sample);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Merge Finished\n\n"); fflush(NULL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MergeHelp()
{
	printf("-Run \"karect -merge [options_list]\" for the merge tool.\n");
	printf("-Available options: (i=integer, c=char, f=file, s=directory)\n");
	printf("-Essential options:\n");
	printf("  \"-inputfile=f\":   Specify an input fasta/fastq file. This option can be repeated for multiple files.\n");
	printf("  \"-mergedfile=f\":  Specify the output file.\n");
	printf("-Basic options:\n");
	printf("  \"-inputdir=s\":    Specify the files directory. Ignored if file paths are complete [Default=.].\n");
	printf("  \"-paired\":        Interlace two paired-end or mate-pairs files.\n");
	printf("  \"-sample=i\":      Reserve only sample/1000 of reads [Default=1000].\n");
	printf("  \"-oneseq\":        Concatenate all sequences into one sequence (works with fasta only).\n");
	printf("  \"-noqualinfo\":    Do not write fastq quality information line (write only +).\n");
	printf("-Advanced options:\n");
	printf("  \"-unknownchar=c\": Specify the character to be put instead of unknown characters\n"
		   "                      (known characters are \"ACGTacgt\") [Default: do not change unknown characters].\n");
	printf("  \"-outfasta\":      Output fasta file [Default: output file type is the same as input file type].\n");
	printf("  \"-outfastq\":      Output fastq file [Default: output file type is the same as input file type].\n");
	printf("  \"-outqual\":       Output separate quality file.\n");
	printf("  \"-adjustqual\":    Do not allow \"@\" to exist as the first quality score for any read.\n");
	printf("  \"-hshrec\":        Concatenate HSHREC output files.\n");
	printf("  \"-quake\":         Concatenate Quake output files.\n");
	printf("  \"-fastaqual\":     Merge fasta and qual files into fastq.\n");
	printf("  \"-replaceinfo\":   Replace info lines of the second file using the ones of the first file.\n");
	printf("  \"-splitpaired\":   Split into two files: fragment, paired.\n");
	printf("  \"-splitpairs\":    Split interlaced pairs into two files: pair1, pair2.\n");
	printf("  \"-outpairids\":    Output pair IDs file to be used by celera.\n");
	printf("-Example:\n"
		   "        ./karect -merge -inputdir=/sra_data -inputfile=SRR001666_1.fasta -inputfile=SRR001666_2.fasta\n"
		   "              -mergedfile=./SRR001666_all.fasta\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
