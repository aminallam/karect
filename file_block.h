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

#ifndef __FILE_BLOCK
#define __FILE_BLOCK

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

char IsEndLineChar[1<<8]; // \r \n
char IsSpaceChar[1<<8]; // space, \t, \v, \f, \r, \n

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GlobalInitFileBlock()
{
	int i;
	for(i=0;i<(1<<8);i++) IsEndLineChar[i]=0;
	IsEndLineChar['\r']=1;
	IsEndLineChar['\n']=1;

	for(i=0;i<(1<<8);i++) IsSpaceChar[i]=0;
	IsSpaceChar[' ']=1;
	IsSpaceChar['\r']=1;
	IsSpaceChar['\n']=1;
	IsSpaceChar['\t']=1;
	IsSpaceChar['\v']=1;
	IsSpaceChar['\f']=1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct StringBlock
{
	char*			str;
	int				size;
	StringBlock* 	next;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct BlockedString
{
	int				block_size;
	long long		total_size;
	StringBlock* 	head;
	StringBlock* 	tail;

	void Initialize(int _block_size)
	{
		block_size=_block_size;
		total_size=0;
		head=tail=AllocateStringBlock();
	}

	void Reset()
	{
		StringBlock* cur=head->next;
		while(cur)
		{
			StringBlock* next=cur->next;
			FreeHeap(cur->str, block_size);
			FreeHeap(cur, sizeof(StringBlock));
			cur=next;
		}
		tail=head;
		head->size=0;
		head->next=0;
		total_size=0;
	}

	void Destroy()
	{
		while(head)
		{
			StringBlock* next=head->next;
			FreeHeap(head->str, block_size);
			FreeHeap(head, sizeof(StringBlock));
			head=next;
		}
	}

	char GetFirstChar()
	{
		if(head->size==0) return 0;
		return head->str[0];
	}

	long long GetSize()
	{
		return total_size;
	}

	StringBlock* AllocateStringBlock()
	{
		StringBlock* block=(StringBlock*)AllocateHeap(sizeof(StringBlock));
		block->str=(char*)AllocateHeap(block_size);
		block->next=0; block->size=0;
		return block;
	}

	long long GetString(char* str, bool append_zero)
	{
		long long cur_size=0;
		StringBlock* cur_node=head;
		while(cur_node)
		{
			StringBlock* next_node=cur_node->next;
			memcpy(str+cur_size, cur_node->str, cur_node->size);
			cur_size+=cur_node->size;
			cur_node=next_node;
		}
		if(append_zero) str[cur_size]=0;
		return cur_size;
	}

	void Append(char* str, long long n)
	{
		long long cur_size=0;

		while(cur_size<n)
		{
			//if(!cur_tail) cur_tail=tail=AllocateStringBlock();
			int rem=block_size-tail->size;
			if(rem>n-cur_size) rem=n-cur_size;

			memcpy(tail->str+tail->size, str+cur_size, rem);
			tail->size+=rem;
			cur_size+=rem;

			if(tail->size==block_size)
			{
				tail->next=AllocateStringBlock();
				tail=tail->next;
			}
		}

		total_size+=n;
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SolFileInOut
{
	FILE* file;
	char* another_block;
	int	  old_cur_pos;
};

void* threadfileout(void* arg)
{
	SolFileInOut* sol=(SolFileInOut*)arg;
	fwrite(sol->another_block, 1, sol->old_cur_pos, sol->file);
	return NULL;
}

void* threadfilein(void* arg)
{
	SolFileInOut* sol=(SolFileInOut*)arg;
	sol->old_cur_pos=fread(sol->another_block, 1, sol->old_cur_pos, sol->file);
	return NULL;
}

struct BufferedOutFile
{
	FILE*		file;
	int			block_size;
	int			cur_pos;
	char*		block;
	char*		another_block;

	bool		thread_running;
	pthread_t	thread_ID;
	void*		exit_status;
	SolFileInOut sol;

	BufferedOutFile() {file=0; block=0; cur_pos=0; block_size=0;}

	void WriteNumBytes(void* _buf, int n)
	{
		int i=0;
		char* buf=(char*)_buf;
		int num_avail=block_size-cur_pos;
		if(num_avail>n) num_avail=n;
		while(i<n)
		{
			memcpy(block+cur_pos, buf+i, num_avail);
			i+=num_avail;
			cur_pos+=num_avail;

			if(cur_pos==block_size)
			{
				Flush();
				num_avail=block_size;
				if(num_avail>n-i) num_avail=n-i;
			}
		}
	}

	void Flush()
	{
		//if(cur_pos) fwrite(block, 1, cur_pos, file);
		//cur_pos=0;
		if(cur_pos)
		{
			if(thread_running) pthread_join(thread_ID, &exit_status);

			thread_running=true;

			char* temp=another_block; another_block=block; block=temp;
			int old_cur_pos=cur_pos; cur_pos=0;

			sol.file=file;
			sol.another_block=another_block;
			sol.old_cur_pos=old_cur_pos;

			pthread_create(&thread_ID, NULL, threadfileout, &sol);
			//fwrite(another_block, 1, old_cur_pos, file);
		}
	}

	void Initialize(FILE* _file, int _block_size)
	{
		file=_file;
		block_size=_block_size;
		block=(char*)malloc(block_size);
		another_block=(char*)malloc(block_size);
		cur_pos=0;
		thread_running=false;
	}

	void Destroy()
	{
		Flush();
		if(thread_running) pthread_join(thread_ID, &exit_status);
		free(block);
		free(another_block);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct BufferedInFile
{
	FILE*		file;
	int			block_size;
	int			cur_pos;
	int			cur_size;
	char*		block;
	char*		another_block;
	bool		finished;

	bool		thread_running;
	pthread_t	thread_ID;
	void*		exit_status;
	SolFileInOut sol;

	BufferedInFile() {file=0; block=0; cur_pos=0; block_size=0; cur_size=0; finished=false;}

	int ReadNumBytes(void* _buf, int n)
	{
		int i=0;
		char* buf=(char*)_buf;
		int num_avail=cur_size-cur_pos;
		if(num_avail>n) num_avail=n;
		while(i<n && (!finished || cur_pos<cur_size))
		{
			if(cur_pos==cur_size)
			{
				ReadBlock();
				if(cur_size<block_size) finished=true;
				num_avail=cur_size;
				if(num_avail>n-i) num_avail=n-i;
			}
			memcpy(buf+i, block+cur_pos, num_avail);
			i+=num_avail;
			cur_pos+=num_avail;
		}
		return i;
	}

	int ReadUntil(void* _buf, char delim) // read delim and discard it, put \0 at end
	{
		int i=0,j;
		char* buf=(char*)_buf;
		int num_avail=cur_size-cur_pos;
		bool delim_found=false;
		for(j=0;j<num_avail;j++) if(block[cur_pos+j]==delim) {num_avail=j; delim_found=true; break;}
		//if(num_avail>n) num_avail=n;
		while(!finished || cur_pos<cur_size)
		{
			if(cur_pos==cur_size)
			{
				ReadBlock();
				if(cur_size<block_size) finished=true;
				num_avail=cur_size;
				for(j=0;j<num_avail;j++) if(block[cur_pos+j]==delim) {num_avail=j; delim_found=true; break;}
				//if(num_avail>n-i) num_avail=n-i;
			}
			if(num_avail>0) memcpy(buf+i, block+cur_pos, num_avail);
			i+=num_avail;
			cur_pos+=num_avail;
			if(delim_found) break;
		}
		buf[i]=0;
		cur_pos++; // skip delimiter
		return i;
	}

	bool ReadLine(BlockedString& str, char exclude_start=0) // read line, skip spaces after it, append line to str if not starting with exculde_start
	{
		int j, num_avail=cur_size-cur_pos;
		bool end_line_found=false;
		for(j=0;j<num_avail;j++) if(IsEndLineChar[(unsigned char)block[cur_pos+j]]) {num_avail=j; end_line_found=true; break;}
		char first_char=0;
		while(!finished || cur_pos<cur_size)
		{
			if(cur_pos==cur_size)
			{
				ReadBlock();
				if(cur_size<block_size) finished=true;
				num_avail=cur_size;
				for(j=0;j<num_avail;j++) if(IsEndLineChar[(unsigned char)block[cur_pos+j]]) {num_avail=j; end_line_found=true; break;}
			}
			if(num_avail>0)
			{
				if(!first_char) first_char=block[cur_pos];
				if(first_char!=exclude_start) str.Append(block+cur_pos, num_avail);
			}
			cur_pos+=num_avail;
			if(end_line_found) break;
		}
		SkipSpaces();
		return end_line_found;
	}

	void SkipSpaces()
	{
		int j, num_avail=cur_size-cur_pos;
		bool non_space_found=false;
		for(j=0;j<num_avail;j++) if(!IsSpaceChar[(unsigned char)block[cur_pos+j]]) {num_avail=j; non_space_found=true; break;}
		while(!finished || cur_pos<cur_size)
		{
			if(cur_pos==cur_size)
			{
				ReadBlock();
				if(cur_size<block_size) finished=true;
				num_avail=cur_size;
				for(j=0;j<num_avail;j++) if(!IsSpaceChar[(unsigned char)block[cur_pos+j]]) {num_avail=j; non_space_found=true; break;}
			}
			cur_pos+=num_avail;
			if(non_space_found) break;
		}
		if(!non_space_found) cur_pos=cur_size;
	}

	void ReadBlock(bool first_time=false)
	{
		//cur_pos=0;
		//cur_size=fread(block, 1, block_size, file);

		if(!first_time)
		{
			if(thread_running) pthread_join(thread_ID, &exit_status);

			char* temp=another_block; another_block=block; block=temp;
			cur_pos=0;
			cur_size=sol.old_cur_pos;
		}

		thread_running=true;

		sol.file=file;
		sol.another_block=another_block;
		sol.old_cur_pos=block_size;

		pthread_create(&thread_ID, NULL, threadfilein, &sol);
	}

	void Initialize(FILE* _file, int _block_size)
	{
		file=_file;
		block_size=_block_size;
		block=(char*)malloc(block_size);
		another_block=(char*)malloc(block_size);
		cur_pos=cur_size=0;
		finished=false;
		thread_running=false;
		ReadBlock(true);
	}

	void Destroy()
	{
		if(thread_running) pthread_join(thread_ID, &exit_status);
		free(block);
		free(another_block);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
