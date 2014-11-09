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

#define _FILE_OFFSET_BITS 64

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <pthread.h>
using namespace std;

#include "memory_utils.h"
#include "file_access.h"
#include "file_block.h"
#include "getreads.h"
#include "stack.h"
#include "edit_dist.h"
#include "suffixfilter.h"
#include "align.h"
#include "evaluate.h"
#include "errorgraph.h"
#include "errorcorrect.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MainHelp()
{
	printf("-Please specify the tool you want to run: (correct-align-eval-merge).\n");
	printf("-Run \"karect -[correct|align|eval|merge]\" to find information about how to run a specific tool.\n");
	printf("  1) correct: a tool to correct assembly reads from fasta/fastq files.\n");
	printf("  2) align:   a tool to align original assembly reads as pre-processing for evaluation.\n");
	printf("  3) eval:    a tool to evaluate assembly reads correction.\n");
	printf("  4) merge:   a tool to concatenate fasta/fastq files.\n");
	printf("-Example: \"./karect -correct\"\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	GlobalInitComp();
	GlobalInitFileBlock();

	if(argc<=1)
	{
		MainHelp();
	}
	else
	{
		if(strcmp(argv[1], "-correct")==0)
		{
			if(argc<=2) CorrectErrorsHelp();
			else CorrectErrorsMain(argc, argv);
		}
		else if(strcmp(argv[1], "-align")==0)
		{
			if(argc<=2) AlignHelp();
			else AlignMain(argc, argv);
		}
		else if(strcmp(argv[1], "-eval")==0)
		{
			if(argc<=2) EvalHelp();
			else EvaluateMain(argc, argv);
		}
		else if(strcmp(argv[1], "-merge")==0)
		{
			if(argc<=2) MergeHelp();
			else MergeFilesMain(argc, argv);
		}
		else
		{
			printf("Unknown tool [%s].\n", argv[1]);
			MainHelp();
		}
	}

	fflush(NULL);
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
