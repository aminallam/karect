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

#ifndef __STACK
#define __STACK

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class StackNode>
struct Stack
{
	StackNode*	head;
	StackNode*	tail;
};

template<class StackNode>
void InitStack(Stack<StackNode>* stack)
{
	stack->head=stack->tail=(StackNode*)AllocateHeap(sizeof(StackNode));
}

template<class StackNode>
void Push(Stack<StackNode>* stack, StackNode* stack_node)
{
	stack_node->prev_node=stack->tail;
	stack->tail=stack_node;
}

template<class StackNode>
StackNode* Pop(Stack<StackNode>* stack)
{
	if(stack->head==stack->tail) return 0;
	StackNode* ret=stack->tail;
	stack->tail=stack->tail->prev_node;
	return ret;
}

template<class StackNode>
void DestroyStack(Stack<StackNode>* stack)
{
	while(stack->head!=stack->tail)
	{
		StackNode* ret=stack->tail;
		stack->tail=stack->tail->prev_node;
		FreeHeap(ret, sizeof(StackNode));
	}
	FreeHeap(stack->head, sizeof(StackNode));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
