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

#include <cmath>
#include "stack.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct PatNodeConn
{
	int		node; // save
	double	prob;
	double	wacc; // save
};

struct PatNode
{
	char				ch; // save
	vector<PatNodeConn>	adj; // save
	union{int path_prev_node; int num_incoming_edges;
	};
	long double path_weight;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ErrorGraph
{
	int				max_str_len;

	int				min_overlap; // The minimum allowed overlap size (try compute it using ac randomization)
	double			min_overlap_percentage;
	double			max_overlap_error_rate; // The maximum allowed error rate in any overlap

	double			truncate_factor;
	double			min_straight_edge_val;
	double			min_straight_edge_percentage;
	double			min_read_weight;
	int				max_mul_len_matches;

	int				match_type;

	int*			T;

	int				num_nodes; // save (also in hamming)
	int				start_node; // save
	int				end_node; // save

	bool			load_mode;
	bool			save_mode;

	vector<PatNode> pat_nodes; // save

	unsigned short	num_hamming_scores; // save hamming (actually this is pat_len)
	double*			hamming_weights; // save hamming

	bool			confirmed; // save both
	bool			excluded; // save both
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
int CompactGetErrorGraphSize(ErrorGraph* eg)
{
	int i, total_size=0;
	Type node_ind;

	total_size+=sizeof(eg->confirmed);
	if(eg->confirmed) return total_size;

	total_size+=sizeof(eg->excluded)+sizeof(eg->num_nodes)+sizeof(node_ind)+sizeof(node_ind);
	total_size+=eg->num_nodes*(sizeof(char)+sizeof(node_ind));

	for(i=0;i<eg->num_nodes;i++)
	{
		int sz=eg->pat_nodes[i].adj.size();
		total_size+=sz*(sizeof(node_ind)+sizeof(float));
	}

	return total_size;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int GetErrorGraphSize(ErrorGraph* eg)
{
	if(eg->match_type==MATCH_SUBSTITUTE_ONLY)
	{
		int total_size=sizeof(eg->confirmed);
		if(eg->confirmed) return total_size;

		total_size+=sizeof(eg->excluded)+sizeof(eg->num_hamming_scores);

		total_size+=sizeof(float)*NUM_VALID_CHARS*eg->num_hamming_scores;
		return total_size;
	}

	if(eg->num_nodes<0xFF) return CompactGetErrorGraphSize<unsigned char>(eg);
	else if(eg->num_nodes<0xFFFF) return CompactGetErrorGraphSize<unsigned short>(eg);
	return CompactGetErrorGraphSize<unsigned int>(eg);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
void CompactSaveErrorGraph(ErrorGraph* eg, char* buf_eg)
{
	int i;
	int cur=0;

	if(eg->match_type==MATCH_SUBSTITUTE_ONLY)
	{
		int i;
		for(i=0;i<eg->num_hamming_scores*NUM_VALID_CHARS;i++)
		{
			Type val=(Type)eg->hamming_weights[i];
			memcpy(buf_eg+cur, &val, sizeof(Type)); cur+=sizeof(Type);
		}
		return;
	}

	Type j;
	Type node_ind;

	node_ind=eg->start_node; memcpy(buf_eg+cur, &node_ind, sizeof(node_ind)); cur+=sizeof(node_ind);
	node_ind=eg->end_node; memcpy(buf_eg+cur, &node_ind, sizeof(node_ind)); cur+=sizeof(node_ind);

	for(i=0;i<eg->num_nodes;i++)
	{
		memcpy(buf_eg+cur, &eg->pat_nodes[i].ch, sizeof(eg->pat_nodes[i].ch)); cur+=sizeof(eg->pat_nodes[i].ch);
		Type sz=eg->pat_nodes[i].adj.size(); memcpy(buf_eg+cur, &sz, sizeof(sz)); cur+=sizeof(sz);
		for(j=0;j<sz;j++)
		{
			node_ind=eg->pat_nodes[i].adj[j].node; memcpy(buf_eg+cur, &node_ind, sizeof(node_ind)); cur+=sizeof(node_ind);
			float wacc=eg->pat_nodes[i].adj[j].wacc; memcpy(buf_eg+cur, &wacc, sizeof(wacc)); cur+=sizeof(wacc);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SaveErrorGraph(ErrorGraph* eg, char* buf_eg)
{
	int cur=0;

	memcpy(buf_eg+cur, &eg->confirmed, sizeof(eg->confirmed)); cur+=sizeof(eg->confirmed);
	if(eg->confirmed) return;

	memcpy(buf_eg+cur, &eg->excluded, sizeof(eg->excluded)); cur+=sizeof(eg->excluded);

	if(eg->match_type==MATCH_SUBSTITUTE_ONLY)
	{
		memcpy(buf_eg+cur, &eg->num_hamming_scores, sizeof(eg->num_hamming_scores)); cur+=sizeof(eg->num_hamming_scores);
		return CompactSaveErrorGraph<float>(eg, buf_eg+cur);
	}

	memcpy(buf_eg+cur, &eg->num_nodes, sizeof(eg->num_nodes)); cur+=sizeof(eg->num_nodes);

	if(eg->num_nodes<0xFF) CompactSaveErrorGraph<unsigned char>(eg, buf_eg+cur);
	else if(eg->num_nodes<0xFFFF) CompactSaveErrorGraph<unsigned short>(eg, buf_eg+cur);
	else CompactSaveErrorGraph<unsigned int>(eg, buf_eg+cur);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Type>
void CompactLoadErrorGraph(ErrorGraph* eg, char* buf_eg)
{
	int i;
	int cur=0;

	if(eg->match_type==MATCH_SUBSTITUTE_ONLY)
	{
		int i;
		for(i=0;i<eg->num_hamming_scores*NUM_VALID_CHARS;i++)
		{
			Type val;
			memcpy(&val, buf_eg+cur, sizeof(Type)); cur+=sizeof(Type);
			eg->hamming_weights[i]=val;
		}
		return;
	}

	Type j;
	Type node_ind;

	memcpy(&node_ind, buf_eg+cur, sizeof(node_ind)); cur+=sizeof(node_ind); eg->start_node=node_ind;
	memcpy(&node_ind, buf_eg+cur, sizeof(node_ind)); cur+=sizeof(node_ind); eg->end_node=node_ind;
	for(i=0;i<eg->num_nodes;i++)
	{
		PatNode pn; eg->pat_nodes.push_back(pn);
		memcpy(&eg->pat_nodes[i].ch, buf_eg+cur, sizeof(eg->pat_nodes[i].ch)); cur+=sizeof(eg->pat_nodes[i].ch);
		Type sz=0; memcpy(&sz, buf_eg+cur, sizeof(sz)); cur+=sizeof(sz); //eg->pat_nodes[i].adj.size();
		for(j=0;j<sz;j++)
		{
			PatNodeConn pnc; eg->pat_nodes[i].adj.push_back(pnc);
			memcpy(&node_ind, buf_eg+cur, sizeof(node_ind)); cur+=sizeof(node_ind); eg->pat_nodes[i].adj[j].node=node_ind;
			float wacc; memcpy(&wacc, buf_eg+cur, sizeof(wacc)); cur+=sizeof(wacc); eg->pat_nodes[i].adj[j].wacc=wacc;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LoadErrorGraph(ErrorGraph* eg, char* buf_eg)
{
	int cur=0;

	memcpy(&eg->confirmed, buf_eg+cur, sizeof(eg->confirmed)); cur+=sizeof(eg->confirmed);
	if(eg->confirmed) return;

	memcpy(&eg->excluded, buf_eg+cur, sizeof(eg->excluded)); cur+=sizeof(eg->excluded);

	if(eg->match_type==MATCH_SUBSTITUTE_ONLY)
	{
		memcpy(&eg->num_hamming_scores, buf_eg+cur, sizeof(eg->num_hamming_scores)); cur+=sizeof(eg->num_hamming_scores);
		return CompactLoadErrorGraph<float>(eg, buf_eg+cur);
	}

	eg->pat_nodes.clear();

	memcpy(&eg->num_nodes, buf_eg+cur, sizeof(eg->num_nodes)); cur+=sizeof(eg->num_nodes);

	if(eg->num_nodes<0xFF) CompactLoadErrorGraph<unsigned char>(eg, buf_eg+cur);
	else if(eg->num_nodes<0xFFFF) CompactLoadErrorGraph<unsigned short>(eg, buf_eg+cur);
	else CompactLoadErrorGraph<unsigned int>(eg, buf_eg+cur);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void InitErrorGraph(ErrorGraph* eg, bool load_mode, bool save_mode, int max_str_len, int min_overlap, double min_overlap_percentage, double max_overlap_error_rate, double truncate_factor, double min_straight_edge_val, double min_straight_edge_percentage, double min_read_weight, int max_mul_len_matches, int match_type)
{
	eg->max_str_len=max_str_len;
	eg->max_overlap_error_rate=max_overlap_error_rate;
	eg->min_overlap=min_overlap;
	eg->min_overlap_percentage=min_overlap_percentage;
	eg->truncate_factor=truncate_factor;
	eg->min_straight_edge_val=min_straight_edge_val;
	eg->min_straight_edge_percentage=min_straight_edge_percentage;
	eg->min_read_weight=min_read_weight;
	eg->max_mul_len_matches=max_mul_len_matches;
	eg->match_type=match_type;

	eg->load_mode=load_mode;
	eg->save_mode=save_mode;

	if(eg->match_type==MATCH_SUBSTITUTE_ONLY)
	{
		eg->hamming_weights=(double*)AllocateHeap(sizeof(double)*eg->max_str_len*NUM_VALID_CHARS);
	}
	else
	{
		eg->T=(int*)AllocateHeap(3*(eg->max_str_len+1)*sizeof(int));
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DestroyErrorGraph(ErrorGraph* eg)
{
	if(eg->match_type==MATCH_SUBSTITUTE_ONLY)
	{
		FreeHeap(eg->hamming_weights, sizeof(double)*eg->max_str_len*NUM_VALID_CHARS);
	}
	else
	{
		FreeHeap(eg->T, 3*(eg->max_str_len+1)*sizeof(int));
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TruncatePath(ErrorGraph* eg, vector<double>& path_w, int& st_path, int& en_path, int& num_truncated)
{
	int i;

	st_path=0;
	en_path=path_w.size();

	double EPS=1e-4;

	if(eg->truncate_factor>EPS)
	{
		int f=1;

		while(f)
		{
			f=0;
			for(i=st_path;i<en_path && en_path>st_path+2;i++)
			{
				if(path_w[i]<=eg->truncate_factor)
				{
					f=1;

					if(i+1-st_path<en_path-i)
					{
						num_truncated+=i+1-st_path;
						st_path=i+1;
					}
					else
					{
						num_truncated+=en_path-i;
						en_path=i;
					}

					break;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void InitializeMatching(ErrorGraph* eg, char* pat, int pat_len, char* pat_quality)
{
	int i;
	eg->pat_nodes.clear();

	// node 0 is start node (does not contain characters)
	// node pat_len+1 is end node (does not contain characters)
	for(i=0;i<pat_len+2;i++)
	{
		PatNode pn;
		eg->pat_nodes.push_back(pn);

		eg->pat_nodes[i].ch=0;
		if(i>0 && i<pat_len+1) eg->pat_nodes[i].ch=pat[i-1];

		if(i<pat_len+1)
		{
			PatNodeConn pnc;
			pnc.node=i+1;
			pnc.wacc=1.0;
			if(pat_quality && i<pat_len) pnc.wacc=CharToProb[(unsigned char)pat_quality[i]];
			eg->pat_nodes[i].adj.push_back(pnc);
		}
	}

	eg->start_node=0;
	eg->end_node=pat_len+1;
	eg->num_nodes=pat_len+2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetTopologicalOrder(ErrorGraph* eg, int* ordered_nodes)
{
	int i,j;

	for(i=0;i<eg->num_nodes;i++) eg->pat_nodes[i].num_incoming_edges=0;

	for(i=0;i<eg->num_nodes;i++)
	{
		for(j=0;j<(int)eg->pat_nodes[i].adj.size();j++)
		{
			eg->pat_nodes[eg->pat_nodes[i].adj[j].node].num_incoming_edges++;
		}
	}

	int num_ordered_nodes=0;

	for(i=0;i<eg->num_nodes;i++)
	{
		if(eg->pat_nodes[i].num_incoming_edges==0) ordered_nodes[num_ordered_nodes++]=i;
	}

	int num_processed_nodes=0;

	while(num_processed_nodes<num_ordered_nodes)
	{
		int cur_node=ordered_nodes[num_processed_nodes++];

		for(j=0;j<(int)eg->pat_nodes[cur_node].adj.size();j++)
		{
			int next_node=eg->pat_nodes[cur_node].adj[j].node;
			eg->pat_nodes[next_node].num_incoming_edges--;
			if(eg->pat_nodes[next_node].num_incoming_edges==0) ordered_nodes[num_ordered_nodes++]=next_node;
		}
	}

	for(i=0;i<eg->num_nodes;i++)
	{
		if(eg->pat_nodes[i].num_incoming_edges!=0)
		{
			printf("Cycle Detected %d %d\n", i, eg->pat_nodes[i].num_incoming_edges);
			fflush(NULL);
		}
	}

	if(num_ordered_nodes!=eg->num_nodes)
	{
		printf("Unexpected Error NumNode Topology %d %d\n", num_ordered_nodes, eg->num_nodes);
		fflush(NULL);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FinalizeMatching(ErrorGraph* eg, int read_len, char* read_quality, char* result_str, int& result_len, char* result_quality, int& num_truncated)
{
	int i,j;

	int* ordered_nodes=new int[eg->num_nodes];

	GetTopologicalOrder(eg, ordered_nodes);

	double EPS=1e-4;

	for(i=0;i<eg->num_nodes;i++)
	{
		eg->pat_nodes[i].path_weight=-1.0;
		eg->pat_nodes[i].path_prev_node=-1;
	}
	eg->pat_nodes[eg->start_node].path_weight=1;

	for(i=0;i<eg->num_nodes;i++)
	{
		int cur_node=ordered_nodes[i];

		double weight_with_n=0;
		double weight_without_n=0;

		double best_score=-1.0;
		double org_best_score=-1.0;

		int num_adj=eg->pat_nodes[cur_node].adj.size();
		for(j=0;j<num_adj;j++)
		{
			PatNodeConn& pnc=eg->pat_nodes[cur_node].adj[j];
			if(cur_node>=eg->start_node && cur_node<eg->end_node && pnc.node==cur_node+1 && eg->pat_nodes[pnc.node].ch!='N') org_best_score=pnc.wacc;
			weight_with_n+=pnc.wacc;
			if(pnc.wacc>best_score) best_score=pnc.wacc;
			if(eg->pat_nodes[pnc.node].ch!='N') weight_without_n+=pnc.wacc;
			pnc.prob=0;
		}

		if(org_best_score>=eg->min_straight_edge_val || best_score<EPS || org_best_score/best_score>=eg->min_straight_edge_percentage)
		{
			if(num_adj) eg->pat_nodes[cur_node].adj[0].prob=1;
		}
		else if(weight_with_n>EPS)
		{
			for(j=0;j<num_adj;j++)
			{
				PatNodeConn& pnc=eg->pat_nodes[cur_node].adj[j];
				if(weight_without_n>EPS){if(eg->pat_nodes[pnc.node].ch!='N') pnc.prob=pnc.wacc/weight_without_n;} else
				pnc.prob=pnc.wacc/weight_with_n;
			}
		}
	}

	long double best_weight=-1.0;

	for(i=0;i<eg->num_nodes;i++)
	{
		int cur_node=ordered_nodes[i];

		for(j=0;j<(int)eg->pat_nodes[cur_node].adj.size();j++)
		{
			PatNodeConn& pnc=eg->pat_nodes[cur_node].adj[j];

			long double new_path_weight=eg->pat_nodes[cur_node].path_weight*pnc.prob;

			if(new_path_weight>eg->pat_nodes[pnc.node].path_weight)
			{
				eg->pat_nodes[pnc.node].path_weight=new_path_weight;
				eg->pat_nodes[pnc.node].path_prev_node=cur_node;

				if(pnc.node==eg->end_node && new_path_weight>best_weight)
				{
					best_weight=new_path_weight;
				}
			}
		}
	}

	vector<char>	path_ch; // the size of path_ch is less than the size of path_w by 1
	vector<double>	path_w;
	vector<int>		path_quality_ind;
	vector<char>	path_quality;

	int cur_node=eg->end_node;

	while(cur_node!=eg->start_node)
	{
		int prev_node=eg->pat_nodes[cur_node].path_prev_node;

		for(i=0;prev_node>=0 && i<(int)eg->pat_nodes[prev_node].adj.size();i++)
		{
			if(eg->pat_nodes[prev_node].adj[i].node==cur_node)
			{
				char ch=eg->pat_nodes[cur_node].ch; if(ch) path_ch.push_back(ch);
				path_w.push_back(eg->pat_nodes[prev_node].adj[i].wacc);
				if(read_quality)
				{
					int q_ind=-1; if(cur_node>eg->start_node && cur_node<eg->end_node) q_ind=cur_node-1;
					if(ch) path_quality_ind.push_back(q_ind); // if(ch) to make sure path_quality_ind.size()==path_ch.size()
				}
				break;
			}
		}

		cur_node=prev_node;
	}

	if(read_quality)
	{
		int pqsz=(int)path_quality_ind.size();
		for(i=0;i<pqsz;i++)
		{
			if(path_quality_ind[i]<0)
			{
				for(j=i;j>=0;j--) if(path_quality_ind[j]>=0) break;
				int en=read_len-1; if(j>=0) en=path_quality_ind[j];

				for(j=i;j<pqsz;j++) if(path_quality_ind[j]>=0) break;
				int st=0; if(j<pqsz) st=path_quality_ind[j];

				int sc=0; for(j=st;j<=en;j++) sc+=CharToQual[(unsigned char)read_quality[j]]; sc/=(en-st+1); path_quality.push_back(QualToChar(sc));
			}
			else
			{
				path_quality.push_back(read_quality[path_quality_ind[i]]);
			}
		}
	}

	int st_path, en_path;

	TruncatePath(eg, path_w, st_path, en_path, num_truncated);

	result_len=en_path-st_path-1;
	if(result_len>eg->max_str_len*2){result_len=0;result_str[0]=0;read_quality[0]=0;return;}
	for(i=result_len-1;i>=0;i--)
	{
		if(i+st_path>=(int)path_ch.size() || (read_quality && path_ch.size()!=path_quality.size())) {printf("Unexpected Error FinalizeMatching\n"); fflush(NULL);}
		result_str[result_len-1-i]=path_ch[i+st_path];
		if(read_quality) result_quality[result_len-1-i]=path_quality[i+st_path];
	}
	result_str[result_len]=0;
	if(read_quality) result_quality[result_len]=0;

	delete[] ordered_nodes;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool IsConfirmed(ErrorGraph* eg, char* pat, int pat_len)
{
	int j;
	for(j=0;j<=pat_len;j++) if(pat[j]=='N' || eg->pat_nodes[j].adj[0].wacc<eg->min_straight_edge_val) break;
	return (j==pat_len+1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int NewNode(ErrorGraph* eg)
{
	eg->num_nodes++;
	if((int)eg->pat_nodes.size()<eg->num_nodes) {PatNode pn; eg->pat_nodes.push_back(pn);}
	return eg->num_nodes-1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int LinkNodes(ErrorGraph* eg, int prev_node, char ch, int cur_node, double w) //, int org_start, double w) // ch = char of destination node
{
	int i;
	if(cur_node==0) {printf("NP Link to Start Node\n"); fflush(NULL);}
	if(prev_node==cur_node) printf("CYCLE %d %d *------------*\n", prev_node, cur_node);

	if(cur_node==-1)
	{
		for(i=0;i<(int)eg->pat_nodes[prev_node].adj.size();i++)
		{
			PatNodeConn& pnc=eg->pat_nodes[prev_node].adj[i];

			if(pnc.node>eg->end_node && eg->pat_nodes[pnc.node].ch==ch)
			{
				pnc.wacc+=w;
				cur_node=pnc.node;
				break;
			}
		}

		if(cur_node==-1)
		{
			cur_node=NewNode(eg);

			PatNodeConn pnc;
			pnc.node=cur_node;
			pnc.wacc=w;

			eg->pat_nodes[cur_node].ch=ch;
			eg->pat_nodes[prev_node].adj.push_back(pnc);
		}
	}
	else
	{
		for(i=0;i<(int)eg->pat_nodes[prev_node].adj.size();i++)
		{
			PatNodeConn& pnc=eg->pat_nodes[prev_node].adj[i];

			if(pnc.node==cur_node)
			{
				if(eg->pat_nodes[pnc.node].ch!=ch) {printf("ERRRRRROR\n"); fflush(NULL);}
				pnc.wacc+=w;
				break;
			}
		}

		if(i==(int)eg->pat_nodes[prev_node].adj.size())
		{
			PatNodeConn pnc;
			pnc.node=cur_node;
			pnc.wacc=w;

			eg->pat_nodes[prev_node].adj.push_back(pnc);
		}
	}

	return cur_node;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct LinkOperation
{
	int a_st, a_en; // a_en=after end
	int type;
	vector<char> chs;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int AddNormalizedMatchToGraph(ErrorGraph* eg, char* a, int na, int first_a, int last_a, vector<ShortNode>& normalized, double w, char* qual)
{
	int i,j;
	int cur_a=first_a;

	int start_node=-1;

	vector<LinkOperation> v;

	for(i=0;i<(int)normalized.size();i++)
	{
		ShortNode node=normalized[i];
		int ca=node.ia;
		char ch=node.ch;

		if(cur_a<ca)
		{
			LinkOperation op;
			op.a_st=cur_a;
			op.a_en=ca;
			op.type=0;
			v.push_back(op);
			cur_a=ca;
			i--;
		}
		else
		{
			LinkOperation op;
			op.a_st=ca;
			op.a_en=ca;
			op.type=1;
			if(node.type==ED_DELETE) {op.a_en++; cur_a++;}
			else op.chs.push_back(ch);

			i++;
			for(;i<(int)normalized.size();i++)
			{
				if(normalized[i].ia>cur_a) break;
				if(normalized[i].type==ED_DELETE) {op.a_en++; cur_a++;}
				else op.chs.push_back(normalized[i].ch);
			}
			i--;
			v.push_back(op);
		}
	}

	if(cur_a<last_a)
	{
		LinkOperation op;
		op.a_st=cur_a;
		op.a_en=last_a;
		op.type=0;
		v.push_back(op);
	}

	int prev_node=first_a;

	if(first_a==0)
	{
		prev_node=eg->start_node;
		start_node=prev_node;
	}
	else
	{
		prev_node=-1;
	}

	int q=0;

	for(i=0;i<(int)v.size();i++)
	{
		if(v[i].type==0)
		{
			for(j=v[i].a_st;j<v[i].a_en;j++)
			{
				if(prev_node==-1) {prev_node=j+1; start_node=prev_node;} //InsertStartNode(eg, prev_node);}
				else
				{
					double cw=w; if(qual) cw*=CharToProb[(unsigned char)qual[q]];
					prev_node=LinkNodes(eg, prev_node, a[j], j+1, cw);
				}
				q++;
			}
		}
		else
		{
			for(j=0;j<(int)v[i].chs.size();j++)
			{
				if(prev_node==-1)
				{
					prev_node=NewNode(eg);
					eg->pat_nodes[prev_node].ch=v[i].chs[j];
					start_node=prev_node;
				}
				else
				{
					double cw=w; if(qual) cw*=CharToProb[(unsigned char)qual[q]];
					prev_node=LinkNodes(eg, prev_node, v[i].chs[j], -1, w);
				}
				q++;
			}
		}
	}

	if(last_a==na)
	{
		LinkNodes(eg, prev_node, 0, eg->end_node, w);
	}

	return start_node;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double GetWeight(double min_read_weight, int num_errors, int cur_allowed_overlap_errors)
{
	double cur_whole_weight=1.0-(1-min_read_weight)*sqrt((double)num_errors/cur_allowed_overlap_errors);
	return cur_whole_weight;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ProcessLocation(ErrorGraph* eg, char* a, int na, char* b, int nb, int first_a, int first_b, int last_a, int last_b, int edit, vector<ShortNode>& v, char* qual)
{
	int i;

	int cur_overlap_size=last_a-first_a;
	int cur_allowed_overlap_errors=eg->max_overlap_error_rate*cur_overlap_size;

	int cur_min_overlap=eg->min_overlap_percentage*na;
	if(cur_min_overlap<eg->min_overlap) cur_min_overlap=eg->min_overlap;

	if(cur_overlap_size-edit<cur_min_overlap || edit>cur_allowed_overlap_errors)
	{
		return -1;
	}

	double cur_whole_weight=GetWeight(eg->min_read_weight, edit, cur_allowed_overlap_errors);

	vector<ShortNode> normalized;

	for(i=0;i<(int)v.size();i++)
	{
		ShortNode& node=v[i];

		if(node.type==ED_SUBSTITUTE)
		{
			ShortNode sn;
			sn.ia=node.ia; sn.ch=node.ch; sn.type=ED_INSERT; normalized.push_back(sn);
			sn.ia=node.ia; sn.ch=a[node.ia]; sn.type=ED_DELETE; normalized.push_back(sn);
		}
		else
		{
			if(node.type==ED_DELETE) node.ch=a[node.ia];
			normalized.push_back(node);
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int last_val=na;

	for(i=(int)normalized.size()-1;i>=0;i--)
	{
		if(normalized[i].type==ED_DELETE)
		{
			int ia=normalized[i].ia;
			int ch=normalized[i].ch;

			ia++; while(ia<last_val && a[ia]==ch) ia++; ia--;
			normalized[i].ia=ia;
			last_val=ia;

			int j=i;

			for(j=i;j<(int)normalized.size()-1;j++)
			{
				if(normalized[j+1].type==ED_DELETE)
					break;

				if(normalized[j].ia>=normalized[j+1].ia)
				{
					ShortNode temp=normalized[j]; normalized[j]=normalized[j+1]; normalized[j+1]=temp;
					normalized[j].ia--;
				}
				else break;
			}
		}
		else if(normalized[i].type==ED_INSERT)
		{
			int ia=normalized[i].ia;
			int ch=normalized[i].ch;

			if(ia<na && a[ia]==ch)
			{
				ia++; while(ia<na && a[ia]==ch) ia++;
				normalized[i].ia=ia;
			}

			int j=i;

			for(j=i;j<(int)normalized.size()-1;j++)
			{
				if(normalized[j+1].type==ED_DELETE)
				{
					if(normalized[j].ia>=normalized[j+1].ia)
					{
						if(normalized[j].ch==normalized[j+1].ch)
						{
							normalized.erase(normalized.begin()+j+1);
							normalized.erase(normalized.begin()+j);
							break;
						}

						ShortNode temp=normalized[j]; normalized[j]=normalized[j+1]; normalized[j+1]=temp;
						if(normalized[j+1].ia<na) normalized[j+1].ia++;

						int ia=normalized[j+1].ia;
						char ch=normalized[j+1].ch;

						if(ia<na && a[ia]==ch)
						{
							ia++; while(ia<na && a[ia]==ch) ia++;
							normalized[j+1].ia=ia;
						}
					}
					else break;
				}
				else
				{
					if(normalized[j].ia>normalized[j+1].ia)
					{
						ShortNode temp=normalized[j]; normalized[j]=normalized[j+1]; normalized[j+1]=temp;
						normalized[j].ia++;
					}
					else break;
				}
			}
		}
	}

	int start_node=AddNormalizedMatchToGraph(eg, a, na, first_a, last_a, normalized, cur_whole_weight, qual);

	return start_node;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FixReadUsingGraph(ErrorGraph* eg, unsigned short* kmer_first, int kmer_size, char* pat, int pat_len, char* pat_quality, vector<PossibleMatch>& possible_matches, char* result_str, int& result_len, char* result_quality, int& num_truncated)
{
	int i;

	if(eg->match_type==MATCH_SUBSTITUTE_ONLY)
	{
		int j,k;

		if(!eg->load_mode)
		{
			eg->confirmed=false;
			eg->excluded=false;

			eg->num_hamming_scores=pat_len;

			for(j=0;j<pat_len;j++)
			{
				for(k=0;k<NUM_VALID_CHARS;k++) eg->hamming_weights[j*NUM_VALID_CHARS+k]=0;

				int cti=CharToInt5Val[(int)pat[j]];
				if(cti<NUM_VALID_CHARS)
				{
					eg->hamming_weights[j*NUM_VALID_CHARS+cti]=1.0;
					if(pat_quality) eg->hamming_weights[j*NUM_VALID_CHARS+cti]=CharToProb[(unsigned char)pat_quality[j]];
				}
			}
		}

		int num=0;

		int cur_min_overlap=eg->min_overlap_percentage*pat_len;
		if(cur_min_overlap<eg->min_overlap) cur_min_overlap=eg->min_overlap;

		for(i=0;i<(int)possible_matches.size() && !eg->confirmed;i++)
		{
			for(j=0;j<pat_len;j++)
			{
				int cti=CharToInt5Val[(int)pat[j]];
				if(cti>=NUM_VALID_CHARS || eg->hamming_weights[j*NUM_VALID_CHARS+cti]<eg->min_straight_edge_val) break;
			}
			if(j==pat_len) {eg->confirmed=true; break;}

			if(!eg->excluded && possible_matches[i].str_len==pat_len && memcmp(pat, possible_matches[i].str, pat_len)==0 &&
					(!pat_quality || !possible_matches[i].qual || memcmp(pat_quality, possible_matches[i].qual, pat_len)==0))
			{
				eg->excluded=true;
				continue;
			}

			int max_overlap_errors=possible_matches[i].max_overlap*eg->max_overlap_error_rate;

			int min_hamming_dist=max_overlap_errors+1;
			double cur_whole_weight=0;

			vector<int> shifts;

			// b starts at a[shift]

			for(k=possible_matches[i].shift1;k<=possible_matches[i].shift2;k++)
			{
				int hd=0;
				int cur_overlap_size=0;

				for(j=0;j<pat_len && hd<min_hamming_dist;j++)
				{
					int ib=j-k;
					if(ib>=0 && ib<possible_matches[i].str_len)
					{
						cur_overlap_size++;
						hd+=(pat[j]!=possible_matches[i].str[ib]);
					}
				}

				int cur_allowed_overlap_errors=eg->max_overlap_error_rate*cur_overlap_size;

				if(cur_overlap_size-hd<cur_min_overlap || hd>cur_allowed_overlap_errors)
					continue;

				if(hd<min_hamming_dist)
				{
					min_hamming_dist=hd;

					shifts.clear();
					shifts.push_back(k);

					double weight=GetWeight(eg->min_read_weight, hd, cur_allowed_overlap_errors);
					cur_whole_weight=weight;
				}
				else if(hd==min_hamming_dist)
				{
					shifts.push_back(k);

					double weight=GetWeight(eg->min_read_weight, hd, cur_allowed_overlap_errors);
					if(weight>cur_whole_weight) cur_whole_weight=weight;
				}
			}

			for(k=0;k<(int)shifts.size();k++)
			{
				int cur_shift=shifts[k];
				double w=cur_whole_weight/shifts.size();

				for(j=0;j<pat_len;j++)
				{
					int ib=j-cur_shift;
					if(ib>=0 && ib<possible_matches[i].str_len)
					{
						int cti=CharToInt5Val[(int)possible_matches[i].str[ib]];
						if(cti<NUM_VALID_CHARS)
						{
							double cw=w; if(possible_matches[i].qual) cw*=CharToProb[(unsigned char)possible_matches[i].qual[ib]];
							eg->hamming_weights[j*NUM_VALID_CHARS+cti]+=cw;
						}
					}
				}

			}

			if(shifts.size()) num++;
		}

		if(!eg->save_mode)
		{
			if(eg->confirmed)
			{
				result_len=pat_len;
				for(i=0;i<result_len;i++) result_str[i]=pat[i];
				result_str[result_len]=0;
				if(pat_quality) strcpy(result_quality, pat_quality);
			}
			else
			{
				vector<char> path_ch;
				vector<double> path_w;

				double EPS=1e-4;

				for(j=0;j<pat_len;j++)
				{
					int cti=CharToInt5Val[(int)pat[j]];
					double best_w=0;
					int best_k=cti;

					for(k=0;k<NUM_VALID_CHARS;k++)
					{
						double w=eg->hamming_weights[j*NUM_VALID_CHARS+k];
						if(k==cti && w>=eg->min_straight_edge_val) {best_k=cti; best_w=w; break;}
						if(w>best_w){best_k=k; best_w=w;}
					}

					if(cti<NUM_VALID_CHARS)
					{
						double org_w=eg->hamming_weights[j*NUM_VALID_CHARS+cti];
						if(best_w<EPS || org_w/best_w>=eg->min_straight_edge_percentage) {best_k=cti; best_w=org_w;}
					}

					path_ch.push_back(IntToChar[best_k]);
					path_w.push_back(best_w);
				}

				int st_path, en_path;
				TruncatePath(eg, path_w, st_path, en_path, num_truncated);

				result_len=en_path-st_path;
				if(result_len>eg->max_str_len*2){result_len=0;result_str[0]=0;return;}
				for(i=0;i<result_len;i++) result_str[i]=path_ch[i+st_path];
				result_str[result_len]=0;
				if(pat_quality) {for(i=0;i<result_len;i++) result_quality[i]=pat_quality[i+st_path]; result_quality[result_len]=0;}
			}
		}

		return;
	}

	if(!eg->load_mode)
	{
		eg->confirmed=false;
		eg->excluded=false;

		InitializeMatching(eg, pat, pat_len, pat_quality);
	}

	char* a=pat;
	int na=pat_len;

	int j;
	int* kmer_next=new int[na];
	for(j=0;j<na;j++) kmer_next[j]=-1;

	unsigned int kmer=GetKmer(a, na, kmer_size);

	for(j=0;j<na-kmer_size+1;j++)
	{
		if(j)
		{
			int new_char=0; if(j+kmer_size-1<na) new_char=CharToInt[(int)a[j+kmer_size-1]];
			kmer=(kmer>>BITS_PER_SYMBOL)|(new_char<<((kmer_size-1)*BITS_PER_SYMBOL));
		}

		kmer_next[j]=kmer_first[kmer]; kmer_next[j]--;
		kmer_first[kmer]=j+1;
	}

	for(i=0;i<(int)possible_matches.size() && !eg->confirmed;i++)
	{
		eg->confirmed=IsConfirmed(eg, pat, pat_len);
		if(eg->confirmed) break;

		if(!eg->excluded && possible_matches[i].str_len==pat_len && memcmp(pat, possible_matches[i].str, pat_len)==0 &&
				(!pat_quality || !possible_matches[i].qual || memcmp(pat_quality, possible_matches[i].qual, pat_len)==0))
		{
			eg->excluded=true;
			continue;
		}

		int max_overlap_errors=possible_matches[i].max_overlap*eg->max_overlap_error_rate;

		char* b=possible_matches[i].str;
		int nb=possible_matches[i].str_len;

		int first_a=-1, first_b=-1;
		int last_a=-1, last_b=-1;
		int edit=-1;
		vector<ShortNode> v;

		int best_sa=0, best_sb=0, best_len=0;

		kmer=GetKmer(b, nb, kmer_size);

		for(j=0;j<nb-kmer_size+1;j++)
		{
			if(j)
			{
				int new_char=0; if(j+kmer_size-1<nb) new_char=CharToInt[(int)b[j+kmer_size-1]];
				kmer=(kmer>>BITS_PER_SYMBOL)|(new_char<<((kmer_size-1)*BITS_PER_SYMBOL));
			}

			int cur=kmer_first[kmer]; cur--;

			while(cur>=0)
			{
				if(j==0 || cur==0 || b[j-1]!=a[cur-1])
				{
					int eb=j+kmer_size, ea=cur+kmer_size;
					while(eb<nb && ea<na && b[eb]==a[ea]) {eb++; ea++;}
					int len=eb-j;

					if(len>best_len) {best_len=len; best_sa=cur, best_sb=j;}
				}
				cur=kmer_next[cur];
			}
		}

		if(best_len==0)
		{
			printf("Bad Loc\n");
			fflush(NULL);
			continue;
		}

		int left_a=best_sa;
		int right_a=best_sa+best_len;

		int left_b=best_sb;
		int right_b=best_sb+best_len;

		int edit_right=0;
		if(na-right_a && nb-right_b)
		{
			edit_right=ApproxEditDist(eg->T, a+right_a, na-right_a, b+right_b, nb-right_b, -ED_INF/3, ED_INF/3, max_overlap_errors, MATCH_EDIT_DIST, false, true, last_a, last_b);
			last_a+=right_a; last_b+=right_b;
		}
		else
		{
			last_a=right_a;
			last_b=right_b;
		}

		int edit_left=0;
		if(left_a && left_b)
		{
			edit_left=ApproxEditDistBackward(eg->T, a, left_a, b, left_b, -ED_INF/3, ED_INF/3, max_overlap_errors, MATCH_EDIT_DIST, false, true, first_a, first_b);
		}
		else
		{
			first_a=left_a;
			first_b=left_b;
		}

		if(edit_right>=0 && edit_left>=0)
		{
			if(edit_left && left_a-first_a && left_b-first_b) ApproxGetOptimalAlignment(eg->T, first_a, a+first_a, left_a-first_a, b+first_b, left_b-first_b, edit_left, MATCH_EDIT_DIST, v);
			if(edit_right && last_a-right_a && last_b-right_b) ApproxGetOptimalAlignment(eg->T, right_a, a+right_a, last_a-right_a, b+right_b, last_b-right_b, edit_right, MATCH_EDIT_DIST, v);

			edit=edit_right+edit_left;

			/////////////////////////////////////////////////////////////////////

			int cur_overlap_size=last_a-first_a; if(last_b-first_b<cur_overlap_size) cur_overlap_size=last_b-first_b;
			int cur_allowed_overlap_errors=eg->max_overlap_error_rate*cur_overlap_size;

			int cur_min_overlap=eg->min_overlap_percentage*na;
			if(cur_min_overlap<eg->min_overlap) cur_min_overlap=eg->min_overlap;

			if(cur_overlap_size-edit>=cur_min_overlap && edit<=cur_allowed_overlap_errors)
			{
				ProcessLocation(eg, a, na, b, nb, first_a, first_b, last_a, last_b, edit, v, possible_matches[i].qual?(possible_matches[i].qual+first_b):0);
			}
		}
	}

	kmer=GetKmer(a, na, kmer_size);

	for(j=0;j<na-kmer_size+1;j++)
	{
		if(j)
		{
			int new_char=0; if(j+kmer_size-1<na) new_char=CharToInt[(int)a[j+kmer_size-1]];
			kmer=(kmer>>BITS_PER_SYMBOL)|(new_char<<((kmer_size-1)*BITS_PER_SYMBOL));
		}

		kmer_first[kmer]=0;
	}

	delete[] kmer_next;

	if(!eg->save_mode)
	{
		if(eg->confirmed)
		{
			result_len=pat_len;
			for(i=0;i<result_len;i++) result_str[i]=pat[i];
			result_str[result_len]=0;
			if(pat_quality) strcpy(result_quality, pat_quality);
		}
		else
		{
			FinalizeMatching(eg, pat_len, pat_quality, result_str, result_len, result_quality, num_truncated);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

