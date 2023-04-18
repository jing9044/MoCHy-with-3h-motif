#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <chrono>

#include "read_data.cpp"
#include "motif_id.cpp"
#include "3h-motif_id_ab1.cpp"

using namespace std;

struct hyperwedge{
	int a, b, C_ab;	
};

inline long long convert_id(int hyperedge_a, int hyperedge_b){
	return hyperedge_a * (1LL << 31) + hyperedge_b;
}

int main(int argc, char *argv[])
{
	clock_t start;
	clock_t run_start;
	int progress;
	
	int sampling_size = stoi(argv[1]);
	string threshold_type = argv[2];

	string graphFile = "dblp_graph.txt";

	cout << "Sampling size: " << sampling_size << endl << endl;

	// Read data
	start = clock();
	vector< vector<int> > node2hyperedge;
	vector< vector<int> > hyperedge2node;
	vector< unordered_set<int> > hyperedge2node_set;
	read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set);

	int V = node2hyperedge.size(), E = hyperedge2node.size();
	cout << "# of nodes: " << V << '\n';
	cout << "# of hyperedges: " << E << '\n';
	cout << "Reading data done: "
		<< (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;

		
	// Adjacency list construction
	start = clock(); run_start = clock();
	hyperedge2node.resize(E); hyperedge2node_set.resize(E);
	vector< vector< pair<int, int> > > hyperedge_adj;
	vector< unordered_map<int, int> > hyperedge_inter;
	hyperedge_adj.resize(E);
	hyperedge_inter.resize(E);
	vector< hyperwedge > W;
	vector<long long> upd_time(E, -1LL);
	
	for(int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++){
		long long l_hyperedge_a = (long long)hyperedge_a;
		for(const int &node: hyperedge2node[hyperedge_a]){
			for(const int &hyperedge_b: node2hyperedge[node]){
				if(hyperedge_b == hyperedge_a) continue;
				if((upd_time[hyperedge_b] >> 31) ^ hyperedge_a){
					upd_time[hyperedge_b] = (l_hyperedge_a << 31) + (long long)hyperedge_adj[hyperedge_b].size();
					hyperedge_adj[hyperedge_b].push_back({hyperedge_a, 0});
				}
				hyperedge_adj[hyperedge_b][(int)(upd_time[hyperedge_b] & 0x7FFFFFFFLL)].second++;
			}
		}
	}

	for (int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++){
		int deg_a = hyperedge_adj[hyperedge_a].size();
		hyperedge_inter[hyperedge_a].rehash(deg_a);
		for (int i = 0; i < deg_a; i++){
			int hyperedge_b = hyperedge_adj[hyperedge_a][i].first;
			int C_ab = hyperedge_adj[hyperedge_a][i].second;
			if (hyperedge_a < hyperedge_b) W.push_back(hyperwedge{hyperedge_a, hyperedge_b, C_ab});
			hyperedge_inter[hyperedge_a].insert({hyperedge_b, C_ab});
		}
	}

	cout << "# of hyperwedges: " << W.size() << "\n";
	cout << "Adjacency list construction done: "
		<< (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;


	// Motif counting via hyperwedge smapling
	start = clock();

	vector<long long> h_motif(30, 0);
	
	vector< vector<long long> > th_motif(27, vector<long long>(40, 0));
    vector<int> th_motif_kind = {0, 6, 8, 8, 12, 16, 24, 6, 12, 16, 32, 20, 40, 8, 24, 40, 40, 3, 6, 8, 16, 10, 20, 4, 12, 20, 20};
	
	vector<int> intersection(E, 0);
	std::fill(upd_time.begin(), upd_time.end(), -1LL);

	random_device rd;
	mt19937 gen(2020);
	uniform_int_distribution<> dist(0, ((int)W.size())-1);

	for (int sample = 0; sample < sampling_size; sample++){
		if (sample % 100000 == 0)
			cout << "Sampling: " << sample << " / " << sampling_size << endl;
		
		int sample_index = dist(gen);
		int hyperedge_a = W[sample_index].a;
		int hyperedge_b = W[sample_index].b;
		int C_ab = W[sample_index].C_ab;
		
		int size_a = (int)hyperedge2node[hyperedge_a].size();
		int size_b = (int)hyperedge2node[hyperedge_b].size();
		int deg_a = (int)hyperedge_adj[hyperedge_a].size();
		int deg_b = (int)hyperedge_adj[hyperedge_b].size();

		upd_time[hyperedge_a] = upd_time[hyperedge_b] = sample;

		int min_ab = hyperedge_a, max_ab = hyperedge_b;
		if (size_a > size_b) min_ab = hyperedge_b, max_ab = hyperedge_a;
		const auto &nodes = hyperedge2node_set[max_ab]; auto it_end = nodes.end(); int cnt = 0;
		for (const int &node: hyperedge2node[min_ab]){ if(nodes.find(node) != it_end) intersection[cnt++] = node;}

		for (int i = 0; i < deg_b; i++){
			int hyperedge_c = hyperedge_adj[hyperedge_b][i].first, C_bc = hyperedge_adj[hyperedge_b][i].second;
			if (upd_time[hyperedge_c] ^ sample){
				upd_time[hyperedge_c] = sample;

				int size_c = (int)hyperedge2node[hyperedge_c].size();
				int C_ca = 0, g_abc = 0;
				C_ca = hyperedge_inter[hyperedge_a][hyperedge_c];
				const auto &nodes = hyperedge2node_set[hyperedge_c]; auto it_end = nodes.end();
				for (int k = 0; k < C_ab; k++){ if(nodes.find(intersection[k]) != it_end) g_abc++; }
				int h_motif_index;
				pair <int, int> th_motif_index;
				if (threshold_type == "none") {
					h_motif_index = get_motif_index_new(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
					h_motif[h_motif_index]++;
				}
				else {
					th_motif_index = get_motif_index_ab1(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
					th_motif[th_motif_index.first][th_motif_index.second]++;
				}
			}
		}

		for (int i = 0; i < deg_a; i++){
			int hyperedge_c = hyperedge_adj[hyperedge_a][i].first, C_ca = hyperedge_adj[hyperedge_a][i].second;
			if (upd_time[hyperedge_c] ^ sample){
				upd_time[hyperedge_c] = sample;
				
				int size_c = (int)hyperedge2node[hyperedge_c].size();
				int C_bc = 0, g_abc = 0;

				int h_motif_index;
				pair <int, int> th_motif_index;
				if (threshold_type == "none") {
					h_motif_index = get_motif_index_new(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
					h_motif[h_motif_index]++;
				}
				else {
					th_motif_index = get_motif_index_ab1(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
					th_motif[th_motif_index.first][th_motif_index.second]++;
				}
			}
		}
	}

	int index = 0;
	if (threshold_type == "none") {
		vector<long double> h_motif_final(30, 0);
		for (int i = 0; i < 30; i++){
			h_motif_final[i] = (long double)h_motif[i];
			h_motif_final[i] *= (long double)W.size() / sampling_size;
			if (20 <= i && i <= 25)
				h_motif_final[i] /= 2.0;
			else
				h_motif_final[i] /= 3.0;
			if (i == 0 || i == 1 || i == 4 || i == 6) continue;
			cout << "h-motif " << ++index << ": " << h_motif_final[i] << endl;
		}
	}
	else {
		vector < vector<long double> > th_motif_final(27, vector<long double>(40, 0));
		for (int i = 1; i < 27; i++){
			for (int j = 0; j < th_motif_kind[i]; j++) {
				th_motif_final[i][j] = (long double)th_motif[i][j];
				th_motif_final[i][j] *= (long double)W.size() / sampling_size;
				if (17 <= i && i <= 22) 
					th_motif_final[i][j] /= 2.0;
				else 
					th_motif_final[i][j] /= 3.0;
				cout << "3h-motif " << ++index << ": " << th_motif_final[i][j] << endl;
			}
		}
	}

	double runtime = (double)(clock() - run_start) / CLOCKS_PER_SEC;

	cout << "\nHypergraph motif counting done: "
		<< (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "Total runtime: " << runtime << endl;
	cout << "-----------------------------------------" << endl << endl;

	W.clear();
	node2hyperedge.clear();
	hyperedge2node.clear();
	hyperedge_adj.clear();
	hyperedge_inter.clear();

	return 0;
}
