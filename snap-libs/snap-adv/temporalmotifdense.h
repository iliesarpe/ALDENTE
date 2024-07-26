#ifndef snap_temporalmotifdense_h
#define snap_temporalmotifdense_h

#include "Snap.h"

#include "temporalmotiftypes.h"

#include "datastructures.h"

#include <vector>
#include <map>
#include <stack>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <thread>
#include <future>
#include <fstream>
#include <sstream>
#include <iterator>
#include <chrono>
#include <boost/heap/fibonacci_heap.hpp>



#define VERBOSE 100
#define THREADS 10
#define SINGLETHREAD 1 
// TODO: no global namespaces

//using namespace std;
using timestamp_t = long long int;

using DataStructures::TEdge;


struct heap_data
{
	long long int node;
	double count;

	heap_data(long long int nodeval, double countval):node(nodeval), count(countval)
	{}
};

struct compare_node
{
	bool operator()(const heap_data& n1, const heap_data& n2) const
	{
		if(n1.count == n2.count)
			return n1.node > n2.node;
		return n1.count > n2.count;
	}
};



class MyStructure{
	private:
		int nodes;
		std::vector<std::vector<long long int>> instance_node_map; // i: instance, vect[i] nodes in instance i
		std::vector<std::set<long long int>> node_instance_map; // i: node, set[i] instances containg node_i
		std::vector<double> node_counts; // Keep for each node the count C_v
		std::vector<bool> nodes_alive; // Keep the nodes that are alive when peeling them
		std::vector<boost::heap::fibonacci_heap<heap_data, boost::heap::compare<compare_node>>::handle_type> node_handles_heap; // Keep for each node "pointer" in heap
		boost::heap::fibonacci_heap<heap_data, boost::heap::compare<compare_node>> heapcounts; // heap to store counts of nodes during peeling
		std::unordered_set<long long int> solution;
		std::unordered_map<std::string, double> subgraphCount;
		std::map<std::string, double> staticCounts;
		std::vector<std::unordered_set<std::string>> nodeSubgraphs;
		std::unordered_map<long long int, long long int> handleMap;
		int alive; // <-- |S|
		double solutionInsts; // k*C_S <-- /|S|
		int iteration; // iteration of randomized batch peeling algorithm
		std::vector<std::vector<double>> familyF; //familyF[i][j] keeps the estimated count of node i on j-th sample
		int weight; //controls weighting function, 1: tau_count, 2: tau_decay
		double decay;
		double min_weight;
		std::vector<long long int> optimalSolution;
	public:
		MyStructure(int node, int weight=1, double decay=1) {
			this->nodes = node;
			for(int i=0; i < node; i++)
			{
				std::set<long long int> set_v;
				node_instance_map.push_back(set_v);
				solution.insert(i);
				std::unordered_set<std::string> uset_v;
				nodeSubgraphs.push_back(uset_v);
			}
			solutionInsts = 0.0;
			node_counts.resize(node, 0);
			nodes_alive.resize(node, true);
			alive = node;
			iteration = 0;
			familyF.clear();
			optimalSolution.clear();
			this->weight=weight;
			this->decay=decay;
			this->min_weight=1.0;
		}
		const std::vector<long long int>& getOptimalSolution()
		{
			return std::cref(optimalSolution);
		}
		void addOptimalNode(long long int node)
		{
			optimalSolution.push_back(node);
		}
		double getDecay()
		{
			return decay;
		}
		void setMinWeight(double w)
		{
			this->min_weight = std::min(min_weight, w);
		}
		double getMinWeight()
		{
			return min_weight;
		}
		void setDecay(double decay)
		{
			this->decay=decay;
			this->weight=2;
			return;
		}
		int getWeight()
		{
			return weight;
		}
		void setNextIteration(int samples)
		{
			iteration = 0;
			familyF.clear();
			for(int i=0; i< nodes; i++)
			{
				std::vector<double> ests(samples, 0.0);
				familyF.push_back(ests);	
			}
			return;
		}
		void increaseIteration()
		{
			iteration++;
			return;
		}
		int getIteration()
		{
			return iteration;
		}
		void updateFamily(long long node, double value, int type=1)
		{
			if(type==1)
			{
				familyF[node][iteration] += value;	
			}
			else if(type==2)
			{
				for(int i=0; i < iteration; i++)
					familyF[node][i] /= value;
			}
			return;
		}
		double getVarianceNode(long long node, int samples)
		{
			double mean = node_counts[node];
			double variance = 0.0;
			for(auto el: familyF[node])
				variance += pow((el-mean), 2);
			variance /= (1.0*samples);
			return variance;
		}
		//returns Instance ID
		long long int addInstance()
		{
			std::vector<long long int> newinst;
			instance_node_map.push_back(newinst);	
			long long int idx_inst = static_cast<long long int>(instance_node_map.size()) -1;
			return idx_inst;
		}
		void addInstance(std::set<long long int>& nodesinst, double weight=1.0)
		{
			//Encode static subgraph by min ids on nodes, update node counts
			std::string key = "";
			bool first = true;
			for(auto& node : nodesinst)
			{
				node_counts[node] +=weight;
				solutionInsts += weight;
				if(first)
				{
					key = key + std::to_string(node);
					first = false;
				}
				else
					key = key + " " + std::to_string(node);
			}		
			for(auto& node : nodesinst)
			{
				nodeSubgraphs[node].insert(key);
			}
			// Keep counts of subgraphs
			if(staticCounts.find(key) == staticCounts.end())
				staticCounts[key] =weight;
			else
				staticCounts[key] = staticCounts[key] +weight;
			return;
		}
		void addSubgraph(std::string key, double w=1.0)
		{
			if(subgraphCount.find(key) == subgraphCount.end())
			{
				subgraphCount[key] = w;
			}
			else
				subgraphCount[key] = subgraphCount[key] + w;
			return;
		}
		const std::unordered_map<std::string, double>& getSubgraphCounts()
		{
			return std::cref(subgraphCount);
		}
		void addNodeToInstance(long long int inst_idx, long long int node)
		{
			instance_node_map[inst_idx].push_back(node);
			node_instance_map[node].insert(inst_idx);
			node_counts[node] += 1.0;
			return;
		}
		void updateNodeEstimate(long long int node, double value, int type=0)
		{
			if(type == 0)
				node_counts[node] += value;
			else if(type ==1)
			{
				node_counts[node] /= value;
			}
			else if(type ==2)
			{
				node_counts[node] = value;
			}
			return;
		}
		double getNodeEstimate(long long int node)
		{
			return node_counts[node];
		}
		void updateEstimate(double value)
		{
				solutionInsts += value;
				return;
		}
		void resetEstimate()
		{
			solutionInsts = 0;
			return;
		}
		int getAlives()
		{
			//return alive;
			return solution.size();
		}
		double getNodeCount(long long int node)
		{
			if(solution.find(node) != solution.end())
				return node_counts[node];
			else
				return -1;
		}
		long long int getMinNodeIdx() //TODO: Change by heap
		{
			return heapcounts.top().node;
		}
		void fixInstances()
		{
			long long int pos = 0;
			for(auto i : solution)
			{
				heap_data nodecount(i, node_counts[i]);
				boost::heap::fibonacci_heap<heap_data, boost::heap::compare<compare_node>>::handle_type handle = heapcounts.push(nodecount);
				node_handles_heap.push_back(handle); // store i'th node handle in i-th position 
				handleMap[i] = pos;
				pos++;
			}
			return;
		}
		std::vector<long long int> getNodesLowerThan(double value)
		{
			std::vector<long long int> answ;
			for(auto& node : solution)
			{
				if(node_counts[node] <= value)
					answ.push_back(node);
			}
			return answ;

		}
		double getCurrentSolutionValue() //TODO: Change by keeping sol
		{
			if (solution.size() == 0)
				return 0;
			return static_cast<double>(solutionInsts/solution.size());
		}
		const std::unordered_set<long long int>& getCurrentSolution() //TODO: Change by keeping sol
		{
			return std::cref(solution);
		}
		void peelNodeApprox(long long int id) 
		{
			solutionInsts -= node_counts[id];
			solution.erase(id);
			return;

		}
		void peelNodeNoTree(long long int id) //TODO: Change by keeping sol
		{
				if(node_counts[id] > 0)
				{
					for(auto& inst : node_instance_map[id]) // this instances should be removed, The node is peeled!
					{
						for(auto& node : instance_node_map[inst]) // This node should not reference anymore the instance
						{
							if(node != id)
							{
								node_instance_map[node].erase(inst);
								node_counts[node]--;
							}
						}

					}
				}
				node_counts[id] = std::numeric_limits<double>::max();
				solution.erase(id);
				return;
		}
		void peelNodeNoTreeEncode(long long int id) //TODO: Change by keeping sol
		{
				if(node_counts[id] > 0)
				{
					solutionInsts -= node_counts[id];
					for(auto& subgraphEncoding : nodeSubgraphs[id]) // this instances should be removed, The node is peeled!
					{
						std::istringstream iss(subgraphEncoding);
						std::vector<std::string> knodes(std::istream_iterator<std::string>{iss},std::istream_iterator<std::string>());
						for(auto& nodeSt : knodes) // This node should not reference anymore the instance
						{
							long long int node = stoi(nodeSt);
							if(node != id)
							{
								nodeSubgraphs[node].erase(subgraphEncoding);
								node_counts[node] -= staticCounts[subgraphEncoding];
								solutionInsts-= staticCounts[subgraphEncoding];
							}
						}

					}
				}
				node_counts[id] = 0;
				solution.erase(id);
				return;
		}
		void peelNode(long long int id) //TODO: Change by keeping sol
		{
				heapcounts.pop();
				solution.erase(id);
				if(node_counts[id] > 0)
				{
					for(auto& inst : node_instance_map[id]) // this instances should be removed, The node is peeled!
					{
						for(auto& node : instance_node_map[inst]) // This node should not reference anymore the instance
						{
							if(node != id) 
							{
								node_instance_map[node].erase(inst);
								auto nodehand = handleMap[node];
								double nodecount = (*(node_handles_heap[nodehand])).count;
								(*(node_handles_heap[nodehand])).count = nodecount -1;
								heapcounts.update(node_handles_heap[nodehand]);
								node_counts[node]--;
							}
						}
					}
				}
				return;
		}
		void peelNodeUseless(long long int id) //TODO: Change by keeping sol
		{
				if(node_counts[id] > 0 )
				{
					std::cout << "ERROR while processing: " << id << std::endl;
				}
				else
					solution.erase(id);
				return;
		}
		void peelNodeEncode(long long int id) //TODO: Change by keeping sol
		{
				solution.erase(id);
				heapcounts.pop();
				if(node_counts[id] > 0)
				{
					solutionInsts -= node_counts[id];
					for(auto& subgraphEncoding : nodeSubgraphs[id]) // this instances should be removed, The node is peeled!
					{
						std::istringstream iss(subgraphEncoding);
						std::vector<std::string> knodes(std::istream_iterator<std::string>{iss},std::istream_iterator<std::string>());
						for(auto& nodeSt : knodes) // This node should not reference anymore the instance
						{
							long long int node = stoi(nodeSt);
							if(node != id) 
							{
								nodeSubgraphs[node].erase(subgraphEncoding);
								auto nodehand = handleMap[node];
								double nodecount = (*(node_handles_heap[nodehand])).count;
								(*(node_handles_heap[nodehand])).count = nodecount -staticCounts[subgraphEncoding];
								heapcounts.update(node_handles_heap[nodehand]);
								node_counts[node]-= staticCounts[subgraphEncoding];
								solutionInsts-= staticCounts[subgraphEncoding];
							}
						}
					}
					nodeSubgraphs[id].clear();
				}
				node_counts[id] = 0;
				return;
		}
		~MyStructure() = default;
};

class MotifProperties{
	public:
		std::vector<int> inDegreeStatic; // i-th position will keep the in degree of the i-th node in the static projected and directed graph of M
		std::vector<int> outDegreeStatic; // i-th position will keep the out degree of the i-th node in the static projected and directed graph of M
		std::vector<int> degreeStatic;// i-th position will keep the degree of the i-th node in the static projected and undirected graph of M
		std::vector<int> inDegreeTemporal; // i-th position will keep the in degree of the i-th node in the multigraph corresponding to the motif
		std::vector<int> outDegreeTemporal; // i-th position will keep the out degree of the i-th node in the multigraph corresponding to the motif
		std::vector<int> matchingOrder; //i-th entry keeps the ordering in which the i-th temporal motif has to be matched
		std::vector<std::pair<int,int>> prev_next; // i-th entry keeps the indices of the edges already matched (only tight ones arew maintained)
		std::vector<int> revMatchingOrder; //i-th entry keeps the ordering in index in matchingOrder
		
		MotifProperties(int nodes, int ell) : inDegreeStatic(nodes, 0), outDegreeStatic(nodes, 0), degreeStatic(nodes, 0), 
										inDegreeTemporal(nodes, 0), outDegreeTemporal(nodes, 0), matchingOrder(ell, 0), prev_next(ell, {-1,-1}),
										revMatchingOrder(ell, -1) {}
		MotifProperties() {}
};

class RunningTimeStats{
	public:
		RunningTimeStats() : timeInitStructures{0.0}, timeRunBacktrackingAlg{0.0} {}
		const double getTimeInitStructures()
		{
			return timeInitStructures;
		}

		const double getTimeRunBacktrackingAlg()
		{
			return timeRunBacktrackingAlg;
		}

		void updateTimeInitStructures(double time)
		{
			timeInitStructures += time;
		}

		void updateTimeRunBacktrackingAlg(double time)
		{
			timeRunBacktrackingAlg += time;
		}
	private:
		double timeInitStructures;
		double timeRunBacktrackingAlg;
};


template<typename R>
  bool is_ready(std::future<R> const& f)
  { return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready; }


struct Partition
{
	unsigned int start;
	long tStart;
	unsigned int end;
	long tEnd;
	bool patch;

	const bool operator<(const Partition& p) const
	{
		if (tStart != p.tStart)
			return tStart < p.tStart;
		if ((tStart == p.tStart) && (tEnd != p.tEnd))
			return tEnd < p.tEnd;
		return false;
	}
};



class TempMotifDense {
public:
	static double get_wall_time();
	struct Datum {
		// Vertex based data structures
		std::vector<int> edgeCount;
    	std::vector<int> mapGM;
    	std::vector<int> mapMG;

		// For each node v_new_id keep v_old_id
		std::vector<int> originalIds;
    	// Edge index based adjacency list
		// Why do we need explicit representation of edges?
		std::vector< std::vector< TEdge > > adj_list;
		std::vector< std::vector< TEdge > > revadj_list;
		std::vector< std::unordered_map<int, std::vector<TEdge> > > adjMap;
		// for each node u keep the list of indexes of (u,v,t) or (v,u,t)
		std::vector< std::vector<int> > node_edge_index_;

		void InitializeStructures(std::vector< TEdge >& edgeList, int Vm) {
			std::unordered_map<int, int> remap;
			int id = 0;
			for (auto e : edgeList) {
				if (!remap.count(e.src)) {
					remap[e.src] = id++;
				}
				if (!remap.count(e.dst)) {
					remap[e.dst] = id++;
				}
			}
			//Keeping original IDs
			originalIds.resize(id, 0);
			for(auto& pair : remap)
			{
				originalIds[pair.second] = pair.first;
			}
			for (auto& e : edgeList) {
				e.src = remap[e.src];
				e.dst = remap[e.dst];
			}
			sort(edgeList.begin(), edgeList.end());

			// Now id is |V| in this set of edges
			for (long long int i = 0; i < (long long int) edgeList.size(); i++) {
				edgeList[i].id = i;
			}

			adj_list.resize(id);
			revadj_list.resize(id);
			adjMap.resize(id);
			for (const auto& edge : edgeList) {
				adj_list[edge.src].push_back(edge);
				revadj_list[edge.dst].push_back(edge);
				adjMap[edge.src][edge.dst].push_back(edge);
			}
			edgeCount.resize(id, 0);
			mapGM.resize(id, -1);
			mapMG.resize(Vm, -1);
			
		} 
	};

	struct DatumV2 {
		// Vertex based data structures
		std::vector<int> edgeCount;
    	std::vector<int> mapGM;
    	std::vector<int> mapMG;
    	// Edge index based adjacency list
		// Why do we need explicit representation of edges?
		std::vector< std::vector< TEdge > > adj_list;
		std::vector< std::vector< TEdge > > revadj_list;
		std::vector< std::unordered_map<int, std::vector<TEdge> > > adjMap;
		// for each node u keep the list of indexes of (u,v,t) or (v,u,t)
		//std::vector< std::vector<int> > node_edge_index_;
		std::vector<int> staticUndirectedDegrees; // i-th node degree undirected static
		std::vector<int> staticDirectedDegreesIn; // i-th node in degree directed static
		std::vector<int> staticDirectedDegreesOut; // i-th node out degree directed static
		std::vector<int> temporalDegreesIn; // i-th node in degree temporal
		std::vector<int> temporalDegreesOut; // i-th node out degree temporal
		// For each node v_new_id keep v_old_id
		std::vector<int> originalIds;

		void InitializeStructures(std::vector< TEdge >& edgeList, int Vm) {
			std::unordered_map<int, int> remap;
			int id = 0;
			for (auto e : edgeList) {
				if (!remap.count(e.src)) {
					remap[e.src] = id++;
				}
				if (!remap.count(e.dst)) {
					remap[e.dst] = id++;
				}
			}

			//Keeping original IDs
			originalIds.resize(id, 0);
			for(auto& pair : remap)
			{
				originalIds[pair.second] = pair.first;
				//std::cout << pair.second << " " << pair.first << '\n';
			}
			//std::cout << "---------------------------\n";

			for (auto& e : edgeList) {
				e.src = remap[e.src];
				e.dst = remap[e.dst];
			}
			sort(edgeList.begin(), edgeList.end());

			// Now id is |V| in this set of edges
			for (int i = 0; i < (int) edgeList.size(); i++) {
				edgeList[i].id = i;
			}

			TNGraph directedGraph = TNGraph(); 
			TUNGraph undirectedGraph = TUNGraph(); 
			// Adding the nodes to the static graphs
			/*
			for(int i{0}; i < id; i++)
			{
					directedGraph.AddNode(id);
					undirectedGraph.AddNode(id);
			}
			*/

			staticUndirectedDegrees.resize(id, 0);
			staticDirectedDegreesIn.resize(id, 0);
			staticDirectedDegreesOut.resize(id, 0);
			temporalDegreesIn.resize(id, 0); 
			temporalDegreesOut.resize(id, 0);

			adj_list.resize(id);
			revadj_list.resize(id);
			adjMap.resize(id);
			for (const auto& edge : edgeList) {
				adj_list[edge.src].push_back(edge);
				revadj_list[edge.dst].push_back(edge);
				adjMap[edge.src][edge.dst].push_back(edge);

				//Updating temporal degrees
				temporalDegreesOut[edge.src] += 1;
				temporalDegreesIn[edge.dst] += 1;

				//Creating graphs undirected and directed by specifying edges
				if(!directedGraph.IsNode(edge.src))
					directedGraph.AddNode(edge.src);
				if(!directedGraph.IsNode(edge.dst))
					directedGraph.AddNode(edge.dst);
				if(!undirectedGraph.IsNode(edge.src))
					undirectedGraph.AddNode(edge.src);
				if(!undirectedGraph.IsNode(edge.dst))
					undirectedGraph.AddNode(edge.dst);
				directedGraph.AddEdge(edge.src, edge.dst);
				undirectedGraph.AddEdge(edge.src, edge.dst);
			}
			// Updating static degrees
			for(int i{0}; i < id; i++)
			{
				staticUndirectedDegrees[i] = undirectedGraph.GetNI(i).GetDeg();
				staticDirectedDegreesIn[i] = directedGraph.GetNI(i).GetInDeg();
				staticDirectedDegreesOut[i] = directedGraph.GetNI(i).GetOutDeg();
			}
			directedGraph.Clr();
			undirectedGraph.Clr();
				
			edgeCount.resize(id, 0);
			mapGM.resize(id, -1);
			mapMG.resize(Vm, -1);
		} 
	};

	TempMotifDense(std::string filenameG, std::string filenameM, bool parallel=false, bool esbool=false);

	static double ExactEnumerateMotifsOrder(
    	const int delta, 
    	const double window, 
    	std::vector<TEdge>* edgesP,
    	const int Vm,
    	std::vector< std::pair<int, int> > edgesM,
	   	const int type,
	   	const long long deltaT, 
		const std::vector< TEdge >& edges_,
	   	const double probability,
	   	std::vector<int>& occurr, const std::vector<Partition>& parts, const std::tuple<const unsigned int,const  unsigned int,const  unsigned int> lims,
		const MotifProperties& motifprops, RunningTimeStats& stats, MyStructure& node_inst) {

    	std::vector<TEdge>& edges = *edgesP;
		auto t1 = get_wall_time();
    	DatumV2 data;
    	data.InitializeStructures(edges, Vm);

    	double res = 0;
    	std::vector<int> eStack;

		t1 = get_wall_time();

		int l = edgesM.size();
		int lastedge{(motifprops.revMatchingOrder)[l-1]};
    	int eG = 0, eM = (motifprops.matchingOrder)[0], uG = -1, vG = -1, uM = -1, vM = -1;
		int matchedEdges = 0, lastMatched = -1;
		double min_w_S{1.0};
		std::vector<long long int> matches(l, -1);// i-th position contains the index of the edges that maps on i-th edge of the motif in the ordering \sigma
    	const long long int INF = std::numeric_limits<long long>::max();
    	long long int _t = INF;
		long long int _t_lb = 0;
		long long int _t_ub = INF;
		unsigned int l3_i = std::get<0>(lims);; // lower index, mid index and upper index for case 3;
		unsigned int m3_i = std::get<1>(lims);; // lower index, mid index and upper index for case 3;
		unsigned int u3_i = std::get<2>(lims);; // lower index, mid index and upper index for case 3;
		int myindex = (int) deltaT; // Used only in case 6
		double decay = node_inst.getDecay();
		int weightingf =node_inst.getWeight();
	   	while (true) {
    		int last = eStack.empty() ? -1 : eStack.back();
    		eG = FindNextMatchOrder(eM, eG, data, _t, edges, last, edgesM, motifprops, _t_lb, _t_ub, matches);
    		TEdge edge = {0, 0, 0, INF};
			// A match is found
    		if (eG < (int) edges.size()) {
				if(matchedEdges == l-1) {
    				// Apply the weight function
    				if (type == 0)
					{
						if(weightingf==1)
							res += 1;
						else
						{
							double instw{0};
							for(int idd{0}; idd<= l-2; idd++)
							{
								int id1 = (l-1 == motifprops.revMatchingOrder[idd+1]) ? eG : eStack[motifprops.revMatchingOrder[idd+1]];
								int id2 = (l-1 == motifprops.revMatchingOrder[idd]) ? eG : eStack[motifprops.revMatchingOrder[idd]];
								instw += pow(std::exp(1.0), -decay*static_cast<double>(edges[id1].tim-edges[id2].tim));
							}
							instw/=(l-1);
							res+=instw;
						}
					}
					else if(type == 2)
					{
						int lastidx = (lastedge == l-1) ? eG : eStack[lastedge];
						double r_u = 1. * ((window*delta) - (edges[lastidx].tim - edges[eStack[0]].tim));
						std::unordered_set<long long int> currnodes;
						for(auto& edgeId : eStack)
						{
							currnodes.insert(data.originalIds[edges[edgeId].src]);
							currnodes.insert(data.originalIds[edges[edgeId].dst]);
						}
						currnodes.insert(data.originalIds[edges[eG].src]);
						currnodes.insert(data.originalIds[edges[eG].dst]);
						double weight{1};
						if(weightingf==2)
						{
							double instw{0};
							for(int idd{0}; idd<= l-2; idd++)
							{
								int id1 = (l-1 == motifprops.revMatchingOrder[idd+1]) ? eG : eStack[motifprops.revMatchingOrder[idd+1]];
								int id2 = (l-1 == motifprops.revMatchingOrder[idd]) ? eG : eStack[motifprops.revMatchingOrder[idd]];
								instw += pow(std::exp(1.0), -decay*static_cast<double>(edges[id1].tim-edges[id2].tim));
							}
							instw/=(l-1);
							weight=instw;
						}
						for(auto& node : currnodes)
						{
							node_inst.updateNodeEstimate(node, weight*deltaT/r_u);
						}
						res += (weight*deltaT/r_u);
					}
					else if(type == 6)
					{
						int lastidx = (lastedge == l-1) ? eG : eStack[lastedge];
						long long tf = edges[eStack[0]].tim;
						long long tl = edges[lastidx].tim;
						bool ispatch = parts[myindex].patch;
						if((myindex > 0) && (myindex < ((int) parts.size()-1) ))
						{
							if(ispatch)
							{
								if(tl <= parts[myindex-1].tEnd)
										;
								else if(tf >= parts[myindex+1].tStart)
										;
								else
										res += 1;
							}
							else
							{
								res +=1;
							}
						}
						else if((myindex == 0))
						{
								res +=1;
						}
						else if((myindex == ((int) parts.size() -1)) && ispatch)
						{
							//Is the last Patch!!
								if(tl <= parts[myindex-1].tEnd)
										;
								else
										res += 1;
						}
						else if((myindex == ((int) parts.size() -1)) && (!ispatch))
						{
								res += 1;	
						}
					}
					else if(type == 7)
					{
						int lastidx = (lastedge == l-1) ? eG : eStack[lastedge];
						long long tf = edges[eStack[0]].tim;
						long long tl = edges[lastidx].tim;
						bool ispatch = parts[myindex].patch;

						double instw{0};
						if(weightingf==1)
							instw=1.0;
						else if(weightingf==2)
						{
							for(int idd{0}; idd<= l-2; idd++)
							{
								int id1 = (l-1 == motifprops.revMatchingOrder[idd+1]) ? eG : eStack[motifprops.revMatchingOrder[idd+1]];
								int id2 = (l-1 == motifprops.revMatchingOrder[idd]) ? eG : eStack[motifprops.revMatchingOrder[idd]];
								instw += pow(std::exp(1.0), -decay*static_cast<double>(edges[id1].tim-edges[id2].tim));
							}
							instw/=(l-1);
							min_w_S = std::min(instw, min_w_S);
						}
						if((myindex > 0) && (myindex < ((int) parts.size()-1) ))
						{
							if(ispatch)
							{
								if(tl <= parts[myindex-1].tEnd)
										;
								else if(tf >= parts[myindex+1].tStart)
										;
								else
								{
									std::set<long long int> currnodes;
									for(auto& edgeId : eStack)
									{
										currnodes.insert(data.originalIds[edges[edgeId].src]);
										currnodes.insert(data.originalIds[edges[edgeId].dst]);
									}
									currnodes.insert(data.originalIds[edges[eG].src]);
									currnodes.insert(data.originalIds[edges[eG].dst]);
									node_inst.addInstance(currnodes,instw);
									res += instw;
								}
							}
							else
							{
								std::set<long long int> currnodes;
								for(auto& edgeId : eStack)
								{
									currnodes.insert(data.originalIds[edges[edgeId].src]);
									currnodes.insert(data.originalIds[edges[edgeId].dst]);
								}
								currnodes.insert(data.originalIds[edges[eG].src]);
								currnodes.insert(data.originalIds[edges[eG].dst]);
								node_inst.addInstance(currnodes,instw);
								res +=instw;
							}
						}
						else if((myindex == 0))
						{
								std::set<long long int> currnodes;
								for(auto& edgeId : eStack)
								{
									currnodes.insert(data.originalIds[edges[edgeId].src]);
									currnodes.insert(data.originalIds[edges[edgeId].dst]);
								}
								currnodes.insert(data.originalIds[edges[eG].src]);
								currnodes.insert(data.originalIds[edges[eG].dst]);
								node_inst.addInstance(currnodes,instw);
								res +=instw;
						}
						else if((myindex == ((int) parts.size() -1)) && ispatch)
						{
								if(tl <= parts[myindex-1].tEnd)
										;
								else
								{
									std::set<long long int> currnodes;
									for(auto& edgeId : eStack)
									{
										currnodes.insert(data.originalIds[edges[edgeId].src]);
										currnodes.insert(data.originalIds[edges[edgeId].dst]);
									}
									currnodes.insert(data.originalIds[edges[eG].src]);
									currnodes.insert(data.originalIds[edges[eG].dst]);
									node_inst.addInstance(currnodes,instw);
									res +=instw;
								}
						}
						else if((myindex == ((int) parts.size() -1)) && (!ispatch))
						{
							std::set<long long int> currnodes;
							for(auto& edgeId : eStack)
							{
								currnodes.insert(data.originalIds[edges[edgeId].src]);
								currnodes.insert(data.originalIds[edges[edgeId].dst]);
							}
							currnodes.insert(data.originalIds[edges[eG].src]);
							currnodes.insert(data.originalIds[edges[eG].dst]);
							node_inst.addInstance(currnodes,instw);
							res += instw;	
						}
					}
					else if(type == 8)
					{
						int lastidx = (lastedge == l-1) ? eG : eStack[lastedge];
						long long tf = edges[eStack[0]].tim;
						long long tl = edges[lastidx].tim;
						bool ispatch = parts[myindex].patch;


						double instw{0};
						if(weightingf==1)
							instw=1.0;
						else if(weightingf==2)
						{
							for(int idd{0}; idd<= l-2; idd++)
							{
								int id1 = (l-1 == motifprops.revMatchingOrder[idd+1]) ? eG : eStack[motifprops.revMatchingOrder[idd+1]];
								int id2 = (l-1 == motifprops.revMatchingOrder[idd]) ? eG : eStack[motifprops.revMatchingOrder[idd]];
								instw += pow(std::exp(1.0), -decay*static_cast<double>(edges[id1].tim-edges[id2].tim));
							}
							instw/=(l-1);
						}
						if((myindex > 0) && (myindex < ((int) parts.size()-1) ))
						{
							if(ispatch)
							{
								if(tl <= parts[myindex-1].tEnd)
										;
								else if(tf >= parts[myindex+1].tStart)
										;
								else
								{
									std::set<long long int> currnodes;
									for(auto& edgeId : eStack)
									{
										currnodes.insert(data.originalIds[edges[edgeId].src]);
										currnodes.insert(data.originalIds[edges[edgeId].dst]);
									}
									currnodes.insert(data.originalIds[edges[eG].src]);
									currnodes.insert(data.originalIds[edges[eG].dst]);
									std::string key = "";
									bool first = true;
									for(auto& node : currnodes)
									{
										if(first)
										{
											key = key + std::to_string(node);
											first = false;
										}
										else
											key = key + " " + std::to_string(node);
									}
									node_inst.addSubgraph(key, instw);
									res += instw;
								}
							}
							else
							{
								std::set<long long int> currnodes;
								for(auto& edgeId : eStack)
								{
									currnodes.insert(data.originalIds[edges[edgeId].src]);
									currnodes.insert(data.originalIds[edges[edgeId].dst]);
								}
								currnodes.insert(data.originalIds[edges[eG].src]);
								currnodes.insert(data.originalIds[edges[eG].dst]);
								std::string key = "";
								bool first = true;
								for(auto& node : currnodes)
								{
									if(first)
									{
										key = key + std::to_string(node);
										first = false;
									}
									else
										key = key + " " + std::to_string(node);
								}
								node_inst.addSubgraph(key, instw);
								res += instw;
							}
						}
						else if((myindex == 0))
						{
								std::set<long long int> currnodes;
								for(auto& edgeId : eStack)
								{
									currnodes.insert(data.originalIds[edges[edgeId].src]);
									currnodes.insert(data.originalIds[edges[edgeId].dst]);
								}
								currnodes.insert(data.originalIds[edges[eG].src]);
								currnodes.insert(data.originalIds[edges[eG].dst]);
								std::string key = "";
								bool first = true;
								for(auto& node : currnodes)
								{
									if(first)
									{
										key = key + std::to_string(node);
										first = false;
									}
									else
										key = key + " " + std::to_string(node);
								}
								node_inst.addSubgraph(key,instw);
								res += instw;
						}
						else if((myindex == ((int) parts.size() -1)) && ispatch)
						{
								if(tl <= parts[myindex-1].tEnd)
										;
								else
								{
									std::set<long long int> currnodes;
									for(auto& edgeId : eStack)
									{
										currnodes.insert(data.originalIds[edges[edgeId].src]);
										currnodes.insert(data.originalIds[edges[edgeId].dst]);
									}
									currnodes.insert(data.originalIds[edges[eG].src]);
									currnodes.insert(data.originalIds[edges[eG].dst]);
									std::string key = "";
									bool first = true;
									for(auto& node : currnodes)
									{
										if(first)
										{
											key = key + std::to_string(node);
											first = false;
										}
										else
											key = key + " " + std::to_string(node);
									}
									node_inst.addSubgraph(key,instw);
									res += instw;
								}
						}
						else if((myindex == ((int) parts.size() -1)) && (!ispatch))
						{
							std::set<long long int> currnodes;
							for(auto& edgeId : eStack)
							{
								currnodes.insert(data.originalIds[edges[edgeId].src]);
								currnodes.insert(data.originalIds[edges[edgeId].dst]);
							}
							currnodes.insert(data.originalIds[edges[eG].src]);
							currnodes.insert(data.originalIds[edges[eG].dst]);
							std::string key = "";
							bool first = true;
							for(auto& node : currnodes)
							{
								if(first)
								{
									key = key + std::to_string(node);
									first = false;
								}
								else
									key = key + " " + std::to_string(node);
							}
							node_inst.addSubgraph(key,instw);
							res += instw;
						}
					}
					if(eG < static_cast<int>(edges.size()))
						_t_lb = static_cast<long long int>(edges[static_cast<int>(eG+1)].tim);
					else
						eG = static_cast<int>(edges.size());
    			} else {
					// Register match
    				edge = edges[eG];
    				uG = edge.src, vG = edge.dst;
    				uM = edgesM[eM].first, vM = edgesM[eM].second;

    				data.mapGM[uG] = uM;
    				data.mapGM[vG] = vM;
    				data.mapMG[uM] = uG;
    				data.mapMG[vM] = vG;

    				data.edgeCount[uG] += 1;
    				data.edgeCount[vG] += 1;
    				eStack.push_back(eG);
					matches[eM] = eG;
					if (eM==0) {
						_t = edge.tim + delta; //Upper bound on the last edge's timestamp
						_t_ub = _t;
						_t_lb = edge.tim; // Lower bound on the time a timestamp can have
					}
					else
					{

						auto pos = (motifprops.prev_next)[(motifprops.matchingOrder)[matchedEdges+1]].second;
						if(pos != -1)
							_t_ub = edges[matches[pos]].tim;
						else
							_t_ub = _t;
						pos = (motifprops.prev_next)[(motifprops.matchingOrder)[matchedEdges+1]].first;
						// Always true since eM> 0 so in the worst case _t_lb = [edges[eStack[0]]].tim
						if(pos != -1)
							_t_lb = edges[matches[pos]].tim;
					}
					matchedEdges += 1;
					eM = (motifprops.matchingOrder)[matchedEdges];
					continue;
    			}
    		}
    		eG += 1;

    		while (eG >= (int) edges.size() || edge.tim > _t_ub) {
    			if (!eStack.empty()) {
					eG = matches[(motifprops.matchingOrder)[matchedEdges-1]] + 1;
    				eStack.pop_back();
				

    				edge = edges[eG - 1];
    				uG = edge.src, vG = edge.dst;
    				uM = edgesM[eM].first, vM = edgesM[eM].second;

    				if (eStack.empty()) { //aka eM = 0 after update
    					_t = INF;
						if(eG < static_cast<int>(edges.size()))
							_t_lb = edges[eG].tim;
						_t_ub = INF;
    				}
    				data.edgeCount[uG] -= 1;
    				data.edgeCount[vG] -= 1;

    				if (data.edgeCount[uG] == 0) {
    					uM = data.mapGM[uG];
    					data.mapMG[uM] = -1;
    					data.mapGM[uG] = -1;
    				}

    				if (data.edgeCount[vG] == 0) {
    					vM = data.mapGM[vG];
    					data.mapMG[vM] = -1;
    					data.mapGM[vG] = -1;
    				}

					
					matches[(motifprops.matchingOrder)[matchedEdges-1]] = -1;
					matchedEdges -= 1;
					eM = (motifprops.matchingOrder)[matchedEdges];
					if(eM != 0)
					{
						auto pos = (motifprops.prev_next)[(motifprops.matchingOrder)[matchedEdges]].second;
						if(pos != -1)
							_t_ub = edges[matches[pos]].tim;
						else
							_t_ub = _t;
						if(eG < static_cast<int>(edges.size()))
						{
							_t_lb = edges[eG].tim;
						}
					}
    			} else {
					node_inst.setMinWeight(min_w_S);
    				delete edgesP;
					return res;
    			}
    		}
    	}
    }

	static double ExactEnumerateMotifsOrderCall(
    	const int delta, 
    	const double window, 
    	std::vector<TEdge>* edgesP,
    	const int Vm,
    	std::vector< std::pair<int, int> > edgesM,
	   	const int type,
	   	const long long deltaT, 
		const std::vector< TEdge >& edges_,
	   	const double probability,
	   	std::vector<int>& occurr,
		const std::tuple<const unsigned int,const unsigned int,const unsigned int> lims,
		const MotifProperties& motifprops, RunningTimeStats& stats, MyStructure& node_inst) {
	std::vector<Partition> fill;
	return ExactEnumerateMotifsOrder(delta, window, edgesP, Vm, edgesM, type, deltaT, std::cref(edges_), probability, std::ref(occurr), std::cref(fill), lims, motifprops, std::ref(stats), node_inst);
	}

	static inline int FindNextMatchOrder(
    	int eM, 
    	int eG, 
    	DatumV2& data, 
    	long long int _t,
    	std::vector<TEdge>& edges,
    	int last,
    	std::vector< std::pair<int, int> >& edgesM,
		const MotifProperties& motifprops,
		long long int lb, long long int ub, std::vector<long long int>& matches) {

    	int uM, vM, uG, vG;
    	uM = edgesM[eM].first;
    	vM = edgesM[eM].second;
    	uG = data.mapMG[uM];
    	vG = data.mapMG[vM];

    	std::vector<TEdge>* S;
    	if (uG >= 0 && vG >= 0) {
    		S = &data.adjMap[uG][vG];
    	} else if (uG >= 0) {
    		S = &data.adj_list[uG];
    	} else if (vG >= 0) {
    		S = &data.revadj_list[vG];
    	} else {
    		S = &edges;
    	}
    	bool small = (S->size() < 16);
    	int head = 0;
		TEdge lb_edge = {0,0,lb,eG};
		
    	if (!small) {
	    	head = std::lower_bound(S->begin(), S->end(), lb_edge) - S->begin();
	    }

    	for (int i = head; i < static_cast<int>(S->size()); i++) {
    		auto edge = (*S)[i];
			// if looking for an edge greater that a one already matched
			// This can even happen in the presence of multiple time stamps with lower bound!
			// if(Ihave a prev && edge.idx < prev.id || edge.tim > ub)
			//std::cerr << edge << std::endl;
			if(edge.tim == lb_edge.tim && (edge.id < eG))
				continue;
			if(((motifprops.prev_next)[eM].first != -1) && ((edge.id < matches[(motifprops.prev_next)[eM].first]) || (edge.tim > ub) || (edge.tim < lb))){
				if (small) continue;
				else break;
			}
			if(eM == 0 && ((edge.tim > ub) || (edge.tim < lb))){
				if (small) continue;
				else break;
			}
    		if (((motifprops.prev_next)[eM].first != -1) && edge.tim <= edges[matches[(motifprops.prev_next)[eM].first]].tim) {
    			continue;
    		}
    		if (((motifprops.prev_next)[eM].second != -1) && edge.tim >= edges[matches[(motifprops.prev_next)[eM].second]].tim) {
    			continue;
    		}

    		int _eG = edge.id;
    		int _uG = edge.src, _vG = edge.dst;
			// Both assigned return
    		if ((uG == _uG) && (vG == _vG)) {
    				return _eG;
    			}
			// Both not assigned
    		else if ((uG < 0 && data.mapGM[_uG] < 0) && (vG < 0 && data.mapGM[_vG] < 0)) {
				if((data.staticUndirectedDegrees[_vG] >= motifprops.degreeStatic[vM]) &&
					(data.staticDirectedDegreesIn[_vG] >= motifprops.inDegreeStatic[vM]) &&
					(data.staticDirectedDegreesOut[_vG] >= motifprops.outDegreeStatic[vM]) &&
					(data.temporalDegreesIn[_vG] >= motifprops.inDegreeTemporal[vM]) &&
					(data.temporalDegreesOut[_vG] >= motifprops.outDegreeTemporal[vM]) &&
					(data.staticUndirectedDegrees[_uG] >= motifprops.degreeStatic[uM]) &&
					(data.staticDirectedDegreesIn[_uG] >= motifprops.inDegreeStatic[uM]) &&
					(data.staticDirectedDegreesOut[_uG] >= motifprops.outDegreeStatic[uM]) &&
					(data.temporalDegreesIn[_uG] >= motifprops.inDegreeTemporal[uM]) &&
					(data.temporalDegreesOut[_uG] >= motifprops.outDegreeTemporal[uM]))
				return _eG;
			}
			// Only one of the nodes already matched
    		else if ((uG < 0 && data.mapGM[_uG] < 0) || (vG < 0 && data.mapGM[_vG] < 0)) {
				std::pair<int, int> tomatch{-1, -1};
				if( (vG < 0 && data.mapGM[_vG] < 0) && (uG == _uG))
					tomatch = std::make_pair(_vG, vM);
				if( (uG < 0 && data.mapGM[_uG] < 0) && (vG == _vG))
					tomatch = std::make_pair(_uG, uM);
				if(tomatch.first==-1)
					continue;
				if((data.staticUndirectedDegrees[tomatch.first] >= motifprops.degreeStatic[tomatch.second]) &&
					(data.staticDirectedDegreesIn[tomatch.first] >= motifprops.inDegreeStatic[tomatch.second]) &&
					(data.staticDirectedDegreesOut[tomatch.first] >= motifprops.outDegreeStatic[tomatch.second]) &&
					(data.temporalDegreesIn[tomatch.first] >= motifprops.inDegreeTemporal[tomatch.second]) &&
					(data.temporalDegreesOut[tomatch.first] >= motifprops.outDegreeTemporal[tomatch.second]) )
					return _eG;
			}

    	}
    	return edges.size();
    }

	static double ProcessSampleOrder(int id, double delta, double c, int type, double deltaT,unsigned int first,unsigned int last, const std::tuple<const unsigned int,const  unsigned int,const  unsigned int>& lims ,
							int Vm_, std::vector<std::pair<int,int>> edgesM, const std::vector< TEdge >& myedges, const MotifProperties& motifprops, RunningTimeStats& stats, MyStructure& node_inst)
	{
        std::vector<TEdge>* tmp = new std::vector<TEdge>();
		tmp->insert(tmp->cbegin(), myedges.cbegin()+first, myedges.cbegin()+last); // Sample goes from [first,last)
		std::vector<int> a;	
		double result = ExactEnumerateMotifsOrderCall(delta, c, tmp, Vm_, edgesM, type, deltaT, std::cref(myedges), 0, std::ref(a), std::cref(lims), motifprops, std::ref(stats), node_inst);
		return result;
	}

	static double ProcessIntervalOrder(int id, int type, const double delta, const int myindex, const std::vector<Partition> parts, const int Vm_,
    	std::vector<std::pair<int,int>> edgesM, const std::vector< TEdge >& edges_, const MotifProperties& motifprops, RunningTimeStats& stats, MyStructure& node_inst)
	{
        std::vector<TEdge>* tmpv = new std::vector<TEdge>();
		tmpv->insert(tmpv->begin(), edges_.begin()+parts[myindex].start, edges_.begin()+parts[myindex].end);
		std::vector<int> a;	
		std::tuple<const unsigned int,const  unsigned int,const  unsigned int> lims;
		double result = ExactEnumerateMotifsOrder(delta, 0, tmpv, Vm_, edgesM, type, myindex, std::cref(edges_), 0, std::ref(a), std::cref(parts), std::cref(lims), motifprops, std::ref(stats), std::ref(node_inst));
		return result;
	}

	void ComputeDelta(long long* delta2, double c, int delta);
	
	double TwoDeltaPatch(int delta, double size, int threads);

	double TwoDeltaPatch(int delta, double size)
	{
		return TwoDeltaPatch(delta, size, THREADS);
	}

	std::vector<Partition> ComputeSplitsAndPatch(int delta, double size);

	int setMotif(std::string filenameM);

	double kApproxDenseMotif(int delta, double decay=2);

	double exactDense(int delta, double decay=2);

	double kApproxDenseMotifBatchPeel(int delta, double epsilon, double decay=2);

	double kApproxDenseMotifBatchPeelSampling(int delta, double epsilon, int samps, long long int seed=2324, double decay=2);

	void ComputeDelta(long long* delta2, double c, int delta, std::vector<TEdge>& tn);

	void fixMotifOrdering();

	std::vector<Partition> ComputeSplitsAndPatch(std::vector<TEdge>& myedges, int delta, double size);

	double kApproxDenseMotifBatchPeelSamplingHybrid(int delta, double epsilon, int samps, long long int seed, double decay, int iters=2);

	int getEll()
	{
		return l_;
	}
	const std::vector<long long int> getOptSol()
	{
		return std::cref(optimalSolution);
	}
	
private:
	static double get_wall_time_inner(){

		struct timeval time;
		if (gettimeofday(&time,NULL)){
			return 0;
		}
		return (double)time.tv_sec + (double)time.tv_usec * .000001;
	}
	std::vector< std::vector<TEdge> > all_components;
	std::vector< TEdge > edges_;
	std::vector<long long int> optimalSolution;
	std::vector< std::pair<int, int> > edgesM_;
    std::default_random_engine generatorLiu;
	MotifProperties* motifprop;
	MotifProperties motifprops;
	int Vm_;
	int l_ = 0;
	int N_ = 0;
	std::mutex mut;
};

#endif  // snap_temporalmotifdense_h
