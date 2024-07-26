#include "utilities.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <queue>

int Utils::LoadEdgeList(const std::string& filenameG, std::vector<TEdge>& edges_, std::vector<size_t>& temporal_degrees_){
	std::ifstream fileG(filenameG);
	if(fileG.is_open())
	{
		std::string line;
		long int ID = 0;
		int self_edges = 0;
		bool first = true;
		long long prev = 0;
		int max_node_id {0};
			while(getline(fileG, line))
			{
				std::istringstream iss(line);
				std::vector<std::string> results(std::istream_iterator<
							std::string>{iss},std::istream_iterator<std::string>());
				if(results.size() != 3)
				{
					std::cerr << "ERROR, dataset is not in format src dst timestamp" << std::endl;
					exit(1);
				}
				else
				{
					int src =stoi(results[0]);
					int dst =stoi(results[1]);
					max_node_id = std::max(std::max(max_node_id, src), dst);
					long long int timestamp = stoll(results[2]);
					if(src != dst)
					{
						TEdge edge = {.src=src, .dst=dst, .tim=timestamp, .id=ID++};
						edges_.push_back(edge);
					}
					else
						self_edges++;
				}
			}
			std::cout << "removed edges is: " << self_edges << std::endl;
			std::cout << "maximum id of a node: " << max_node_id << std::endl;
	}
	fileG.close();
	return 0;
}

int Utils::LoadTemporalMotif(const std::string& filenameM, std::vector<std::pair<int, int>>& edgesM_, int& Vm_, int& l_){

	std::ifstream fileM(filenameM);
	if(fileM.is_open())
	{
		std::string line;
		int max_nodes = 0;
		while(getline(fileM, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> results(std::istream_iterator<
						std::string>{iss},std::istream_iterator<std::string>());
			if(results.size() != 3)
			{
				std::cerr << "ERROR, motif is not in format src dst order" << std::endl;
				exit(1);
			}
			else
			{
				// Assume the edges are already ordered
				int src =stoi(results[0]);
				int dst =stoi(results[1]);
				if(src != dst)
				{
      				edgesM_.push_back(std::make_pair(src, dst));
	  				l_++; // counting the number of edges
					max_nodes = std::max(max_nodes, src + 1);
					max_nodes = std::max(max_nodes, dst + 1);
				}
			}
		}
		Vm_ = max_nodes;
	}

	fileM.close();

	return 0;
}

int Utils::LoadEdgeList(const std::string& filenameG, TVec< THash<TInt, TVec<TInt64> >> &temporal_data, std::vector<size_t>& temporal_degrees_, long long int& totedges, timestamp_t& timespan){
	std::ifstream fileG(filenameG);
	timestamp_t t1;
	timestamp_t tm;
	const char* DELIM = " "; // TODO: Possibly add custom delims
	char currentString[150]; // Bound on the number of characters on a line: 149
	if(fileG.is_open())
	{
		//std::string line;
		long int ID = 0;
		int self_edges = 0;
		bool first = true;
		long long prev = 0;
		int max_node_id {0};
			//while(getline(fileG, line))
			while (fileG.getline(currentString, 150, '\n'))
			{
				char* chunk = strtok(currentString, DELIM);
				int src = atoi(chunk);
				chunk = strtok(NULL, DELIM);
				int dst = atoi(chunk);
				max_node_id = std::max(std::max(max_node_id, src), dst);
				chunk = strtok(NULL, DELIM);
				long long int timestamp = atoll(chunk); 
				chunk = strtok(NULL, DELIM);
				if(src != dst)
				{
					temporal_data[src](dst).Add(timestamp);
					//Trial
					temporal_degrees_[src]++;
					temporal_degrees_[dst]++;
					totedges++;
					tm = timestamp;
					if(first)
					{
						first = false;
						t1 = timestamp;
					}
				}
				else
					self_edges++;
				delete chunk;
				/*
				std::istringstream iss(line);
				std::vector<std::string> results(std::istream_iterator<
							std::string>{iss},std::istream_iterator<std::string>());
				if(results.size() != 3)
				{
					std::cerr << "ERROR, dataset is not in format src dst timestamp" << std::endl;
					continue;
				}
				else
				{
					int src =stoi(results[0]);
					int dst =stoi(results[1]);
					max_node_id = std::max(std::max(max_node_id, src), dst);
					long long int timestamp = stoll(results[2]);
					if(src != dst)
					{
						temporal_data[src](dst).Add(timestamp);
						//Trial
						temporal_degrees_[src]++;
						temporal_degrees_[dst]++;
						totedges++;
						tm = timestamp;
						if(first)
						{
							first = false;
							t1 = timestamp;
						}
					}
					else
						self_edges++;
				}
				*/
			}
			std::cout << "removed edges is: " << self_edges << std::endl;
			std::cout << "maximum id of a node: " << max_node_id << std::endl;
	}
	//delete DELIM;
	fileG.close();
	timespan = tm -t1;
	return 0;
}

int Utils::RemapNodes(std::vector<TEdge>& edgeList){
	std::unordered_map<node_t, node_t> remap;
	node_t id = 0;
	for (auto e : edgeList) {
		if (!remap.count(e.src)) {
			remap[e.src] = id++;
		}
		if (!remap.count(e.dst)) {
			remap[e.dst] = id++;
		}
	}
	for (auto& e : edgeList) {
		e.src = remap[e.src];
		e.dst = remap[e.dst];
	}
	return 0;
}

int Utils::PrintGraph(std::vector<TEdge>& edgeList){
	for(auto& edge : edgeList)
		std::cout << edge.src << " " << edge.dst << " " << edge.tim << '\n';
	return 0;
}


int Utils::kDFSFromNode(PUNGraph graph, const int depthlimit, const int startingNode, std::unordered_set<staticEdge_t, DataStructures::hash_pair>& subgraph){
	subgraph.clear();	
	std::unordered_map<int,int> node_dist;
	std::stack<int> tovisit;
	// Insert first node
	tovisit.push(startingNode);
	node_dist[startingNode] = 0;
	while(!tovisit.empty())
	{
		int currNode = tovisit.top();
		tovisit.pop();
		int currentDist{ node_dist[currNode] };
		//std::cout << "visiting: " << currNode << " with distance "<< currentDist << '\n';
		if(currentDist == depthlimit) { continue; }
		TUNGraph::TNodeI nodeIcurrent = graph->GetNI(currNode);
		// Scan the neighborhood
		int nextDist{ currentDist+1 };
		for(int i=0; i < nodeIcurrent.GetDeg(); i++)
		{
			int neighbor = nodeIcurrent.GetNbrNId(i);
			// Update datastructures 
			staticEdge_t currEdge;
			if(currNode < neighbor)
				currEdge = std::make_pair<node_t, node_t>(currNode, neighbor);
			else
		   		currEdge = std::make_pair<node_t, node_t>(neighbor, currNode);
			//std::cout << "src: " << currEdge.first << " dst: " << currEdge.second << '\n';
			if ((subgraph.find(currEdge)) == subgraph.end()) // Not Found
			{
				auto neigmap = node_dist.find(neighbor);
				if( (neigmap != node_dist.end()) && (neigmap->second == depthlimit) )
						;
				else
					subgraph.insert(currEdge);
				/*
				else
				{
					if((neigmap->second != depthlimit) || )
						subgraph.insert(currEdge);
				}
				*/

			}
			if((node_dist.find(neighbor)) == node_dist.end()) // Not already visited
			{
				tovisit.push(neighbor);
				node_dist[neighbor] = nextDist;
			}
		}

	}

	return (static_cast<int>(node_dist.size()));
}

int Utils::kBFSFromNode(PUNGraph graph, const int depthlimit, const int startingNode, std::unordered_set<staticEdge_t, DataStructures::hash_pair>& subgraph){
	subgraph.clear();	
	std::unordered_map<int,int> node_dist;
	std::queue<int> tovisit;
	// Insert first node
	tovisit.push(startingNode);
	node_dist[startingNode] = 0;
	while(!tovisit.empty())
	{
		int currNode = tovisit.front();
		tovisit.pop();
		int currentDist{ node_dist[currNode] };
		//std::cout << "visiting: " << currNode << " with distance "<< currentDist << '\n';
		//if(currentDist == depthlimit) { continue; }
		TUNGraph::TNodeI nodeIcurrent = graph->GetNI(currNode);
		// Scan the neighborhood
		int nextDist{ currentDist+1 };
		for(int i=0; i < nodeIcurrent.GetDeg(); i++)
		{
			int neighbor = nodeIcurrent.GetNbrNId(i);
			// Update datastructures 
			staticEdge_t currEdge;
			if(currNode < neighbor)
				currEdge = std::make_pair<node_t, node_t>(currNode, neighbor);
			else
		   		currEdge = std::make_pair<node_t, node_t>(neighbor, currNode);
			//std::cout << "src: " << currEdge.first << " dst: " << currEdge.second << '\n';
			if ((subgraph.find(currEdge)) == subgraph.end()) // Edge has not been considered
			{
				auto neigmap = node_dist.find(neighbor); // is the neighbor already visited or in list?
				if( (neigmap != node_dist.end()) && (neigmap->second <= depthlimit) )
					subgraph.insert(currEdge);
			}
			if(currentDist != depthlimit)
			{
				if((node_dist.find(neighbor)) == node_dist.end()) // Not already visited
				{
					tovisit.push(neighbor);
					node_dist[neighbor] = nextDist;
				}
			}
		}

	}

	return (static_cast<int>(node_dist.size()));
}

int Utils::GetTrianglesFromNode(PUNGraph graph, const int startingNode, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles){
	TUNGraph::TNodeI sNode = graph->GetNI(startingNode);
	//int bag{0};
	for(int i=0; i < sNode.GetDeg(); i++)
	{
		int firstNeighbor{ sNode.GetNbrNId(i) };
		TUNGraph::TNodeI secondNode = graph->GetNI(firstNeighbor);

		//TODO: improve this by cycling over the min degree e.g. as it is now or the original node if it has a lower degree
		for(int j=0; j < secondNode.GetDeg(); j++)
		{
			int secondNeighbor{ secondNode.GetNbrNId(j) };
			if( (secondNeighbor != startingNode) && (firstNeighbor < secondNeighbor) && (graph->IsEdge(startingNode, secondNeighbor)) )
			{
				//if(bag < 1000)
				//{
						std::vector<std::pair<node_t, node_t>> triangle;
						triangle.push_back(std::make_pair<node_t, node_t>(startingNode, firstNeighbor) );
						triangle.push_back(std::make_pair<node_t, node_t>(firstNeighbor, secondNeighbor) );
						triangle.push_back(std::make_pair<node_t, node_t>(startingNode, secondNeighbor) );
						triangles.push_back(triangle);
						//bag++;
				//}
				//else
					//return 0;

				/*
				triangles.push_back(std::vector<std::pair<node_t, node_t>>{ std::make_pair<node_t, node_t>(startingNode, firstNeighbor) } );
				triangles.back().push_back(std::make_pair<node_t, node_t>(firstNeighbor, secondNeighbor));
				triangles.back().push_back(std::make_pair<node_t, node_t>(startingNode, secondNeighbor));
				*/
			}
		}
	}

	return 0;
}

int Utils::GetEdgesFromNode(PUNGraph graph, const int startingNode, std::vector<std::vector<std::pair<node_t, node_t> > >& edges){
	TUNGraph::TNodeI sNode = graph->GetNI(startingNode);
	for(int i=0; i < sNode.GetDeg(); i++)
	{
		int firstNeighbor{ sNode.GetNbrNId(i) };
		std::vector<std::pair<node_t, node_t>> edge;
		edge.push_back(std::make_pair(startingNode, firstNeighbor));
		edges.push_back(edge);
	}
}


int Utils::GetDirectedTrianglesFromNode(PNGraph graph, const int startingNode, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles, bool isCycle){
	TNGraph::TNodeI sNode = graph->GetNI(startingNode);
	if(isCycle)
	{
		;
	}
	else // not a cycle triangle e.g. u -> v, u -> w, w -> v
	{
		for(int i=0; i < sNode.GetOutDeg(); i++)
		{
			int firstNeighbor{ sNode.GetOutNId(i) };
			TNGraph::TNodeI secondNode = graph->GetNI(firstNeighbor);
	
			for(int j=0; j < secondNode.GetOutDeg(); j++)
			{
				int secondNeighbor{ secondNode.GetOutNId(j) };
				if( (secondNeighbor != startingNode) && (graph->IsEdge(startingNode, secondNeighbor)) )
				{
					std::vector<std::pair<node_t, node_t>> triangle;
					triangle.push_back(std::make_pair<node_t, node_t>(startingNode, firstNeighbor) );
					triangle.push_back(std::make_pair<node_t, node_t>(firstNeighbor, secondNeighbor) );
					triangle.push_back(std::make_pair<node_t, node_t>(startingNode, secondNeighbor) );
					triangles.push_back(triangle);
	
				}
			}
		}
	}
	
	

	return 0;
}



int Utils::GetSquaresFromNode(PUNGraph graph, const int startingNode, std::vector<std::vector<std::pair<node_t, node_t> > >& squares)
{
	TUNGraph::TNodeI sNode = graph->GetNI(startingNode);
	for(int i=0; i < sNode.GetDeg(); i++)
	{
		int firstNeighbor{ sNode.GetNbrNId(i) };
		TUNGraph::TNodeI secondNode = graph->GetNI(firstNeighbor);

		//TODO: improve this by cycling over the min degree e.g. as it is now or the original node if it has a lower degree
		//for(int j=0; j < secondNode.GetDeg(); j++)
		for(int j=0; j < sNode.GetDeg(); j++)
		{
			/*
			if( (secondNeighbor != startingNode) && (firstNeighbor < secondNeighbor) && (graph->IsEdge(startingNode, secondNeighbor)) )
			{
					std::vector<std::pair<node_t, node_t>> square;
					square.push_back(std::make_pair<node_t, node_t>(startingNode, firstNeighbor) );
					square.push_back(std::make_pair<node_t, node_t>(firstNeighbor, secondNeighbor) );
					square.push_back(std::make_pair<node_t, node_t>(startingNode, secondNeighbor) );
					squares.push_back(square);

			}
			*/
			int secondNeighbor{ sNode.GetNbrNId(j) };
			if(secondNeighbor >= firstNeighbor)
				continue;
			TUNGraph::TNodeI thirdNode = graph->GetNI(secondNeighbor);
			for(int k=0; k < thirdNode.GetDeg(); k++)
			{
				int thirdNeighbor{ thirdNode.GetNbrNId(k) };
				if((thirdNeighbor == startingNode) || (thirdNeighbor == firstNeighbor))
					continue;
				if(graph->IsEdge(firstNeighbor,  thirdNeighbor)) // 0 -> 1, 1-> 2, 2->3, 3->0 i.e. 4-path identified
				{
					std::vector<std::pair<node_t, node_t>> square;
					/*
					square.push_back(std::make_pair<node_t, node_t>(startingNode, firstNeighbor) );
					square.push_back(std::make_pair<node_t, node_t>(firstNeighbor, secondNeighbor) );
					square.push_back(std::make_pair<node_t, node_t>(secondNeighbor, thirdNeighbor) );
					square.push_back(std::make_pair<node_t, node_t>(thirdNeighbor, startingNode) );
					*/
					square.push_back(std::make_pair<node_t, node_t>(startingNode, firstNeighbor) );
					square.push_back(std::make_pair<node_t, node_t>(firstNeighbor, thirdNeighbor) );
					square.push_back(std::make_pair<node_t, node_t>(secondNeighbor, thirdNeighbor) );
					square.push_back(std::make_pair<node_t, node_t>(secondNeighbor, startingNode) );
					squares.push_back(square);
				}
			}
		}
	}

	return 0;

}



int Utils::GetTrianglesFromEdge(const PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles)
{
	TUNGraph::TNodeI firstNode = graph->GetNI(sampledEdge.first);
	TUNGraph::TNodeI secondNode = graph->GetNI(sampledEdge.second);
	//bool first = (firstNode.GetDeg() < secondNode.GetDeg()) ? true : false;
	TUNGraph::TNodeI cycleNode = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.first) : graph->GetNI(sampledEdge.second);
	TUNGraph::TNodeI toMatch = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.second) : graph->GetNI(sampledEdge.first);

	for(int i=0; i < cycleNode.GetDeg(); i++)
	{
		int firstNeighbor{ cycleNode.GetNbrNId(i) };

		if( (firstNeighbor != toMatch.GetId()) && (graph->IsEdge(firstNeighbor, toMatch.GetId())) )
		{
			std::vector<std::pair<node_t, node_t>> triangle;
			//triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), toMatch.GetId()) );
			triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), firstNeighbor) );
			triangle.push_back(std::make_pair<node_t, node_t>(firstNeighbor, toMatch.GetId()) );
			triangles.push_back(triangle);

		}
	}

	return 0;

}

int Utils::GetG6FromEdge(const PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles)
{
	TUNGraph::TNodeI firstNode = graph->GetNI(sampledEdge.first);
	TUNGraph::TNodeI secondNode = graph->GetNI(sampledEdge.second);
	//bool first = (firstNode.GetDeg() < secondNode.GetDeg()) ? true : false;
	TUNGraph::TNodeI cycleNode = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.first) : graph->GetNI(sampledEdge.second);
	TUNGraph::TNodeI toMatch = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.second) : graph->GetNI(sampledEdge.first);

	for(int i=0; i < cycleNode.GetDeg(); i++)
	{
		int firstNeighbor{ cycleNode.GetNbrNId(i) };

		if( (firstNeighbor != toMatch.GetId()) && (graph->IsEdge(firstNeighbor, toMatch.GetId())) )
		{
			// cycle first node to search for single edges
			for(int fn=0; fn < cycleNode.GetDeg();fn++)
			{
				int lastmatch{ cycleNode.GetNbrNId(fn) };
				if((lastmatch != firstNeighbor) && (lastmatch != toMatch.GetId()))
				{
					std::vector<std::pair<node_t, node_t>> g6;
					//triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), toMatch.GetId()) );
					g6.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), firstNeighbor) );
					g6.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), lastmatch) );
					g6.push_back(std::make_pair<node_t, node_t>(firstNeighbor, toMatch.GetId()) );
					triangles.push_back(g6);
				}
			}
			// cycle second node to search for single edges
			for(int fn=0; fn < toMatch.GetDeg();fn++)
			{
				int lastmatch{ toMatch.GetNbrNId(fn) };
				if((lastmatch != firstNeighbor) && (lastmatch != cycleNode.GetId()))
				{
					std::vector<std::pair<node_t, node_t>> g6;
					//triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), toMatch.GetId()) );
					g6.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), firstNeighbor) );
					g6.push_back(std::make_pair<node_t, node_t>(toMatch.GetId(), lastmatch) );
					g6.push_back(std::make_pair<node_t, node_t>(firstNeighbor, toMatch.GetId()) );
					triangles.push_back(g6);
				}
			}
			// cycle second node to search for single edges
			TUNGraph::TNodeI lastnode = graph->GetNI(firstNeighbor);
			for(int fn=0; fn < lastnode.GetDeg();fn++)
			{
				int lastmatch{ lastnode.GetNbrNId(fn) };
				if((lastmatch != cycleNode.GetId()) && (lastmatch != toMatch.GetId()))
				{
					std::vector<std::pair<node_t, node_t>> g6;
					//triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), toMatch.GetId()) );
					g6.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), firstNeighbor) );
					g6.push_back(std::make_pair<node_t, node_t>(firstNeighbor, lastmatch) );
					g6.push_back(std::make_pair<node_t, node_t>(firstNeighbor, toMatch.GetId()) );
					triangles.push_back(g6);
				}
			}

		}
	}

	return 0;

}

int Utils::GetTrianglesFromEdgeThreads(const PUNGraph& graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles)
{
	TUNGraph::TNodeI firstNode = graph->GetNI(sampledEdge.first);
	TUNGraph::TNodeI secondNode = graph->GetNI(sampledEdge.second);
	//bool first = (firstNode.GetDeg() < secondNode.GetDeg()) ? true : false;
	TUNGraph::TNodeI cycleNode = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.first) : graph->GetNI(sampledEdge.second);
	TUNGraph::TNodeI toMatch = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.second) : graph->GetNI(sampledEdge.first);

	for(int i=0; i < cycleNode.GetDeg(); i++)
	{
		int firstNeighbor{ cycleNode.GetNbrNId(i) };

		if( (firstNeighbor != toMatch.GetId()) && (graph->IsEdge(firstNeighbor, toMatch.GetId())) )
		{
			std::vector<std::pair<node_t, node_t>> triangle;
			//triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), toMatch.GetId()) );
			triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), firstNeighbor) );
			triangle.push_back(std::make_pair<node_t, node_t>(firstNeighbor, toMatch.GetId()) );
			triangles.push_back(triangle);

		}
	}

	return 0;

}


int Utils::GetWedgesFromEdge(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& wedges)
{
	TUNGraph::TNodeI firstNode = graph->GetNI(sampledEdge.first);
	TUNGraph::TNodeI secondNode = graph->GetNI(sampledEdge.second);

	for(int i=0; i < firstNode.GetDeg(); i++)
	{
		int firstNeighbor{ firstNode.GetNbrNId(i) };

		if( (firstNeighbor != sampledEdge.second))
		{
			std::vector<std::pair<node_t, node_t>> wedge;
			//wedge.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(sampledEdge.first), static_cast<node_t>(sampledEdge.second)) );
			wedge.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(sampledEdge.first), firstNeighbor) );
			wedges.push_back(wedge);

		}
	}

	for(int i=0; i < secondNode.GetDeg(); i++)
	{
		int firstNeighbor{ secondNode.GetNbrNId(i) };

		if( (firstNeighbor != sampledEdge.first))
		{
			std::vector<std::pair<node_t, node_t>> wedge;
			//wedge.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(sampledEdge.first), static_cast<node_t>(sampledEdge.second)) );
			wedge.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(sampledEdge.second), firstNeighbor) );
			wedges.push_back(wedge);

		}
	}
	return 0;

}


int Utils::GetSquaresFromEdge(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& squares)
{
	node_t src = static_cast<node_t>(sampledEdge.first);
	node_t dst = static_cast<node_t>(sampledEdge.second);
	TUNGraph::TNodeI firstNode = graph->GetNI(src);
	TUNGraph::TNodeI secondNode = graph->GetNI(dst);

	for(int i=0; i < firstNode.GetDeg(); i++)
	{
		int firstNeighbor{ firstNode.GetNbrNId(i) };
		TUNGraph::TNodeI firstNeighborNode = graph->GetNI(firstNeighbor);

		if(firstNeighbor != dst)
		{
			for(int j=0; j < firstNeighborNode.GetDeg(); j++)
			{
				int secondNeighbor{ firstNeighborNode.GetNbrNId(j) };

				if( (firstNeighbor < secondNeighbor) && (secondNeighbor != dst) && (secondNeighbor != src) && (graph->IsEdge(secondNeighbor, dst)) )
				{
					std::vector<std::pair<node_t, node_t>> square;
					//square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(firstNeighbor)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(secondNeighbor)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(secondNeighbor), static_cast<node_t>(dst)) );
					squares.push_back(square);
				}
			}
		}
	}

	for(int i=0; i < secondNode.GetDeg(); i++)
	{
		int firstNeighbor{ secondNode.GetNbrNId(i) };
		TUNGraph::TNodeI firstNeighborNode = graph->GetNI(firstNeighbor);

		if(firstNeighbor != src)
		{
			for(int j=0; j < firstNeighborNode.GetDeg(); j++)
			{
				int secondNeighbor{ firstNeighborNode.GetNbrNId(j) };

				if( (firstNeighbor < secondNeighbor) && (secondNeighbor != src) && (secondNeighbor != dst) && (graph->IsEdge(secondNeighbor, src)) )
				{
					std::vector<std::pair<node_t, node_t>> square;
					//square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(dst), static_cast<node_t>(firstNeighbor)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(secondNeighbor)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(secondNeighbor), static_cast<node_t>(src)) );
					squares.push_back(square);
				}
			}
		}
	}
	return 0;


}

int Utils::GetSquaresFromEdgeMin(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& squares)
{
	node_t src = static_cast<node_t>(sampledEdge.first);
	node_t dst = static_cast<node_t>(sampledEdge.second);
	TUNGraph::TNodeI firstNode = graph->GetNI(src);
	TUNGraph::TNodeI secondNode = graph->GetNI(dst);
	bool minFirst = (firstNode.GetDeg() < secondNode.GetDeg());


	if(minFirst)
	{
		for(int i=0; i < firstNode.GetDeg(); i++)
		{
			int firstNeighbor{ firstNode.GetNbrNId(i) };
			TUNGraph::TNodeI firstNeighborNode = graph->GetNI(firstNeighbor);

			if(firstNeighbor != dst)
			{
				for(int j=0; j < firstNeighborNode.GetDeg(); j++)
				{
					int secondNeighbor{ firstNeighborNode.GetNbrNId(j) };

					if( (secondNeighbor != dst) && (secondNeighbor != src) && (graph->IsEdge(secondNeighbor, dst)) )
					{
						std::vector<std::pair<node_t, node_t>> square;
						//square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(firstNeighbor)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(secondNeighbor)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(secondNeighbor), static_cast<node_t>(dst)) );
						squares.push_back(square);
					}
				}
			}
		}
	}
	else
	{
		for(int i=0; i < secondNode.GetDeg(); i++)
		{
			int firstNeighbor{ secondNode.GetNbrNId(i) };
			TUNGraph::TNodeI firstNeighborNode = graph->GetNI(firstNeighbor);

			if(firstNeighbor != src)
			{
				for(int j=0; j < firstNeighborNode.GetDeg(); j++)
				{
					int secondNeighbor{ firstNeighborNode.GetNbrNId(j) };

					if((secondNeighbor != src) && (secondNeighbor != dst) && (graph->IsEdge(secondNeighbor, src)) )
					{
						std::vector<std::pair<node_t, node_t>> square;
						//square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(dst), static_cast<node_t>(firstNeighbor)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(secondNeighbor)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(secondNeighbor), static_cast<node_t>(src)) );
						squares.push_back(square);
					}
				}
			}
		}
	}
	return 0;
}

double Utils::GetSampledWedge(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& squares, TVec< THash<TInt, TInt64 >>& mapEdgeIdx, const std::vector<size_t>& weights, std::mt19937& generator, int alpha)
{
	TUNGraph::TNodeI firstNode = graph->GetNI(sampledEdge.first);
	TUNGraph::TNodeI secondNode = graph->GetNI(sampledEdge.second);
	node_t src = static_cast<node_t>(sampledEdge.first);
	node_t dst = static_cast<node_t>(sampledEdge.second);

	//Uniform sampling among candidates!!
	long long degSrc = firstNode.GetDeg();
	long long degDst = secondNode.GetDeg(); 
	if( (degSrc + degDst -2) == 0)
		return -1;
	std::vector<long long> massEdges(degSrc+degDst);
	std::fill(massEdges.begin(), massEdges.end(), 1);
	std::discrete_distribution<long long int> sampledist(massEdges.begin(), massEdges.end());
	int samps = alpha;
	int currsamps = 0;
	while(currsamps < samps)
	{
		long long edgeIdx = sampledist(generator);

		if(edgeIdx < degSrc)
		{
			//Check first Neighborhood
			int firstNeighbor{ firstNode.GetNbrNId(edgeIdx) };
			if(firstNeighbor != dst)
			{
				std::vector<std::pair<node_t, node_t>> square;
				square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
				square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(firstNeighbor)) );
				squares.push_back(square);

				currsamps++;
			}
		}
		else
		{
			//Check dst neighborhood
			int firstNeighbor{ secondNode.GetNbrNId(edgeIdx-degSrc) };
			if(firstNeighbor != src)
			{
				std::vector<std::pair<node_t, node_t>> square;
				square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
				square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(dst), static_cast<node_t>(firstNeighbor)) );
				squares.push_back(square);

				currsamps++;
			}
		}
	}


	/*
	//Importance sampling among candidates!!
	std::vector<size_t> massEdges;
	std::vector<std::pair<node_t, node_t>> edgesCandidate;

	
	for(int i=0; i < firstNode.GetDeg(); i++)
	{
		int firstNeighbor{ firstNode.GetNbrNId(i) };

		if( (firstNeighbor != sampledEdge.second))
		{
			if(src < firstNeighbor)
			{
				long long idxEdge = mapEdgeIdx[src].GetDat(firstNeighbor);
				massEdges.push_back(weights[idxEdge]);
				edgesCandidate.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(firstNeighbor)) );
			}
			else if(src > firstNeighbor)
			{
				long long idxEdge = mapEdgeIdx[firstNeighbor].GetDat(src);
				massEdges.push_back(weights[idxEdge]);
				edgesCandidate.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(src)) );
			}

		}
	}

	for(int i=0; i < secondNode.GetDeg(); i++)
	{
		int firstNeighbor{ secondNode.GetNbrNId(i) };

		if( (firstNeighbor != sampledEdge.first))
		{
			if(dst < firstNeighbor)
			{
				long long idxEdge = mapEdgeIdx[dst].GetDat(firstNeighbor);
				massEdges.push_back(weights[idxEdge]);
				edgesCandidate.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(dst), static_cast<node_t>(firstNeighbor)) );
			}
			else if(dst > firstNeighbor)
			{
				long long idxEdge = mapEdgeIdx[firstNeighbor].GetDat(dst);
				massEdges.push_back(weights[idxEdge]);
				edgesCandidate.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(dst)) );
			}
		}
	}

	std::discrete_distribution<long long int> sampledist(massEdges.begin(), massEdges.end());
	if(massEdges.size() < 1)
		return -1;
	int samps { 10 };
	//int samps { 1 };
	for(int i = 0; i < samps; i++)
	{
		long long edgeIdx = sampledist(generator);
		std::pair<node_t, node_t> secondEdge = edgesCandidate[edgeIdx];

		std::vector<std::pair<node_t, node_t>> square;
		square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
		square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(secondEdge.first), static_cast<node_t>(secondEdge.second)) );
		squares.push_back(square);
	}
	*/

	//return sampledist.probabilities()[edgeIdx];
	return samps;
}
