#include "Snap.h"

#include "temporalmotifdense.h"
#include <random>
#include "ctpl.h"

#include <vector>
#include <algorithm>
#include <iomanip>
#include <queue>
#include "./graph.h"

#define WEIGHTING 1


TempMotifDense::TempMotifDense(std::string filenameG, std::string filenameM, bool parallel, bool esbool)
{
    int nodes = 0;
    std::ifstream fileG(filenameG);
    if(fileG.is_open())
    {
        std::string line;
        long int ID = 0;
        int self_edges = 0;
        bool first = true;
        long long prev = 0;
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
                long long int timestamp = stoll(results[2]);
                if(src > nodes)
                {
                    if(src > nodes) nodes = src;
                }
                if(dst > nodes)
                {
                    if(dst > nodes) nodes = dst;
                }
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
    }
    N_ = nodes +1;
    std::cout << "Nodes in graph: " << N_ << std::endl;
    fileG.close();
    // Enable if edges need to be sorted, we assume data is already preprocessed!
    //sort(edges_.begin(), edges_.end());

	// Loading the motif
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
	fixMotifOrdering();
}

void TempMotifDense::fixMotifOrdering()
{
	
	// Build motif graph
	auto start = get_wall_time();
	TNGraph directed_motif; 
	TUNGraph undirected_motif; 

	MotifProperties prop = MotifProperties(Vm_, l_);
	motifprop = &prop;


	// Keep track of the node degrees
	for(auto& edge : edgesM_)
	{
		if(!directed_motif.IsNode(edge.first))
			directed_motif.AddNode(edge.first);
		if(!directed_motif.IsNode(edge.second))
			directed_motif.AddNode(edge.second);
		if(!undirected_motif.IsNode(edge.first))
			undirected_motif.AddNode(edge.first);
		if(!undirected_motif.IsNode(edge.second))
			undirected_motif.AddNode(edge.second);
		directed_motif.AddEdge(edge.first, edge.second);
		undirected_motif.AddEdge(edge.first, edge.second);
		(motifprop->inDegreeTemporal)[edge.second] += 1;
		(motifprop->outDegreeTemporal)[edge.first] += 1;
	}

	// Assuming the nodes of the motif are labelled from 0 to Vm_-1
	// This may also be done with iterator on the nodes of the graphs
	for(int i{0}; i< Vm_; i++)
	{
		auto directedNode = directed_motif.GetNI(i);
		auto undirectedNode = undirected_motif.GetNI(i);
		(motifprop->inDegreeStatic)[i] = directedNode.GetInDeg();
		(motifprop->outDegreeStatic)[i] = directedNode.GetOutDeg();
		(motifprop->degreeStatic)[i] = undirectedNode.GetDeg();
	}

	// Matching the first edge always first, then expand by keeping the motif connected and by matching last first
	std::vector<bool> alreadyVisited(l_, false);
	// key: distance from the already explored component, std::vector<int> indices in edgesM_ of the edges with that distance
	//std::unordered_map<int, std::vector<int>> priorityEdges;
	// priority: -1 not to be considred, 0, mathcing an existing edge, 1 reverse edge from one already selected, 2 distance 1 from existing component
	std::vector<int> priority(l_, 3);

	std::queue<int> toExplore;
	// Mathching the first edge always first
	toExplore.push(0);

	int matched = 0;
	while(!toExplore.empty())
	{
		int indexMatching = toExplore.front();
		toExplore.pop();
		auto currentEdge = edgesM_[indexMatching];
		(motifprop->matchingOrder)[matched++] = indexMatching;
		alreadyVisited[indexMatching] = true;
		//std::cout << "current edge: " << currentEdge.first << " " << currentEdge.second << "\n";
		for(int j{0}; j< l_; j++)
		{
			if(alreadyVisited[j])
				continue;
			int distance = 3;
			if((edgesM_[j].first == currentEdge.first) && (edgesM_[j].second == currentEdge.second)) 
				distance = 0;
			else if((edgesM_[j].first == currentEdge.second) && (edgesM_[j].second == currentEdge.first))
				distance = 1;
			else if((edgesM_[j].first == currentEdge.first) || (edgesM_[j].second == currentEdge.second))
				distance = 2;
			else if((edgesM_[j].first == currentEdge.second) || (edgesM_[j].second == currentEdge.first))
				distance = 2;
			priority[j] = std::min(distance, priority[j]);
		}
		int indexmin = -1;
		int distmin = 4;
		//std::cout << "candidate edges and distance\n";
		for(int j{0}; j< l_; j++)
		{
			if(alreadyVisited[j])
				continue;
			if(distmin >= priority[j])
			{
				distmin = priority[j];
				indexmin=j;
			}
		}
		if(indexmin!=-1)
			toExplore.push(indexmin);
	}


    /*
	for(int i{0}; i< Vm_; i++)
	{
		std::cout << "Node: " << i << " of the motif, has the following degrees\n";
		std::cout << "Temporal: (in, out): " << (motifprop->inDegreeTemporal)[i] << " , " << (motifprop->outDegreeTemporal)[i] << '\n';
		std::cout << "Directed static: (in, out): " << (motifprop->inDegreeStatic)[i] << " , " << (motifprop->outDegreeStatic)[i] << '\n';
		std::cout << "Undirected statuc: (in,out): " << (motifprop->degreeStatic)[i] << '\n';
	}
    */
	for(int i{0}; i < l_; i++)
	{
		auto curredge{ (motifprop->matchingOrder)[i] };
		auto myprev{-1};
		auto mynext{-1};
		for(int j{0}; j < i; j++)
		{
			if(i != j)
			{
				if((motifprop->matchingOrder)[i] > (motifprop->matchingOrder)[j])
				{
					if(myprev ==-1)
						myprev = (motifprop->matchingOrder)[j];
					else
						myprev = std::max(myprev, (motifprop->matchingOrder)[j]);
				}
				else if((motifprop->matchingOrder)[i] < (motifprop->matchingOrder)[j])
				{
					if(mynext ==-1)
						mynext = (motifprop->matchingOrder)[j];
					else
						mynext = std::min(mynext, (motifprop->matchingOrder)[j]);
				}
			}
		}
		//std::cout << " curredge: " << curredge << " prev,next: " << myprev << " " << mynext << std::endl;
		(motifprop->prev_next)[(motifprop->matchingOrder)[i]] = std::make_pair(myprev, mynext);

		//std::cout << "i: " << (motifprop->matchingOrder)[i]+1 << " edge is matched as: " << i+1 << "\n";
		(motifprop->revMatchingOrder)[(motifprop->matchingOrder)[i]] = i;

	}

    /*
	for(int i{0}; i < l_; i++)
	{
		std::cout << "i: " << (motifprop->matchingOrder)[i] << " prev: " << (motifprop->prev_next)[i].first << " next: " << (motifprop->prev_next)[i].second<< '\n';
		std::cout << "i: revmatching " << (motifprop->revMatchingOrder)[i] << "\n";
	}
    */
	motifprops = prop;
	auto end = get_wall_time();
	std::cout << "Time to compute motif ordering: " << (end-start) << '\n';

	// Form ordering enforcing last first and connectedness

}

int TempMotifDense::setMotif(std::string filenameM)
{
	std::ifstream fileM(filenameM);
	edgesM_.clear();
	l_ = 0;
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
					//G.motif(src, dst);
	  				l_++; // counting the number of edges
					max_nodes = std::max(max_nodes, src + 1);
					max_nodes = std::max(max_nodes, dst + 1);
				}
			}
		}
		Vm_ = max_nodes;
	}

	fileM.close();

	fixMotifOrdering();

	return 0;
}


double TempMotifDense::get_wall_time()
{
	struct timeval time;
	if (gettimeofday(&time,NULL)){
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void TempMotifDense::ComputeDelta(long long* delta2, double c, int delta, std::vector<TEdge>& tn)
{
		if(tn.size() - l_ < 0)
		{
			std::cerr << "Graph is too small" << std::endl;
			*delta2= 0;
		}
		*delta2 = floor(tn[tn.size()-l_-1].tim - tn[l_-1].tim + (c * delta));
}

double TempMotifDense::kApproxDenseMotifBatchPeelSampling(int delta, double epsilon, int samps, long long int seed, double decay)
{

		double c{1.25};
		//std::cout << " c is: " << c << ", delta is: " << delta << '\n';
		std::random_device rd;
		std::mt19937 generator(seed);
		long long delta2 = 0;
		auto s = samps; 

		std::vector<TEdge> currtn = edges_;
		MyStructure node_insts(N_);
		if(decay <1)
			node_insts.setDecay(decay);
		double bestSolutionValue = 0;
		std::unordered_set<long long int> sol;
		bool firstIteration = {true};
		bool F{false};
		auto startt = get_wall_time();
		double cM=0;
		int its{0};
		double epsApprox={0.95};
		double eta={0.1};
		while(node_insts.getAlives() >= Vm_)
		{

			//std::cout << "Nodes alive: " << node_insts.getAlives() << '\n';
			delta2 =0;
			ComputeDelta(&delta2, c, delta, std::ref(currtn));
			//int boundsamps = (delta2/((c-1.0)*delta)-1)*1/((1+epsApprox)*log(1+epsApprox) - epsApprox)*log(2/(eta/pow(2, its+2)));
			//std::cout << "Bound Samples at iteration i: " << its+1 << " is: " << boundsamps << '\n';
			if(delta2 == 0) {return -1.0;}
			
			long long min_t =  floor(currtn[l_-1].tim - (c*delta));
			long long max_t = currtn[currtn.size()-1-l_].tim; 

			std::uniform_int_distribution<long long> randomTime(min_t, max_t);
			

			double approx = 0;
			std::vector<double> probs;

			s = samps;
			//std::cout << "samps used is: " << s << '\n';
			long long new_samps = 0;
			double st{ get_wall_time() };
			
			int id=0;
			double timeInitStructs{0};
			RunningTimeStats stats = RunningTimeStats();
			while(new_samps < s)
			{	
				long long t_lower = randomTime(generator);
				long long t_upper = ceil(t_lower + (c * delta));
				long long bound = ceil(t_lower);
				TEdge target = {0,0, bound, 0};
				long long b1 = ceil(t_upper);
				TEdge targ = {0, 0, b1, 0};
				unsigned int first = lower_bound(currtn.begin(), currtn.end(), target) - currtn.begin();
				unsigned int last = upper_bound(currtn.begin()+first, currtn.end(), targ) - currtn.begin();
				new_samps++;
				if((last-first) < l_)
					continue;


				std::tuple<const unsigned int,const  unsigned int,const  unsigned int> lims{0, 0, 0};
				double est = ProcessSampleOrder(id, delta, c, 2, delta2, first, last,std::cref(lims),
								Vm_, edgesM_, std::cref(currtn), std::cref(motifprops), std::ref(stats), std::ref(node_insts));
				approx+= est;
				node_insts.increaseIteration();
			}
			auto solution = node_insts.getCurrentSolution();
			double sup=0.0;
			for(auto& node : solution)
			{
				node_insts.updateNodeEstimate(node, s, 1);
				double up = node_insts.getNodeEstimate(node);
				node_insts.updateEstimate(up);
			}

			auto curr_sol = node_insts.getCurrentSolutionValue()/Vm_;
			if(curr_sol > bestSolutionValue)
			{
				bestSolutionValue = curr_sol;
				sol= node_insts.getCurrentSolution();
				if(firstIteration)
				{
					cM = node_insts.getCurrentSolutionValue() * solution.size()/Vm_;
					firstIteration = false;
					if(abs(cM-0.00001) < 0)
						return 0;
				}
			}
			
			auto curr_sol_tmp = node_insts.getCurrentSolutionValue()/Vm_;
			auto small_nodes= node_insts.getNodesLowerThan((Vm_ + Vm_*epsilon) * curr_sol_tmp);
			std::unordered_set<long long int> toRemove;
			for(auto& node : small_nodes)
			{
				node_insts.peelNodeApprox(node);
				toRemove.insert(node);
			}
			std::vector<TEdge> newgraph;
			for(auto edge : currtn)
			{
				if((toRemove.find(edge.src) == toRemove.end()) && (toRemove.find(edge.dst) == toRemove.end()))
					newgraph.push_back(edge);
			}
			for(auto& node : node_insts.getCurrentSolution())
			{
				node_insts.updateNodeEstimate(node, 0, 2);
				node_insts.resetEstimate();
			}
			currtn = newgraph;
			if(F)
			{
				F=false;
				edges_.clear();
			}
			
			its++;
		}

		std::cout << "Inner Time enumerate: " << get_wall_time() - startt << " s\n";
		//Evaluate solution
		std::unordered_set<long long int> solutionSet;
		for(auto& node: sol)
		{
			solutionSet.insert(node);
		}
		std::vector<TEdge> newgraph;
		for(auto edge : edges_)
		{
			if((solutionSet.find(edge.src) != solutionSet.end()) && (solutionSet.find(edge.dst) != solutionSet.end()))
				newgraph.push_back(edge);
		}
		std::tuple<const unsigned int,const  unsigned int,const  unsigned int> lims{0, 0, 0};
		RunningTimeStats stats = RunningTimeStats();
		double insts = ProcessSampleOrder(0, delta, c, 0, delta2, 0, newgraph.size(),std::cref(lims),
						Vm_, edgesM_, std::cref(newgraph), std::cref(motifprops), std::ref(stats), std::ref(node_insts));

	std::cout << "SOLUTION: Current solution value is: " << bestSolutionValue << " and has: "<< sol.size() << " nodes accounting for: "<< bestSolutionValue*sol.size() << " instances, fraction: " << (bestSolutionValue*sol.size())/cM << '\n'; 
	std::cout << "SOLUTIONEX: Current solution value is: " << insts/sol.size() << " and has: "<< sol.size() << " nodes accounting for: "<< insts << " instances, fraction: " << (insts)/cM << '\n'; 
	std::cout << "SSSIZE: "<< samps << " SOLUTIONEX: Current solution value is: " << insts/sol.size() << " and has: "<< sol.size() << " nodes accounting for: "<< insts << " instances, fraction: " << (insts)/cM << '\n'; 
	for(auto& node : sol)
	{
		std::cout << "Node: " << node << '\n';
	}

	return cM;
}


double TempMotifDense::kApproxDenseMotifBatchPeelSamplingHybrid(int delta, double epsilon, int samps, long long int seed, double decay, int iters)
{

		double c{1.25};
		std::cout << " c is: " << c << ", delta is: " << delta << '\n';
		std::random_device rd;
		std::mt19937 generator(seed);
		long long delta2 = 0;
		auto s = samps; 

		std::vector<TEdge> currtn = edges_;
		MyStructure node_insts(N_);
		if(decay <1)
			node_insts.setDecay(decay);
		double bestSolutionValue = 0;
		std::unordered_set<long long int> sol;
		bool firstIteration = {true};
		bool F{false};
		auto startt = get_wall_time();
		double cM=0;
		int its{0};
		double epsApprox={0.95};
		double eta={0.1};
		
		std::vector<bool> removed(N_, false);
		std::vector<bool> best_removed_sol{removed};
		long long int cnt_useless = 0;
		bool isBatchSol{true};
		while(its < iters)
		{

			//std::cout << "Nodes alive: " << node_insts.getAlives() << " that induce edges: " << currtn.size() << '\n';
			delta2 =0;
			ComputeDelta(&delta2, c, delta, std::ref(currtn));
			//int boundsamps = (delta2/((c-1.0)*delta)-1)*1/((1+epsApprox)*log(1+epsApprox) - epsApprox)*log(2/(eta/pow(2, its+2)));
			//std::cout << "Bound Samples at iteration i: " << its+1 << " is: " << boundsamps << '\n';
			if(delta2 == 0) {return -1.0;}
			
			long long min_t =  floor(currtn[l_-1].tim - (c*delta));
			long long max_t = currtn[currtn.size()-1-l_].tim; 

			std::uniform_int_distribution<long long> randomTime(min_t, max_t);
			

			double approx = 0;
			std::vector<double> probs;

			s = samps;
			std::cout << "samps used is: " << s << '\n';
			long long new_samps = 0;
			double st{ get_wall_time() };
			
			int id=0;
			double timeInitStructs{0};
			RunningTimeStats stats = RunningTimeStats();
			while(new_samps < s)
			{	
				long long t_lower = randomTime(generator);
				long long t_upper = ceil(t_lower + (c * delta));
				long long bound = ceil(t_lower);
				TEdge target = {0,0, bound, 0};
				long long b1 = ceil(t_upper);
				TEdge targ = {0, 0, b1, 0};
				unsigned int first = lower_bound(currtn.begin(), currtn.end(), target) - currtn.begin();
				unsigned int last = upper_bound(currtn.begin()+first, currtn.end(), targ) - currtn.begin();
				new_samps++;
				if((last-first) < l_)
					continue;


				std::tuple<const unsigned int,const  unsigned int,const  unsigned int> lims{0, 0, 0};
				double est = ProcessSampleOrder(id, delta, c, 2, delta2, first, last,std::cref(lims),
								Vm_, edgesM_, std::cref(currtn), std::cref(motifprops), std::ref(stats), std::ref(node_insts));
				approx+= est;
				node_insts.increaseIteration();
			}
			auto solution = node_insts.getCurrentSolution();
			double sup=0.0;
			for(auto& node : solution)
			{
				node_insts.updateNodeEstimate(node, s, 1);
				double up = node_insts.getNodeEstimate(node);
				node_insts.updateEstimate(up);
			}

			auto curr_sol = node_insts.getCurrentSolutionValue()/Vm_;
			if(curr_sol > bestSolutionValue)
			{
				bestSolutionValue = curr_sol;
				sol= node_insts.getCurrentSolution();
				best_removed_sol = removed;
				if(firstIteration)
				{
					cM = node_insts.getCurrentSolutionValue() * solution.size()/Vm_;
					firstIteration = false;
					if(abs(cM-0.00001) < 0)
						return 0;
				}
			}
			
			auto curr_sol_tmp = node_insts.getCurrentSolutionValue()/Vm_;
			auto small_nodes= node_insts.getNodesLowerThan((Vm_ + Vm_*epsilon) * curr_sol_tmp);
			std::unordered_set<long long int> toRemove;
			for(auto& node : small_nodes)
			{
				node_insts.peelNodeApprox(node);
				toRemove.insert(node);
				removed[node] = true;
				cnt_useless++;
			}
			std::vector<TEdge> newgraph;
			for(auto edge : currtn)
			{
				if((toRemove.find(edge.src) == toRemove.end()) && (toRemove.find(edge.dst) == toRemove.end()))
					newgraph.push_back(edge);
			}
			for(auto& node : node_insts.getCurrentSolution())
			{
				node_insts.updateNodeEstimate(node, 0, 2);
				node_insts.resetEstimate();
			}
			currtn = newgraph;
			if(F)
			{
				F=false;
				edges_.clear();
			}
			its++;
		}

		//std::cout << "Time randomized-batch: " << get_wall_time() - startt << " s\n";
		//std::cout << "Nodes alive after random peeling: " << node_insts.getAlives() << " with edges: " << currtn.size() << '\n';
		//std::cout << "Estimated instances on the whole graph cM: " << cM << " \n";
		std::vector<Partition> partitions = ComputeSplitsAndPatch(std::ref(currtn), delta, 2);

		// Starting Greedy peeling on the remaining nodes!

	auto init = get_wall_time();
	int iterations = (int) partitions.size();
	auto nodeAlives = node_insts.getCurrentSolution();
	RunningTimeStats stats = RunningTimeStats();
	cM=0;
	for(int i=0; i < iterations; i++)
	{
      	cM += ProcessIntervalOrder(i, 7, delta, i, std::cref(partitions),
							Vm_, edgesM_, std::cref(currtn), std::cref(motifprops), std::ref(stats), std::ref(node_insts));
	}
	//std::cout << "Actual unstances in the remaining network after random peeling: " << cM << " \n";
	auto timeEnum = get_wall_time() - init;
	init = get_wall_time();
	double minW =node_insts.getMinWeight();
	//std::cout << "Min weight ins: " << minW << '\n';
	auto useless_nodes= node_insts.getNodesLowerThan(minW/3);
	for(auto& node : useless_nodes)
	{
		node_insts.peelNodeUseless(node);
		removed[node] = true;
		cnt_useless++;
	}

	node_insts.fixInstances();
	double solStart = node_insts.getCurrentSolutionValue()/Vm_;
	if(solStart > bestSolutionValue)
	{
		bestSolutionValue = solStart;
		best_removed_sol = removed;
		isBatchSol = false;
	}
	std::vector<bool> removed_sol{removed};
	auto timeBuild = get_wall_time() - init;
	//std::cout << "Current solution value is: " << solStart << '\n'; 
	init = get_wall_time();
		
	for(int i=(N_-cnt_useless); i>Vm_; i--)
	{
		auto min_idx= node_insts.getMinNodeIdx();
		auto curr_sol_tmp = node_insts.getCurrentSolutionValue()/Vm_;
		node_insts.peelNodeEncode(min_idx);
		removed[min_idx] = true;
		auto curr_sol = node_insts.getCurrentSolutionValue()/Vm_;
		if(curr_sol > bestSolutionValue)
		{
			bestSolutionValue= curr_sol;
			best_removed_sol = removed;
			isBatchSol = false;
		}
	}
	std::unordered_set<long long int> solution;
	for(long long int i=0; i< best_removed_sol.size(); i++)
	{
		if(!best_removed_sol[i])
		{
			solution.insert(i);
		}
	}
	auto timePeel = get_wall_time() - init;
	if(isBatchSol)
		std::cout << "Solution is from randomized batch\n";
	else
		std::cout << "Solution is from greedy peeling\n";
	std::cout << "Time enumerate: " << timeEnum << " Time build: " << timeBuild << " Time peel: " << timePeel << '\n';
	std::cout << "SOLUTION: Current solution value is: " << bestSolutionValue << " and has: "<< solution.size() << " nodes accounting for: "<< bestSolutionValue*solution.size() << " instances, fraction: " << (bestSolutionValue*solution.size())/cM << '\n'; 
	std::unordered_set<long long int> solutionSet;
	for(auto& node : solution)
	{
		std::cout << "Node: " << node << '\n';
		solutionSet.insert(node);
	}
	std::vector<TEdge> solgraph;
	for(auto edge : edges_)
	{
		if((solutionSet.find(edge.src) != solutionSet.end()) && (solutionSet.find(edge.dst) != solutionSet.end()))
		{
			solgraph.push_back(edge);
		}
	}
	sort(solgraph.begin(), solgraph.end());
	std::tuple<const unsigned int,const  unsigned int,const  unsigned int> lims{0, 0, 0};
	double est = ProcessSampleOrder(0, delta, 2, 0, 0, 0, solgraph.size(),std::cref(lims),
								Vm_, edgesM_, std::cref(solgraph), std::cref(motifprops), std::ref(stats), std::ref(node_insts));
	double top = est/(1.0*solutionSet.size());
	int precision = std::numeric_limits<double>::digits10;
	std::cout << "Solution value: " << std::setprecision(precision) << top << " insts: " << est << '\n';

	return cM;
}

double TempMotifDense::kApproxDenseMotifBatchPeel(int delta, double epsilon, double decay)
{
	std::vector<Partition> partitions = ComputeSplitsAndPatch(delta, 2);
	
	MyStructure node_insts(N_);
	if(decay <1)
		node_insts.setDecay(decay);
	RunningTimeStats stats = RunningTimeStats();
	auto init = get_wall_time();
	int iterations = (int) partitions.size();
	double cM = 0;
	for(int i=0; i < iterations; i++)
	{
      	cM += ProcessIntervalOrder(i, 7, delta, i, std::cref(partitions),
							Vm_, edgesM_, std::cref(edges_), std::cref(motifprops), std::ref(stats), std::ref(node_insts));
	}
	auto timeEnum = get_wall_time() - init;
	double sol = node_insts.getCurrentSolutionValue()/Vm_;
	if(abs(sol-0.000001) < 0)
		return 0;
	std::vector<bool> removed(N_, false);
	std::vector<bool> removed_sol;
	std::cout << "Current solution value is: " << sol << '\n'; 
	init = get_wall_time();
	int its{0};
	while(node_insts.getAlives() > Vm_)
	{
		auto curr_sol_tmp = node_insts.getCurrentSolutionValue()/Vm_;
		auto small_nodes= node_insts.getNodesLowerThan((Vm_ + Vm_*epsilon) * curr_sol_tmp);
		for(auto& node : small_nodes)
		{
			node_insts.peelNodeNoTreeEncode(node);
			removed[node] = true;
		}
		auto curr_sol = node_insts.getCurrentSolutionValue()/Vm_;
		if(curr_sol > sol)
		{
			sol = curr_sol;
			removed_sol = removed;
		}
		its++;
	}
	std::unordered_set<long long int> solution;
	for(long long int i=0; i< removed_sol.size(); i++)
	{
		if(!removed_sol[i])
			solution.insert(i);
	}
	auto timePeel = get_wall_time() - init;
	std::cout << "Time enumerate: " << timeEnum << " Time peel: " << timePeel << '\n';
	std::cout << "SOLUTION: Current solution value is: " << sol << " and has: "<< solution.size() << " nodes accounting for: "<< sol*solution.size() << " instances, fraction: " << (sol*solution.size())/cM << '\n'; 
	std::cout << "Iterations: " << its << '\n';
	for(auto& node : solution)
	{
		std::cout << "Node: " << node << '\n';
	}

	return cM;
}

double TempMotifDense::kApproxDenseMotif(int delta, double decay)
{
	std::vector<Partition> partitions = ComputeSplitsAndPatch(delta, 2);
	
	MyStructure node_insts(N_);
	if(decay <1)
		node_insts.setDecay(decay);
	RunningTimeStats stats = RunningTimeStats();
	auto init = get_wall_time();
	int iterations = (int) partitions.size();
	double cM = 0;
	for(int i=0; i < iterations; i++)
	{
      	cM += ProcessIntervalOrder(i, 7, delta, i, std::cref(partitions),
							Vm_, edgesM_, std::cref(edges_), std::cref(motifprops), std::ref(stats), std::ref(node_insts));
	}
	auto timeEnum = get_wall_time() - init;
	init = get_wall_time();
	double minW =node_insts.getMinWeight();
	std::cout << "Min weight ins: " << minW << '\n';
	auto useless_nodes= node_insts.getNodesLowerThan(minW/3);
	std::vector<bool> removed(N_, false);
	long long int cnt_useless = 0;
	for(auto& node : useless_nodes)
	{
		node_insts.peelNodeUseless(node);
		removed[node] = true;
		cnt_useless++;
	}
	node_insts.fixInstances();
	double sol = node_insts.getCurrentSolutionValue()/Vm_;
	auto timeBuild = get_wall_time() - init;
	std::cout << "Current solution value is: " << sol << '\n'; 
	init = get_wall_time();
	std::vector<bool> removed_sol{removed};
		
	for(int i=(N_-cnt_useless); i>Vm_; i--)
	{
		auto min_idx= node_insts.getMinNodeIdx();
		auto curr_sol_tmp = node_insts.getCurrentSolutionValue()/Vm_;
		node_insts.peelNodeEncode(min_idx);
		removed[min_idx] = true;
		auto curr_sol = node_insts.getCurrentSolutionValue()/Vm_;
		if(curr_sol > sol)
		{
			sol = curr_sol;
			removed_sol = removed;
		}
	}
	std::unordered_set<long long int> solution;
	for(long long int i=0; i< removed_sol.size(); i++)
	{
		if(!removed_sol[i])
		{
			solution.insert(i);
		}
	}
	auto timePeel = get_wall_time() - init;
	std::cout << "Time enumerate: " << timeEnum << " Time build: " << timeBuild << " Time peel: " << timePeel << '\n';
	std::cout << "SOLUTION: Current solution value is: " << sol << " and has: "<< solution.size() << " nodes accounting for: "<< sol*solution.size() << " instances, fraction: " << (sol*solution.size())/cM << '\n'; 
	std::unordered_set<long long int> solutionSet;
	for(auto& node : solution)
	{
		std::cout << "Node: " << node << '\n';
		solutionSet.insert(node);
	}
	std::vector<TEdge> solgraph;
	for(auto edge : edges_)
	{
		if((solutionSet.find(edge.src) != solutionSet.end()) && (solutionSet.find(edge.dst) != solutionSet.end()))
		{
			solgraph.push_back(edge);
		}
	}
	sort(solgraph.begin(), solgraph.end());
	std::tuple<const unsigned int,const  unsigned int,const  unsigned int> lims{0, 0, 0};
	double est = ProcessSampleOrder(0, delta, 2, 0, 0, 0, solgraph.size(),std::cref(lims),
								Vm_, edgesM_, std::cref(solgraph), std::cref(motifprops), std::ref(stats), std::ref(node_insts));
	double top = est/(1.0*solutionSet.size());
	int precision = std::numeric_limits<double>::digits10;
	std::cout << "Solution value: " << std::setprecision(precision) << top << " insts: " << est << '\n';

	return cM;
}

double TempMotifDense::TwoDeltaPatch(int delta, double size, int threads)
{
	std::vector<Partition> partitions = ComputeSplitsAndPatch(delta, size);
	
	ctpl::thread_pool p(threads);
	std::vector<std::future<double>> futures;
	int iterations = (int) partitions.size();
	RunningTimeStats stats = RunningTimeStats();
	MyStructure node_insts(N_);
	for(int i=0; i < iterations; i++)
	{
      	futures.push_back(p.push(ProcessIntervalOrder, 6, delta, i, std::cref(partitions),
							Vm_, edgesM_, std::cref(edges_), std::cref(motifprops), std::ref(stats), std::ref(node_insts)));
	}

	double cM = 0;
	for (int i = 0; i < (int) futures.size(); i++) 
	{
		double part_count = futures[i].get();
        cM += part_count; 
	}

	return cM;
}

double TempMotifDense::exactDense(int delta, double decay)
{
	

	//std::cout << "Funciton called!" << std::endl;

	std::vector<Partition> partitions = ComputeSplitsAndPatch(delta, 2);

	MyStructure node_insts(N_);
	if(decay <1)
		node_insts.setDecay(decay);
	std::cout << "decay: " << node_insts.getDecay() << '\n';
	RunningTimeStats stats = RunningTimeStats();
	auto init = get_wall_time();
	int iterations = (int) partitions.size();
	double cM = 0;
	for(int i=0; i < iterations; i++)
	{
      	cM += ProcessIntervalOrder(i, 8, delta, i, std::cref(partitions),
							Vm_, edgesM_, std::cref(edges_), std::cref(motifprops), std::ref(stats), std::ref(node_insts));
	}
	std::unordered_map<std::string, double> subcounts = node_insts.getSubgraphCounts();
	std::cout << "Runing time to compute motifs: " << get_wall_time() - init << std::fixed << '\n';
	
	// Build flow graph
	// 1. Add Nodes
	init = get_wall_time();
	typedef Graph<double, double ,double> GraphType;
	GraphType *g = new GraphType(/*estimated # of nodes*/ N_ + subcounts.size() , /*estimated # of edges*/ subcounts.size()*Vm_); 
	std::cout << "Flow-Network will have: " << N_ + subcounts.size() << " nodes and edges: " << subcounts.size()*Vm_ << std::endl;
	double lower = 0.0;
	double upper = cM;
	double zeta {(upper+lower)/2.0};
	std::cout << "zeta: " << zeta << '\n';
	long long int nodes_added = 0;
	for(int i=0; i < N_; i++)
	{
		auto nodeid = g -> add_node(); //Node i
		g -> add_tweights( nodeid,   /* capacities */  0.0, zeta );
		nodes_added++;
	}
	//std::cout << "First nodes added" << std::endl;

	
	double edges_added = 0.0;
	for(auto& el : subcounts)
	{
		g -> add_node(); //Node corresponds to subgraph
		g -> add_tweights( nodes_added,   /* capacities */  static_cast<double>(el.second), 0.0 );
		std::istringstream iss(el.first);
		std::vector<std::string> results(std::istream_iterator<std::string>{iss},std::istream_iterator<std::string>());
		for(auto& nn : results)
		{
			int node = stoi(nn);
			g -> add_edge( nodes_added, node,    /* capacities */ std::numeric_limits<double>::max(), 0.0 );
			edges_added +=1.0;
		}
		nodes_added++;

	}
	std::cout << "Runing time to add variables on flow network: " << get_wall_time() - init << std::fixed << '\n';
	//std::cout << "Second nodes added" << std::endl;
	std::unordered_set<int> solution;
	double variation = 0;
	std::cout << "Starting with zeta: " << zeta << std::endl;
	std::cout << "upper is: " << upper << " lower is: " << lower << std::endl;
	std::cout << "Nodes added in the network: " << nodes_added << " edges added: " << edges_added << std::endl;
	double threshold = (N_*(N_-1.0));
	//std::cout << " nodes are: " << N_ << " and: " << N_-1 << " threshold nodes: " << threshold << std::endl;
	init = get_wall_time();
	while((threshold*(upper-lower)) >= (1.0))
	{
		double flow = g -> maxflow();
		int sizeT =0;
		std::unordered_set<int> newsolution;
		for(int i=0; i < N_; i++)
		{
			if (g->what_segment(i) == 0)
			{
				newsolution.insert(i);
			}
			else
				sizeT++;
		}
		if(sizeT == N_)
		{
			upper = zeta;
		}
		else
		{
			lower = zeta;
			solution = newsolution;
		}
		double oldz = zeta;
		zeta = (upper+lower)/2.0;
		variation = zeta - oldz;
		double update{-oldz};
			
		for(int i=0; i < N_; i++)
		{
			g -> add_tweights( i,   /* capacities */  0.0, variation );
		}

	}
	std::cout << "Runing time to compute exact solution: " << get_wall_time() - init << std::fixed << '\n';
	std::vector<TEdge> solgraph;
	std::unordered_set<long long int> solutionSet;
	this->optimalSolution.clear();
	for(auto& node: solution)
	{
		std::cout << "Node: " << node << '\n';
		solutionSet.insert(node);
		this->optimalSolution.push_back(node);
	}
	for(auto edge : edges_)
	{
		if((solutionSet.find(edge.src) != solutionSet.end()) && (solutionSet.find(edge.dst) != solutionSet.end()))
			solgraph.push_back(edge);
	}
	std::sort(solgraph.begin(), solgraph.end());
	std::tuple<const unsigned int,const  unsigned int,const  unsigned int> lims{0, 0, 0};
	double est = ProcessSampleOrder(0, delta, 2, 0, 0, 0, solgraph.size(),std::cref(lims),
								Vm_, edgesM_, std::cref(solgraph), std::cref(motifprops), std::ref(stats), std::ref(node_insts));
	std::cout << "Solution value: " << est/solutionSet.size() << " insts: " << est << " size: " << solutionSet.size() << '\n';
	return cM;
}

std::vector<Partition> TempMotifDense::ComputeSplitsAndPatch(int delta, double size)
{
	bool done = false;
	unsigned int prev = 0;
	std::vector<Partition> parts;
	while(!done)
	{
		
		unsigned int begin = prev;
		long long tbegin = edges_[begin].tim;
		long long ttarg = ceil(tbegin + delta*size);
		TEdge target = {0,0,ttarg, 0};
		unsigned int last = std::upper_bound(edges_.begin(), edges_.end(), target) - edges_.begin();
		long long tlast = edges_[last-1].tim;
		if((last-prev) >= l_)
		{
			Partition part = {.start=begin, .tStart=tbegin, .end=last, .tEnd=tlast, .patch=false};
			parts.push_back(part);
		}
		else if((last-prev) < l_)
		{
				;
		}
		
		if(last == edges_.size())
		{
			done = true;
			continue;
		}

		long long tlowdelta = ceil(ttarg - delta);
		long long tupdelta = ceil(ttarg + delta);
		target.tim = tlowdelta;
		unsigned int idxlow = std::lower_bound(edges_.begin(), edges_.end(), target) - edges_.begin();
		target.tim = tupdelta;
		unsigned int idxup = std::upper_bound(edges_.begin(), edges_.end(), target) - edges_.begin();
		if((idxup-idxlow) < l_)
				;
		else
		{
			Partition patch = {.start=idxlow, .tStart=edges_[idxlow].tim, .end=idxup, .tEnd=edges_[idxup-1].tim, .patch=true};
			parts.push_back(patch);
		}

		prev = last;
		if(prev ==  edges_.size())
			done = true;
	}

	return parts;
}

std::vector<Partition> TempMotifDense::ComputeSplitsAndPatch(std::vector<TEdge>& myedges, int delta, double size)
{
	bool done = false;
	unsigned int prev = 0;
	std::vector<Partition> parts;
	while(!done)
	{
		
		unsigned int begin = prev;
		long long tbegin = myedges[begin].tim;
		long long ttarg = ceil(tbegin + delta*size);
		TEdge target = {0,0,ttarg, 0};
		unsigned int last = std::upper_bound(myedges.begin(), myedges.end(), target) - myedges.begin();
		long long tlast = myedges[last-1].tim;
		if((last-prev) >= l_)
		{
			Partition part = {.start=begin, .tStart=tbegin, .end=last, .tEnd=tlast, .patch=false};
			parts.push_back(part);
		}
		else if((last-prev) < l_)
		{
				;
		}
		
		if(last == myedges.size())
		{
			done = true;
			continue;
		}

		long long tlowdelta = ceil(ttarg - delta);
		long long tupdelta = ceil(ttarg + delta);
		target.tim = tlowdelta;
		unsigned int idxlow = std::lower_bound(myedges.begin(), myedges.end(), target) - myedges.begin();
		target.tim = tupdelta;
		unsigned int idxup = std::upper_bound(myedges.begin(), myedges.end(), target) - myedges.begin();
		if((idxup-idxlow) < l_)
				;
		else
		{
			Partition patch = {.start=idxlow, .tStart=myedges[idxlow].tim, .end=idxup, .tEnd=myedges[idxup-1].tim, .patch=true};
			parts.push_back(patch);
		}

		prev = last;
		if(prev ==  myedges.size())
			done = true;
	}

	return parts;
}
