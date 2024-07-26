#include "temporalmotifdense.h"
#include <nlohmann/json.hpp>

#include <omp.h>

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <dirent.h>
#include <mutex>

//For memmory stuff
#include <stdlib.h>

#define MEM_ENV_VAR "MAXMEM_MB"

void setmemlimit();
using json = nlohmann::json;
//

void setmemlimit()
{
	struct rlimit memlimit;
	long bytes;

	if(getenv(MEM_ENV_VAR)!=NULL)
	{
		bytes = atol(getenv(MEM_ENV_VAR))*(1024*1024);
		memlimit.rlim_cur = bytes;
		memlimit.rlim_max = bytes;
		setrlimit(RLIMIT_AS, &memlimit);
		//setrlimit(RLIMIT_STACK, &memlimit);
	}
}

double get_wall_time(){

    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int main(int argc, char* argv[]) {
	const rlim_t kStackSize = RLIM_INFINITY;   // min stack size = 16 MB
    struct rlimit rl;
    int status;

	setmemlimit();

	status = getrlimit(RLIMIT_STACK, &rl);
    if (status == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            status = setrlimit(RLIMIT_STACK, &rl);
            if (status != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", status);
            }
        }
    }
	

  Env = TEnv(argc, argv, TNotify::StdNotify);

  const TStr temporal_graph_filename =
    Env.GetIfArgPrefixStr("-i:", "simple-example.txt",
			  "Input directed temporal graph file");
  const TStr motif_graph_filename =
    Env.GetIfArgPrefixStr("-m:", "simple-motif.txt",
        "Input directed motif graph file");
  const TFlt delta =
    Env.GetIfArgPrefixFlt("-delta:", 4096, "Time window delta");
  const TFlt c =
    Env.GetIfArgPrefixFlt("-c:", 30, "Window size multiplier");
  const TFlt decay =
    Env.GetIfArgPrefixFlt("-lambda:", 2.0, "Exponential decay over the network usually << 1");
  const TInt type =
    Env.GetIfArgPrefixInt("-type:", 0, "Which alg?");
  const int samp=
    Env.GetIfArgPrefixInt("-samp:", 100, "number of samples to be used in Embedded sampling algorithm");
  const int hyb=
    Env.GetIfArgPrefixInt("-ithyb:", 2, "number of iterations for hybrid peeling algorithm");
  const TStr motif_folder =
    Env.GetIfArgPrefixStr("-fol:", "../motifs/",
        "Folder of the motifs to be used");
    const TFlt epsilon =
    Env.GetIfArgPrefixFlt("-epsilon:", 0.01, "Threshold for batch peel");
  
  int num_threads = {1};
  bool ises { 0 };
  double st = get_wall_time(); 
  TempMotifDense tmc(temporal_graph_filename.CStr(), motif_graph_filename.CStr(), false, ises);

  Env.PrepArgs(TStr::Fmt("Temporalmotifs. build: %s, %s. Time: %s",
       __TIME__, __DATE__, TExeTm::GetCurTm()));  
	printf("Time to load the dataset: %lfs\n\n", get_wall_time() - st);
    	srand(time(nullptr));

    double init;
    double end;
    double frac = 0;
    int samples = 0;
    double a=0;
	
    // COUNTING/ENUMERATING TEMPORAL MOTIF INSTANCES WITH PARALLEL 2-DELTA PATCH
	if(type == 0)
	{
		init = get_wall_time();
		double count = tmc.TwoDeltaPatch(delta, 2, num_threads);
		std::cout <<"--Exact counting/enum delta-inst-- run time : "<< get_wall_time() - init << "s" << std::endl;
		std::cout.precision(2);
		std::cout <<"--Exact counting/enum delta-inst-- instances counted: "<< std::fixed << count << std::endl;
	}
    // 1/k-APPROX PEELING GREEDY ALGORITHM
	if(type == 1)
	{
		init = get_wall_time();
		double count = tmc.kApproxDenseMotif(delta, decay);
		std::cout <<"--k-Approx-- run time : "<< get_wall_time() - init << "s" << std::endl;
		std::cout.precision(2);
		std::cout <<"--k-Approx-- instances counted: "<< std::fixed << count << std::endl;
	}
    // 1/(k(1+\xi))-APPROX BATCH PEELING ALGORTIHM
	if(type == 2)
	{
		init = get_wall_time();
		double count = tmc.kApproxDenseMotifBatchPeel(delta, epsilon, decay);
		std::cout <<"--Batch-peel-- run time : "<< get_wall_time() - init << "s" << std::endl;
		std::cout.precision(2);
		std::cout <<"--Batch-peel-- instances counted: "<< std::fixed << count << std::endl;
	}
    // ALDENTE-ProbPeel
	if(type == 3)
	{
		long long int seeds[5]{1234, 48593, 293934, 38484750, 274540193}; // Execution on five different random seeds
		double count = 0;
		for(int i=0; i < 5; i++)
		{
			init = get_wall_time();
			count = tmc.kApproxDenseMotifBatchPeelSampling(delta, epsilon, samp, seeds[i], decay);
			std::cout <<"--ALDENTE-ProbPeel-- run time : "<< get_wall_time() - init << "s" << std::endl;
			std::cout.precision(2);
			std::cout <<"--ALDENTE-ProbPeel-- instances counted: "<< std::fixed << count << std::endl;
		}
	}
    // EXACT-FLOW BASED ALGORITHM
	if(type == 4)
	{
		init = get_wall_time();
		double count = tmc.exactDense(delta, decay);
		std::cout <<"--Exact-FlowBased-- run time : "<< get_wall_time() - init << "s" << std::endl;
		std::cout.precision(2);
		std::cout <<"--Exact-FlowBased-- instances counted: "<< std::fixed << count << std::endl;
		auto optSol = tmc.getOptSol();
		json j = {
		{"Nodes", optSol}};
		std::cout << "JSONSOL " << j << '\n';
	}
    // ALDENTE-HybridPeel
	if(type == 5)
	{
		long long int seeds[10]{1234, 48593, 293934, 38484750, 274540193, 839, 2739494, 52547839, 19209, 4995}; 
		double count = 0;
		for(int i=0; i < 5; i++)// Execution on five different random seeds
		{
			init = get_wall_time();
			count = tmc.kApproxDenseMotifBatchPeelSamplingHybrid(delta, epsilon, samp, seeds[i], decay, hyb);
			std::cout <<"--2DeltaPatchTime-- run time : "<< get_wall_time() - init << "s" << std::endl;
			std::cout.precision(2);
			std::cout <<"--2DeltaPatchSol-- instances counted: "<< std::fixed << count << std::endl;
		}
	}
    // TEST PROB-PEEL CONVERGENCE SAMPLE SIZE, see paper
	if(type == 6)
	{
		long long int seeds[5]{1234, 48593, 293934, 3884750, 2540193};
		double count = 0;
		std::vector<double> mysamps;
		int uppersamp{20*samp};
		int step=(uppersamp-samp)/30;
		for(int ii{samp}; ii<=uppersamp; ii+=step)
			mysamps.push_back(ii);
		int j=0;
		for(auto ss : mysamps)
		{
			for(int i=0; i < 5; i++)
			{
				init = get_wall_time();
				count = tmc.kApproxDenseMotifBatchPeelSampling(delta, epsilon, ss, seeds[i]+(j*10));
				std::cout <<"--2DeltaPatchTime-- run time : "<< get_wall_time() - init << "s" << std::endl;
				std::cout.precision(2);
				std::cout <<"--2DeltaPatchSol-- instances counted: "<< std::fixed << count << std::endl;
			}
		}
	}
    // TEST SENSITIVITY FOR BATCH-PEEL ACCORDING TO VARIOUS XI VALUES, SEE PAPER
	if(type == 7)
	{
		double eps_max = epsilon;
		double eps_init = 0.001;
		double step = 0.005;

		std::vector<double> epses;
		while(eps_init <= eps_max)
		{
			epses.push_back(eps_init);
			eps_init += step;
		}
		epses.push_back(eps_max);
		std::cout.precision(5);
		for(auto eps : epses)
		{
				std::cout << "EPSVAL: " << std::setprecision(5)<< eps << '\n';
				init = get_wall_time();
				double count = tmc.kApproxDenseMotifBatchPeel(delta, eps);
				std::cout <<"--2DeltaPatchTime-- run time : "<< get_wall_time() - init << "s" << std::endl;
				std::cout.precision(2);
				std::cout <<"--2DeltaPatchSol-- instances counted: "<< std::fixed << count << std::endl;
		}
	}
	
  return 0;
}
