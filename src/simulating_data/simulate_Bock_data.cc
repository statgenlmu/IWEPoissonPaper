
#include <iostream>
#include <vector>
#include <random>
#include <math.h> 
#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <map>

#include <ctime>


#include "back_tree.hh"

#include "mutationclass.hh"

#include "parameters.hh"

#include "cpgislandcontainer.hh"

#include "functionobjects.hh"

using namespace std;


mt19937 re;


int main(int argc , char * argv[])
{


 
  unsigned int seed=clock();
  unsigned int second_seed=clock();
  seed*=second_seed;
  
  string cmp;
  for (int i=1;i<argc;++i)
    {
      cmp=argv[i];
    if (cmp=="-seed")
      {
	
	seed=atoi(argv[i+1]);}
    if(cmp=="-reuse_seed")
      {
	string seedstr;
	cout << "reuse seed" << endl;
	ifstream ifs ("seed.txt");
	ifs >> seedstr;
	
	seed =stoi(seedstr);
      }
       
    }

  
  re.seed(seed);

 
  ofstream seed_stream;
  
  seed_stream.open("seed.txt");
  seed_stream<< seed;
  seed_stream.close();

  uniform_real_distribution<double> sigma(0,1);

  vector <back_tree*> node_container;
 
  Parameter_Struct params;
  set_parameters_from_command_line(params, node_container, argc ,  argv);
  back_tree HSC;
  back_tree MPP1;
  back_tree MPP2;
  back_tree CLP;
  back_tree CD4;
  back_tree CD8;
  back_tree B_cell;
  back_tree CMP;
  back_tree MEP;
  back_tree Eryth;
  back_tree GMP;
  back_tree Granu;
  back_tree Mono;
	
  if (node_container.size()==0)
    {
    
      HSC.set_name("HSC");
      node_container.push_back(&HSC);
      HSC.t_coal=6;
      HSC.branch_length=exp(-5);

     
      MPP1.set_name("MPP1");
      HSC.add_to_children(&MPP1);
      node_container.push_back(&MPP1);
      MPP1.t_coal=5;
      MPP1.branch_length=exp(-5);

   
      MPP2.set_name("MPP2");
      MPP1.add_to_children(&MPP2);
      node_container.push_back(&MPP2);
      MPP2.t_coal=4;
      MPP2.branch_length=exp(-5);

      
      CLP.set_name("CLP");
      MPP2.add_to_children(&CLP);
      node_container.push_back(&CLP);
      CLP.t_coal=3;
      CLP.branch_length=exp(-5);
  

    
      CD4.set_name("CD4");
      CLP.add_to_children(&CD4);
      node_container.push_back(&CD4);
      CD4.t_coal=2;
      CD4.branch_length=exp(-5);

    
      CD8.set_name("CD8");
      CLP.add_to_children(&CD8);
      node_container.push_back(&CD8);
      CD8.t_coal=2;
      CD8.branch_length=exp(-5);
  
    
      B_cell.set_name("B_cell");
      CLP.add_to_children(&B_cell);
      node_container.push_back(&B_cell);
      B_cell.t_coal=2;
      B_cell.branch_length=exp(-5);
  
     
      CMP.set_name("CMP");
      MPP2.add_to_children(&CMP);
      node_container.push_back(&CMP);
      CMP.t_coal=3;
      CMP.branch_length=exp(-5);
  
    
      MEP.set_name("MEP");
      CMP.add_to_children(&MEP);
      node_container.push_back(&MEP);
      MEP.t_coal=2;
      MEP.branch_length=exp(-5);
  
      Eryth.set_name("Eryth");
      MEP.add_to_children(&Eryth);
      node_container.push_back(&Eryth);
      Eryth.t_coal=1;
      Eryth.branch_length=exp(-5);
  
     
      GMP.set_name("GMP");
      CMP.add_to_children(&GMP);
      node_container.push_back(&GMP);
      GMP.t_coal=2;
      GMP.branch_length=exp(-5);
  
    
      Granu.set_name("Gran");
      GMP.add_to_children(&Granu);
      node_container.push_back(&Granu);
      Granu.t_coal=1;
      Granu.branch_length=exp(-5);

    
      Mono.set_name("Mono");
      GMP.add_to_children(&Mono);
      node_container.push_back(&Mono);
      Mono.t_coal=1;
      Mono.branch_length=exp(-5);

    }
	  

  
  vector<double> log_branch_lengths={6,-5,-4,-0.4,-2,-3,-6,-5,-2,-3,-2,-4,0};
  if (params.log_branches.size()>0)
    {
      log_branch_lengths=params.log_branches;
    }
      
  for (unsigned i=0; i<node_container.size(); ++i)
    node_container[i]->branch_length=exp(log_branch_lengths[i]);
  
  vector<cpg_island> cpg_islands;
 
  
  DistObjects dist;

  LogNormalProposal proposal_for_branch_length(1);
  LogNormalPriorDensity length_prior;
  dist.proposal_for_branch_length=&proposal_for_branch_length;
  dist.length_prior=&length_prior;

  LogNormalProposal proposal_for_gamma_shape(1);
  dist.gamma_rate_construction= GammaRateConstruction(exp(params.log_shape),3);
  auto rate_acc=dist.gamma_rate_construction.get_rate_acc(re);
  LogNormalPriorDensity gamma_shape_prior(1,2);
  dist.proposal_for_gamma_shape=&proposal_for_gamma_shape;
  dist.gamma_shape_prior=&gamma_shape_prior;

  LogNormalProposal pcrateproposal(1);
  LogNormalPriorDensity pcrate_prior;

  dist.pcrateproposal=&pcrateproposal;
  dist.pcrate_prior=&pcrate_prior;

  vector<vector<double>> meth_freq;
  meth_freq.push_back({0.8,0.15,0.05}); 
  double pcrate=exp(params.log_pc_rate);

  ofstream target_param_stream("target_params.target");
  target_param_stream<<"log shape "<<log(dist.gamma_rate_construction.get_alpha())<<endl;
  target_param_stream<<"rate acc ";
  for (auto &i : rate_acc)
    target_param_stream<<i<<" ";
  target_param_stream<<endl;
   target_param_stream<<"log branch lengths ";
  for (auto &i : log_branch_lengths)
    target_param_stream<<i<<" ";
  target_param_stream<<endl;
  vector<unsigned> island_sizes={43,754,2703,321,5430,43,32,31,410,430,120,210,2,60,800,5,430,500,850,4,34,32,2,7426,3074,343,2,5,8,4,23,32};

  if (params.island_sizes.size()>0)
    {
      island_sizes=params.island_sizes;
    }
 
  
  
  vector<vector<int>> site_factors;
  uniform_int_distribution<int> which_factor_dist(0,rate_acc.size()-1);

  double p_inv=params.p_inv;
  target_param_stream<<"p_inv "<<p_inv<<endl;
  target_param_stream<<"log pcrate "<<log(pcrate)<<endl;

  bernoulli_distribution invariant_site(p_inv);
  
  for (const auto &i : island_sizes)
    {
      site_factors.push_back(vector<int>(i,0));
      for (unsigned j=0; j<i; j++)
	{
	  if (invariant_site(re))
	    site_factors.back()[j]=-1;
	  else
	    site_factors.back()[j]=which_factor_dist(re);
	}
    }
  node_container[0]->initialize_island_meth_changes(rate_acc, dist, pcrate, island_sizes);
  target_param_stream<<"n_meth_events ";
  for (auto node : node_container) target_param_stream<<node->get_n_meth_events()<<" ";
  target_param_stream<<endl;
  
    
  node_container[0]->simulate_methylation( island_sizes, site_factors);
  node_container[0]->print_cpg_islands();
  ofstream probability_change_stream;
  probability_change_stream.open("target_p_change_info.txt");
  node_container[0]->print_transition_steps( probability_change_stream);
}
