

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

#include "MCMC_step.hh"

#include "functionobjects.hh"

using namespace std;


mt19937 re;


int main(int argc , char * argv[])
{

  clock_t clock_begin = clock();
 
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
 
  back_tree HSC;
  HSC.set_name("HSC");
  node_container.push_back(&HSC);
 

  back_tree MPP1;
  MPP1.set_name("MPP1");
  HSC.add_to_children(&MPP1);
  node_container.push_back(&MPP1);
 
  

  back_tree MPP2;
  MPP2.set_name("MPP2");
  MPP1.add_to_children(&MPP2);
  node_container.push_back(&MPP2);
  

  back_tree CLP;
  CLP.set_name("CLP");
  MPP2.add_to_children(&CLP);
  node_container.push_back(&CLP);
 
  
  

  back_tree CD4;
  CD4.set_name("CD4");
  CLP.add_to_children(&CD4);
  node_container.push_back(&CD4);
 
 

  back_tree CD8;
  CD8.set_name("CD8");
  CLP.add_to_children(&CD8);
  node_container.push_back(&CD8);
  
  
  back_tree B_cell;
  B_cell.set_name("B_cell");
  CLP.add_to_children(&B_cell);
  node_container.push_back(&B_cell);
 
  
  
  back_tree CMP;
  CMP.set_name("CMP");
  MPP2.add_to_children(&CMP);
  node_container.push_back(&CMP);
 
 
  
  back_tree MEP;
  MEP.set_name("MEP");
  CMP.add_to_children(&MEP);
  node_container.push_back(&MEP);
 
 
  
  back_tree Eryth;
  Eryth.set_name("Eryth");
  MEP.add_to_children(&Eryth);
  node_container.push_back(&Eryth);
  
  
  back_tree GMP;
  GMP.set_name("GMP");
  CMP.add_to_children(&GMP);
  node_container.push_back(&GMP);
 
  
  back_tree Granu;
  Granu.set_name("Gran");
  GMP.add_to_children(&Granu);
  node_container.push_back(&Granu);


  back_tree Mono;
  Mono.set_name("Mono");
  GMP.add_to_children(&Mono);
  node_container.push_back(&Mono);


  vector<back_tree*> nc;

  set_starting_tree(argc,argv,nc);
  if(nc.size()>0) node_container=nc;


  
  vector<double> log_branch_lengths;
  normal_distribution<double> branch_prior_dist(-2,1);
  for (unsigned i=0; i<node_container.size(); ++i)
    log_branch_lengths.push_back(-2);    
  for (unsigned i=0; i<node_container.size(); ++i)
    node_container[i]->branch_length=exp(log_branch_lengths[i]);

  vector<cpg_island> cpg_islands;
  /*
   *I/O reads values from system and writes them into param
   */
  set_methylation_states_from_files( node_container);
 

  DistObjects dist;

  BoundedNormalProposal ip_prop_small(0.0,1.0, 0.1);
  BoundedNormalProposal ip_prop_large(0.0,1.0, 0.6);
  MultiProposal ip_prop({&ip_prop_small,&ip_prop_large});
  
  dist.invariant_probability_poposal=&ip_prop;
  UniformPriorDensity ip_prior;
  dist.invariant_probability_prior=&ip_prior;
  
  LogNormalProposal proposal_for_branch_length_small(0.2);
  LogNormalProposal proposal_for_branch_length_large(1.5);
  MultiProposal proposal_for_branch_length({&proposal_for_branch_length_small,&proposal_for_branch_length_large});
  
  
  LogNormalPriorDensity length_prior(-2,1);

  
  dist.proposal_for_branch_length=&proposal_for_branch_length;
  dist.length_prior=&length_prior;


  LogNormalProposal proposal_for_gamma_shape_small(0.2);
  LogNormalProposal proposal_for_gamma_shape_large(3.5);
  MultiProposal proposal_for_gamma_shape({&proposal_for_gamma_shape_small,&proposal_for_gamma_shape_large});

  double prob_inv=0.5;
  double log_shape=2;
  double log_pc_rate=-2;
  set_start_params(argc, argv, prob_inv, log_shape,log_pc_rate);
  double pcrate=exp(log_pc_rate);
  
  dist.gamma_rate_construction= GammaRateConstruction(exp(log_shape),3);
  auto rate_acc=dist.gamma_rate_construction.get_rate_acc();
  LogNormalPriorDensity gamma_shape_prior(2,1);
  dist.proposal_for_gamma_shape=&proposal_for_gamma_shape;
  dist.gamma_shape_prior=&gamma_shape_prior;

  LogNormalProposal pcrateproposal_small(0.2);
  LogNormalProposal pcrateproposal_large(1.5);
  MultiProposal pcrateproposal({&pcrateproposal_small,&pcrateproposal_large});
  LogNormalPriorDensity pcrate_prior(-2,1);

  dist.pcrateproposal=&pcrateproposal;
  dist.pcrate_prior=&pcrate_prior;

  vector<vector<vector<double>>> meth_freq=node_container[0]->get_empirical_meth_freq();
 

  
  for (auto i: node_container) i->initialize_island_meth_changes(rate_acc,   meth_freq);

  if(use_transition_solution(argc, argv)) read_probability_change_from_file(node_container, rate_acc);
  
  
  ofstream rate_acc_stream("rate_acc.param");
  ofstream n_iter_stream("n_iter.param");
  ofstream branch_length_stream("log_branch_length.param");
  ofstream pcrate_stream("log_pcrate.param");
  ofstream computing_time;
  computing_time.open("computing_time.txt");
  ofstream log_shape_stream("log_shape.param");
  ofstream n_pc_stream("npc.param");
  ofstream likelihoodstream("likelihood.param");
  ofstream invpstream("invp.param");
  ofstream n_affected_island("n_affected_island.param");
  ofstream probability_change_stream;
  probability_change_stream.open("p_change_info.param");

  ofstream sample_params;
  sample_params.open("sample_params.param");
  auto lhd=node_container[0]->llhd_of_transition_with_fixed_ends_whole_tree_all_islands(rate_acc.size(),prob_inv);
  auto old_lhd=node_container[0]->llhd_of_transition_with_fixed_ends_whole_tree_all_islands(rate_acc.size(),prob_inv);
  int choice=0;
  vector<bool> which_step(6, true);
  which_step[1]=true;
  which_step = disable_modes(argc,argv,which_step);
  double logMeth=0;
 
  for (int i=0; i< 5e7; i++)
    {
   
      

      if (i%1000==0)
	{
	  old_lhd=lhd;
	  lhd=node_container[0]->llhd_of_transition_with_fixed_ends_whole_tree_all_islands(rate_acc.size(),prob_inv);
       
	  for (auto &i: rate_acc) rate_acc_stream<<i<<" ";
	  rate_acc_stream<<endl;
	  n_iter_stream<<i<<endl;
	  for (auto &node: node_container) branch_length_stream<<log(node->branch_length)<<" ";
	  for (auto &node: node_container) n_pc_stream<<node->get_n_meth_events()<<" ";
	  for (auto &node: node_container) n_affected_island<<node->get_n_affected_islands()<<" ";
	  n_affected_island<<endl;
	  branch_length_stream<<endl;
	  n_pc_stream<<endl;
	  pcrate_stream<<log(pcrate)<<endl;
	  clock_t clock_end = clock();
	  double elapsed_time=double(clock_end - clock_begin) / CLOCKS_PER_SEC;
	  computing_time<<elapsed_time<<endl;
	  node_container[0]->print_transition_steps( probability_change_stream);
	  probability_change_stream<<endl;

	  log_shape_stream<<log(dist.gamma_rate_construction.get_alpha())<<endl;

	  logMeth=0;
	  int n_islands=node_container[0]->get_n_islands();
	  int n_meth_events=0;
	  double tree_length=0;
	  for (auto node: node_container)
	    {
	      if (node->parent!=0)
		{
		  n_meth_events+=node->get_n_meth_events();
		  tree_length+=node->branch_length;
		}
	    }
	  logMeth=log10poisdensity(pcrate*n_islands*tree_length, n_meth_events);

	  likelihoodstream<<lhd.first<<"|"<<lhd.second<<" "<<logMeth<<endl;
	  invpstream<<prob_inv<<endl;


	  switch(choice)
	    {
	    case 0:
	      sample_params<<"sample branch length "<<endl;
    	      break;
	    case 1:
	      sample_params<<"sample rate accell "<<endl;     
	      break;
	    case 2:
	      sample_params<<"sample pc rate "<<endl;
    	      break;
	    case 3:
	      sample_params<<"sample pc event "<<endl;
	      break;
	    case 4:
	      sample_params<<"Sample new invariant probability "<<endl;
	      break;
	    case 5:
	      sample_params << "Sample new beginning "<<endl;
	      break;
	    default: break;
	    }
	}


	
	
      choice=MCMC_step(node_container, dist, rate_acc, pcrate,prob_inv,which_step,re);

    }

  sample_params.close();
  rate_acc_stream.close();
  n_iter_stream.close();
  computing_time.close();
  log_shape_stream.close();
  probability_change_stream.close();
  invpstream.close();

 
  cout<<"Finished MCMC "<<endl;

  
}
