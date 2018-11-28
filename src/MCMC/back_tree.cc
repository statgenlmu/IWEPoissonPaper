
/* Copyright (C) 2018 Konrad Grosser and Dirk Metzler */
/* This file is part of IWEinference.
   
   IWEinference is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   IWEinference is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
    
   You should have received a copy of the GNU General Public License
   along with IWEinference.  If not, see <https://www.gnu.org/licenses/>. */

#include "mutationclass.hh"
#include "back_tree.hh"
#include <iostream>
#include <vector>
#include <deque>
#include <math.h>
#include <random>
#include <algorithm>
#include <map>
#include <memory>
#include <boost/math/distributions/beta.hpp>
#include <fstream>
#include "felsenstein_pruning.hh"
#include "functionobjects.hh"
#include "MCMC_step.hh"
 
#include <boost/algorithm/string.hpp>

using namespace std;
using boost::math::beta_distribution;

class Event;


back_tree::back_tree()
{
  parent=nullptr;
}

back_tree::~back_tree()
{
}



/*
 * This function parses files of the format CPG_ISLANDS_celltype.txt for the celltype name 
 * and initializes states for MCMC compuation.
 */

void back_tree::set_methylation_states_from_files()
{
  string file_name="CPG_ISLANDS_";
  file_name+=name+".txt";
  ifstream file(file_name);
  string current_cpgi;
  while ( getline(file, current_cpgi )) {
      if (current_cpgi=="") cerr<<"Warning: Empty line in CPG_ISLAND file"<<endl;
      else {
	  cpg_island cgi(current_cpgi);
	  //constructor of cpg_island takes a string to store all relevant sites, code in mutationclass.cc
	  cpg_container.push_back(cgi);
      }
  }
}


/*
 * This function parses files of the format CPG_ISLANDS_celltype.txt for the celltype name 
 * and initializes states for MCMC compuation.
 * file path is the path the files are in, for example ../testrun/
 */

void back_tree::set_methylation_states_from_files(const string & file_path)
{
  string file_name=file_path+"CPG_ISLANDS_";
  file_name+=name+".txt";
  ifstream file(file_name);
  string current_cpgi;
  vector <cpg_island> cpgis;
  
  while ( getline(file, current_cpgi ))
    {
      if (current_cpgi=="") cerr<<"Warning: Empty line in CPG_ISLAND file"<<endl;
      else
	{
	  cpg_island cgi(current_cpgi);
	  //constructor of cpg_island takes a string to store all relevant sites, code in mutationclass.cc
	  cpgis.push_back(cgi);
	}
    }
  cpg_container=cpgis;
}




/*
 * The following function takes a pair of probability and scale factor a 
 * and probability and scalefactor B and returns their product with appropriate scale factor. 
 * Scalefactor repreresents the power of -10 that needs to be multiplied to the probability to 
 * get the actual value.
 */


const pair <double, int> llhdmulti(const pair <double, int>  A, const pair <double, int>  B)
{
  double prod=A.first*B.first;
  if (prod>1.0000000001) cerr<< "llhdmulti Warning: likelihoods should never be bigger than 1 in the first factor "<<endl; 
  unsigned scale=A.second+B.second;
  if (prod<1e-100)
    {
      prod*=1e+100; 
      ++scale;
    }
  return {prod,scale};
}



/*
 * Function for determining the likelihood of transiton of a specific site and rate 
 * with fixed endpoint at the current node and begninning at the node above.
 *
 * Input:
 * Index of rate r
 * index of site s
 * index of island i
 *
 * output: pair of probability and scale factor
 */
 

pair<double, int> back_tree::llhd_of_single_site_transition_with_fixed_ends(const int i, const int r, const int s)
{
  pair<double, int> output={1.0,0};
  if (parent==0) return output; //transition at the root node is not well defined
  int x=parent->cpg_container[i].states[s];
  if (x!=-1) //case of state bein available
    {
      output.first=transition[i]->get_part_L(s,r)[x];
      output.second=transition[i]->get_scalefactor();
    }
  return output;
}



/*
 * The following function checks wheter a site changes along the tree or stays invariant.
 * it returns true if there are no detectable changes, false otherwise.
 *
 * The input argument state is the state the site was in in cells that have been checked so far. It is -1 at initialization or if states
 * so far have been NAs. Else 0,1,2 corresponding to U,P,M.
 * i stands for the index of the island, s for the index of the site.
 */


bool back_tree::possible_inv(const int i, const int s, int & state)
{
  if (state==-1) state=cpg_container[i].states[s];
  if (state!=cpg_container[i].states[s]&&cpg_container[i].states[s]!=-1)return false;
  for (auto node: children)if (!node->possible_inv(i,s,state)) return false;
  return true;
}



/*
 * Function for determining the likelihood of transition of a specific 
 * site and rate along the whole tree with a given rate. Should be called at the 
 * tree root.
 *
 * Input:
 * Index of rate r
 * index of site s
 * index of island i
 *
 * output: pair of probability and scale factor
 */

const pair<double, int> back_tree::llhd_of_single_site_transition_whole_tree(const int i, const int r, const int s)
{
  pair<double,int> output;
  if (parent==0) output=make_pair(1.0,0);
  else output=llhd_of_single_site_transition_with_fixed_ends(i,r,s);
  for (auto iter =children.begin(); iter!=children.end(); ++iter)
    {
      auto newlhd=(*iter)->llhd_of_single_site_transition_whole_tree(i,r,s);
      output=llhdmulti(output,newlhd);
    }
  return output;
}



/*
 * Function for determining the likelihood of transition of a specific 
 * island along the whole tree with a given rate. Should be called at the 
 * tree root.
 *
 * Input:
 * index of island i
 * probability of sites being invariant prob_inv
 * n_rates: number of rate acceleration parameters.
 *
 * output: pair of probability and scale factor
 */

const pair<double, int> back_tree::llhd_of_island_transition_whole_tree(const int i, const int n_rates, const double prob_inv)
{
  pair<double,int> output={1,0};
  for (unsigned s=0; s < cpg_container[i].states.size(); ++s)
    {
      vector<pair<double, int>> likelihoods_of_site(n_rates, {1.0,0});
      for (int r=0; r<n_rates; ++r) likelihoods_of_site[r]=llhd_of_single_site_transition_whole_tree(i, r, s);
      auto mean_lhd_of_site=meanlhd(likelihoods_of_site);
      int state=-1;
      const bool poss_inv=possible_inv(i,s,state);
      mean_lhd_of_site.first*=(1-prob_inv);
      if(poss_inv)
	{
	  mean_lhd_of_site.first=mean_lhd_of_site.first*(mean_lhd_of_site.second==0)+prob_inv;
	  mean_lhd_of_site.second=0;
	}
      output=llhdmulti(output,mean_lhd_of_site);
    }
  return output;
}


/*
 * Function for determining the combined likelihoods of transition of all sites averaged over 
 * rate factors of all islands along the whole tree. Should be called at the 
 * tree root.
 *
 * Input:
 * nrates: # of rate acceleration factors
 * prob_inv is the probability that the site is invariant
 *
 * output: pair of probability and scale factor
 */

pair <double, int>  back_tree::llhd_of_transition_with_fixed_ends_whole_tree_all_islands(const int n_rates, const double prob_inv)
{
  pair<double,int> output={1,0};
  for (unsigned i=0; i<transition.size(); i++)
    {
      for (unsigned s=0; s < cpg_container[i].states.size(); ++s)
	{
	  vector<pair<double, int>> likelihoods_of_site(n_rates, {1.0,0});
	  for (int r=0; r<n_rates; ++r) likelihoods_of_site[r]=llhd_of_single_site_transition_whole_tree(i, r, s);
	  auto mean_lhd_of_site=meanlhd(likelihoods_of_site);
	  int state=-1;
	  bool poss_inv=possible_inv(i,s,state);
	  mean_lhd_of_site.first*=(1-prob_inv); 
	  if(poss_inv)
	    {
	      mean_lhd_of_site.first=mean_lhd_of_site.first*(mean_lhd_of_site.second==0)+prob_inv;
	      mean_lhd_of_site.second=0;
	    }
	  output=llhdmulti(output,mean_lhd_of_site);
	}
    }
  return output;
}



/*
 * When in the root the depednence of root states upon methylation probabilities in the root is not captured
 * by felsenstein pruning lhds. This function calculates that lhd.
 *
 *
 * Parameters:
 * index: is the index of the island.
 * meth_freq: are the state probabilities in that island.
 *
 * Returns the lhd of the states in the island coming from the probabilities.
 *
 */
pair<double,int> back_tree::lhd_of_root_states(const int index, const vector<double> & meth_freq)
{
  pair<double,int> lhd={1,0};
  for (const auto & s:cpg_container[index].states)
    {
      if(s!=-1) lhd.first*=meth_freq[s];
      if (lhd.first<1e-100)
	{
	  lhd.first*=1e100;
	  lhd.second++;
	}
    }
  return lhd;
}




/*
 * Fucntion that stretches all time differences between this node and its parent stored in
 * the island_meth classes.
 */
void back_tree::stretch_branch(double const factor)
{
  for (const auto & i : transition) i->stretch_branch_length(factor);
}



void back_tree::set_rate_acc_in_tree(vector<double> & rate_acc)
{
  for (auto & i : transition) i->set_rate_acc(rate_acc);
  for (auto i =children.begin(); i!=children.end(); ++i) (*i)->set_rate_acc_in_tree(rate_acc);
}


const unsigned back_tree::get_n_meth_events()
{
  unsigned n=0;
  for (auto & i : transition) n+=i->get_n_intervals()-1;
  return n;
}


/*
   * Islands can have multiple events on them along a branch - it is often intersting 
   * to know wheter they merely cluster on one branch or are spread out, so I give out 
   * the the number of islands carrying an event as well.
   */
 
const unsigned back_tree::get_n_affected_islands()
{
  unsigned n=0;
  for (auto & i : transition) n+=(i->get_n_intervals()>1);
  return n;
}

const unsigned back_tree::get_n_meth_events_whole_tree()
{
  unsigned n=0;
  if(parent!=0) n=get_n_meth_events();
  for (auto & i :children) n+=i->get_n_meth_events_whole_tree();
  return n;
}

const double back_tree::get_complete_tree_length()
{
  double l=0;
  if (parent!=0) l=branch_length;
  for (auto & i :children)
    l+=i->get_complete_tree_length();
  return l;
}
  


  /* If a new methylation event happens, methylation probabilities
   * at the beginnings of daughter nodes 
   * must be adjusted iteratively, until those nodes have a methylation change. The following
   * function does that.
   *
   * Input:
   * meth_freq is the methylation probability that is inserted
   * index is the index of the island where the methylation event is inserted
   */

void back_tree::set_meth_freqs_in_island_until_change(const vector<double> & meth_freq, const unsigned index)
{
  transition[index]->set_first_meth_freq(meth_freq);
  if (transition[index]->get_n_intervals()==1)
    for (auto i =children.begin(); i!=children.end(); ++i) (*i)->set_meth_freqs_in_island_until_change(meth_freq,index);
}

const unsigned back_tree::get_n_islands() {return transition.size();}
const unsigned back_tree::get_n_meth_events_at_island(unsigned i)
{
  return transition[i]->get_n_intervals()-1;
}
 
   
void back_tree::create_new_event(back_tree* tree_root,  DistObjects & dists, const double nmethdensity, const vector<double> & rate_acc, const unsigned index, const double prob_inv, mt19937 & re)
{
    const double log1e100 = 230.2585; // ln(10^100)
 
    pair<double, int>old_lhd=tree_root->llhd_of_island_transition_whole_tree(index, rate_acc.size(), prob_inv);
  

  vector<double> prev_change_time=transition[index]->get_change_time();
  uniform_real_distribution<double> insertion_time_distribution(0, prev_change_time.back());
  const double new_time=insertion_time_distribution(re);
  //stores the change time we had previously so that we can revert if a jump is not accepted
  vector<vector<double>> prev_methyl_freq=transition[index]->get_methyl_freq();
  //stores the frequencies we had previously so that we can revert if a jump is not accepted
  auto prev_daughter=transition[index]->get_daughter();
  
  auto old_probabilities=transition[index]->get_last_meth_freq();
  auto new_probabilities=dists.dirichlet111(re);

  transition[index]->plus_one_meth_change(new_probabilities,new_time);
  transition[index]->update_all_transprobs();
  transition[index]->update_partial_llhs();

  auto last_freq=transition[index]->get_last_meth_freq();
  for (auto i =children.begin(); i!=children.end(); ++i)(*i)->set_meth_freqs_in_island_until_change(last_freq,index);

  pair<double, int>new_lhd=tree_root->llhd_of_island_transition_whole_tree(index, rate_acc.size(), prob_inv);
  double acceptance_probability=min(1.0,nmethdensity*new_lhd.first/old_lhd.first * exp(-log1e100*(new_lhd.second-old_lhd.second)));

  bernoulli_distribution acceptance(acceptance_probability);
  bool decision=acceptance(re);

  if (decision) return;
  for (auto i =children.begin(); i!=children.end(); ++i) (*i)->set_meth_freqs_in_island_until_change(old_probabilities,index);

  delete transition[index];
  transition[index]=new island_meth_changes( prev_change_time,  prev_daughter,prev_methyl_freq);
  transition[index]->set_rate_acc(rate_acc);
  
}

  
void back_tree::delete_new_event(back_tree* tree_root,  DistObjects & dists, const double nmethdensity, const vector<double> & rate_acc, unsigned index, const double prob_inv,mt19937 & re)
{
  const double log1e100 = 230.2585;
  pair<double, int> old_lhd=tree_root->llhd_of_island_transition_whole_tree(index, rate_acc.size(), prob_inv);
  const vector<double> prev_change_time=transition[index]->get_change_time();
  //stores the change time we had previously so that we can revert if a jump is not accepted
  const vector<vector<double>> prev_methyl_freq=transition[index]->get_methyl_freq();
  //stores the frequencies we had previously so that we can revert if a jump is not accepted
  auto prev_daughter=transition[index]->get_daughter();
  auto old_probabilities=transition[index]->get_last_meth_freq();
  auto new_probabilities=dists.dirichlet111(re);
  uniform_int_distribution<int> event_del_dist(1, transition[index]->get_n_intervals()-1);
  const unsigned event_index=event_del_dist(re);
  transition[index]-> minus_meth_change(event_index);

  transition[index]-> update_all_transprobs();
  transition[index]-> update_partial_llhs();

  auto last_freq=transition[index]->get_last_meth_freq();
  for (auto i =children.begin(); i!=children.end(); ++i) (*i)->set_meth_freqs_in_island_until_change(last_freq,index);
  pair<double, int>new_lhd=tree_root->llhd_of_island_transition_whole_tree(index, rate_acc.size(), prob_inv); 
  bernoulli_distribution acceptance(min(1.0,nmethdensity*new_lhd.first/old_lhd.first * exp(-log1e100*(new_lhd.second-old_lhd.second))));
  bool decision=acceptance(re);
  if (decision) return;
  for (auto i =children.begin(); i!=children.end(); ++i)(*i)->set_meth_freqs_in_island_until_change(old_probabilities,index);
  delete transition[index];
  transition[index]=new island_meth_changes(prev_change_time, prev_daughter,prev_methyl_freq);
  transition[index]->set_rate_acc(rate_acc);
}


/*
 * This function finds the emprical frequency of U, P and M states in a node and stores them in a vector that can be used 
 * to initialize them in the whole tree as a starting point.
 *
 * No Input.
 *
 * Output: A vector output[i][j][k] where i constitutes the index of the cpg island, [j] is a place holder for frequency changes
 * and [k] denotes the relative frequency of state k in 0, 1, 2
 */ 

vector<vector<vector<double>>> back_tree::get_empirical_meth_freq()
{
  vector<vector<vector<double>>> output;
  for (unsigned i=0; i <cpg_container.size(); ++i)
    {
      vector<double> counts(3,0.1);
      double sum=0.3;
      for (unsigned j=0; j <cpg_container[i].states.size(); ++j)
	{
	  if (cpg_container[i].states[j]!=-1)
	    {
	      counts[cpg_container[i].states[j]]++;
	      sum++;
	    }
	}
      for (auto & j: counts) j/=sum;
      vector<vector<double>> longit={counts};
      output.push_back(longit);
    }
  return output;
}
  
/*
 * In the beginning of the MCMC inference we have methylation changes
 * along the tree initialized at some value. This function does this.
 */

void back_tree::initialize_island_meth_changes(vector<double> & rate_acc, vector<vector<vector<double>>> & meth_freq)
{
  
  for (unsigned i=0; i <cpg_container.size(); ++i)
    {
      vector<double> change_time={branch_length};     
      vector<vector<vector<double>>>  part_L;
      for (auto &j : cpg_container[i].states_as_llhood)	part_L.push_back(vector<vector<double>> (rate_acc.size(), j));
      auto daughter_pll=new island_partial_likelihoods_numscale( rate_acc, 0, part_L);
      auto meth_changes=new island_meth_changes(change_time , daughter_pll,meth_freq[i]); 
      transition.push_back(meth_changes);
    }
  
}



 /*
   * The following function prints the information about transition probability changes in internal
   * nodes to a stream.
   * The format here is:
   *
   * name is printed, whihc contains tells us the Newick representation of the subtree under the node.
   * In the following line we get the methylation island whose transition probablities we will
   * see in the form "methylation island " followed by its number.
   * After each methylation island line we get a line representing its changes,
   * of the form (probability, time interval for this probability). here the order is temporally going
   * from the parent node of the subnode in question downward.
   * This is done for all methylation islands.
   */
   
void back_tree:: print_transition_steps(std::ofstream & outputstream)
{
  outputstream<<name<< endl;
  for (unsigned i=0;i<transition.size(); i++)
    {
      if (parent==0)
	{
	  const vector<double> c_times=transition[i]->get_change_time();
	  outputstream<<"methylation island "<<i<< endl;
      	  const vector<vector<double>>  meth_freq=transition[i]->get_methyl_freq();
	  for (unsigned j=0;j<c_times.size(); j++)
	    {
	      outputstream<<" "<<c_times[j]<<"("<<meth_freq[j][0]<<"|"<<meth_freq[j][1]<<"|"<<meth_freq[j][2] <<")";
	    }
	  outputstream<<endl;
	}
      else if (transition[i]->get_n_intervals()>1)
	{
	  const vector<double> c_times=transition[i]->get_change_time();
	  outputstream<<"methylation island "<<i<< endl;
      	  const vector<vector<double>>  meth_freq=transition[i]->get_methyl_freq();
	  for (unsigned j=0;j<c_times.size(); j++)
	    {
	      outputstream<<" "<<c_times[j]<<"("<<meth_freq[j][0]<<"|"<<meth_freq[j][1]<<"|"<<meth_freq[j][2] <<")";
	    }
	  outputstream<<endl;
	}
    }
  for (auto node:children) node->print_transition_steps(outputstream);
}



/*
 * The following functiontakes a string of the form t(d_1|d_2|d_3) repeated and converts 
 * it into a pointer to a meth_freq class inserted at index i into the transitions. 
 * For this to work, transtion should allready be initialized with initialize meth_freq.
 */

void back_tree::set_transition_at_index_by_string(string  s, const vector<double> & rate_acc, const int i)
{
  boost::trim_if(s, boost::is_any_of("\t "));
  vector<string> event_strings;
  boost::split(event_strings, s, boost::is_any_of("  "));
  vector<double> t_change;
  vector<vector<double>> meth_freqs;
  for (auto & es: event_strings)
    {
      //cout << es << endl;
      if(es.size()>0)
	{
	  int index=0;
	  string time_string;
	  for ( index=0; es[index]!='('; ++index) time_string.push_back(es[index]);
	  t_change.push_back(stof(time_string));
	  index++;
	  vector<double> event_probs(3,0.0);
	  for (int j=0; j<3; ++j)
	    {
	      string prob_string;
	      while (es[index]!=')'&&es[index]!='|')
		{
		  prob_string.push_back(es[index]);
		  index++;
		}
	      index++;
	      event_probs[j]=stof(prob_string);
	    }
	  meth_freqs.push_back(event_probs);
	}
    }
  auto prev_daughter=transition[i]->get_daughter();
  branch_length=t_change.back();
  delete transition[i];
  transition[i]=new island_meth_changes(t_change, prev_daughter,meth_freqs);
  transition[i]->set_rate_acc(rate_acc);
  transition[i]-> update_all_transprobs();
  transition[i]-> update_partial_llhs();
}
      
	  
	      
/*
 * The following function samples a new first event and then accepts or rejects it according to MCMC.
 *
 * Parameters: 
 * DistObjects & dists: is a struct containing the different distributions in use throughout the code.
 * vector<double> & rate_acc : Are the rate acceleration factors used to model rate heterogeneity. 
 *                             This is a vector of length 3 that was sampled using the means of 
 *                             the quantiles of a Gamma distribution.
 * int index : index of the island we are potentially changing. This index is chosen uniformly over the
 *             number of all islands.
 * double prob_inv : rate of invariant sites
 * double pcrate: rate of IWEs
 * mt19937 & re: RNG in use throughout the simulation.
 */

void back_tree::sample_new_beginning(DistObjects & dists, const vector<double> & rate_acc, const int index, const double prob_inv, const double pcrate, mt19937 & re)
{

  //store the old probabilities. transition is a vector of pointers to island_meth_changes objects. These objects are defined in felsenstein_pruning.hh
  auto old_probabilities=(transition[index]->get_methyl_freq())[0];
  auto new_probabilities=dists.dirichlet111(re);
   
  //Calculate the old lhd:
  pair<double, int> old_lhd=llhd_of_island_transition_whole_tree(index, rate_acc.size(), prob_inv);
  pair<double, int> old_root_lhd=lhd_of_root_states(index, old_probabilities);
 
  old_lhd=llhdmulti(old_lhd, old_root_lhd);

  //new probabilities are set:
  transition[index]-> set_first_meth_freq(new_probabilities);

  //These probabilities have to be set in children of the node as well
  for (auto i =children.begin(); i!=children.end(); ++i) (*i)->set_meth_freqs_in_island_until_change(new_probabilities,index);

  //calculate lhd and determine acceptance:
  pair<double, int> new_lhd=llhd_of_island_transition_whole_tree(index, rate_acc.size(), prob_inv);
  pair<double, int> new_root_lhd=lhd_of_root_states(index, new_probabilities);
  new_lhd=llhdmulti(new_lhd, new_root_lhd);

  const double llhr=log10llhdratio(new_lhd,old_lhd);
  double acceptance_probability=min(1.0,pow(10,llhr));

  bernoulli_distribution acceptance(acceptance_probability);
  bool decision=acceptance(re);

 
  //If acceptance return as is, else revert tree to initial state.
  if (decision) return;
  for (auto i =children.begin(); i!=children.end(); ++i) (*i)->set_meth_freqs_in_island_until_change(old_probabilities,index);
  transition[index]-> set_first_meth_freq(old_probabilities);
}


/*
 * The folllowing function returns a vector ith length of n_islands and entries equal to the number
 * of events in the respective island
 * Vector returned is double so that it can be more easily used for probability calcs.
 *
 */

const vector<double> back_tree::get_event_number_in_islands()
{
  vector<double> output(transition.size(),0);
  for (unsigned i=0; i<transition.size(); ++i) output[i]=transition[i]->get_n_intervals()-1;
  return output;
}


/*
 * Parses the nodecontainer and returns the number of events per node as entries in a vector for probabiltiy caclualtion.
 */
const vector<double> get_event_numbers_in_nodes(const vector<back_tree*> & node_container)
{
  vector <double> output(node_container.size(),0);
  for (unsigned i=0; i<output.size(); ++i) output[i]=node_container[i]->get_n_meth_events();
  return output;
}
  


/*
 * Parses the nodecontainer and returns the branch lengths as entries in a vector for probabiltiy caclualtion.
 */
const vector<double> get_branch_lengths_nodes(const vector<back_tree*> & node_container)
{
  vector<double> output(node_container.size(),0);
  for (unsigned i=1; i<output.size();++i) output[i]=node_container[i]->branch_length;
  return output;
}
  
