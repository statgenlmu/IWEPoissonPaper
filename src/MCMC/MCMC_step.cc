
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

#include <iostream>
#include <fstream>
#include <random>
#include "MCMC_step.hh"
#include "back_tree.hh"
#include "functionobjects.hh"
#include <boost/math/distributions/beta.hpp>
#include <cmath>

using namespace std;
using boost::math::beta_distribution;


const double log10faculty(const int n)
{
  return n>0 ? log10(n)+log10faculty(n-1) : 0;
}


const double log10poisdensity(const double l, const int n)
{
  if (l==0&&n>0)
    {
      cerr<<"0 probability in log10poisdensity "<<endl;
      return -pow(10, 50)*n;
    }
  return n*log10(l)-l*log10(exp(1))-log10faculty(n);
}

/*
 * Returns the log 10 ratio a/b.
 * Here a and b are pairs representing probabilities, consiting of a leading term and a scale factor,
 * the actual value is of the probability represented by a is a.first * 10^(-100*a.second)
 */
const double log10llhdratio(const pair<double, int> & a,const pair<double, int> & b)
{
  return log10(a.first/b.first)-100*a.second+100*b.second;
}

/*
 * Returns the ratio c=a/b as a pair of double and int. The value represented by this pair is  c.first * 10^(-100*c.second)
 * Here a and b are pairs representing probabilities, consiting of a leading term and a scale factor,
 * the actual value is of the probability represented by a is a.first * 10^(-100*a.second)
 */
const pair<double,int> lhd_ratio(const pair<double, int> & a,const pair<double, int> & b)
{
  if (a.first==0||b.first==0) cerr<<"lhd_ratio: Likelihood is 0!."<<endl;
  return make_pair(a.first/b.first, a.second-b.second);
}


const pair <double, int> meanlhd(const vector<pair <double, int>> & llhds)
{
  int minscale=0;
  for (auto & j:llhds) if(j.first!=0) minscale=j.second;
  for (auto & j:llhds) if(j.first!=0) minscale=min(j.second, minscale);
  pair <double, int> meanlld={0,minscale};
  for (auto & j:llhds) meanlld.first+=(j.first/(pow(10,100*(j.second-minscale))*llhds.size()));
  return meanlld;
}

const double log10llhdratio(const vector<pair <double, int>> &  a,const vector<pair <double, int>> &  b)
{
  return log10llhdratio(meanlhd(a),meanlhd(b));
}



/*
 * The following function is used to sample new branch lengths. It picks a branch at random, 
 * suggests a new cranch legnth, and accepts the change according to MH probability. 
 * The new branch length is suggested according to the proposal distribution proposal_for_branch_length
 * contained in dists.

 * Parameters:
 * node_container: contains the pointers to the nodes in the tree.
 * dists: struct containing references to various probability distributions in use.
 *        Contains the proposal distribution.
 * pcrate: rate with which IWEs occur, measured in IWE per island per unit length.
 * prob_inv: Probability of sites to remain invariant.
 *
 * There is no output, howver if the change is accepted a tree node changes in length.
 */





void sample_new_branch_length(const vector<back_tree *> & node_container,const DistObjects & dists,vector<double> & rate_acc, const double pcrate, const double prob_inv, mt19937 & re)
{
  const double log10e100 = 230.2585;
  //choose branch
  uniform_int_distribution<int> which_branch(0,node_container.size()-2);
  int branch=which_branch(re);
  back_tree* node=node_container[branch+1];
  //sample new length
  double current_branch_length=node->branch_length;
  double new_branch_length=(*dists.proposal_for_branch_length)(current_branch_length, re);
  // calculate lklhd before legnth change
  pair<double,int> old_lhd=node_container[0]->llhd_of_transition_with_fixed_ends_whole_tree_all_islands(rate_acc.size(), prob_inv);
  //insert in tree
  node->branch_length=new_branch_length;
  node->stretch_branch(new_branch_length/current_branch_length);
  pair<double,int> new_lhd=node_container[0]->llhd_of_transition_with_fixed_ends_whole_tree_all_islands(rate_acc.size(),prob_inv);
 
  //determine acceptance
  double ln_prior_ratio=(*dists.length_prior).logratio(new_branch_length,current_branch_length);

  int n_islands=node->get_n_islands();
  int n_meth_events=node->get_n_meth_events();
  double ln_iwe_density_ratio=n_meth_events*log(new_branch_length/current_branch_length)
    -n_islands*pcrate*(new_branch_length-current_branch_length);

  double acceptance_prob=min(1.0,(new_lhd.first/old_lhd.first)*exp(ln_iwe_density_ratio+ln_prior_ratio-log10e100*(new_lhd.second-old_lhd.second)));
  bernoulli_distribution acceptance(acceptance_prob);
  //revert back if necessary
  if (acceptance(re)) return;
  node->stretch_branch(current_branch_length/new_branch_length);
  node->branch_length=current_branch_length;


}



/*
 * The following function is used to sample new iwe rates. It picks the current iwe rate and 
 * suggests a new one based on the distribution: pcrateproposal. It is neccessary for the 
 * proposal distribution to be symmetric.
 *
 * Parameters:
 * node_container: contains the pointers to the nodes in the tree.
 * dists: struct containing references to various probability distributions in use.
          It contains the proposal for new pc_rates.
 * pcrate: rate with which IWEs occur, measured in IWE per island per unit length.
 * prob_inv: Probaiblity of sites to remain invariant.
 *
 * There is no output, however if the change is accepted a tree node changes in length.
 */

 void sample_new_pcrate(const vector<back_tree *> & node_container,const DistObjects & dists,const vector<double> & rate_acc,  double & pcrate, mt19937 & re)
 {
   double new_pcrate=(*dists.pcrateproposal)(pcrate, re);
   double ln_meth_density_ratio=0;

   for (auto node: node_container)
     {
       if (node->parent!=0)
 	{
 	  int n_islands=node->get_n_islands();
 	  int n_meth_events=node->get_n_meth_events();
 	  double branch_length=node->branch_length;
 	  ln_meth_density_ratio+=n_meth_events*log(new_pcrate/pcrate)-n_islands*branch_length*(new_pcrate-pcrate);
 	}
     }
  
   double logpriorratio=((*dists.pcrate_prior).logratio(new_pcrate,pcrate));
   bernoulli_distribution acceptance(min(1.0,exp(ln_meth_density_ratio + logpriorratio))); 
   if(acceptance(re)) pcrate=new_pcrate;  
 }


/* 
 * The following function is used to insert or delete IWEs in the tree. It flips a coin to decide 
 * which of the two is going to happen and then picks a transition and island at random. 
 * If there are zero evets at the chosen location and deletion is chosen, nothing happens. Else
 * a specialized function for either deletion or creation of an IWE according to MH probabilities
 * is called.

 *Parameters:
 * node_container: contains the pointers to the nodes in the tree.
 * dists: struct containing references to various probability distributions in use.
 * pcrate: rate with which IWEs occur, measured in IWE per island per unit length.
 * prob_inv: Probaiblity of sites to remain invariant.
 * rate_acc: is a vector of rate acceleration factors used for likelihood calculation.
 */

void sample_new_pcevent(const vector<back_tree *> & node_container, DistObjects & dists,const vector<double> & rate_acc, double & pcrate, const double prob_inv, mt19937 & re)
{
  uniform_int_distribution<int> which_branch(0,node_container.size()-2);
  back_tree* node=node_container[which_branch(re)+1];
  back_tree* tree_root=node_container[0];
  //vector<double> n_events_in_island=node->get_event_number_in_islands();
  //discrete_distribution<int> which_island(n_events_in_island.begin(),n_events_in_island.end());
  uniform_int_distribution<int> which_island_uniform(0, node->get_n_islands()-1);
  unsigned index=which_island_uniform(re);
  bool zero_events=(node->get_n_meth_events_at_island(index))==0;
  bernoulli_distribution creation_or_deletion(0.5);
  if (creation_or_deletion(re))
    {
      double nmethdensity=(node->branch_length)*pcrate/(node->get_n_meth_events_at_island(index)+1);
      node->create_new_event(tree_root, dists, nmethdensity, rate_acc, index, prob_inv,re);
      return;
    }
  if (!zero_events)
    {
      double nmethdensity=(node->get_n_meth_events_at_island(index))/((node->branch_length)*pcrate);
      node->delete_new_event(tree_root, dists, nmethdensity, rate_acc, index,prob_inv,  re);
    }
}

/* 
 * The following function is used to insert or delete IWEs in the tree. It flips a coin to decide 
 * which of the two is going to happen and then picks a transition and island at random. 
 * If there are zero evets at the chosen location and deletion is chosen, nothing happens. Else
 * a specialized function for either deletion or creation of an IWE according to MH probabilities
 * is called.

 *Parameters:
 * node_container: contains the pointers to the nodes in the tree.
 * dists: struct containing references to various probability distributions in use.
 * pcrate: rate with which IWEs occur, measured in IWE per island per unit length.
 * prob_inv: Probability of sites to remain invariant.
 * rate_acc: is a vector of rate acceleration factors used for likelihood calculation.
 */

void sample_new_pcevent_tree(const vector<back_tree *> & node_container, DistObjects & dists,const vector<double> & rate_acc, const double pcrate, const double prob_inv, mt19937 & re)
{
  back_tree* tree_root=node_container[0];
  const double tree_length=tree_root->get_subtree_length();
  const double n_meth_events=tree_root->get_subtree_nevents();
  const bool zero_events=n_meth_events==0;
  bernoulli_distribution creation_or_deletion(0.5);
 
  if (creation_or_deletion(re))
    {
      const vector<double> branch_lengths=get_branch_lengths_nodes(node_container);
      discrete_distribution<int> which_branch_by_length(branch_lengths.begin(), branch_lengths.end());
      back_tree* node=node_container[which_branch_by_length(re)];
      const int n_islands= node->get_n_islands();
      uniform_int_distribution<int> which_island_uniform(0,n_islands-1);
      const unsigned index=which_island_uniform(re);
      const double nmethdensity=n_islands*tree_length*pcrate/(n_meth_events+1);
      node->create_new_event(tree_root, dists, nmethdensity, rate_acc, index, prob_inv,re);
      return;
    }
  if (!zero_events)
    {
       
      const vector<double> branch_weights=get_event_numbers_in_nodes(node_container);

      discrete_distribution<int> which_branch_event_based(branch_weights.begin(), branch_weights.end());
      const int branch_index=which_branch_event_based(re);
      back_tree* node=node_container[branch_index];
      const int n_islands= node->get_n_islands();
      vector<double> n_events_in_island=node->get_event_number_in_islands();
      discrete_distribution<int> which_island_event_based(n_events_in_island.begin(),n_events_in_island.end());
    
      const unsigned index=which_island_event_based(re);
      const double nmethdensity=n_meth_events/(tree_length*pcrate*n_islands);
      node->delete_new_event(tree_root, dists, nmethdensity, rate_acc, index,prob_inv,  re);
    }
}


/* 
 * The following function is used to change the shape parameter of the Gamma distribution underlying 
 * rate heterogenuity. It samples a new porposed shape parameter, calculates respective rates from

 *Parameters:
 * node_container: contains the pointers to the nodes in the tree.
 * dists: struct containing references to various probability distributions in use.
          It contains the proposal for new shape parameters and a method to calculate 
          new rate heterogenuity.
 * pcrate: rate with which IWEs occur, measured in IWE per island per unit length.
 * prob_inv: Probaiblity of sites to remain invariant.
 * rate_acc: is a vector of rate acceleration factors used for likelihood calculation.
 */

void sample_new_rate_accell(const vector<back_tree *> & node_container,DistObjects & dists,vector<double> & rate_acc,const double prob_inv,mt19937 & re) //TODO: Dont give rate_acc to function, everything is in gamma_rate_cosntruction
{
  const double log1e100 = 230.2585; // ln of 10^100
  back_tree* node=node_container[0];
  const double old_alpha=dists.gamma_rate_construction.get_alpha();
  const double new_alpha=(*dists.proposal_for_gamma_shape)(old_alpha, re);
  dists.gamma_rate_construction.set_alpha(new_alpha);
  auto new_rates=dists.gamma_rate_construction.get_rate_acc();
  // calculate lklhd before rate change
  auto old_lhd=node->llhd_of_transition_with_fixed_ends_whole_tree_all_islands(rate_acc.size(), prob_inv); //TODO: lhd does not need rate_acc
  //insert in tree
  node->set_rate_acc_in_tree(new_rates);
  auto new_lhd=node->llhd_of_transition_with_fixed_ends_whole_tree_all_islands(rate_acc.size(),prob_inv);
  //determine acceptance
  const double ln_prior_ratio=(*dists.gamma_shape_prior).logratio(new_alpha,old_alpha);
  const double acceptance_probability=min(1.0,(new_lhd.first/old_lhd.first)*exp(ln_prior_ratio-log1e100*(new_lhd.second-old_lhd.second)));
  bernoulli_distribution acceptance(acceptance_probability);
  bool decision=acceptance(re);
  //revert back if necessary
  if (decision)
    {
      rate_acc=new_rates;
      return;
    }
  dists.gamma_rate_construction.set_alpha(old_alpha);
  node->set_rate_acc_in_tree(rate_acc);
}


/*
 * This takes a reference to the first node, picks an island at random and proposes changing the first entry of the island. It accepts or rejects this proposal acocrding to MCMC.
 *
 */
void sample_new_beginning(back_tree* first_node,  DistObjects & dists, const vector<double> & rate_acc, const double prob_inv, const double pcrate, mt19937 & re)
{
  uniform_int_distribution<int> which_island(0,first_node->get_n_islands()-1);
  unsigned index=which_island(re);
  first_node->sample_new_beginning( dists,rate_acc,index,prob_inv, pcrate, re);
}


/*
 * The following function makes a MCMC jump where invariant probability is changed, if accepted
 */
void sample_new_invariant_prop(back_tree* tree_root, DistObjects & dists,const vector<double> & rate_acc, double & prob_inv, mt19937 & re)
{
  const double log1e100 = 230.2585; // ln(10^100)
  auto old_lhd=tree_root->llhd_of_transition_with_fixed_ends_whole_tree_all_islands(rate_acc.size(),prob_inv);
  const double new_prob_inv=(*dists.invariant_probability_poposal)(prob_inv,re);
  auto new_lhd=tree_root->llhd_of_transition_with_fixed_ends_whole_tree_all_islands(rate_acc.size(),new_prob_inv);
  const double logpriorratio=(*dists.invariant_probability_prior).logratio(new_prob_inv,prob_inv); 
  const double acceptance_probability=min(1.0,new_lhd.first/old_lhd.first * exp(-(new_lhd.second-old_lhd.second) * log1e100 + logpriorratio));
 
  bernoulli_distribution acceptance(acceptance_probability);
  const bool decision=acceptance(re);
  //accept if necessary
  if (decision) prob_inv=new_prob_inv;
}
  

/*
 * This function decides which MCMC step is going to be executed and then calls the function to execute this MCMC step.
 * It has the additional option which_steps. This is a boolean vector with length 6. If an entry is 0 the corresponding MCMC step is not executed.

 * Indeces correspond to:
 * 0: sampling new branch legnth
 * 1: sampling new rate acceleration 
 * 2: sampling new IWE rate
 * 3: Sample a new IWE or delete an existing one
 * 4: sampling new invariant rate
 * 5 Sample new first IWE

 */

int MCMC_step(vector<back_tree *> & node_container,DistObjects & dists,vector<double> & rate_acc,double & pcrate, double & prob_inv, const vector<bool> & which_steps,mt19937 & re)
  {

    vector<double> relative_weights={0.1,0.01,0.1,1, 0.01, 0.3};
    discrete_distribution<int> what_jump(relative_weights.begin(), relative_weights.end());
    int choice=what_jump(re);
    if (!which_steps[choice]) return choice;
    //do the respective step
  switch(choice)
    {
    case 0:
      ////cout<<"sample branch length "<<endl;
      sample_new_branch_length(node_container,dists,rate_acc,pcrate,prob_inv, re);
      break;
    case 1:
      ////cout<<"sample rate accell "<<endl;
      sample_new_rate_accell(node_container,dists,rate_acc,prob_inv,re);
      break;
    case 2:
      ////cout<<"sample pc rate "<<endl;
      sample_new_pcrate( node_container,dists,rate_acc, pcrate, re);
      break;
    case 3:
      ////cout<<"sample pc event "<<endl;
      sample_new_pcevent_tree( node_container,dists,rate_acc, pcrate,prob_inv, re);
      break;
    case 4:
      ////cout<<"Sample new invariant probability "<<endl;
      sample_new_invariant_prop(node_container[0],  dists, rate_acc, prob_inv, re);
      break;
    case 5:
      sample_new_beginning(node_container[0], dists, rate_acc, prob_inv, pcrate,  re);
      break;
    default: break;
    }
  return choice; 
  
  }

