   /* This file is part of IWEinferrence.

    IWEinferrence is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IWEinferrence is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with IWEinferrence.  If not, see <https://www.gnu.org/licenses/>. */


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
#include "cpgislandcontainer.hh"
#include "../MCMC/felsenstein_pruning.hh"
#include "functionobjects.hh"
#include <regex>
#include <boost/algorithm/string.hpp>


using namespace std;
using boost::math::beta_distribution;

class Event;
extern mt19937 re;

back_tree::back_tree()
{
  parent=nullptr;
  t_coal=0;
  number_neutral=0;
  alive=true;
  mark=false;
}

back_tree::~back_tree()
{
 
  parent=0;
}




t_d::t_d(vector <cpg_island> &cpg_islands)
{
  for (unsigned i=0; i<cpg_islands.size();++i)
    {
      vector <double> a(1,cpg_islands[i].r);
      all_rates.push_back(a);
    }
}

Cpgi_Data::Cpgi_Data(){}
Cpgi_Data::Cpgi_Data(vector <cpg_island> & cpg_islands)
  {
    methylation_rate=0;
    total_rate=0;
    for (unsigned i=0; i< cpg_islands.size(); ++i)
      {
	methylation_rate+=cpg_islands[i].m_rate;
	m_rates.push_back(cpg_islands[i].m_rate);
	c_rates.push_back(cpg_islands[i].change_rate);
	total_rate+=cpg_islands[i].m_rate+cpg_islands[i].change_rate;
       }
  }
 


pair <double, int> llhdmulti(const pair <double, int>  A, const pair <double, int>  B)
{
  return {A.first*B.first,A.second+B.second};
}

/* 
 * Samples probability change events along ree then calculates 
 * the resulting transition probabilities and stores them in transition_probs 
 */
void back_tree::initialize_island_meth_changes(vector<double> & rate_acc, DistObjects dist, double pcrate, const vector<unsigned> & island_sizes)
{
  transition.clear();
  
  for (unsigned i=0; i <island_sizes.size(); ++i)
    {
      vector<double> change_time(1, branch_length);
      vector<vector<double>> meth_freq;
      if (parent!=0)
	{
	  vector<vector<double>> parent_meth_freq=parent->meth_freqs[i];
	  meth_freq={parent_meth_freq.back()};
	}
      else meth_freq={dist.dirichlet111(re)};
      uniform_real_distribution<double> location_of_changes(0, branch_length);
      poisson_distribution<int> number_of_events(branch_length*pcrate);
      int n_events=parent==0?0:number_of_events(re);
      for (int j=0; j<n_events; ++j)
	{
	  change_time.push_back(location_of_changes(re));
	  meth_freq.push_back(dist.dirichlet111(re));
	}
      sort(change_time.begin(), change_time.end());
	  
      vector<vector<vector<double>>>  part_L;
      for (unsigned j=0; j <3; ++j)
	{
	  part_L.push_back(vector<vector<double>> (rate_acc.size(), {j==0,j==1,j==2}));
	}
	  

      auto daughter_pll=new island_partial_likelihoods_numscale( rate_acc, 0, part_L);
      auto meth_changes=new island_meth_changes(change_time , daughter_pll,meth_freq);
      
      transition.push_back(meth_changes);

      vector<vector<vector<double>>> transition_matrix; //first vecotr for different rate_accels
      for (unsigned j=0; j<rate_acc.size(); ++j)
	{
	  vector<vector<double>> tm(3,vector<double>(3,0));
	    
	  for (int k=0;k<3;++k)
	    for(int l=0;l<3;++l)
	      {
		tm[k][l]=(meth_changes->get_part_L(l,j))[k];
		//cout<<"tm["<<k<<"]["<<l<<"] is " <<tm[k][l]<<endl;
	      }
	  transition_matrix.push_back(tm);
	  
	}
      transition_probs.push_back(transition_matrix);
      meth_freqs.push_back(meth_freq);
				 
    }
  for (auto node: children)
    node->initialize_island_meth_changes(rate_acc, dist, pcrate, island_sizes);
}

/* 
 * Samples probability change events along tree then calculates 
 * the resulting transition probabilities and stores them in transition_probs 
 * Needs to call set_transtion_matrices to get the actual transition probabilities
 */
void back_tree::initialize_island_meth_changes_empty_transitions(vector<double> & rate_acc, const vector<unsigned> & island_sizes, const vector<vector<double>>& start_freqs)
{
  
  for (unsigned i=0; i <island_sizes.size(); ++i)
    {
      vector<double> change_time(1, branch_length);
      vector<vector<double>> meth_freq;
      if (parent!=0)
	{
	  vector<vector<double>> parent_meth_freq=parent->meth_freqs[i];
	  meth_freq={parent_meth_freq.back()};
	}
      else meth_freq={start_freqs[i]};
    
	  
      vector<vector<vector<double>>>  part_L;
      for (unsigned j=0; j <3; ++j)
	{
	  part_L.push_back(vector<vector<double>> (rate_acc.size(), {j==0,j==1,j==2}));
	}
	  

      auto daughter_pll=new island_partial_likelihoods_numscale( rate_acc, 0, part_L);
      auto meth_changes=new island_meth_changes(change_time , daughter_pll,meth_freq);
      
      transition.push_back(meth_changes);

     
      meth_freqs.push_back(meth_freq);
				 
    }
  for (auto node: children)
    node->initialize_island_meth_changes_empty_transitions(rate_acc, island_sizes, start_freqs);
}



/*
 * The following function iterates through thte tree and sets transition matrices from the transiton. 
 * Invoke after custom transitons are inserted.
 */
void back_tree::set_transition_matrices(vector<double> & rate_acc, const vector<unsigned> & island_sizes)
{
  int index=0;
  for (const auto & i:island_sizes)
    {
      vector<vector<vector<double>>> transition_matrix; //first vecotr for different rate_accels
      for (unsigned j=0; j<rate_acc.size(); ++j)
	{
	  vector<vector<double>> tm(3,vector<double>(3,0));
	    
	  for (int k=0;k<3;++k)
	    for(int l=0;l<3;++l)
	      {
		tm[k][l]=(transition[i]->get_part_L(l,j))[k];
		//cout<<"tm["<<k<<"]["<<l<<"] is " <<tm[k][l]<<endl;
	      }
	  transition_matrix.push_back(tm);
	  
	}
      transition_probs.push_back(transition_matrix);
      index++;
    }
  for (auto node: children)
    node->set_transition_matrices(rate_acc, island_sizes);
}




/*
 * The following uses the tree with allready sampled probability change events to set the cpg_containers to contain methylation states sampled with probabilities erived from the sampled probability change events stored in transition_probs*/
void back_tree::simulate_methylation(const vector<unsigned> & island_sizes, const vector<vector<int>> & site_factors)
 {
   if (parent==0)
     {
       for (unsigned i=0; i<island_sizes.size(); ++i)
	 {
	   auto meth_start=transition[i]->get_methyl_freq(0);
	   discrete_distribution<int> initial_state_dist(meth_start.begin(), meth_start.end());
	   cpg_island cpgi;
	   for (unsigned j=0; j<island_sizes[i]; ++j) cpgi.states.push_back(initial_state_dist(re));	     
	   cpg_container.push_back(cpgi);
	 }
     }
   else
     {
       for (unsigned i=0; i<island_sizes.size(); ++i)
	 {
	   cpg_island cpgi;
	   for (unsigned j=0; j<island_sizes[i]; ++j)
	     {
	       int parent_state=parent->cpg_container[i].states[j];
	       if (site_factors[i][j]==-1) cpgi.states.push_back(parent_state);
	       else
		 {
		   auto meth_prob=transition_probs[i][site_factors[i][j]][parent_state];
		   discrete_distribution<int> state_dist(meth_prob.begin(), meth_prob.end());
		   cpgi.states.push_back(state_dist(re));
		 }
	     }
	   cpg_container.push_back(cpgi);
	 }
     }
   for (auto node: children) node->simulate_methylation(island_sizes, site_factors);
 }

const unsigned back_tree::get_n_meth_events()
{
  unsigned n=0;
  for (auto & i : transition)
    n+=i->get_n_intervals()-1;
  return n;
}


void back_tree::set_meth_freqs_in_island_until_change(const vector<double> & meth_freq, unsigned index)
{
  transition[index]->set_first_meth_freq(meth_freq);
  if (transition[index]->get_n_intervals()==1)
    for (auto i =children.begin(); i!=children.end(); ++i)
      (*i)->set_meth_freqs_in_island_until_change(meth_freq,index);
}

const unsigned back_tree::get_n_islands() {return transition.size();}
const unsigned back_tree::get_n_meth_events_at_island(unsigned i) {
  return transition[i]->get_n_intervals()-1;
}
 
   
/*
 * This function prints all methylation states in a vector of
 * cpg islands into a file with the name of the node indicated 
 * in the form of CPG_ISLANDS_[name]_txt. The format is: Each row 
 * contains the methylation states of an island.
 *
 */
void back_tree::print_cpg_islands()
{
 
  string file_name="CPG_ISLANDS_"+name+".txt";
  ofstream island_stream;
  island_stream.open(file_name);
  int index=0;
  //cout<<"cpgislands # "<<cpg_container.size()<<endl;
  for (auto i : cpg_container)
    {
      index++;
      for (auto j : i.states)
	{
	  island_stream<<j;
	}
      island_stream<<endl;
    }
  island_stream.close();

  for (auto node:children)
    node->print_cpg_islands();
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
     
	  vector<double> c_times=transition[i]->get_change_time();
	  outputstream<<"methylation island "<<i<< endl;
	  
	  vector<vector<double>>  meth_freq=transition[i]->get_methyl_freq();
	  for (unsigned j=0;j<c_times.size(); j++)
	    {
	      outputstream<<" "<<c_times[j]<<"("<<meth_freq[j][0]<<"|"<<meth_freq[j][1]<<"|"<<meth_freq[j][2] <<") ";
	    }
	  outputstream<<endl;
    }
  for (auto node:children)
    node->print_transition_steps(outputstream);
}




bool is_transition_string(string &s)
{
  regex r("(\\ ?\\d+(\\.\\d+)?{\\d+(\\.\\d+)?\\|\\d+(\\.\\d+)?\\|\\d+(\\.\\d+)?})+");
  return regex_match(s,r);
}

bool is_island_index (string &s)
{
  regex r("methylation\\ island\\ \\d+");
  return regex_match(s,r);
}

int get_island_index_from_string(string &s)
{
  vector<string> sub_strings;
  boost::split(sub_strings, s, boost::is_any_of(" "));
  return stoi(sub_strings.back());
}

	      
