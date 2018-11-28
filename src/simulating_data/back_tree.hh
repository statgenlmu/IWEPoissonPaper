
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

#ifndef BACKCLASS_H
#define BACKCLASS_H

#include <iostream>
#include <memory>
#include <list>
#include <map>
#include <deque>
#include <random>
#include <cassert>
#include "cpgislandcontainer.hh"
#include <fstream>

using namespace std;

class dyn_pop;
class mutation_type;
struct cpg_island;
class Event;
class Cpg_Island_Container;
class island_meth_changes;
struct DistObjects;

struct t_d
{
  t_d(){};
  t_d(vector <cpg_island> &cpg_islands);
  
  vector<vector <double>> all_rates; //cpg isle , then iteration
  vector <double> times;


  
  /*
   *This function stores info when a cpgi changes probability for new methylation states, 
   *which states these are and which cpgi this is and the time gap
   */
  void push_back(double time_gap, double rate, int w_cpgi)
  {
     times.push_back(time_gap);
     for (unsigned iter=0; iter <all_rates.size(); ++iter)
       all_rates[iter].push_back(all_rates[iter].back());
     all_rates[w_cpgi][all_rates[w_cpgi].size()-1]=rate;
  }
};




class back_tree
{
private:
  bool alive;
  string name;
  double gw_alpha;
  double gw_beta;
  vector<vector<vector<double>>> meth_freqs;
  //gives the length of the branch above it
  vector<cpg_island> cpg_container;
  vector<island_meth_changes*> transition;
  vector<vector<vector<vector<double>>>> transition_probs;
  /*
   * Contains the transition probabilities of all islands conditional on rate parameters.
   * transition_probs[a][b][c][d] is the probability for a site at island a to transtion
   * from methylation state c to the methylation state d if the rate acceleration factor is rate_acc[b]
   */
public:
  back_tree();
  ~back_tree();

  double branch_length;
  void initialize_island_meth_changes(vector<double> & rate_acc, DistObjects dist, double pcrate, const vector<unsigned> & island_sizes);
  void simulate_methylation(const vector<unsigned> & island_sizes, const vector<vector<int>> & site_factors);
  pair <double, int>  llhd_of_transition_with_fixed_ends(const int index_of_island, const vector<double> & rate_acc);
  pair <double, int>  llhd_of_transition_with_fixed_ends_all_islands(const vector<double> & rate_acc);
  pair <double, int>  llhd_of_transition_with_fixed_ends_whole_tree(const vector<double> & rate_acc);
  pair <double, int>  llhd_of_transition_with_fixed_ends_until_meth_change(const int index_of_island, const vector<double> & rate_acc);

  void stretch_branch(double const factor);
  void print_cpg_islands();

  void set_rate_acc_in_tree(vector<double> & rate_acc);
 
  const unsigned get_n_islands();
  const unsigned get_n_meth_events();

  const unsigned get_n_meth_events_at_island(unsigned i);
 
  
  void set_name(string n)
  {
    name=n;
  }
  void add_to_children(back_tree* a)
  {
    a->parent=this;    
    children.push_back(a);
  }
  string get_name()
  {
    return name;
  }
  void set_gw_alpha(double a)
  {
    gw_alpha=a;
  }
  double get_gw_alpha(void)
  {
    return gw_alpha;
  }
  void set_gw_beta(double a)
  {
    gw_beta=a;
  }
  double get_gw_beta(void)
  {
    return gw_beta;
  }
  
  double t_coal;
 
  int number_neutral;
  bool mark;

  back_tree * parent;
  vector<back_tree *> children;
  
  vector <cpg_island> cpg_data;//stores information about the cpgs contained here

  void set_meth_freqs_in_island_until_change(const vector<double> & meth_freq, unsigned index);
  void create_new_event(DistObjects & dists,vector<double> & rate_acc, double & pcrate, unsigned index, mt19937 & re);
  void delete_new_event(DistObjects & dists,vector<double> & rate_acc, double & pcrate, unsigned index, mt19937 & re);
    
  /*
   * This function parses files of the format CPG_ISLANDS_celltype.txt for the celltype name 
   * and initializes states for MCMC compuation.
   */

  void set_methylation_states_from_files();
  
  /*
   *This function prints all methylation states in a vector of
   *cpg islands into a file with the name of the node together with 
   *_methylation.csv. The format is: Each row contains a state, 
   *consisting of two numbers. The fist is the cpg island, the second is
   *the methylation state as a boolean
   */
  void print_cpg_islands(vector <cpg_island> & cpg_islands);
    
  
  void set_mutations(double neutral_rate, double t_last_coal, double gw_rate, vector <cpg_island> cpg_islands, mt19937 & re); //Goes through the tree and inserts a number of neutral mutations

  /*
   * The following function goes through the tree and sets neutral mutations methylations and 
   * changes in methylation probabilities along the tree. Different from version without step, in
   * so far as it does not use continous time, but time steps for improved speed.
   */ 
  void set_mutations_step(double neutral_rate, double time_step, double t_last_coal, double gw_rate, Cpg_Island_Container cpg_islands, mt19937 & re);

  
  
  /*
   * The following function sets methylations that happened between this node and its ancestor
   */
  void set_methylations(double t_last_coal, double gw_rate, vector <cpg_island> & cpg_islands, mt19937 & re);

  /*
   * The following function sets methylation that happened between this node and its ancestor. 
   * It does not do so continously, but after each time step you draw their quantity according 
   * to a poisson distribution.
   */
  void set_methylations_step(double time_step, double t_last_coal, double gw_rate, Cpg_Island_Container & cpg_islands, mt19937 & re);

  
  /*
   *The following function decides wheter the next even is a probability change or
   *methylation and then does that.
   */
  void methylation_or_p_event(Cpgi_Data & cpgi_data,vector<cpg_island> & cpg_islands, double time_gap, mt19937 & re);

  
  /*
   * The following function calculates the likelyhood of this node transitioning from its parent node.
   */

  double likelyhood_to_transition_from_parent_node(const int n_quantiles);

  /*
   * The following function iterates through a tree where methylation data is known at every node
   * and calculates its likelyhood given the the gw distribution parameters and the times of
   * the difference. For this purpose, not single methylation states are used but cpg islands are
   * classified according to islands and the sorted into corresponding quantiles.
   */

    double combined_likelyhood_of_tree_with_known_nodes_using_quantiles(const int n_quantiles);
  
  vector<vector<vector<double>>> tree_likelyhood();

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
   
  void print_transition_steps(std::ofstream & outputstream);
  
  /*
   * The following function iterates through thte tree and sets transition matrices from the transiton. 
   * Invoke after custom transitons are inserted.
   */
  void set_transition_matrices(vector<double> & rate_acc, const vector<unsigned> & island_sizes);

  /* 
 * Samples probability change events along ree then calculates 
 * the resulting transition probabilities and stores them in transition_probs 
 * Needs to call set_transtion_matrices to get the actual transition probabilities
 */
  void initialize_island_meth_changes_empty_transitions(vector<double> & rate_acc, const vector<unsigned> & island_sizes, const vector<vector<double>>& start_freqs);

};


vector <vector <double>> matrix_multi (vector <vector <double>> a, vector <vector <double>>b);

vector <vector < double>> transition_probability(vector<cpg_island> cpg_data, t_d transition_data, int isle);

/*
 * Given parameters of a beta distribution, the following function sets
 * the probabilities of methylation states in a vector of cpg_islands
 * equal to values drawn from that distribution.
 */
void initialize_cpg_island_probabilities(double gw_alpha, double gw_beta, vector<cpg_island> & cpg_islands, mt19937 & re);


/*
 * Given indeces of two quantiles, a time and a of quantile change rate,
 * the following function gives the probability that quantile change happens. 
 * This is done by giving the repective entry of the exponential of the rate matrix.
 */
double t_matrix_entry_JC(double m_rate, double time, int i, int j, int n_quantiles);


#endif





  

  





  
