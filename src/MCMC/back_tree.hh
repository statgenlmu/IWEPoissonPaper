
#ifndef BACKCLASS_H
#define BACKCLASS_H

#include <iostream>
#include <memory>
#include <list>
#include <map>
#include <deque>
#include <random>
#include <cassert>


using namespace std;

class dyn_pop;
class mutation_type;
struct cpg_island;
class Event;
class Cpg_Island_Container;
class island_meth_changes;
struct DistObjects;



/*
 * THe back_tree class stores the tree in a graph structure with parent and child nodes. 
 * It features various ways to calculate different likelihoods that are needed in the program 
 * and stores information about the transitions on the branches as well as the cpg island data 
 * in the nodes.
 */


class back_tree
{
private:

  friend class likelihoodstests;//for testing with the cxx testsuite
  string name;
  
  //gives the length of the branch above it
   
  vector<island_meth_changes*> transition;
  // this vector contains the addresses of the island meth changes associated
  //with specific islands. For more detail about these functions consult the felsenstain_pruning.hh
  vector <cpg_island> cpg_container;//stores information about the cpgs contained here.
  //Specifically contains the states the cpg islands are in.  struct cpg islands is in mutationclass.hh
  
public:
  back_tree();
  ~back_tree();

  double branch_length;
  back_tree * parent;
  //contains the address of the parent node
  vector<back_tree *> children;
  //contains the addresses of all the daughter nodes of a node

  /*
*************Likelihood functions************************
*/

  
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
 

/*
 * The following function checks wheter a site changes along the tree or stays invariant.
 * it returns true if there are no detectable changes, false otherwise.
 *
 * The input argument state is the state the site was in in cells that have been checked so far. It is -1 at initialization or if states
 * so far have been NAs. Else 0,1,2 corresponding to U,P,M.
 * i stands for the index of the island, s for the index of the site.
 */


  bool possible_inv(const int i, const int s, int & state);
  
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
  pair<double,int> lhd_of_root_states(const int index, const vector<double> & meth_freq);


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
 

  pair<double, int> llhd_of_single_site_transition_with_fixed_ends(const int i, const int r, const int s);
  
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

  const pair<double, int> llhd_of_single_site_transition_whole_tree(const int i, const int r, const int s);



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

  const pair<double, int> llhd_of_island_transition_whole_tree(const int i, const int n_rates, const double prob_inv);


  
/*
 * Function for determining the combined likelihoods of transition of all sites averaged over 
 * rate factors of all islands along the whole tree. Should be called at the 
 * tree root.
 *
 * Input:
 * n_rates: # rate acceleration factors
 * prob_inv is the probability that the site is invariant
 *
 * output: pair of probability and scale factor
 */

  pair <double, int>  llhd_of_transition_with_fixed_ends_whole_tree_all_islands(const int n_rates, const double prob_inv);


   /*
*************End Likelihood functions************************
*/



  
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
   * This is done for all methylation islands with changes.
   */
   
  void print_transition_steps(std::ofstream & outputstream);


  
  /*
   * Fucntion that stretches all time differences between this node and its parent stored in
   * the island_meth classes.
   */
  void stretch_branch(double const factor);
  void set_rate_acc_in_tree(vector<double> & rate_acc);
  const unsigned get_n_islands();
  const unsigned get_n_meth_events();

  /*
   * The following function returns the length of the subtree bellow the current node.
   */
  const double get_subtree_length()
  {
    double T=0;
    for (auto node : children)T+=node->branch_length+node->get_subtree_length();
    return T;
  }

  
  /*
   * The following function returns the length of the subtree bellow the current node.
   */
  const int get_subtree_nevents()
  {
    int n=0;
    for (auto node : children) n+=node->get_n_meth_events()+node->get_subtree_nevents();
    return n;
  }

  /*
   * Islands can have multiple events on them along a branch - it is often intersting 
   * to know wheter they merely cluster on one branch or are spread out, so I give out 
   * the the number of islands carrying an event as well.
   */
 
  const unsigned get_n_affected_islands();
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
  const string get_name()
  {
    return name;
  }

  /* If a new methylation event happens, methylation probabilities
   * at the beginnings of daughterr nodes 
   * must be adjusted iteratively, until those nodes have a methylation change. The following
   * function does that.
   *
   * Input:
   * meth_freq is the methylation probability that is inserted
   * index is the index of the island where the methylation event is inserted
   */
  void set_meth_freqs_in_island_until_change(const vector<double> & meth_freq, const unsigned index);

  /*
   * Te following two functions respectively insert a methylation change event or delete one from an island and then use the acceptance probability derived from Metropolis Hastings to accept or reject this lcaim. 
   */
 
  void create_new_event(back_tree* tree_root,  DistObjects & dists,const double nmethdensity, const vector<double> & rate_acc, unsigned index, const double prob_inv, mt19937 & re);
  void delete_new_event(back_tree* tree_root,  DistObjects & dists,const double nmethdensity, const vector<double> & rate_acc, unsigned index, const double prob_inv, mt19937 & re);
    
  /*
   * This function parses files of the format CPG_ISLANDS_celltype.txt for the celltype name 
   * and initializes states for MCMC compuation.
   */

  void set_methylation_states_from_files();
  void set_methylation_states_from_files(const string & file_path);

  
/*
 * This function finds the emprical frequency of U, P and M states in a node and stores them in a vector that can be used 
 * to initialize them in the whole tree as a starting point.
 *
 * No Input.
 *
 * Output: A vector output[i][j][k] where i constitutes the index of the cpg island, [j] is a place holder for frequency changes
 * and [k] denotes the relative frequency of state k in 0, 1, 2
 */ 

  vector<vector<vector<double>>> get_empirical_meth_freq();
  void initialize_island_meth_changes(vector<double> & rate_acc, vector<vector<vector<double>>> & meth_freq);

 
  const unsigned get_n_meth_events_whole_tree();
  const double get_complete_tree_length();

  /*
   * The following functiontakes a string of the form t(d_1|d_2|d_3) repeated and converts 
   * it into a pointer to a meth_freq class inserted at index i into the transitions. 
   * For this to work, transtion should allready be initialized with initialize meth_freq.
   */

  void set_transition_at_index_by_string(string  s, const vector<double> & rate_acc, const int i);

  
/*
 * The folllowing function returns a vector ith length of n_islands and entries equal to the number
 * of events in the respective island*
 *
 */

  const vector<double> get_event_number_in_islands();
  	      
/*
 * The following function samples a new first event and then accepts or rejects it according to MCMC.
 */

  void sample_new_beginning(DistObjects & dists, const vector<double> & rate_acc, const int index, const double prob_inv,  const double pcrate, mt19937 & re);

};



/*
 * The following function takes a pair of probability and scale factor a 
 * and probability and scalefactor B and returns their product with appropriate scale factor. 
 * Scalefactor repreresents the power of -10 that needs to be multiplied to the probability to 
 * get the actual value.
 */

const pair <double, int> llhdmulti(const pair <double, int>  A, const pair <double, int>  B);



/*
 * Parses the nodecontainer and returns the number of events in each node as entries in a vector for probabiltiy caclualtion.
 */
const vector<double> get_event_numbers_in_nodes(const vector<back_tree*> & node_container);


/*
 * Parses the nodecontainer and returns the branch lengths as entries in a vector for probabiltiy caclualtion.
 */
const vector<double> get_branch_lengths_nodes(const vector<back_tree*> & node_container);
  


#endif





  

  





  
