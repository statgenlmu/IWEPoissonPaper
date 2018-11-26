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
#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <cassert>
#include <regex>
#include "mutationclass.hh"
#include "back_tree.hh"
#include <fstream>

using namespace std;

class back_tree;

void set_cpg_island_sizes_default(vector<cpg_island> & cpg_islands)
{
  for (unsigned i=0; i<45000; i++)
    {
      cpg_island a;
      a.states=vector<int>(20,0);
      cpg_islands.push_back(a);
    }
}

/*
 * The following function sets the starting composition of the cpg_islands.
 * It searches argv for the number of  cpg islands marked by the string -n_islands.
 * If those are not found then
 * default is applied. If they are found, it searches argv for the the string -island_sizes
 * and then parses the following numbers for the sizes of the respective cpg isladns which
 *  are then initialized.
 */

void set_cpg_island_sizes(vector<cpg_island> & cpg_islands, int argc , char ** argv)
{

  bool found_n_islands=false;
  unsigned n_islands=0;
  for (int i=0; i<argc-1; i++)
    {
      string comp="-n_islands";
      if (argv[i]==comp)
	{
	  n_islands=atoi(argv[i+1]);
	  found_n_islands=true;
	}
    }
  if (!found_n_islands)
    {
      cerr<<"Warning: CPG island sizes not set. Initialized 45000 CPG islands with a size of 20 sites as default"<<endl;
      set_cpg_island_sizes_default(cpg_islands);
      return;
    }

  
  for (int i=0; i<argc-1; i++)
    {
      string comp="-island_sizes";
      if (argv[i]==comp)
	{
	  for (unsigned j=i+1; j <= n_islands; ++j)
	    {
	      cpg_island a;
	      a.states=vector<int>(atoi(argv[j]));
	      cpg_islands.push_back(a);
	    }
	  break;
	}
    }
}



/*
 * Checks wheter a given string is a number and returns true,
 * if yes, false else
 */

bool string_is_number(string a)
{
  regex e ("^-?\\d*\\.?\\d+");
  if (regex_match (a,e))
    {
      return true;
    }
  return false;
}

/*
 * This function parses files of the format CPG_ISLANDS_celltype.txt for each celltype name 
 * in the tree and initializes states for MCMC compuation.
 */

void set_methylation_states_from_files(vector<back_tree*> & node_container)
{
  for (auto i:node_container) i->set_methylation_states_from_files();
}


 

/***********************Building a tree from its newick representation************************/

/*
 * Finds the name of a root of a sub tree by parsing a Newick string s
 * For example if given the string (a,b)Alex:5.0 it returns the string 'Alex'
 */

string find_newick_name(string s)
{
  unsigned n;
  string r;
  for ( n=s.size(); n>0&&s[n-1]!=')';--n);
  for (; n<s.size()&&s[n]!=':';++n) if (s[n]!=';') r.push_back(s[n]);
  return r;
}



/*
 * Finds the distance of an edge adjacent to the root by parsing a Newick string s
 */

double find_newick_distance(string s)
{
  unsigned n;
  string r;
  if (s.find(":")==string::npos) return 0;
  for ( n=s.size(); n>0&&s[n-1]!=':';--n);
  for (; n<s.size();++n) r.push_back(s[n]);
  if (r=="") return 0;
  if (! string_is_number(r))
    cerr<<"find_newick_distance: Branch length is NAN"<<endl;
  return stof(r);
}



/* 
 * The following function parses a Newick string for subtree strings.
 * For example for the string ((a,B):7.0,sdh,(sa,sb))Alex:5.0 it would return
 * {"(a,B):7.0","sdh","(sa,sb)"}
 */

vector<string> find_sub_tree_strings(string s)
{
  unsigned int nesting=1;
  vector<string> output;
  string substring;
  for (int i=1; nesting!=0; ++i)
    {
      if (nesting==1&&((s[i]==',')||(s[i]==')')))
	{
	  output.push_back(substring);
	  substring.clear();
	}
      else substring.push_back(s[i]);
      nesting=nesting + (s[i]=='(') - (s[i]==')');
    }
  return output;
}
    
/*
 * The following function parses a string s in Newick format and builds a binary tree from it and then 
 * returns the root.
 * The node_container then contains all the nodes, with the node_container[0] being the root when called
 * with the complete string.
 */
back_tree* tree_from_string(string s, vector<back_tree*> & node_container)
{
  //adjust s if necessary
  if (s.size()==0)
    return 0;
  if (s.back()==';')
    s.pop_back();

  //create return node
  back_tree * sub_tree_root(new back_tree); //just sub_tree, because the function is recursive
 
  sub_tree_root->set_name(find_newick_name(s));
  sub_tree_root->branch_length=find_newick_distance(s);
  node_container.push_back(sub_tree_root);

  if (s.find(")") != string::npos)
    {
      vector<string> sub_tree_strings=find_sub_tree_strings(s);
      for (auto i:sub_tree_strings)
	sub_tree_root->add_to_children(tree_from_string(i, node_container));
    }
  
  return sub_tree_root;
}




/*
 * The following function looks for a tree in Newick format in the file tree.sol. 
 * If it finds this tree it fills node_container with references to the nodes so that
 * node_container[0] is the root of the tree.
 */

void initialize_tree_from_file(vector<back_tree*> & node_container)
{
  ifstream tree_stream("tree.sol");
  string tree_string;
  if (!getline(tree_stream, tree_string))
    {
      cerr<<"No tree found, \"tree.sol\" empty or not existent. Default configuration was used."<<endl;
      return;
    }
 tree_from_string(tree_string,node_container);
}


/******************Reading in data about methylation events to resume computation after stop************/

bool is_methylation_island_index(string& s)
{
  regex e("methylation island +\\d+");
  return regex_match(s,e);
}

bool is_transition_frequency_data(string & s)
{
  regex e("(\\s\\d+\\.\\d+\\(\\d+\\.\\d+\\|\\d+\\.\\d+\\|\\d+\\.\\d+\\)\s+)+");
  // e("(\\d+\\.\\d+\\(\\d+\\.\\d+\\|\\d+\\.\\d+\\|\\d+\\.\\d+\\)  )*(\\d+\\.\\d+\\(\\d+\\.\\d+\\|\\d+\\.\\d+\\|\\d+\\.\\d+\\))");
  return regex_match(s,e);
}

const int find_island_index(string & s)
{
  int i=s.size()-1;
  for ( i=s.size()-1; s[i]>='0' && s[i]<='9'; --i);
  string m="0";
  string h= s.substr(i+1, s.size()-i-1);
  return stoi(h);
}


/*
 *
 * The following function reads an allready existing
 * collection of estimated methylation frequency change events into the tree from a file called
 * p_change_info.sol. It can be used to restart a cancelled MCMC run in order to debug.
 *
 * When this function is invoked, there should allready be initialized transtions in tthe tree.
 */

void read_probability_change_from_file(vector <back_tree*> & node_container, const vector<double> & rate_acc)
{
  ifstream pc_stream("p_change_info.sol");
  string name_of_cell;
  string next_line;
  back_tree* current_node;
  int current_island_index;
  int kind_of_line=0;
  while (getline(pc_stream,next_line))
    {

      kind_of_line=2;
      for (auto &node : node_container)
	{
	  if(node->get_name()==next_line)
	    {
	      current_node=node;
	      kind_of_line=0;
	    }
	}
     
      if (is_methylation_island_index(next_line))
	kind_of_line=1;
    
   
      switch (kind_of_line)
	{
	case -1: cerr<<"read_probability_change_from_file: Could not parse line in p_change_info.sol"<<endl;
	  cerr<<"String read "<<next_line<<endl;
	  break;
	case 0:
	  break;
	case 1:
	  current_island_index=find_island_index(next_line);
	  break;
	case 2:
	  current_node->set_transition_at_index_by_string(next_line, rate_acc, current_island_index);
	  break;
	}
    }
	  
	  
	
}




/**************************************************Command Line Input*********************************/

/*
 * The following funciton searches the array argv for the string comp and returns the string following it turned into a number,
 * if com is found. Else it returns default_param.
 */
const double set_double_param(const int argc, char ** argv, const  string & comp, const double default_param)
{
  double param=default_param;
  for (int i=0; i<argc; i++)
    if (argv[i]==comp)
      {
	if (i==argc-1) cerr<<"parameter missing before end of input while using "<<comp<<endl;
	string params=argv[i+1];
	if (!string_is_number(params)) cerr<<"parameter is NAN while setting "<<comp<<endl;
	param=stof(params);
	break;
      }
  return param; 
}

/*The following funciton is used to take start parameters from the command line.
 *
 * Current paramters supported:
 * invp: To set the parameter for the fraction of invariant sites add -set_invp [number] while calling methalytion_simulator
 * log_shape: To set the parameter for the fraction of invariant sites add -set_log_shape [number] while calling methalytion_simulator
 * log_pcrate: Use -set_log_pcrate [number] to set this parameter.

 * Parameters:
 * argv: Contains the command line arguments.
 * argc: length of argv.
 * invp: is the invariant site rate in use in maind.
 * log_shape: is the log_shape parameter in use in main.
 * log_pcrate: is the log_pcrate parameter in use in main.
 */
void set_start_params(const int argc,char ** argv, double & invp, double & log_shape, double & log_pcrate)
{
  string comp="-set_invp";
  invp=set_double_param(argc, argv, comp, invp);

  comp="-set_log_shape";
  log_shape=set_double_param(argc, argv, comp, log_shape);

  comp="-set_log_pcrate";
  log_pcrate=set_double_param(argc, argv, comp, log_pcrate);
}



/*
 * The following function reads the command line input to find whether the default tree should be used or a custom tree.
 * If a custom tree should be used this tree is built from a Newick string from the command line. Tree is stored in
 * node container and first node in node container is the root of the tree.
 *
 * Format is "-use_custom_tree" followed by the tree in Newick format.
 *
 */

void set_starting_tree(const int argc, char ** argv, vector<back_tree*> & node_container)
{
  string comp="-use_custom_tree";
   for (int i=0; i<argc; i++)
    if (argv[i]==comp)
      {
	if (i==argc-1)
	  {
	    cerr<<"parameter missing before end of input while using "<<comp<<endl;
	    exit(1);
	  }
	string Newick=argv[i+1];
	tree_from_string(Newick,node_container);
	return;
      }  
}



/*
 * The following function reads in from the command line to determine whether certain paramaters should
 * not be estimated but set at their starting value.
 *
 * Format is "-disable_modes" followed by a string of unsigned ints  corresponding to a string
 * of modes to be disabled.
 *
 * 0 corresponds to sampling of new log branch lengths.
 * 1 corresponds to sampling of new log shape params.
 * 2 corresponds to sampling of new log IWE rates.
 * 3 corresponds to sampling of new events.
 * 4 corresponds to sampling of new invariant probabilities.
 * 5 corresponds to sampling of new probabilites in the root.
 */

const vector<bool> disable_modes(const int argc, char ** argv, const vector<bool> & default_modes)
{
  vector<bool> modes=default_modes;
  string comp="-disable_modes";
  for (int i=0; i<argc; i++)
    if (argv[i]==comp)
      {
	if (i==argc-1) cerr<<"parameter missing before end of input while using "<<comp<<endl;
	++i;
	string new_mode_string=argv[i];
	while (string_is_number(new_mode_string)&&i<argc)
	  {
	    int new_mode=stoi(new_mode_string);
	    if (new_mode<0 || new_mode>=modes.size()) cerr<< "Warning, no mode: " << new_mode<< endl;
	    modes[new_mode]=false;
	    ++i;
	    if (i<argc) new_mode_string=argv[i];	    
	  }
      }
  return modes;
}

/*
 * The following function determines on whether IWE transitions are read from P_change_info.sol.
 * It searches the command line input for -use_custom_transitions and if found returns True, else false.
 */

const bool use_transition_solution(const int argc, char ** argv)
{
  string comp="-use_custom_transitions";
  for(int i=0;i<argc;++i) if (argv[i]==comp) return true;
  return false;
}
