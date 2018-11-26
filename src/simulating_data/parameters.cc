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


#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <cassert>
#include <regex>
#include "back_tree.hh"
#include "parameters.hh"

using namespace std;


back_tree* tree_from_string(string s, vector<back_tree*> & node_container, vector<double> & log_branch_lengths);
bool isint(string s)
{
  std::regex e ("\\d+");
  if (std::regex_match (s,e)) return true;
  return false;
}


bool isdouble(string s)
{
  std::regex e ("-?\\d+(\\.\\d+)?");
  if (std::regex_match (s,e)) return true;
  return false;
}



/*
 * take the parameter struct and fill it from the command line.
 * The command line should contain
 * -island_sizes  followed by unsinged island sizes
 * -log_branch_length followed by 13 log branch lengths which are double
 * -log shape, follwoed by a singledouble variable
 * -log_pc_rate followed by a single double variable
 */
void set_parameters_from_command_line(Parameter_Struct & params, vector<back_tree*> & node_container,int argc , char ** argv)
{

  for (int i=1; i< argc;)
    {
      string comp=argv[i];

      if (comp=="-island_sizes")
	{
	  i++;
	  if (i==argc-1) cerr<<"Error reading from commandline: island_sizes are too short."<<endl;
	  string input=argv[i];
	  if (!isint(input)) cerr<<"Error reading from commandline: island_sizes are too short or not unsinged"<<endl;
	  while (isint(input)&&i!=argc)
	    {
	      params.island_sizes.push_back(stoi(input));
	      i++;
	      if (i<argc)
		input=argv[i];
	    }
	}
       if (comp=="-log_branch_lengths")
	{
	  i++;
	  if (i>=argc-12) cerr<<"Error reading from commandline: log branch ends abruptly."<<endl;
	  string input=argv[i];
	  if (!isdouble(input)) cerr<<"Error reading from commandline: log branch are too short."<<endl;
	  while (isdouble(input)&&i!=argc)
	    {
	      params.log_branches.push_back(stof(input));
	      i++;
	      if (i<argc)
		input=argv[i];
	    }
	}
       if (comp=="-log_shape")
	{
	  i++;
	  if (i>=argc) cerr<<"Error reading from commandline: log shape ends abrubtly"<<endl;
	  string input=argv[i];
	  if (!isdouble(input)) cerr<<"Error reading from commandline: log shape is not double"<<endl;
	  params.log_shape=stof(input);
	  i++;
	}
       if (comp=="-log_pc_rate")
	 {
	   i++;
	   if (i>=argc) cerr<<"Error reading from commandline: pc rate ends abrubtly"<<endl;
	   string input=argv[i];
	   if (!isdouble(input)) cerr<<"Error reading from commandline: pc rate is not double"<<endl;
	   params.log_pc_rate=stof(input);
	   i++;
	 }
       if (comp=="-p_inv")
	  {
	   i++;
	   if (i>=argc) cerr<<"Error reading from commandline: p invariant ends abrubtly"<<endl;
	   string input=argv[i];
	   if (!isdouble(input)) cerr<<"Error reading from commandline: fraction of invariant sites is not double"<<endl;
	   params.p_inv=stof(input);
	   i++;
	 }
       if (comp=="-use_custom_tree")
	 {
	   i++;
	   if (i>=argc) cerr<<"Error reading from commandline: custom tree ends abrubtly"<<endl;
	   string input=argv[i];
	   if (params.log_branches.size()>0)
	     {
	       cerr<<"Warning, cannot set branch lengths in Newick and as sepparate command line input. Please specify in Newick."<<endl;
	       exit(1);
	     }
	   tree_from_string(input,  node_container, params.log_branches);
	   //set_tree_from_Newick(node_container,input);
	   i++;
	     
	 }
	 
    }
      
	  
  
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
  if (s.find(":")==string::npos)
    {
      cerr<<"find_newick_distance: Branch length is missing"<<endl;
      exit(1);
    }
  for ( n=s.size(); n>0&&s[n-1]!=':';--n);
  for (; n<s.size();++n) r.push_back(s[n]);
  if (r=="")
    {
      cerr<<"find_newick_distance: Branch length is missing"<<endl;
      exit(1);
    }
  if (! isdouble(r))
    {
      cerr<<"find_newick_distance: Branch length is NAN"<<endl;
      exit(1);
    }
  cout<<r<<endl;
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
back_tree* tree_from_string(string s, vector<back_tree*> & node_container, vector<double> & log_branch_lengths)
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
  log_branch_lengths.push_back(log(find_newick_distance(s)));
  node_container.push_back(sub_tree_root);

  if (s.find(")") != string::npos)
    {
      vector<string> sub_tree_strings=find_sub_tree_strings(s);
      for (auto i:sub_tree_strings)
	sub_tree_root->add_to_children(tree_from_string(i, node_container, log_branch_lengths));
    }
  
  return sub_tree_root;
}




/*
 * The following function looks for a tree in Newick format in the file tree.sol. 
 * If it finds this tree it fills node_container with references to the nodes so that
 * node_container[0] is the root of the tree.
 */

void initialize_tree_from_file(vector<back_tree*> & node_container, vector<double> & log_branches)
{
  ifstream tree_stream("tree.sol");
  string tree_string;
  if (!getline(tree_stream, tree_string))
    {
      cerr<<"No tree found, \"tree.sol\" empty or not existent. Default configuration was used."<<endl;
      return;
    }
  tree_from_string(tree_string,node_container, log_branches);
}
