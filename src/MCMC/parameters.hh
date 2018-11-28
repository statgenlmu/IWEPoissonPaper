
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

#ifndef PARAMETERS_H
#define PARAMETERS_H


/*
 * The following function looks for a tree in Newick format in the file tree.sol. 
 * If it finds this tree it fills node_container with references to the nodes so that
 * node_container[0] is the root of the tree.
 */

void initialize_tree_from_file(vector<back_tree*> & node_container);



/*
 * This function parses files of the format CPG_ISLANDS_celltype.txt for each celltype name 
 * in the tree and initializes states for MCMC compuation.
 */

void set_methylation_states_from_files(vector<back_tree*> & node_container);

/*
 * Finds the name of a root of a sub tree by parsing a Newick string s
 */

std::string find_newick_name(std::string s);

/*
 * Finds the distance of an edge adjacent to the root by parsing a Newick string s
 */

double find_newick_distance(string s);

/* 
 * The following function parses a Newick string for two subtree strings.
 */

vector<string> find_sub_tree_strings(string s);

/*
 * The following function parses a string s in Newick format and builds a binary tree from it and then 
 * returns the root.
 * t_coals are saved not as absolute times but distances of the edge above to the root. 
 * This is then converted to absoluete heights later.
 */
back_tree * tree_from_string(string s,vector<back_tree*> & node_container);


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
void set_start_params(const int argc, char ** argv, double & invp, double & log_shape, double & log_pcrate);


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

const vector<bool> disable_modes(const int argc,char ** argv, const vector<bool> & default_modes);



/*
 *
 * The following function reads an allready existing
 * collection of estimated methylation frequency change events into the tree from a file called
 * p_change_info.sol. It can be used to restart a cancelled MCMC run in order to debug.
 *
 * When this function is invoked, there should allready be initialized transtions in tthe tree.
 */

void read_probability_change_from_file(vector <back_tree*> & node_container, const vector<double> & rate_acc);


/*
 * The following function determines on whether IWE transitions are read from P_change_info.sol.
 * It searches the command line input for -use_custom_transitions and if found returns True, else false.
 */

const bool use_transition_solution(const int argc, char ** argv);


/*
 * The following function reads the command line input to find whether the default tree should be used or a custom tree.
 * If a custom tree should be used this tree is built from a Newick string from the command line. Tree is stored in
 * node container and first node in node container is the root of the tree.
 *
 * Format is "-use_custom_tree" followed by the tree in Newick format.
 *
 */

void set_starting_tree(const int argc, char ** argv, vector<back_tree*> & node_container);


#endif

