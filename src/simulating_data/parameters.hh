

#ifndef PARAMETERS_H
#define PARAMETERS_H


/*
 * Struct containing the parameters of the simulation set by the user.
 */

struct Parameter_Struct
{
  Parameter_Struct():log_pc_rate(-8), p_inv(0) {}; 
  vector<unsigned> island_sizes;
  vector<double> log_branches;
  double log_shape;
  double log_pc_rate;
  double p_inv;
};


/*
 * take the parameter struct and fill it from the command line.
 * The command line should contain
 * -island_sizes  followed by unsinged island sizes
 * -log_branch_length followed by 13 log branch lengths which are double
 * -log shape, follwoed by a singledouble variable
 * -log_pc_rate followed by a single double variable
 */
void set_parameters_from_command_line(Parameter_Struct & params, vector<back_tree*> & node_container,int argc , char ** argv);

#endif

