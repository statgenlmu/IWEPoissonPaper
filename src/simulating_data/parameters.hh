
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

