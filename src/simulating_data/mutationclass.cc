
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
#include "mutationclass.hh"
#include <string>
#include <random>
#include <boost/math/distributions/beta.hpp>

using namespace std;
using boost::math::beta_distribution;
mutation_type::mutation_type()
{
  severity=0;
  fitness_effect_h=0.0;
  fitness_effect_f=0.0;
  name = " ";

}

mutation_type::~mutation_type()
{
  
}
  


const std::string mutation_type::get_name()
{
  return name;
}

const int mutation_type::get_severity()
  {
    return severity;
  }

const double mutation_type::get_fitness_effect_f()
  {
    return fitness_effect_f;
  }

const double mutation_type::get_fitness_effect_h()
  {
    return fitness_effect_h;
  }

const double mutation_type::get_rate()
{
  return rate;
}

void mutation_type::set_name(std::string s)
{
  name=s;
}


void mutation_type::set_severity(int s)
 {
   if( (s<0)||(s>2))
     std::cout << "Problem: No such mutation severity" << std::endl;
   else
     severity=s;

 }

void mutation_type::set_fitness_effect_f(double f)
 {
   fitness_effect_f=f;

 }

void mutation_type::set_fitness_effect_h(double h)
 {
   fitness_effect_h=h;

 }

void mutation_type::set_rate(double r)
{
  rate=r;
}


void cpg::site_rate_change()
{
  base_m_rate*=m_rate_change;
}



cpg_island::cpg_island(){}
cpg_island::~cpg_island(){}

cpg_island::cpg_island(string & s)
{
  for (unsigned i=0; i< s.size(); ++i)
    {
      if (s[i]=='0') states.push_back(0);
      if (s[i]=='1') states.push_back(1);
      if (s[i]=='2') states.push_back(2);
      if (s[i]=='?') states.push_back(NAN);
      if (s[i]=='N') states.push_back(NAN);
    }

  for (unsigned i=0; i< s.size(); ++i)
    {
      if (s[i]=='0') states_as_llhood.push_back({1,0,0});
      if (s[i]=='1') states_as_llhood.push_back({0,1,0});
      if (s[i]=='2') states_as_llhood.push_back({0,0,1});
      if (s[i]=='?') states_as_llhood.push_back({1,1,1});
      if (s[i]=='N') states_as_llhood.push_back({1,1,1});
    }

  
}




/*
 *The following function picks one site at random.
 *Then it changes the methylation state at this site with probability r
 *to methylated and with probability 1-r to unmethylated
 */

void cpg_island::methylation_event(mt19937 & re)//TODO:Unit test
{
  uniform_int_distribution <int> which_cpg (0, states.size());
  int w_cpg=which_cpg(re);
  bernoulli_distribution state_change(r);
  states[w_cpg]=state_change(re);
}

/*
 *The following function takes the parameters of the global beta distribution and changes 
 *r, the probability of a respective methylation state at an individual CPG after a 
 *substitution event at that CPG, by drawing from it
 */

void cpg_island::probability_change_event(double gw_alpha, double gw_beta, mt19937 & re) //TODO:Unit test
{
  beta_distribution <>  new_m_density(gw_alpha,gw_beta);
  uniform_real_distribution <double> rfu(0.0,1.0);
  double rand_from_uni=rfu(re);
  
  double new_fraction=quantile(new_m_density, rand_from_uni);	   
  r=new_fraction;
}

/*
 * The following takes two g_w rate parameters and determines which quantile of the 
 * genome wide rate the methylation states belong to. Quantiles are indexed and indexes
 * range from 0 to n_quantiles minus 1. The function returns the index.
 */
int cpg_island::which_quantile(const double gw_alpha, const double gw_beta, const int n_quantiles)
{
  int n_meth=0;
  for (auto i : states) n_meth+=i;
  double r_meth=n_meth;
  r_meth/=states.size();
  beta_distribution <>  m_density(gw_alpha,gw_beta);
  return min(double(n_quantiles-1),quantile(m_density, r_meth)*n_quantiles);
}
