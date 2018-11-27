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


#ifndef MCMC_STEP_H
#define MCMC_STEP_H

using namespace std;

class back_tree;
struct DistObjects;


int MCMC_step(vector<back_tree *> & node_container,DistObjects & dists,vector<double> & rate_acc,double & pcrate, double & prob_inv,const vector<bool> & which_steps,mt19937 & re);


const double log10faculty(const int n);
const double log10poisdensity(const double l, const int n);
const double log10llhdratio(const pair<double, int> & a,const pair<double, int> & b);
const double log10llhdratio(const vector<pair <double, int>> &  a,const vector<pair <double, int>> &  b);
const pair <double, int> meanlhd(const vector<pair <double, int>> & llhds);

/*
 * Returns the ratio c=a/b as a pair of double and int. The value represented by this pair is  c.first * 10^(-100*c.second)
 * Here a and b are pairs representing probabilities, consiting of a leading term and a scale factor,
 * the actual value is of the probability represented by a is a.first * 10^(-100*a.second)
 */
const pair<double,int> lhd_ratio(const pair<double, int> & a,const pair<double, int> & b);


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

void sample_new_pcevent_tree(const vector<back_tree *> & node_container, DistObjects & dists,const vector<double> & rate_acc, const double pcrate, const double prob_inv, mt19937 & re);


#endif
