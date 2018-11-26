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

#ifndef FUNCTIONOBJECT_H
#define FUNCTIONOBJECT_H

#include <boost/math/distributions.hpp>
#include <random>

using namespace std;

class RealDistributionObject
{
public:
  virtual  double operator()(mt19937 &re)=0;
  virtual const double density( double x)=0;
  
};



class RealProposalDistributionObject
{
public:
  virtual  double operator()(const double current_value,mt19937 &re)=0;
};

class LogNormalProposal: public RealProposalDistributionObject
{
private:
  double standard_deviation;
  normal_distribution<double> distribution;
public:
  LogNormalProposal(double sd): standard_deviation(sd), distribution(0,sd) {};
  double operator()(const double current_rate, mt19937 & re);
 
};

/*
 * The following Proposal distribution uses several proposal distributions at once by 
 * having a fixed rate for each of them and using them with a probability proprtional to the
 * rate.
 */
class MultiProposal: public RealProposalDistributionObject
{
private:
  vector<RealProposalDistributionObject *> proposal_objects;
  vector<double> proposal_rates;
  discrete_distribution<int> which_dist;
public:
  // If no rates are specified it uses all distributions with egal frequency.
  MultiProposal(const vector<RealProposalDistributionObject*>  proposal_objects):proposal_objects(proposal_objects),proposal_rates(proposal_objects.size(), 1.0), which_dist(proposal_rates.begin(), proposal_rates.end()) {};

  MultiProposal(const vector<RealProposalDistributionObject*>  proposal_objects, const vector<double> proposal_rates):proposal_objects(proposal_objects),proposal_rates(proposal_rates), which_dist(proposal_rates.begin(), proposal_rates.end()) {};
   
  double operator()(const double current_value, mt19937 & re)
  {
    return (*proposal_objects[which_dist(re)])(current_value, re); 
  }
 
};

class BoundedNormalProposal: public RealProposalDistributionObject
{
private:
  double lower_bound;
  double upper_bound;
  double standard_deviation;
  normal_distribution<double> distribution;  
public:
  BoundedNormalProposal(double lb, double ub, double sd ): lower_bound(lb), upper_bound(ub), standard_deviation(sd), distribution(0,sd) {};
  double operator()(const double current_rate, mt19937 & re);
 
};

class BetaDistribution: public RealDistributionObject
{
  double alpha;
  double beta;
public:  
  BetaDistribution(double a, double b): alpha(a), beta(b) {};
  double operator()(mt19937 & re);

  /*
   * The following funciton creates a new proposal beta distribution by 
   * perturbing alpha with RealProposalDistributionObject 
   * a and beta with RealProposalDistributionObject b
   */
  BetaDistribution propose_new(RealProposalDistributionObject &a,RealProposalDistributionObject &b, mt19937 & re);

  /*
   * The following is the density function.
   */
  const double density(double x);
};


class RealFunctionObject
{
public:
  virtual  double operator()(const double current_value)=0;
  virtual  double logratio(const double a, const double b)=0;
};

class LogNormalPriorDensity : public  RealFunctionObject
{
  boost::math::normal_distribution<> dist;
  double m;
  double sigma;
public:
  LogNormalPriorDensity (): dist(0,1), m(0), sigma(1){}
  LogNormalPriorDensity (const double mean, const double sd): dist(mean,sd), m(mean), sigma(sd){}
  double operator()(const double x);
  double logratio(const double new_param, const double old_param)
  {
    return ((log(old_param)-m)*(log(old_param)-m)-(log(new_param)-m)*(log(new_param)-m))/(2*sigma*sigma);
  }
};


class UniformPriorDensity : public  RealFunctionObject
{
public:
  double operator()(const double x) {return 1.0;}
  double logratio(const double a, const double b) {return 0;}
};


class Dirichlet111DistributionObject
{
  exponential_distribution<double> helper;
public:
  Dirichlet111DistributionObject(): helper(1) {};
 
  const vector<double> operator()(mt19937 & re)
  {
    vector<double> output(3,0.0);
    double sum=0;
    for (unsigned i=0; i< 3; i++)
      {
	output[i]=helper(re);
	sum+=output[i];
      }
    for (auto &i :output)
      i/=sum;
    return output;
  }
};


/*
 * This object stores the current value of the log_shape parameter, alpha, and can be use to produce rate_accs from this value via get_rate_acc.
 */
class GammaRateConstruction
{
  friend class functionobject_tests;
  double alpha;
  boost::math::gamma_distribution<> gd;
  boost::math::gamma_distribution<> helper_gd;
  int n_rates;
public:
  GammaRateConstruction(double a, int n): alpha(a),  gd(a), helper_gd(a+1), n_rates(n) {};
  void set_alpha(double a)
  {
    alpha=max(pow(10, -7),min(a, 1000.0));
    gd= boost::math::gamma_distribution<>(max(pow(10, -7),min(a, 1000.0)));
    helper_gd= boost::math::gamma_distribution<>(max(pow(10, -7),min(a, 1000.0))+1);
  }
  const double get_alpha() {return alpha;}
  vector <double> get_rate_acc();
};




struct DistObjects
{
  DistObjects(): gamma_rate_construction(1,3){};
  RealProposalDistributionObject* proposal_for_branch_length;
  RealFunctionObject* length_prior;

  RealProposalDistributionObject* proposal_for_gamma_shape;
  GammaRateConstruction gamma_rate_construction;
  RealFunctionObject* gamma_shape_prior;

  RealProposalDistributionObject*pcrateproposal;
  RealFunctionObject* pcrate_prior;

  RealProposalDistributionObject*invariant_probability_poposal;
  RealFunctionObject* invariant_probability_prior;

  Dirichlet111DistributionObject dirichlet111;
  
};

#endif
