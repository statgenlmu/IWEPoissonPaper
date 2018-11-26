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
  double logratio(const double a, const double b)
  {
    return ((log(a)-m)*(log(a)-m)-(log(b)-m)*(log(b)-m))/(sigma*sigma);
  }
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


class GammaRateConstruction
{

  double alpha;
  boost::math::gamma_distribution<> gd;
  boost::math::gamma_distribution<> helper_gd;
  int n_rates;
  
public:

  
  GammaRateConstruction(double a, int n): alpha(max(pow(10, -7),min(a, 1000.0))),  gd(alpha), helper_gd(alpha+1), n_rates(n) {};

  void set_alpha(double a) {
    alpha=max(pow(10, -7),min(a, 1000.0));
    gd= boost::math::gamma_distribution<>(max(pow(10, -7),min(a, 1000.0)));
    helper_gd= boost::math::gamma_distribution<>(max(pow(10, -7),min(a, 1000.0))+1);
  }

  const double get_alpha() {return alpha;}

  vector <double> get_rate_acc(mt19937 & re);
 
  


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

  Dirichlet111DistributionObject dirichlet111;
  
};

#endif
