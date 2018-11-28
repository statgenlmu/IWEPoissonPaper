 
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

#include<iostream>
#include "functionobjects.hh"
//#include <boost/math/distributions/gamma_distribution.hpp>

double BetaDistribution::operator()(mt19937 & re)
{
  boost::math::beta_distribution<> distr(alpha, beta);
  uniform_real_distribution<double> uni_dist(0,1);
   
  return boost::math::quantile(distr,uni_dist(re));
}

BetaDistribution BetaDistribution::propose_new(RealProposalDistributionObject &a,RealProposalDistributionObject &b, mt19937 & re)
{
  return BetaDistribution(a(alpha,re),b(beta,re));
}

const double BetaDistribution::density(double x)
{
  boost::math::beta_distribution<> distr(alpha, beta);
  return boost::math::pdf(distr,x);
}


double LogNormalProposal::operator()(const double current, mt19937 & re)
{
  if (current>1e100||current<1e-100)
    {
      cerr<<"LogNormalProposal::operator(): Proposal is getting very small or large"<<endl;
      return  log(current)>0 ? pow(10,100):pow(10,-100);
    }
  return exp(distribution(re)+log(current));
}

double BoundedNormalProposal::operator()(const double current, mt19937 & re)
{
  double prop=current+distribution(re);
  while (prop<lower_bound || prop>upper_bound) //mirror the proposition until it is within bounds
    {
      if (prop>upper_bound) prop=2*upper_bound-prop;
      if (prop<lower_bound) prop=2*lower_bound-prop;
    }
      
  return prop;
}


double LogNormalPriorDensity::operator()(const double x)
{
  return boost::math::pdf(dist,log(x));
}



vector <double> GammaRateConstruction::get_rate_acc()
{

  vector<double>rate_acc(n_rates,0);
    
 
  double upperq=1.0/n_rates;

  double prev_limit=0;
  double next_limit=boost::math::quantile(gd,upperq);

  for (int i=0; i<n_rates;++i)
    {
      // Note for the next step that due to the form of the odf, increasing the form parameter by 1
      // exactly calculates the mean over an interval. helper gd is exactly increased by one.
      rate_acc[i]=boost::math::cdf(helper_gd, next_limit)-boost::math::cdf(helper_gd, prev_limit);
      upperq=min(0.999999,(1.0*(i+2))/n_rates);
      prev_limit=next_limit;
      next_limit=boost::math::quantile(gd,upperq);
    }

  double sum=0;

  for (const auto &i : rate_acc) sum+=i;
  for (auto &i : rate_acc) i*=n_rates/sum;

  return rate_acc;   
}
