#include <iostream>
#include "mutationclass.hh"
#include <string>
#include <random>
#include <boost/math/distributions/beta.hpp>

using namespace std;
using boost::math::beta_distribution;



cpg_island::cpg_island(){}
cpg_island::~cpg_island(){}

cpg_island::cpg_island(string & s)
{
  for (unsigned i=0; i< s.size(); ++i)
    {
      if (s[i]=='0') states.push_back(0);
      if (s[i]=='1') states.push_back(1);
      if (s[i]=='2') states.push_back(2);
      if (s[i]=='?') states.push_back(-1);
      if (s[i]=='N') states.push_back(-1);
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


