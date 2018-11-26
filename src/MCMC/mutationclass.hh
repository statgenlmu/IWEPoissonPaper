#ifndef MUTATION_TYPE_H
#define MUTATION_TYPE_H

#include <string>
#include <vector>
#include <random>

using namespace std;


class cpg_island
{
public:
  cpg_island();
  cpg_island(string& s);
  ~cpg_island();

  /*
   * states contains the methylation states the island is in. The possible states are
   * unmethylated partially methylated and methylated. Their numeric representation is 
   * 0, 1 , 2 respectively. NAs stand for missing data.
   */

  std::vector <int> states;



  /*
   * states_as_llhd contains the states as in above, states_as_llhd[i][j] is the states
   * of site i expressed as a llhd, with 
   * methylated meaning the vector (0,0,1)
   * partially methylated meaning the vector (0,1,0)
   * unmethylated meaning the vector (1,0,0)
   * NA meaning the vector (1,1,1)
   */
  std::vector<vector<double>> states_as_llhood;
  
  double change_rate;
  double m_rate;
  double r;//probability of methylation/non methylation

  const double size()
  {
    return states.size();
  }




  const vector<vector<double>> & get_llhd()
  {
    return states_as_llhood;
  }

  
  /*
   *The following function picks one site at random.
   *Then it changes the methylation state at this site with probability r
   *to methylated and with probability 1-r to unmethylated
   */
  void methylation_event(mt19937 & re);

  /*
   *The following function takes the parameters of the global beta distribution and changes 
   *r by drawing from it
   */
  void probability_change_event(double gw_alpha, double gw_beta, mt19937 & re);

  /*
   * The following takes two g_w rate parameters and determines which quantile of the 
   * genome wide rate the methylation states belong to. Quantiles range from 0 to n_quantiles minus 1.
   */
  int which_quantile(const double gw_alpha, const double gw_beta, const int n_quantiles);
  

  
};
#endif
