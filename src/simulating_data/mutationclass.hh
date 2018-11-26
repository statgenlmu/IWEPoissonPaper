   /* This file is part of IWEinferrence.

    IWEinferrence is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IWEinferrence is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with IWEinferrence.  If not, see <https://www.gnu.org/licenses/>. */

MUTATION_TYPE_H
#define MUTATION_TYPE_H

#include <string>
#include <vector>
#include <random>

using namespace std;

class mutation_type
{
public:
  mutation_type(std::string n,  
		int sev, 
		double f_f,
		double f_h,
		double r): name(n), severity(sev), fitness_effect_f(f_f), fitness_effect_h(f_h), rate(r){}
  mutation_type();
  
  ~mutation_type();

  const std::string get_name(void);
  const int get_severity(void);
  const double get_fitness_effect_f(void);
  const double get_fitness_effect_h(void);
  const double get_rate(void);

  void set_severity(int);
  void set_fitness_effect_f(double);
  void set_fitness_effect_h(double);
  void set_rate(double);
  void set_name(std::string);


  std::string name;  
  int severity; // 0 for neutral, 1 for secondary driver, 2 for primary driver
  double fitness_effect_f;
  double fitness_effect_h;
  double rate;
  double CPG_influence;

};


class cpg
{
public:

  bool fitness_affect;//denotes wheter the cpg site may have an effect on fitness
  double base_m_rate; //denotes baserate of methylation
  double m_rate_change; //denotes changing of rates via methylation

  void site_rate_change();

};


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
