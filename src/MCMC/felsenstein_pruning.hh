#ifndef FELSENSTEIN_H
#define FELSENSTEIN_H
#include<iostream>

using namespace std;

class island_partial_likelihoods_numscale{
  // partial likelihoods of the states
  // unmethylated, partially methylated and methylated
  // for the CpG sites of an island
private :
  vector<double> rate_acc; // the possible rate acceleration factors
  // for the rate of methylation changes at CpG sites
  unsigned scalefactor;
  vector<vector<vector<double>>> part_L;
  // part_L[s][r][x]/(10^(100*scalefactor)) is the partial likelihood
  // of state x (0: unmethylated, 1: partially methylated, 2: methylated)
  // at position s within the island, for the case that the rate
  // acceleration factor at that site is rate_acc[r];
public :
  island_partial_likelihoods_numscale(const vector<double> & rate_acc, unsigned scf,
				      const vector<vector<vector<double>>> & part_L):
    rate_acc(rate_acc), part_L(part_L){scalefactor=scf;}
  const vector<double> & get_part_L(unsigned s, unsigned r) const {return part_L[s][r];}
  void set_part_L(unsigned s, unsigned r, vector<double> v) {part_L[s][r] = v;}
  double get_scalefactor() {return scalefactor;}
  void set_scalefactor(unsigned s) {scalefactor=s;}
  void set_rate_acc(const vector<double> & rate_a) {rate_acc=rate_a;}
  void rescale(); // increases scalefactor by one and adapts oll part_L
  void rescale_if_necessary(); // rescales if maximum part_L is below 10^(-100)
  vector<double> get_rate_acc() const {return rate_acc;}
  unsigned nsites() const {return part_L.size();} // number of site in CpG island
  unsigned nratecat() const {return rate_acc.size();} // number of rate categories
};

class island_meth_changes {
  // when and how do methylation frequencies change in
  // a CpG island on a branch of the genealogy
private :  
  bool up_to_date_flag; // are all the partial likelihoods up to date?
  vector<double> change_time; // TODO: RECFACTOR CHANGE_TIME TO LIST<VECTOR<DOUBLE>> FOR BETTER EFFICIENCY
  // time points on a branch in which events take place that
  // change methylation frequencies on the island.
  // last entry (which may be the only one) is total length.
  // The time of the k-th IWE is change_time[k-1] and
  // the last change_time is the branch length.
  vector<island_partial_likelihoods_numscale*> pll;  // TODO: RECFACTOR PPL TO LIST<ISLAND_PARTIAL_LIKELIHOODS_NUMSCALE*> FOR BETTER EFFICIENCY
  // note that pll.back()==daughter;
  // The partial likelihoods in the mother node are *pll[0] and those in the k-th IWE are *pll[k]
  // Note that pll.size()==change_time()+1.
  island_partial_likelihoods_numscale * daughter; 
  vector<vector<double>> methyl_freq; // TODO: RECFACTOR PPL TO LIST<VECTOR<DOUBLE>> FOR BETTER EFFICIENCY
  // probabilities of unmethylated, partially methylated and methylated
  // (in this order)
  // in the beginning of the branch (methyl_freq[0]) 
  // after change at time change_time.
  // The probs starting at the k-th IWE 
  // The probs methyl_freq[k] hold from change_time[k-1] to change_time[k].
  // The vectors methyl_freq and change_time have the same size, and this size is n+1 if n is the number of IWEs.
  vector<vector<vector<vector<double>>>*> transprob; 
  // (*transprob[k])[g][i][j] will be the transition prob from state i to state j
  // assuming rate category g between IWEs k and k+1,
  // or between the parent node and IWE 1 if k=0.
  void update_transprobs(unsigned k);
  // update *transprob[k]
  void update_transprobs_alt(unsigned k);
  // update *transprob[k]
  // (alternative implementation)
public :
  island_meth_changes(const vector<double> & change_time, island_partial_likelihoods_numscale * daughter_pll,  
		      const vector<vector<double>> & methyl_freq);
  ~island_meth_changes(){
      for(unsigned i=0; i+1<pll.size(); ++i) delete pll[i]; 
      for(unsigned i=0; i<transprob.size(); ++i) delete transprob[i];
  }
  island_meth_changes(const island_meth_changes & imc){
      cerr << "copy constructor for island_meth_changes not yet implemented\n";
      exit(1);
  };
  island_meth_changes & operator=(const island_meth_changes & imc){
      cerr << "operator= for island_meth_changes not yet implemented\n";
      exit(1);
      return *this;
  };
  void update_all_transprobs() {for(unsigned i=0; i<transprob.size(); ++i) update_transprobs(i);}
  const vector<double> get_change_time() const {return change_time;}
  double get_change_time(unsigned i) const {return change_time[i];}
  unsigned get_n_intervals() const {return change_time.size();}
  const vector<vector<double>> & get_methyl_freq() const {return methyl_freq;}
  const vector<double> & get_methyl_freq(unsigned i) const {return methyl_freq[i];}
  const vector<double> & get_part_L(unsigned s, unsigned r) const {return pll[0]->get_part_L(s,r);}
  const vector<double> & get_part_L(unsigned k, unsigned s, unsigned r) const {return pll[k]->get_part_L(s,r);}
  double get_transprob(unsigned k, unsigned g, unsigned i, unsigned j) const {return (*transprob[k])[g][i][j];}
  island_partial_likelihoods_numscale * get_daughter() {return daughter;}
  void stretch_branch_length(const double factor);
  //stretches all change times by the constant factor factor.
  // Note that partial likelihoods are not automatically updated.
  void plus_one_meth_change(const vector<double> & meth_freq, const double r);
  // add a new island-wide methylation event at time r.
  void minus_meth_change(const unsigned k);
  // remove the k-th meth change event k. (counting starts with 1 here)
  void update_partial_llhs();
  // for a CpG island calculates contribution to partial likelihoods from a branch that is characterized in imc.
  // pll are the partial likelihoods at the daughter node at the other end at the branch
  void update_partial_llhs(unsigned n);
  // for a CpG island calculates contribution to partial likelihoods from a branch that is characterized in imc.
  // pll are the partial likelihoods at an intermediate node, where the n-th rate change on the branch take place.
  // pll[0] is the table of partial likelihoods at the parent node.
  // Note that it is assumed that the transition probabilities are up to date for the next node.
  // Note also that, if the rates between nodes n-1 and n have changed, first the transition probabilities between
  // nodes n and n+1 needs to be updated.
  void accept_proposal(island_meth_changes * prop, bool copydaughter);
  // The proposal *prop is accepted. pll[1] to pll[pll.size()-2] will then point to the 
  // previous pll tables on prop. The previous pll tables of *this are deleted and the
  // *prop will not have pointers its previous tables for global rate changes.
  // If copydaughter==true, the partial likelihoods at the daughter node are copied from *prop,
  // otherwise they are unchanged. In the latter case it is assumed (but not checked!) that they were the same anyway.
  void copy_to_proposal(island_meth_changes * prop) const;
  // *prop will have the same strucure as *this but the partial likelihood
  // at the daughter node will not be changed. 
  // The partial likelihoods of *prop are not automatically updated.
  // If needed, this must be done with prop->partial_llh()
  // It is assumed that the pll tables that already exist in *prop have 
  // the same dimensions as those in *this.
  const double get_scalefactor() const {return pll[0]->get_scalefactor();}
  const vector<double> get_last_meth_freq()
  {
    return methyl_freq.back();
  }
  void set_first_meth_freq(const vector<double> & new_f) {
    methyl_freq[0]=new_f;
    update_transprobs(0);
    if (get_n_intervals()>1) {
      update_transprobs(1);
      update_partial_llhs(1);
    } else {
      update_partial_llhs(0);
    }
  }
  void set_rate_acc(const vector<double> & rate_a) {
    for (auto &i: pll) i->set_rate_acc(rate_a);
    daughter->set_rate_acc(rate_a);
    update_all_transprobs();
    update_partial_llhs();
    is_up_to_date();
  }
  bool up_to_date() const {return up_to_date_flag;}
  void to_be_updated() {up_to_date_flag=false;}
  void is_up_to_date() {up_to_date_flag=true;}
};



class island_partial_logL{
// logarithms of partial likelihoods of the states
// unmethylated, partially methylated and methylated
// for the CpG sites of an island
private :
    vector<double> rate_acc; // the possible rate acceleration factors
    // for the rate of methylation changes at CpG sites
    vector<vector<vector<double>>> log_part_L;
    // log_part_L[s][r][x] is the log of the partial likelihood
    // of state x (0: unmethylated, 1: partially methylated, 2: methylated)
    // at position s within the island, for the case that the rate
    // acceleration factor at that site is rate_acc[r];
public :
    island_partial_logL(const vector<double> & rate_acc, const vector<vector<vector<double>>> & log_part_L) :
	rate_acc(rate_acc), log_part_L(log_part_L){};
    island_partial_logL(const vector<double> & rate_acc) :
	rate_acc(rate_acc) {};
    vector<double> get_log_part_L(unsigned s, unsigned r) const {return log_part_L[s][r];}
    vector<double> get_rate_acc() const {return rate_acc;}
    unsigned nsites() const {return log_part_L.size();} // number of site in CpG island
    unsigned nratecat() const {return rate_acc.size();} // number of rate categories
    void set_log_part_L(const vector<vector<vector<double>>> & v) {log_part_L=v;}
};


island_partial_logL island_branch_partial_logL(const island_partial_logL & pll, const island_meth_changes & imc);
// for a CpG island calculates contribution to log partial likelihoods from a branch that is characterized in imc.
// pll are the partial likelihoods at the daughter node at the other end at the branch

island_partial_logL island_branch_partial_logL(const island_partial_logL & pll, const island_meth_changes & imc, unsigned n);
// for a CpG island calculates contribution to log partial likelihoods from a branch that is characterized in imc.
// pll are the partial likelihoods at an intermediate node, where the n-th rate change on the branch take place.

#endif


