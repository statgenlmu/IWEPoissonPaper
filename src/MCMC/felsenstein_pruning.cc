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

#include <cmath>
#include <vector>
#include <utility>
#include <random>
#include <algorithm>
#include <iostream>
#include "felsenstein_pruning.hh"


extern mt19937 re;
// --------------------- class island_meth_changes ----------------------------

// private:


void island_meth_changes::update_transprobs(unsigned k) {
    // update *transprob[k]
    // This implementation is presumably faster than the other one but
    // cannot be generalized to parent-dependent dynamics.
    if(k < 0 || k>=transprob.size()) {
	cerr << "k out of range in call of void update_transprobs(k)\n";
	exit(1);
    }
    if(k==0) {
	// no immediate transitions allowed in parent node
	vector<double> erate((*transprob[0]).size());
	const double t=get_change_time(0);
	const vector<double> rate_acc(pll[0]->get_rate_acc());
	for(unsigned i=0; i<erate.size();++i) {
	    erate[i] = exp(-rate_acc[i]*t);
	}  // FOR NOW ONLY F81 TYPE OF MODEL IS IMPLEMENTED. NEED FUNCTION FOR MATRIX EXPONENTIAL TO GENERALIZE.
	for(unsigned g=0; g<transprob[0]->size(); ++g) {
	    for(unsigned i=0; i<3; ++i) {
		for(unsigned j=0; j<3; ++j) {
		    (*transprob[0])[g][i][j] = (1-erate[g])*get_methyl_freq(0)[j];
		}
		(*transprob[0])[g][i][i] += erate[g];
	    }
	}
    } else {
	vector<double> erate(transprob[k]->size());
	const double t=get_change_time(k)-get_change_time(k-1);
	const vector<double> rate_acc(pll[k]->get_rate_acc());
	for(unsigned i=0; i<erate.size();++i) {
	    erate[i] = exp(-rate_acc[i]*t);
	}  // FOR NOW ONLY F81 TYPE OF MODEL IS IMPLEMENTED. NEED FUNCTION FOR MATRIX EXPONENTIAL TO GENERALIZE.
	
	vector<double> freqdiff(3); // changes in state frequencies
	for(unsigned i=0; i<3; ++i) freqdiff[i]=get_methyl_freq(k)[i]-get_methyl_freq(k-1)[i];
	for(unsigned g=0; g<pll[k]->nratecat(); ++g)
	    for(unsigned i=0; i<3; ++i) {
		for(unsigned j=0; j<3; ++j) {
		    (*transprob[k])[g][i][j] = (1-erate[g])*get_methyl_freq(k)[j];
		    if(freqdiff[i]<0 && freqdiff[j]>0) {
			// Here we use that the dynamics on the branch is (similar to F81) parent-independent.
                        // If a state changes on the branch, the instantaneous change in the node is not relevant.
                        // That is, the probability of a change from state x to state y is increased by the
                        // probability that there is no "type-forget" event on the branch times the probability
                        // that there is a change in the node from x to y.
			// Note that if the probability of state i becomes smaller at an island-wide methylation
			// change event, each site in state i has a probability of freqdiff[i]/f_i to change to one of
			// the other states, where f_i := get_methyl_freq(n-1)[i] is here the frequency of i state sites.
			// If both other states increase their probability, the sum of their differences
			// must be -freqdiff[i], and the probability that a changed state becomes j is then the fraction 
			// freqdiff[j]/(-freqdiff[i]). Thus, -freqdiff[i] cancels and the probability that a site that
			// was in i becomes j is just freqdiff[j]/f_i.
			// If, however, only the probability of j increases, all changed sites that were i before become j,
                        // and the probability to go from i to j in this way is just -freqdiff[i]/f_i.
                        // All this can be hidden by another single-site methylation change, so the increase
			// in the probability to change from i to j contains the factor erate[g], which is the probability
			// that no single-site methylation change affects the end of the next branch segment.
			(*transprob[k])[g][i][j] += erate[g]*(freqdiff[j] < -freqdiff[i] ? freqdiff[j] : -freqdiff[i])/get_methyl_freq(k-1)[i];
		    }
		}
		double f=erate[g];
		if(freqdiff[i]<0) f *= (1+freqdiff[i]/get_methyl_freq(k-1)[i]);
		(*transprob[k])[g][i][i] += f;
	    } 
    }
}


void island_meth_changes::update_transprobs_alt(unsigned k) {
    // update *transprob[k]
    if(k < 0 || k>=transprob.size()) {
	cerr << "k out of range in call of void update_transprobs(k)\n";
	exit(1);
    }
    if(k==0) {
	// no immediate transitions allowed in parent node
	vector<double> erate((*transprob[0]).size());
	const double t=get_change_time(0);
	const vector<double> rate_acc(pll[0]->get_rate_acc());
	for(unsigned i=0; i<erate.size();++i) {
	    erate[i] = exp(-rate_acc[i]*t);
	}  // FOR NOW ONLY F81 TYPE OF MODEL IS IMPLEMENTED. NEED FUNCTION FOR MATRIX EXPONENTIAL TO GENERALIZE.
	for(unsigned g=0; g<transprob[0]->size(); ++g) {
	    for(unsigned i=0; i<3; ++i) {
		for(unsigned j=0; j<3; ++j) {
		    (*transprob[0])[g][i][j] = (1-erate[g])*get_methyl_freq(0)[j];
		}
		(*transprob[0])[g][i][i] += erate[g];
	    }
	}
    } else {
	vector<double> erate(transprob[k]->size());
	const double t=get_change_time(k)-get_change_time(k-1);
	const vector<double> rate_acc(pll[k]->get_rate_acc());
	for(unsigned i=0; i<erate.size();++i) {
	    erate[i] = exp(-rate_acc[i]*t);
	}  // FOR NOW ONLY F81 TYPE OF MODEL IS IMPLEMENTED. NEED FUNCTION FOR MATRIX EXPONENTIAL TO GENERALIZE.
	
	vector<double> freqdiff(3); // changes in state frequencies
	for(unsigned i=0; i<3; ++i) freqdiff[i]=get_methyl_freq(k)[i]-get_methyl_freq(k-1)[i];
	for(unsigned g=0; g<pll[k]->nratecat(); ++g) {
	    vector<vector<double>> M(3,vector<double>(3,0)), P(3,vector<double>(3));
	    // M is transition matrix within IWE, P transition matrix to next IWE; result will be M*P
	    // We start with P:
	    for(unsigned i=0; i<3; ++i) {
		for(unsigned j=0; j<3; ++j) {
		    P[i][j] = (1-erate[g])*get_methyl_freq(k)[j];
		}
		P[i][i] += erate[g];
	    }
	    // Now for Matrix M:
	    bool onlyonegreater=false;
	    unsigned i,j,h; // i will be the focal one that is different from the others
	    if(freqdiff[0]>0) {
		if(freqdiff[1]>0) {
		    i=2; j=0; h=1;
		} else
		{  // freqdiff[1]<=0
		    if(freqdiff[2]>0) {
			i=1; j=2; h=0;
		    } else
		    { // freqdiff[2]<=0
			onlyonegreater=true;
			i=0; j=1; h=2; 
		    }
		}
	    } else { // freqdiff[0]<=0
		  if(freqdiff[1]>0) {
			if(freqdiff[2]>0) {
			     i=0; j=1; h=2;
		        } else {
			     onlyonegreater=true;
			     i=1; j=2; h=0;
		        }
	         } else { // freqdiff[1]<=0
			onlyonegreater=true;
			i=2; j=0; h=1;
		    }
	    }
	    if(onlyonegreater) {// freqdiff[i]>0, freqdiff[j]<=0, freqdiff[h]<=0 
		 M[i][i] = 1.0;
		 M[j][i] = -freqdiff[j]/get_methyl_freq(k-1)[j];
		 M[j][j] = get_methyl_freq(k)[j]/get_methyl_freq(k-1)[j];
		 M[h][i] = -freqdiff[h]/get_methyl_freq(k-1)[h];
		 M[h][h] = get_methyl_freq(k)[h]/get_methyl_freq(k-1)[h];
	    } else { // freqdiff[i]<=0, freqdiff[j]>0, freqdiff[h]>0
		 M[i][i] = get_methyl_freq(k)[i]/get_methyl_freq(k-1)[i];
		 M[i][j] = freqdiff[j]/get_methyl_freq(k-1)[i];
		 M[i][h] = freqdiff[h]/get_methyl_freq(k-1)[i];
		 M[j][j] = 1.0;
		 M[h][h] = 1.0;
	    }
	    // now the matrix product M*P:
	    for(unsigned i=0; i<3; ++i) {
	        for(unsigned j=0; j<3; ++j) {
		   (*transprob[k])[g][i][j]=0.0;
		   for(unsigned h=0; h<3; ++h) (*transprob[k])[g][i][j] += M[i][h]*P[h][j];
		}
	    }
	}
    }
}


// public:

island_meth_changes::island_meth_changes(const vector<double> & c_time, island_partial_likelihoods_numscale * daughter_pll,  
					 const vector<vector<double>> & methyl_freq) : 
    change_time(c_time), pll(0), methyl_freq(methyl_freq), transprob(c_time.size()) {
    if(c_time.size()!=methyl_freq.size()) {
	cerr << "error: c_time and methyl_freq must have the same length when given to " 
	     << "island_meth_changes::island_meth_changes const vector<double> & c_time, "
	     << "island_partial_likelihoods_numscale * daughter_pll, "
	     << "const vector<vector<double>> & methyl_freq)\n";
    }
    for(unsigned i=0; i<change_time.size(); ++i) {
	island_partial_likelihoods_numscale * p = new island_partial_likelihoods_numscale(*daughter_pll);
	pll.push_back(p);
	transprob[i]=new vector<vector<vector<double>>>(daughter_pll->get_rate_acc().size(),vector<vector<double>>(3,vector<double>(3)));
	update_transprobs(i);
    }
    pll.push_back(daughter_pll);
    daughter=daughter_pll;
    update_partial_llhs();
}

void island_meth_changes::stretch_branch_length(const double factor)
{
  for (auto & j:change_time)
    j*=factor;
  update_all_transprobs();
  update_partial_llhs(); 
  is_up_to_date();
}


void island_meth_changes::plus_one_meth_change(const vector<double> & meth_freq, const double r) {
    // add a new island-wide methylation event at time r.
    
    const unsigned n=methyl_freq.size()-1; // n is the current number of IWEs
    unsigned k=1;
    while(change_time[k-1]<r && k<n+2) ++k;
    if(k==n+2) {
	cerr << "error in void island_meth_changes::plus_one_meth_change(const vector<double> & meth_freq, const double r):"
	     << " r is beyond branch length\n";
	exit(1);
    }
    // r now is between change_time[k-2] and change_time[k-1], and will be new change_time[k-1].
    // That is, the newly inserted IWE will be the k-th IWE
    if(!(change_time[k-1]>r)) {
	cerr << "warning in void island_meth_changes::plus_one_meth_change(const vector<double> & meth_freq, const double r):"
	     << " new island-wide methylation events at exact same time as existing one or branch end\n";
    }
    methyl_freq.push_back(methyl_freq[n]);
    transprob.push_back(transprob[n]);
    // cout << n << "      " << change_time.size() << endl;
    change_time.push_back(change_time[n]); 
    // methyl_freq, change_time and transprob should now have size n+2
    for(unsigned j=n; j>=k; --j) {    // TODO: REFACTOR VECTORS TO LISTS TO GET RID OF THIS LOOP
	methyl_freq[j+1]=methyl_freq[j];
	transprob[j+1]=transprob[j];
	change_time[j]=change_time[j-1];
    }
    // if(transprob.size()!=n+2) {
    // 	cerr << "transprob.size() is " <<  transprob.size() << " but should be " << n+2 << endl;
    // 	exit(1);
    // }
    methyl_freq[k]=meth_freq;
    change_time[k-1]=r;
    pll.push_back(daughter);
    for(unsigned j=pll.size()-2; j>k; --j) {   // TODO: REFACTOR VECTORS TO LISTS TO GET RID OF THIS LOOP
	pll[j]=pll[j-1];
    }
    pll[k]=new island_partial_likelihoods_numscale(*daughter);

    transprob[k]=new vector<vector<vector<double>>>(daughter->get_rate_acc().size(), vector<vector<double>>(3,vector<double>(3)));
    update_transprobs(k);
    if(k>0) update_transprobs(k-1);
    if(k+1<transprob.size()) update_transprobs(k+1);
    update_partial_llhs(k < n ? k+1 : n+1);
}

	
void island_meth_changes::minus_meth_change(const unsigned k) {
    // remove the k-th meth change event. (counting starts with 1 here)
    const unsigned n=methyl_freq.size()-1; // n is the current number of IWEs
    if(k>0 && k<=n) {
	delete transprob[k];
	for(unsigned i=k; i<n; ++i) {   // TODO: REFACTOR VECTORS TO LISTS TO GET RID OF THIS LOOP
	    methyl_freq[i]=methyl_freq[i+1];
	    transprob[i]=transprob[i+1];
	    change_time[i-1]=change_time[i];
	}
	change_time[n-1]=change_time[n];
	methyl_freq.resize(n);
	transprob.resize(n);
	change_time.resize(n);
	if(pll.size()!=n+2) {cerr << "wrong size of pll in island_meth_changes::minus_meth_change\n"; exit(1);}
	delete pll[k];
	for(unsigned i=k; i<n; ++i) { // TODO: REFACTOR VECTORS TO LISTS TO GET RID OF THIS LOOP
	    pll[i]=pll[i+1];
	}
	pll[n]=daughter;
	pll.resize(n+1);
	// number of IWEs is now n-1
	if(k>0) update_transprobs(k-1);
	if(k<n) update_transprobs(k);
	update_partial_llhs(k < n ? k : n-1);
    } else {
      cerr << "island_meth_changes::minus_meth_change(const unsigned): attempt to delete no-existing island-wide methylation event"<<endl;
    }
}

void island_meth_changes::update_partial_llhs() {
// for a CpG island calculates contribution to partial likelihoods from a branch that is characterized in imc.
// pll are the partial likelihoods at the daughter node at the other end at the branch
    update_partial_llhs(get_n_intervals()-1); 
}

void island_meth_changes::update_partial_llhs(unsigned n) {
// for a CpG island calculates contribution to partial likelihoods from a branch that is characterized in imc.
// pll are the partial likelihoods at an intermediate node, where the n-th rate change on the branch take place.
// That is, for n=1 the first IWE is updated.
// Note that it is assumed that the transition probabilities are up to date for the next node.
// Note also that, if the rates between nodes n-1 and n have changed, first the transition probabilities between
// nodes n and n+1 needs to be updated.

    if(n==0) { // this calculates the partial likelihoods for the parent node, assuming partial likelihoods
               // and also transition probabilities for all IWEs are up to date.
	pll[0]->set_scalefactor(pll[1]->get_scalefactor());
	bool rescale = true;
	for(unsigned s=0; s<pll[0]->nsites(); ++s) 
	    for(unsigned g=0; g<pll[0]->nratecat(); ++g) {
		vector<double> v(3,0);
		for(unsigned i=0; i<3; ++i) {	
		    for(unsigned j=0; j<3; ++j) {
			v[i] +=  (*transprob[0])[g][i][j] * (pll[1] -> get_part_L(s,g)[j]);
		    }
		    if(v[i]>1e-100) rescale=false;
		}
		pll[0] -> set_part_L(s,g,v);
	    }
	if(rescale) pll[0]->rescale();
	is_up_to_date();
    } else {
	if(n+1>=pll.size()) {
	    cerr << "too large n in void island_meth_changes::partial_llh(unsigned n)" << endl;
	    exit(1);
	}
	pll[n]->set_scalefactor(pll[n+1]->get_scalefactor());
	bool rescale = true;
	for(unsigned s=0; s<pll[n]->nsites(); ++s) 
	    for(unsigned g=0; g<pll[n]->nratecat(); ++g) {
		vector<double> v(3,0);
		for(unsigned i=0; i<3; ++i) {	
		    for(unsigned j=0; j<3; ++j) {
			v[i] +=  (*transprob[n])[g][i][j] * pll[n+1] -> get_part_L(s,g)[j];
		    }
		    if(v[i]>1e-100) rescale=false;
		}
		pll[n] -> set_part_L(s,g,v);
	    }
	if(rescale) pll[n]->rescale();
	update_partial_llhs(n-1);
    }
}

void island_meth_changes::accept_proposal(island_meth_changes * prop, bool copydaughter=false) {
    // The proposal *prop is accepted. pll[1] to pll[pll.size()-2] will then point to the 
    // previous pll tables on prop. The previous pll tables of *this are deleted and the
    // *prop will not have pointers its previous tables for global rate changes.
    // If copydaughter==true, the partial likelihoods at the daughter node are copied from *prop,
    // otherwise they are unchanged. In the latter case it is assumed (but not checked!) that they were the same anyway.
    change_time = prop->change_time;
    methyl_freq = prop->methyl_freq;
    daughter=pll.back(); // just to make sure...
    for(unsigned i=0; i+1 < pll.size(); ++i) {
	delete pll[i];
    }
    pll = prop->pll;
    pll.resize(prop->pll.size());
    for(unsigned i=0; i+1 < pll.size(); ++i) {
	pll[i] = prop->pll[i];
    }
    pll.back()=daughter;
    if(copydaughter) *daughter = *(prop->pll.back());
    is_up_to_date();
    prop->to_be_updated();
    while(prop->change_time.size() > 1) prop->minus_meth_change(1);
}

void island_meth_changes::copy_to_proposal(island_meth_changes * prop) const {
    // *prop will have the same strucure as *this but the partial likelihood
    // at the daughter node will not be changed. 
    // The partial likelihoods of *prop are not automatically updated.
    // If needed, this must be done with prop->partial_llh().
    // It is assumed that the pll tables that already exist in *prop have 
    // the same dimensions as those in *this.
    prop->change_time=change_time;
    prop->methyl_freq=methyl_freq;
    int sdiff = prop->pll.size() - pll.size();
    if(sdiff > 0) {
	unsigned n = pll.size()-1, m = prop->pll.size()-1;
	for(unsigned i=n; i < m; ++i) delete prop->pll[i];
	prop->pll.resize(pll.size());
	prop->pll.back()=prop->daughter;
    } else if(sdiff < 0) {
	unsigned n = prop->pll.size()-1;
	prop->pll.resize(pll.size());
	for(unsigned i=n; i+1 < pll.size(); ++i) 
	    prop->pll[i] = new island_partial_likelihoods_numscale(*pll[i]);
	prop->pll.back()=prop->daughter;
    }
    prop->to_be_updated();
}

// -------------- class island_partial_likelihoods_numscale ------------

//public:
void island_partial_likelihoods_numscale::rescale() {
    ++scalefactor;
    for(auto i=part_L.begin(); i!=part_L.end(); ++i) 
	for(auto j= i->begin(); j!=i->end(); ++j) 
	    for(auto k= j->begin(); k!=j->end(); ++k) 
		(*k)*=1e100;
}

void island_partial_likelihoods_numscale::rescale_if_necessary() {
    bool necessary=true;
    for(auto i=part_L.begin(); i!=part_L.end(); ++i) 
	for(auto j= i->begin(); j!=i->end(); ++j) 
	    for(auto k= j->begin(); k!=j->end(); ++k) 
		if ((*k) > 1e-100)
		    necessary=false;
    if(necessary) rescale();
}

// --------------------- (no class) -------------------------------------------

inline double addlog(double a, double b) {
    // returns log(exp(a)+exp(b)) without calculating exp(a) or exp(b) as they might be beyond double range
    return a<b ? b+log(1+exp(a-b)) : a+log(1+exp(b-a));
}

island_partial_logL island_branch_partial_logL(const island_partial_logL & pll, const island_meth_changes & imc) {
// for a CpG island calculates contribution to log partial likelihoods from a branch that is characterized in imc.
// pll are the partial likelihoods at the daughter node at the other end at the branch
    return island_branch_partial_logL(pll, imc, imc.get_n_intervals()-1); 
}

island_partial_logL island_branch_partial_logL(const island_partial_logL & pll, const island_meth_changes & imc, unsigned n) {
// for a CpG island calculates contribution to log partial likelihoods from a branch that is characterized in imc.
// pll are the partial likelihoods at an intermediate node, where the n-th rate change on the branch take place.
    
    vector<vector<vector<double>>> log_part_L(pll.nsites(),vector<vector<double>>(pll.nratecat(),vector<double>(3,1)));
    // this will be in the return object

    island_partial_logL ppll(pll.get_rate_acc());

    if(n==0) { // this is the direct intermediate node after the parent node
	vector<vector<vector<double>>> logtransprob(pll.nratecat(),vector<vector<double>>(3,vector<double>(3))); 
	// logtransprob[g][i][j] will be the log of the transition prob from state i to state j assuming rate category g
	vector<double> erate(pll.nratecat()), ler(pll.nratecat());
	const double t=imc.get_change_time(0);
	const vector<double> rate_acc(pll.get_rate_acc());
	for(unsigned i=0; i<pll.nratecat();++i) {
	    erate[i] = exp(-rate_acc[i]*t);
	    ler[i] = log(1-erate[i]);
	}  // FOR NOW ONLY F81 TYPE OF MODEL IS IMPLEMENTED. NEED FUNCTION FOR MATRIX EXPONENTIAL TO GENERALIZE. MAYBE USE ARMADILLO LIBRARY
	vector<double> logfreq(3); // to be the log frequencies of the states
	for(unsigned i=0; i<3; ++i) logfreq[i]=log(imc.get_methyl_freq(0)[i]); 
	for(unsigned g=0; g<pll.nratecat(); ++g)
	    for(unsigned i=0; i<3; ++i) {
		for(unsigned j=0; j<3; ++j) {
		    logtransprob[g][i][j] = ler[g]+logfreq[j];
		}
		logtransprob[g][i][i] = log(exp(logtransprob[g][i][i])+erate[g]);		
	    } 
	for(unsigned s=0; s<pll.nsites(); ++s) 
	    for(unsigned g=0; g<pll.nratecat(); ++g) {
		vector<double> lfpl(pll.get_log_part_L(s,g));
		for(unsigned i=0; i<3; ++i) {
		    log_part_L[s][g][i] = logtransprob[g][i][0] + lfpl[0];
		    for(unsigned j=1; j<3; ++j) {
			log_part_L[s][g][i] = addlog(log_part_L[s][g][i], logtransprob[g][i][j] + lfpl[j]);
		    }
		}
	    }
	ppll.set_log_part_L(log_part_L);
    } else {
	vector<vector<vector<double>>> logtransprob(pll.nratecat(),vector<vector<double>>(3,vector<double>(3))); 
	// logtransprob[g][i][j] will be the log of the transition prob from state i to state j assuming rate category g
	vector<double> erate(pll.nratecat()), ler(pll.nratecat());
	const double t=imc.get_change_time(n)-imc.get_change_time(n-1);
	const vector<double> rate_acc(pll.get_rate_acc());
	for(unsigned i=0; i<pll.nratecat();++i) {
	    erate[i] = exp(-rate_acc[i]*t);
	    ler[i] = log(1-erate[i]);
	}  // FOR NOW ONLY F81 TYPE OF MODEL IS IMPLEMENTED. NEED FUNCTION FOR MATRIX EXPONENTIAL TO GENERALIZE. MAYBE USE ARMADILLO LIBRARY.
	vector<double> logfreq(3); // to be the log frequencies of the states
	for(unsigned i=0; i<3; ++i) logfreq[i]=log(imc.get_methyl_freq(n)[i]); 
	vector<double> freqdiff(3); // changes in state frequencies
	for(unsigned i=0; i<3; ++i) freqdiff[i]=imc.get_methyl_freq(n)[i]-imc.get_methyl_freq(n-1)[i]; 
	for(unsigned g=0; g<pll.nratecat(); ++g)
	    for(unsigned i=0; i<3; ++i) {
		for(unsigned j=0; j<3; ++j) {
		    logtransprob[g][i][j] = ler[g]+logfreq[j];
		    if(freqdiff[i]<0 && freqdiff[j]>0) {
			// note here that if the probability of state i becomes smaller at an island-wide methylation
			// change event, each site in state i has a probability of freqdiff[i]/f_i to change to one of
			// the other states, where f_i := get_methyl_freq(n-1)[i] is here the frequency of i state sites.
			// If both other states increase their probability, the sum of their differences
			// must be -freqdiff[i], and the probability that a changed state becomes j is then the fraction 
			// freqdiff[j]/(-freqdiff[i]). Thus, -freqdiff[i] cancels and the probability that a site that
			// was in i becomes j is just freqdiff[j]/f_i.
			// If, however, only the probability of j increases, all changed sites that were i before become j,
                        // and the probability to go from i to j in this way is just -freqdiff[i]/f_i.
                        // All this can be hidden by another single-site methylation change, so the increase
			// in the probability to change from i to j contains the factor erate[g], which is the probability
			// that no single-site methylation change affects the end of the next branch segment.
			double f=erate[g]*( freqdiff[j] < -freqdiff[i] ? freqdiff[j] : -freqdiff[i] )/imc.get_methyl_freq(n-1)[i];
			logtransprob[g][i][j] = addlog(logtransprob[g][i][j],log(f));
		    }
		}
		double f=erate[g];
		if(freqdiff[i]<0) f *= (1+freqdiff[i]/imc.get_methyl_freq(n-1)[i]);
		logtransprob[g][i][i] = log(exp(logtransprob[g][i][i])+f);
	    } 
	for(unsigned s=0; s<pll.nsites(); ++s) 
	    for(unsigned g=0; g<pll.nratecat(); ++g) {
		vector<double> lfpl(pll.get_log_part_L(s,g));
		for(unsigned i=0; i<3; ++i) {	
		    log_part_L[s][g][i] = logtransprob[g][i][0] + lfpl[0];
		    for(unsigned j=1; j<3; ++j) {
			log_part_L[s][g][i] = addlog(log_part_L[s][g][i], logtransprob[g][i][j] + lfpl[j]);
		    }
		}
	    }
	ppll.set_log_part_L(log_part_L);
	ppll = island_branch_partial_logL(ppll ,imc, n-1);
    }
    return ppll;
}

