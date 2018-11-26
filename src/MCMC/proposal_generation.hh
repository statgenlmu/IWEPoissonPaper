
   /* This file is part of methylation_simulator.

    methylation_simulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    methylation_simulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with methylation_simulator.  If not, see <https://www.gnu.org/licenses/>. */


/*
 * The following struct contains information about the state the Markov chain was
 * in before the jump, or may be in after a jump.
 * Since always only one node of the tree is modified, we just store the
 * address to the tree, the likelyhood and the parameters that can be subject to change.
 */
struct Markov_State_Info
{
  double likelyhood;
  back_tree * address_of_node;

  double p_change_rate;
  double gw_alpha;
  double gw_beta;

  Markov_State_Info(){};
  Markov_State_Info(back_tree * node, double lklhood): address_of_node(node), likelyhood(lklyhood)
  {
    gw_alpha=node->gw_alpha;
    gw_beta=node->gw_beta;
    p_change_rate= (node->cpg_container->get_c_rate());
  }
  /*
   * The following function accesses the node stored in address_of_node
   * of the tree and replaces p_change_rate and gw_rates in the tree
   */
  void insert_itself_tree()
  {
    node->set_gw_alpha(gw_alpha);
    node->set_gw_beta(gw_beta);
    node->set_c_rate(p_change_rate);
  }
  
};


/*
 * The following function generates a proposal for
 * a new shape parameter from an old one by adding a normally
 * distributed random variable with mean zero and sd 0.5
 */

double propose_new_shape_parameter(double old_shape, mt19937 & re)
{
  normal_distribution<double> proposal_dist(0, 0.5);
  double new_shape=old_shape+proposal_dist(re);
  while (new_shape <0) new_shape=old_shape+proposal_dist(re);
  return new_shape;
}

/*
 * The following function takes a Markov state and generates information to
 * implement a proposal.
 */

//TODO:Clarify proposal distributions with Dirk.
Markov_State_Info generate_proposal_state(Markov_State_Info & old_state, mt19937 & re)
{
  // Copy node address and likelyhood
  Markov_State_Info new_state=old_state;
  //decide wheter you want to modify change rate or distribution
  Bernoulli_distribution mod(0.5);
  //produce the proposal
  if (mod(re))
    {
      normal_distribution<double> proposal_dist(0, 0.5);
      new_state.p_change_rate=exp(log(new_state.p_change_rate)+proposal_dist(re));
      return new_state;
    }
  new_state.gw_alpha=propose_new_shape_parameter(new_state.gw_alpha, re);
  new_state.gw_beta=propose_new_shape_parameter(new_state.gw_beta, re);
  return new_state;  
}

/*
 * The following function checks wheter a given proposal should be accepted or
 * rejected based on likelyhood ratios and proposal_probability ratios.
 */

bool accept_proposal(double probability_ratios, double likelyhood_ratios, mt19937 & re);
{
  Bernoulli_distribution acceptance(min(1,probability_ratios*likely_hood_ratios));
  return acceptance(re);
}

void MCMC_step(vector<back_tree *> & node_container, const int n_quantiles, mt19937 & re)
  {
    //find_likelyhood
    double likelyhood=node_container[0]->combined_likelyhood_of_tree_with_known_nodes_using_quantiles(n_quantiles);
    //find node to change
    uniform_int_distribution next_node(0, node_container.size());
    back_tree * node=next_node(re);
    //store information about the current state
    Markov_State_Info old_state(node, likelyhood);
    //generate proposal state
    Markov_State_Info new_state=generate_proposal(old_state, re);
    //insert proposal in tree
    new_state.insert_itself_tree();
    //calculate_proposal_probability_ratios
    double proposal_probability_ratio=1; //One at the moment because proposals are symmetric. TODO:Clarify proposal distributions with Dirk
    //find likelyhood of new tree
    double new_likely_hood=node_container[0]->combined_likelyhood_of_tree_with_known_nodes_using_quantiles(n_quantiles);
    //accept or reject if reject reverttree
    if (!accept_proposal( proposal_probability_ratio, new_likely_hood/likelyhood, re))
      old_state.insert_itself_tree();
  }
      
