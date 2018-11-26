
#include <iostream>
#include <vector>
#include "info_classes.hh"



bool ancestry_compare(clone_info a,clone_info b,int k)
{
  for (int i=0; i<k;++i)
    if (a.mutation[i]!=b.mutation[i]||
	a.gene[i]!=b.gene[i]||
	a.fitness_effect_f[i]!=b.fitness_effect_h[i]||
	a.fitness_effect_h[i]!=b.fitness_effect_h[i])
      return false;
  return true;
}

	  
clone_info::clone_info()
{

  total_size=0;

  
}
