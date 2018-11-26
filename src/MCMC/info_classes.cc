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
