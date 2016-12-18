/* Routine for evaluating population members  */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#include <iostream> // awi

# include "global.h"
# include "rand.h"

extern int gind_idx; // awi

/* Routine to evaluate objective function values and constraints for a population */
void evaluate_pop (population *pop)
{
  int i;
  for (i=0; i<popsize; i++)
    {
      gind_idx = i; // awi
      evaluate_ind (&(pop->ind[i]));
      //std::cout << "constr_violation :: " << pop->ind[i].constr_violation << std::endl;
    }
  return;
}

/* Routine to evaluate objective function values and constraints for an individual */
void evaluate_ind (individual *ind)
{
  int j;
  //test_problem (ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
  test_problem (ind->xreal, ind->xbin, ind->gene, ind->obj, 
		ind->constr, ind->data, ind->hit, ind->reached,
		ind->proximity, ind->lifetime, ind->dist,
		ind->lc, ind->fc, ind->rc); // awi
  if (ncon==0)
    {
      ind->constr_violation = 0.0;
    }
  else
    {
      ind->constr_violation = 0.0;
      for (j=0; j<ncon; j++)
        {
	  if (ind->constr[j]<0.0)
            {
	      ind->constr_violation += ind->constr[j];
            }
        }
    }
  return;
}
