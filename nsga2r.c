/* NSGA-II routine (implementation of the 'main' function) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>

# include "global.h"
# include "rand.h"
#include "awi.h"

# include <iostream>
# include <fstream>
#include <sstream>

# include <limits>

int nreal;
int nbin;
int nobj;
int ncon;
int popsize;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;
double *min_binvar;
double *max_binvar;
int bitlength;
int choice;
int obj1;
int obj2;
int obj3;
int angle1;
int angle2;

double gbest_x = 0.0, gbest_y = 0.0;
int ggen_idx = 1, gind_idx = 0, gbest_gen = 1, gbest_ind = 0;
size_t gbest_i = 0;
std::ofstream gsummary_file("summary.txt", std::ios::out);

void hack();

int gport = 6665;

void write_data(population * pop, int gen_idx, std::ostream & os) {
  for(size_t i = 0; i < popsize; ++i) {
    os << "# " << gen_idx << "," << i << "\n";
    for(size_t r = 0; r < pop->ind[i].data.size(); ++r) {
      os << r << "\t";
      for(size_t c = 0; c < pop->ind[i].data[r].size(); ++c) {
	os << pop->ind[i].data[r][c] << "\t";
      }
      os << "\n";
    }
    os << "\n";
  }
  os << "\n";
}


void write_metrics(population * pop, int gen_idx, std::ostream & os) {

  os << "# " << gen_idx << "\n";
  for(size_t i = 0; i < popsize; ++i) {
    os << i << "\t" 
       << pop->ind[i].hit << "\t"
       << pop->ind[i].reached << "\t"
       << pop->ind[i].proximity << "\t"
       << pop->ind[i].lifetime << "\t"
       << pop->ind[i].dist << "\t"
       << pop->ind[i].lc << "\t"
       << pop->ind[i].fc << "\t"
       << pop->ind[i].rc << "\n";
  }
  os << "\n";
}
      

/* void write_feasible_data(population * pop, int gen_idx, std::ostream & os) { */
/*   for(size_t i = 0; i < popsize; ++i) { */
/*     if(pop->ind[i].constr_violation == 0.0) { */
/*       os << "# " << gen_idx << "," << i << "\n"; */
/*       for(size_t r = 0; r < pop->ind[i].data.size(); ++r) { */
/* 	os << r << "\t"; */
/* 	for(size_t c = 0; c < pop->ind[i].data[r].size(); ++c) { */
/* 	  os << pop->ind[i].data[r][c] << "\t"; */
/* 	} */
/* 	os << "\n"; */
/*       } */
/*       os << "\n"; */
/*     } */
/*   } */
/*   os << "\n"; */
/* } */

void write_feasible_data(population * pop, int gen_idx, std::ostream & os) {
  for(size_t i = 0; i < popsize; ++i) {
    os << "# " << gen_idx << "," << i << "\n";
    if(pop->ind[i].constr_violation == 0.0) {
      for(size_t r = 0; r < pop->ind[i].data.size(); ++r) {
	os << r << "\t";
	for(size_t c = 0; c < pop->ind[i].data[r].size(); ++c) {
	  os << pop->ind[i].data[r][c] << "\t";
	}
	os << "\n";
      }
    } else {
      os << 0.0 << "\t?\n";
    }
    os << "\n";
  }
  os << "\n";
}


int main (int argc, char **argv)
{
  int i;
  FILE *fpt1;
  FILE *fpt2;
  FILE *fpt3;
  FILE *fpt4;
  FILE *fpt5;
  FILE *gp;
  population *parent_pop;
  population *child_pop;
  population *mixed_pop;

  std::ios::fmtflags fmtflags = std::cout.flags();
  std::cout.setf(std::ios::showpos | std::ios::showpoint);
  
  //std::cout.precision(6);
  //std::cout.precision(std::numeric_limits<double>::digits10);

  std::ofstream robo_file("r.txt", std::ios::out);
  robo_file.setf(std::ios::showpos | std::ios::showpoint);

  std::ofstream metrics_file("m.txt", std::ios::out);
  metrics_file.setf(std::ios::showpos | std::ios::showpoint);

  std::ofstream feas_file("f.txt", std::ios::out);
  std::ofstream nhits_file("nhits.txt", std::ios::out);
  std::ofstream nreaches_file("nreaches.txt", std::ios::out);
  size_t nreaches = 0;
  size_t nhits = 0;

  if (argc<2)
    {
      printf("\n Usage ./nsga2r random_seed \n");
      fflush(stdout);
      exit(1);
    }
  seed = (double)atof(argv[1]);
  if (seed<=0.0 || seed>=1.0)
    {
      printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
      fflush(stdout);
      exit(1);
    }

  {// awhan  
    gport = atoi(argv[2]);
  }
  

  fpt1 = fopen("initial_pop.out","w");
  fpt2 = fopen("final_pop.out","w");
  fpt3 = fopen("best_pop.out","w");
  fpt4 = fopen("all_pop.out","w");
  fpt5 = fopen("params.out","w");
  fprintf(fpt1,"# This file contains the data of initial population\n");
  fprintf(fpt2,"# This file contains the data of final population\n");
  fprintf(fpt3,"# This file contains the data of final feasible population (if found)\n");
  fprintf(fpt4,"# This file contains the data of all generations\n");
  fprintf(fpt5,"# This file contains information about inputs as read by the program\n");
#ifdef VERBOSE
  printf("\n Enter the problem relevant and algorithm relevant parameters ... ");
  printf("\n Enter the population size (a multiple of 4) : ");
  fflush(stdout);
#endif
  scanf("%d",&popsize);
  if (popsize<4 || (popsize%4)!= 0)
    {
      printf("\n population size read is : %d",popsize);
      printf("\n Wrong population size entered, hence exiting \n");
      fflush(stdout);
      exit (1);
    }
#ifdef VERBOSE
  printf("\n Enter the number of generations : ");
#endif
  scanf("%d",&ngen);
  if (ngen<1)
    {
      printf("\n number of generations read is : %d",ngen);
      printf("\n Wrong nuber of generations entered, hence exiting \n");
      fflush(stdout);
      exit (1);
    }
#ifdef VERBOSE
  printf("\n Enter the number of objectives : ");
  fflush(stdout);
#endif
  scanf("%d",&nobj);
  if (nobj<1)
    {
      printf("\n number of objectives entered is : %d",nobj);
      printf("\n Wrong number of objectives entered, hence exiting \n");
      exit (1);
    }
#ifdef VERBOSE
  printf("\n Enter the number of constraints : ");
  fflush(stdout);
#endif
  scanf("%d",&ncon);
  if (ncon<0)
    {
      printf("\n number of constraints entered is : %d",ncon);
      printf("\n Wrong number of constraints enetered, hence exiting \n");
      exit (1);
    }
#ifdef VERBOSE
  printf("\n Enter the number of real variables : ");
  fflush(stdout);
#endif
  scanf("%d",&nreal);
  if (nreal<0)
    {
      printf("\n number of real variables entered is : %d",nreal);
      printf("\n Wrong number of variables entered, hence exiting \n");
      exit (1);
    }
  if (nreal != 0)
    {
      min_realvar = (double *)malloc(nreal*sizeof(double));
      max_realvar = (double *)malloc(nreal*sizeof(double));
      for (i=0; i<nreal; i++)
        {
#ifdef VERBOSE
	  printf ("\n Enter the lower limit of real variable %d : ",i+1);
#endif
	  scanf ("%lf",&min_realvar[i]);
#ifdef VERBOSE
	  printf ("\n Enter the upper limit of real variable %d : ",i+1);
#endif
	  scanf ("%lf",&max_realvar[i]);
	  if (max_realvar[i] <= min_realvar[i])
            {
	      printf("\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n");
	      exit(1);
            }
        }
#ifdef VERBOSE
      printf ("\n Enter the probability of crossover of real variable (0.6-1.0) : ");
#endif
      scanf ("%lf",&pcross_real);
      if (pcross_real<0.0 || pcross_real>1.0)
        {
	  printf("\n Probability of crossover entered is : %e",pcross_real);
	  printf("\n Entered value of probability of crossover of real variables is out of bounds, hence exiting \n");
	  exit (1);
        }
#ifdef VERBOSE
      printf ("\n Enter the probablity of mutation of real variables (1/nreal) : ");
#endif
      scanf ("%lf",&pmut_real);
      if (pmut_real<0.0 || pmut_real>1.0)
        {
	  printf("\n Probability of mutation entered is : %e",pmut_real);
	  printf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
	  exit (1);
        }
#ifdef VERBOSE
      printf ("\n Enter the value of distribution index for crossover (5-20): ");
#endif
      scanf ("%lf",&eta_c);
      if (eta_c<=0)
        {
	  printf("\n The value entered is : %e",eta_c);
	  printf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
	  exit (1);
        }
#ifdef VERBOSE
      printf ("\n Enter the value of distribution index for mutation (5-50): ");
#endif
      scanf ("%lf",&eta_m);
      if (eta_m<=0)
        {
	  printf("\n The value entered is : %e",eta_m);
	  printf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
	  exit (1);
        }
    }
#ifdef VERBOSE
  printf("\n Enter the number of binary variables : ");
#endif
  scanf("%d",&nbin);
  if (nbin<0)
    {
      printf ("\n number of binary variables entered is : %d",nbin);
      printf ("\n Wrong number of binary variables entered, hence exiting \n");
      exit(1);
    }
  if (nbin != 0)
    {
      nbits = (int *)malloc(nbin*sizeof(int));
      min_binvar = (double *)malloc(nbin*sizeof(double));
      max_binvar = (double *)malloc(nbin*sizeof(double));
      for (i=0; i<nbin; i++)
        {
	  printf ("\n Enter the number of bits for binary variable %d : ",i+1);
	  scanf ("%d",&nbits[i]);
	  if (nbits[i] < 1)
            {
	      printf("\n Wrong number of bits for binary variable entered, hence exiting");
	      exit(1);
            }
	  printf ("\n Enter the lower limit of binary variable %d : ",i+1);
	  scanf ("%lf",&min_binvar[i]);
	  printf ("\n Enter the upper limit of binary variable %d : ",i+1);
	  scanf ("%lf",&max_binvar[i]);
	  if (max_binvar[i] <= min_binvar[i])
            {
	      printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
	      exit(1);
            }
        }
#ifdef VERBOSE
      printf ("\n Enter the probability of crossover of binary variable (0.6-1.0): ");
#endif
      scanf ("%lf",&pcross_bin);
      if (pcross_bin<0.0 || pcross_bin>1.0)
        {
	  printf("\n Probability of crossover entered is : %e",pcross_bin);
	  printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
	  exit (1);
        }
      printf ("\n Enter the probability of mutation of binary variables (1/nbits): ");
      scanf ("%lf",&pmut_bin);
      if (pmut_bin<0.0 || pmut_bin>1.0)
        {
	  printf("\n Probability of mutation entered is : %e",pmut_bin);
	  printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
	  exit (1);
        }
    }
  if (nreal==0 && nbin==0)
    {
      printf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
      exit(1);
    }
  choice=0;
#ifdef VERBOSE
  printf("\n Do you want to use gnuplot to display the results realtime (0 for NO) (1 for yes) : ");
#endif
  scanf("%d",&choice);
  if (choice!=0 && choice!=1)
    {
      printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
      exit(1);
    }
  if (choice==1)
    {
      gp = popen(GNUPLOT_COMMAND,"w");
      if (gp==NULL)
        {
	  printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
	  printf("\n Edit the string to suit your system configuration and rerun the program\n");
	  exit(1);
        }
      if (nobj==2)
        {
	  printf("\n Enter the objective for X axis display : ");
	  scanf("%d",&obj1);
	  if (obj1<1 || obj1>nobj)
            {
	      printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
	      exit(1);
            }
	  printf("\n Enter the objective for Y axis display : ");
	  scanf("%d",&obj2);
	  if (obj2<1 || obj2>nobj)
            {
	      printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
	      exit(1);
            }
	  obj3 = -1;
        }
      else
        {
	  printf("\n #obj > 2, 2D display or a 3D display ?, enter 2 for 2D and 3 for 3D :");
	  scanf("%d",&choice);
	  if (choice!=2 && choice!=3)
            {
	      printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
	      exit(1);
            }
	  if (choice==2)
            {
	      printf("\n Enter the objective for X axis display : ");
	      scanf("%d",&obj1);
	      if (obj1<1 || obj1>nobj)
                {
		  printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
		  exit(1);
                }
	      printf("\n Enter the objective for Y axis display : ");
	      scanf("%d",&obj2);
	      if (obj2<1 || obj2>nobj)
                {
		  printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
		  exit(1);
                }
	      obj3 = -1;
            }
	  else
            {
	      printf("\n Enter the objective for X axis display : ");
	      scanf("%d",&obj1);
	      if (obj1<1 || obj1>nobj)
                {
		  printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
		  exit(1);
                }
	      printf("\n Enter the objective for Y axis display : ");
	      scanf("%d",&obj2);
	      if (obj2<1 || obj2>nobj)
                {
		  printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
		  exit(1);
                }
	      printf("\n Enter the objective for Z axis display : ");
	      scanf("%d",&obj3);
	      if (obj3<1 || obj3>nobj)
                {
		  printf("\n Wrong value of Z objective entered, value entered was %d\n",obj3);
		  exit(1);
                }
	      printf("\n You have chosen 3D display, hence location of eye required \n");
	      printf("\n Enter the first angle (an integer in the range 0-180) (if not known, enter 60) :");
	      scanf("%d",&angle1);
	      if (angle1<0 || angle1>180)
                {
		  printf("\n Wrong value for first angle entered, hence exiting \n");
		  exit(1);
                }
	      printf("\n Enter the second angle (an integer in the range 0-360) (if not known, enter 30) :");
	      scanf("%d",&angle2);
	      if (angle2<0 || angle2>360)
                {
		  printf("\n Wrong value for second angle entered, hence exiting \n");
		  exit(1);
                }
            }
        }
    }
#ifdef VERBOSE
  printf("\n Input data successfully entered, now performing initialization \n");
#endif

  fprintf(fpt5,"\n Population size = %d",popsize);
  fprintf(fpt5,"\n Number of generations = %d",ngen);
  fprintf(fpt5,"\n Number of objective functions = %d",nobj);
  fprintf(fpt5,"\n Number of constraints = %d",ncon);
  fprintf(fpt5,"\n Number of real variables = %d",nreal);
  if (nreal!=0)
    {
      for (i=0; i<nreal; i++)
        {
	  fprintf(fpt5,"\n Lower limit of real variable %d = %e",i+1,min_realvar[i]);
	  fprintf(fpt5,"\n Upper limit of real variable %d = %e",i+1,max_realvar[i]);
        }
      fprintf(fpt5,"\n Probability of crossover of real variable = %e",pcross_real);
      fprintf(fpt5,"\n Probability of mutation of real variable = %e",pmut_real);
      fprintf(fpt5,"\n Distribution index for crossover = %e",eta_c);
      fprintf(fpt5,"\n Distribution index for mutation = %e",eta_m);
    }
  fprintf(fpt5,"\n Number of binary variables = %d",nbin);
  if (nbin!=0)
    {
      for (i=0; i<nbin; i++)
        {
	  fprintf(fpt5,"\n Number of bits for binary variable %d = %d",i+1,nbits[i]);
	  fprintf(fpt5,"\n Lower limit of binary variable %d = %e",i+1,min_binvar[i]);
	  fprintf(fpt5,"\n Upper limit of binary variable %d = %e",i+1,max_binvar[i]);
        }
      fprintf(fpt5,"\n Probability of crossover of binary variable = %e",pcross_bin);
      fprintf(fpt5,"\n Probability of mutation of binary variable = %e",pmut_bin);
    }
  fprintf(fpt5,"\n Seed for random number generator = %e",seed);
  bitlength = 0;
  if (nbin!=0)
    {
      for (i=0; i<nbin; i++)
        {
	  bitlength += nbits[i];
        }
    }
  fprintf(fpt1,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
  fprintf(fpt2,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
  fprintf(fpt3,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
  fprintf(fpt4,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
  nbinmut = 0;
  nrealmut = 0;
  nbincross = 0;
  nrealcross = 0;
  parent_pop = (population *)malloc(sizeof(population));
  child_pop = (population *)malloc(sizeof(population));
  mixed_pop = (population *)malloc(sizeof(population));
  allocate_memory_pop (parent_pop, popsize);
  allocate_memory_pop (child_pop, popsize);
  allocate_memory_pop (mixed_pop, 2*popsize);
  randomize();
  initialize_pop (parent_pop);
#ifdef VERBOSE
  printf("\n Initialization done, now performing first generation\n");
#endif
  decode_pop(parent_pop);

  {// awhan
    std::stringstream ss;
    ss << gport;
    std::string command =  "player -p " + ss.str() + " simple.cfg > /dev/null 2>&1 &";
    system(command.c_str());
    system("sleep 5");
    //system("sleep 10"); // extra time to enable trails
    //hack();
  }

  evaluate_pop (parent_pop);

  assign_rank_and_crowding_distance (parent_pop);
  report_pop (parent_pop, fpt1);
  fprintf(fpt4,"# gen = 1\n");
  report_pop(parent_pop,fpt4);

  {// awi
    write_data(parent_pop, 1, robo_file);
    write_feasible_data(parent_pop, 1, feas_file);
    write_metrics(parent_pop, 1, metrics_file);
    nreaches = num_reaches(parent_pop);
    nreaches_file << 1 << "\t" << nreaches << "\n";
    nhits = num_hits(parent_pop);
    nhits_file << 1 << "\t" << nhits << "\n";
  }

  printf("\n ------------------------------gen = 1");
  std::cout << " num feas = " << nreaches
	    << " nhits = " << nhits << std::endl;

  fflush(stdout);
  if (choice!=0)    onthefly_display (parent_pop,gp,1);
  fflush(fpt1);
  fflush(fpt2);
  fflush(fpt3);
  fflush(fpt4);
  fflush(fpt5);
  sleep(1);
  for (i=2; i<=ngen; i++){
    ggen_idx = i;
    selection (parent_pop, child_pop);
    mutation_pop (child_pop);
    decode_pop(child_pop);

    evaluate_pop(child_pop);
    merge (parent_pop, child_pop, mixed_pop);
    fill_nondominated_sort (mixed_pop, parent_pop);
    //show_data(parent_pop);
    /* Comment following four lines if information for all
       generations is not desired, it will speed up the execution */
    fprintf(fpt4,"# gen = %d\n",i);
    report_pop(parent_pop,fpt4);
    write_data(parent_pop, i, robo_file);// awi
    write_feasible_data(parent_pop, i, feas_file);// awi
    write_metrics(parent_pop, i, metrics_file);
    fflush(fpt4);
    //if (choice!=0)    onthefly_display (parent_pop,gp,i);
    printf("\n ------------------------------gen = %d",i);
    nreaches = num_reaches(parent_pop);
    nreaches_file << i << "\t" << nreaches << "\n";
    nhits = num_hits(parent_pop);
    nhits_file << i << "\t" << nhits << "\n";
    std::cout << " num feas = " << nreaches
	      << " nhits = " << nhits << std::endl;
  }
  printf("\n Generations finished, now reporting solutions");
  report_pop(parent_pop,fpt2);
  report_feasible(parent_pop,fpt3);
  if (nreal!=0)
    {
      fprintf(fpt5,"\n Number of crossover of real variable = %d",nrealcross);
      fprintf(fpt5,"\n Number of mutation of real variable = %d",nrealmut);
    }
  if (nbin!=0)
    {
      fprintf(fpt5,"\n Number of crossover of binary variable = %d",nbincross);
      fprintf(fpt5,"\n Number of mutation of binary variable = %d",nbinmut);
    }
  fflush(stdout);
  fflush(fpt1);
  fflush(fpt2);
  fflush(fpt3);
  fflush(fpt4);
  fflush(fpt5);
  fclose(fpt1);
  fclose(fpt2);
  fclose(fpt3);
  fclose(fpt4);
  fclose(fpt5);
  if (choice!=0)
    {
      pclose(gp);
    }
  if (nreal!=0)
    {
      free (min_realvar);
      free (max_realvar);
    }
  if (nbin!=0)
    {
      free (min_binvar);
      free (max_binvar);
      free (nbits);
    }
  deallocate_memory_pop (parent_pop, popsize);
  deallocate_memory_pop (child_pop, popsize);
  deallocate_memory_pop (mixed_pop, 2*popsize);
  free (parent_pop);
  free (child_pop);
  free (mixed_pop);
  printf("\n Routine successfully exited \n");

  {// awhan
    
    std::cout << "gbest_x = " << gbest_x << std::endl;
    std::cout << "gbest_y = " << gbest_y << std::endl;
    std::cout << "gbest_gen = " << gbest_gen << std::endl;
    std::cout << "gbest_ind = " << gbest_ind << std::endl;


    gsummary_file << "gbest_x = " << gbest_x << std::endl;
    gsummary_file << "gbest_y = " << gbest_y << std::endl;
    gsummary_file << "gbest_gen = " << gbest_gen << std::endl;
    gsummary_file << "gbest_ind = " << gbest_ind << std::endl;
  }

  std::cout.flags(fmtflags);
  return 0;
}
