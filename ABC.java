package binMeta;

import java.io.*;
public class ABC extends binMeta {

	/* Control Parameters of ABC algorithm */
	int NP; /* The number of colony size (employed bees+onlooker bees) */
	int FoodNumber; /* The number of food sources equals the half of the colony size */
	int limit; /* A food source which could not be improved through "limit" trials is abandoned
				 * by its employed bee */
	int maxCycle; /* The number of cycles for foraging {a stopping criteria} */

	/* Problem specific variables */
	int D; /* The number of parameters of the problem to be optimized */
	double Lbvec[]; /* lower bound of the parameters. */
	double Ubvec[]; /* upper bound of the parameters. */
	
	double f[]=new double[FoodNumber];        /*f is a vector holding objective function values associated with food sources */
	double fitness[]=new double[FoodNumber];      /*fitness is a vector holding fitness (quality) values associated with food sources*/
	double trial[]=new double[FoodNumber];         /*trial is a vector holding trial numbers through which solutions can not be improved*/
	double prob[]=new double[FoodNumber];          /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
	
	double GlobalMin;        /*Optimum solution obtained by ABC algorithm*/
	double GlobalParams[]=new double[D];              /*Parameters of the optimum solution*/
	double r; /*a random number in the range [0,1)*/

	// ABC constructor
	public ABC(Data foodSources, Objective obj, long maxTime, int iNP, int ilimit, int imaxCycle, double[] iLbvec,double[] iUbvec) {
		try {
			  this.Lbvec=iLbvec;
			  this.Ubvec=iUbvec;
			this.NP = iNP;
			this.FoodNumber = NP / 2;
			this.limit = ilimit;
			this.maxCycle = imaxCycle;
			String msg = "Impossible to create ABC object: ";
			if (maxTime <= 0)
				throw new Exception(msg + "the maximum execution time is 0 or even negative");
			this.maxTime = maxTime;
			if (foodSources == null)
				throw new Exception(msg + "the reference to the starting point is null");
			this.solution = foodSources;
			if (obj == null)
				throw new Exception(msg + "the reference to the objective is null");
			this.obj = obj;
			this.objValue = this.obj.value(this.solution);
			this.metaName = "ABC";
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	  //Initializing food resources,fitness values, global optimum parameters 
    void initialize() {
    	
    	Data Da = new Data(this.solution);
    }

	/* Fitness function */
	double CalculateFitness(double fun) {
		double result = 0;
		if (fun >= 0) {
			result = 1 / (fun + 1);
		} else {

			result = 1 + Math.abs(fun);
		}
		return result;
	}
	
	//Boundary correction
		 double[] simpleboundss(double s[])
		 {
		   for(int i=0;i<D;i++)
		   {if(s[i]<Lbvec[i])
		    {s[i]=Lbvec[i];}
		    if(s[i]>Ubvec[i])
		    {s[i]=Ubvec[i];}
		   }	 
		   return s;	 
		 }

	void SendEmployedBees() {
		int i, j;
		/* Employed Bee Phase */
		for (i = 0; i < FoodNumber; i++) {
			/* The parameter to be changed is determined randomly */
			r = ((double) Math.random() * 32767 / ((double) (32767) + (double) (1)));
			param2change = (int) (r * D);

			/*
			 * A randomly chosen solution is used in producing a mutant solution of the
			 * solution i
			 */
			r = ((double) Math.random() * 32767 / ((double) (32767) + (double) (1)));
			neighbour = (int) (r * FoodNumber);

			/* Randomly selected solution must be different from the solution i */
			// while(neighbour==i)
			// {
			// r = ( (double)Math.random()*32767 / ((double)(32767)+(double)(1)) );
			// neighbour=(int)(r*FoodNumber);
			// }
			for (j = 0; j < D; j++)
				solutions[j] = Foods[i][j];

			/* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
			r = ((double) Math.random() * 32767 / ((double) (32767) + (double) (1)));
			solutions[param2change] = Foods[i][param2change]
					+ (Foods[i][param2change] - Foods[neighbour][param2change]) * (r - 0.5) * 2;

			/*
			 * if generated parameter value is out of boundaries, it is shifted onto the
			 * boundaries
			 */
			
			simplebounds(solutions);
			
	 		//ObjValSol = calculateFunction(solution);
			FitnessSol = CalculateFitness(ObjValSol);

			/*
			 * a greedy selection is applied between the current solution i and its mutant
			 */
			if (FitnessSol > fitness[i]) {

				/*
				 * If the mutant solution is better than the current solution i, replace the
				 * solution with the mutant and reset the trial counter of solution i
				 */
				trial[i] = 0;
				for (j = 0; j < D; j++)
					Foods[i][j] = solutions[j];
				f[i] = ObjValSol;
				fitness[i] = FitnessSol;
			} else { /* if the solution i can not be improved, increase its trial counter */
				trial[i] = trial[i] + 1;
			}

		}

		/* end of employed bee phase */

	}

	@Override
	public void optimize() // by ABC
	{

	}
}
