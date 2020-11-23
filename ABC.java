package binMeta;

import java.io.*;

public class ABC extends binMeta{
	
	/* Control Parameters of ABC algorithm*/
	int NP; /* The number of colony size (employed bees+onlooker bees)*/
	int FoodNumber; /*The number of food sources equals the half of the colony size*/
	int limit;  /*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
	int maxCycle; /*The number of cycles for foraging {a stopping criteria}*/

	/* Problem specific variables*/
	int D; /*The number of parameters of the problem to be optimized*/
	double lb; /*lower bound of the parameters. */
	double ub; /*upper bound of the parameters. */
	
	  // ABC constructor
	   public ABC(Data foodSources, Objective obj, long maxTime, int iNP, int ilimit, int imaxCycle, double iLb, double iUb)
	   {
	      try
	      {
	    	 this.NP = iNP;
	  		 this.FoodNumber = NP / 2;
	  		 this.limit = ilimit;
	  		 this.maxCycle = imaxCycle;
	         String msg = "Impossible to create ABC object: ";
	         if (maxTime <= 0) throw new Exception(msg + "the maximum execution time is 0 or even negative");
	         this.maxTime = maxTime;
	         if (foodSources == null) throw new Exception(msg + "the reference to the starting point is null");
	         this.solution = foodSources;
	         if (obj == null) throw new Exception(msg + "the reference to the objective is null");
	         this.obj = obj;
	         this.objValue = this.obj.value(this.solution);
	         this.metaName = "ABC";
	      }
	      catch (Exception e)
	      {
	         e.printStackTrace();
	         System.exit(1);
	      }
	   }
	 
	   /*
	// Initializing food resources,fitness values, global optimum parameters
		void initialize() {
			for (int i = 0; i < FoodNumber; i++) {
				for (int j = 0; j < D; j++) {
					Foods[i][j] = Lbvec[j] + ((Ubvec[j] - Lbvec[j]) * Math.random());
				}
				Objval[i] = ff.func(Foods[i]);
			}
			Fitness = Fitcalc(Objval);

			double[] dep = getminval_index(Objval);
			fmin = dep[0];
			indexfmin = (int) (dep[1]);
			Globalmin = Objval[indexfmin];
			for (int i = 0; i < D; i++) {
				GlobalParams[i] = Foods[indexfmin][i];
			}
		}*/
	   
	 @Override
	   public void optimize()  // by ABC
	   {
		 
	   }
}