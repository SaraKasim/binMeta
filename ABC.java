package binMeta;

public class ABC extends binMeta {

	/* Control Parameters of ABC algorithm */
	int NP; /* The number of colony size (employed bees+onlooker bees) */
	int FoodNumber; /* The number of food sources equals the half of the colony size */
	int limit; /*
				 * A food source which could not be improved through "limit" trials is abandoned
				 * by its employed bee
				 */
	int maxCycle; /* The number of cycles for foraging {a stopping criteria} */

	/* Problem specific variables */
	int D; /* The number of parameters of the problem to be optimized */
	int lb; /* lower bound of the parameters. */
	int ub; // upper bound of the parameters.

	int runtime = 30; /* Algorithm can be run many times in order to see its robustness */

	Data Foods[][];
	double f[];
	double fitness[];
	int trial[];
	Data solutions[];

	double GlobalMin; // Optimum solution obtained by ABC algorithm */
	Data GlobalParams[]; /* Parameters of the optimum solution */
	double ObjValSol; // Objective function value of new solution
	double FitnessSol; // Fitness value of new solution
	int neighbour, param2change; 
	double prob[]; /*
					 * prob is a vector holding probabilities of food sources (solutions) to be
					 * chosen
					 */

	// ABC constructor
	public ABC(Data startPoint, Objective obj, int iNP, int ilimit, int imaxCycle, int iLb, int iUb) {
		this.lb = iLb;
		this.ub = iUb;
		this.NP = iNP;
		this.limit = ilimit;
		FoodNumber = iNP / 2;
		D = 1;
		this.maxCycle = imaxCycle;
		f = new double[FoodNumber];
		fitness = new double[FoodNumber];
		trial = new int[FoodNumber];
		prob = new double[FoodNumber];
		GlobalParams = new Data[D];
		Foods = new Data[FoodNumber][D];
		solutions = new Data[D];
		try {
			String msg = "Impossible to create RandomWalk object: ";
			if (startPoint == null)
				throw new Exception(msg + "the reference to the starting point is null");
			this.solution = startPoint;
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

	/* The best food source is memorized */
	void MemorizeBestSource() {
		
		int i, j;

		for (i = 0; i < FoodNumber; i++) {
			if (f[i] < GlobalMin) {
				GlobalMin = f[i];
				for (j = 0; j < D; j++) {
					GlobalParams[j] = Foods[i][j];
					if (this.objValue > obj.value(GlobalParams[j])) {
						this.objValue = obj.value(GlobalParams[j]);
						this.solution = GlobalParams[j];
					}

				}
			}
		}
	}

	void init(int index) {
		
		int j;
		Data solD = new Data(this.solution);
		for (j = 0; j < D; j++) {

			// int test=(int)(Math.random()*(ub-lb)+lb);
			Foods[index][j] = solD.randomSelectInNeighbour(solution.numberOfBits());
			// Foods[index][j] = new Data(test);
			// System.out.println("value " + test);
			solutions[j] = Foods[index][j];
		}
		
		  f[index] = obj.value(solutions[0]);
		  fitness[index] = CalculateFitness(f[index]);
		  trial[index] = 0;
		 
	}

	/* All food sources are initialized */
	void initial() {
		
		int i;
		for (i = 0; i < FoodNumber; i++) {
			init(i);
		}
		GlobalMin = f[0];
		for (i = 0; i < D; i++)
			GlobalParams[i] = Foods[0][i];
	}

	void SendEmployedBees() throws InterruptedException {
		Data solD = new Data(this.solution);
		int i, j;
		/* Employed Bee Phase */
		for (i = 0; i < FoodNumber; i++) {
			/* The parameter to be changed is determined randomly */
			param2change = (int) Math.floor(D * Math.random());// creates random integer that is element of [0,D]
			/*
			 * A randomly chosen solution is used in producing a mutant solution of the
			 * solution i
			 */
			neighbour = (int) Math.floor(FoodNumber * Math.random());// creates random integer that is element of [0,FS]

			/* Randomly selected solution must be different from the solution i */
			while (neighbour == i)// neighbour values must be different from i
			{
				neighbour = (int) Math.floor(FoodNumber * Math.random());
			}

			for (j = 0; j < D; j++) {
				solutions[j] = Foods[i][j];
			}

			/* generate a new random solution */

			solutions[param2change] = solD.randomSelectInNeighbour(solution.numberOfBits());

			/*
			 * if generated parameter value is out of boundaries, it is shifted onto the
			 * boundaries
			 */

			ObjValSol = obj.value(solutions[0]);
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
				if (this.objValue > obj.value(solutions[0])) {
					this.solution = solutions[0];
					this.objValue = ObjValSol;
				}

				f[i] = ObjValSol;
				fitness[i] = FitnessSol;
			} else { /* if the solution i can not be improved, increase its trial counter */
				trial[i] = trial[i] + 1;
			}

		}

		/* end of employed bee phase */

	}

	/*
	 * A food source is chosen with the probability which is proportioal to its
	 * quality
	 */
	/* Different schemes can be used to calculate the probability values */
	/* For example prob(i)=fitness(i)/sum(fitness) */
	/* or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b */
	/*
	 * probability values are calculated by using fitness values and normalized by
	 * dividing maximum fitness value
	 */
	void CalculateProbabilities() {
		int i;
		double maxfit;
		maxfit = fitness[0];
		for (i = 1; i < FoodNumber; i++) {
			if (fitness[i] > maxfit)
				maxfit = fitness[i];
		}

		for (i = 0; i < FoodNumber; i++) {
			prob[i] = (0.9 * (fitness[i] / maxfit)) + 0.1;
		}

	}

	void SendOnlookerBees() throws InterruptedException {
		Data solD = new Data(this.solution);
		int i, j, t;
		i = 0;
		t = 0;
		/* onlooker Bee Phase */
		while (t < FoodNumber) {

			if (Math.random() < prob[i]) /* choose a food source depending on its probability to be chosen */
			{
				t++;

				/* The parameter to be changed is determined randomly */
				param2change = (int) Math.floor(D * Math.random());// creates random number that is element of [0,D]

				/*
				 * A randomly chosen solution is used in producing a mutant solution of the
				 * solution i
				 */
				neighbour = (int) Math.floor(FoodNumber * Math.random());// creates random number that is element of
																			// [0,FoodNumber]

				/* Randomly selected solution must be different from the solution i */
				while (neighbour == i) {
					neighbour = (int) Math.floor(FoodNumber * Math.random());
				}
				for (j = 0; j < D; j++) {
					solutions[j] = Foods[i][j];
				}

				/* randomly select new solution */
				solutions[param2change] = solD.randomSelectInNeighbour(solution.numberOfBits());
				/*
				 * if generated parameter value is out of boundaries, it is shifted onto the
				 * boundaries
				 */
				ObjValSol = obj.value(solutions[0]);
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
					if (this.objValue > obj.value(solutions[0])) {
						this.solution = solutions[0];
						this.objValue = ObjValSol;
					}

					f[i] = ObjValSol;
					fitness[i] = FitnessSol;
				} else { /* if the solution i can not be improved, increase its trial counter */
					trial[i] = trial[i] + 1;
				}
			} /* if */
			i++;
			if (i == FoodNumber)
				i = 0;
		} /* while */

		/* end of onlooker bee phase */
	}

	/*
	 * determine the food sources whose trial counter exceeds the "limit" value. In
	 * Basic ABC, only one scout is allowed to occur in each cycle
	 */
	void SendScoutBees() throws InterruptedException {
		int maxtrialindex, i;
		maxtrialindex = 0;
		for (i = 1; i < FoodNumber; i++) {
			if (trial[i] > trial[maxtrialindex])
				maxtrialindex = i;
		}
		if (trial[maxtrialindex] >= limit) {
			init(maxtrialindex);
		}
	}

	@Override
	public void optimize() {// by ABC
		int iter = 0;
		int run = 0;

		for (run = 0; run < runtime; run++) {
			initial();
			MemorizeBestSource();
			for (iter = 0; iter < maxCycle; iter++) {
				try {
					SendEmployedBees();
					CalculateProbabilities();
					SendOnlookerBees();
					MemorizeBestSource();
					SendScoutBees();
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
				}
			}
		}
	}

	// ThreadEm for Employed bee phase
	static class ThreadEm extends Thread {

		ABC abc;

		// overriding the run method
		@Override
		public void run() {

			try {
				//System.out.println("Thread start " + Thread.currentThread().getName());

				for (int run = 0; run < abc.runtime; run++) {
					abc.initial();
					abc.MemorizeBestSource();

					for (int iter = 0; iter < abc.maxCycle; iter++) {
						abc.SendEmployedBees();
					}
				}
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
			}

		}

	}

	// Thread for onlooker phase
	static class ThreadOn extends Thread {

		ABC abc;

		// overriding the run method
		@Override
		public void run() {
			//System.out.println("Thread start " + Thread.currentThread().getName());

			try {
				for (int run = 0; run < abc.runtime; run++) {
					abc.initial();
					abc.MemorizeBestSource();

					for (int iter = 0; iter < abc.maxCycle; iter++) {

						abc.CalculateProbabilities();
						;
						abc.SendOnlookerBees();
					}
				}
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
			}

		}

	}

	// Thread for scout phase
	static class ThreadSc extends Thread {

		ABC abc;

		// overriding the run method
		@Override
		public void run() {

			try {
				//System.out.println("Thread start " + Thread.currentThread().getName());

				for (int run = 0; run < abc.runtime; run++) {
					abc.initial();
					abc.MemorizeBestSource();

					for (int iter = 0; iter < abc.maxCycle; iter++) {

						abc.MemorizeBestSource();
						abc.SendScoutBees();
					}
				}
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
			}

		}

	}

	// main
	public static void main(String[] args) {
		// BitCounter
		int n = 50;
		int swarmSize = 6;
		int lim = 10;
		int cycles = 5;
		int lb = 1;
		int ub = 1048575;
		Objective obj = new BitCounter(n);
		Data D = obj.solutionSample();
		ABC abc = new ABC(D, obj, swarmSize, lim, cycles, lb, ub);
		System.out.println(abc);
		System.out.println("starting point : " + abc.getSolution());
		System.out.println("optimizing ...");

		ThreadEm thEm = new ThreadEm();
		ThreadOn thOn = new ThreadOn();
		ThreadSc thSc = new ThreadSc();

		thEm.abc = abc;
		thOn.abc = abc;
		thSc.abc = abc;

		thEm.start();
		thOn.start();
		thSc.start();
		
		System.out.println(abc);
		System.out.println("multi threaded solution : " + abc.getSolution());
		System.out.println();

		abc.optimize();
		System.out.println(abc);
		System.out.println("solution : " + abc.getSolution());
		System.out.println();

		// Fermat
		int exp = 2;
		int ndigits = 10;
		obj = new Fermat(exp, ndigits);
		Data DF = obj.solutionSample();
		ABC abcF = new ABC(DF, obj, swarmSize, lim, cycles, lb, ub);
		System.out.println(abcF);
		System.out.println("starting point : " + abcF.getSolution());
		System.out.println("optimizing ...");
		abcF.optimize();
		System.out.println(abcF);
		System.out.println("solution : " + abcF.getSolution());
		Data x = new Data(abcF.solution, 0, ndigits - 1);
		Data y = new Data(abcF.solution, ndigits, 2 * ndigits - 1);
		Data z = new Data(abcF.solution, 2 * ndigits, 3 * ndigits - 1);
		System.out.print(
				"equivalent to the equation : " + x.posLongValue() + "^" + exp + " + " + y.posLongValue() + "^" + exp);
		if (abcF.objValue == 0.0)
			System.out.print(" == ");
		else
			System.out.print(" ?= ");
		System.out.println(z.posLongValue() + "^" + exp);
		System.out.println();

		// ColorPartition
		n = 4;
		int m = 14;
		ColorPartition cp = new ColorPartition(n, m);
		Data DCP = cp.solutionSample();
		ABC abcCP = new ABC(DCP, cp, swarmSize, lim, cycles, lb, ub);
		System.out.println(abcCP);
		System.out.println("starting point : " + abcCP.getSolution());
		System.out.println("optimizing ...");
		abcCP.optimize();
		System.out.println(abcCP);
		System.out.println("solution : " + abcCP.getSolution());
		cp.value(abcCP.solution);
		System.out.println("corresponding to the matrix :\n" + cp.show());

	}

}

