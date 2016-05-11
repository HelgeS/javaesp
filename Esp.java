//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: Esp.java,v $
//////////////////////////////////////////////////////////////////////////
import java.util.*;
import java.io.*;

public class Esp {
	private static final int FEEDFORWARD = 0;
	private static final int SIMPLERECURRENT = 1;
	private static final int SECONDORDER = 2;
	private static final int FULLYRECURRENT = 3;
	
//	public static int NUM_IN;          	// number of input from environment
//	public static int NUM_INPUTS;      	// number of input units; includes recurrent input
//	public static int NUM_OUTPUTS;     	// number of output units 
//	private static int GENE_SIZE;        	// size is chromosomes
	private static int maxGenerations = 10000;
	private static boolean SKIP = false;	// skip recombination;
	private static double MUT_RATE = 0.4;	// mutation rate in %
	private static int STAG = 20;			// the length of the performance queue.
	private static boolean PLOTTING = false;
	private static boolean DELTA = false;
	public static boolean MIN = false;
	private static long seed;
	private static int seedNet;
	private static boolean analyze;
	private RandomCauchy rndCauchy = new RandomCauchy();
	public Network bestNetwork;//make private after debug
	public Network phaseBest;  
	private int subPopulationCount;	// The number of subpopulations in the gene pool
	private int subPopulationSize;
	private int numTrials;
	private int generation;
	//  espParams params;
	private int netType;
	private List perfQ = new Vector();	// double
	private Environment environ;
	public SubPopulation[] subpops;	// SubPopulation
 			// neurons 
	// constructor
	public Esp() {
	}
	
	public void init( Environment e ) {
		environ = e;
		generation = 0;          // start at generation 0
		numTrials = subPopulationSize * 10;  // # of trials ~10/neuron
		
		bestNetwork = genNetType(); // "              "
		bestNetwork.create(); 
		
		phaseBest = genNetType();   // generate new net 
		phaseBest.create();               // create net  
		
		subpops = new SubPopulation[subPopulationCount];
		for( int i = 0; i < subpops.length; ++i) {  // construct subPopulationCount # of subPops
    		subpops[i] = new SubPopulation(subPopulationSize, bestNetwork);
		}
	}

//----------------------------------------------------------------------
// create the subpopulations of neurons, initializing them to 
//   random weights
	public void create() { // creates a random population of sub-populations
		for( int i=0; i<subpops.length; ++i ) {
    		subpops[i].create(); // create each subpop
    	}
	}
	
//////////////////////////////////////////////////////////////////////
//
//  evolve is the main genetic function.  The subpopulations are first
//  evaluated in the given task.  Then for each subpop the neurons are
//  ranked (sorted by fitness) and recombined using crossover.  
//  Mating is allowed only between the top 'numBreed' members of 
//  the same subpop. The neurons are then mutated.
//
//////////////////////////////////////////////////////////////////////
// cycles: the number of generations to evolve
	public void evolve(int cycles) {
		while( ( SignalHandler.gInterrupt == false ) && ( ++generation < cycles ) ) {
			evalPop();                      	//evaluate neurons
			if( SKIP ) {
				SKIP = false;					// skip recombination if we have just 
			} else {            				//    perturb the population
				for(int j = 0; j < subpops.length; ++j) { 
					subpops[j].qsortNeurons();
		//      	for(int k = 0; k < subpops[j].numNeurons; ++k) 
		//      		System.out.print("%f ", subpops[j].getNeuron(k).fitness);
		//      	System.out.println();
					subpops[j].recombine();
				}
					      
				//printNeurons();
						
				// mutate population 
				for( int j = 0; j < subpops.length; ++j ) {
					subpops[j].mutate(MUT_RATE);
				}
			}
		}
		endEvolution();
	}
	
//////////////////////////////////////////////////////////////////////
//
//  eval_pop evaluates entire population to get the average fitness of
//  each neuron.  This function is called by evolve and should put the
//  fitness of each neuron in the neuron's fitness field
//  (getNeuron(i).fitness).
//
//////////////////////////////////////////////////////////////////////
	public void evalPop() {
		evalReset();   // reset fitness and test values of each neuron
		performEval(); // evaluate the networks
		average();     // normalize neuron fitnesses
	}
	
	public void findChampion() {
		double fitness, best=-999999;
		StringBuffer subname = new StringBuffer("000best.bin");
		DataInputStream fptr;
		SubPopulation net;
		  
		for( int l = 0; l < 1000; ++l ) {
			char j = (char)(l/100);
			subname.setCharAt(0, (char)(j + 48));
			j = (char)(l%100);
			subname.setCharAt(1, (char)(j/10 + 48));
			subname.setCharAt(2, (char)(j%10 + 48));
			try {
				fptr = new DataInputStream( new FileInputStream( subname.toString() ) );
				fptr.close();
				net = new SubPopulation(subname.toString());
    			fitness = 0;
			      
    			// evaluate network
    			//ccfitness += state.apply(net);
    			if( fitness > best) {
        			best = fitness;
        		}
				PrintWriter pw = new PrintWriter( new FileOutputStream( "report.txt", true ) );
				pw.write(l+" ");
    			// Aluminum stuff 
    			//  fprintf(fptr, "%f %f",
    			//      ((totdeltasum + deltamaara) / nsulatus) / (MAXSUBS + 1),
    			//     deltamaara / nsulatus);
    			//for( j = 0; j < MAXSUBS; j++)
    			//	fprintf(fptr, " %f", deltasum[j] / nsulatus);
    			//fprintf(fptr, "\n");

    			pw.close();
			} catch( IOException e ) {
				System.out.println(e);
			}
		}
		  
		try {
			PrintWriter pw = new PrintWriter( new FileOutputStream( "analyze.out", true ) );
			pw.println(1/best);
			pw.close();
		} catch( IOException e ) {
			System.out.println(e);
		}
		System.out.println(1/best);		
	}  
	
//----------------------------------------------------------------------
//
	public void loadSeedNet(String filename) {
		netType = SECONDORDER;  //change this to read the type from the filename;
		subPopulationCount = bestNetwork.load(filename);
		System.out.println("NUM_POPS " + subPopulationCount);
		System.out.println("GENE_SIZE " + bestNetwork.getGeneSize());
	}
	
//----------------------------------------------------------------------
// add a new subbpop. In effect: add a unit to the networks.
	public void addSubPop() {
		SubPopulation newSp;    //pointer to new subpop.
		addConnections();
		 
		for(int i=0; i < subpops.length; ++i) { //freeze best network? false=yes
    		subpops[i].evolvable = true;    // true = no.
    	}
		++subPopulationCount;             // inc # subpops
//		setupNetDimensions();  //adjust NUM_INPUTS, GENE_SIZE
		phaseBest.addNeuron();  // add a neuron to the phasebest and best net.
		//  phaseBest.fitness = 0;
		bestNetwork.addNeuron(); 
		newSp = new SubPopulation(subPopulationSize, bestNetwork); //constuct and create
		newSp.create();                         //  new subpop

//		subpops.add( newSp );		// put its pointer in Esp's vector of sp pointers.
		SubPopulation[] temp = new SubPopulation[subpops.length+1];
		System.arraycopy( subpops, 0, temp, 0, subpops.length );
		temp[ temp.length-1 ] = newSp;
		subpops = temp;
		
		System.out.println( "Adding SUBPOP " + (subPopulationCount-1) );
	}
	
	public void endEvolution() {
		printStats();
		System.out.println( "BYE!" );
		System.exit(0);
	}


//	private Esp(final Esp) // just for safety
//	private Esp(final Esp)

//////////////////////////////////////////////////////////////////////
//
//  evaluation stage.  Create numTrials networks, each containing
//  ZETA neurons selected randomly from the population.  Each network
//  is evaluated and each participating neuron receives the fitness
//  evaluations of each network it parcipates in.
//
//////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------
// output a new network of the appropriate type
	private Network genNetType() {
		Network rc = null;
		switch (netType) {
			case FEEDFORWARD: {
    			rc = ( new FeedForwardNetwork( environ.getInputSize(), subPopulationCount, environ.getOutputSize()) );
    		} break;
    		case FULLYRECURRENT : {
    			rc = ( new FullyRecurrentNetwork( environ.getInputSize(), subPopulationCount, environ.getOutputSize()) );
    		} break;
			case SIMPLERECURRENT : {
    			rc = ( new SimpleRecurrentNetwork( environ.getInputSize(), subPopulationCount, environ.getOutputSize()) );
    		} break;
			case SECONDORDER : {
    			rc = ( new SecondOrderRecurrentNetwork( environ.getInputSize(), subPopulationCount, environ.getOutputSize()) );
    		} break;
    		default: {
    			rc = ( new FullyRecurrentNetwork(0, 0, 0) );
    		} break;
		}
		return rc;
	}
	
//----------------------------------------------------------------------
// reset fitness vals
	private void evalReset() {
		for( int i = 0; i < subpops.length; ++i) {
    		subpops[i].evalReset();
    	}
	}
	
//--------------------------------------------------------------------
// evaluate the networks on the task
	int evaluations = 0;
	private void performEval() {
		Network net = genNetType(); 
		Network bestNet = genNetType(); 
		//bestNet.create();
		for( int i = 0; i < numTrials; ++i ) {
			// find random subpopulation
			net.setRandomNetwork( subpops );
			net.incrementTests();
				    
			++evaluations;
			//evaluate the network
			net.fitness = (double) environ.evalNet(net);
				     
			//net.printWeight(ALL_NETS);

			//    System.out.println("fit %d  %f",i, net.fitness);
			if( MIN ) { 
				if( net.fitness < bestNet.fitness) {
					bestNet = (Network)net.clone();
				}
			} else if(net.fitness > bestNet.fitness) { 
				bestNet = (Network)net.clone();
			}

			net.addFitness(); // add network fitness to its neurons 
		} 
		if(MIN) { 
			if(bestNet.fitness < phaseBest.fitness) {
				phaseBest = (Network)bestNet.clone(); 
			}
		} else if(bestNet.fitness > phaseBest.fitness) { 
			phaseBest = (Network)bestNet.clone(); 
		}
			  
		perfQ.add(0, new Double(phaseBest.fitness));

		if( MIN ) {  
			if( bestNet.fitness <= environ.getMaxFitness() ) {
				phaseBest.fitness = 1000000;
				perfQ.clear();
				environ.nextTask();
			}
		} else if(bestNet.fitness >= environ.getMaxFitness() ) {
			phaseBest.fitness = 0.0;
			perfQ.clear();
			environ.nextTask();
		}

		// if performance stagnates, do something
		if( perfQ.size() >= STAG) {
			if(MIN) {
				if(bestNet.fitness >= ((Double)perfQ.get(STAG-1)).doubleValue()) {
					handleStagnation();
				}
			} else if(bestNet.fitness <= ((Double)perfQ.get(STAG-1)).doubleValue()) {
				handleStagnation();
			}
		}

		//System.out.println("%f %f" + perfQ.front() + perfQ[STAG-1] );
		System.out.println( "gen " + generation + ": best " + bestNet.fitness + ", task best " + phaseBest.fitness 
			//+ ", overall best " + bestNetwork.fitness 
			);
	}
	
//----------------------------------------------------------------------
// get average fitness level.  
	private void average() {
		for( int i = 0; i < subpops.length; ++i)  {
    		subpops[i].average();
    	}
	}
	
/////////////////////////////////////////////////////////////////////
//
// ESP I/O fns
	private void printNeurons() {
		String filename = "./"+seed/*getpid()*/+"_"+subPopulationCount+"_"+generation+".pop";

		try {
			PrintWriter fptr = new PrintWriter( new FileOutputStream( filename ) );
			for(int i=0; i < subpops.length; ++i) {
    			subpops[i].printWeight(fptr);
			}
			fptr.close();
		} catch( IOException e ) {
    		System.out.println(" Error - cannot open " + filename);
			System.out.println( e );
			System.exit(1);
		}

	}
	
//----------------------------------------------------------------------
// make a decision about what to do when performace stagnates
	private void handleStagnation() {
  perfQ.clear();   
  
  if( phaseBest.fitness == bestNetwork.fitness ) {
    if( removeSubPop( phaseBest ) ) {
      bestNetwork = (Network)phaseBest.clone();
      newDeltaPhase();
    }
    else
      addSubPop();
    environ.simplifyTask();
    phaseBest.fitness = 0.0;
  } else {
    bestNetwork = (Network)phaseBest.clone();
    //sprintf(filename, "./NETS/net%lu_%d_%d", getpid(), subpops.length, generation);
    //nets[b].saveText(filename);
    newDeltaPhase();    
  }
}
	
//----------------------------------------------------------------------
// Start a new delta phase 
	private void newDeltaPhase() {
		System.out.println( "DELTA started");
		SKIP = true;
		for(int k=0; k < subpops.length; ++k) {
    		subpops[k].deltify( bestNetwork.getNeuron(k) );
		}
	}
	
//----------------------------------------------------------------------
// opposite of addSubPop.
	private boolean removeSubPop(int sp) {
		boolean rc = false;
		if( sp >= 0) {
    		removeConnections(sp);
    		--subPopulationCount;
//    		setupNetDimensions();
    		phaseBest.removeNeuron(sp);
    		bestNetwork.removeNeuron(sp);
//    		delete subpops[sp];
//    		subpops.erase(subpops.begin()+sp);      
			SubPopulation[] temp = new SubPopulation[subpops.length-1];
			System.arraycopy( subpops, 0, temp, 0, sp );
			System.arraycopy( subpops, sp+1, temp, sp, temp.length-sp );
			subpops = temp;
		      
    		System.out.println( "Remove SUBPOP " + sp);
    		System.out.println( "Now " + subpops.length + " subpops");
    		rc = true;
    	} 
    	return rc;
	}
	
//----------------------------------------------------------------------
	private boolean removeSubPop(Network net) {
		int sp = net.lesion(environ);  //see which (if any) subpop can be removed
		return ( removeSubPop(sp) );  // remove it.
	}

	
//----------------------------------------------------------------------
// add the appropriate connection wieghts that are needed to add
//   a unit to each type of network.
	private void addConnections() {
		int i,j;
		switch (netType) {
		case FEEDFORWARD:
    		break;
		case FULLYRECURRENT:
		case SIMPLERECURRENT:
    		for(i=0; i < subpops.length; ++i)  // add connection to neurons in spops
    			subpops[i].addConnection(environ.getInputSize() + subPopulationCount);
    		break;
		case SECONDORDER:
    		for(i=0; i < subpops.length; ++i)  // add connection to neurons in spops
    			for(j=1; j < environ.getInputSize()+1; ++j)  // add connection to neurons in spops
					subpops[i].addConnection(subpops.length*j+j-1);
    		break;
		}
	}
	
//----------------------------------------------------------------------
// opposite of addConnections.
	private void removeConnections(int sp) {
  int i,j;
  switch (netType) {
  case FEEDFORWARD:
    break;
  case FULLYRECURRENT:
  case SIMPLERECURRENT:
    for(i=0; i < subpops.length; ++i)  // remove connection to neurons in spops
      subpops[i].removeConnection(environ.getInputSize() + subPopulationCount-1);
    break;
  case SECONDORDER:
    for(i=0; i < subpops.length; ++i) 
      for(j=1; j < environ.getInputSize()+1; ++j)  // remove connection to neurons in spops
	subpops[i].removeConnection((subpops.length-1)*j);
  }
}

	
	private void printStats() {
		System.out.println();
		System.out.println("Total number of network evaluations : "+ (generation * numTrials));
	}
	
	
	static void parseArgs( String[] args, Esp esp ) {
		// Command-line defaults.
		final int DEFAULT_NUM_POPS = 10; 
		final int DEFAULT_POP_SIZE = 100;
		  
		// Set defaults.
		esp.subPopulationCount = DEFAULT_NUM_POPS;  
		esp.subPopulationSize = DEFAULT_POP_SIZE;
		esp.netType = FEEDFORWARD;
		esp.seed = System.currentTimeMillis();//getpid();
		
		for( int argc = 0; argc < args.length; argc++ ) {
			if( args[argc].equals( "--plot" ) ) {
				PLOTTING = true;
			} else if(  args[argc].equals( "--delta" ) ) {
				DELTA = true;
			} else if(  args[argc].equals( "--minimize" ) ) {
				MIN = true;
			} else if(  args[argc].equals( "--stag" ) ) {
				STAG = Integer.valueOf( args[argc+1] ).intValue();
				argc++;
	    		System.out.println("STAG : "+ STAG);
			} else if(  args[argc].equals( "--mutation" ) ) {
				MUT_RATE = Double.valueOf( args[argc+1] ).doubleValue();
				argc++;
	    		System.out.println("MUT_RATE : "+ MUT_RATE);
			} else if(  args[argc].equals( "-z" ) ) {
				esp.subPopulationCount = Integer.valueOf(args[argc+1]).intValue();
				argc++;				
			} else if(  args[argc].equals( "-n" ) ) {
				esp.subPopulationSize = Integer.valueOf(args[argc+1]).intValue();
				argc++;
			} else if(  args[argc].equals( "-s" ) ) {
				esp.seed = Integer.valueOf(args[argc+1]).intValue();
				argc++;
			} else if(  args[argc].equals( "-t" ) ) {
				esp.netType = Integer.valueOf(args[argc+1]).intValue();
				argc++;

			//Additions by Jacob Schrum, August 2007
			} else if(  args[argc].equals( "--maxGens" ) ) {
				maxGenerations = Integer.valueOf( args[argc+1] ).intValue();
				argc++;
	    			System.out.println("MAX GENERATIONS : "+ maxGenerations );
            		} else if(  args[argc].equals( "-toFile" ) ) {
                	
				try{
                		    System.setOut( new PrintStream( new FileOutputStream(args[argc+1]) ) );
                		}catch(Exception e){
				    System.out.println("Must provide valid filename.");
                		    System.out.println(e);
                		}
                		argc++;
			//End Additions by Jacob Schrum


			} else if(  args[argc].equals( "-x" ) ) {
				esp.analyze = true;
			} else if(  args[argc].equals( "-f" ) ) {
				esp.seedNet = Integer.valueOf(args[argc+1]).intValue();
				argc++;
			} else if(  args[argc].equals( "-h" ) ) {
    			System.out.println( "Usage: " + args[0] + " <options>");
    			System.out.println( "options:");
    			System.out.println( "  -z Z : Z number of hidden units.");
    			System.out.println( "  -n N : N neurons in each subpop");
    			System.out.println( "  -s S : Initial integer seed S.");
    			System.out.println( "  -t T : Type of Network.");
    			System.out.println( "         0 = Feed Forward, 1 = Simple Recurrent, 2 = 2ndOrder Recurrent");
    			System.out.println( "         3 = Fully Recurrent");
    			//System.out.println( "  -f F : Seed network. ");
    			System.out.println( "  --mutation M : Mutation rate in percent (0,1). ");
    			System.out.println( "  --stag STAG : Stagnation threshold, number of generations. ");
    			System.out.println( "  --minimize : Minimize fitness. ");
    			//      System.out.println( "  -x   : analyze nets.");

			//Additions by Jacob Schrum, August 2007
			System.out.println( "  --maxGens MG : Maximum generations to run for. ");
			System.out.println( "  --toFile FN : Output to file FN instead of console. ");
			//End Additions by Jacob Schrum

    			System.out.println( "  -h   : this message.");
    			System.out.println( 
    				"Default: " + args[0]
					+ " -z " + DEFAULT_NUM_POPS
					+ " -n " + DEFAULT_POP_SIZE
					+ " -s getpid()"
					+ " -t " + FEEDFORWARD
					+ " --mutation 0.4"
					+ " --stag 20"
					+ " --maxGens " + maxGenerations );
    			System.out.println();
    			System.exit(0);
			}
		}		
	}
 
	static void checkArgs (int subPopulationCount, int popSize) {
		if(MIN) {
			System.out.println( "MINIMIZING");
		}
	}

	public static void main(String[] args)	{
		// Parse command-line arguments.
		int seedNet = 0;
		int netType;
		boolean analyze = false;
		
		synchronized( SignalHandler.class ) {
			SignalHandler.setUpCtrlC();
			Esp esp = new Esp(); //construct the ESP.
			parseArgs (args, esp);
			checkArgs (esp.subPopulationCount, esp.subPopulationSize);
			RandomSingleton.getInstance().setSeed(esp.seed);
			CartPole cartpole = new CartPole(0); // this is our experiment.
			esp.init(cartpole);
			  
	//		Esp esp = new Esp(subPopulationCount, subPopulationCount, cartpole, netType); //construct the ESP.

			if( analyze ) { 
    			esp.findChampion();
    			System.exit(1);
			}				
			esp.create();                     // create the neurons.
			esp.evolve(maxGenerations);       // evolve them.
		}
	}
}
