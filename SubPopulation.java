//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: SubPopulation.java,v $
//////////////////////////////////////////////////////////////////////////
import java.util.*;
import java.io.*;

public class SubPopulation /*implements Cloneable*/ {
	//protected:
	private Neuron[] pop;  // Neuron
	private int geneCount;

	private static RandomCauchy rndCauchy = new RandomCauchy();
	
	public SubPopulation(int populationCount, Network net) {
		this.geneCount = net.getGeneSize();
		pop = new Neuron[populationCount];
		evolvable = true;
		numBreed = (int) populationCount/4;   //we leave subPop::create to 
		// actually create the neurons.  This way neurons are not created
		// when we initially create a Network
		// for (int i = 0; i < getGeneCount(); ++i)
		// pop[i] = new neuron();
	}
/*	
	public Object clone() {
		try {
			SubPopulation rc = (SubPopulation)super.clone();
			rc.pop = (Neuron[])pop.clone();
			for( int i=0; i < rc.pop.length; i++ ) {
				rc.pop[i] = (Neuron)rc.pop[i].clone();
			}
			return rc;
		} catch( CloneNotSupportedException e ) {
			return null;
		}
	}
*/	
	public SubPopulation( String fname ) {
		try {
			DataInputStream fptr = new DataInputStream( new FileInputStream( new File( fname ) ) );
			int populationSize = fptr.readInt();
			geneCount = fptr.readInt();
			pop = new Neuron[populationSize];
			for( int i=0; i < pop.length; i++ ) {
				pop[i] = new Neuron(getGeneCount());
			}
			evolvable = true;
			numBreed = (int) populationSize/4;   //we leave subPop::create to 
			
			fptr.close();
		} catch( IOException e ) {
    		System.out.println(" Error - cannot open " + fname);
			System.out.println( e );
			System.exit(1);
		}
	}
	
	public int getPopulationSize() {
		return pop.length;
	}
	
	public int getGeneCount() {
		return geneCount;
	}
/*	
	public void addNeuron() {
		Neuron[] temp = new Neuron[ pop.length+1 ];
		for( int i = 0; i < pop.length; i++ ) {
			temp[i] = pop[i];
		}
		temp[pop.length] = new Neuron( getGeneCount() );
		pop = temp;
	}
	
	public void removeNeuron( int index ) {
		Neuron[] temp = new Neuron[ pop.length-1 ];
		for( int i = 0; i < index; i++ ) {
			temp[i] = pop[i];
		}
		for( int i = index+1; i < pop.length; i++ ) {
			temp[i-1] = pop[i];
		}
		pop = temp;
	}
	
	private void setNeuron( int index, Neuron n ) {
		pop[index] = (Neuron)n.clone();
	}
*/	
	public Neuron getNeuron( int index ) {
		return pop[index];
	}
	
	//  SubPopulation(const SubPopulation &s);
//---------------------------------------------------------------------
// create the neurons, initial their weights, and put them in the subpop.
	public void create() { // creates a random subpopulation of neurons
		if(evolvable) {
			for (int i = 0; i < pop.length; ++i) {
				pop[i] = new Neuron(getGeneCount());
				pop[i].create();
			}
		}
	}
	
//----------------------------------------------------------------------
// reset fitness and test vals of all neurons
	public void evalReset() {
		for(int i=0; i < pop.length; ++i)  {
    		pop[i].fitness = 0;
    		pop[i].tests = 0;
		}
	}

//----------------------------------------------------------------------
// select a neuron at random
	public Neuron selectNeuron() {
		int random = Math.abs(RandomSingleton.getInstance().nextInt());
		int length = pop.length;
		int index = random % length;
		return pop[ index ];
	}
	
//----------------------------------------------------------------------
// normalize the neuron fitnesses 
	public void average() {
		for(int i = 0; i < pop.length; ++i) {
    		if( pop[i].tests != 0 ) {
	  			pop[i].fitness = pop[i].fitness/pop[i].tests;
    		} else {
	    		pop[i].fitness = 0;
    		}
		}
	}
	
//----------------------------------------------------------------------
// sort the neurons in each subpop using quicksort.
	private static final Comparator minimize_fit = new MinimizeFit();
	private static final Comparator maximize_fit = new MaximizeFit();
	public void qsortNeurons() {
		if( Esp.MIN ) {
			Arrays.sort( pop, minimize_fit );
		} else {
			Arrays.sort( pop, maximize_fit );
    	}
	}
	
//----------------------------------------------------------------------
// recombine neurons with members of their subpop using crossover.
	public void recombine() {
		for (int i = 0; i < numBreed; ++i) {
			int mate = findMate(i);
    		crossover(pop[i].weight, pop[mate].weight, 
	    		pop[pop.length-(1+i*2)].weight, 
	    		pop[pop.length-(2+i*2)].weight);
	    }
	}
	
//----------------------------------------------------------------------
// mutate half of the neurons with cauchy noise.
	public void mutate(double rate) {
		for (int i = numBreed*2 ; i < pop.length; ++i)  {
			if( RandomSingleton.getInstance().nextDouble() < rate) {
    			pop[i].weight[Math.abs(RandomSingleton.getInstance().nextInt())%getGeneCount()] += rndCauchy.randomFunction( 0.3 );
    		}
    	}
	}
	
//---------------------------------------------------------------------
// used to perform "delta-coding" like 
	public void deltify(Neuron bestNeuron) {
		//  neuron tmp = *pop[0];
		for(int i= 0; i < pop.length; ++i) {
			pop[i].perturb( bestNeuron ); 	// make each neuron a perturbation of the
												// neuron in the best network that 
												// corresponds to that subpop
			//delete pop[i];
			//if(i > (pop.length/2))
			//  pop[i] = bestNeuron.perturb();    
			//else
			//pop[i] = tmp.perturb();
		}
	}
	
//	public void save(String fname);
// add a weight at position 'locus' for all neuron in the subpop.
	public void addConnection(int locus) {
		// TODO: Make more efficient
    	for(int i=0; i < pop.length; ++i) {
    		pop[i].weight = doubleArrayInsert(pop[i].weight, 1.0, locus);
		}
	}
	
//----------------------------------------------------------------------
// opposite of addConnection.
	public void removeConnection(int locus) {
		if( locus < pop.length ) {
    		for( int i = 0; i < pop.length; ++i ) {
    			pop[i].weight = doubleArrayRemove( pop[i].weight, locus );
    		}
    	}
	}	
	
	public void saveBin(String fname) {
		try {
			DataOutputStream dos = new DataOutputStream( new FileOutputStream( fname ) );
			dos.writeInt( pop.length );
			dos.writeInt( getGeneCount() );

			for( int i = 0; i < pop.length; ++i ) {
    			for( int j = 0; j < getGeneCount(); ++j ) {
    				dos.writeDouble(pop[i].weight[j]);
    			}
    		}
    		dos.close();
		} catch( IOException e ) {
			System.out.println();
    		System.out.print("Error - cannot open " + fname);
    		System.exit(1);
		}
	} 

	public void saveText(String fname) {
		try {
			PrintWriter pw = new PrintWriter( new FileOutputStream( fname ) );
			pw.println( pop.length );
			pw.println( getGeneCount() );
			for( int i = 0; i < pop.length; ++i ) {
    			for( int j = 0; j < getGeneCount(); ++j ) {
    				pw.print( pop[i].weight[j] );
    			}
    			pw.println();
			}
			pw.close();
		} catch( IOException e ) {
			System.out.println();
    		System.out.print(" Error - cannot open " + fname);
    		System.exit(1);
		}
	}
	
	public int load(String filename) {
		try {
			BufferedReader br = new BufferedReader( new InputStreamReader( new FileInputStream( filename ) ) ) ;
			StreamTokenizer st = new StreamTokenizer( br );
			st.nextToken();
			pop = new Neuron[(Integer.valueOf(st.toString())).intValue()];
			geneCount = (Integer.valueOf(st.toString())).intValue();
			for( int i = 0; i < pop.length; ++i ) {
				pop[i] = new Neuron(getGeneCount());
    			for( int j = 0; j < getGeneCount(); ++j ) {
					st.nextToken();
    				pop[i].weight[j] = Double.valueOf(st.toString()).doubleValue();
    			}
			}
			br.close();
		} catch( IOException e ) {
			System.out.println();
    		System.out.print(" Error - cannot open " + filename);
    		System.exit(1);
		}
		return pop.length;
	}

	public void print() {
		int wSize = pop[0].weight.length;
		System.out.println( "GENE_SIZE " + wSize );
		for(int i=0; i < pop.length; ++i){
    		for(int j=0; j < wSize; ++j) {
    			System.out.print( pop[i].weight[j]  );
    		}
			System.out.println();
		}
	}
	
	public void printWeight( PrintWriter pw ) {
		int wSize = pop[0].weight.length;
		for( int i=0; i < pop.length; ++i ){
    		for( int j=0; j < wSize; ++j ) {
    			pw.print( pop[i].weight[j]+" " );
    		}
    		//    fprintf(file,"\n");
    		pw.close();
		}
		pw.println();
	}
	
//	public void printDelta(File file);
	  
	//private:
	public boolean evolvable;

	private int numBreed;  //number of neurons to be mated in SubPopulation
//----------------------------------------------------------------------
// 1-point crossover
	private void crossover(	final double[] parent1, final double[] parent2, double[] child1, double[] child2) {
		//find crossover point
		int cross1 = Math.abs(RandomSingleton.getInstance().nextInt())%getGeneCount();
		System.arraycopy( parent2, 0, child1, 0, child1.length );
		System.arraycopy( parent1, 0, child2, 0, child2.length );
		
		for( int i = 0; i < cross1; i++ ) {
			double temp = child1[i];
			child1[i] = child2[i];
			child2[i] = temp;
		}
	}

//----------------------------------------------------------------------
// randomly find a mate in the same subpop.
	private int findMate(int num) {
		if( num == 0 ) {
    		return (Math.abs(RandomSingleton.getInstance().nextInt())%numBreed);
		} else {
    		return (Math.abs(RandomSingleton.getInstance().nextInt())%num);
    	}
	}
	
	private double[] doubleArrayInsert( double[] array, double value, int position ) {
		double[] rc = new double[ array.length+1 ];
		System.arraycopy( array, 0, rc, 0, position );
//		for( int i = 0; i < position; i++ ) {
//			rc[i] = array[i];
//		}
		rc[position] = value;
		System.arraycopy( array, position, rc, position+1, array.length-position );
//		for( int i = position; i < array.length; i++ ) {
//			rc[i+1] = array[i];
//		}
		return rc;
	}
		
	private double[] doubleArrayRemove( double[] array, int position ) {
		double[] rc = new double[ array.length-1 ];
		System.arraycopy( array, 0, rc, 0, position );
//		for( int i = 0; i < position; i++ ) {
//			rc[i] = array[i];
//		}
		System.arraycopy( array, position+1, rc, position, array.length-position-1 );
//		for( int i = position+1; i < array.length; i++ ) {
//			rc[i-1] = array[i];
//		}
		return rc;
	}
}