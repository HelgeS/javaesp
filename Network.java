//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: Network.java,v $
//////////////////////////////////////////////////////////////////////////
import java.io.*;
import java.util.Arrays;

public abstract class Network implements Cloneable {
	private int inputCount;
	private int outputCount;
	protected double[] activation;
	public double fitness;
	final static private double LESION_THRESHOLD = 0.9;
	private Neuron[] neurons;  // Neuron
//	private int geneCount;
	
	//Functions added by Jacob Schrum, August 2007
	public double[] getActivation() { return activation; }
	public Neuron[] getNeurons() { return neurons; }
	//End functions by Jacob Schrum

	static double sigmoid(double x, double slope ) {
		return (1/(1+Math.exp(-(slope * x)) ) );
	}
 
	static double sigmoid(double x ) {
		return sigmoid( x, 1.0 );
	} 
	  
	// , getGeneSize()
	public Network( int inputCount, int neuronCount, int outputCount ) {
		this.inputCount = inputCount;
		this.outputCount = outputCount;
		activation = new double[neuronCount]; 
		if( Esp.MIN ) {
    		fitness = 10000000;
		} else { 
    		fitness = 0.0; 
		}
		neurons = new Neuron[neuronCount];
		for( int i=0; i < neurons.length; i++ ) {
			neurons[i] = new Neuron(getGeneSize());
		}
	}
	
	public Object clone() {
		try {
			Network rc = (Network)super.clone();
			rc.activation = (double[])rc.activation.clone();
			rc.neurons = (Neuron[])rc.neurons.clone();
			for( int i=0; i < rc.neurons.length; i++ ) {
				rc.neurons[i] = (Neuron)rc.neurons[i].clone();
			}
			return rc;
		} catch( CloneNotSupportedException e ) {
			return null;
		}
	}
	
	public void copy( Network net ) {
		
	}
/*	  
	public Network( Network n ) {
		inputCount = n.inputCount;
		outputCount = n.outputCount;
		activation = new double[n.activation.length];
		for( int i = 0; i < activation.length; i++ ) {
			activation[i] = n.activation[i];
		}
		fitness = n.fitness;
		neurons = new Neuron[n.getNeuronCount()];
		for( int i = 0; i < n.getNeuronCount(); ++i ) {
    		setNeuron( i, n.getNeuron(i) );
    	}
	}
*/	  
	//  Network(const Network &n) : subPop(n){;}
	public void resetActivation() {
		Arrays.fill( activation, 0.0 );
//		for( int i = 0; i < activation.length; ++i ) {
//   		activation[i] = 0.0;
//    	}
	}
/*	  
	public void setNeuron( Neuron n, int position ) {
		setNeuron(position, n);
	}
	  
	public void setNetwork( Network n ) {
		fitness = n.fitness;
		for(int i = 0; i < getNeuronCount(); ++i) {
    		setNeuron(i, n.getNeuron(i));
    	}
	}
*/
	public void addFitness() {
		for(int i = 0; i < getNeuronCount(); ++i) {
    		getNeuron(i).fitness += fitness;
    	}
	}

	public void incrementTests() {
		for(int i = 0; i < getNeuronCount(); ++i) {
    		++getNeuron(i).tests;
    	}
	}
	  
	public int lesion( Environment e ) {
		int sp = -1;
		double lesionFitness, max = 0.0;
		fitness = e.evalNet(this);
		System.out.println("UNlesioned : " + fitness);
		for( int i = 0; i < getNeuronCount(); ++i) {
    		getNeuron(i).lesioned = true;
    		lesionFitness = e.evalNet(this);
    		System.out.println("lesion " + i + " : " + lesionFitness);
    		getNeuron(i).lesioned = false;
    		if(lesionFitness > max) {
    			max = lesionFitness;
    			sp = i;
    		}
		}
		if( max < (fitness * Network.LESION_THRESHOLD) ) {
			sp = -1;
    	}
    	return sp;
	}
	
	public int getInputCount() {
		return inputCount;
	}
	
	public int getOutputCount() {
		return outputCount;
	}
	
	public int getPopulationSize() {
		return neurons.length;
	}
	
//	public int getGeneCount() {
//		return geneCount;
//	}
	
	public void addNeuron() {
		Neuron[] temp = new Neuron[ neurons.length+1 ];
		System.arraycopy( neurons, 0, temp, 0, neurons.length );
		temp[neurons.length] = new Neuron(getGeneSize());
		neurons = temp;
	}
	
	public void removeNeuron( int index ) {
		Neuron[] temp = new Neuron[ neurons.length-1 ];
		System.arraycopy( neurons, 0, temp, 0, index );
		System.arraycopy( neurons, index+1, temp, index, neurons.length-index-1 );
		neurons = temp;
	}
/*	
	private void setNeuron( int index, Neuron n ) {
		neurons[index] = (Neuron)n.clone();
	}
*/	
	public Neuron getNeuron( int index ) {
		return neurons[index];
	}
	
	public int getNeuronCount() {
		return neurons.length;
	}

	public void create() { // creates a random subpopulation of neurons
		for (int i = 0; i < neurons.length; ++i) {
			neurons[i] = new Neuron(getGeneSize());
			neurons[i].create();
		}
	}
	
	public int load(String filename) {
		try {
			BufferedReader br = new BufferedReader( new InputStreamReader( new FileInputStream( filename ) ) ) ;
			StreamTokenizer st = new StreamTokenizer( br );
			st.nextToken();
			neurons = new Neuron[(Integer.valueOf(st.toString())).intValue()];
			int geneCount = (Integer.valueOf(st.toString())).intValue();
			for( int i = 0; i < neurons.length; ++i ) {
				neurons[i] = new Neuron(geneCount);
    			for( int j = 0; j < geneCount; ++j ) {
					st.nextToken();
    				neurons[i].weight[j] = Double.valueOf(st.toString()).doubleValue();
    			}
			}
			br.close();
		} catch( IOException e ) {
			System.out.println();
    		System.out.print(" Error - cannot open " + filename);
    		System.exit(1);
		}
		return neurons.length;
	}
	
	public void setRandomNetwork( SubPopulation[] pops ) {
		for( int j = 0; j < pops.length; ++j ) {
			neurons[j] = pops[j].selectNeuron();
//			setNeuron( j, pops[j].selectNeuron() );
		}
	}
	
	public void saveText(String fname) {
		try {
			PrintWriter pw = new PrintWriter( new FileOutputStream( fname ) );
			pw.println( neurons.length );
			pw.println( getGeneSize() );
			for( int i = 0; i < neurons.length; ++i ) {
    			for( int j = 0; j < getGeneSize(); ++j ) {
    				pw.print( neurons[i].weight[j] );
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
	
	public void addConnection(int locus) {
		// TODO: Make more efficient
    	for(int i=0; i < neurons.length; ++i) {
    		neurons[i].weight = doubleArrayInsert(neurons[i].weight, 1.0, locus);
		}
	}
	
//----------------------------------------------------------------------
// opposite of addConnection.
	public void removeConnection(int locus) {
		if( locus < neurons.length ) {
    		for( int i = 0; i < neurons.length; ++i ) {
    			neurons[i].weight = doubleArrayRemove( neurons[i].weight, locus );
    		}
    	}
	}	

	public abstract int getGeneSize();
	public abstract int getTotalInputs();	// returns sum of environment and recurrent inputs
	public abstract void activate( double[] input, double[] output);
	public abstract void save(String filename);
  //  public void print();
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