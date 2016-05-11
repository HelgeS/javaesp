//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: SecondOrderRecurrentNework.java,v $
//////////////////////////////////////////////////////////////////////////
public class SecondOrderRecurrentNetwork extends Network {
	private int outOffset;
	private double[][] newWeights;
	
	public SecondOrderRecurrentNetwork( int envInputCount, int neuronCount, int outputCount ) {
		super( envInputCount, neuronCount, outputCount );
		outOffset = getInputCount() * neuronCount;
		newWeights = new double[neuronCount][outOffset];
	}
	
	public Object clone() {
		SecondOrderRecurrentNetwork rc = (SecondOrderRecurrentNetwork)super.clone();
		rc.newWeights = (double[][])rc.newWeights.clone();
		for( int i=0; i < rc.newWeights.length; i++ ) {
			rc.newWeights[i] = (double[])rc.newWeights[i].clone();
		}
		return rc;
	}
/*	
	public SecondOrderRecurrentNetwork( SecondOrderRecurrentNetwork net ) {
		super( net );
		outOffset = net.outOffset;
		newWeights = new double[net.getPopulationSize()][outOffset];
		for( int i = 0; i < newWeights.length; i++ ) {
			for( int j = 0; j < newWeights[i].length; j++ ) {
				newWeights[i][j] = net.newWeights[i][j];
			}
		}
	}
*/	
	public void addNeuron() {
		for( int i=1; i < getInputCount()+1; ++i )  // add connection to neurons in spops
    		addConnection( getNeuronCount()*i + i - 1 );
		super.addNeuron();
	}
 
	public void removeNeuron(int sp) {
		for( int i=1; i < getInputCount()+1; ++i ) {  // remove connection to neurons in spops
    		removeConnection(getNeuronCount()*i);
		}
		super.removeNeuron(sp);
	}
	  
	public void activate( double[] input, double[] output) {
		// calculate new weights wij 
		for( int i = 0; i < getNeuronCount(); ++i )
    		if(!getNeuron(i).lesioned) {
			//Modification by Jacob Schrum, August 2007
    			for( int j = 0; j < getInputCount(); ++j ) {
					(newWeights[i])[j] = 0.0;
					for( int k = 0; k < activation.length; ++k ) { 
						(newWeights[i])[j] += activation[k] * getNeuron(i).weight[j*activation.length+k];
					}
    			}    
    		}
		// activate hidden layer
		for( int i = 0; i < activation.length; ++i ) {
    		activation[i] = 0.0;
    		if(!getNeuron(i).lesioned) {
    			for( int j = 0; j < getInputCount(); ++j ) {
					activation[i] += (newWeights[i])[j] * input[j];
					//printf("%f\n", activation[i]);
    			}
    			activation[i] = sigmoid( activation[i] );//, fabs(getNeuron(i).weight[0]/6.0));
    		}
		}
		// for( i = 0; i < getOutputCount();++i) 
		//output[i] = activation[i];
		for( int i = 0; i < getOutputCount();++i) {
    		output[i] = 0.0;
    		for( int j = 0; j < activation.length;++j) {
    			output[i] += activation[j] * getNeuron(j).weight[outOffset+i];
    		}
    		output[i] = sigmoid( output[i] );
    		//printf("%f\n", output[0]); //.weight[j]);
		} 
	}
	  
	public void save(String filename) {
		saveText( filename + "_2ndOrder" );
	}

	public int getTotalInputs() {
		return getInputCount();
	}

	public int getGeneSize() {
		return getInputCount() * getNeuronCount() + getOutputCount();
	}
}