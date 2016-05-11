//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: FeedForwardNetwork.java,v $
//////////////////////////////////////////////////////////////////////////
public class FeedForwardNetwork extends Network {
	public FeedForwardNetwork( int inputCount, int neuronCount, int outputCount ) {
		super( inputCount, neuronCount, outputCount );
	}
	
	public void activate( double[] input, double[] output) {
		// evaluate hidden/output layer 
		for( int i = 0; i < getNeuronCount(); ++i )
		{  //for each hidden unit
			Neuron n = getNeuron(i);

    			if( !n.lesioned )
			{
    				activation[i] = 0.0;
    				for( int j = 0; j < getInputCount(); ++j ) 
				{
					activation[i] += n.weight[j] * input[j];
	    			}
    				activation[i] = sigmoid( activation[i] );
    			}
		}

		for( int i = 0; i < output.length; ++i ) 
		{
    			output[i] = 0.0;
    			for( int j = 0; j < getNeuronCount(); ++j ) 
			{
    				output[i] += activation[j] * getNeuron(j).weight[getInputCount()+i];
    			}
    			output[i] = sigmoid( output[i] );
		}  
	}
	  
	public void save( String filename ) {
		saveText( filename + "_FF" );
	}

	public int getTotalInputs() {
		return getInputCount();
	}

	public int getGeneSize() {
		return getInputCount() + getOutputCount();
	}
}