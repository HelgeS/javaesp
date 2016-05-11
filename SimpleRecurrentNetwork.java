//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: SimpleRecurrentNetwork.java,v $
//////////////////////////////////////////////////////////////////////////
public class SimpleRecurrentNetwork extends Network {
	public SimpleRecurrentNetwork( int inputCount, int neuronCount, int outputCount ) {
		super( inputCount, neuronCount, outputCount );
	}
/*	
	public SimpleRecurrentNetwork( SimpleRecurrentNetwork net ) {
		super( net );
	}
*/	
/*	
	public void addNeuron() {
		addConnection( getTotalInputs()-1 );
		super.addNeuron();
	}
 
	public void removeNeuron(int sp) {
		removeConnection( getInputCount()+sp );
		super.removeNeuron(sp);
	}
*/	  
	private double[] tmp;

	public void activate( double[] input, double[] output) {
		/* evaluate hidden/output layer */
		tmp = new double[getInputCount() + activation.length];
		System.arraycopy( input, 0, tmp, 0, getInputCount() );
		System.arraycopy( activation, 0, tmp, getInputCount(), activation.length );

		//System.arraycopy( activation, 0, input, getInputCount(), activation.length );

		for( int i = 0; i < activation.length; ++i ) 
		{  /*for each hidden unit*/
    			activation[i] = 0.0;
    			if(!getNeuron(i).lesioned)
			{
    				for( int j = 0; j < getTotalInputs(); ++j ) 
				{
					activation[i] += getNeuron(i).weight[j] * tmp[j];
			
    				}
    				activation[i] = sigmoid( activation[i] ); 
    			}
		}

		for( int i = 0; i < output.length; ++i ) 
		{
    			output[i] = 0.0;
    			for( int j = 0; j < activation.length; ++j ) 
			{
    				output[i] += activation[j] * getNeuron(j).weight[getTotalInputs()+i];
    			}
    			output[i] = sigmoid( output[i] );	
		}  
	}
	  
	public void save(String filename) {
		saveText( filename + "_SRN" );
	}

	public int getTotalInputs() {
		return getInputCount() + getNeuronCount();
	}

	public int getGeneSize() {
		return getTotalInputs() + getOutputCount();
	}
}