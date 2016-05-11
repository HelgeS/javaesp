//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: Neuron.java,v $
//////////////////////////////////////////////////////////////////////////
public class Neuron implements Cloneable {
	public double[] weight;
//	public double[] delta;
	public double fitness;               //neuron's fitness value
	public int tests;                 //holds the total # of networks participated in
	public double activation;         //neuron's activation level 
	public boolean lesioned = false;
	
	//constructor
	public Neuron(int numNeurons) {
		weight = new double[numNeurons];	//fill weight vector with zeros
	}
	
	public Object clone() {
		try {
			Neuron rc = (Neuron)super.clone();
			rc.weight = (double[])rc.weight.clone();
			return rc;
		} catch( CloneNotSupportedException e ) {
			return null;
		}
	}
/*	
	public Neuron(final Neuron n) {
		weight = new double[n.weight.length];	
		for( int i=0; i < weight.length; i++ ) {
			weight[i] = n.weight[i];
		}
//		delta = new double[n.delta.length];	
//		for( int i=0; i < delta.length; i++ ) {
//			delta[i] = n.delta[i];
//		}
		fitness = n.fitness;
		tests = n.tests;
		activation = n.activation;
		lesioned = n.lesioned;
	}
*/	
	public void create() {
		for( int i = 0;i < weight.length; ++i ) {
    		weight[i] = (RandomSingleton.getInstance().nextDouble() * 12.0) - 6.0;
    	}
	}
	
	RandomCauchy rndCauchy = new RandomCauchy();
	public void perturb( Neuron n ) {
		perturb( n, rndCauchy, 0.3 );
	}
	
	public void perturb(Neuron n, RandomFunction randFn, double coeff) {
		for( int i = 0; i < weight.length; ++i ) {
    		weight[i] = n.weight[i] + randFn.randomFunction(coeff);
    	}
	}
}