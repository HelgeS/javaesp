//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: Environment.java,v $
//////////////////////////////////////////////////////////////////////////
public interface Environment {
/*
	public int inputSize;
	public int outputSize;
	public double maxFitness;
*/	
	public double evalNet( Network net );
	public void nextTask();
	public void simplifyTask();
	public int getInputSize();
	public int getOutputSize();
	public double getMaxFitness();
	public void setupInput(double[] input);
}