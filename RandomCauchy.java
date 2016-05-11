//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: RandomCauchy.java,v $
//////////////////////////////////////////////////////////////////////////
public class RandomCauchy implements RandomFunction {
    public double randomFunction( double wtrange ) {
		double u = 0.5, Cauchy_cut = 10.0;
		 
		while (u == 0.5) {
			u = RandomSingleton.getInstance().nextDouble();
		}
		u = wtrange * Math.tan(u * Math.PI);
		if(Math.abs(u) > Cauchy_cut) {
			return randomFunction(wtrange);
		} else {
			return u;
		}
    }

}