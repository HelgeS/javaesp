//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: RandomSingleton.java,v $
//////////////////////////////////////////////////////////////////////////
import java.util.Random;

class RandomSingleton {
	private static Random m_random = new Random();
	
	public static Random getInstance() {
		return m_random;
	}
}

