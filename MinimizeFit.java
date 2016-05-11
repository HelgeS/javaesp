//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: MinimizeFit.java,v $
//////////////////////////////////////////////////////////////////////////
import java.util.Comparator;

public class MinimizeFit implements Comparator {
    public int compare(Object x, Object y) {
    	if( x instanceof Neuron && y instanceof Neuron ) {
    		if( ((Neuron)x).fitness < ((Neuron)y).fitness ) {
    			return -1;
    		} else if ( ((Neuron)x).fitness > ((Neuron)y).fitness ) {
    			return 1;
    		} else {
        		return 0;
        	}
        } else {
        	throw new ClassCastException();
        }
    }
    
    public boolean equals( Object obj ) {
    	return this.getClass() == obj.getClass();
    }
}