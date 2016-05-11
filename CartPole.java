//////////////////////////////////////////////////////////////////////
//   Alan Oursland  
// 
//   ESP Java implementation.
//
//	 This code is a direct port of Faustino Gomez's ESP code
//   that was code started from Daniel Polani's translation of 
//   Dave Moriarty's original SANE code from C to C++.
//   
//   $Log: CartPole.java,v $
//////////////////////////////////////////////////////////////////////////
public class CartPole implements Environment {
	public static final double MUP      =        0.000002;
	public static final double MUC      =        0.0005;
	public static final double GRAVITY  = -9.8;
	public static final double MASSCART = 1.0;
	public static final double MASSPOLE_1 = 0.1;
	public static		double MASSPOLE_2 = 0.01;
	public static final double LENGTH_1 = 0.5;		  /* actually half the pole's length */
	public static		double LENGTH_2 = 0.05;
	public static final double FORCE_MAG = 10.0;
	public static final double TAU = 0.01;		  //seconds between state updates 

	public static final double one_degree = 0.0174532;	/* 2pi/360 */
	public static final double six_degrees = 0.1047192;
	public static final double twelve_degrees = 0.2094384;
	public static final double fifteen_degrees = 0.2617993;
	public static final double thirty_six_degrees = 0.628329;
	public static final double fifty_degrees = 0.87266;

	public static final boolean MARKOV = false; // markov (full state) of non-markov (no-velocities).

	private static final int inputSize = 7;
	private static final int outputSize = 1;
	private static final double maxFitness = 1000;
	private static final boolean RK4 = true;
	private static final double EULER_TAU = (TAU/4);
	
	public static final double MIN_INC = 0.001;
	private double state[] = new double[6];
	private double POLE_INC = 0.05;
	private double MASS_INC = 0.01;
	
	public CartPole( int index ) {
	}
	
	private boolean first_time = true;
	public void init( boolean randomize )
	{
		/*if (randomize) {
    		state[0] = (lrand48()%4800)/1000.0 - 2.4;
    		state[1] = (lrand48()%2000)/1000.0 - 1;
    		state[2] = (lrand48()%400)/1000.0 - 0.2;
    		state[3] = (lrand48()%400)/1000.0 - 0.2;
    		state[4] = (lrand48()%3000)/1000.0 - 1.5;
    		state[5] = (lrand48()%3000)/1000.0 - 1.5;
		}
		else {*/
		state[0] = state[1] = state[3] = state[4] = state[5] = 0;
		state[2] = 0.07; // one_degree;
    		//}
		if(first_time){
    		System.out.println("Initial Long pole angle = " + state[2]);
    		System.out.println("Initial Short pole length = " + LENGTH_2);
    		first_time = false;
		}
	}


	private static final double ML_1 = LENGTH_1 * MASSPOLE_1;
	private static double ML_2 = LENGTH_2 * MASSPOLE_2;
	void step( double action, double[] st, double[] derivs ) {
    	double force =  (action - 0.5) * FORCE_MAG * 2;
    	double costheta_1 = Math.cos(st[2]);
    	double sintheta_1 = Math.sin(st[2]);
    	double gsintheta_1 = GRAVITY * sintheta_1;
    	double costheta_2 = Math.cos(st[4]);
    	double sintheta_2 = Math.sin(st[4]);
    	double gsintheta_2 = GRAVITY * sintheta_2;
	    
    	double temp_1 = MUP * st[3] / ML_1;
    	double temp_2 = MUP * st[5] / ML_2;
    	double fi_1 = (ML_1 * st[3] * st[3] * sintheta_1) + (0.75 * MASSPOLE_1 * costheta_1 * (temp_1 + gsintheta_1));
    	double fi_2 = (ML_2 * st[5] * st[5] * sintheta_2) +	(0.75 * MASSPOLE_2 * costheta_2 * (temp_2 + gsintheta_2));
    	double mi_1 = MASSPOLE_1 * (1 - (0.75 * Math.pow(costheta_1,2)));
    	double mi_2 = MASSPOLE_2 * (1 - (0.75 * Math.pow(costheta_2,2)));
	    
	    
    	derivs[1] = (force + fi_1 + fi_2) / (mi_1 + mi_2 + MASSCART);
	    
    	derivs[3] = -0.75 * (derivs[1] * costheta_1 + gsintheta_1 + temp_1) / LENGTH_1;
    	derivs[5] = -0.75 * (derivs[1] * costheta_2 + gsintheta_2 + temp_2) / LENGTH_2;
	}
	
	void performAction( double[] output ) {
		int i;
		double[]  dydx = new double[6];
		 
		/*random start state for long pole*/
		/*state[2]= drand48();   */
		  
		    
		/*--- Apply action to the simulated cart-pole ---*/

		if(RK4) {
    		for(i=0;i<2;++i){
    		dydx[0] = state[1];
    		dydx[2] = state[3];
    		dydx[4] = state[5];
    		step(output[0],state,dydx);
    		rk4(output[0],state,dydx,state);
    		}
		} else {
    		for(i=0;i<8;++i){
    			step(output[0],state,dydx);
    			state[0] += EULER_TAU * dydx[0];
    			state[1] += EULER_TAU * dydx[1];
    			state[2] += EULER_TAU * dydx[2];
    			state[3] += EULER_TAU * dydx[3];
    			state[4] += EULER_TAU * dydx[4];
    			state[5] += EULER_TAU * dydx[5];
    		}
		}
	}
	
	void rk4( double f, double[] y, double[] dydx, double[] yout ) {
		int i;

		double hh;
		double h6;
		double[] dym = new double[6];
		double[] dyt = new double[6];
		double[] yt = new double[6];


		hh=TAU*0.5;
		h6=TAU/6.0;
		for (i=0;i<=5;i++) yt[i]=y[i]+hh*dydx[i];
		step(f,yt,dyt);
		dyt[0] = yt[1];
		dyt[2] = yt[3];
		dyt[4] = yt[5];
		for (i=0;i<=5;i++) yt[i]=y[i]+hh*dyt[i];
		step(f,yt,dym);
		dym[0] = yt[1];
		dym[2] = yt[3];
		dym[4] = yt[5];
		for (i=0;i<=5;i++) {
			yt[i]=y[i]+TAU*dym[i];
			dym[i] += dyt[i];
		}
		step(f,yt,dyt);
		dyt[0] = yt[1];
		dyt[2] = yt[3];
		dyt[4] = yt[5];
		for (i=0;i<=5;i++)
			yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	}
	
	boolean outsideBounds() {
		final double failureAngle = thirty_six_degrees; 

		return 
    		state[0] < -2.4              || 
    		state[0] > 2.4               || 
    		state[2] < -failureAngle     ||
    		state[2] > failureAngle      ||
    		state[4] < -failureAngle     ||
    		state[4] > failureAngle;  
	}

	public double evalNet( Network net ) {
		int steps=0;
		double[] input = new double[net.getTotalInputs()];
		double[] output = new double[getOutputSize()];
		  
		//init(randomize);		// restart at some point


		init(false);

		net.resetActivation();
		while (steps++ < maxFitness) {
    		setupInput(input);
    		net.activate(input, output);
    		performAction(output);
    		if (outsideBounds())	// if failure
    		break;			// stop it now
		}
		return (double) steps;
	}
	
	public void nextTask() {
		LENGTH_2 += POLE_INC;   /* LENGTH_2 * INCREASE;   */
		MASSPOLE_2 += MASS_INC; /* MASSPOLE_2 * INCREASE; */
		ML_2 = LENGTH_2 * MASSPOLE_2;
		//  ++new_task;
		System.out.println("#Pole Length "+LENGTH_2);
	}
	
	public void simplifyTask() {
		if(POLE_INC > MIN_INC) {
    		POLE_INC = POLE_INC/2;
    		MASS_INC = MASS_INC/2;
    		LENGTH_2 -= POLE_INC;
    		MASSPOLE_2 -= MASS_INC;
			ML_2 = LENGTH_2 * MASSPOLE_2;
			System.out.println("#SIMPLIFY");
    		System.out.println("#Pole Length "+LENGTH_2);
		} else {
    		System.out.println("#NO TASK CHANGE");
    	}
	}
	
	public int getInputSize() {
		return inputSize;
	}
	
	public int getOutputSize() {
		return outputSize;
	}
	
	public double getMaxFitness() {
		return maxFitness;
	}
	
	public void setupInput(double[] input) {
		if(MARKOV){      
    		input[0] = state[0] / 4.8;
    		input[1] = state[1] /2;
    		input[2] = state[2]  / 0.52;
    		input[3] = state[3] /2;
    		input[4] = state[4] / 0.52;
    		input[5] = state[5] /2;
    		input[6] = .5;
		}
		else {
    		input[0] = state[0] / 4.8;
    		input[1] = 0.0;
    		input[2] = state[2]  / 0.52;
    		input[3] = 0.0;
    		input[4] = state[4] / 0.52;
    		input[5] = 0.0;
    		input[6] = .5;
		    
		}
	}

}