package beast.app;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.special.*;

/**
 * Truncated normal distribution, mostly from
 * http://en.wikipedia.org/wiki/Truncated_normal_distribution
 * However, that page has some errors, which have been corrected.
 * 
 * @author mao
 *
 */
public class TruncatedNormal extends AbstractRealDistribution {

	private static final long serialVersionUID = 3279850685379757823L;

	protected final double mu;
	protected final double sigma;	
	protected final double cdf_a;
	protected final double cdf_b;
	protected final double Z;
	
	protected final double alpha;
	protected final double beta;
	protected final double aa;
	protected final double bb;
	protected final org.apache.commons.math3.distribution.NormalDistribution norm;//=new org.apache.commons.math3.distribution.NormalDistribution(mean, sd);
	
    public static double sampleUpgraded(double m, double s, double l, double u)
    {
    	double x=sampleUpgraded((l-m)/s,(u-m)/s);
    	
    	double Z=(m+(s*x));
    	return Z;
    }
    
    public static double ntail(double l, double u)
    {
    	RandomGenerator rng=new Well19937c();
    	double c=(l*l)/2.0;
    	double doubleU=(u*u)/2.0;

    	double n=1.0;
    	double f= java.lang.Math.expm1(c-doubleU);
    	double x;//=c-Math.log(1+rng.nextDouble()*f);
    	double rndval;//=rng.nextDouble();
    	
    	do
    	{
    		x=c-Math.log(1+rng.nextDouble()*f);
    		rndval=rng.nextDouble();
    	}
    	while((rndval*rndval*x)>=c);
    	
    	double result=Math.sqrt(2*x);

    	return result;
    }
    
    public static double trnd(double l, double u)
    {
    	org.apache.commons.math3.distribution.NormalDistribution normG = new org.apache.commons.math3.distribution.NormalDistribution();
    	double x;
    	
    	do
    	{
    		x=normG.sample();
    	}while(x<l || x>u);
    	
    	return x;
    }
    
    public static double tn(double l, double u)
    {
    	RandomGenerator rng=new Well19937c();
    	double tol=2.0;
    	double result =0.0;
    	
    	if(Math.abs(u-l)>tol)
    	{
    		result=trnd(l,u);
    	}
    	else
    	{
    	    //tl=l(I); tu=u(I); pl=erfc(tl/sqrt(2))/2; pu=erfc(tu/sqrt(2))/2;
    	    double pl=org.apache.commons.math3.special.Erf.erfc(l/Math.sqrt(2));
    	    double pu=org.apache.commons.math3.special.Erf.erfc(u/Math.sqrt(2));
    	    result=Math.sqrt(2)*org.apache.commons.math3.special.Erf.erfcInv(2*(pl-(pl-pu)*rng.nextDouble()));
    	    //x(I)=sqrt(2)*erfcinv(2*(pl-(pl-pu).*rand(size(tl))));
    	}
    	
    	return result;
    }
    
    public static double sampleUpgraded(double l, double u)
    {
    	double a=.66; // threshold
    	double result=0.0;
    	
    	if(l>=a) // a<=l<u
    	{
    		result=ntail(l, u);
    	}
    	else
    	{
    		if(u<=-a) // l<u<=-a
    		{
        	    //tl=-u(J); tu=-l(J); x(J)=-ntail(tl,tu);
        	    result=-ntail(-u,-l);
    		}
    		else
    		{// otherwise
    			result=tn(l, u);
    		}
    	}
    	
    	return result;
    }

    public TruncatedNormal(RandomGenerator rng, double mean, double sd, double lb, double ub) {
		super(rng);
		
		if( mean == Double.NaN || sd == Double.NaN ||
				lb == Double.NaN || ub == Double.NaN ) 
			throw new IllegalArgumentException("Cannot take NaN as an argument");				
		
		norm=new org.apache.commons.math3.distribution.NormalDistribution(mean, sd);

		mu = mean;
		sigma = sd;
		cdf_a = norm.cumulativeProbability((lb - mu)/sigma);//NormalDist.cdf01((lb - mu)/sigma);
		cdf_b = norm.cumulativeProbability((ub - mu)/sigma);//NormalDist.cdf01((ub - mu)/sigma);
		Z = cdf_b - cdf_a;
		
		if(Z==0)
		{
			int dbg=0;
		}
		
		alpha = (lb - mu) / sigma;
		beta = (ub - mu) / sigma;
		this.aa = lb;
		this.bb = ub;
	}
	
	protected TruncatedNormal(double mean, double sd, double lb, double ub, double cdf_a, double cdf_b) {
		super(new Well19937c());
		
		norm=new org.apache.commons.math3.distribution.NormalDistribution(mean, sd);

		this.mu = mean;
		this.sigma = sd;
		this.cdf_a = cdf_a;
		this.cdf_b = cdf_b;
		this.Z = cdf_b - cdf_a;
		
		this.alpha = (lb - mu) / sigma;
		this.beta = (ub - mu) / sigma;
		this.aa = lb;
		this.bb = ub;
	}
	
	public TruncatedNormal(double mean, double sd, double lb, double ub) { 
		this(new Well19937c(), mean, sd, lb, ub);
	}
	
	@Override
	public double probability(double x) {		
		return 0;
	}

	@Override
	public double density(double x) {
		if( x <= aa || x >= bb ) return 0;

        return norm.density((x - mu)/sigma) / (sigma * Z);
		// return NormalDist.density01((x - mu)/sigma) / (sigma * Z);
	}

	@Override
	public double cumulativeProbability(double x) {
		if( x <= aa ) return 0;
		else if( x >= bb ) return 1;
			
		double u = norm.cumulativeProbability((x - mu)/sigma);//NormalDist.cdf01((x - mu)/sigma);
		return (u - cdf_a) / Z;
	}

	@Override
	public double inverseCumulativeProbability(double p) throws OutOfRangeException {		
		if (p < 0.0 || p > 1.0) throw new OutOfRangeException(p, 0, 1);        
		
		double val = p * Z + cdf_a;
		
		return mu + sigma * norm.inverseCumulativeProbability(val);//NormalDist.inverseF01(val);
	}

	@Override
	public double getNumericalMean() {
		double phi_a = norm.density(alpha);// NormalDist.density01(alpha);
		double phi_b = norm.density(beta); //NormalDist.density01(beta);
		
		return mu + (phi_a - phi_b) * sigma / Z;
	}

	@Override
	public double getNumericalVariance() {
		//double phi_a = NormalDist.density01(alpha);
		//double phi_b = NormalDist.density01(beta);
		double phi_a = norm.density(alpha);// NormalDist.density01(alpha);
		double phi_b = norm.density(beta); //NormalDist.density01(beta);
		double sq = ( phi_a - phi_b ) / Z;
		double br = 1 + ( alpha * phi_a - beta * phi_b ) / Z - sq * sq;  
		return sigma * sigma * br;
	}

	@Override
	public double sample() {			
		double val = super.random.nextDouble() * Z + cdf_a;
		double retVal = mu + sigma * norm.inverseCumulativeProbability(val);
		if(Double.isNaN(retVal) || (!Double.isFinite(retVal)))
		{
			int sullo=0;
		}
		return mu + sigma * norm.inverseCumulativeProbability(val); //NormalDist.inverseF01(val);
	}

	@Override
	public double getSupportLowerBound() {		
		return aa;
	}

	@Override
	public double getSupportUpperBound() {		
		return bb;
	}

	@Override
	public boolean isSupportLowerBoundInclusive() {		
		return false;
	}

	@Override
	public boolean isSupportUpperBoundInclusive() {		
		return false;
	}

	@Override
	public boolean isSupportConnected() {		
		return true;
	}

}
