package simGWAS;

import java.security.SecureRandom;

public class RandomNum{
	final int IM1 = 2147483563;
	final int IM2 = 2147483399;
	final double AM = (1.0/IM1);
	final int IMM1 = (IM1-1);
	final int  IA1 = 40014;
	final int IA2  = 40692;
	final int IQ1 = 53668;
	final int IQ2 = 52774;
	final int IR1 = 12211;
	final int IR2 = 3791;
	final int  NTAB = 32;
	final double  NDIV = (1+(double)IMM1/(double)NTAB);
	final  double EPS = 1.2e-12;
	final  double  RNMX = (1.0-EPS);
	long randseed;
	long idum;
	SecureRandom generator;
	 long idum2,iy;
	 long[] iv;
	 
	 public RandomNum(){
		 idum2 = 123456789;
		 iy = 0;
		 iv = new long[NTAB];
		 generator = new SecureRandom();
	 }
	 public synchronized double randomDouble(){
		 return ran2();
	 }
	 public  void setRngSeed(long rseed){
		 randseed = rseed;
	 }
	 public  long getRandSeed(){
		return randseed;
	 }
	 public  long seedRNG(){
		seedRandom(generator.nextInt());
		return randseed;
	 }
	 public  void seedRandom(long aSeed){
		long newSeed = aSeed;
		newSeed = Math.abs(newSeed);
		randseed = newSeed;
		newSeed *= (-1);
		idum = newSeed;
		ran2();	
	}
	double ran2(){
	        int j;
	        long k;
	        double temp;
	
	        if (idum <= 0) {
	                if (-(idum) < 1) idum=1;
	                else idum = -(idum);
	                idum2=(idum);
	                for (j=NTAB+7;j>=0;j--) {
	                        k=(idum)/IQ1;
	                        idum=IA1*(idum-k*IQ1)-k*IR1;
	                        if (idum < 0) idum += IM1;
	                        if (j < NTAB) iv[j] = idum;
	                }
	                iy=iv[0];
	        }
	        k=(idum)/IQ1;
	        idum=IA1*(idum-k*IQ1)-k*IR1;
	        if (idum < 0) idum += IM1;
	        k=idum2/IQ2;
	        idum2=IA2*(idum2-k*IQ2)-k*IR2;
	        if (idum2 < 0) idum2 += IM2;
	        j=(int) ((double)iy/NDIV);
	        if (j == NTAB) j--;
	        iy=iv[j]-idum2;
	        iv[j] = idum;
	        if (iy < 1) iy += IMM1;
	        temp = AM*iy;
	        if (temp > RNMX) return RNMX;
	        else return temp;
	}
}
