package simGWAS;

import java.util.Random;

public class RecombWorker {
	recombStructure recombs;//populated using user-defined recombination map
	int length = ParamParser.chromLength;
	double recombRate;
	
	private class recombStructure{
		int startBase;
		double rate;
		private double cumulative;
		recombStructure next;
		public recombStructure(int aStartBase,double aRate, double aCumulative,recombStructure aNext){
			startBase = aStartBase;
			rate = aRate;
			next = aNext;
		}
		public double getCumulative() {
			return cumulative;
		}
		public void setCumulative(double cumulative) {
			this.cumulative = cumulative;
		}
	}
	
	public void addRecombSiteLL(int aStart,double aRate){
		recombStructure  newRecomb = new recombStructure(aStart,aRate,0.0,null);
		recombStructure tempRecomb = recombs;
		if(recombs == null){
			recombs = newRecomb;
			return;
		}
		else if (recombs.startBase>aStart){
			newRecomb.next = recombs;
			recombs = newRecomb;
			return;			
		}
		else{
			while (tempRecomb.next !=null){
				if(tempRecomb.next.startBase > aStart){
					newRecomb.next = tempRecomb.next;
					tempRecomb.next = newRecomb; 
					return;
				}
				tempRecomb = tempRecomb.next;
			}
			newRecomb.next = null;
			tempRecomb.next = newRecomb;
		}		
	}

	public recombStructure getRecomb(){
		return recombs;
	}
	
	public void recomb_calc_r(){
		  double rr = 0;
		  recombStructure temprecomb = recombs;

		  while (temprecomb != null && temprecomb.startBase < length) {
		    if (temprecomb.next == null || 
			temprecomb.next.startBase > length) {
		      rr += (length - temprecomb.startBase) * 
			temprecomb.rate;
		    }
		    else
		      rr += (temprecomb.next.startBase - 
			     temprecomb.startBase) * temprecomb.rate;
				
		    temprecomb.setCumulative(rr);
		    temprecomb = temprecomb.next;

		  }
		  recombRate = rr;
	}
	
	public double pickRecombLoc(){
		  double temp;
		  double temp1;
		  recombStructure temprecomb = recombs;
		  double rr = 0;
		  double end;
				       
		  /* choosing location.. */
		  Random random = new Random();
		  temp1 = random.nextDouble();
		  temp = (double) (temp1 * (recombRate));

		  while (temprecomb != null && temprecomb.startBase < length && rr < temp) {
			    if (temprecomb.next == null || temprecomb.next.startBase > length) {
			      rr += (length - temprecomb.startBase) * temprecomb.rate;
			    }
			    else {
			      rr += (temprecomb.next.startBase - temprecomb.startBase) * temprecomb.rate;			
			    }
			    if (rr < temp) {
			      temprecomb = temprecomb.next;
			    }
		  }
		  if (temprecomb.next == null || temprecomb.next.startBase > length)
			  end = length;
		  else end = temprecomb.next.startBase;
		  return (double)((int) (end - (rr - temp) / temprecomb.rate)) / length;
	}
}
