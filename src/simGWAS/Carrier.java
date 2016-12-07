package simGWAS;

public class Carrier extends NonCarrier {
	public Carrier(){
		super();
	}
	public Carrier(int[][] haps, boolean chr0WithCV, boolean chr1WithCV){
		super();
		String hap = "";
		for(int i = 0; i < haps[0].length; i++)
			hap += haps[0][i];
		this.recordHap(0, ((chr0WithCV)?'1':'0') + hap);
		hap = "";
		for(int i = 0; i < haps[1].length; i++)
			hap += haps[1][i];
		this.recordHap(1, ((chr1WithCV)?'1':'0')+hap);
	}
	public boolean isCarrier() {
		return (chroms[0][0]=='1' || chroms[1][0]=='1');
	}
	public boolean doesChromWithCV(int id) {
		if(id!=0 && id!=1){
			System.out.println("Wrong parameters.\n");
			return false;
		}
		else if(id==0) return (chroms[0][0]=='1');
		else return (chroms[1][0]=='1');
	}
	public boolean isHomozygoteWithCV() {
		return (chroms[0][0]=='1' && chroms[1][0]=='1');
	}
}
