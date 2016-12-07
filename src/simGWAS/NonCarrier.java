package simGWAS;

public class NonCarrier {
	protected char[][] chroms;
	public NonCarrier(){
		chroms = new char[2][];
	}
	public NonCarrier(int[][] haps){
		this();
		String hap = "";
		for(int j = 0; j < haps[0].length; j++)
			hap += haps[0][j];
		chroms[0] = hap.toCharArray();
		hap = "";
		for(int j = 0; j < haps[1].length; j++)
			hap += haps[0][j];
		chroms[1] = hap.toCharArray();
	}
	public void recordHap(int chrID, String hap){
		this.chroms[chrID] = hap.toCharArray();
	}
	public char[] getSingleChrom(int id){
		if(id!=0 && id!=1){
			System.out.println("Wrong parameters.\n");
			return null;
		}
		return chroms[id];
	}
	public char[][] getBothChroms(){
		return chroms;
	}
}
