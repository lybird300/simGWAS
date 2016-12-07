package simGWAS;

public class Subject {
	private String ID;
	private int chromA;
	private int chromB;
	private int carrierOrNot;//0-non-carrier, 1-carrier
	public Subject(String id, Integer chrom1, Integer chrom2, int carryStatus) {
		this.ID = id;
		this.chromA = chrom1.intValue();
		this.chromB = chrom2.intValue();
		carrierOrNot = carryStatus;
	}
	public int isCarrier() {
		return carrierOrNot;
	}
	public void setCarrier(int carried) {
		this.carrierOrNot = carried;
	}
	public String getID() {
		return ID;
	}
	public void setID(String id) {
		ID = id;
	}
	public int getChromA() {
		return chromA;
	}
	public void setChromA(int chromA) {
		this.chromA = chromA;
	}
	public int getChromB() {
		return chromB;
	}
	public void setChromB(int chromB) {
		this.chromB = chromB;
	}
}

