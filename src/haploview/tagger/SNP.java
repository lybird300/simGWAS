package haploview.tagger;

//import java.util.Vector;
//import java.util.Comparator;
import java.util.HashSet;

public class SNP extends VariantSequence {
    private int position;
    private double MAF;
    private HashSet<SNP> LDList;
	//Alleles are represented by 1(mutated but not necessarily minor allele) and 2 (original but not necessarily major allele) in Dan's data set;
	//In this function, they are further converted to 2 and 3 respectively (as "C" and "G" in Haploview.tagger)
    private byte minor;
    private byte major;


    public SNP(int l,double maf,byte minor,byte major) {
        position =l;
        MAF = maf;
        LDList = new HashSet<SNP>();
        this.minor = minor;
        this.major = major;
    }  

    public byte getMinor(){
        return minor;
    }

    public byte getMajor(){
        return major;
    }

    public int compareTo(Object o) {
    	SNP s = (SNP)o;
        if(this.equals(s)) {
            return 0;
        } else if(this.position == s.position) {
            return 0;
        } else {
            return this.position > s.position ? 1 : -1;
        }
    }

    public boolean equals(Object o) {
        if (o instanceof SNP){
            SNP s = (SNP)o;
            if(position == s.position) {
                return true;
            }else {
                return false;
            }
        }
        return false;
    }

    public String getName() {
        return String.valueOf(position);
    }

    public double getMAF() {
        return MAF;
    }

    public int getPosition() {
        return position;
    }

    public void addToLDList(SNP s){
        LDList.add(s);
    }

    public HashSet<SNP> getLDList() {
        return LDList;
    }

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + position;
		return result;
	}
 }
