package haploview.tagger;

//import java.util.Hashtable;

public class LocusCorrelation {
    Allele a;
    double rsq;
    double dprime;

    public LocusCorrelation(Allele a, double rsq, double dprime) {
        this.a = a;
        this.rsq = rsq;
        this.dprime = dprime;
    }

    public Allele getAllele(){
        return a;
    }

    public double getRsq(){
        return rsq;
    }
    
    public double getDPrime(){
        return dprime;
    }
}
