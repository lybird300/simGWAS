package haploview.tagger;


public class PairwiseLinkage{

    private double dprime, lod, r2/*, ci_low, ci_high*/;
    //private double[] freqs;

    PairwiseLinkage(double d, double l, double r/*, double lo, double hi, double[] f*/){
        dprime = d;
        lod = l;
        r2 = r;
        /*ci_low = lo;
        ci_high = hi;
        freqs = f;*/
    }

    public double getDPrime(){
        return dprime;
    }

    public double getLOD(){
        return lod;
    }

    public double getRSquared(){
        return r2;
    }

    /*public double getConfidenceLow(){
        return ci_low;
    }

    public double getConfidenceHigh(){
        return ci_high;
    }

    public double[] getFreqs(){
        return freqs;
    }

    public String toString(){
        return new String(dprime + "\t" + lod + "\t" + r2 + "\t" + ci_low + "\t" + ci_high);
    }*/
}
