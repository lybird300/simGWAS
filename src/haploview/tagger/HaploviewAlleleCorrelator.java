package haploview.tagger;

//import java.util.ArrayList;
import java.util.Hashtable;
//import java.util.List;
//import java.util.HashSet;
//import java.util.Iterator;

public class HaploviewAlleleCorrelator{
    private Hashtable<VariantSequence, Integer> indicesByVarSeq;
    private DPrimeTable dpTable;
    private HaploData theData;
    //private Haplotype[] phasedCache;
    //private Hashtable<VariantSequence, Integer> phasedCacheIndicesByVarSeq;
    //private Hashtable<Comparison, LocusCorrelation> lcByComparison;

    public HaploviewAlleleCorrelator(Hashtable<VariantSequence, Integer> indicesByVarSeq2, HaploData hd) {
        this.indicesByVarSeq = indicesByVarSeq2;
        dpTable = HaploData.dpTable;
        theData = hd;
        //lcByComparison = new Hashtable<Comparison, LocusCorrelation>();
    }

    public LocusCorrelation getCorrelation(VariantSequence v1, VariantSequence v2) {
        if(v1.getName().equalsIgnoreCase(v2.getName())) {
            return new LocusCorrelation(null,1,1);
        }

        //if (v1 instanceof SNP && v2 instanceof SNP){
            //we are comparing two snps
            int v1Index = (indicesByVarSeq.get(v1)).intValue();
            int v2Index = (indicesByVarSeq.get(v2)).intValue();
            double rsq, dprime;
            PairwiseLinkage pl;
            if (v1Index > v2Index){
                pl = dpTable.getLDStats(v2Index,v1Index);
            }else{
                pl = dpTable.getLDStats(v1Index,v2Index);
            }

            if(pl == null) {
                rsq = 0;
                dprime=0;
            } else {
                rsq = pl.getRSquared();
                dprime = pl.getDPrime();
            }
			LocusCorrelation lc = new LocusCorrelation(null,rsq,dprime);
            return lc;
        //}
        /*else{
            //we are comparing a snp vs. a block
            SNP theSNP;
            Block theBlock;
            if (v1 instanceof SNP){
                theSNP = (SNP) v1;
                theBlock = (Block) v2;
            }else{
                theSNP = (SNP) v2;
                theBlock = (Block) v1;
            }

            Comparison c = new Comparison(theSNP, theBlock);
            if (lcByComparison.containsKey(c)){
                return (LocusCorrelation) lcByComparison.get(c);
            }

            
            Allele curBestAllele = null;
            double curBestRsq = 0;
            if(HaploData.isPhasedData()){
                double[] rsquareds;
                int snpid = (indicesByVarSeq.get(theSNP)).intValue();
                int[] block = null;
                if(theBlock.getMarkerCount() == 2){
                    block = new int[2];
                    block[0] = (indicesByVarSeq.get(theBlock.getSNP(0))).intValue();
                    block[1] = (indicesByVarSeq.get(theBlock.getSNP(1))).intValue();
                    HaploData.RSquared rsquaredWrapper = theData.getPhasedRSquared(snpid, block);
                    rsquareds = rsquaredWrapper.getRsquareds();

                    for(int i=0;i<rsquareds.length;i++){
                        if(rsquareds[i] > curBestRsq) {
                            curBestRsq = rsquareds[i];
                            String allele = "";

                            if(i < 2){
                                allele += SNP.getMajor();
                            }else{
                                allele += SNP.getMinor();
                            }
                            if(i % 2 == 0){
                                allele += SNP.getMajor();
                            }else{
                                allele += SNP.getMinor();
                            }
                            curBestAllele = new Allele(theBlock,allele);
                        }

                    }
                }else if (theBlock.getMarkerCount() == 3){
                    block = new int[3];
                    block[0] = (indicesByVarSeq.get(theBlock.getSNP(0))).intValue();
                    block[1] = (indicesByVarSeq.get(theBlock.getSNP(1))).intValue();
                    block[2] = (indicesByVarSeq.get(theBlock.getSNP(2))).intValue();
                    HaploData.RSquared rsquaredWrapper = theData.getPhasedRSquared(snpid, block);
                    rsquareds = rsquaredWrapper.getRsquareds();

                    for(int i=0;i<rsquareds.length;i++){
                        if(rsquareds[i] > curBestRsq) {
                            curBestRsq = rsquareds[i];
                            String allele = "";
                            
                            if(i < 4 ){
                                allele += SNP.getMajor();
                            }else{
                                allele += SNP.getMinor();
                            }
                            if(i % 4 < 2){
                                allele += SNP.getMajor();
                            }else{
                                allele += SNP.getMinor();
                            }
                            if(i % 2 == 0){
                                allele += SNP.getMajor();
                            }else{
                                allele += SNP.getMinor();
                            }

                            curBestAllele = new Allele(theBlock,allele);
                        }
                    }
                }

            }else{
                int[][] genos = new int[phasedCache.length][theBlock.getMarkerCount()+1];
                for (int i = 0; i < phasedCache.length; i++){
                    //create a temporary set of mini hap genotypes with theSNP as the first marker and theBlock's markers as the rest
                    genos[i][0] = phasedCache[i].getGeno()[phasedCacheIndicesByVarSeq.get(theSNP).intValue()];
                    for (int j = 1; j < theBlock.getMarkerCount()+1; j++){
                        genos[i][j] = phasedCache[i].getGeno()[phasedCacheIndicesByVarSeq.get((theBlock.getSNP(j-1))).intValue()];
                    }
                }
                for (int i = 0; i < genos.length; i++){
                    double aa=0,ab=0,bb=0,ba=0;
                    for (int j = 0; j < genos.length; j++){
                        if (genos[j][0] == genos[0][0]){
                            if(!sameHap(genos[i], genos[j])){
                                ab += phasedCache[j].getPercentage();
                            }else{
                                aa += phasedCache[j].getPercentage();
                            }
                        }else{
                            if(!sameHap(genos[i], genos[j])){
                                bb += phasedCache[j].getPercentage();
                            }else{
                                ba += phasedCache[j].getPercentage();
                            }
                        }
                    }
                    //p is snp's freq, q is hap's freq
                    double p = aa+ab;
                    double q = ba+aa;
                    //round to 5 decimal places.
                    double rsq = Util.roundDouble(Math.pow((aa*bb - ab*ba),2)/(p*(1-p)*q*(1-q)),3);
                    if (rsq > curBestRsq){
                        StringBuffer sb = new StringBuffer();
                        for (int j = 1; j < genos[i].length; j++){
                            sb.append(genos[i][j]);
                        }
                        curBestAllele = new Allele(theBlock,sb.toString());
                        curBestRsq = rsq;
                    }
                }

            }
            LocusCorrelation lc = new LocusCorrelation(curBestAllele, curBestRsq);
            lcByComparison.put(new Comparison(theSNP, theBlock),lc);

            return (lc);
        }*/
    }

    /*private boolean sameHap(int[] a, int[] b) {
        if(a == null || b == null) {
            throw new NullPointerException("blah");
        }
        if(a.length != b.length) {
            return false;
        }
        for(int i=1;i<a.length;i++) {
            if(a[i] != b[i]) {
                return false;
            }
        }
        return true;
    }*/

    /*public void phaseAndCache(HashSet<SNP> snpList){
        phasedCacheIndicesByVarSeq = new Hashtable<VariantSequence, Integer>();
        if(!HaploData.isPhasedData()){
                //create a pseudo block of the SNP to be correlated and the markers from the multi-test
        	Integer[] blockArray = new Integer[snpList.size()];
            Iterator<SNP> itr = snpList.iterator();
            for (int i = 0; i < blockArray.length; i++){
                SNP n = itr.next();
                blockArray[i] = indicesByVarSeq.get(n).intValue();
                phasedCacheIndicesByVarSeq.put(n, new Integer(i));
            }
            List<Integer[]> victor = new ArrayList<Integer[]>();
            victor.add(blockArray);
            phasedCache = theData.generateHaplotypes(victor, true)[0];
        }
    }*/

    private class Comparison{
        private SNP s;
        private Block b;

        private Comparison(SNP s, Block b){
            this.s = s;
            this.b = b;
        }

        public boolean equals(Object o){
            Comparison c = (Comparison) o;
            return (s.equals(c.s) && b.equals(c.b));
        }

        public int hashCode(){
            return s.hashCode() + b.hashCode();
        }
    }

}
