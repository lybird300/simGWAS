package haploview.tagger;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;


public class TaggerController {
    private Tagger tagger;
    private ArrayList<TagSequence> results;
    private boolean taggingCompleted = false;
    public List<TagThread> threadList = new ArrayList<TagThread>();
    //private Hashtable<String, SNP> snpHash;//Hashtable is synchronized and thread-safe (comparing to HashMap)

    public TaggerController(HaploData hd, List<SNP> allSNPs, int aggressionLevel, int maxDist, int maxNumTags, boolean findTags){
        List<SNP> taggerSNPs = new ArrayList<SNP>();
        
        Hashtable<VariantSequence, Integer> indicesByVarSeq = new Hashtable<VariantSequence, Integer>();

        for(int i=0;i<allSNPs.size();i++) {
            SNP s = allSNPs.get(i);
            taggerSNPs.add(s);
            indicesByVarSeq.put(s, new Integer(i));
        }

        for (int i = 0; i < allSNPs.size(); i++){
        	SNP taggerSNP = taggerSNPs.get(i);
            for (int j = 1; j < HaploData.dpTable.getLength(i); j++){
                PairwiseLinkage pl = HaploData.dpTable.getLDStats(i,j+i);
                if (pl != null && pl.getLOD() >= Tagger.DEFAULT_LOD_CUTOFF){
                	if(j+i < allSNPs.size()){//equivalent to: if (indicesByVarSeq.containsValue(new Integer(j+i)))
                        SNP ldsnp = (SNP) taggerSNPs.get(j+i);
                        taggerSNP.addToLDList(ldsnp);
                        ldsnp.addToLDList(taggerSNP);
                    }
                }
            }
        }

        HaploviewAlleleCorrelator hac = new HaploviewAlleleCorrelator(indicesByVarSeq, hd);
        tagger = new Tagger(taggerSNPs, taggerSNPs, hac, maxDist, maxNumTags);//change default values at Tagger.java L53
    }

    public void runTagger() {
    	/*while(Thread.activeCount() < 32 && taggingCompleted==false){
    		TagThread tagThread = new TagThread(tagger);   		
            taggingCompleted = false;
    		tagThread.start();
    		threadList.add(tagThread);
    	}*/
    	results = tagger.findTags();
    	taggingFinished();
    }

    public int getTaggedSoFar() {
        return tagger.taggedSoFar;
    }

    public int getUntaggableCount() {
        return tagger.getUntaggableCount();
    }

    public List<SNP> getForceIncludeds(){
        return tagger.getForceInclude();
    }


    private void taggingFinished() {
        taggingCompleted = true;
    }

    public boolean isTaggingCompleted() {
        return taggingCompleted;
    }

    public ArrayList<TagSequence> getResults() {
        return results;
    }

    public int getNumTagSNPs(){
        return tagger.getTagSNPs().size();
    }


    public class TagThread extends Thread{
        Tagger tagger;
        public TagThread(Tagger t) {
            tagger = t;
        }
        public void run() {
            results = tagger.findTags();
            taggingFinished();
        }
    }

    public double getMeanRSq(){
        return tagger.getMeanRSq();
    }

    public int getPercentCaptured(){
        return tagger.getPercentCaptured();
    }
    
    public ArrayList<SNP> getTagSNPs(){
    	return tagger.getTagSNPs();
    }
}
