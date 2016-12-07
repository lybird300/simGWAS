package haploview.tagger;

import java.util.*;
//import java.util.concurrent.TimeUnit;

public class Tagger {
    public static final double DEFAULT_RSQ_CUTOFF = 0.8;
    public static final double DEFAULT_LOD_CUTOFF = 3.0;
    public static final short DEFAULT_MIN_DISTANCE = 0;
    //private static final int PAIRWISE_ONLY = 0;
    //private static final int AGGRESSIVE_DUPLE = 1;
    //private static final int AGGRESSIVE_TRIPLE = 2;
    public static final int DEFAULT_MAXDIST = 500000;
    //public static final int DEFAULT_MAXDIST = 100000;
    private static final int DEFAULT_MAXNUMTAGS = 6000;
    public static final double DEFAULT_MIN_DESIGNSCORE = 0;

    //vector of SNP objects, which contains every SNP (tags and non-tags)
    private List<SNP> snps;

    //vector just those SNPs which should be captured.
    private List<SNP> capture;

    //vector of SNPs which must be included in the set of Tags
    //no object may be present in forceInclude and forceExclude concurrently.
    private List<SNP> forceInclude = new ArrayList<SNP>();


    //private Hashtable designScores = null;

    private HaploviewAlleleCorrelator alleleCorrelator;
    private double meanRSq;
    private int percentCapped;
    private double minRSquared;
    private double minDPrime = 0.5;
    //private int aggression;
    private int maxNumTags;
    private long maxComparisonDistance;
    //ArrayList of Tag objects determined by the most recent call to findTags()
    private ArrayList<TagSequence> tags;
    //vector of sites which are not tagged by anything in the tags vector
    private ArrayList<SNP> untagged;

    public int taggedSoFar;

    public Tagger(List<SNP> allSNPs, List<SNP> capture, HaploviewAlleleCorrelator ac){
        this(allSNPs, capture, ac,DEFAULT_RSQ_CUTOFF/*,PAIRWISE_ONLY*/,DEFAULT_MAXDIST,DEFAULT_MAXNUMTAGS,false);
    }
    
    public Tagger(List<SNP> allSNPs, List<SNP> capture, HaploviewAlleleCorrelator ac, int maxDist, int maxNumTags){
        this(allSNPs, capture, ac,DEFAULT_RSQ_CUTOFF/*,PAIRWISE_ONLY*/,maxDist,maxNumTags,false);
    }

    private Tagger(List<SNP> allSNPs, List<SNP> capture, HaploviewAlleleCorrelator ac, double rsqCut,
                  /*int aggressionLevel,*/ long maxCompDist, int maxNumTags, boolean printAllTags){
        minRSquared = rsqCut;
        //aggression = aggressionLevel;
        this.maxNumTags = maxNumTags;
        if(maxCompDist < 0 ) {
            maxComparisonDistance = DEFAULT_MAXDIST;
        } else {
            maxComparisonDistance = maxCompDist;
        }

        if(allSNPs != null) {
            snps = allSNPs;
        } else {
            snps = new ArrayList<SNP>();
        }
        
        if (capture != null){
            this.capture = capture;
        }else{
            this.capture = new ArrayList<SNP>();
        }

        alleleCorrelator = ac;


        for(int i=0;i<snps.size();i++) {
            VariantSequence curVarSeq = (VariantSequence)snps.get(i);
            if(curVarSeq.getTagComparator() == null) {
                TagRSquaredComparator trsc = new TagRSquaredComparator(curVarSeq);
                curVarSeq.setTagComparator(trsc);
            }
        }
    }

    private double getPairwiseCompRsq(VariantSequence a, VariantSequence b){
        return getPairwiseComp(a,b).getRsq();
    }
    
    private double getPairwiseCompDprime(VariantSequence a, VariantSequence b){
    	return alleleCorrelator.getCorrelation(a,b).getDPrime();
    }

    private LocusCorrelation getPairwiseComp(VariantSequence a, VariantSequence b) {
        return alleleCorrelator.getCorrelation(a,b);
    }

    /**
     * This method finds Tags for the SNPs in the snps vector, and returns a vector of Tag objects
     */
    public ArrayList<TagSequence> findTags() {
        tags = new ArrayList<TagSequence>();
        untagged = new ArrayList<SNP>();
        taggedSoFar = 0;

        //potentialTagsHash stores the PotentialTag objects keyed on the corresponding sequences
        Hashtable<VariantSequence, PotentialTag> potentialTagByVarSeq = new Hashtable<VariantSequence, PotentialTag>();
        PotentialTagComparator ptcomp = new PotentialTagComparator();
        SNP currentVarSeq;

        //create SequenceComparison objects for each potential Tag, and
        //add any comparisons which have an r-squared greater than the minimum.
        int ctPassRsqThresh = 0;
        int ctPassDistThresh = 0;
        for(int i=0;i<snps.size();i++) {
            currentVarSeq = snps.get(i);
            PotentialTag tempPT = new PotentialTag(currentVarSeq);
            for(int j=0;j<capture.size();j++) {
            	if(Math.abs(currentVarSeq.getPosition() - capture.get(j).getPosition()) <= maxComparisonDistance)  {
            		ctPassDistThresh++;
            		double r2 = getPairwiseCompRsq(currentVarSeq, capture.get(j));
            		//double dprime = getPairwiseCompDprime(currentVarSeq, capture.get(j));
            		if(/*r2 > 1 && dprime < 0)||*/r2 >= minRSquared /*&& r2 <= 1*/){
                    	ctPassRsqThresh++;
                        tempPT.addTagged(capture.get(j));
                    }
                }
            	else
            		break;
            }
            potentialTagByVarSeq.put(currentVarSeq,tempPT);
        }
        System.out.println(ctPassDistThresh + " SNPs passed the distance threshold.\n" + ctPassRsqThresh + " SNPs passed the rsq threshold.\n");

        List<SNP> sitesToCapture = new ArrayList<SNP>();
        sitesToCapture.addAll(capture);
        capture = null;

        //HaploText.logger.debug("snps to tag: " + sitesToCapture.size());

        int countTagged = 0;

            //loop until all snps are tagged
        while(sitesToCapture.size() > 0) {
                ArrayList<PotentialTag> potentialTags = new ArrayList<PotentialTag>(potentialTagByVarSeq.values());
                if(potentialTags.size() == 0) {
                    //we still have sites left to capture, but we have no more available tags.
                    //this should only happen if the sites remaining in sitesToCapture were specifically
                    //excluded from being tags. Since we can't add any more tags, break out of the loop.
                    break;
                }

                //sorts the array of potential tags according to the number of untagged sites they can tag.
                //the last element is the one which tags the most untagged sites, so we choose that as our next tag.
                Collections.sort(potentialTags, ptcomp);
                PotentialTag currentBestTag = potentialTags.get(potentialTags.size()-1);

                if (currentBestTag.taggedCount() > 0){
                    //there are no tags which can be chosen -- this should only happen if you've force
                    //excluded something which you want to capture, and none of the available tags can
                    //capture any of the requested alleles
                    HashSet<VariantSequence> newlyTagged = addTag(currentBestTag,potentialTagByVarSeq);
                    countTagged += newlyTagged.size();
                    sitesToCapture.removeAll(newlyTagged);
                    sitesToCapture.remove(currentBestTag.sequence);
                }else{
                    potentialTagByVarSeq.remove(currentBestTag.sequence);
                }
        }
        taggedSoFar = countTagged;

        if(sitesToCapture.size() > 0) {
            //any sites left in sitesToCapture could not be tagged, so we add them all to the untagged ArrayList
            untagged.addAll(sitesToCapture);
        }

        System.out.println("tagged " + countTagged + " SNPS using " + tags.size() +" tags" );
        System.out.println(untagged.size() + " SNPs could not be tagged.\n");

        /*if (aggression != PAIRWISE_ONLY){
            //peelback starting with the worst tag (i.e. the one that tags the fewest other snps.
            ArrayList<TagSequence> tags2BPeeled = new ArrayList<TagSequence>();
            tags2BPeeled.addAll(tags);
            Collections.reverse(tags2BPeeled);
            long startTime = System.currentTimeMillis();
            peelBack(tags2BPeeled);
            long stopTime = System.currentTimeMillis();
			System.out.println("Peeling back takes approximately " + TimeUnit.MILLISECONDS.toMinutes(stopTime - startTime) + " minutes.\n");
        }*/

        //we've done the best we can. now we check to see if there's a limit to the
        //num of tags we're allowed to choose.
        if (maxNumTags > 0){
            //if so we need to chuck out the extras. figure out the utility of each tagSNP
            //i.e. how many SNPs for which they and their combos are the only tags

            while (getTagSNPs().size() > maxNumTags){
                ArrayList<SNP> tagSNPs = getTagSNPs();
                potentialTagByVarSeq = new Hashtable<VariantSequence, PotentialTag>();
                Hashtable<PotentialTag, TagSequence> tagSeqByPotentialTag = new Hashtable<PotentialTag, TagSequence>();
                //account for stuff tagged by snps themselves
                for (int i = 0; i < tagSNPs.size(); i++){
                    TagSequence ts = new TagSequence(tagSNPs.get(i));
                    PotentialTag pt = new PotentialTag(ts.getSequence());
                    pt.addTagged(ts.getTagged());
                    potentialTagByVarSeq.put(ts.getSequence(),pt);
                    tagSeqByPotentialTag.put(pt,ts);
                }
                tagSNPs = null;
                //go through all pt's and add their utilities as members of combos
                List<TagSequence> tagHaps = getTagHaplotypes();
                for (int i = 0; i < tagHaps.size(); i++){
                    TagSequence ts = (TagSequence) tagHaps.get(i);
                    Block b = (Block) ts.getSequence();
                    for (int j = 0; j < b.getSnps().size(); j++){
                        (potentialTagByVarSeq.get(b.getSNP(j))).addTagged(ts.getTagged());
                    }
                }

                //now perform the steps of sorting and peeling              
                ArrayList<PotentialTag> potTagVec = new ArrayList<PotentialTag>();
                potTagVec.addAll(potentialTagByVarSeq.values());
                Collections.sort(potTagVec,ptcomp);
                //Now all tags (either a SNP or a block) are ordered (ascending) by the number of others they can tag

                int count = 0;
                PotentialTag dumpedPT = potTagVec.get(count);//dumptedPT is a potential tag to dump
                /*while (forceInclude.contains(dumpedPT.sequence)){
                    count++;
                    dumpedPT = (PotentialTag)potTagVec.get(count);
                }*/
                TagSequence dumpedTS = tagSeqByPotentialTag.get(dumpedPT);//get the tag of dumptedPT
                ArrayList<VariantSequence> taggedByCurTag = dumpedTS.getTagged();
                int j = 0;
                for (; j < taggedByCurTag.size(); j++){
                    //note for everything tagged by this guy that they're no longer tagged by him
                    VariantSequence vs =  (VariantSequence)taggedByCurTag.get(j);
                    vs.removeTag(dumpedTS);
                    if (vs.getTags().size() == 0){
                        taggedSoFar--;
                    }
                }
                if(j > 0) tagHaps = getTagHaplotypes();
                for (int i = 0; i < tagHaps.size(); i++){
                    TagSequence ts = (TagSequence) tagHaps.get(i);
                    Block b = (Block) ts.getSequence();
                    boolean isContain = false;
                    String snp1 = dumpedTS.getSequence().getName();
                    Iterator<SNP> itr = b.getSnps().iterator();
                    while(itr.hasNext()) {
                    	String snp2 = (itr.next()).getName();
                    	if(snp1.equalsIgnoreCase(snp2)){
                    		isContain = true;
                    		break;
                    	}
                    }
                    if (isContain){
                        //this hap tag is now defunct because it was comprised in part by dumpedTS
                        ArrayList<VariantSequence> taggedByHap = ts.getTagged();
                        for (int k = 0; k < taggedByHap.size(); k++){
                            VariantSequence vs =  (VariantSequence)taggedByHap.get(k);
                            vs.removeTag(ts);
                            if (vs.getTags().size() == 0){
                                taggedSoFar--;
                            }
                        }
                        tags.remove(ts);
                    }
                }
                if(!tags.remove(dumpedTS)){
                	System.out.println("Not removed.");
                }
            }
        }

        return new ArrayList<TagSequence>(tags);
    }

    /*private void peelBack(ArrayList<TagSequence> tagsToBePeeled){
        Hashtable<Allele, TagSequence> blockTagsByAllele = new Hashtable<Allele, TagSequence>();
        HashSet<VariantSequence> snpsInBlockTags = new HashSet<VariantSequence>();

        //HaploText.logger.debug("starting peelback. untagged.size() = " + untagged.size());

        ArrayList<VariantSequence> availTagSNPs = new ArrayList<VariantSequence>();
        for (int j = 0; j < tags.size(); j++){
            availTagSNPs.add(tags.get(j).getSequence());
        }

        ListIterator<SNP> uitr = untagged.listIterator();

        //try to tag things that weren't taggable in pairwise with haps
        while(uitr.hasNext()) {
            SNP curSnp = uitr.next();
            HashSet<SNP> comprehensiveBlock = new HashSet<SNP>();

            comprehensiveBlock.add(curSnp);
            HashSet<SNP> victor = curSnp.getLDList();
            victor.retainAll(availTagSNPs);
            comprehensiveBlock.addAll(victor);

            alleleCorrelator.phaseAndCache(comprehensiveBlock);

            LocusCorrelation bestPredictor = null;
            ArrayList<TagSequence> availableTags = new ArrayList<TagSequence>();
            availableTags.addAll(tags);
            ArrayList<Block> potentialTests = generateTests(curSnp, availableTags);
            for (int j = 0; j < potentialTests.size(); j++){
                LocusCorrelation lc = getPairwiseComp(potentialTests.get(j), curSnp);
                if (lc.getRsq() >= minRSquared){
                    if (bestPredictor != null){
                        if (lc.getRsq() > bestPredictor.getRsq()){
                            bestPredictor = lc;
                        }
                    }else{
                        bestPredictor= lc;
                    }
                }
            }

            if(bestPredictor != null) {
                Allele bpAllele = bestPredictor.getAllele();
                snpsInBlockTags.addAll(((Block)bpAllele.getLocus()).getSnps());
                if (blockTagsByAllele.containsKey(bpAllele)){
                    TagSequence ts = blockTagsByAllele.get(bpAllele);
                    ts.addTagged(curSnp);
                }else{
                    TagSequence ts = new TagSequence(bpAllele);
                    ts.addTagged(curSnp);
                    tags.add(ts);
                    blockTagsByAllele.put(bpAllele,ts);
                }
                uitr.remove();
                //note that we've caught another SNP
                taggedSoFar++;
            }
        }

        //HaploText.logger.debug("finished attempt at pairwise untaggables. untagged.size() = " + untagged.size());
        int peelct = 0;
        for (int i = 0; i < tagsToBePeeled.size(); i++){     	
            TagSequence curTag = tagsToBePeeled.get(i);
            if (forceInclude.contains(curTag.getSequence()) ||
                    snpsInBlockTags.contains(curTag.getSequence())){
                continue;
            }
            ArrayList<VariantSequence> taggedByCurTag = curTag.getTagged();

            //a hashset that contains all snps tagged by curtag
            //and all tag snps in LD with any of them
            HashSet<SNP> comprehensiveBlock = new HashSet<SNP>();
            availTagSNPs = new ArrayList<VariantSequence>();
            for (int j = 0; j < tags.size(); j++){
                availTagSNPs.add(tags.get(j).getSequence());
            }
            availTagSNPs.remove(curTag.getSequence());
            for (int j = 0; j < taggedByCurTag.size(); j++) {
                SNP snp = (SNP) taggedByCurTag.get(j);
                comprehensiveBlock.add(snp);
                HashSet<SNP> victor = snp.getLDList();
                victor.retainAll(availTagSNPs);
                comprehensiveBlock.addAll(victor);
            }
            alleleCorrelator.phaseAndCache(comprehensiveBlock);

            Hashtable<SNP, LocusCorrelation> bestPredictor = new Hashtable<SNP, LocusCorrelation>();
            boolean peelSuccessful = true;
            for (int k = 0; k < taggedByCurTag.size(); k++){
                //look to see if we can find a predictor for each thing curTag tags
            	SNP thisTaggable = (SNP) taggedByCurTag.get(k);
                ArrayList<TagSequence> victor = new ArrayList<TagSequence>();
                victor.addAll(tags);
                victor.remove(curTag);
                ArrayList<Block> potentialTests = generateTests(thisTaggable, victor);
                for (int j = 0; j < potentialTests.size(); j++){
                    LocusCorrelation lc = getPairwiseComp(potentialTests.get(j), thisTaggable);
                    if (lc.getRsq() >= minRSquared){
                        if (bestPredictor.containsKey(thisTaggable)){
                            if (lc.getRsq() >
                                    (bestPredictor.get(thisTaggable)).getRsq()){
                                bestPredictor.put(thisTaggable,lc);
                            }
                        }else{
                            bestPredictor.put(thisTaggable,lc);
                        }
                    }
                }
                if (thisTaggable.getTags().size() == 1 && !bestPredictor.containsKey(thisTaggable)){
                    peelSuccessful = false;
                    break;
                }
            }
            if (peelSuccessful){
                for (int k = 0; k < taggedByCurTag.size(); k++){
                    SNP thisTaggable = (SNP) taggedByCurTag.get(k);
                    //if more than one snp is tagged by the same
                    if (bestPredictor.containsKey(thisTaggable)){
                        Allele bpAllele = (bestPredictor.get(thisTaggable)).getAllele();
                        snpsInBlockTags.addAll(((Block)bpAllele.getLocus()).getSnps());
                        if (blockTagsByAllele.containsKey(bpAllele)){
                            TagSequence ts = blockTagsByAllele.get(bpAllele);
                            ts.addTagged(thisTaggable);
                        }else{
                            TagSequence ts = new TagSequence(bpAllele);
                            ts.addTagged(thisTaggable);
                            tags.add(ts);
                            blockTagsByAllele.put(bpAllele,ts);
                        }
                    }
                    thisTaggable.removeTag(curTag);
                }
                boolean success = tags.remove(curTag);
                if(success) peelct++;
            }
        }

        //this removes multimarker tags that are not the best tag for any alleles
         for (int i = 0; i < tags.size(); i++){
            if (tags.get(i).getBestTagged().size() == 0){
                TagSequence ts = tags.get(i);
                ArrayList<VariantSequence> taggedByHap = ts.getTagged();
                for (int j = 0; j < taggedByHap.size(); j++){
                    VariantSequence vs = taggedByHap.get(j);
                    vs.removeTag(ts);
                    if (vs.getTags().size() == 0){
                       //this should never happen!
                    }
                }
                //HaploText.logger.debug(((TagSequence)tags.get(i)).getSequence().getName() + ": " + ((TagSequence)tags.get(i)).getTagged().size() + " " + ((TagSequence)tags.get(i)).getBestTagged().size());
                //HaploText.logger.debug("Removing the above tag since it isn't the best tag for anything.");
                if(tags.remove(i)!=null){
                	peelct++;
                	i--;
                }
            }
         }
         System.out.println(peelct + " tags were peeled off.\n");
    }
    
    //returns all duples and triples (2 or 3-marker blocks) from availTags which are in LD with SNP s, and with each other
    private ArrayList<Block> generateTests(SNP s, ArrayList<TagSequence> availTags){      
        HashSet<SNP> tagsInLD = new HashSet<SNP>();
        tagsInLD.addAll(s.getLDList());
        
        ArrayList<VariantSequence> availTagSNPs = new ArrayList<VariantSequence>();
        for (int i = 0; i < availTags.size(); i++){
            availTagSNPs.add(availTags.get(i).getSequence());
        }
        tagsInLD.retainAll(availTagSNPs);
        HashSet<Block> tests = new HashSet<Block>();
        Iterator<SNP> lditr = tagsInLD.iterator();
        while (lditr.hasNext()){
            SNP curTag = lditr.next(); 
            HashSet<SNP> hs = new HashSet<SNP>();
            hs.addAll(curTag.getLDList());
            hs.retainAll(tagsInLD);
            ArrayList<SNP> victor = new ArrayList<SNP>();
            victor.addAll(hs);
            if (aggression == AGGRESSIVE_DUPLE || aggression == AGGRESSIVE_TRIPLE){
                //2 marker blocks
                for (int i = 0; i < victor.size(); i++) {
                    ArrayList<SNP> block = new ArrayList<SNP>();
                    block.add(curTag);
                    block.add(victor.get(i));
                    tests.add(new Block(block));
                    if (aggression == AGGRESSIVE_TRIPLE){
                        //3 marker blocks
                        for(int j=i+1;j<victor.size();j++) {
                            //make sure these two snps are in LD with each other
                            if (victor.get(i).getLDList().contains(victor.get(j))){
                                ArrayList<SNP> block2 = new ArrayList<SNP>();
                                block2.addAll(block);
                                block2.add(victor.get(j));
                                tests.add(new Block(block2));
                            }
                        }
                    }
                }
            }
        }

        return new ArrayList<Block>(tests);
    }*/

	private HashSet<VariantSequence> addTag(PotentialTag theTag,Hashtable<VariantSequence, PotentialTag> potentialTagHash) {
        potentialTagHash.remove(theTag.sequence);
        //newlyTagged contains alleles which were not tagged by anything in the set of tags before,
        //and are now tagged by theTag.
        HashSet<VariantSequence> newlyTagged = theTag.tagged;

        TagSequence tagSeq = new TagSequence(theTag.sequence);
        tags.add(tagSeq);

        Iterator<VariantSequence> itr = potentialTagHash.keySet().iterator();
        ArrayList<VariantSequence> toRemove = new ArrayList<VariantSequence>();
        //iterate through the list of available tags, and remove the newly tagged alleles from
        //the list of alleles that each PotentialTag can tag. (since we want to choose our next tag
        // according to which will tag the most untagged alleles )
        while(itr.hasNext()) {
            PotentialTag pt = (PotentialTag) potentialTagHash.get(itr.next());
            pt.removeTagged(newlyTagged);
            //if a PotentialTag cannot tag any other uncaptured sites, then we want to remove it from contention,
            //unless its sequence still needs to be captured.
            if(pt.taggedCount() == 0){
                toRemove.add(pt.sequence);
            }
        }

        for(int i=0;i<toRemove.size();i++) {
            potentialTagHash.remove(toRemove.get(i));
        }

        //loop through the list of alleles the newly added tag can capture, and
        //add them to the TagSequence object.
        //we add all the alleles the tag can capture, _not_ just the newly tagged alleles.
        Iterator<VariantSequence> ptitr = theTag.allTagged.iterator();
        while(ptitr.hasNext()) {
            tagSeq.addTagged((VariantSequence)ptitr.next());
        }

        return newlyTagged;
    }

    class PotentialTag {
        VariantSequence sequence;
        // tagged contains the sequences which this sequence can tag, which are not yet tagged
        //(this is used in the while loop in findTags() )
        HashSet<VariantSequence> tagged;
        //allTagged contains all sequences that this sequence can tag, regardless of what tags have already been chosen
        HashSet<VariantSequence> allTagged;

        public PotentialTag(VariantSequence s) {
            sequence = s;
            tagged = new HashSet<VariantSequence>();
            allTagged = new HashSet<VariantSequence>();
        }

        public void addTagged(VariantSequence vs) {
            tagged.add(vs);
            allTagged.add(vs);
        }

        public void addTagged(Collection<VariantSequence> c){
            tagged.addAll(c);
            allTagged.addAll(c);
        }

        public void removeTagged(VariantSequence vs) {
            tagged.remove(vs);
        }

        public void removeTagged(Collection<VariantSequence> c) {
            tagged.removeAll(c);
        }

        //this returns the number of sequences that havent been tagged this sequence can tag
        public int taggedCount() {
            return tagged.size();
        }

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime * result
					+ ((sequence == null) ? 0 : sequence.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (obj != null && obj instanceof PotentialTag) {
				PotentialTag other = (PotentialTag) obj;
				if (sequence != null){
					if(other.sequence == null)
						return false;
					else if(sequence.equals(other.sequence)){
						return true;
					}
				}
			}
			return false;
		}

		private Tagger getOuterType() {
			return Tagger.this;
		}
    }

    

    public ArrayList<TagSequence> getTags() {
        return tags;
    }

    public ArrayList<SNP> getTagSNPs(){
        ArrayList<SNP> res = new ArrayList<SNP>();
        Iterator<TagSequence> itr = tags.iterator();
        while (itr.hasNext()){
            TagSequence t = (TagSequence) itr.next();
            if (t.getSequence() instanceof SNP){
                res.add((SNP)t.getSequence());
            }
        }
        return res;
    }

    public List<TagSequence> getTagHaplotypes(){
        List<TagSequence> res = new ArrayList<TagSequence>();
        Iterator<TagSequence> itr = tags.iterator();
        while (itr.hasNext()){
            TagSequence t = (TagSequence) itr.next();
            if (t.getSequence() instanceof Block){
                res.add(t);
            }
        }
        return res;
    }

    public List<SNP> getForceInclude() {
        return forceInclude;
    }

    public int getUntaggableCount() {
        return untagged.size();
    }

    /*public ArrayList saveResults(File outFile) throws IOException {
    	ArrayList res = getTagSNPs();
        BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));

        bw.write("#captured " + taggedSoFar + " of " + capture.size() +" alleles at r^2 >= " + minRSquared);
        bw.newLine();
        bw.write("#captured " + percentCapped + " percent of alleles with mean r^2 of " + Util.roundDouble(meanRSq, 3));
        bw.newLine();
        bw.write("#using " + res.size() + " Tag SNPs in " + tags.size() + " tests.");
        bw.newLine();

        bw.write("Allele\tBest Test\tr^2 w/test");
        bw.newLine();
        for (int i = 0; i < capture.size(); i++) {
            StringBuffer line = new StringBuffer();
            SNP snp = (SNP) capture.get(i);
            line.append(snp.getName()).append("\t");
            TagSequence theTag = snp.getBestTag();
            if(theTag != null) {
                line.append(theTag.getName()).append("\t");
                line.append(getPairwiseCompRsq(snp,theTag.getSequence())).append("\t");
            }
            bw.write(line.toString());
            bw.newLine();
        }

        bw.newLine();

        bw.write("Test\tAlleles Captured");
        bw.newLine();
        for(int i=0;i<tags.size();i++) {
            StringBuffer line = new StringBuffer();
            TagSequence theTag = (TagSequence) tags.get(i);
            line.append(theTag.getName()).append("\t");
            List<VariantSequence> tagged = null;
            if (printAllTags){
                tagged = theTag.getTagged();
            }else{
                tagged = theTag.getBestTagged();
            }
            for (int j = 0; j < tagged.size(); j++) {
                VariantSequence varSeq = (VariantSequence) tagged.get(j);
                if(j !=0){
                    line.append(",");
                }
                line.append(varSeq.getName());
            }
            bw.write(line.toString());
            bw.newLine();
        }

        bw.close();
        return res;
    }*/

    class PotentialTagComparator implements Comparator<PotentialTag> {
        public int compare(PotentialTag o1, PotentialTag o2) {
            return o1.taggedCount() -  o2.taggedCount();
        }
    }

    class TagRSquaredComparator implements Comparator<TagSequence> {
        VariantSequence seq;

        public TagRSquaredComparator(VariantSequence s) {
            seq = s;
        }

        public int compare(TagSequence o1, TagSequence o2) {
            //if one of the compared tags actually is this sequence, always promote
            //it to the front (i.e. a SNP should always pick itself as its own best tag
            //if possible).
            if (seq.equals(o1.getSequence())){
                return 1;
            }else if (seq.equals(o2.getSequence())){
                return -1;
            }

            if(getPairwiseCompRsq(seq,o1.getSequence()) ==
                    getPairwiseCompRsq(seq, o2.getSequence())) {
                return 0;
            } else if (getPairwiseCompRsq(seq, o1.getSequence()) >
                    getPairwiseCompRsq(seq, o2.getSequence())) {
                return 1;
            } else {
                return -1;
            }
        }
    }

    public double getMeanRSq() {
        return meanRSq;
    }

    public int getPercentCaptured() {
        return percentCapped;
    }
}
