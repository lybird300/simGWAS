package haploview.tagger;

import java.util.ArrayList;
import java.util.List;


public class TagSequence {
    private Allele allele;
    private ArrayList<VariantSequence> tagged;//A list of sites tagged by this object
    private VariantSequence sequence;

    public TagSequence(Allele a) {
        allele = a;
        sequence = a.getLocus();
        tagged = new ArrayList<VariantSequence>();
    }

    public TagSequence(VariantSequence snp){
        tagged = new ArrayList<VariantSequence>();
        sequence = snp;
        allele = null;
    }
    
    /**
     * @param t The site to add to the list of sites tagged by this object
     */
    public void addTagged(VariantSequence t) {
        tagged.add(t);
        t.addTag(this);
    }


    public ArrayList<VariantSequence> getTagged() {
        return tagged;
    }

    public Allele getAllele() {
        return allele;
    }

    public VariantSequence getSequence(){
        return sequence;
    }

    public String getName(){
        if (allele == null){
            return sequence.getName();
        }else{
            return sequence.getName() + " : " + allele.getGenotypeString();
        }
    }

    public String getTestName(){
        if (allele == null){
            return sequence.getName();
        }else{
            return sequence.getName() + "\t" + allele.getTestFileFormat();
        }
    }

    public List<VariantSequence> getBestTagged() {
        List<VariantSequence> result = new ArrayList<VariantSequence>();

        for (int i = 0; i < tagged.size(); i++) {
        	VariantSequence taggable = tagged.get(i);
            if(taggable.getBestTag() == this) {
                result.add(taggable);
            }
        }
        return result;
    }

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((allele == null) ? 0 : allele.hashCode());
		result = prime * result
				+ ((sequence == null) ? 0 : sequence.hashCode());
		result = prime * result + ((tagged == null) ? 0 : tagged.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof TagSequence){
			TagSequence other = (TagSequence) o;
	    	if (sequence != null) {
	    		if (other.sequence == null) {
	    				return false;
	    		} else if (sequence.equals(other.sequence)) {
	    			return true;
	    		}
	        }
		}
		return false;
	}
}
