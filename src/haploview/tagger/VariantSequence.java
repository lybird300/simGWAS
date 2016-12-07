package haploview.tagger;

import java.util.*;

public abstract class VariantSequence{

	private List<TagSequence> tags;
    private Comparator<TagSequence> tagComparator;

    public VariantSequence() {
        tags = new ArrayList<TagSequence>();
    }
    
    /*public VariantSequence(SNP s){
    	this();
    	variants.add(s);
    }*/

    /**
     * Adds t to the list of sites that tag this object.
     * @param t Tag object which tags this object 
     */
    public void addTag(TagSequence t) {
        tags.add(t);
    }

    public void removeTag(TagSequence t){
        tags.remove(t);
    }
    
    public abstract boolean equals(Object o);
    public abstract int hashCode();

    /**
     * returns a list of the Tags which capture this site
     * @return List list of Tag objects
     */
    public List<TagSequence> getTags(){
        return tags;
    }

    public void setTagComparator(Comparator<TagSequence> tagComparator) {
        this.tagComparator = tagComparator;
    }

    public Comparator<TagSequence> getTagComparator() {
        return tagComparator;
    }

    public abstract String getName();

    /**
     * Find and then return the "best" tag for this VariantSequence.
     * Uses the Comparator tagComparator to decide the "best" tag.
     * If no comparator has been set, then an arbitrary tag is returned from the Vector tags.
     * Best is defined by the particular implementation
     * @return Tag.  
     */
    public TagSequence getBestTag() {
        if(tags!= null && tags.size()>0) {
            //if a comparator has been set, then use it to sort the tags. the "best" tag will be the
            //last one in the sorted vector.
            //if there isn't a comparator, we have no basis for choosing which tag is "best",
            //so we just return whatever is currently in the last spot in the
            if(tagComparator != null) {
                Collections.sort(tags, tagComparator);
            }
            return tags.get(tags.size()-1);
        }
        return null;
    }
}
