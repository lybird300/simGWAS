package haploview.tagger;

import java.util.ArrayList;
import java.util.Iterator;

public class Block extends VariantSequence{
    //this will have tagger SNP objects
    private ArrayList<SNP> snps;

    public Block(ArrayList<SNP> block) {
        this.snps = block;
    }

    public boolean equals(Object o) {
        if (o instanceof Block){
            Block b = (Block)o;
            if(b.snps.size() != snps.size()) {
                return false;
            }
            for (int i = 0; i < snps.size(); i++){
                if (!b.snps.contains(snps.get(i))){
                    return false;
                }
            }
            return true;
        }
        return false;
    }

    public int hashCode() {
        int h = 0;
        Iterator<SNP> itr = snps.iterator();
        while(itr.hasNext()) {
            Object obj = itr.next();
            h += obj.hashCode();
        }
        return h;
    }

    public String getName() {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < snps.size()-1; i++){
            sb.append((snps.get(i)).getName()).append(",");
        }
        sb.append((snps.get(snps.size()-1)).getName());

        return sb.toString();
    }

    public int getMarkerCount() {
        return snps.size();
    }

    public SNP getSNP(int i) {
        return snps.get(i);
    }

    public String toString() {
        return getName();
    }

    public ArrayList<SNP> getSnps() {
        return snps;
    }
}
