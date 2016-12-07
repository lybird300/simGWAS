package haploview.tagger;

import java.util.List;

public class DPrimeTable {
    private PairwiseLinkage[][] theTable;

    public DPrimeTable(int numMarkers){
        theTable = new PairwiseLinkage[numMarkers][];
    }

    /*public void addMarker(Vector marker, int pos){
        theTable[pos] = (PairwiseLinkage[]) marker.toArray(new PairwiseLinkage[0]);
    }*/
    
    public void addMarker(List<PairwiseLinkage> dpTemp, int pos){
    	if(dpTemp == null)
    		theTable[pos]= null;
    	else if(dpTemp.size() == 0)
    		theTable[pos]= null;
    	else
    		theTable[pos] = (PairwiseLinkage[]) dpTemp.toArray(new PairwiseLinkage[dpTemp.size()]);
    }

    public PairwiseLinkage getLDStats(int pos1, int pos2){
        //we need to convert the input of an absolute position into the relative position
        //to index into the DP array. here we jump through the additional hoop of un-filtering the input
        //numbers
        int x = pos1;
        int y = pos2 - x - 1;
        if (x < theTable.length-1 && theTable[x] != null){
            if (y < theTable[x].length){
                return theTable[x][y];
            }else{
                return null;
            }
        }else{
            return null;
        }
    }

    public int getLength(int x){
        if (x >= theTable.length-1 || theTable[x]== null) return 0;
        else return theTable[x].length;
    }
    
    public void setValue(int i, int j, PairwiseLinkage pl){
    	theTable[i][j] = pl;
    }
}