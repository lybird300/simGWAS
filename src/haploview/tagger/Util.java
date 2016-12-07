package haploview.tagger;

import java.io.File;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.FieldPosition;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Map.Entry;

import simGWAS.FileHandler;
import simGWAS.ParamParser;

public abstract class Util {

    public static String formatPValue(double pval){
        DecimalFormat df;
        //java truly sucks for simply restricting the number of sigfigs but still
        //using scientific notation when appropriate
        if (pval < 0.0001){
            df = new DecimalFormat("0.0000E0", new DecimalFormatSymbols(Locale.US));
        }else{
            df = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.US));
        }
        String formattedNumber =  df.format(pval, new StringBuffer(), new FieldPosition(NumberFormat.INTEGER_FIELD)).toString();
        return formattedNumber;
    }

    public static double roundDouble (double d, int places){
        double factor = Math.pow(10,places);
        return Math.rint(d*factor)/factor;
    }
    
	//Get indices of specific marker SNPs on selected chromsomes.
	public static Map<Integer, int[]> gatherMutMarkersfromDistFiles(Integer[] chroms, Integer[] markerIndex) {		
		Arrays.sort(chroms);
		Arrays.sort(markerIndex);
		Map<Integer, int[]> chromsToMarkerSNPs = new HashMap<Integer, int[]>();
		Map<String, ArrayList<Integer>> chromRecord = getOriginDataLocation(ParamParser.mutInput, chroms);	
		for(Entry<String, ArrayList<Integer>> entry: chromRecord.entrySet()){
			FileHandler readMarkerData = new FileHandler(ParamParser.dataOutputDir, entry.getKey());
			readMarkerData.openForRead();
			readMarkerData.nextLine();//first line is the header
			ArrayList<Integer> records = entry.getValue();
			int maxRecIndex = Collections.max(records);
			Integer[] recordsInThisFile = records.toArray(new Integer[records.size()]);
			Arrays.sort(recordsInThisFile);		
			String nextLine = null;
			for(int recordCt = 0; (nextLine = readMarkerData.nextLine()) != null && recordCt <= maxRecIndex; recordCt++){
				if(Arrays.binarySearch(recordsInThisFile, recordCt) < 0)
					continue;//this record is not what we want
				String[] s = nextLine.split("\\s+");					
				int chromID = Integer.parseInt(s[0]);//In out.pos-1_# files, s[0] is chromID
				if(Arrays.binarySearch(chroms, chromID) < 0)
					System.out.println("There is something wrong with you strategy in reading distributed files.\n");
				int[] mutations = new int[s.length-1];
				for(int j = 0; j < s.length-1; j++)
					mutations[j] = Integer.valueOf(s[j+1]);
				Arrays.sort(mutations);
				int[] markers = new int[markerIndex.length];
				for(int k = 0; k < markerIndex.length; k++){
					if(Arrays.binarySearch(mutations, markerIndex[k]) >= 0)
						markers[k] = 1;
					else
						markers[k] = 2;
				}
				chromsToMarkerSNPs.put(chromID, markers);				
			}				
			readMarkerData.closeReader();
		}
		return chromsToMarkerSNPs;
	}
	
	public static int[] gatherMutMarkersfromDistFiles(int chromID, Integer[] markerIndex) {		
		Arrays.sort(markerIndex);
		int[] markers = new int[markerIndex.length];
		int fileSize = ParamParser.sampleSize/ParamParser.inputFileNum;
		String fileName = ParamParser.mutInput + "_" + (chromID/fileSize) + ".gz";
		FileHandler readMarkerData = new FileHandler(ParamParser.dataOutputDir, fileName);
		readMarkerData.openForRead();
		readMarkerData.nextLine();//first line is the header		
		String nextLine = null;
		while((nextLine = readMarkerData.nextLine()) != null){
			String[] s = nextLine.split("\\s+");					
			if(Integer.parseInt(s[0]) != chromID)//In out.pos-1_# files, s[0] is chromID
				continue;
			int[] mutations = new int[s.length-1];
			for(int j = 0; j < s.length-1; j++)
				mutations[j] = Integer.valueOf(s[j+1]);
			Arrays.sort(mutations);
			
			for(int k = 0; k < markerIndex.length; k++){
				if(Arrays.binarySearch(mutations, markerIndex[k]) >= 0)
					markers[k] = 1;
				else
					markers[k] = 2;
			}
			break;
		}
		readMarkerData.closeReader();
		return markers;
	}
	
	public static char[] getFullHap(int chromID) {		
		char[] markers = new char[ParamParser.num_mut];
		int fileSize = ParamParser.sampleSize/ParamParser.inputFileNum;
		String fileName = ParamParser.mutInput + "_" + (chromID/fileSize) + ".gz";
		FileHandler readMarkerData = new FileHandler(ParamParser.dataOutputDir, fileName);
		readMarkerData.openForRead();
		readMarkerData.nextLine();//first line is the header		
		String nextLine = null;
		while((nextLine = readMarkerData.nextLine()) != null){
			String[] s = nextLine.split("\\s+");					
			if(Integer.parseInt(s[0]) != chromID)//In out.pos-1_# files, s[0] is chromID
				continue;
			int[] mutations = new int[s.length-1];
			for(int j = 0; j < s.length-1; j++)
				mutations[j] = Integer.valueOf(s[j+1]);
			Arrays.sort(mutations);
			
			for(int k = 0; k < markers.length; k++){
				if(Arrays.binarySearch(mutations, k) >= 0)
					markers[k] = '1';
				else
					markers[k] = '2';
			}
			break;
		}
		readMarkerData.closeReader();
		return markers;
	}
	
	public static Map<Integer, String> getHapsfromDistFiles(Integer[] chroms, Integer[] markerIndex) {		
		Arrays.sort(chroms);
		if(markerIndex!=null) Arrays.sort(markerIndex);
		Map<Integer, String> chromsToMarkerSNPs = new HashMap<Integer, String>();
		Map<String, ArrayList<Integer>> chromRecord = getOriginDataLocation(ParamParser.mutInput, chroms);	
		for(Entry<String, ArrayList<Integer>> entry: chromRecord.entrySet()){
			FileHandler readMarkerData = new FileHandler(ParamParser.dataOutputDir, entry.getKey());
			readMarkerData.openForRead();
			readMarkerData.nextLine();//first line is the header
			ArrayList<Integer> records = entry.getValue();
			int maxRecIndex = Collections.max(records);
			Integer[] recordsInThisFile = records.toArray(new Integer[records.size()]);
			Arrays.sort(recordsInThisFile);		
			String nextLine = null;
			for(int recordCt = 0; (nextLine = readMarkerData.nextLine()) != null && recordCt <= maxRecIndex; recordCt++){
				if(Arrays.binarySearch(recordsInThisFile, recordCt) < 0)
					continue;//this record is not what we want
				String[] s = nextLine.split("\\s+");					
				int chromID = Integer.parseInt(s[0]);//In out.pos-1_# files, s[0] is chromID
				if(Arrays.binarySearch(chroms, chromID) < 0)
					System.out.println("There is something wrong with you strategy in reading distributed files.\n");
				int[] mutations = new int[s.length-1];
				for(int j = 0; j < s.length-1; j++)
					mutations[j] = Integer.valueOf(s[j+1]);
				Arrays.sort(mutations);
				int length, compareValue;
				if(markerIndex==null) length = ParamParser.num_mut;
				else length = markerIndex.length; 
				char[] markers = new char[length];
				for(int k = 0; k < length; k++){
					if(markerIndex==null) compareValue = k;
					else compareValue = markerIndex[k];
					if(Arrays.binarySearch(mutations, compareValue) >= 0)
						markers[k] = '1';
					else
						markers[k] = '2';
				}
				chromsToMarkerSNPs.put(chromID, new String(markers));				
			}				
			readMarkerData.closeReader();
		}
		return chromsToMarkerSNPs;
	}
    
	/**
	 * Returns a pseudo-random number between min and max, inclusive.
	 * The difference between min and max can be at most
	 * <code>Integer.MAX_VALUE - 1</code>.
	 *
	 * @param min Minimum value
	 * @param max Maximum value.  Must be greater than min.
	 * @return Integer between min and max, inclusive.
	 * @see java.util.Random#nextInt(int)
	 */
	public static int randInt(int min, int max) {
		Random randomGenerator = new Random();
	    // nextInt is normally exclusive of the top value, so add 1 to make it inclusive
	    return randomGenerator.nextInt((max - min) + 1) + min;
	}
	
	/**
	 * Find out where the information of each chromosome exists and store the results in a map.
	 * Each entry of this map is <name of the file, line in the file>
	 * @param filePrefix
	 * @param usedChroms
	 * @return
	 */
	public static Map<String, ArrayList<Integer>> getOriginDataLocation(String filePrefix, Integer[] usedChroms){
		Map<String, ArrayList<Integer>> fileNameToChromIDs = new HashMap<String, ArrayList<Integer>>();
		if(usedChroms==null){
			for(int j = 0; j< ParamParser.inputFileNum; j++)
				fileNameToChromIDs.put(filePrefix + "_" + j + ".gz", new ArrayList<Integer>());
		}
		else{
			for(int i = 0; i<usedChroms.length; i++){
				int fileSize = ParamParser.sampleSize/ParamParser.inputFileNum;
				String fileName = filePrefix + "_" + (usedChroms[i]/fileSize) + ".gz";
				if(fileNameToChromIDs.get(fileName)==null){
					fileNameToChromIDs.put(fileName, new ArrayList<Integer>());	
				}
				fileNameToChromIDs.get(fileName).add(usedChroms[i]%fileSize);	
			}
		}
		return fileNameToChromIDs;	
	}
	
	//get snp folder
	public static File returnSNPDir(int snp, double maf){
     	String dir_2nd_name = "SNP_" + snp + "_";
      	if(maf >= ParamParser.lo_maf[1])
      		dir_2nd_name += ParamParser.lo_maf[1]+ "_" + ParamParser.hi_maf[1];
      	else if(maf < ParamParser.lo_maf[1])
      		dir_2nd_name += ParamParser.lo_maf[0]+ "_" + ParamParser.hi_maf[0];
      	if(dir_2nd_name.charAt(dir_2nd_name.length()-1)=='0')//if the directory name ends with "0.0010", change it to "0.001"
      		dir_2nd_name = dir_2nd_name.substring(0, dir_2nd_name.length()-1);
		File dir_2nd = new File(ParamParser.dataOutputDir, dir_2nd_name);
		return dir_2nd;
	}
	
	public static Integer[] getSNPRealPosition(Integer[] SNPIndex) {
		Arrays.sort(SNPIndex);
		Integer[] snpPos = new Integer[SNPIndex.length];
		//Read the absolute positions of all SNPs from the "out.snp-1.gz" file
		FileHandler readSNPPos = new FileHandler(ParamParser.dataOutputDir, ParamParser.posInput);
		readSNPPos.openForRead();
		String nextLine = null;
		int i = 0;
		for(int recordCt = 0;(nextLine = readSNPPos.nextLine()) != null; recordCt++)
			if(Arrays.binarySearch(SNPIndex, recordCt) >= 0)
				snpPos[i++] = Integer.valueOf(nextLine);
		readSNPPos.closeReader();
		return snpPos;
	}
}
