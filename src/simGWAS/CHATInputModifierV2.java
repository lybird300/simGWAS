package simGWAS;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import haploview.tagger.Util;

public class CHATInputModifierV2 extends ChatInputGenerator {
	private int cv;
	private int repID;
	RecombWorker recomb = new RecombWorker();
	//private Set<Integer> genToGenerate = new HashSet<Integer>();
	private int gID;
	private Map<Integer, int[]> chromToMarkers = null;
	private List<Integer> chromsWithoutCV = null;
	private Map<Integer, Integer> allSNPPosToIDMap = null;

	public CHATInputModifierV2(int cv, Map<Integer, Integer> tagSNPs, File panelSpecificChroms, int caseIndex, int ctrlIndex, int replicate, 
			int simChromName, String recombFile, Map<Integer, Integer> allSNPPosToID) {
		super(tagSNPs, panelSpecificChroms, caseIndex, ctrlIndex, simChromName);
		repID = replicate;
		this.cv = cv;
		processRecombMap(recombFile);
		this.gID = simChromName - ParamParser.founderChrom;
		allSNPPosToIDMap = allSNPPosToID;
	}

	public void run(Collection<Integer> chromsWithSNP, int grrIndex, int bdrIndex, double snpPopMAF, String planGWAS, 
			Integer[] missingGenoPos, Integer[] snpToBeRandShuffle) {
		FileHandler writePEDFile = null;
		FileHandler recordRefCtrlCarriers = null;
		FileHandler writeControlRef = null;
		File ctrlcarrierRecord = null;
		Integer[] refCtrls = null;
		String chrLabel = ParamParser.padNumberAsString(simChromName,2);
		FileHandler readPanelSpecificChroms = new FileHandler(chromFile);
		readPanelSpecificChroms.openForRead();
		readPanelSpecificChroms.nextLine();//first line is the header	
		String nextLine = null;
		chromToMarkers = new HashMap<Integer, int[]>();
		while((nextLine = readPanelSpecificChroms.nextLine()) != null){
			String[] s = nextLine.split(" ");					
			if(Integer.valueOf(s[0]) != repID) continue;
			Integer[] panelSpecificChroms = new Integer[s.length-1];
			for(int j = 1; j < s.length; j++)
				panelSpecificChroms[j-1] = Integer.valueOf(s[j]);
			chromToMarkers = Util.gatherMutMarkersfromDistFiles(panelSpecificChroms, tag2MIndex);
			s = null;
			break;
		}
		readPanelSpecificChroms.closeReader();
		
		chromsWithoutCV = new ArrayList<Integer>();
	  	for(int k=0; k<ParamParser.sampleSize; k++){
	  		if(!chromsWithSNP.contains((Integer)k))
	  			chromsWithoutCV.add((Integer)k);
	  	}	
		
		int caseCt = ParamParser.num_ca[caseIndex];
		int ctrlCt = ParamParser.num_co[caseIndex];
		int numOfSubjects = caseCt + ctrlCt;
		int[][] genotypes = new int[tag20MIndex.length][numOfSubjects];
		
		File dir_4th = new File(panelFolderName + File.separator + "ChatInput_tagSNP" + File.separator + "Rep" + repID + File.separator + "ChromosomeSpecificData");
		if(!dir_4th.exists()) dir_4th.mkdirs();		
		
		if(planGWAS.equalsIgnoreCase("t")){
			//create Mach input the PED file
			File dir2_4th = new File(panelFolderName, "Imputation_MACHInput");
			if(!dir2_4th.exists()) dir2_4th.mkdirs();
			File machInput = new File(dir2_4th, "GenotypeChrom_" + chrLabel + "_Rep" + repID + ".ped.gz");
			writePEDFile = new FileHandler(machInput);
			writePEDFile.openForWrite(true);
			
			//create an all-ctrl hap ref for running minimac
			List<Integer> allCtrls = new ArrayList<Integer>();
			for(int k=0; k<ctrlCt; k++)
				allCtrls.add((Integer)k);
	  		Collections.shuffle(allCtrls);
			Set<Integer> selectedCtrls = new HashSet<Integer>();
			for(int i = 0; selectedCtrls.size() < ParamParser.num_refSub; i++)
				selectedCtrls.add(allCtrls.get(i));
			refCtrls = selectedCtrls.toArray(new Integer[selectedCtrls.size()]);
			Arrays.sort(refCtrls);
			allCtrls = null;
			selectedCtrls = null;			
			dir2_4th = new File(panelFolderName, "Imputation_MiniMacOutput_RefCtrl");
			if(!dir2_4th.exists()) dir2_4th.mkdirs();
			ctrlcarrierRecord = new File(dir2_4th, "refCtrlCarriers_" + chrLabel + repID + ".txt");
			File ctrlRef = new File(dir2_4th, "refHap_allControls" + chrLabel + "_Rep"  + repID + ".gz");
			recordRefCtrlCarriers = new FileHandler(ctrlcarrierRecord);
			recordRefCtrlCarriers.openForWrite(true);
			writeControlRef = new FileHandler(ctrlRef);
			writeControlRef.openForWrite(true);
		}
		int num_recombPerChrom = (int)(ParamParser.num_forwardGen[gID]*ParamParser.prob_recombPerChrom);
		FileHandler readSubjectFile = new FileHandler(panelFolderName, "CaseControl_Rep" + repID + ".csv");
		readSubjectFile.openForRead();
		readSubjectFile.nextLine();//first line is the header
		nextLine = null;
		int subjectCt = 0;
		int refsubCt = 0;
		while((nextLine = readSubjectFile.nextLine()) != null && subjectCt < numOfSubjects){
			boolean pickAsRef = false;
			String[] s = nextLine.split(",");
			Integer[] subChromID = new Integer[2];
			subChromID[0] = Integer.parseInt(s[1]);
			subChromID[1] = Integer.parseInt(s[2]);
			List<String> refHaps = null;
			String genotypeForMachPED = "";	
			if((planGWAS.equalsIgnoreCase("t")) && (refsubCt < ParamParser.num_refSub) && 
					(Integer.valueOf(s[3]) == 0) && (Arrays.binarySearch(refCtrls, subjectCt - caseCt) >= 0)){
				pickAsRef = true;
				refHaps = new ArrayList<String>(2);
				for(int i = 0; i < 2; i++){
					refHaps.add("S_" + refsubCt + "->S_" + refsubCt + " HAP" + (i+1) + " ");
				}
				refsubCt++;
			}
			int[][] childGeno = new int[2][];
			for(int i = 0; i < 2; i++)
				childGeno[i] = generateChildHap(cv, subChromID[i], chromsWithSNP.contains(subChromID[i]), num_recombPerChrom, refHaps, i);
			for(int j = 0; j < tag20MIndex.length; j++){
				genotypes[j][subjectCt] = childGeno[0][j] - 1 + childGeno[1][j] - 1;//Alleles are represented by 1(minor) and 2 (major) in Dan's data set; here they are converted to 0 and 1 respectively.
				if(planGWAS.equalsIgnoreCase("t")){
					genotypeForMachPED += " " + childGeno[0][j] + "/" + childGeno[1][j];
				}
			}
			if(planGWAS.equalsIgnoreCase("t")){
				writePEDFile.writeString("S_" + subjectCt + " S_" + subjectCt + " 0 0 2" + genotypeForMachPED);
				if(pickAsRef && refHaps != null){
					for(int i = 0; i < 2; i++)
						writeControlRef.writeString(refHaps.get(i));
					if(Integer.valueOf(s[4]) == 1)//this control is a cv carrier
						recordRefCtrlCarriers.writeString(nextLine);
					refHaps.clear();
				}
			}
			subjectCt++;
		}
		if(planGWAS.equalsIgnoreCase("t")){
			writePEDFile.closeWriter();
			writeControlRef.closeWriter();
			recordRefCtrlCarriers.closeWriter();
			if(ctrlcarrierRecord.length()<=0) ctrlcarrierRecord.delete();
			allSNPPosToIDMap.clear();
		}
		String subjectLine = "Case0";
		String dxLine = "2";
		String sexLine = "2";//2-Female, 1-Male, 0-unknown
		String includeLine = "1";	
		for(int i = 1; i < caseCt; i++){
			subjectLine += ",Case" + i;
			dxLine += ",2";
			sexLine += ",2";
			includeLine += ",1";
		}
		for(int i = 0; i < ctrlCt; i++){
			subjectLine += ",Control" + i;
			dxLine += ",1";
			sexLine += ",2";
			includeLine += ",1";
		}
		//for(Integer gID : genToGenerate){
			String chatInputFileName = "GenotypeChrom_" + chrLabel + ".gz";		
			FileHandler writeCHATinput = new FileHandler(dir_4th.getAbsolutePath(), chatInputFileName);
			writeCHATinput.openForWrite(true);
			writeCHATinput.writeString(subjectLine);
			writeCHATinput.writeString(dxLine);
			writeCHATinput.writeString(sexLine);
			writeCHATinput.writeString(includeLine);
			for(int j = 0; j < tag20MIndex.length; j++){
				String aMarker = simChromName + "," +  tag20MIndex[j] + ",rs" + j + ",-9999.00000000,1,";
				//String aMarker = simChrom + "," +  markerIndex[j] + ",rs" + j + ",-9999.00000000,1,";//for testing
				for(int k=0; k < numOfSubjects; k++)
					aMarker += String.valueOf(genotypes[j][k]);
				writeCHATinput.writeString(aMarker);	
			}
			writeCHATinput.preventFutureWriting();
			writeCHATinput.closeWriter();
		//}			
	}	

	private int[] generateChildHap(int cv, int hapChromID, boolean carryCV, int num_recombPerChrom, List<String> refHaps, int hapID) {
		int[] childHap = chromToMarkers.get(hapChromID);
		if(refHaps != null){
			char[] hapSeq = Util.getFullHap(hapChromID);
			refHaps.add(hapID, refHaps.remove(hapID) + String.valueOf(hapSeq));
		}
		if(num_recombPerChrom==0)	
			return childHap;
		Set<Integer> recombIndex = new HashSet<Integer>();
		boolean[] segWithCV = new boolean[num_recombPerChrom + 1];
		recombIndex.add(0);//so that the first element of the sorted list is always 0, the minimal index
		recombIndex.add(tag2MIndex.length-1);//so that the last element of the sorted list is always the maximal index
		while(recombIndex.size() < num_recombPerChrom + 2)
			recombIndex.add(getClosestTagIndex(Math.round(recomb.pickRecombLoc()*ParamParser.chromLength)));
		List<Integer> sortRecombIndex = new ArrayList<Integer>(recombIndex);
		Collections.sort(sortRecombIndex);
		//mark the recomb-site-delineated segment that carries the cv for both carriers and non-carriers; this segment will be replaced in non-carriers but not in carriers
		boolean findSegWithCV = false;
		for(int i = 0; i < sortRecombIndex.size()-1 && !findSegWithCV; i++){
			if(cv >= tag20MIndex[sortRecombIndex.get(i)] && cv < tag20MIndex[sortRecombIndex.get(i+1)]){
				segWithCV[i] = true;
				findSegWithCV = true;
				break;
			}
		}
		if(!findSegWithCV)//cv falls exactly at the last tag SNP; may happen but extremely low probability
			segWithCV[segWithCV.length-1] = true;
		for(int j = 0; j < segWithCV.length; j++){
			if(carryCV && segWithCV[j]) continue;
			int segStartID = sortRecombIndex.get(j);
			int segEndID = sortRecombIndex.get(j+1);
			if(segWithCV[j])//At this point, this condition is true only when the hap came from a non-carrier, so we need to replace the segment with that from a sample chrom without cv
				childHap = replaceSegWithStartWithoutEnd(childHap, segStartID, segEndID, chromsWithoutCV, refHaps, hapID);
			else//Otherwise, the segment does not contain the cv, so we can replace it with that from any sample chrom in the population
				childHap = replaceSegWithStartWithoutEnd(childHap, segStartID, segEndID, null, refHaps, hapID);
		}
		return childHap;
	}

	private int[] replaceSegWithStartWithoutEnd(int[] childHap, Integer segStartTagIndex, Integer segEndTagIndex, List<Integer> chromsWithoutCV,
			List<String> refHaps, int hapID) {
		int sourceChromID;
		if(chromsWithoutCV == null)
			sourceChromID = Util.randInt(0,ParamParser.sampleSize-1);
		else
			sourceChromID = chromsWithoutCV.get(Util.randInt(0,chromsWithoutCV.size()-1));
		int[] sourceTagHap = Util.gatherMutMarkersfromDistFiles(sourceChromID, tag2MIndex);
		 
		for(int i = segStartTagIndex; i < segEndTagIndex; i++)
			childHap[i] = sourceTagHap[i];
		
		if(refHaps!=null){
			char[] sourceFullHap = Util.getFullHap(sourceChromID);
			String[] tokens = refHaps.get(hapID).split("\\s+");
			char[] hapSeq = tokens[2].toCharArray();
			for(int j = allSNPPosToIDMap.get(tag20MIndex[segStartTagIndex]); j < allSNPPosToIDMap.get(tag20MIndex[segEndTagIndex]); j++)
				hapSeq[j] = sourceFullHap[j];
			refHaps.remove(hapID);
			refHaps.add(hapID, tokens[0] + " " + tokens[1] + " " + String.valueOf(hapSeq));
			sourceFullHap = null;
			hapSeq = null;
		}
		sourceTagHap = null;	
		return childHap;
	}

	
	private int getClosestTagIndex(long recombPos) {
		long minDistance = Math.abs(tag20MIndex[0] - recombPos);
		int idx = 0;
		for(int i = 1; i < tag20MIndex.length; i++){
			long newDistance = Math.abs(tag20MIndex[i] - recombPos);
			if(newDistance < minDistance){
				idx = i;
				minDistance = newDistance;
			}
		}
		return idx;
	}	
	
	private void processRecombMap(String line) {
		FileHandler recombMapFile = new FileHandler(line.trim());
		recombMapFile.openForRead();
		String nextLine = null;
		while((nextLine = recombMapFile.nextLine())!=null){
			String[] result = nextLine.split("\\s+");
			int start = Integer.parseInt(result[0]);
			double rate = Double.parseDouble(result[1]);
			recomb.addRecombSiteLL(start, rate);
		}
		recombMapFile.closeReader();
		recomb.recomb_calc_r();	
	}
}
