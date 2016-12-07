package simGWAS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;
import java.util.Map.Entry;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

import haploview.tagger.HaploData;
import haploview.tagger.SNP;
import haploview.tagger.TaggerController;
import haploview.tagger.Util;

public class CaseControlSampler /*implements Runnable*/ {
	private int snp;
	private double snpMaf;
	private int caseIndex;
	private int ctrlIndex;
	private int grrIndex;
	private int bdrIndex;
	private Collection<Integer> chromsWithCV;
	
	public CaseControlSampler(int aSNP, double aMAF, int caseIndex, int ctrlIndex, int grrIndex, int bdrIndex, Collection<Integer> snpToChrom) {
		this(caseIndex,ctrlIndex, grrIndex, bdrIndex, aMAF);
		snp = aSNP;
		chromsWithCV = snpToChrom;
	}

	public CaseControlSampler(int caseIndex, int ctrlIndex, int grrIndex, int bdrIndex, double snpPopMAF){
		snpMaf = snpPopMAF;
		this.caseIndex = caseIndex;
		this.ctrlIndex = ctrlIndex;
		this.grrIndex = grrIndex;
		this.bdrIndex = bdrIndex;
	}
	
	/*public void createAllControlRef(int replicateID){		
		File dir_2nd = Util.returnSNPDir(snp, snpMaf);
		String dir_3rd_name = "Case" + ParamParser.num_ca[ccIndex] + "_Control" + ParamParser.num_co[ccIndex] + "_GRR" + ParamParser.grr[grrIndex];
		File dir_3rd = new File(dir_2nd, dir_3rd_name);
		List<Integer> allCtrls = new ArrayList<Integer>();
		for(int k=0; k<ParamParser.num_co[ccIndex]; k++)
			allCtrls.add((Integer)k);
  		Collections.shuffle(allCtrls);
		Set<Integer> selectedCtrls = new HashSet<Integer>();
		for(int i = 0; selectedCtrls.size() < ParamParser.num_refSub; i++)
			selectedCtrls.add(allCtrls.get(i));
		Integer[] refCtrls = selectedCtrls.toArray(new Integer[selectedCtrls.size()]);
		Arrays.sort(refCtrls);
		allCtrls = null;
		selectedCtrls = null;
		FileHandler recordRefCtrlCarriers = new FileHandler(dir_3rd.getAbsolutePath()+ File.separator + "Imputation_MiniMacOutput_RefCtrl" + File.separator + "refCtrlCarriers.txt");
		recordRefCtrlCarriers.openForWrite(true);
		FileHandler readSubjectFile = new FileHandler(dir_3rd.getAbsolutePath(), "CaseControl_Rep" + replicateID + ".csv");
		readSubjectFile.openForRead();
		readSubjectFile.nextLine();//first line is the header
		Set<Integer> refHaps = new HashSet<Integer>();
		String nextLine = null;
		int ctrlCt = 0;
		boolean refCtrlIncludeCarrier = false;
		while((nextLine = readSubjectFile.nextLine()) != null){
			String[] s = nextLine.split(",");
			if(s[0].startsWith("Case")) continue;
			if(Arrays.binarySearch(refCtrls, ctrlCt) >= 0){
				refHaps.add(Integer.valueOf(s[1]));
				refHaps.add(Integer.valueOf(s[2]));
				if(Integer.valueOf(s[3]) == 1)
					System.out.println("Something's wrong. This subject is a case.\n");
				if(Integer.valueOf(s[4]) == 1){//this control is a cv carrier
					recordRefCtrlCarriers.writeString(nextLine);
					refCtrlIncludeCarrier = true;
				}
			}
			ctrlCt++;
		}
		readSubjectFile.closeReader();
		File ctrlcarrierRecord = recordRefCtrlCarriers.getFile();
		recordRefCtrlCarriers.closeWriter();
		if(!refCtrlIncludeCarrier) ctrlcarrierRecord.delete();
		Map<Integer, String> hapList = Util.getHapsfromDistFiles(refHaps.toArray(new Integer[refHaps.size()]), null);
		FileHandler writeControlRef = new FileHandler(dir_3rd.getAbsolutePath()+ File.separator + "Imputation_MiniMacOutput_RefCtrl" + File.separator + "refHap_allControls_Rep" + replicateID + ".gz");
		writeControlRef.openForWrite(true);
		int sampleHapCt = 0;
		for(Entry<Integer, String> hap : hapList.entrySet()){			
			String ref = "S_" + sampleHapCt/2 + "->S_" + sampleHapCt/2;
			if(sampleHapCt < 2)
				ref += (sampleHapCt==0)?" HAP1 ":" HAP2 ";
			else
				ref += (sampleHapCt%2==0)?" HAP1 ":" HAP2 ";
			ref += hap.getValue();
			writeControlRef.writeString(ref);
			sampleHapCt++;
			ref = null;
		}
		writeControlRef.closeWriter();
	}*/

	public void run() {
      	// At the 2nd level directory, create a folder for the current SNP and name the folder as:
		// SNP_pos<its position>_freqLo<lower frequency bound>_freqHi<higher frequency bound>
		File dir_2nd = Util.returnSNPDir(snp, snpMaf);
		if(!dir_2nd.exists()) dir_2nd.mkdirs();
		
		/*At the 3rd level directory, create a folder for the current case-control panel and name the folder as:
		 * Case<numOfCase>_Control<numOfControl>_grr<valueOfGRR>*/
		String dir_3rd_name = "Case" + ParamParser.num_ca[caseIndex] + "_Control" + ParamParser.num_co[ctrlIndex] + 
				"_GRR" + ParamParser.grr[grrIndex] + "_BDR" + ParamParser.bdr[bdrIndex];
		File dir_3rd = new File(dir_2nd, dir_3rd_name);
		if(!dir_3rd.exists()) dir_3rd.mkdirs();
		
		/* The following code outputs all data related to the first replicate of a case-control panel.
		 * The chromosomes selected for different replicates of the case-control panel is stored in a file called "ChromsUsed***.txt"
		 * inside the folder for that specific panel and will be used later to figure out the SNP markers for that specific panel. */		
		File panelChromsFile = new File(dir_3rd, "PanelUsedChroms.txt");
		FileHandler recordPanelSpecificChroms = null;
		if(!panelChromsFile.exists()){
			recordPanelSpecificChroms = new FileHandler(panelChromsFile);
			recordPanelSpecificChroms.openForWrite(true);
			recordPanelSpecificChroms.writeString("ReplicateID\tIDs of selected chromosomes");
		}
		else{
			recordPanelSpecificChroms = new FileHandler(panelChromsFile);
			recordPanelSpecificChroms.allowFutureWriting();
			recordPanelSpecificChroms.openForWrite(true);		
		}
				
  		//create case-control panel data for one specific causal variant				
  		//We will use the chromosome sample WITH replacement to assemble cases and controls
		List<Integer> chromsWithMut = new ArrayList<Integer>(chromsWithCV);
		List<Integer> chromsWithoutMut = new ArrayList<Integer>();
  		for(int k=0; k<ParamParser.sampleSize; k++){
  			if(!chromsWithMut.contains((Integer)k))
  				chromsWithoutMut.add((Integer)k);
  		}
  		
  		List<Subject> controlList = new ArrayList<Subject>();
	  	List<Subject> caseList = new ArrayList<Subject>();
	  	Set<Integer> usedChroms = new HashSet<Integer>();
	  	String usedChromSet;
  		File firstRep = new File(dir_3rd.getAbsolutePath(), "CaseControl_Rep0.csv");
  		if(!firstRep.exists()){//we may have simulated some replicates already
	  		assembCaseControlPanel(chromsWithMut, chromsWithoutMut, caseList, controlList, usedChroms);
 			usedChromSet = "";
			for(Integer chrom : usedChroms)
				usedChromSet += " " + chrom;
			/*  Here we choose to construct a string (usedChromSet) instead of using usedChroms.toString() because the latter
			 *  would output nasty things like [1, 2, 3]. The string representation consists of a list of the array's elements,
			 *  enclosed in square brackets ("[]"). Adjacent elements are separated by the characters ", " (a comma followed by a space).*/
			outputPanelData(dir_3rd.getAbsolutePath(), 0, caseList, controlList);
			recordPanelSpecificChroms.writeString("0" + usedChromSet);	
  		}
		
		//output the rest of replicate data
		for(int j=1; j<ParamParser.replica; j++){
			File curRep = new File(dir_3rd.getAbsolutePath(), "CaseControl_Rep" + j + ".csv");
			if(!curRep.exists()){
				caseList.clear();
				controlList.clear();
				usedChroms.clear();
				assembCaseControlPanel(chromsWithMut, chromsWithoutMut, caseList, controlList, usedChroms);
				usedChromSet = "";
				for(Integer chrom : usedChroms)
					usedChromSet += " " + chrom;
		    	outputPanelData(dir_3rd.getAbsolutePath(), j, caseList, controlList);
		      	recordPanelSpecificChroms.writeString((j) + usedChromSet);
			}
		}
		recordPanelSpecificChroms.preventFutureWriting();
		recordPanelSpecificChroms.closeWriter();
	}
	
	private static void outputPanelData(String directory, int replicaID, List<Subject> caseList, List<Subject> controlList) {
		FileHandler outputSubject = new FileHandler(directory, "CaseControl_Rep" + replicaID + ".csv");
		outputSubject.openForWrite(true);
		outputSubject.writeSubjectFile(caseList, controlList);
		//String fName = outputSubject.getFilename();//must capture the file before close the writer; otherwise it would be nullified
		outputSubject.preventFutureWriting();
		outputSubject.closeWriter();
		//return fName;
	}
	
	public int getNextIndex(List<Integer> listToCheck, int curIndex){
		if(curIndex == listToCheck.size()){
			Collections.shuffle(listToCheck);
			return 0;
		}
		return ++curIndex;
	}
	
	/**
	 * Controlled subjects are defined as explicitly chosen for not having the disease under study
	 * @param chromsWithMut
	 * @param chromsWithoutMut
	 * @param usedChroms 
	 * @return
	 */
	private void assembCaseControlPanel(List<Integer>chromsWithMut, List<Integer>chromsWithoutMut, 
			List<Subject> caseList, List<Subject> controlList, Set<Integer> usedChroms){
		Random randomGenerator = new Random();
		int num_caseAndCarrier = (int)(ParamParser.num_ca[caseIndex]*ParamParser.bdr[bdrIndex]);
		//int num_caseNotCarrier = ParamParser.num_ca[ccIndex] - num_caseAndCarrier;
		int num_ctrlAndCarrier = (randomGenerator.nextDouble() < (ParamParser.num_co[ctrlIndex]*ParamParser.bdr[bdrIndex]/ParamParser.grr[grrIndex]))?1:0;
		int num_ctrlNotCarrier = ParamParser.num_co[ctrlIndex] - num_ctrlAndCarrier;
		
		//calculate the proportion of homozygotes in cases that are carriers
		double prop_homoInCarrier = snpMaf*snpMaf/(snpMaf*snpMaf + 2*snpMaf*(1-snpMaf));
		
		int caseCt = 0;
		int controlCt = 0;
		int numOfChromWithMut = chromsWithMut.size();
		int numOfChromWithoutMut = chromsWithoutMut.size();
		
		Collections.shuffle(chromsWithMut);
		Collections.shuffle(chromsWithoutMut);
		//int curMutChromIndex = -1;
		//int curWildChromIndex = -1;
		Subject newOne = null;
		randomGenerator = new Random();
		//select cases from carriers WITH REPLACEMENT (sometimes when the causal variant is very rare there may be only one or two chromosomes having it)	
		while(caseCt < num_caseAndCarrier){
			if(randomGenerator.nextDouble() < prop_homoInCarrier){
				newOne = new Subject("Case" + caseCt, chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)),
						chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)), 1);
				/*newOne = new Subject("Case" + caseCt, chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)),
						chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)), 1);*/
				caseList.add(newOne);
			}
			else{
				newOne =new Subject("Case" + caseCt, chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)),
						chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)), 1);
				/*newOne =new Subject("Case" + caseCt, chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)),
						chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)), 1);*/
				caseList.add(newOne);	
			}
			caseCt++;
			usedChroms.add(new Integer(newOne.getChromA()));
			usedChroms.add(new Integer(newOne.getChromB()));
		}
		
		//select cases from non-carriers WITH REPLACEMENT
		while(caseCt < ParamParser.num_ca[caseIndex]){
			newOne = new Subject("Case" + caseCt, chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)),
					chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)), 0);
			/*newOne = new Subject("Case" + caseCt, chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)),
					chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)), 0);*/
			caseList.add(newOne);
			caseCt++;
			usedChroms.add(new Integer(newOne.getChromA()));
			usedChroms.add(new Integer(newOne.getChromB()));
		}
		
		//select controls from non-carriers WITH REPLACEMENT		
		while(controlCt < num_ctrlNotCarrier){
			newOne = new Subject("Control" + controlCt, chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)),
				chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)), 0);
			/*newOne = new Subject("Control" + controlCt, chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)),
					chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)), 0);*/
			controlList.add(newOne);
			controlCt++;
			usedChroms.add(new Integer(newOne.getChromA()));
			usedChroms.add(new Integer(newOne.getChromB()));
		}
		
		//select controls from carriers WITH REPLACEMENT		
		while(controlCt < ParamParser.num_co[ctrlIndex]){
			if(randomGenerator.nextDouble() < prop_homoInCarrier){
				newOne = new Subject("Control" + controlCt, chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)),
						chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)), 1);
				/*newOne = new Subject("Control" + controlCt, chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)),
						chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)), 1);*/
				controlList.add(newOne);
			}
			else {
				newOne = new Subject("Control" + controlCt, chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)),
						chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)), 1);
				/*newOne = new Subject("Control" + controlCt, chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)),
						chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)), 1);*/
				controlList.add(newOne);	
			}
			controlCt++;
			usedChroms.add(new Integer(newOne.getChromA()));
			usedChroms.add(new Integer(newOne.getChromB()));
		}		
	}
	
	//For the sake of simplicity, no requirements on the number of homozygotes
	public void assembCaseControlPanel(List<Carrier> carriers, List<NonCarrier> nonCarriers, List<char[][]> caseList, List<char[][]> controlList, 
			int num_caseAndCarrier, int num_ctrlNotCarrier){		
		int caseCt = 0;
		int controlCt = 0;	
		Collections.shuffle(carriers);
		Collections.shuffle(nonCarriers);
		int curCarrierIndex = 0;
		int curNonCarrierIndex = 0;
		while(caseCt < num_caseAndCarrier){
			if(curCarrierIndex == carriers.size()){
				Collections.shuffle(carriers);
				curCarrierIndex = 0;
			}
			char[][] diploid = new char[2][ParamParser.num_mark];//must be a different object
			System.arraycopy(carriers.get(curCarrierIndex).getSingleChrom(0), 1, diploid[0], 0, ParamParser.num_mark);
			System.arraycopy(carriers.get(curCarrierIndex).getSingleChrom(1), 1, diploid[1], 0, ParamParser.num_mark);
			caseList.add(diploid);
			curCarrierIndex++;
			caseCt++;
		}
		while(caseCt < ParamParser.num_ca[caseIndex]){
			if(curNonCarrierIndex == nonCarriers.size()){
				Collections.shuffle(nonCarriers);
				curNonCarrierIndex = 0;
			}
			char[][] diploid = new char[2][ParamParser.num_mark];
			diploid[0] = nonCarriers.get(curNonCarrierIndex).getSingleChrom(0).clone();
			diploid[1] = nonCarriers.get(curNonCarrierIndex).getSingleChrom(1).clone();
			caseList.add(diploid);
			curNonCarrierIndex++;
			caseCt++;
		}				
		while(controlCt < num_ctrlNotCarrier){
			if(curNonCarrierIndex == nonCarriers.size()){
				Collections.shuffle(nonCarriers);
				curNonCarrierIndex = 0;
			}
			char[][] diploid = new char[2][ParamParser.num_mark];
			diploid[0] = nonCarriers.get(curNonCarrierIndex).getSingleChrom(0).clone();
			diploid[1] = nonCarriers.get(curNonCarrierIndex).getSingleChrom(1).clone();
			controlList.add(diploid);
			curNonCarrierIndex++;
			controlCt++;
		}				
		while(controlCt < ParamParser.num_co[ctrlIndex]){
			if(curCarrierIndex == carriers.size()){
				Collections.shuffle(carriers);
				curCarrierIndex = 0;
			}
			char[][] diploid = new char[2][ParamParser.num_mark];
			System.arraycopy(carriers.get(curCarrierIndex).getSingleChrom(0), 1, diploid[0], 0, ParamParser.num_mark);
			System.arraycopy(carriers.get(curCarrierIndex).getSingleChrom(1), 1, diploid[1], 0, ParamParser.num_mark);
			controlList.add(diploid);
			curCarrierIndex++;
			controlCt++;
		}		
	}
}
