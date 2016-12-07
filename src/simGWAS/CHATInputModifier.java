package simGWAS;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import haploview.tagger.Util;

public class CHATInputModifier extends ChatInputGenerator {
	private int cv;
	private int repID;
	private int curCarrierIndex = 0;
	private int curNonCarrierIndex = 0;
	private List<Carrier> carriers = null;
	private List<NonCarrier> nonCarriers = null;
	RecombWorker recomb = new RecombWorker();

	public CHATInputModifier(int cv, Map<Integer, Integer> tagSNPs, File panelSpecificChroms, int caseIndex, int ctrlIndex, int replicate, int simChrom, String recombFile) {
		super(tagSNPs, panelSpecificChroms, caseIndex, ctrlIndex, simChrom);
		repID = replicate;
		this.cv = cv;
		processRecombMap(recombFile);
	}
	
	public int getCurCarrierIndex(){
		if(curCarrierIndex == carriers.size()){
			Collections.shuffle(carriers);
			curCarrierIndex = 0;
		}
		return curCarrierIndex++;
	}
	
	public int getCurNonCarrierIndex(){
		if(curNonCarrierIndex == nonCarriers.size()){
			Collections.shuffle(nonCarriers);
			curNonCarrierIndex = 0;
		}
		return curNonCarrierIndex++;
	}

	public void run(Collection<Integer> chromWithSNP, int num_Gen, int grrIndex, int bdrIndex, double snpPopMAF) {
		FileHandler readPanelSpecificChroms = new FileHandler(chromFile);
		readPanelSpecificChroms.openForRead();
		readPanelSpecificChroms.nextLine();//first line is the header
		
		File dir_4th = new File(panelFolderName + File.separator + "ChatInput_tagSNP" + File.separator + "Rep" + repID);
		if(!dir_4th.exists()) dir_4th.mkdirs();
		
		String nextLine;
		Map<Integer, int[]> chromToMarkers = new HashMap<Integer, int[]>();
		while((nextLine = readPanelSpecificChroms.nextLine()) != null){
			String[] s = nextLine.split(" ");					
			int replicateID = Integer.valueOf(s[0]);
			if(replicateID != repID) continue;
			Integer[] panelSpecificChroms = new Integer[s.length-1];
			for(int j = 1; j < s.length; j++)
				panelSpecificChroms[j-1] = Integer.valueOf(s[j]);
			chromToMarkers = Util.gatherMutMarkersfromDistFiles(panelSpecificChroms, tag2MIndex);
			s = null;
			break;
		}
		readPanelSpecificChroms.closeReader();
		FileHandler readSubjectFile = new FileHandler(panelFolderName, "CaseControl_Rep" + repID + ".csv");
		readSubjectFile.openForRead();
		readSubjectFile.nextLine();//first line is the header
		List<Carrier> newCarriers = null; 
		List<NonCarrier> newNonCarriers = null;
		int cvChromCt = 0;
		int newCVChromCt = 0;
		int outputGen = 1;
		for(int i = 0; i < num_Gen; i++){
			if(i==0){
				carriers = new ArrayList<Carrier>();
				nonCarriers = new ArrayList<NonCarrier>();
				nextLine = null;
				while((nextLine = readSubjectFile.nextLine()) != null){
					String[] s = nextLine.split(",");
					boolean chr0WithCV = false;
					boolean chr1WithCV = false;
					int chromA = Integer.parseInt(s[1]);
					if(chromWithSNP.contains(chromA)){
						chr0WithCV = true;
						cvChromCt++;
					}
					int chromB = Integer.parseInt(s[2]);
					if(chromWithSNP.contains(chromB)){
						chr1WithCV = true;
						cvChromCt++;
					}
					int[][] diploid = new int[2][];
					diploid[0] = chromToMarkers.get(chromA);
					diploid[1] = chromToMarkers.get(chromB);
					if(Integer.parseInt(s[4])==1){//Carrier
						Carrier aSub = new Carrier(diploid, chr0WithCV, chr1WithCV);
						carriers.add(aSub);
					}
					else{
						NonCarrier aSub = new NonCarrier(diploid);
						nonCarriers.add(aSub);
					}
				}
			}
			else{
				carriers = newCarriers;
				nonCarriers = newNonCarriers;
				cvChromCt = newCVChromCt;
				newCVChromCt = 0;
			}		
			Collections.shuffle(carriers);
			Collections.shuffle(nonCarriers);
			newCarriers = new ArrayList<Carrier>();
			newNonCarriers = new ArrayList<NonCarrier>();
			int new_num_caseAndCarrier = (int)(ParamParser.num_ca[caseIndex]*ParamParser.bdr[bdrIndex]);
			int new_num_ctrlAndCarrier = ((new Random()).nextDouble() < ParamParser.num_co[ctrlIndex]*ParamParser.bdr[bdrIndex]/ParamParser.grr[grrIndex])?1:0;
			int new_num_ctrlNotCarrier = ParamParser.num_co[ctrlIndex] - new_num_ctrlAndCarrier;
			double prop_homoInCarrier = snpPopMAF*snpPopMAF/(snpPopMAF*snpPopMAF + 2*snpPopMAF*(1-snpPopMAF));
			while(newCVChromCt < cvChromCt){
				Carrier child = generateCarrierOffSpring(carriers.get(getCurCarrierIndex()), carriers.get(getCurCarrierIndex()), cv);
				if(child.isHomozygoteWithCV()){					
					if((new Random()).nextDouble() < prop_homoInCarrier){
						newCarriers.add(child);
						newCVChromCt += 2;
					}
					else continue;
				}
				else if(child.isCarrier()){
					newCarriers.add(child);
					newCVChromCt++;
				}
			}
			int nonCarrierCt = ParamParser.num_ca[caseIndex] + ParamParser.num_co[ctrlIndex] - newCarriers.size();
			while(newNonCarriers.size() < nonCarrierCt)
				newNonCarriers.add(generateNonCarrierOffSpring(nonCarriers.get(getCurNonCarrierIndex()), nonCarriers.get(getCurNonCarrierIndex())));
			CaseControlSampler cc = new CaseControlSampler(this.caseIndex, this.ctrlIndex, grrIndex, bdrIndex,snpPopMAF);
			List<char[][]> cases = new ArrayList<char[][]>();
			List<char[][]> controls = new ArrayList<char[][]>();
			cc.assembCaseControlPanel(newCarriers, newNonCarriers, cases, controls, new_num_caseAndCarrier, new_num_ctrlNotCarrier);
			if((i+1)%20==0) outputCHATInputFile(dir_4th.getAbsolutePath(), outputGen++, cases, controls);
		}	
	}
	
	private void outputCHATInputFile(String absolutePath, int forwardGen, List<char[][]> cases, List<char[][]> controls){    
		    String simChrom = String.valueOf(forwardGen);//use a fake chromosome name
			String outputFileName = "GenotypeChrom_" + padNumberAsString(forwardGen,2) + ".gz";
			File chatInputFile = new File(absolutePath + File.separator + "ChromosomeSpecificData" + File.separator + outputFileName);
			int numOfSubjects = cases.size() + controls.size();
			if(chatInputFile.exists() && sanityCheck(tag20MIndex, chatInputFile, cases.size() + controls.size(), true)) return;
			else if(chatInputFile.exists())
				chatInputFile.delete();
			String subjectLine = "Case0";
			String dxLine = "2";
			String sexLine = "0";//2-Female, 1-Male, 0-unknown
			String includeLine = "1";
			int numOfMarkers = tag20MIndex.length;
			int[][] genotypes = new int[numOfMarkers][numOfSubjects];
			for(int i = 0; i < cases.size(); i++){
				if(i > 0){
					subjectLine += ",Case" + i;
					dxLine += ",2";
					sexLine += ",0";
					includeLine += ",1";
				}
				for(int j = 0; j < numOfMarkers; j++){
					genotypes[j][i] = Integer.valueOf(String.valueOf(cases.get(i)[0][j])) - 1 + Integer.valueOf(String.valueOf(cases.get(i)[1][j])) - 1;//Alleles are represented by 1(minor) and 2 (major) in Dan's data set; here they are converted to 0 and 1 respectively.
					/*System.out.println("case" + i + "chr0:" + String.valueOf(cases.get(i)[0]) + "\n");
					System.out.println("case" + i + "chr1:" + String.valueOf(cases.get(i)[1]) + "\n");
					if(cases.get(i)[0][j] != cases.get(i)[1][j]){
						System.out.println("case" + i + "chr0 at position " + j + ": " + cases.get(i)[0][j] + "\n");
						System.out.println("case" + i + "chr1 at position " + j + ": " + cases.get(i)[1][j] + "\n");
						System.out.println("case" + i + "genotype at position " + j + ": " + genotypes[j][i] + "\n");
					}*/
				}
			}
			for(int i = 0; i < controls.size(); i++){
				subjectLine += ",Control" + i;
				dxLine += ",1";
				sexLine += ",0";
				includeLine += ",1";
				for(int j = 0; j < numOfMarkers; j++){					
					genotypes[j][i+cases.size()] = Integer.valueOf(String.valueOf(controls.get(i)[0][j])) - 1 + Integer.valueOf(String.valueOf(controls.get(i)[1][j])) - 1;
					/*System.out.println("control" + i + "chr0:" + String.valueOf(controls.get(i)[0]) + "\n");
					System.out.println("control" + i + "chr1:" + String.valueOf(controls.get(i)[1]) + "\n");
					if(controls.get(i)[0][j] != controls.get(i)[1][j]){
						System.out.println("case" + i + "chr0 at position " + j + ": " + controls.get(i)[0][j] + "\n");
						System.out.println("case" + i + "chr1 at position " + j + ": " + controls.get(i)[1][j] + "\n");
						System.out.println("case" + i + "genotype at position " + j + ": " + genotypes[j][i+cases.size()] + "\n");
					}*/
				}
			}
		
			FileHandler writeCHATinput = new FileHandler(chatInputFile);
			writeCHATinput.openForWrite(true);
			writeCHATinput.writeString(subjectLine);
			writeCHATinput.writeString(dxLine);
			writeCHATinput.writeString(sexLine);
			writeCHATinput.writeString(includeLine);
			for(int j = 0; j < numOfMarkers; j++){
				String aMarker = simChrom + "," +  tag20MIndex[j] + ",rs" + j + ",-9999.00000000,1,";
				//String aMarker = simChrom + "," +  markerIndex[j] + ",rs" + j + ",-9999.00000000,1,";//for testing
				for(int k=0; k < numOfSubjects; k++)
					aMarker += String.valueOf(genotypes[j][k]);
				writeCHATinput.writeString(aMarker);	
			}
			writeCHATinput.preventFutureWriting();
			writeCHATinput.closeWriter();
	}

	//assume there is only one child
	public Carrier generateCarrierOffSpring(Carrier parent_p, Carrier parent_m, int cvPos){
		Carrier child = new Carrier();	
		//get paternal recombination result
		child.recordHap(0, generateParentChr(parent_p, cvPos));
		//get maternal recombination result
		child.recordHap(1, generateParentChr(parent_m, cvPos));
		return child;
	}
	
	//assume there is only one child
	public NonCarrier generateNonCarrierOffSpring(NonCarrier parent_p, NonCarrier parent_m){
		NonCarrier child = new NonCarrier();
		//get paternal recombination result
		child.recordHap(0, generateParentChr(parent_p));
		//get maternal recombination result
		child.recordHap(1, generateParentChr(parent_m));
		return child;
	}
	
	private String generateParentChr(NonCarrier oneParent) {
		Random randomGenerator = new Random();
		if(randomGenerator.nextDouble() < 0.5){
			if(randomGenerator.nextDouble() < 0.5)
				return String.valueOf(oneParent.getSingleChrom(0));
			else
				return String.valueOf(oneParent.getSingleChrom(1));
		}
		else{
			long recombPos = Math.round(recomb.pickRecombLoc()*ParamParser.chromLength);
			int recombIndex = getClosestTagIndex(recombPos);
			char[][] recombinedChroms = new char[2][ParamParser.num_mark];
			for(int i = 0; i < recombIndex; i++){
				recombinedChroms[0][i] = oneParent.getSingleChrom(0)[i];
				recombinedChroms[1][i] = oneParent.getSingleChrom(1)[i];
			}
			for(int i = recombIndex; i < ParamParser.num_mark; i++){
				recombinedChroms[0][i] = oneParent.getSingleChrom(1)[i];
				recombinedChroms[1][i] = oneParent.getSingleChrom(0)[i];
			}
			if(randomGenerator.nextDouble() < 0.5)
				return String.valueOf(recombinedChroms[0]);
			else
				return String.valueOf(recombinedChroms[1]);			
		}
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

	private String generateParentChr(Carrier oneParent, int cvPos) {
		Random randomGenerator = new Random();
		if(randomGenerator.nextDouble() < 0.5){
			if(randomGenerator.nextDouble() < 0.5){
				return String.valueOf(oneParent.getSingleChrom(0));
			}
			else{
				return String.valueOf(oneParent.getSingleChrom(1));
			}
		}
		else{
			char[] haplotypeToChild = null;
			long recombPos = Math.round(recomb.pickRecombLoc()*ParamParser.chromLength);
			int recombIndex = getClosestTagIndex(recombPos);
			char[][] recombinedChroms = new char[2][ParamParser.num_mark];
			for(int i = 0; i < recombIndex; i++){
				recombinedChroms[0][i] = oneParent.getSingleChrom(0)[i+1];
				recombinedChroms[1][i] = oneParent.getSingleChrom(1)[i+1];
			}
			for(int i = recombIndex; i < ParamParser.num_mark; i++){
				recombinedChroms[0][i] = oneParent.getSingleChrom(1)[i+1];
				recombinedChroms[1][i] = oneParent.getSingleChrom(0)[i+1];
			}
			if(randomGenerator.nextDouble() < 0.5){
				haplotypeToChild = recombinedChroms[0];
				if((oneParent.doesChromWithCV(0) && cvPos < recombPos) || (oneParent.doesChromWithCV(1) && cvPos >= recombPos))
					return "1" + String.valueOf(haplotypeToChild);
				else return "0" + String.valueOf(haplotypeToChild);
			}
			else{
				haplotypeToChild = recombinedChroms[1];
				if((oneParent.doesChromWithCV(1) && cvPos < recombPos) || (oneParent.doesChromWithCV(0) && cvPos >= recombPos))
					return "1" + String.valueOf(haplotypeToChild);
				else return "0" + String.valueOf(haplotypeToChild);
			}
		}
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
	
    public static String padNumberAsString(int val, int padLength) {
        StringBuffer outString = new StringBuffer();
        String inputString = String.valueOf(val);
        for (int x = inputString.length() + 1; x <= padLength; x++)
            outString.append("0");
        outString.append(inputString);
        return outString.toString();
    }
}
