package simGWAS;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

import haploview.tagger.*;

/**
 * A tool for creating a human disease GWAS data set using DNA sequences generated from a coalescent simulation
 * !!!!!!!!!!!!!!!!!!!!Alert: this version of code can only work on a largemem node, as it directly inputs coJava hap file
 * @author Yuan Lin
 */
public class simGWASmain {
	private static ParamParser params;
    
    //Mapping a SNP to the frequency of its minor allele; entry format: <index of a SNP regarding the all-SNP array (stored in out.snp-1.gz), its MAF regarding a specific panel> 
    private static Map<Integer, Double> minorAlleleFreq = new HashMap<Integer, Double>();
    
    //Mapping a SNP to chromosomes with mutation on that SNP; entry format: <index of a SNP regarding the all-SNP array (stored in out.snp-1.gz), chromID>
    private static Map<Integer, Collection<Integer>> snpToChrom = new HashMap<Integer, Collection<Integer>>();
    
	
	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length == 0){
            System.out.println("\nPlease indicate name of the .properties file that provides default values for program parameters.\n");
            System.exit(0);
		}
        else{
	       	// Obtain and assign parameter values
	       	params = new ParamParser(args[0]);
	       	
	       	//create the 1st level directory that holds all data output	
	    	String dir_1st_name = ParamParser.dataOutputDir /*+ File.separator + "simGWASNewData"*/;
	    	File dir_1st = new File(dir_1st_name);
	       	if(!dir_1st.exists()) dir_1st.mkdirs();
	       	
	       	//If the file that contains mutated positions (mutInput) is unavailable, use the code below to derive it from raw sequence data (hapInput) in parallel
/*	       	byte[][] hapArray = new byte[ParamParser.sampleSize][ParamParser.num_mut];
       		FileHandler readHapInput = new FileHandler(ParamParser.dataInputDir, ParamParser.hapInput);
       		readHapInput.openForRead();       			
       		String input = null;
       		int sampleCt = 0;
   			while((input = readHapInput.nextLine())!=null && sampleCt < ParamParser.sampleSize){
   				String[] s = input.split("\\s+");
   				if(s.length-2 != ParamParser.num_mut){
   					System.out.println("Sample " + sampleCt + "has incomplete sequence data.\n");
   					System.exit(0);
   				}
   				for(int i = 2; i < s.length; i++)
   					hapArray[sampleCt][i-2] = Byte.parseByte(s[i]);
   				sampleCt++;
   			}
   			if(sampleCt < ParamParser.sampleSize){
				System.out.println("There are only " + (sampleCt-1) + "samples.\n");
				System.exit(0);
			}
			readHapInput.closeReader();
			
			int processorAvailable = Runtime.getRuntime().availableProcessors();
			ExecutorService executor = Executors.newFixedThreadPool(processorAvailable);
			for(int k = 0; k < ParamParser.inputFileNum; k++){
				int start = ParamParser.sampleSize/ParamParser.inputFileNum*k;
				int end = start + ParamParser.sampleSize/ParamParser.inputFileNum - 1;
				Runnable worker = new SnpPosWriter(start,end, hapArray,dir_1st_name,ParamParser.mutInput + "_" + k + ".gz");
				executor.execute(worker);
			}
			executor.shutdown();
			while(!executor.isTerminated()){}
			System.out.println("Finish creating out.pos-# files.\n");
*/	       	
	       	/*try{
	       		String input = null;
	       		int doneLine = 0;
	       		File mutations = new File(dir_1st_name, params.mutInput);
	       		FileHandler processMutInput = new FileHandler(mutations);
	       		if(!mutations.exists()){
	       			processMutInput.openForWrite(true);
	       			processMutInput.writeString("ChromID\tMutation Position.....");//output header
				}
	       		else{//already exists; go to the end of this file; check whether it is complete; if not return the last chromID; start a new line
	       			processMutInput.openForRead();
	       			doneLine = 0;
	       			processMutInput.nextLine();//header
	       			while((input = processMutInput.nextLine())!=null){
	       				doneLine++;
       				}
	       			processMutInput.openForWrite(true);
	       			processMutInput.writer.newLine();
	       		}
	       		FileHandler readHapInput = new FileHandler(params.dataInputDir, params.hapInput);
	       		readHapInput.openForRead();       			
	       		input = null;
       			int readLine = 0;
       			while((input = readHapInput.nextLine())!=null && readLine <= 16333){
       				if(readLine < doneLine){
       					readLine++;
       					continue;
       				}
       				Scanner s = new Scanner(input);
       				s.useDelimiter("\\s+");
       				String output = s.next();//capture chrom ID
       				s.next();//skip population ID
       				int ct=0;
       				while(s.hasNext()){
       					if(s.nextInt()==1)
       						output += "\t" + ct;
       					ct++;
       				}
       				s.close();
       				processMutInput.writeString(output);
       			}
				readHapInput.closeReader();
				processMutInput.closeWriter();
	       	}catch(Exception e){
				e.printStackTrace();
			}*/
		
	       	//check selected casual variants (at the first level directory)
			File cvFile = new File(dir_1st_name, "selectedCausalVariants.txt");
			FileHandler handleCausalVar = null;
			//If the file does not exist, it means causal variants have not been selected yet. Go ahead selecting them and write relevant 
	      	//information into the file. The minorAlleleFreq map and the snpToChrom map will be created during the process. 	     
			if(!cvFile.exists()){
				handleCausalVar = new FileHandler(cvFile);
				handleCausalVar.openForWrite(true);
				handleCausalVar.writeString("CV Index\tMAF\tChromosomes with this CV");
				for(int i = 0; i < ParamParser.lo_maf.length; i++){
					selectCausalVariants(ParamParser.segmentStart, ParamParser.segmentEnd, ParamParser.lo_maf[i], ParamParser.hi_maf[i], ParamParser.numOfCVPerMAFRange);
					for(Entry<Integer, Double> var : minorAlleleFreq.entrySet()){
				      	int aSNP = var.getKey();
				      	List<Integer> chromsWithMut = new ArrayList<Integer>(snpToChrom.get(aSNP));
				      	String chromLine = "";
				      	for(Integer chrom : chromsWithMut)
				      		chromLine += "\t" + chrom;
			      		handleCausalVar.writeString(String.valueOf(aSNP)+ "\t" + String.valueOf(var.getValue()) + chromLine);
					}
	      		}
	      		handleCausalVar.closeWriter();
			}
			
	      	// If the file already exists, it means causal variants have been selected and we are working in a recovery mode.
			// Then read their info from the file to populate the minorAlleleFreq map and the snpToChrom map. 
	      	else{
	      		handleCausalVar = new FileHandler(cvFile);
	      		handleCausalVar.openForRead();
	      		handleCausalVar.nextLine();//header
				String cvRecord = null;
				while((cvRecord = handleCausalVar.nextLine())!=null){
					Scanner cv = new Scanner(cvRecord);
					cv.useDelimiter("\t");
					int aSNP = Integer.parseInt(cv.next());
					minorAlleleFreq.put(aSNP, Double.parseDouble(cv.next()));
					while(cv.hasNext()){
						Collection<Integer> chromsWithSpecificSNPMut = snpToChrom.get(aSNP);			
						if (chromsWithSpecificSNPMut == null) {
							chromsWithSpecificSNPMut = new HashSet<Integer>();
							snpToChrom.put(aSNP, chromsWithSpecificSNPMut);
						}
						chromsWithSpecificSNPMut.add(Integer.parseInt(cv.next()));
					}
					cv.close();
				}
				handleCausalVar.closeReader();
	      	}
	      	
			//int processorAvailable = Runtime.getRuntime().availableProcessors();
			//ExecutorService executor = Executors.newFixedThreadPool(processorAvailable);
			int aSNP = Integer.valueOf(args[1]);
/*			int caseIndex = Integer.valueOf(args[2]);
			int ctrlIndex = Integer.valueOf(args[3]);
			int grrIndex = Integer.valueOf(args[4]);
			int bdrIndex = Integer.valueOf(args[5]); 
			int replicate = Integer.valueOf(args[6]);
			String recombMapFile = args[7];
			Collection<Integer> chromWithSNP = snpToChrom.get(aSNP);*/
			double aMAF = minorAlleleFreq.get(aSNP);			
			minorAlleleFreq = null;
			snpToChrom = null;
			
/*			CaseControlSampler sampler = new CaseControlSampler(aSNP, aMAF, caseIndex, ctrlIndex, grrIndex, bdrIndex, chromWithSNP);
			sampler.run();*/
		    //System.out.println("Finish creating cases and controls.\n");
		    
		    //!!!!!!!!First launch one job per cv to get a set of tags. All scenarios of a cv uses the same set of tags
		    File dir_2nd = Util.returnSNPDir(aSNP, aMAF);
		    File tagFile = new File(dir_2nd.getAbsolutePath(), "markerSNPs_tag.txt");
		    Integer[] tagPos;
		    tagPos = gatherTagSNPs(tagFile, aSNP, null);
		    //System.out.println("Finish selecting tag SNPs.\n");
	
			//Create an absolute position to index map from the "out.snp-1.gz" file
			Map<Integer, Integer> posToIndex = new HashMap<Integer, Integer>();
			FileHandler readSNPPos = new FileHandler(ParamParser.dataOutputDir, ParamParser.posInput);
			readSNPPos.openForRead();
			String nextLine = null;	
			int index = 0;
			while((nextLine = readSNPPos.nextLine()) != null){
				posToIndex.put(Integer.valueOf(nextLine), index);
				index++;
			}
			readSNPPos.closeReader();
		    Map<Integer, Integer> tagSNPs = new HashMap<Integer, Integer>();
		    Map<Integer, Integer> mutCtOfTagSNP = new HashMap<Integer, Integer>();
			for(Integer tag20MPos : tagPos){
				int tag2MIndex = posToIndex.get(tag20MPos);
				tagSNPs.put(tag20MPos, tag2MIndex);
				mutCtOfTagSNP.put(tag20MPos,0);
			}
			
			//The following code obtain minor allele frequency of all tag SNPs based on the entire population
			for(int i = 0; i < ParamParser.inputFileNum; i++){
				FileHandler readMutData = new FileHandler(ParamParser.dataOutputDir, ParamParser.mutInput + "_" + i + ".gz");
				readMutData.openForRead();
				readMutData.nextLine();//first line is the header
				String currLine = null;
				while ((currLine = readMutData.nextLine()) != null){
					String[] s = currLine.split("\\s+");
					int[] mutations = new int[s.length-1];
					for(int j = 1; j < s.length; j++)//In out.pos-1_# files, s[0] is chromID
						mutations[j-1] = Integer.valueOf(s[j]);
					Arrays.sort(mutations);
					for(Integer tag20MPos : tagPos){
						if(Arrays.binarySearch(mutations, tagSNPs.get(tag20MPos)) >= 0){
							mutCtOfTagSNP.put(tag20MPos, (mutCtOfTagSNP.get(tag20MPos)+1));
							//System.out.println(tag20MPos + " with " + mutCtOfTagSNP.get(tag20MPos));
						}
					}
				}				
				readMutData.closeReader();
			}
 			File freqFile = new File(dir_2nd,"tagSNPMAF.txt");
  			if(freqFile.exists()) freqFile.delete();
  			FileHandler writeAlleleFreq = new FileHandler(freqFile);
  			writeAlleleFreq.openForWrite(true);
  			writeAlleleFreq.writeString("MarkerPos\tMinorAlleleFreq");  			
  			for(Entry<Integer, Integer> entry : mutCtOfTagSNP.entrySet()){
  				double mutAlleleFreq = ((double)entry.getValue())/ParamParser.sampleSize;
  				writeAlleleFreq.writeString(entry.getKey() + "\t" + ((mutAlleleFreq>0.5)?(1-mutAlleleFreq):mutAlleleFreq));
  			}
  			writeAlleleFreq.closeWriter();
			
/*			String dir_3rd_name = "Case" + ParamParser.num_ca[caseIndex] + "_Control" + ParamParser.num_co[ctrlIndex] + 
					"_GRR" + ParamParser.grr[grrIndex] + "_BDR" + ParamParser.bdr[bdrIndex];
			File dir_3rd = new File(dir_2nd.getAbsolutePath(), dir_3rd_name);
			File panelSpecificChroms = new File(dir_3rd, "PanelUsedChroms.txt");
			
			Collections.shuffle(tagSNPList);
			Integer[] snpWithMissingData = new Integer[(int) (ParamParser.num_mark*ParamParser.missingRate)];
			FileHandler writeMissingSNPFile = new FileHandler(dir_3rd_name, "SNPwithMissingGeno.txt");
			writeMissingSNPFile.openForWrite(true);
			int i = 0;
			for(; i < snpWithMissingData.length; i++){
				snpWithMissingData[i] = tagSNPList.get(i);
				writeMissingSNPFile.writeString(Integer.toString(snpWithMissingData[i]));
			}
			writeMissingSNPFile.closeWriter();
			Integer[] snpWithRandError = new Integer[(int) (ParamParser.num_mark*ParamParser.errorRate)];
			FileHandler writeRandErrSNPFile = new FileHandler(dir_3rd_name, "SNPwithRandGenoErr.txt");
			writeRandErrSNPFile.openForWrite(true);
			for(int j = 0; j < snpWithRandError.length; j++){
				snpWithRandError[j] = tagSNPList.get(i+j);
				writeRandErrSNPFile.writeString(Integer.toString(snpWithRandError[j]));
			}
			writeRandErrSNPFile.closeWriter();
			tagSNPList = null;
			//ChatInputGenerator chatDataWriter = new ChatInputGenerator(tagSNPs, panelSpecificChroms, caseIndex, ctrlIndex, ParamParser.founderChrom);
			//chatDataWriter.run(args[9]);
			//System.out.println("Finish generating CHAT inputs.\n");
			CHATInputModifierV2 chatInputModifier = new CHATInputModifierV2(aSNP, tagSNPs, panelSpecificChroms, caseIndex, ctrlIndex, 
					replicate, Integer.valueOf(args[8]), recombMapFile, posToIndex);
			tagSNPs.clear();
			chatInputModifier.run(chromWithSNP, grrIndex, bdrIndex, aMAF, args[9], snpWithMissingData, snpWithRandError);
*/        }
	}
	
	private static Map<Integer, Integer> gatherTagSNPs(String parentDirName, Map<Integer, Integer> posToIndex) throws IOException {
		Map<Integer, Integer> tags = new HashMap<Integer, Integer>();
	    FileHandler handleTagSNPs = new FileHandler(parentDirName, "markerSNPs_tag.txt");
	    handleTagSNPs.openForRead();
	    handleTagSNPs.nextLine();//header line
	    String currentLine = null;
	    while((currentLine = handleTagSNPs.nextLine()) != null){
	    	int tag20MPos = Integer.parseInt(currentLine);
	    	tags.put(tag20MPos, posToIndex.get(tag20MPos));	
	    }
		handleTagSNPs.closeReader();
		return tags;
	}
	
	private static Integer[] gatherTag20MPos(String parentDirName) throws IOException {
		List<Integer> tags = new ArrayList<Integer>();
	    FileHandler handleTagSNPs = new FileHandler(parentDirName, "markerSNPs_tag.txt");
	    handleTagSNPs.openForRead();
	    handleTagSNPs.nextLine();//header line
	    String currentLine = null;
	    while((currentLine = handleTagSNPs.nextLine()) != null){
	    	tags.add(Integer.parseInt(currentLine));
	    }
		handleTagSNPs.closeReader();
		Integer[] tagPos = tags.toArray(new Integer[tags.size()]);
		return tagPos;
	}
	
	private static Integer[] gatherTagSNPs(File tagFile, int aSNP, List<Integer> chromsWithoutMut) throws IOException {
		List<Integer> tags = new ArrayList<Integer>();
	    FileHandler handleTagSNPs = null;
	    //If the file does not exist, tag SNPs have not been generated yet. Go ahead generating them and write information in this file.
		if(chromsWithoutMut!=null){
			handleTagSNPs = new FileHandler(tagFile);
			tags = selectTagSNPs(chromsWithoutMut);
   			handleTagSNPs.openForWrite(true);
   			handleTagSNPs.writeString("Absolute Positions of Tag SNPs.....");//output header	
   			for(Integer tag : tags)
   				handleTagSNPs.writeString(String.valueOf(tag));
   			handleTagSNPs.writer.flush();
   			handleTagSNPs.closeWriter();
		}
		//If the file already exists, we are probably in the recovery mode. Get tag SNPs from the file
		//Note that these are the positions of tag SNPs on a 20Mb chromosome
		else{
			handleTagSNPs = new FileHandler(tagFile);
		    handleTagSNPs.openForRead();
		    handleTagSNPs.nextLine();//header line
		    String currentLine = null;
		    while((currentLine = handleTagSNPs.nextLine()) != null){
		    	tags.add(Integer.parseInt(currentLine));
		    }
   			handleTagSNPs.closeReader();
   		}
		Integer[] tagArray = tags.toArray(new Integer[tags.size()]);
		return tagArray;
	}
	
	/**
	 * @param parentDirName 
	 * @param cauVarID 
	 * @param controlChroms
	 * @param tagSNPs stores the tag snps selected for each causal variant
	 * @throws IOException 
	 */
	private static List<Integer> selectTagSNPs(List<Integer> controlChroms) throws IOException {
		Integer[] snpPos = null;
	    ArrayList<Integer> tagSNPs = new ArrayList<Integer>();
		HaploData hd = new HaploData();
		//snpPos = hd.inputChromMarkerData(parentDirName, (int)(params.markerMAFTresh*100));
		//if(snpPos == null){			
			//Randomly select 1000 controls WITH REPLACEMENT
			List<Subject> subjectList = new ArrayList<Subject>();
			Set<Integer> chromSet = new HashSet<Integer>();
			Collections.shuffle(controlChroms);
			int numOfControlChroms = controlChroms.size();
			for(int i = 1; i <= 1000; i++){		
				Subject unSub  =  new Subject("S" + i, controlChroms.get(Util.randInt(0,numOfControlChroms-1)),
					controlChroms.get(Util.randInt(0,numOfControlChroms-1)), 0);
				subjectList.add(unSub);
				chromSet.add(new Integer(unSub.getChromA()));
				chromSet.add(new Integer(unSub.getChromB()));
			}		
			Integer[] chromArray = chromSet.toArray(new Integer[chromSet.size()]);
			Map<Integer, Double> markerFreq = selectComVariants(chromArray, ParamParser.markerMAFTresh);
			Integer[] markerIndex = markerFreq.keySet().toArray(new Integer[markerFreq.size()]);
			Map<Integer, int[]> chromToMarkerSNPs = Util.gatherMutMarkersfromDistFiles(chromArray, markerIndex);
			snpPos = Util.getSNPRealPosition(markerIndex);
			//The following function combines HaploData.linkageToChrom() and HaploData.prepareMarkerInput()
			hd.setupChromAndMarker(markerFreq, markerIndex, chromToMarkerSNPs, snpPos, subjectList);
			//hd.saveChromMarkerData(parentDirName, (int)(params.markerMAFTresh*100));
	       	markerFreq.clear();
        	markerFreq = null;
        	markerIndex = null;
        	chromToMarkerSNPs.clear();
        	chromToMarkerSNPs = null;
        	subjectList.clear();
        	subjectList = null;
		//}	
			//File dumpDprimeFile = new File(parentDirName + File.separator + "HapTagger_maf" + (int)(params.markerMAFTresh*100) + ".dprime");
			//if(!dumpDprimeFile.exists()){
		        //long startTime = System.currentTimeMillis();
			hd.generateDPrimeTable(snpPos.length);
				//long stopTime = System.currentTimeMillis();
				//System.out.println("Generating the D prime table takes about " + TimeUnit.MILLISECONDS.toMinutes(stopTime - startTime) + " minute(s).\n");
				//hd.saveDprimeToText(dumpDprimeFile, snpPos);
			//}
			//else
				//hd.inputDPrimeTable(snpPos.length,dumpDprimeFile);
			TaggerController taggerExecutor = new TaggerController(hd, HaploData.markers, 1, ParamParser.chromLength-1, ParamParser.num_mark, true);
			taggerExecutor.runTagger();
			if(taggerExecutor.isTaggingCompleted()){
				ArrayList<SNP> tagList = taggerExecutor.getTagSNPs();
				for(SNP aTag : tagList)
					tagSNPs.add(aTag.getPosition());//should do a safe long-to-int conversion
			}
		if(tagSNPs.size() > 1)
			Collections.sort(tagSNPs);
		return tagSNPs;
	}
	
	private static Map<Integer, Double> selectComVariants(Integer[] panelSpecificChroms, double freq_lowerbound) {
		Arrays.sort(panelSpecificChroms);
	    //Mapping a SNP to the frequency of its minor allele; entry format: <index of a SNP regarding the all-SNP array (stored in out.snp-1.gz), its MAF regarding a specific panel> 
	    Map<Integer, Double> snpToFreq = new HashMap<Integer, Double>();
		Map<String, ArrayList<Integer>> chromRecord = Util.getOriginDataLocation(ParamParser.mutInput, panelSpecificChroms);
		for(Entry<String, ArrayList<Integer>> entry: chromRecord.entrySet()){			
			List<Integer> usedChromsInThisFile = entry.getValue();
			int maxRecIndex = Collections.max(usedChromsInThisFile);
			Integer[] recLocation = usedChromsInThisFile.toArray(new Integer[usedChromsInThisFile.size()]);
			Arrays.sort(recLocation);//sort for binary search
			FileHandler getMutationData = new FileHandler(ParamParser.dataOutputDir, entry.getKey());
			getMutationData.openForRead();
			getMutationData.nextLine();//first line is the header
			String nextLine = null;
			Integer mutPos;
			String[] s = null;
			for(int recordCt = 0; (nextLine = getMutationData.nextLine()) != null && recordCt <= maxRecIndex; recordCt++){
				if(recLocation != null && Arrays.binarySearch(recLocation, recordCt) < 0 )
					continue;//this record is not what we want
				s = nextLine.split("\\s+");
				int chromID = Integer.parseInt(s[0]);
				if(Arrays.binarySearch(panelSpecificChroms, chromID) < 0)
					System.out.println("There is something wrong with you strategy in reading distributed files.\n");
				for(int j=1; j<s.length; j++){
					mutPos = Integer.parseInt(s[j]);				
					//update the minorAlleleFreq map
					if(snpToFreq.get(mutPos)!=null){
						double currentValue = snpToFreq.get(mutPos).doubleValue();
						snpToFreq.put(mutPos, currentValue + 1/(double)panelSpecificChroms.length);
					}
					else
						snpToFreq.put(mutPos, 1/(double)panelSpecificChroms.length);
				}
			}
			getMutationData.closeReader();
		}
		
		/*System.out.println("In this data set, MAF ranges from " + Collections.min(minorAlleleFreq.values())
				+ " to " + Collections.max(minorAlleleFreq.values()) + ".\n");*/
		List<Integer> allSNPs = new ArrayList<Integer>(snpToFreq.keySet());
	    Collections.shuffle(allSNPs);
	    
	    Map<Integer, Double> selectedSNPs_maf = new HashMap<Integer, Double>();
		for (int i = 0; i < allSNPs.size(); i++){
	      	int pos = allSNPs.get(i);
	      	double freq = snpToFreq.get(pos);			
			if (freq < freq_lowerbound || (1 - freq) < freq_lowerbound)
				continue;//continue when freq == 1 to avoid monomorphic (i.e., non-polymorphic) sites (GWAS quality control)
			selectedSNPs_maf.put(pos, freq);
		}
		System.out.println("This data set has " + snpToFreq.size() + " common variants with minor allele frequency >= " + ParamParser.markerMAFTresh + ".\n");
		snpToFreq = null;
		return selectedSNPs_maf;
	}
	
	private static void selectCausalVariants(int pos_start, int pos_end, double lo_maf, double hi_maf, int maxNum) {	
		Map<String, ArrayList<Integer>> chromRecord = Util.getOriginDataLocation(ParamParser.mutInput, null);
		int ctReadPosFile = 0;
		for(Entry<String, ArrayList<Integer>> entry: chromRecord.entrySet()){
			ctReadPosFile++;
			if(ctReadPosFile==ParamParser.inputFileNum)
				System.out.println("Last pos file.\n");
			FileHandler getMutationData = new FileHandler(ParamParser.dataOutputDir, entry.getKey());
			getMutationData.openForRead();
			getMutationData.nextLine();//first line is the header			
			String nextLine = null;
			Integer mutInx;
			String[] s = null;
			while((nextLine = getMutationData.nextLine()) != null){
				if(nextLine=="" || nextLine=="\n") continue;
				s = nextLine.split("\\s+");
				int chromID = Integer.parseInt(s[0]);
				for(int j=1; j<s.length; j++){
					mutInx = Integer.parseInt(s[j]);
					boolean exceedFreqThresh = false;
					
					//update the minorAlleleFreq map
					if(minorAlleleFreq.get(mutInx)!=null){
						double currentValue = minorAlleleFreq.get(mutInx).doubleValue() + 1/(double)ParamParser.sampleSize;
						minorAlleleFreq.put(mutInx, currentValue);						
						if(currentValue > hi_maf)
							exceedFreqThresh = true;
					}
					else
						minorAlleleFreq.put(mutInx, 1/(double)ParamParser.sampleSize);
					
					if(!exceedFreqThresh){//update the snpToChrom map only when exceedFreqThresh == false
						Collection<Integer> chromsWithSpecificSNPMut = snpToChrom.get(mutInx);			
						if (chromsWithSpecificSNPMut == null) {
							//initiate the set of chromosomes with mutation on a specific snp, if it does not exist
							//The Set data structure is used to avoid redundancy
							chromsWithSpecificSNPMut = new HashSet<Integer>();//
							snpToChrom.put(mutInx, chromsWithSpecificSNPMut);
						}
						chromsWithSpecificSNPMut.add(chromID);		
					}
					else if(snpToChrom.get(mutInx) != null)
						snpToChrom.remove(mutInx);
				}
				
			}
			getMutationData.closeReader();
		}		
		/*System.out.println("In this data set, MAF ranges from " + Collections.min(minorAlleleFreq.values())
				+ " to " + Collections.max(minorAlleleFreq.values()) + ".\n");*/
		List<Integer> allSNPs = new ArrayList<Integer>(minorAlleleFreq.keySet());
	    Collections.shuffle(allSNPs);
	    //System.out.println("no\n");
		Integer[] snpPos = new Integer[ParamParser.num_mut];
		//Read the absolute positions of all SNPs from the "out.snp-1.gz" file
		FileHandler readSNPPos = new FileHandler(ParamParser.dataOutputDir, ParamParser.posInput);
		readSNPPos.openForRead();
		readSNPPos.nextLine();//header
		String nextLine = null;
		int j = 0;
		while((nextLine = readSNPPos.nextLine()) != null)
			snpPos[j++] = Integer.valueOf(nextLine);
		readSNPPos.closeReader();
	    
	    Map<Integer, Double> selectedSNPs_maf = new HashMap<Integer, Double>();
	    Map<Integer, Collection<Integer>> selectedSNPs_chroms = new HashMap<Integer, Collection<Integer>>();
		for (int i = 0; i < allSNPs.size(); i++){
	      	int snpIndex = allSNPs.get(i);
	      	double freq = minorAlleleFreq.get(snpIndex);
	      	int realPos = snpPos[snpIndex];
			if (realPos < pos_start || realPos > pos_end || freq < lo_maf || freq > hi_maf)
				continue;
			Collection<Integer> c = snpToChrom.get(snpIndex);
			if(c == null)
				continue;
			else{
	      		selectedSNPs_maf.put(realPos, freq);
	      		selectedSNPs_chroms.put(realPos, c);
			}
			if(selectedSNPs_maf.size() >= maxNum) break;
		}
		minorAlleleFreq = selectedSNPs_maf;
		snpToChrom = selectedSNPs_chroms;
	}
	
	


}
