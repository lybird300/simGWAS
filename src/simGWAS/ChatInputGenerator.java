package simGWAS;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;
import java.util.Scanner;

import org.apache.commons.io.FilenameUtils;

import haploview.tagger.Util;

public class ChatInputGenerator /*implements Runnable*/ {
	protected Integer[] tag2MIndex;
	protected Integer[] tag20MIndex;
	protected File chromFile;
	protected int caseIndex;
	protected int ctrlIndex;
	protected String panelFolderName;
	protected int simChromName;

	public ChatInputGenerator(Map<Integer, Integer> tag20MPosTo2MPos, File f, int caseIndex, int ctrlIndex, int chromID) {
		int num_tags = tag20MPosTo2MPos.size();
		if(num_tags!=ParamParser.num_mark)
			System.out.println("Insufficient markers.\n");
		tag2MIndex = tag20MPosTo2MPos.values().toArray(new Integer[num_tags]);
		tag20MIndex = tag20MPosTo2MPos.keySet().toArray(new Integer[num_tags]);
		Arrays.sort(tag2MIndex);
		Arrays.sort(tag20MIndex);
		chromFile = f;
		panelFolderName = FilenameUtils.getFullPathNoEndSeparator(chromFile.getAbsolutePath());
		this.caseIndex = caseIndex;
		this.ctrlIndex = ctrlIndex;
		this.simChromName = chromID;
	}

	public ChatInputGenerator(Integer[] tagPos, String folderName, int caseIndex, int ctrlIndex, int chromID, int firstRep) {
		int num_tags = tagPos.length;
		if(num_tags!=ParamParser.num_mark)
			System.out.println("Insufficient markers.\n");
		tag2MIndex = null;
		tag20MIndex = tagPos;
		Arrays.sort(tag20MIndex);
		panelFolderName = folderName;
		this.caseIndex = caseIndex;
		this.ctrlIndex = ctrlIndex;
		this.simChromName = chromID;
	}

	public void run() {
		try {
			if(tag2MIndex==null){
				checkCHATChromSpecFile();
			}
			else{
				generateCHATChromSpecFile(ParamParser.num_ca[caseIndex] + ParamParser.num_co[ctrlIndex]);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}

	private void checkCHATChromSpecFile() {			
		File dir_4th = new File(panelFolderName, "ChatInput_tagSNP");
		if(!dir_4th.exists()){
    		System.out.println(dir_4th.getAbsolutePath() + " does not exist.\n");
			return;
		}
		String outputFileName = null;
		for(int replicateID = 0; replicateID < ParamParser.replica; replicateID++){
			outputFileName = "GenotypeChrom_" + ParamParser.padNumberAsString(simChromName,2) + "_Rep" + replicateID + ".gz";
			File chatInputFile = new File(dir_4th.getAbsoluteFile(), outputFileName);
			if(!chatInputFile.exists()){
	    		System.out.println(chatInputFile.getAbsolutePath() + " does not exist.\n");
				return;
			}
			else
				sanityCheck(tag20MIndex, chatInputFile, ParamParser.num_ca[caseIndex] + ParamParser.num_co[ctrlIndex], true);
		}
	}

	private void generateCHATChromSpecFile(int numOfSubjects) throws IOException {
		FileHandler readPanelSpecificChroms = new FileHandler(chromFile);
		readPanelSpecificChroms.openForRead();
		readPanelSpecificChroms.nextLine();//first line is the header
		
		File dir_4th = new File(panelFolderName, "ChatInput_tagSNP");
		if(!dir_4th.exists()) dir_4th.mkdirs();
		
		String nextLine;
		String[] s = null;
		while((nextLine = readPanelSpecificChroms.nextLine()) != null){
			s = nextLine.split(" ");								
			int replicateID = Integer.valueOf(s[0]);
			File dir_5th = new File(dir_4th.getAbsolutePath() + File.separator + "Rep" + replicateID + File.separator + "ChromosomeSpecificData");
			if(!dir_5th.exists()) dir_5th.mkdirs();
			File chatInputFile = new File(dir_5th, "GenotypeChrom_" + ParamParser.padNumberAsString(simChromName,2) + ".gz");
			if(chatInputFile.exists() && sanityCheck(tag20MIndex, chatInputFile, numOfSubjects, true)) continue;
			else{
				//populate the array of selected chromIDs
				Integer[] panelSpecificChroms = new Integer[s.length-1];
				for(int j = 1; j < s.length; j++)
					panelSpecificChroms[j-1] = Integer.valueOf(s[j]);
				
				//If tag SNPs are not provided, randomly select certain number (indicated by params.num_mark)
				//of common variants (MAF >= 5% indicated by params.maf_mark) as markers
				/*if(markerIndex==null){		
					Map<Integer, Double> markerMaf = selectComVariants(panelSpecificChroms, params.markerMAFTresh);
					markerIndex = markerMaf.keySet().toArray(new Integer[params.num_mark]);		
				}*/
				
				//The function below returns a map of <selected chromID, mutated markers>
				//Recording only mutated markers (labeled as "1" while ancestral markers as "2") rather than all of them
				//allows us to reconstruct raw sequence data relatively easily, as the former is much fewer than the latter.
				Map<Integer, int[]> chromToMarkers = Util.gatherMutMarkersfromDistFiles(panelSpecificChroms, tag2MIndex);			
				outputCHATInputFile(chatInputFile, replicateID, chromToMarkers, numOfSubjects);
			}
		}
		s = null;
		readPanelSpecificChroms.closeReader();
	}
	
	/**
	 * Simulate input files in the "ChromosomeSpecificData" folder. File format refers to
	 * (1)chat-prep/org.renci.chat.interpreter.MakeChromSpecificChatFromLgenRunnable.run()
	 * (2)chat-common/org.renci.chat.common.DataSet.readLgenDataSet(...)
	 * (3)chat-common/org.renci.chat.common.DataSet.writeDataSet(File)
	 * chromToHaplotype stores all selected haplotypes (chromosome segments based on markers).
	 * The original position of each SNP in the haplotype is stored in the Integer array of "markers"
	 * @throws IOException 
	 */
	protected void outputCHATInputFile(File chatInputFile, int replicateID, Map<Integer, int[]> chromToMarkerSNPs, int numOfSubjects){	    
		if(chatInputFile.exists()){
			//System.out.println("Check " + chatSetFile.getAbsolutePath() + "\n");
			//return;
			chatInputFile.delete();
		}
		FileHandler readSubjectFile = new FileHandler(panelFolderName, "CaseControl_Rep" + replicateID + ".csv");
		readSubjectFile.openForRead();
		readSubjectFile.nextLine();//first line is the header	

		String subjectLine = "";
		String dxLine = "";
		
		int numOfMarkers = tag20MIndex.length;
		int[][] genotypes = new int[numOfMarkers][numOfSubjects]; 
		 
		int subjectCt = 0;
		String nextLine = null;
		while((nextLine = readSubjectFile.nextLine()) != null && subjectCt < numOfSubjects){
			Scanner s = new Scanner(nextLine);
			s.useDelimiter(",");
			if(subjectCt > 0){
				subjectLine += ",";
				dxLine += ",";
			}
			subjectLine += s.next();//get subject ID
			int[] markerDataA = chromToMarkerSNPs.get(Integer.parseInt(s.next()));
			int[] markerDataB = chromToMarkerSNPs.get(Integer.parseInt(s.next()));
			
			for(int i = 0; i < numOfMarkers; i++)					
				genotypes[i][subjectCt] = markerDataA[i] - 1 + markerDataB[i] - 1;//Alleles are represented by 1(minor) and 2 (major) in Dan's data set; here they are converted to 0 and 1 respectively.			
			dxLine += s.next().equals("0")? "1":"2";//unaffected (0->1); affect (1->2)
			
			subjectCt++;
			s.close();
		}
		if(subjectCt < numOfSubjects)
			System.out.println("The subject file " + readSubjectFile.getAbsolutePath() + " has been corrupted.\n");

		String sexLine = "";
		String includeLine = "";
		for(int i=0; i<subjectCt; i++){
			if(i>0){
				sexLine += ',';
				includeLine += ',';
			}
			sexLine += '2'; //2-Female, 1-Male, 0-unknown
			includeLine += '1';
		}
	
		FileHandler writeCHATinput = new FileHandler(chatInputFile);
		writeCHATinput.openForWrite(true);
		writeCHATinput.writeString(subjectLine);
		writeCHATinput.writeString(dxLine);
		writeCHATinput.writeString(sexLine);
		writeCHATinput.writeString(includeLine);
		for(int j = 0; j < numOfMarkers; j++){
			String aMarker = simChromName + "," +  tag20MIndex[j] + ",rs" + j + ",-9999.00000000,1,";
			//String aMarker = simChrom + "," +  markerIndex[j] + ",rs" + j + ",-9999.00000000,1,";//for testing
			for(int k=0; k < numOfSubjects; k++)
				aMarker += String.valueOf(genotypes[j][k]);
			writeCHATinput.writeString(aMarker);	
		}
		writeCHATinput.preventFutureWriting();
		writeCHATinput.closeWriter();
	}
	
	protected boolean sanityCheck(Integer[] markerPos, File chatInputFile, int numOfSubjects, boolean detailedCheck) {
		//double fileSizeThreshInKB = (numOfSubjects==2000)?11000:35000;
		//if(chatInputFile.length()/1024.0 < fileSizeThreshInKB) return false;
		/*else*/ if(detailedCheck){
			FileHandler readChatInputFile = new FileHandler(chatInputFile);
			readChatInputFile.openForRead();
			String nextLine = null;
			int rowCt = 1;
			for(;rowCt <= 4; rowCt++){
				String[] s = lineSplit(readChatInputFile.nextLine(), ',');
				if(s.length != numOfSubjects){
					System.out.println(chatInputFile.getAbsolutePath() + "\nRow " + rowCt + "does not have enough data.\n");
					return false;
				}
			}
			int markerIndex = 0;
			while((nextLine = readChatInputFile.nextLine())!=null){
				String[] ff = lineSplit(nextLine, ',');
			    if(ff.length != 6) return false;
			    if(ff[1].equalsIgnoreCase("null")){
			    	System.out.println(chatInputFile.getAbsolutePath() + "One marker position appears to be null.\n");
			    	return false;
			    }
			    if(Integer.parseInt(ff[1]) != markerPos[markerIndex++]){
			    	System.out.println(chatInputFile.getAbsolutePath() + "One marker position appears to be wrong.\n");
			    	return false;
			    }
			    else{
			    	char[] genotype = ff[5].trim().toCharArray();
			    	if(genotype.length != numOfSubjects){
			    		System.out.println(chatInputFile.getAbsolutePath() + "\n" + ff[2] + "does not have enough data.\n");
						return false;
			    	}
			    }
			    rowCt++;
			}
			if(rowCt < ParamParser.num_mark + 4){
				System.out.println(chatInputFile.getAbsolutePath() + "\nOnly have" + (rowCt-4) + "marker records.\n");
				return false;
			}
		}
		return true;
	}
	public String[] lineSplit(final String input, char separator) {
		int pieces = 0;

		// First we count how many pieces we will need to store
		// ( = separators + 1)
		int position = -1;
		do {
			pieces++;
			position = input.indexOf(separator, position + 1);
		} while (position != -1);

		// Then we allocate memory
		final String[] result = new String[pieces];

		// And start cutting and copying the pieces.
		int previousposition = 0;
		int currentposition = input.indexOf(separator);
		int piece = 0;
		final int lastpiece = pieces - 1;
		while (piece < lastpiece) {
			result[piece++] = input.substring(previousposition, currentposition);
			previousposition = currentposition + 1;
			currentposition = input.indexOf(separator, previousposition);
		}
		result[piece] = input.substring(previousposition);

		return result;
	}


}
