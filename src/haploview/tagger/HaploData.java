package haploview.tagger;

import java.io.*;
import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;

import simGWAS.FileHandler;
import simGWAS.Subject;

public class HaploData{

    private byte[][] chromosomes;
    //private Haplotype[][] haplotypes;
    public static DPrimeTable dpTable;
    public static List<SNP> markers = new ArrayList<SNP>();
    
	public Integer[] inputChromMarkerData(String parentDirName, int mafThresh) throws FileNotFoundException {
		File chromFile = new File(parentDirName + File.separator + "HapTagger_maf" + mafThresh + ".chrom");
		File markerFile = new File(parentDirName + File.separator + "HapTagger_maf" + mafThresh + ".marker");
		if(!chromFile.exists() || !markerFile.exists()) return null;
		List<byte[]> ret = new ArrayList<byte[]>();
		Scanner inputChroms = new Scanner(chromFile);
		while (inputChroms.hasNextLine()) {
			 // read a line, and turn it into the characters
			 String[] oneLine = inputChroms.nextLine().split(" ");
			 byte[] intLine = new byte[oneLine.length];
			 for(int i=0; i < intLine.length; i++)
				 intLine[i] = Byte.parseByte(oneLine[i].trim());
			 ret.add(intLine);
		}
		inputChroms.close();
		chromosomes = new byte[ret.size()][];
		for(int i=0; i < ret.size(); i++){
			byte[] row = ret.get(i);
			chromosomes[i] = row;
		}

		FileHandler inputMarkers = new FileHandler(markerFile);
		inputMarkers.openForRead();
		inputMarkers.nextLine();//header
		String markerLine;
		List<Integer> positions = new ArrayList<Integer>();
	    while((markerLine = inputMarkers.nextLine())!=null){
			Scanner s = new Scanner(markerLine).useDelimiter(" ");
			int realPos = Integer.parseInt(s.next());
			positions.add(realPos);
			double mutfreq = Double.parseDouble(s.next());
			if(mutfreq <= 0.5)
				markers.add(new SNP(realPos, mutfreq,(byte)2,(byte)3));
			else
				markers.add(new SNP(realPos, 1-mutfreq,(byte)3,(byte)2));
			s.close();
	    }
	    inputMarkers.closeReader();
	    return positions.toArray(new Integer[positions.size()]);
	}

	public void saveChromMarkerData(String parentDirName, int mafThresh) throws IOException {
		File chromFile = new File(parentDirName + File.separator + "HapTagger_maf" + mafThresh + ".chrom");
		FileHandler outputChroms = new FileHandler(chromFile);
		outputChroms.openForWrite(true);
		for(int i = 0; i < chromosomes.length; i++){
			for(int j = 0; j < chromosomes[i].length; j++){
				outputChroms.writer.write(chromosomes[i][j] + " ");
			}
			outputChroms.writer.write("\n");
		}
		outputChroms.closeWriter();
		File markerFile = new File(parentDirName + File.separator + "HapTagger_maf" + mafThresh + ".marker");
		FileHandler outputMarkers = new FileHandler(markerFile);
		outputMarkers.openForWrite(true);
		outputMarkers.writeString("realPos\tmaf");
        for (int i = 0; i < markers.size(); i++){
        	SNP marker = markers.get(i);
        	int realPos = marker.getPosition();
			double maf = marker.getMAF();
			outputMarkers.writeString(realPos + " " + maf);
        }
        outputMarkers.closeWriter();
	}
    
	/*private class ParallelChromAndMarkerSetter implements Runnable{
    	private int from; //inclusive
    	private int to; //exclusive
    	private Map<Integer, Double> markerMaf = null;
		private Map<Integer, int[]> chromToMarkerSNPs = null;
		private Integer[] snpPos = null;
		private List<Subject> subjectList = null;
		private Integer[] markerIndex = null;
    	
    	private ParallelChromAndMarkerSetter(int fromIndex, int toIndex, Map<Integer, Double> markerMaf,
			Map<Integer, int[]> chromToMarkerSNPs, Integer[] snpPos2, List<Subject> subjectList, Integer markerIndex[]) {
    		super();
            this.from = fromIndex;
            this.to = toIndex;
            this.markerMaf = markerMaf;
            this.chromToMarkerSNPs = chromToMarkerSNPs;
            this.snpPos = snpPos2;
            this.subjectList = subjectList;
            this.markerIndex = markerIndex;
        }

		@Override
		public void run() {
			System.out.println(Thread.currentThread() + " set up markers from Sample " + this.from  + " (inclusive) to Sample " + this.to + " (exclusive)\n");
			int numMarkers = snpPos.length;
			for(int k = this.from*2; k < this.to*2; k+=2){
				Subject s = this.subjectList.get(k/2);
				int[] haplotypeA = this.chromToMarkerSNPs.get(s.getChromA());
				int[] haplotypeB = this.chromToMarkerSNPs.get(s.getChromB());
				for(int i = 0; i < numMarkers; i++){			
					byte thisMarkerA = (byte)(haplotypeA[i]+2);
	                byte thisMarkerB = (byte)(haplotypeB[i]+2);
	                if (thisMarkerA == thisMarkerB){
	                	chromosomes[k][i] = thisMarkerA;
	                	chromosomes[k+1][i] = thisMarkerB;
		            }else{
		            	chromosomes[k][i] = (byte)(4+thisMarkerA);
		            	chromosomes[k+1][i] = (byte)(4+thisMarkerB);
		            }			
					int realPos = this.snpPos[i];
					int relaPos = this.markerIndex[i];
					double maf = this.markerMaf.get(relaPos);
					SNP marker = new SNP(String.valueOf(realPos), realPos, maf);
					markers.add(marker);
				}
			}		
		}
	}
	
	public void setupChromAndMarker(Map<Integer, Double> markerMaf, Integer[] markerIndex,
			Map<Integer, int[]> chromToMarkerSNPs, Integer[] snpPos, List<Subject> subjectList){
		int numMarkers = snpPos.length;
		int numChroms = subjectList.size()*2;
		chromosomes = new byte[numChroms][numMarkers];	
		int numOfThreads = Math.min(Runtime.getRuntime().availableProcessors(),32);
        ExecutorService service = Executors.newFixedThreadPool(numOfThreads);
        int toIndex = -1;
	   	for (int i = 0; i < numOfThreads; i++) {
    		int itemsInCurThread = getThreadItemCount(numChroms, i, numOfThreads);
    		int fromIndex = Math.max(itemsInCurThread * i, toIndex);
    		toIndex = fromIndex + itemsInCurThread;
    		ParallelChromAndMarkerSetter setter = new ParallelChromAndMarkerSetter(fromIndex, toIndex, 
    				markerMaf, chromToMarkerSNPs, snpPos, subjectList, markerIndex);  		
    		service.execute(setter);
    	}
         try {
        	service.shutdown();
 			service.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
 		} catch (InterruptedException ignore) {			
 		}
         finally {
        	 markerMaf.clear();
        	 markerMaf = null;
        	 markerIndex = null;
        	 chromToMarkerSNPs.clear();
        	 chromToMarkerSNPs = null;
        	 snpPos = null;
        	 subjectList.clear();
        	 subjectList = null;
        	    if (!service.isTerminated()) {
        	        System.err.println("cancel non-finished tasks");
        	    }
        	    service.shutdownNow();
         }
    }*/
	
	public void setupChromAndMarker(Map<Integer, Double> markerMaf, Integer[] markerIndex,
			Map<Integer, int[]> chromToMarkerSNPs, Integer[] snpPos, List<Subject> subjectList){
		int numMarkers = snpPos.length;
		int numChroms = subjectList.size()*2;
		chromosomes = new byte[numChroms][numMarkers];
		for(int k = 0; k < numChroms; k+=2){
			Subject s = subjectList.get(k/2);
			int[] haplotypeA = chromToMarkerSNPs.get(s.getChromA());
			int[] haplotypeB = chromToMarkerSNPs.get(s.getChromB());
			for(int i = 0; i < numMarkers; i++){			
				byte thisMarkerA = (byte)(haplotypeA[i]+1);
                byte thisMarkerB = (byte)(haplotypeB[i]+1);
                if (thisMarkerA == thisMarkerB){
                	chromosomes[k][i] = thisMarkerA;
                	chromosomes[k+1][i] = thisMarkerB;
	            }else{
	            	chromosomes[k][i] = (byte)(4+thisMarkerA);
	            	chromosomes[k+1][i] = (byte)(4+thisMarkerB);
	            }			
			}
		}
		markers = new ArrayList<SNP>();
		for(int i = 0; i < numMarkers; i++){
			int realPos = snpPos[i];
			double mutfreq = markerMaf.get(markerIndex[i]);
			if(mutfreq <= 0.5)
				markers.add(new SNP(realPos, mutfreq,(byte)2,(byte)3));
			else
				markers.add(new SNP(realPos, 1-mutfreq,(byte)3,(byte)2));
		}
	}
    
    private class ParallelDPrimeTableGenerator implements Runnable{
    	final private int from; //inclusive
    	final private int to; //exclusive
    	final private DPrimeTable subTable;
    	final private int chromSize;
    	
    	private ParallelDPrimeTableGenerator(int fromIndex, int toIndex, DPrimeTable dpTable, int chromSize) {
    		super();
            this.from = fromIndex;
            this.to = toIndex;
            this.subTable = dpTable;
            this.chromSize = chromSize;
        }

		public void run() {
			System.out.println(Thread.currentThread() + " calculate D prime on: " + this.from  + "(inclusive) to " + this.to + "(exclusive)\n");
    		List<PairwiseLinkage> dpTemp = null;
    		for (int pos1 = from; pos1 < to; pos1++) {
    			dpTemp = new ArrayList<PairwiseLinkage>();
                for (int pos2 = pos1 + 1; pos2 < chromSize; pos2++){
                    //if the markers are too far apart don't try to compare them
                    long sep = markers.get(pos2).getPosition() - markers.get(pos1).getPosition();
                    if (sep <= Tagger.DEFAULT_MAXDIST){
                    	PairwiseLinkage p = computeDPrime(pos1,pos2);
                        if(p!=null){
                        	dpTemp.add(p);
                        }
                        //else
                        	//System.out.println("null pairwise linkage result.\n");
                    }
                    else
                    	break;              
                }
                subTable.addMarker(dpTemp,pos1);
                dpTemp = null;
		   }
    		dpTemp = null;
		}
    }
    
    private int getThreadItemCount(int numOfItems, int threadID, int numOfThreads){
    	return (int) (numOfItems / numOfThreads) + (((numOfItems % numOfThreads) >= (threadID + 1)) ? 1 : 0);
    }

    public void generateDPrimeTable(int chromSize){
        //calculating D prime requires the number of each possible 2 marker haplotype in the dataset
        dpTable = new DPrimeTable(chromSize);
        int numOfThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService service = Executors.newFixedThreadPool(numOfThreads);
    	//ArrayList<Thread> threadList = new ArrayList<Thread>(numOfThreads);
        int toIndex = -1;
    	for (int i = 0; i < numOfThreads; i++) {
    		int itemsInCurThread = getThreadItemCount(chromSize - 1, i, numOfThreads);
    		int fromIndex = Math.max(itemsInCurThread * i, toIndex);
    		toIndex = fromIndex + itemsInCurThread;
    		ParallelDPrimeTableGenerator dpGenerator = new ParallelDPrimeTableGenerator(fromIndex, toIndex, dpTable, chromSize);  		
    		service.execute(dpGenerator);
    	}
         try {
        	service.shutdown();
 			service.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
 		} catch (InterruptedException ignore) {			
 		}
         finally {
        	    if (!service.isTerminated()) {
        	        System.err.println("cancel non-finished tasks");
        	    }
        	    service.shutdownNow();
         }
    }
    
    /*public void generateDPrimeTable(int chromSize){
    	//calculating D prime requires the number of each possible 2 marker haplotype in the dataset
        dpTable = new DPrimeTable(chromSize);
		List<PairwiseLinkage> dpTemp = null;
		for (int pos1 = 0; pos1 < chromSize - 1; pos1++) {
			dpTemp = new ArrayList<PairwiseLinkage>();
            for (int pos2 = pos1 + 1; pos2 < chromSize; pos2++){
                //if the markers are too far apart don't try to compare them
                long sep = markers.get(pos2).getPosition() - markers.get(pos1).getPosition();
                if (sep <= Tagger.DEFAULT_MAXDIST){
                	PairwiseLinkage p = computeDPrime(pos1,pos2);
                    if(p!=null){
                    	dpTemp.add(p);
                    }
                    //else
                    	//System.out.println("null pairwise linkage result.\n");
                }
                else
                	break;              
            }
            dpTable.addMarker(dpTemp,pos1);
            dpTemp = null;
	   }
		dpTemp = null;
    }*/
    
	public void inputDPrimeTable(int chromSize, File dumpDprimeFile) {
		dpTable = new DPrimeTable(chromSize);
		FileHandler handleTaggerDpTable = new FileHandler(dumpDprimeFile);
		handleTaggerDpTable.openForRead();
		handleTaggerDpTable.nextLine();//header
		String curMarkerPosi = null;
		List<PairwiseLinkage> dpTemp = null;
        String dpTableLine = null;
        int pos = 0;
        while((dpTableLine = handleTaggerDpTable.nextLine())!=null){
			Scanner dpEntry = new Scanner(dpTableLine);
			dpEntry.useDelimiter("\t");
			String markerPosi = dpEntry.next();
			if(curMarkerPosi == null ||!markerPosi.equalsIgnoreCase(curMarkerPosi)){
				curMarkerPosi = markerPosi;
				if(dpTemp != null){
					dpTable.addMarker(dpTemp, pos);
					pos++;
				}
				dpTemp= new ArrayList<PairwiseLinkage>();
			}
			String markerPosj = dpEntry.next();
			if(!markerPosj.equalsIgnoreCase("null")){
				double dprime = Double.parseDouble(dpEntry.next());
				double lod = Double.parseDouble(dpEntry.next());
				double r2 = Double.parseDouble(dpEntry.next());
				//double ci_low = Double.parseDouble(dpEntry.next());
				//double ci_high = Double.parseDouble(dpEntry.next());
				PairwiseLinkage pl = new PairwiseLinkage(dprime, lod, r2/*, ci_low, ci_high, null*/);
				dpTemp.add(pl);
			}
			else{
				dpTemp = null;
				dpTable.addMarker(null, pos);
			}
			dpEntry.close();
        }
		handleTaggerDpTable.closeReader();
	}
	
    /*public Haplotype[][] generateHaplotypes(List<Integer[]> victor, boolean storeEM){
        Haplotype[][] rawHaplotypes = new Haplotype[victor.size()][];

        for (int k = 0; k < victor.size(); k++){
            Integer[] preFiltBlock = victor.get(k);
            Integer[] theBlock;

            Integer[] selectedMarkers = new Integer[0];
            int[] equivClass = new int[0];
            if (preFiltBlock.length > 30){
                equivClass = new int[preFiltBlock.length];
                int classCounter = 0;
                for (int x = 0; x < preFiltBlock.length; x++){
                    int marker1 = preFiltBlock[x];

                    //already been lumped into an equivalency class
                    if (equivClass[x] != 0){
                        continue;
                    }

                    //start a new equivalency class for this SNP
                    classCounter ++;
                    equivClass[x] = classCounter;

                    for (int y = x+1; y < preFiltBlock.length; y++){
                        int marker2 = preFiltBlock[y];
                        if (marker1 > marker2){
                            int tmp = marker1; marker1 = marker2; marker2 = tmp;
                        }
                        if ( dpTable.getLDStats(marker1,marker2) != null
                                && dpTable.getLDStats(marker1,marker2).getRSquared() == 1.0){
                            //these two SNPs are redundant
                            equivClass[y] = classCounter;
                        }
                    }
                }

                //parse equivalency classes
                selectedMarkers = new Integer[classCounter];
                for (int x = 0; x < selectedMarkers.length; x++){
                    selectedMarkers[x] = -1;
                }
                for (int x = 0; x < classCounter; x++){
                    for (int y = 0; y < equivClass.length; y++){
                        if (equivClass[y] == x+1){
                            selectedMarkers[x] = preFiltBlock[y];
                        }
                    }
                }

                theBlock = selectedMarkers;
            }else{
                theBlock = preFiltBlock;
            }

            //kirby patch
            EM theEM = new EM(chromosomes);
            theEM.doEM(theBlock);

            Haplotype[] tempArray = new Haplotype[theEM.numHaplos()];
            int[][] returnedHaplos = theEM.getHaplotypes();
            double[] returnedFreqs = theEM.getFrequencies();
            for (int i = 0; i < theEM.numHaplos(); i++){
                int[] genos = new int[returnedHaplos[i].length];
                for (int j = 0; j < genos.length; j++){
                    if (returnedHaplos[i][j] == 1){
                        markers.get(theBlock[j]);
						genos[j] = SNP.getMajor();
                    }else{
                        markers.get(theBlock[j]);
						if (SNP.getMinor() == 0){
                            genos[j] = 8;
                        }else{
                            markers.get(theBlock[j]);
							genos[j] = SNP.getMinor();
                        }
                    }
                }

                if (selectedMarkers.length > 0){
                    //we need to reassemble the haplotypes
                    Hashtable<Integer, Integer> hapsHash = new Hashtable<Integer, Integer>();
                    //add to hash all the genotypes we phased
                    for (int q = 0; q < genos.length; q++){
                        hapsHash.put(new Integer(theBlock[q]), new Integer(genos[q]));
                    }
                    //now add all the genotypes we didn't bother phasing, based on
                    //which marker they are identical to
                    for (int q = 0; q < equivClass.length; q++){
                        int currentClass = equivClass[q]-1;
                        if (selectedMarkers[currentClass] == preFiltBlock[q]){
                            //we alredy added the phased genotypes above
                            continue;
                        }
                        int indexIntoBlock=0;
                        for (int x = 0; x < theBlock.length; x++){
                            if (theBlock[x] == selectedMarkers[currentClass]){
                                indexIntoBlock = x;
                                break;
                            }
                        }
                        //this (somewhat laboriously) reconstructs whether to add the minor or major allele
                        //for markers with MAF close to 0.50 we can't use major/minor alleles to match
                        //'em up 'cause these might change given missing data
                        boolean success = false;
                        if (markers.get(selectedMarkers[currentClass]).getMAF() > 0.4){
                            for (int z = 0; z < chromosomes.length; z++){
                                byte[] thisChrom = chromosomes[z];
                                byte[] nextChrom = chromosomes[++z];
                                int theGeno = (int)thisChrom[selectedMarkers[currentClass]];
                                int nextGeno = (int)nextChrom[selectedMarkers[currentClass]];
                                if (theGeno == nextGeno && theGeno == genos[indexIntoBlock]
                                        && (int)thisChrom[preFiltBlock[q]] != 0){
                                    hapsHash.put(new Integer(preFiltBlock[q]),
                                            new Integer(thisChrom[preFiltBlock[q]]));
                                    success = true;
                                    break;
                                }
                            }
                        }

                        //either we didn't use careful counting or it didn't work due to missing data
                        if(!success){
                            markers.get(selectedMarkers[currentClass]);
							if (SNP.getMajor() ==
                                    genos[indexIntoBlock]){
                                hapsHash.put(new Integer(preFiltBlock[q]),
                                        new Integer(SNP.getMajor()));
                            }else{
                                hapsHash.put(new Integer(preFiltBlock[q]),
                                        new Integer(SNP.getMinor()));
                            }
                        }
                    }
                    genos = new int[preFiltBlock.length];
                    for (int q = 0; q < preFiltBlock.length; q++){
                        genos[q] = ((Integer)hapsHash.get(new Integer(preFiltBlock[q]))).intValue();
                    }
                }

                if(storeEM) {
                    tempArray[i] = new Haplotype(genos, returnedFreqs[i], preFiltBlock, theEM);
                }else {
                    tempArray[i] = new Haplotype(genos, returnedFreqs[i], preFiltBlock, null);
                }
            }
            //make the rawHaplotypes array only large enough to hold haps
            //which pass threshold above
            rawHaplotypes[k] = new Haplotype[theEM.numHaplos()];
            for (int z = 0; z < theEM.numHaplos(); z++){
                rawHaplotypes[k][z] = tempArray[z];
            }
        }

        return rawHaplotypes;
    }*/


    //this method computes the dPrime value for the pair of markers pos1,pos2.
    //the method assumes that the given pair are not too far apart (ie the distance
    //between them is less than maximum distance).
    public PairwiseLinkage computeDPrime(int pos1, int pos2){
        
    	int AA = 0;
        int AB = 1;
        int BA = 2;
        int BB = 3;
        double TOLERANCE = 0.00000001;
        double LN10 = Math.log(10.0);
        int unknownDH=-1;
        int total_chroms=-1;
        double const_prob=-1.0;
        double[] known = new double[5];
        double[] numHaps = new double[4];
        double[] probHaps = new double[4];
    	
        int doublehet = 0;
        int[][] twoMarkerHaplos = new int[3][3];

        for (int i = 0; i < twoMarkerHaplos.length; i++){
            for (int j = 0; j < twoMarkerHaplos[i].length; j++){
                twoMarkerHaplos[i][j] = 0;
            }
        }
        //The created matrix is used to calculate haplotype frequency and it looks like below, 
        //where A and B represents major and minor alleles respectively
        //*********pos2****
        //			1	2
        //pos1	1   AA	AB
        //		2	BA	BB
        //*****************

        //check for non-polymorphic (i.e., monomorphic) markers
        if (markers.get(pos1).getMAF() == 0 || markers.get(pos2).getMAF() == 0){
            return null;
        }

        int[] marker1num = new int[5]; int[] marker2num = new int[5];

        marker1num[0]=0;
		marker1num[markers.get(pos1).getMajor()]=1;
		marker1num[markers.get(pos1).getMinor()]=2;
        marker2num[0]=0;
		marker2num[markers.get(pos2).getMajor()]=1;
		marker2num[markers.get(pos2).getMinor()]=2;

        byte a1,a2,b1,b2;
        //iterate through all chromosomes in dataset
        for (int i = 0; i < chromosomes.length; i++){
            //assign alleles for each of a pair of chromosomes at a marker to four variables
                a1 = chromosomes[i][pos1];
                a2 = chromosomes[i][pos2];
                b1 = chromosomes[++i][pos1];
                b2 = chromosomes[i][pos2];

                if (a1 == 0 || a2 == 0 || b1 == 0 || b2 == 0){
                    System.out.println("Missing data...\n");
                } else if (((a1 >= 5 || b1 >= 5) && (a2 >= 5 || b2 >= 5)) || (a1 >= 5 && !(a2 == b2)) || (a2 >= 5 && !(a1 == b1))){
                    doublehet++;//both markers are heterozygotes regarding this specific pair of chromosomes
                    //find doublehets and resolved haplotypes
                } else if (a1 >= 5 || b1 >= 5){//only pos1 is heterozygote
                    twoMarkerHaplos[1][marker2num[a2]]++;
                    twoMarkerHaplos[2][marker2num[a2]]++;
                } else if (a2 >= 5 || b2 >= 5){//only pos2 is heterozygote
                    twoMarkerHaplos[marker1num[a1]][1]++;
                    twoMarkerHaplos[marker1num[a1]][2]++;
                } else {//both markers are homozygotes
                    twoMarkerHaplos[marker1num[a1]][marker2num[a2]]++;
                    twoMarkerHaplos[marker1num[b1]][marker2num[b2]]++;
                }
        }
        //another monomorphic marker check
        int r1, r2, c1, c2;
        r1 = twoMarkerHaplos[1][1] + twoMarkerHaplos[1][2];
        r2 = twoMarkerHaplos[2][1] + twoMarkerHaplos[2][2];
        c1 = twoMarkerHaplos[1][1] + twoMarkerHaplos[2][1];
        c2 = twoMarkerHaplos[1][2] + twoMarkerHaplos[2][2];
        if ( (r1==0 || r2==0 || c1==0 || c2==0) && doublehet == 0){
            return  new PairwiseLinkage(1,0,0/*,0,0,new double[0]*/);
        }

        //compute D Prime for this pair of markers.
        int /*i,*/count;
        //int j,k,itmp;
        //int low_i = 0;
        //int high_i = 0;
        double loglike, oldloglike;// meand, mean2d, sd;
        //double tmp;//g,h,m,tmp,r;
        double num, denom1, denom2, denom, dprime;//, real_dprime;
        double pA1, pB1, pA2, pB2, loglike1, loglike0, rsq;
        //double tmpAA, tmpAB, tmpBA, tmpBB, dpr;// tmp2AA, tmp2AB, tmp2BA, tmp2BB;
        //double total_prob, sum_prob;
        //double lsurface[] = new double[101];

        //store arguments in externals and compute allele frequencies
        known[AA]=twoMarkerHaplos[1][1];//numbers of haplotype AA that we are sure about
        known[AB]=twoMarkerHaplos[1][2];//numbers of haplotype AB that we are sure about
        known[BA]=twoMarkerHaplos[2][1];//numbers of haplotype BA that we are sure about
        known[BB]=twoMarkerHaplos[2][2];//numbers of haplotype BB that we are sure about
        //store unphasing counts: we know both positions (1 and 2) have major allele (A) and minor allel (B) in the sample,
        //but not sure whether it is A1B2/B1A2 or B1A2/A1B2. Only the double heterozygote [AB][AB] results in ambiguous reconstruction,
        //so we'll count the obligates then tack on the [AB][AB] for clarity.
        unknownDH=doublehet;        
        total_chroms= (int)(known[AA]+known[AB]+known[BA]+known[BB]+(2*unknownDH));
        pA1 = (known[AA]+known[AB]+unknownDH) / (double) total_chroms;//estimated major allele frequency of SNP at pos1 (inflate major allele frequency)
        pB1 = 1.0-pA1;//estimated minor allele frequency of SNP at pos1
        pA2 = (known[AA]+known[BA]+unknownDH) / (double) total_chroms;//estimated major allele frequency of SNP at pos2
        pB2 = 1.0-pA2;//estimated minor allele frequency of SNP at pos2
        //const_prob = 0.1;

        //set initial conditions
        if (const_prob < 0.00) {
            probHaps[AA]=pA1*pA2;//expected frequency of haplotype AA if the two SNPs are linkage equilibrium
            probHaps[AB]=pA1*pB2;//expected frequency of haplotype AB if the two SNPs are linkage equilibrium
            probHaps[BA]=pB1*pA2;//expected frequency of haplotype BA if the two SNPs are linkage equilibrium
            probHaps[BB]=pB1*pB2;//expected frequency of haplotype BB if the two SNPs are linkage equilibrium
        } else {
            probHaps[AA]=const_prob;
            probHaps[AB]=const_prob;
            probHaps[BA]=const_prob;
            probHaps[BB]=const_prob;;

            /*so that the first count step will produce an initial estimate without inferences (this should
            be closer and therefore speedier than assuming they are all at equal frequency) */
            //former count_haps(0)
            numHaps[AA] = known[AA];
            numHaps[AB] = known[AB];
            numHaps[BA] = known[BA];
            numHaps[BB] = known[BB];
            
            //estimate_p();
            double total= numHaps[AA]+numHaps[AB]+numHaps[BA]+numHaps[BB]+(4.0*const_prob);
            probHaps[AA]=(numHaps[AA]+const_prob)/total; if (probHaps[AA] < 1e-100) probHaps[AA]=1e-100;
            probHaps[AB]=(numHaps[AB]+const_prob)/total; if (probHaps[AB] < 1e-100) probHaps[AB]=1e-100;
            probHaps[BA]=(numHaps[BA]+const_prob)/total; if (probHaps[BA] < 1e-100) probHaps[BA]=1e-100;
            probHaps[BB]=(numHaps[BB]+const_prob)/total; if (probHaps[BB] < 1e-100) probHaps[BB]=1e-100;
            
        }

        //now we have an initial reasonable guess at p, with which we can start the EM
        const_prob=0.0;
        count=1; loglike=-999999999.0;
        do {
            oldloglike=loglike;
            //count_haps(count);
            numHaps[AA] += unknownDH* (probHaps[AA]*probHaps[BB])/((probHaps[AA]*probHaps[BB])+(probHaps[AB]*probHaps[BA]));
            numHaps[BB] += unknownDH* (probHaps[AA]*probHaps[BB])/((probHaps[AA]*probHaps[BB])+(probHaps[AB]*probHaps[BA]));
            numHaps[AB] += unknownDH* (probHaps[AB]*probHaps[BA])/((probHaps[AA]*probHaps[BB])+(probHaps[AB]*probHaps[BA]));
            numHaps[BA] += unknownDH* (probHaps[AB]*probHaps[BA])/((probHaps[AA]*probHaps[BB])+(probHaps[AB]*probHaps[BA]));            
            
            loglike = (known[AA]*Math.log(probHaps[AA]) + known[AB]*Math.log(probHaps[AB]) + known[BA]*Math.log(probHaps[BA]) + known[BB]*Math.log(probHaps[BB]))/LN10 + ((double)unknownDH*Math.log(probHaps[AA]*probHaps[BB] + probHaps[AB]*probHaps[BA]))/LN10;
            if (Math.abs(loglike-oldloglike) < TOLERANCE) break;
            
            //estimate_p();
            double total= numHaps[AA]+numHaps[AB]+numHaps[BA]+numHaps[BB]+(4.0*const_prob);
            probHaps[AA]=(numHaps[AA]+const_prob)/total; if (probHaps[AA] < 1e-10) probHaps[AA]=1e-10;
            probHaps[AB]=(numHaps[AB]+const_prob)/total; if (probHaps[AB] < 1e-10) probHaps[AB]=1e-10;
            probHaps[BA]=(numHaps[BA]+const_prob)/total; if (probHaps[BA] < 1e-10) probHaps[BA]=1e-10;
            probHaps[BB]=(numHaps[BB]+const_prob)/total; if (probHaps[BB] < 1e-10) probHaps[BB]=1e-10;
            
            count++;
        } while(count < 1000);
        /* in reality I've never seen it need more than 10 or so iterations
        to converge so this is really here just to keep it from running off into eternity */

        loglike1 = (known[AA]*Math.log(probHaps[AA]) + known[AB]*Math.log(probHaps[AB]) + known[BA]*Math.log(probHaps[BA]) + known[BB]*Math.log(probHaps[BB]) + (double)unknownDH*Math.log(probHaps[AA]*probHaps[BB] + probHaps[AB]*probHaps[BA]))/LN10;
        loglike0 = (known[AA]*Math.log(pA1*pA2) + known[AB]*Math.log(pA1*pB2) + known[BA]*Math.log(pB1*pA2) + known[BB]*Math.log(pB1*pB2) + (double)unknownDH*Math.log(2*pA1*pA2*pB1*pB2))/LN10;

        num = probHaps[AA]*probHaps[BB] - probHaps[AB]*probHaps[BA];
        /*if(probHaps[AA] <= 1e-100 || probHaps[AB] <= 1e-100 || probHaps[BA] <= 1e-100 || probHaps[BB] <= 1e-100){
        	System.out.println("unknownDH = " + unknownDH + " AA=" + known[AA] + " BA=" + known[BA] + " BB=" + known[BB] + " AB=" + known[AB] + "\n" +
        			"Prob[AA]=" + probHaps[AA] + " Prob[AB]=" + probHaps[AB] + " Prob[BA]=" + probHaps[BA] + " Prob[BB]=" + probHaps[BB] + ".\n");
        }*/
        /*if (num < 0) {
            //flip matrix so we get the positive D'
            //flip AA with AB and BA with BB
            tmp=probHaps[AA]; probHaps[AA]=probHaps[AB]; probHaps[AB]=tmp;
            tmp=probHaps[BB]; probHaps[BB]=probHaps[BA]; probHaps[BA]=tmp;
            //flip frequency of second allele
            //done in this slightly asinine way because of a compiler bugz0r in the dec-alpha version of java
            //which causes it to try to parallelize the swapping operations and mis-schedules them
            pA2 = pA2 + pB2;
            pB2 = pA2 - pB2;
            pA2 = pA2 - pB2;
            //flip counts in the same fashion as p's
            tmp=numHaps[AA]; numHaps[AA]=numHaps[AB]; numHaps[AB]=tmp;
            tmp=numHaps[BB]; numHaps[BB]=numHaps[BA]; numHaps[BA]=tmp;
            //num has now undergone a sign change
            num = probHaps[AA]*probHaps[BB] - probHaps[AB]*probHaps[BA];
            if(num < 0){
            	System.out.println("negative D prime.\n");
            }
            //flip known array for likelihood computation
            tmp=known[AA]; known[AA]=known[AB]; known[AB]=tmp;
            tmp=known[BB]; known[BB]=known[BA]; known[BA]=tmp;
        }*/
        if(num>=0){
	        denom1 = (probHaps[AA]+probHaps[BA])*(probHaps[BA]+probHaps[BB]);
	        denom2 = (probHaps[AA]+probHaps[AB])*(probHaps[AB]+probHaps[BB]);
        }
        else{
	        denom1 = (probHaps[AB]+probHaps[BB])*(probHaps[BA]+probHaps[BB]);
	        denom2 = (probHaps[AA]+probHaps[AB])*(probHaps[AA]+probHaps[BA]);
        }
        denom = (denom1 <= denom2)? denom1:denom2;
        dprime = num/denom;

        /* add computation of r^2 = (D^2)/p(1-p)q(1-q) */
        rsq = num*num/(pA1*pB1*pA2*pB2);
        if(rsq > 1)/*&& (probHaps[AA] >= 0.001) && (probHaps[AB] >= 0.001) && (probHaps[BA] >= 0.001) && (probHaps[BB] >= 0.001)*/{
        	System.out.println("RSQ=" + rsq + " d=" + num + "\n pA1=" + pA1 + " pA2=" + pA2 + " pB1=" + pB1 + " pB2=" + pB2 + ".\n" +
        		"Prob[AA]=" + probHaps[AA] + " Prob[AB]=" + probHaps[AB] + " Prob[BA]=" + probHaps[BA] + " Prob[BB]=" + probHaps[BB] + ".\n");
        }

        //real_dprime=dprime;

        /*for (i=0; i<=100; i++) {
            dpr = (double)i*0.01;
            tmpAA = dpr*denom + pA1*pA2;
            tmpAB = pA1-tmpAA;
            tmpBA = pA2-tmpAA;
            tmpBB = pB1-tmpBA;
            if (i==100) {
                //one value will be 0
                if (tmpAA < 1e-10) tmpAA=1e-10;
                if (tmpAB < 1e-10) tmpAB=1e-10;
                if (tmpBA < 1e-10) tmpBA=1e-10;
                if (tmpBB < 1e-10) tmpBB=1e-10;
            }
            lsurface[i] = (known[AA]*Math.log(tmpAA) + known[AB]*Math.log(tmpAB) + known[BA]*Math.log(tmpBA) + known[BB]*Math.log(tmpBB) + (double)unknownDH*Math.log(tmpAA*tmpBB + tmpAB*tmpBA))/LN10;
        }*/

        /* Confidence bounds #2 - used in Gabriel et al (2002) - translate into posterior dist of D' -
        assumes a flat prior dist. of D' - someday we may be able to make
        this even more clever by adjusting given the distribution of observed
        D' values for any given distance after some large scale studies are complete */

        /*total_prob=sum_prob=0.0;

        for (i=0; i<=100; i++) {
            lsurface[i] -= loglike1;
            lsurface[i] = Math.pow(10.0,lsurface[i]);
            total_prob += lsurface[i];
        }

        for (i=0; i<=100; i++) {
            sum_prob += lsurface[i];
            if (sum_prob > 0.05*total_prob &&
                    sum_prob-lsurface[i] < 0.05*total_prob) {
                low_i = i-1;
                break;
            }
        }

        sum_prob=0.0;
        for (i=100; i>=0; i--) {
            sum_prob += lsurface[i];
            if (sum_prob > 0.05*total_prob &&
                    sum_prob-lsurface[i] < 0.05*total_prob) {
                high_i = i+1;
                break;
            }
        }
        if (high_i > 100){ high_i = 100; }


        double[] freqarray = {probHaps[AA], probHaps[AB], probHaps[BB], probHaps[BA]};*/

        return new PairwiseLinkage(Util.roundDouble(dprime,3), Util.roundDouble((loglike1-loglike0),2),
                Util.roundDouble(rsq,3)/*, ((double)low_i/100.0), ((double)high_i/100.0), freqarray*/);
    }

    public void saveDprimeToText(File dumpDprimeFile, Integer[] snpPos) throws IOException{
        FileWriter saveDprimeWriter = new FileWriter(dumpDprimeFile);
        //we use this LinkedList to store the dprime computations for a window of 5 markers
        //in either direction in the for loop farther down.
        //a tInt value is calculated for each marker which requires the dprime calculations
        //for the 5 markers before and 5 after the current marker, so we store these values while we need
        //them to avoid unnecesary recomputation.
        //long dist;
        saveDprimeWriter.write("L1\tL2\tD'\tLOD\tr^2\n");

        PairwiseLinkage currComp = null;

        for (int i = 0; i < snpPos.length; i++){
        	long posi = snpPos[i];
            for (int j = i+1; j < snpPos.length; j++){
                currComp = dpTable.getLDStats(i,j);                
                if(currComp != null) {
                    long posj = snpPos[j];
                    //dist =  posj - posi;
                    saveDprimeWriter.write(posi + "\t" + posj + "\t" + currComp.getDPrime() + "\t" + currComp.getLOD() + "\t" + 
                    		currComp.getRSquared() /*+ "\t" + currComp.getConfidenceLow() + "\t" + currComp.getConfidenceHigh() + "\t" + dist*/ + "\n");
                }
                else{
                	if(j==i+1)
                		saveDprimeWriter.write(posi + "\tnull\n");
                	break;
                }
            }
        }
        saveDprimeWriter.close();
    }

    /*public Haplotype[][] getHaplotypes() {
        return haplotypes;
    }*/

    /*public class RSquared {
        private double[] rsquareds;
        private double[] conditionalProbs;

        public RSquared(double[] rsquareds, double[] conditionalProbs) {
            this.rsquareds = rsquareds;
            this.conditionalProbs = conditionalProbs;
        }

        public double[] getRsquareds() {
            return rsquareds;
        }

        public double[] getConditionalProbs() {
            return conditionalProbs;
        }
    }*/

    /*public RSquared getPhasedRSquared(int snp, int[] block){

        double rsquareds[] = null;
        double conditionalProbs[] = null;
        double alleleCounts[][] = null;
        int maxIndex =0;
        int[] multiMarkerHaplos = null;
        boolean monomorphic = false;

        if(block.length == 2){
            multiMarkerHaplos = new int[8];

            int pos1 = snp;
            int pos2 = block[0];
            int pos3 = block[1];

            int[] marker1num = new int[5]; int[] marker2num = new int[5]; int[] marker3num = new int[5];

            markers.get(pos1);
			marker1num[SNP.getMajor()]=0;
            markers.get(pos1);
			marker1num[SNP.getMinor()]=4;
            markers.get(pos2);
			marker2num[SNP.getMajor()]=0;
            markers.get(pos2);
			marker2num[SNP.getMinor()]=2;
            markers.get(pos3);
			marker3num[SNP.getMajor()]=0;
            markers.get(pos3);
			marker3num[SNP.getMinor()]=1;

            alleleCounts = new double[3][5];

            byte a1,a2,a3,b1,b2,b3;
   
            for (int i = 0; i < chromosomes.length; i++){

                    a1 = chromosomes[i][pos1];
                    a2 = chromosomes[i][pos2];
                    a3 = chromosomes[i][pos3];
                    b1 = chromosomes[++i][pos1];
                    b2 = chromosomes[i][pos2];
                    b3 = chromosomes[i][pos3];

                    multiMarkerHaplos[marker1num[a1] + marker2num[a2] + marker3num[a3]]++;
                    multiMarkerHaplos[marker1num[b1] + marker2num[b2] + marker3num[b3]]++;
                    alleleCounts[0][a1]++;
                    alleleCounts[0][b1]++;
                    alleleCounts[1][a2]++;
                    alleleCounts[1][b2]++;
                    alleleCounts[2][a3]++;
                    alleleCounts[2][b3]++;

            }
            markers.get(pos1);
			markers.get(pos1);
			markers.get(pos2);
			markers.get(pos2);
			markers.get(pos3);
			markers.get(pos3);
			//check for any monomorphic SNPs
            if(alleleCounts[0][SNP.getMajor()] == 0 || alleleCounts[0][SNP.getMinor()] == 0
                    || alleleCounts[1][SNP.getMajor()] == 0 || alleleCounts[1][SNP.getMinor()] == 0
                    || alleleCounts[2][SNP.getMajor()] == 0 || alleleCounts[2][SNP.getMinor()] == 0){
                monomorphic = true;
            }
            maxIndex = 4;
        }else if (block.length == 3){
            multiMarkerHaplos = new int[16];

              int pos1 = snp;
            int pos2 = block[0];
            int pos3 = block[1];
            int pos4 = block[2];

            int[] marker1num = new int[5];
            int[] marker2num = new int[5];
            int[] marker3num = new int[5];
            int[] marker4num = new int[5];

			marker1num[SNP.getMinor()]=8;
			marker2num[SNP.getMinor()]=4;
			marker3num[SNP.getMinor()]=2;
			marker4num[SNP.getMinor()]=1;

            alleleCounts = new double[4][5];

            byte a1,a2,a3,a4,b1,b2,b3,b4;
            //iterate through all chromosomes in dataset
            for (int i = 0; i < chromosomes.length; i++){

                    a1 = chromosomes[i][pos1];
                    a2 = chromosomes[i][pos2];
                    a3 = chromosomes[i][pos3];
                    a4 = chromosomes[i][pos4];
                    b1 = chromosomes[++i][pos1];
                    b2 = chromosomes[i][pos2];
                    b3 = chromosomes[i][pos3];
                    b4 = chromosomes[i][pos4];

                    multiMarkerHaplos[marker1num[a1] + marker2num[a2] + marker3num[a3] + marker4num[a4]]++;
                    multiMarkerHaplos[marker1num[b1] + marker2num[b2] + marker3num[b3] + marker4num[b4]]++;

                    alleleCounts[0][a1]++;
                    alleleCounts[0][b1]++;
                    alleleCounts[1][a2]++;
                    alleleCounts[1][b2]++;
                    alleleCounts[2][a3]++;
                    alleleCounts[2][b3]++;
                    alleleCounts[3][a4]++;
                    alleleCounts[3][b4]++;
            }
            markers.get(pos1);
			markers.get(pos1);
			markers.get(pos2);
			markers.get(pos2);
			markers.get(pos3);
			markers.get(pos3);
			markers.get(pos4);
			markers.get(pos4);
			if(alleleCounts[0][SNP.getMajor()] == 0 || alleleCounts[0][SNP.getMinor()] == 0
                    || alleleCounts[1][SNP.getMajor()] == 0 || alleleCounts[1][SNP.getMinor()] == 0
                    || alleleCounts[2][SNP.getMajor()] == 0 || alleleCounts[2][SNP.getMinor()] == 0
                    || alleleCounts[3][SNP.getMajor()] == 0 || alleleCounts[3][SNP.getMinor()] == 0){
                monomorphic = true;
            }
            maxIndex =8;
        }
        //the rest of the code is the same for 2 and 3 marker blocks
        rsquareds = new double[maxIndex];
        conditionalProbs = new double[maxIndex];

        if(monomorphic){
            Arrays.fill(rsquareds,0);
            Arrays.fill(conditionalProbs,0);
            return new RSquared(rsquareds,conditionalProbs);
        }

        int totalChroms=0;

        for(int i = 0;i < multiMarkerHaplos.length;i++){
            totalChroms += multiMarkerHaplos[i];
        }

        double[] freqs = new double[multiMarkerHaplos.length];
        for(int i=0;i<freqs.length;i++){
            freqs[i] = multiMarkerHaplos[i]/(double)totalChroms;
        }

        double p=0;
        for(int i=0;i< maxIndex; i++){
            p += freqs[i];
        }

        for(int i =0;i< maxIndex;i++){
            //calculate r^2
            double aa = freqs[i];
            double ab = p - freqs[i];
            double ba = freqs[i+maxIndex];
            double bb = (1-p) - freqs[i+maxIndex];

            double q = ba + aa;
            double c = aa*bb - ab*ba;
            rsquareds[i] = Util.roundDouble((c*c)/(p*(1-p)*q*(1-q)),3);

            //calculate conditional prob (ie P(snp | hap))
            conditionalProbs[i] = freqs[i]/(freqs[i] + freqs[i+maxIndex]);
        }

        return new RSquared(rsquareds,conditionalProbs);
    }*/


    /*public static boolean isPhasedData() {
        return phasedData;
    }*/
}
