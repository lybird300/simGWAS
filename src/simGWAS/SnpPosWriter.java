package simGWAS;

import java.io.File;

public class SnpPosWriter implements Runnable {
	private int firstSample;
	private int lastSample;
	private byte[][] hapArray;
	private String outDir;
	private String outFile;
	//private String inDir;
	//private String inFile;
    
    SnpPosWriter(int start, int end, byte[][] haps, String outDir, String outFile/*, String dataInputDir, String hapInput*/){
    	super();
    	firstSample = start;
    	lastSample = end;
    	hapArray = haps;
    	this.outDir = outDir;
    	this.outFile = outFile;
    	//inDir = dataInputDir;
    	//inFile = hapInput;
    }
    
	public void run() {
		try {
			int doneLine = 0;
			File snpPosFile = new File(outDir, outFile);
	   		FileHandler processMutInput = null;
       		if(!snpPosFile.exists()){
       			processMutInput = new FileHandler(snpPosFile);
       			processMutInput.openForWrite(true);
       			processMutInput.writeString("ChromID\tMutation Position.....");//output header
			}
       		else{//already exists; go to the end of this file; check whether it is complete; if not return the last chromID; start a new line
       			processMutInput = new FileHandler(snpPosFile);
       			processMutInput.openForRead();
       			processMutInput.nextLine();//header
       			while(processMutInput.nextLine()!=null){
       				doneLine++;
   				}
       			processMutInput.closeReader();
       			processMutInput = new FileHandler(snpPosFile);
       			processMutInput.openForWrite(true);
       			//processMutInput.writer.newLine();//This will create a new line between existing records and new appended records
       		}
       		/*FileHandler readHapInput = new FileHandler(inDir, inFile);
       		readHapInput.openForRead();       			
       		String input = null;
   			int readLine = 0;
   			while((input = readHapInput.nextLine())!=null && readLine <= lastSample){
   				if(readLine < firstSample){
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
	   		String output = null;
	   		for(int i = firstSample; i <= lastSample; i++){
	   			if(doneLine > 0){
	   				doneLine--;
	   				continue;
	   			}
	   			output = String.valueOf(i);
	   			for(int j = 0; j < hapArray[0].length; j++){		
					if(hapArray[i][j]==(byte)1)
						output += "\t" + j;
				}
	   			processMutInput.writeString(output);
	   			output = null;
	   		}
			processMutInput.closeWriter();
			output = null;
			Thread.sleep(1000);
		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
