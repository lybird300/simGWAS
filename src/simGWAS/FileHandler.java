package simGWAS;

import java.io.FileOutputStream;
import java.util.List;
import java.util.zip.GZIPOutputStream;
import java.io.OutputStreamWriter;
import java.io.RandomAccessFile;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.BufferedWriter;
import java.util.zip.GZIPInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;

public class FileHandler {
	
	public BufferedReader reader = null;
    
    public BufferedWriter writer = null;

    private String filename = null;
    
    private String directory = null;
    
    private boolean useGzip = false;
    
    private boolean newFile = false;
    
    private File f = null;    

    public FileHandler(String _d, String _f){
    	try{
    		filename = _f;
    		directory = _d;
    		f = new File(directory + File.separator + filename);
    		newFile = f.createNewFile();
       		String absolutePath = f.getAbsolutePath();
    		directory = absolutePath.substring(0,absolutePath.lastIndexOf(File.separator));//exclude the separator in the end
    		filename = f.getName();
    		useGzip = isGZipped(f, newFile);
    	}catch (IOException e) {
    		System.out.println(directory + filename + " file could not be opened");
    	} 
    }
    
    public void preventFutureWriting(){
    	f.setReadOnly();
    }
    
    public void allowFutureWriting(){
    	f.setWritable(true);
    }
    
    /**
     * Checks if a file is gzipped.
     */
    public static boolean isGZipped(File f, boolean newFile) {
    	if(newFile)
    		return f.getName().endsWith(".gz");
    	else{
		     int magic = 0;
		     try {
		      RandomAccessFile raf = new RandomAccessFile(f, "r");
		      magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
		      raf.close();
		     } catch (Throwable e) {
		      e.printStackTrace(System.err);
		     }
		     return magic == GZIPInputStream.GZIP_MAGIC;
    	}
    }
    
    public FileHandler(String completeFileName){
    	try{
    		f = new File(completeFileName);
     		newFile = f.createNewFile();   		
    		String absolutePath = f.getAbsolutePath();
    		directory = absolutePath.substring(0,absolutePath.lastIndexOf(File.separator));//exclude the separator in the end
    		filename = f.getName();
    		useGzip = isGZipped(f, newFile);
    	}catch (IOException e) {
    		System.out.println(directory + filename + " file could not be opened");
    	} 
    }
    
    public FileHandler(File file){
    	try{
    		f = file;
    		newFile = f.createNewFile();
       		String absolutePath = f.getAbsolutePath();
    		directory = absolutePath.substring(0,absolutePath.lastIndexOf(File.separator));//exclude the separator in the end
    		filename = f.getName();
    		useGzip = isGZipped(f, newFile);
    	}catch (IOException e) {
    		System.out.println(f.getAbsolutePath() + " could not be opened");
    	} 
    }
    
    public boolean isNewFile(){
    	return newFile;
    }
    
    public String getAbsolutePath(){
    	return f.getAbsolutePath();
    }
    
    public void openForRead(){   	
    	try{
    		FileInputStream fis = new FileInputStream(f);
    		if(useGzip){   			
    			GZIPInputStream gzipStream = new GZIPInputStream(fis);          	   
    			reader = new BufferedReader(new InputStreamReader(gzipStream, Charset.forName("UTF-8").newDecoder()));    			
    		}
    		else{    			
    			reader = new BufferedReader(new InputStreamReader(fis, Charset.forName("UTF-8").newDecoder()));
    		}
    	} catch (FileNotFoundException e) {
    		System.out.println(directory + filename + " file could not be found");
    	} catch (IOException e) {
    		System.out.println(directory + filename + " file could not be opened");
    	}    	
    }
    
    public void closeReader() {
        try {
           	reader.close();
            reader = null;
            f = null;
        } catch (IOException e) {
            System.out.println("Trouble closing buffered file ");
            e.printStackTrace();
            System.exit(0);
        }
    }
    
    public void closeWriter() {
        try {
        	writer.flush();
            writer.close();
            f = null;
        } catch (IOException e) {
            System.out.println("Trouble closing buffered file ");
            e.printStackTrace();
        }
    }

    public String nextLine() {
        try {
            return reader.readLine();
        } catch (IOException e) {
            System.out.println("An I/O error has occured ");
            e.printStackTrace();
            System.exit(0);
        }
        return null;
    }
   
    public void open(String path) {
        directory = path;
        this.openForWrite(false);
    }
    
    public void openWithAppend(String path) {
        directory = path;
        this.openForWrite(true);
    }
    
    public void openForWrite(boolean append) {
        try {
        	FileOutputStream fos = new FileOutputStream(f, append);
        	OutputStreamWriter outWriter = null;
            if(useGzip){            	
            	GZIPOutputStream gzipStream = new GZIPOutputStream(fos);
            	outWriter = new OutputStreamWriter(gzipStream, Charset.forName("UTF-8").newEncoder());
            }
            else{
            	outWriter = new OutputStreamWriter(fos, Charset.forName("UTF-8").newEncoder());
            }
            writer = new BufferedWriter(outWriter);                  
        } catch (IOException e) {
            System.out.println(filename + " file could not be opened");
            e.printStackTrace();
            System.exit(-1);
        }
    }

    public void writeString(String lineToWrite) {
        try {
            writer.append(lineToWrite);
            writer.newLine();
        } catch (IOException e) {
            System.out.println("An I/O error has occured ");
            System.exit(-1);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

	public void writeString(String writeString, boolean echo) {
		writeString(writeString);
		if(echo)System.out.println(writeString);		
	}
	
	public void writeSubjectFile(List<Subject> cases, List<Subject> controls){
		try{
			//output the head row the first time we write into this file
			if(newFile){
			    writer.append("SubjectID");
			    writer.append(',');
			    writer.append("ChromA_ID");
			    writer.append(',');
			    writer.append("ChromB_ID");
			    writer.append(',');
			    writer.append("CaseOrNot");//0--control(unaffected), 1--case(affected)
			    writer.append(',');
			    writer.append("CarrierOrNot");//0--non-carrier, 1--carrier
			    writer.newLine();
			}
			for(Subject sub : cases){
				writer.append(sub.getID());
				writer.append(',');
				writer.append(String.valueOf(sub.getChromA()));
				writer.append(',');
				writer.append(String.valueOf(sub.getChromB()));
				writer.append(',');
				writer.append(String.valueOf(1));
				writer.append(',');
				writer.append(String.valueOf(sub.isCarrier()));
				writer.newLine();
			}
			for(Subject sub : controls){
				writer.append(sub.getID());
				writer.append(',');
				writer.append(String.valueOf(sub.getChromA()));
				writer.append(',');
				writer.append(String.valueOf(sub.getChromB()));
				writer.append(',');
				writer.append(String.valueOf(0));
				writer.append(',');
				writer.append(String.valueOf(sub.isCarrier()));
				writer.newLine();
			}
		}catch (IOException e) {
            System.out.println("An I/O error has occured ");
            System.exit(-1);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
	}

	public String getFilename() {
		return filename;
	}

	public void setFilename(String filename) {
		this.filename = filename;
	}
	public String getDirectory() {
		return directory;
	}

	public void setDirectory(String d) {
		this.directory = d;
	}

	public boolean canWrite() {
		return f.canWrite();
	}

	public File getFile() {
		return f;
	}
}
