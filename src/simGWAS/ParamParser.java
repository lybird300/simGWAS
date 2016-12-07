package simGWAS;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;

public class ParamParser {
	public static String propertyFile;
	public static Properties props = null;
	//public static String workDir;
	public static String dataOutputDir;
	public static String dataInputDir;
	public static String hapInput;
	public static String mutInput;
	public static String posInput;
	public static int inputFileNum;
	public static double[] lo_maf = new double[2];
	public static double[] hi_maf = new double[2];
	public static double markerMAFTresh;
	public static double[] grr;
	public static double[] bdr;
	public static int[] num_ca;
	public static int[] num_co;
	public static int replica;
	public static int sampleSize;
	public static int num_mut;
	public static int numOfCVPerMAFRange;
	public static int segmentStart;
	public static int segmentEnd;
	public static int chromLength;
	public static int num_mark;
	public static int num_refSub;
	public static int[] num_forwardGen;
	public static double prob_recombPerChrom;
	public static int founderChrom;
	public static double missingRate;
	public static double errorRate;
	
	ParamParser(String arg){
		ParamParser.propertyFile = arg;
		setParamValue();
	}
	private static void setParamValue(){
		//if (args.length == 1) {
	        System.out.println("********** load the property file **********\n" + propertyFile + "\n********************");
	        File currentPropertiesFile = new File(propertyFile);
	        Properties props = new Properties();
	        try {
	            props.load(new FileInputStream(currentPropertiesFile));
	        } catch (FileNotFoundException e) {
	            e.printStackTrace();
	            System.exit(-1);
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
	        getParamValueFromPropertyFile(props);	        	
	        //focus on the 1/4 segment in the middle of the haplotype
	        segmentStart = (int) (0 + chromLength*(1.0/2.0-1.0/8.0));
	        segmentEnd = (int) (segmentStart + chromLength/4.0);    				        	
		//}
	}
	private static void getParamValueFromPropertyFile(Properties props) {
		//workDir=System.getProperty("user.dir");
		dataOutputDir = props.getProperty("dataOutputDir");
		dataInputDir = props.getProperty("dataInputDir");
		hapInput = props.getProperty("InputHapFileName");
		mutInput = props.getProperty("InputMutFileName");
		posInput = props.getProperty("InputPosFileName");
		inputFileNum = Integer.parseInt(props.getProperty("NumberOfInputFiles"));
		String[] maf_range = props.getProperty("MinAlleleFreqRange").split(",");
		lo_maf[0] = Double.parseDouble(maf_range[0]);
		lo_maf[1] = Double.parseDouble(maf_range[2]);
		hi_maf[0] = Double.parseDouble(maf_range[1]);
		hi_maf[1] = Double.parseDouble(maf_range[3]);
		markerMAFTresh = Double.parseDouble(props.getProperty("MinAlleleFreqOfMarkers"));
		num_mark = Integer.parseInt(props.getProperty("NumberOfMarkers"));
		String[] grr_string = props.getProperty("GenoRelativeRisk").split(",");
		grr = new double[grr_string.length];
		for(int i = 0; i<grr_string.length; i++)
			grr[i] = Double.parseDouble(grr_string[i]);
		String[] bdr_string = props.getProperty("BaseDiseaseRisk").split(",");
		bdr = new double[bdr_string.length];
		for(int i = 0; i<bdr_string.length; i++)
			bdr[i] = Double.parseDouble(bdr_string[i]);
		String[] num_ca_string = props.getProperty("NumberOfCases").split(",");
		num_ca = new int[num_ca_string.length];
		for(int j = 0; j<num_ca_string.length; j++)
			num_ca[j] = Integer.parseInt(num_ca_string[j]);
		//num_ca = Integer.parseInt(props.getProperty("NumberOfCases"));
		String[] num_co_string = props.getProperty("NumberOfControls").split(",");
		num_co = new int[num_co_string.length];
		for(int j = 0; j<num_co_string.length; j++)
			num_co[j] = Integer.parseInt(num_co_string[j]);
		replica = Integer.parseInt(props.getProperty("NumberOfReplicates"));
		sampleSize = Integer.parseInt(props.getProperty("SampleSize"));
		num_mut = Integer.parseInt(props.getProperty("NumberOfMutations"));
		numOfCVPerMAFRange = Integer.parseInt(props.getProperty("NumberOfCVPerRange"));
		chromLength = Integer.parseInt(props.getProperty("LengthOfChrom"));
		num_refSub = Integer.parseInt(props.getProperty("NumberOfSubsInRefPanel"));
		String[] gen_string = props.getProperty("NumOfForwardGenerations").split(",");
		num_forwardGen = new int[gen_string.length];
		for(int k = 0; k<gen_string.length; k++)
			num_forwardGen[k] = Integer.parseInt(gen_string[k]);
		prob_recombPerChrom = Double.parseDouble(props.getProperty("ProbOfRecombPerChromPerGen"));
		founderChrom = Integer.parseInt(props.getProperty("FounderChrom"));
		//missingRate = Double.parseDouble(props.getProperty("MissingRatio"));
		//errorRate = Double.parseDouble(props.getProperty("ErrorRate"));
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
