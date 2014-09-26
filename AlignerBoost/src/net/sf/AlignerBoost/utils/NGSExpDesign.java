/**
 * A class to represent the NGS experimental design file information for a given library
 */
package net.sf.AlignerBoost.utils;
import java.io.*;
import java.util.*;
import java.util.regex.*;

/**
 * @author Qi Zheng
 * @version 1.1.2
 * @since 1.1
 */
public class NGSExpDesign {
	/**
	 * A static method to create a List of NGSExpDesign objects from a File
	 * @param designFileName
	 * @return a List of NGSExpDesign objects
	 * @throws IOException
	 * @throws IllegalArgumentException
	 */
	public static List<NGSExpDesign> createNGSExpDesignsFromFile(String designFileName) throws IOException, IllegalArgumentException {
		// initiate lists
		List<NGSExpDesign> designLs = new ArrayList<NGSExpDesign>();
		optNames = new ArrayList<String>();
		
		// read expDesign file
		BufferedReader in = new BufferedReader(new FileReader(designFileName));
		String line = null;
		while((line = in.readLine()) != null) {
			if(line.startsWith("#")) { // a header line
				Matcher match = globalPat.matcher(line);
				if(match.find()) { // a global opt line
					switch(match.group(1)) {
					case "MAX_PROC":
						MAX_PROC = Integer.parseInt(match.group(2));
					default:
							throw new IllegalArgumentException("Uknown global option '" + match.group(1) + "' found at\n" + line);
					}
				}
				match = localPat.matcher(line); // try local opts
				if(match.find()) // a local opt line
					optNames.add(match.group(1));
			}
			else { // a opt-value line
				String[] optValues = line.split("\t");
				if(optNames.size() != optValues.length)
					throw new IllegalArgumentException("Incorrect number of value fields found at\n$line\nExpect " +
							optNames.size() + " fields but found " + optValues.length);
				NGSExpDesign design = new NGSExpDesign();
				// construct NGSExpDesign objects with paired opt names and values
				for(int i = 0; i < optNames.size(); i++) {
					String name = optNames.get(i);
					String value = optValues[i];
					switch(name) {
					case "libname":
						design.libName = value;
						break;
					case "read_file":
						design.readFile = value;
						break;
					case "read_len":
						design.readLen = Integer.parseInt(value);
						break;
					case "ascii_offset":
						design.asciiOffset = Integer.parseInt(value);
						break;
					case "is_paired":
						design.isPaired = value.toUpperCase().equals("YES");
						break;
					case "mate_file":
						design.mateFile = value;
						break;
					case "strand_type":
						design.strandType = Integer.parseInt(value);
						break;
					case "has_spliced":
						design.hasSpliced = value.toUpperCase().equals("YES");
						break;
					case "do_trim":
						design.doTrim = value.toUpperCase().equals("YES");
						break;
					case "3_adapter_seq":
						design.adapterSeq3 = value;
						break;
					case "5_adapter_seq":
						design.adapterSeq5 = value;
						break;
					case "trim_mis":
						design.trimMis = Float.parseFloat(value) / 100;
						break;
					case "min_trim":
						design.minTrim = Integer.parseInt(value);
						break;
					case "do_NR":
						design.doNR = value.toUpperCase().equals("YES");
						break;
					case "aligner":
						design.aligner = value;
						break;
					case "seed_len":
						design.seedLen = Integer.parseInt(value);
						break;
					case "seed_mis":
						design.seedMis = Float.parseFloat(value) / 100;
						break;
					case "all_mis":
						design.allMis = Float.parseFloat(value) / 100;
						break;
					case "all_indel":
						design.allIndel = Float.parseFloat(value) / 100;
						break;
					case "min_insert":
						design.minInsert = Integer.parseInt(value);
						break;
					case "max_hit":
						design.maxHit = Integer.parseInt(value);
						break;
					case "min_frag_len":
						design.minFragLen = Integer.parseInt(value);
						break;
					case "max_frag_len":
						design.maxFragLen = Integer.parseInt(value);
						break;
					case "max_div":
						design.maxDiv = Float.parseFloat(value) / 100;
						break;
					case "max_best":
						design.maxBest = Integer.parseInt(value);
						break;
					case "max_report":
						design.maxReport = Integer.parseInt(value);
						break;
					case "ref_genome":
						design.refGenome = value;
						break;
					case "ref_index":
						design.refIndex = value;
						break;
					case "transcriptome_GFF":
						design.transcriptomeGFF = value;
						break;
					case "transcriptome_index":
						design.transcriptomeIndex = value;
						break;
					case "other_aligner_opts":
						design.otherAlignerOpts = value;
						break;
					default:
						throw new IllegalArgumentException("Unknown per-lib option '" + name + "' found at\n" + line);							
					}
				} // end of each value field
				// fill default value for this libOpts
				design.setDefaultOptValues();
				designLs.add(design);
			} // end of each value line (design)
			in.close();
		}
		return designLs;
	}
	
	/**
	 * @return the MAX_PROC
	 */
	public static int getMAX_PROC() {
		return MAX_PROC;
	}

	/**
	 * @param MAX_PROC the mAX_PROC to set
	 */
	public static void setMAX_PROC(int max) {
		MAX_PROC = max;
	}

	/**
	 * @return the libName
	 */
	public String getLibName() {
		return libName;
	}

	/**
	 * @return the readFile
	 */
	public String getReadFile() {
		return readFile;
	}

	/**
	 * @return the readLen
	 */
	public int getReadLen() {
		return readLen;
	}

	/**
	 * @return the asciiOffset
	 */
	public int getAsciiOffset() {
		return asciiOffset;
	}

	/**
	 * @return the isPaired
	 */
	public boolean isPaired() {
		return isPaired;
	}

	/**
	 * @return the mateFile
	 */
	public String getMateFile() {
		return mateFile;
	}

	/**
	 * @return the strandType
	 */
	public int getStrandType() {
		return strandType;
	}

	/**
	 * @return the hasSpliced
	 */
	public boolean hasSpliced() {
		return hasSpliced;
	}

	/**
	 * @return the doTrim
	 */
	public boolean doTrim() {
		return doTrim;
	}

	/**
	 * @return the adapterSeq3
	 */
	public String getAdapterSeq3() {
		return adapterSeq3;
	}

	/**
	 * @return the adapterSeq5
	 */
	public String getAdapterSeq5() {
		return adapterSeq5;
	}

	/**
	 * @return the trimMis
	 */
	public float getTrimMis() {
		return trimMis;
	}

	/**
	 * @return the minTrim
	 */
	public int getMinTrim() {
		return minTrim;
	}

	/**
	 * @return the doNR
	 */
	public boolean doNR() {
		return doNR;
	}

	/**
	 * @return the aligner
	 */
	public String getAligner() {
		return aligner;
	}

	/**
	 * @return the seedLen
	 */
	public int getSeedLen() {
		return seedLen;
	}

	/**
	 * @return the seedMis
	 */
	public float getSeedMis() {
		return seedMis;
	}

	/**
	 * @return the allMis
	 */
	public float getAllMis() {
		return allMis;
	}

	/**
	 * @return the allIndel
	 */
	public float getAllIndel() {
		return allIndel;
	}

	/**
	 * @return the minInsert
	 */
	public int getMinInsert() {
		return minInsert;
	}

	/**
	 * @return the maxHit
	 */
	public int getMaxHit() {
		return maxHit;
	}

	/**
	 * @return the minFragLen
	 */
	public int getMinFragLen() {
		return minFragLen;
	}

	/**
	 * @return the maxFragLen
	 */
	public int getMaxFragLen() {
		return maxFragLen;
	}

	/**
	 * @return the otherAlignerOpts
	 */
	public String getOtherAlignerOpts() {
		return otherAlignerOpts;
	}

	/**
	 * @return the maxDiv
	 */
	public float getMaxDiv() {
		return maxDiv;
	}

	/**
	 * @return the maxBest
	 */
	public int getMaxBest() {
		return maxBest;
	}

	/**
	 * @return the maxReport
	 */
	public int getMaxReport() {
		return maxReport;
	}

	/**
	 * @return the refGenome
	 */
	public String getRefGenome() {
		return refGenome;
	}

	/**
	 * @return the refIndex
	 */
	public String getRefIndex() {
		return refIndex;
	}

	/**
	 * @return the transcriptomeGFF
	 */
	public String getTranscriptomeGFF() {
		return transcriptomeGFF;
	}

	/**
	 * @return the transcriptomeIndex
	 */
	public String getTranscriptomeIndex() {
		return transcriptomeIndex;
	}

	/**
	 * @return the optNames
	 */
	public static List<String> getOptNames() {
		return optNames;
	}

	private void setDefaultOptValues() {
		if(readFile.equals("NA"))
			readFile = libName + ".fastq";
		if(refIndex.equals("NA"))
			refIndex = refGenome;
		if(transcriptomeIndex.equals("NA"))
			transcriptomeIndex = "transcriptome/" + transcriptomeGFF.replaceFirst("(?i:\\.gff)", "");		
		if(!doTrim && !hasSpliced)
			minInsert = readLen;
	}

	
	// global options
	public static int MAX_PROC = 6; // maximum processors to use
	
	// basic library options
	private String libName;
	private String readFile;
	private int readLen;
	private int asciiOffset = 33;
	private boolean isPaired;
	private String mateFile;
	private int strandType;
	private boolean hasSpliced;
	// trimming options
	private boolean doTrim;
	private String adapterSeq3;
	private String adapterSeq5;
	private float trimMis = 0.1f;
	private int minTrim = 10;
	// NR options
	private boolean doNR = true;
	// mapping options
	private String aligner;
	private int seedLen = 25;
	private float seedMis = 0.04f;
	private float allMis = 0.06f;
	private float allIndel = 0;
	private int minInsert = 15;
	private int maxHit = 10;
	// mate options
	private int minFragLen = 0;
	private int maxFragLen = 500;
	// user-specified options
	private String otherAlignerOpts = "";
	// best-stratum options
	private float maxDiv = 0.04f;
	private int maxBest = 0;
	private int maxReport = 0;
	// reference genome options
	private String refGenome;
	private String refIndex;
	// transcriptome options
	private String transcriptomeGFF;
	private String transcriptomeIndex;
	
	private static List<String> optNames;
	private static Pattern globalPat = Pattern.compile("^## (\\w+)=(.*):");
	private static Pattern localPat = Pattern.compile("^### (\\w+):");
}
