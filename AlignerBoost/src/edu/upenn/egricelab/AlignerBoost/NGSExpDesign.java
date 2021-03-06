/*******************************************************************************
 *     This file is part of AlignerBoost, a generalized software toolkit to boost
 *     the NextGen sequencing (NGS) aligner precision and sensitivity.
 *     Copyright (C) 2015  Qi Zheng
 *
 *     AlignerBoost is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     AlignerBoost is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
/**
 * A class to represent the NGS experimental design file information for a given library
 */
package edu.upenn.egricelab.AlignerBoost;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Qi Zheng
 * @version 1.2
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
	public static List<NGSExpDesign> createNGSExpDesignsFromFile(String designFileName)
			throws IOException, IllegalArgumentException {
		// initiate lists
		List<NGSExpDesign> designLs = new ArrayList<NGSExpDesign>();
		optNames = new ArrayList<String>();
		
		// read expDesign file
		BufferedReader in = new BufferedReader(new FileReader(designFileName));
		String line = null;
		while((line = in.readLine()) != null) {
			if(line.isEmpty())
				continue;
			if(line.startsWith("#")) { // a header line
				Matcher match = globalPat.matcher(line);
				if(match.find()) { // a global opt line
					String name = match.group(1);
					switch(name) {
					case "VERSION":
						VERSION = match.group(2);
						if(cmpVersion(VERSION, progVer) > 0)
							throw new IllegalArgumentException("Config file version is newer than the progVersion, please download the latest " + progName);
						break;
					case "MAX_PROC":
						MAX_PROC = Integer.parseInt(match.group(2));
						break;
					case "INIT_MEM":
						INIT_MEM = match.group(2);
						break;
					case "MAX_MEM":
						MAX_MEM = match.group(2);
						break;
					case "SH_PATH":
						SH_PATH = match.group(2);
						break;
					case "PROJECT_DIR":
						PROJECT_DIR = match.group(2).replace("/$", "");
						break;
					case "WORK_DIR":
						WORK_DIR = match.group(2).replace("/$", "");
						break;
					default:
						in.close();
						throw new IllegalArgumentException("Uknown global option '" + match.group(1) + "' found at\n" + line);
					}
				}
				match = localPat.matcher(line); // try local opts
				if(match.find()) // a local opt line
					optNames.add(match.group(1));
			}
			else { // a opt-value line
				String[] optValues = line.split("\t");
				if(optNames.size() < optValues.length) {
					in.close();
					throw new IllegalArgumentException("Incorrect number of value fields found at" + newLine + line + newLine + "Expect at most " +
							optNames.size() + " fields but found " + optValues.length);
				}
				NGSExpDesign design = new NGSExpDesign();
				// construct NGSExpDesign objects with paired opt names and values
				for(int i = 0; i < optValues.length; i++) {
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
						if(!(design.asciiOffset == 0 || design.asciiOffset == 33 || design.asciiOffset == 64))
							throw new IllegalArgumentException("option 'ascii_offset' must be 33 or 64 or 0 (not recommended)");
						break;
					case "is_paired":
						design.isPaired = value.toUpperCase().equals("YES");
						break;
					case "mate_file":
						design.mateFile = value;
						break;
					case "mate_len":
						design.mateLen = Integer.parseInt(value);
						if(design.mateLen == 0)
							design.mateLen = design.readLen; // default mate_len == read_len
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
					case "trim_prog":
						design.trimProg = value.toLowerCase();
						break;
					case "3_adapter_seq":
						design.adapterSeq3 = value;
						break;
					case "5_adapter_seq":
						design.adapterSeq5 = value;
						break;
					case "trim_mis":
						design.trimMis = Float.parseFloat(value);
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
						design.seedMis = Float.parseFloat(value);
						break;
					case "seed_indel":
						design.seedIndel = Float.parseFloat(value);
						break;
					case "all_mis":
						design.allMis = Float.parseFloat(value);
						break;
					case "all_indel":
						design.allIndel = Float.parseFloat(value);
						break;
					case "min_align_rate":
						design.setMinAlignRate(Float.parseFloat(value));
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
					case "min_mapQ":
						design.minMapQ = Integer.parseInt(value);
						break;
					case "max_best":
						design.maxBest = Integer.parseInt(value);
						break;
					case "max_report":
						design.maxReport = Integer.parseInt(value);
						break;
					case "best_mode":
						if(value.equals("best")) {
							design.maxBest = 0;
							design.maxReport = 1;
						}
						else if(value.equals("best-only")) {
							design.maxBest = 1;
							design.maxReport = 1;
						}
						else if(value.equals("none")) { } // do nothing
						else
							throw new IllegalArgumentException("Unknown best_mode option '" + value + "' found at\n" + line);
						break;
					case "sort_method":
						design.sortMethod = value;
						break;
					case "known_SNP_file":
						if(!value.equals("NA"))
							design.knownSnpFile = value;
						break;
					case "ref_genome":
						design.refGenome = value;
						break;
					case "ref_index":
						design.refIndex = value;
						break;
					case "transcriptome_GFF":
						if(!value.equals("NA"))
							design.transcriptomeGFF = value;
						break;
					case "transcriptome_index":
						design.transcriptomeIndex = value;
						break;
					case "other_aligner_opts":
						design.otherAlignerOpts = value;
						break;
					case "other_filter_opts":
						design.otherFilterOpts = value;
						break;
					default:
						throw new IllegalArgumentException("Unknown per-lib option '" + name + "' found at\n" + line);							
					}
				} // end of each value field
				// fill default value for this libOpts
				design.setDefaultOptValues();
				designLs.add(design);
			} // end of each value line (design)
		}
		in.close();
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
	 * @return the iNIT_MEM
	 */
	public static String getINIT_MEM() {
		return INIT_MEM;
	}

	/**
	 * @param iNIT_MEM the iNIT_MEM to set
	 */
	public static void setINIT_MEM(String iNIT_MEM) {
		INIT_MEM = iNIT_MEM;
	}

	/**
	 * @return the mAX_MEM
	 */
	public static String getMAX_MEM() {
		return MAX_MEM;
	}

	/**
	 * @param mAX_MEM the mAX_MEM to set
	 */
	public static void setMAX_MEM(String maxMEM) {
		MAX_MEM = maxMEM;
	}

	/**
	 * @return the sH_PATH
	 */
	public static String getSH_PATH() {
		return SH_PATH;
	}

	/**
	 * @param sH_PATH the sH_PATH to set
	 */
	public static void setSH_PATH(String path) {
		SH_PATH = path;
	}

	/**
	 * @return the PROJECT_DIR
	 */
	public static String getPROJECT_DIR() {
		return PROJECT_DIR;
	}

	/**
	 * @param dir the PROJECT_DIR to set
	 */
	public static void setPROJECT_DIR(String dir) {
		PROJECT_DIR = dir;
	}

	/**
	 * @return the WORK_DIR
	 */
	public static String getWORK_DIR() {
		return WORK_DIR;
	}

	/**
	 * @param dir the WORK_DIR to set
	 */
	public static void setWORK_DIR(String dir) {
		WORK_DIR = dir;
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
	 * @return the mateLen
	 */
	public int getMateLen() {
		return mateLen;
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
	 * @return  the trimProg
	 */
	public String getTrimProg() {
		return trimProg;
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
	 * @return the minAlignRate
	 */
	public float getMinAlignRate() {
		return minAlignRate;
	}

	/**
	 * @param minAlignRate the minAlignRate to set
	 */
	public void setMinAlignRate(float minAlignRate) {
		// check minAlignRate
		if(minAlignRate < 0)
			minAlignRate = !hasSpliced ? DEFAULT_DNA_ALIGN_RATE :
				isRNAAligner(aligner) ? DEFAULT_RNA_ALIGN_RATE_RNA_ALIGNER : DEFAULT_RNA_ALIGN_RATE_DNA_ALIGNER;
		this.minAlignRate = minAlignRate;
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
	 * @return the otherFilterOpts
	 */
	public String getOtherFilterrOpts() {
		return otherFilterOpts;
	}
	
	/**
	 * @return the minMapQ
	 */
	public int getMinMapQ() {
		return minMapQ;
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
	 * @return the sortMethod
	 */
	public String getSortMethod() {
		return sortMethod;
	}

	/**
	 * @return the knownSnpFile
	 */
	public String getKnownSnpFile() {
		return knownSnpFile;
	}

	/**
	 * @param knownSnpFile the knownSnpFile to set
	 */
	public void setKnown_SNP_file(String knownSnpFile) {
		this.knownSnpFile = knownSnpFile;
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
		if(transcriptomeGFF != null && transcriptomeIndex.equals("NA"))
			transcriptomeIndex = "transcriptome/" + transcriptomeGFF.replaceFirst("(?i:\\.gff)", "");		
	}

	/**
	 * get trimmed read file name
	 * @return  trimmed read filename
	 */
	public String getTrimmedReadFileName() {
		if(!doTrim)
			return getReadFile();
		String fn = !isPaired ? libName + "_trimmed.fastq" : libName + "_trimmed_1.fastq";
		return WORK_DIR.equals(".") ? fn : WORK_DIR + "/" + fn;
	}

	/**
	 * get trimmed mate file name
	 * @return  trimmed mate filename
	 */
	public String getTrimmedMateFileName() {
		if(!doTrim)
			return getMateFile();
		String fn = !isPaired ? "" : libName + "_trimmed_2.fastq";
		return WORK_DIR.equals(".") ? fn : WORK_DIR + "/" + fn;
	}
	
	/**
	 * get NR read file name
	 * @return  NR read filename
	 */
	public String getNRReadFileName() {
		if(!doNR)
			return getTrimmedReadFileName();
		String fn = !isPaired ? libName + "_NR.fastq" : libName + "_NR_1.fastq";
		return !WORK_DIR.equals(".") ? WORK_DIR + "/" + fn : fn;
	}
	
	/**
	 * get NR mate file name
	 * @return  NR mate filename
	 */
	public String getNRMateFileName() {
		if(!doNR)
			return getTrimmedMateFileName();
		String fn = !isPaired ? "" : libName + "_NR_2.fastq";
		return WORK_DIR.equals(".") ? fn : WORK_DIR + "/" + fn;
	}
	
	/**
	 * get raw alignment bam file name for this library
	 * @return  BAM filename
	 */
	public String getAlignRawFileName() {
		String fn = libName + "_" + refGenome + "_raw.bam"; // always bam output
		return WORK_DIR.equals(".") ? fn : WORK_DIR + "/" + fn;
	}

	/**
	 * get filtered & best-stratum selected alignment bam file name for this library
	 * @return  BAM filename
	 */
	public String getAlignFilteredFileName() {
		String fn = !sortMethod.equals("coordinate") ? libName + "_" + refGenome + "_filtered.bam" : 
			libName + "_" + refGenome + "_filtered_sorted.bam"; // always bam output
		return PROJECT_DIR.equals(".") ? fn : PROJECT_DIR + "/" + fn;
	}
	
	/**
	 * Test whether a given aligner support spliced alignment
	 * @param aligner  aligner name to test
	 * @return true  if this is a known RNA aligner
	 */
	public static boolean isRNAAligner(String aligner) {
		return aligner.equals("tophat2") || aligner.equals("STAR");
	}
	
	/**
	 * Test whether a aligner support controling of the max-hit
	 * @param aligner aligner  aligner name
	 * @return true  if this aligner is known to support max-hit
	 */
	public static boolean supportMaxHit(String aligner) {
		switch(aligner) {
		case "bowtie": case "bowtie2":
			return true;
		case "novoalign": case "seqalto":
			return true;
		case "tophat1": case "tophat2":
			return true;
		case "STAR":
			return true;
		case "bwa": case "bwa-mem": case "bwa-aln":
			return false;
		case "bwa-sw":
			return true;
		default:
			return true;
		}
	}
	// global options
	static String VERSION = EnvConstants.progVer;
	static int MAX_PROC = 6; // maximum processors to use
	static String INIT_MEM = "4G";
	static String MAX_MEM = "16G";
	static String SH_PATH = "/bin/sh";
	static String PROJECT_DIR = ".";
	static String WORK_DIR = ".";
	static final int MIN_UNIQ_INSERT = 15;
	static final float DEFAULT_DNA_ALIGN_RATE = 0.9f;
	static final float DEFAULT_RNA_ALIGN_RATE_RNA_ALIGNER = 0.8f;
	static final float DEFAULT_RNA_ALIGN_RATE_DNA_ALIGNER = 0.5f;
	
	// basic library options
	String libName;
	String readFile;
	int readLen;
	int asciiOffset = 33;
	boolean isPaired;
	String mateFile;
	int mateLen;
	int strandType;
	boolean hasSpliced;
	// trimming options
	boolean doTrim;
	String trimProg = "cutadapt"; // default program
	String adapterSeq3;
	String adapterSeq5;
	float trimMis = 10; // default trimming mismatch %
	int minTrim = 10;
	// NR options
	boolean doNR = true;
	// mapping options
	String aligner;
	int seedLen = 25;
	float seedMis = 4;
	float seedIndel = 0;
	float allMis = 6;
	float allIndel = 0;
	//int minInsert = 15;
	float minAlignRate = -1;
	int maxHit = 10;
	// mate options
	int minFragLen = 0;
	int maxFragLen = 600;
	// user-specified options
	String otherAlignerOpts = "";
	String otherFilterOpts = "";
	// best-stratum options
	int minMapQ = 0;
	int maxBest = 0;
	int maxReport = 0;
	String sortMethod = "none";
	String knownSnpFile;
	// reference genome options
	String refGenome;
	String refIndex;
	// transcriptome options
	String transcriptomeGFF;
	String transcriptomeIndex;
	
	private static List<String> optNames;
	private static Pattern globalPat = Pattern.compile("^### (\\w+)=(.*):");
	private static Pattern localPat = Pattern.compile("^## (\\w+):");
}
