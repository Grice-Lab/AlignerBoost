/**
 * 
 */
package net.sf.AlignerBoost;

import static net.sf.AlignerBoost.EnvConstants.newLine;
import static net.sf.AlignerBoost.EnvConstants.progFile;
import static net.sf.AlignerBoost.NGSExpDesign.MAX_PROC;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * @author Qi Zheng
 *
 */
public class PrepareMapCmd {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// parse options
		if(args.length != 2) {
			printUsage();
			return;
		}
		inFile = args[0];
		outFile = args[1];
		
		BufferedWriter out = null;
		try {
			// Get configs
			configs = NGSExpDesign.createNGSExpDesignsFromFile(inFile);
			shPath = NGSExpDesign.SH_PATH;
			
			out = new BufferedWriter(new FileWriter(outFile));
			out.write("#!" + shPath + newLine + newLine);
			// process each lib conf
			for(NGSExpDesign conf : configs) {
				if(conf.refGenome.equals("NA"))
					continue; // ignore if no ref genome
				// set max_insert for trimmed and untrimmed reads separately
				int maxInsertTR = conf.doTrim ? conf.readLen - conf.minTrim + 1 : conf.readLen;
				int maxInsertUN = conf.readLen;
				int seedNMis = (int) Math.floor(conf.seedLen * conf.seedMis);
				// determine maximum #mismatches
				int maxNMisTR, maxNMisUN;
				if(conf.hasSpliced && conf.aligner.equals("bowtie")) { // has spliced reads and using the non-sw bowtie aligner
					maxNMisTR = (int) Math.floor(conf.allMis * conf.minInsert) + maxInsertUN - conf.minInsert;
					maxNMisUN = (int) Math.floor(conf.allMis * conf.minInsert) + maxInsertUN - conf.minInsert;
				}
				else { // no spliced reads or using SW aligner
					maxNMisTR = (int) Math.floor(conf.allMis * maxInsertTR);
					maxNMisUN = (int) Math.floor(conf.allMis * maxInsertUN);
				}
				// prepare mapping cmd
				String cmd = "";
				String prog = "";
				String outFn = conf.getAlignRawFileName();
				String readIn = conf.getNRReadFileName(); // getNRReadFileName is automatic set by doTrim and doNR
				String mateIn = conf.getNRMateFileName();

				// set average qual score, for bowtie only
				int avgQ = conf.doNR ? FASTA_DEFAULT_Q : FASTQ_DEFAULT_Q;
				
				switch(conf.aligner){
				case "bowtie":
					prog = "bowtie";
					String inType = conf.doNR ? " -f " : " -q ";
					String qual = " --phred" + conf.asciiOffset + "-quals ";
					int e = Math.max(maxNMisTR,  maxNMisUN) * avgQ;
					String hit = " -k " + conf.maxHit;
					if(conf.maxHit == 1)
						hit += " --best";
					if(seedNMis > MAX_BOWTIE_SEED_NMIS)
						seedNMis = MAX_BOWTIE_SEED_NMIS;
					String quiet = conf.doTrim && conf.isPaired ? " --quiet " : " ";
					String inFn = !conf.isPaired ? " " + readIn : " -1 " + readIn + " -2 "+ mateIn;
					cmd = prog + inType + qual + " -n " + seedNMis + " -e " + e + " -l " + conf.seedLen + " -y " + hit +
							" --nomaqround -p " + MAX_PROC + " --sam" + quiet + conf.otherAlignerOpts +
							" " + conf.refIndex + inFn + " | samtools view -S -b -o " + outFn + " -";
					break;
				case "bowtie2":
				    prog = "bowtie2";
				    inType = conf.doNR ? " -f " : " -q ";
				    qual = " --phred" + conf.asciiOffset + " ";
				    String mode = " --local ";
				    hit = conf.maxHit > 1 ? " -k " + conf.maxHit : ""; // use default mode if max_hit == 1
				    if(seedNMis > MAX_BOWTIE2_SEED_NMIS)
				    	seedNMis = MAX_BOWTIE2_SEED_NMIS;
				    if(conf.seedLen > MAX_BOWTIE2_SEED_LEN)
				    	conf.seedLen = MAX_BOWTIE2_SEED_LEN;
				    quiet = conf.doTrim && conf.isPaired ? " --quiet " : " ";
				    inFn = !conf.isPaired ? " -U " + readIn : " -1 " + readIn + " -2 " + mateIn;
				    cmd = prog + inType + mode + qual + " -N " + seedNMis + " -L " + conf.seedLen + hit +
				    		" -p " + MAX_PROC + quiet + conf.otherAlignerOpts +
				    		" -x " + conf.refIndex + inFn + " | samtools view -S -b -o " + outFn + " -";
				    break;
				case "bwa-mem": case "bwa":
				    prog = "bwa mem";
				    int minScore = conf.minInsert * BWA_MEM_MATCH_SCORE - (int) Math.floor(conf.minInsert * conf.allMis) * BWA_MEM_MISMATCH_PENALTY;
				    if(minScore < 0)
				    	minScore = 0;
				    inFn = !conf.isPaired ? " " + readIn : " " + readIn + " " + mateIn;
				    cmd = prog + " -t " + MAX_PROC + " -k " + conf.seedLen + " -T " + minScore + " " + conf.otherAlignerOpts + 
				    		" -a " + conf.refIndex + inFn + " | samtools view -S -b -o " + outFn + " -";
				    break;
				case "bwa-sw":
				    prog = "bwa bwasw";
				    minScore = conf.minInsert - (int) Math.floor(conf.minInsert * conf.allMis * BWA_SW_MISMATCH_PENALTY / BWA_SW_MATCH_SCORE);
				    int zBest = conf.maxHit > BWA_SW_MAX_Z_BEST ? BWA_SW_MAX_Z_BEST : conf.maxHit;
				    inFn = !conf.isPaired ? " " + readIn : " " + readIn + " " + mateIn;
				    cmd = prog + " -t " + MAX_PROC + " -T " + minScore + " -z " + zBest + " " + conf.otherAlignerOpts +
				    		" " + conf.refIndex + inFn + " | samtools view -S -b -o " + outFn + " -";
				    break;
				case "bwa-aln":
				    prog = "bwa aln";
				    float maxEdit = conf.allMis + conf.allIndel;
				    String qualFormat = conf.asciiOffset == 33 ? " " : " -I ";
				    inFn = readIn;
				    if(!conf.isPaired) {
				      String saiOut = conf.libName + "_" + conf.refGenome + ".sai";
				      cmd = prog + " -n " + maxEdit + " -l " + conf.seedLen + " -k " + seedNMis + " -t " + MAX_PROC +
				    		  qualFormat + conf.otherAlignerOpts + " " +
				    		  conf.refIndex + " " + inFn + " > " + saiOut + newLine;
				      // format sai file to bam file
				      cmd += "  bwa samse " + conf.refIndex + " " + saiOut + " " + inFn + " | samtools view -S -b -o " + outFn + " -";
				      break;
				    }
				    else {
				      String saiOut1 = conf.libName + "_" + conf.refGenome + ".1.sai";
				      String saiOut2 = conf.libName + "_" + conf.refGenome + ".2.sai";
				      cmd = prog + " -n " + maxEdit + " -l " + conf.seedLen + " -k " + seedNMis + " -t " + MAX_PROC +
				    		  qualFormat + " " + conf.otherAlignerOpts + " " +
				    		  conf.refIndex + " " + readIn + " > " + saiOut1 + newLine;
				      cmd += prog + " -n " + maxEdit + " -l " + conf.seedLen + " -k " + seedNMis + " -t " + MAX_PROC +
				    		  qualFormat + " " + conf.otherAlignerOpts + " " +
				    		  conf.refIndex + " " + mateIn + " > " + saiOut2 + newLine;
				      // format paired sai files to bam file
				      cmd += "  bwa sampe " + conf.refIndex + " " + saiOut1 + " " + saiOut2 + " " +
				    		  readIn + " " + mateIn + " | samtools view -S -b -o " + outFn + " -";
				      break;
				    }
				case "tophat1": case "tophat2":
				    prog = "tophat2";
				    // fix options
				    if(conf.aligner.equals("tophat1")) {
				    		if(seedNMis > MAX_BOWTIE_SEED_NMIS)
				    			seedNMis = MAX_BOWTIE_SEED_NMIS;
				    }
				    else {
				    	if(seedNMis > MAX_BOWTIE2_SEED_NMIS)
				    		seedNMis = MAX_BOWTIE2_SEED_NMIS;
				    	if(conf.seedLen > MAX_BOWTIE2_SEED_LEN)
				    		conf.seedLen = MAX_BOWTIE2_SEED_LEN;
				    }
				    int readNMis = Math.max(maxNMisTR, maxNMisUN);
				    int readNGap = (int) Math.floor(conf.allIndel * conf.readLen);
				    int readEdit = readNMis + readNGap;
				    int maxIns = readNGap;
				    int maxDel = readNGap;
				    qual = conf.asciiOffset == 33 ? " --solexa-quals " : " --solexa1.3-quals ";
				    String libType;
				    if(conf.strandType == 0)
				    	libType = " fr-unstranded ";
				    else if(conf.strandType == 1)
				    	libType = " fr-firststrand ";
				    else if(conf.strandType == 2)
				      libType = " fr-secondstrand ";
				    else
				    	throw new IllegalArgumentException("Invalid option 'strand_type': " + conf.strandType +
				    			", must be one of 0, 1 or 2");
				    int segLen = conf.doTrim ? (conf.readLen - conf.minTrim + 1) / 2 : conf.readLen / 2;
				    // prevent too short/large segments
				    if(segLen < MIN_TOPHAT_SEG_LEN)
				    	segLen = MIN_TOPHAT_SEG_LEN;
				    if(segLen > conf.seedLen)
				    	segLen = conf.seedLen;
				    int segNMis = (int) Math.floor(conf.seedMis * segLen);
				    String juncSearch = " ";
				    if(!conf.hasSpliced) {
				      System.err.println("Warning: not recommended to use tophat2 for non-spliced read mapping for lib: " +
				    		  conf.libName);
				      juncSearch = " --no-novel-juncs --no-gtf-juncs --no-coverage-search --no-novel-indels ";
				    }
				    inFn = !conf.isPaired ? " " + readIn : " " + readIn + " " + mateIn;
				    String dir = conf.libName + "_" + conf.refGenome + "_" + conf.aligner;
				    if(conf.aligner.equals("tophat1"))
				      cmd = prog + " --bowtie1 -N " + readNMis + " --read-gap-length " + readNGap + " --read-edit-dist " + readEdit +
				      " -g " + conf.maxHit + " -x " + conf.maxHit +
				      " --max-insertion-length " + maxIns + " --max-deletion-length " + maxDel + qual + "--library-type" + libType +
				      "-p " + MAX_PROC + " -G " + conf.transcriptomeGFF + " --transcriptome-index " + conf.transcriptomeIndex +
				      " --segment-length " + segLen + " --segment-mismatches " + segNMis + juncSearch +
				      " --no-sort-bam " + conf.otherAlignerOpts + " -o " + dir + " " + conf.refIndex + inFn + newLine;
				    else
				      cmd = prog + " -N " + readNMis + " --read-gap-length " + readNGap + " --read-edit-dist " + readEdit +
				      " -g " + conf.maxHit + " -x " + conf.maxHit +
				      " --max-insertion-length " + maxIns + " --max-deletion-length " + maxDel + qual + "--library-type" + libType +
				      "-p " + MAX_PROC + " -G " + conf.transcriptomeGFF + " --transcriptome-index " + conf.transcriptomeIndex +
				      " --segment-length " + segLen + " --segment-mismatches " + segNMis +
				      "--b2-N " + seedNMis + "--b2-L " + conf.seedLen + juncSearch +
				      " --no-sort-bam " + conf.otherAlignerOpts + " -o " + dir + " " + conf.refIndex + inFn + newLine;
				    // move and rename tophat2 result out
				    cmd += "mv " + dir + "/accepted_hits.bam " + outFn;
				    break;
				case "STAR":
				    prog = "STAR";
				    float minScoreRate = 1 - conf.allMis - conf.allIndel * STAR_INDEL_PENALTY / STAR_MATCH_SCORE;
				    float minMatchRate = minScoreRate;
				    int divRange = (int) Math.floor(conf.maxDiv * conf.readLen);
				    int maxNMis = Math.max(maxNMisTR, maxNMisUN);
				    float maxMisRate = conf.allMis;
				    inFn = !conf.isPaired ? readIn : readIn + " " + mateIn;
				    cmd = prog + " --genomeDir " + conf.refIndex + " --readFilesIn " + inFn + " --runThreadN " + MAX_PROC +
				    		" --outFilterScoreMin " + minScoreRate + " --outFilterMatchNminOverLread " + minMatchRate +
				    		" --outFilterMultimapScoreRange " + divRange + " --outFilterMultimapNmax " + conf.maxHit +
				    		" --outFilterMismatchNmax " + maxNMis + " --outFilterMismatchNoverLmax " + maxMisRate +
				    		" " + conf.otherAlignerOpts + " --outStd SAM - | samtools view -S -b -o " + outFn + " -";
				    break;
				default:
					throw new IllegalArgumentException("Unknown aligner '" + conf.aligner + "'");
				}
				// print cmd
				if(!(new File(outFn)).exists())
					out.write(cmd + newLine);
				else {
					System.err.println("Alignment output file(s) already exist, won't override");
					cmd = cmd.replace("\n", "\n#");
					out.write("#" + cmd + newLine);
				}
			} // end each config
			// chmod of the outFn
			(new File(outFile)).setExecutable(true);
		}
		catch(IOException e) {
			System.err.println("Error: " + e.getMessage());
		}
		catch(IllegalArgumentException e) {
			System.err.println("Error: " + e.getMessage());
		}
		finally {
			try {
				if(out != null)
					out.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void printUsage() {
		System.err.println("Usage:    java -jar " + progFile + " prepare align <EXPERIMENT-CONFIG-INFILE> <BASH-OUTFILE>");
	}

	private static String shPath;
	private static String inFile;
	private static String outFile;
	private static List<NGSExpDesign> configs;
	
	public static final int MAX_BOWTIE_SEED_NMIS = 3;
	public static final int MAX_BOWTIE2_SEED_NMIS = 1;
	public static final int MAX_BOWTIE2_SEED_LEN = 32;
	public static final int BWA_MEM_MATCH_SCORE = 1;
	public static final int BWA_MEM_MISMATCH_PENALTY = 4;
	public static final int BWA_SW_MATCH_SCORE = 1;
	public static final int BWA_SW_MISMATCH_PENALTY = 3;
	public static final int BWA_SW_MAX_Z_BEST = 10;
	public static final int STAR_MATCH_SCORE = 1;
	public static final int STAR_INDEL_PENALTY = 2;
	private static final int FASTA_DEFAULT_Q = 40;
	private static final int FASTQ_DEFAULT_Q = 20;
	private static final int MIN_TOPHAT_SEG_LEN = 25;
}