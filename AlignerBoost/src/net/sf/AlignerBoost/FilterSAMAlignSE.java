package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.*;

import java.io.*;
import java.util.*;

import net.sf.AlignerBoost.utils.ProcessStatusTask;
import net.sf.AlignerBoost.utils.Stats;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.GroupOrder;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.variant.vcf.VCFFileReader;

/** Filter SAM/BAM single-end (SE) alignments as well as do best-stratum selection to remove too divergent hits
 * @author Qi Zheng
 * @version 1.2
 * @since 1.1
 */
public class FilterSAMAlignSE {
	public static void main(String[] args) {
		if(args.length == 0) {
			printUsage();
			return;
		}
		try {
			parseOptions(args);
		}
		catch(IllegalArgumentException e) {
			System.err.println("Error: " + e.getMessage());
			printUsage();
			return;
		}

		if(verbose > 0) {
			// Start the processMonitor
			processMonitor = new Timer();
			// Start the ProcessStatusTask
			statusTask = new ProcessStatusTask();
			// Schedule to show the status every 1 second
			processMonitor.scheduleAtFixedRate(statusTask, 0, statusFreq);
		}
		
		// Read in chrList, if specified
		if(chrFile != null) {
			chrFilter = new HashSet<String>();
			try {
				BufferedReader chrFilterIn = new BufferedReader(new FileReader(chrFile));
				String chr = null;
				while((chr = chrFilterIn.readLine()) != null)
					chrFilter.add(chr);
				chrFilterIn.close();
				if(verbose > 0)
					System.err.println("Only looking at alignments on " + chrFilter.size() + " specified chromosomes");
			}
			catch(IOException e) {
				System.err.println("Error: " + e.getMessage());
				return;
			}
		}
		
		// Read in known SNP file, if specified
		if(knownSnpFile != null) {
			if(verbose > 0)
				System.err.println("Checking known SNPs from user specified VCF file");
			knownVCF = new VCFFileReader(new File(knownSnpFile));
		}
		
		SamReaderFactory readerFac = SamReaderFactory.makeDefault();
		SAMFileWriterFactory writerFac = new SAMFileWriterFactory();
		if(!isSilent)
			readerFac.validationStringency(ValidationStringency.LENIENT); // use LENIENT stringency
		else
			readerFac.validationStringency(ValidationStringency.SILENT); // use SILENT stringency

		SamReader in = readerFac.open(new File(inFile));
		SAMFileHeader inHeader = in.getFileHeader();
		if(inHeader.getGroupOrder() == GroupOrder.reference || inHeader.getSortOrder() == SortOrder.coordinate)
			System.err.println("Warning: Input file '" + inFile + "' might be sorted by coordinate and cannot be correctly processed!");

		SAMFileHeader header = inHeader.clone(); // copy the inFile header as outFile header
		//System.err.println(inFile + " groupOrder: " + in.getFileHeader().getGroupOrder() + " sortOrder: " + in.getFileHeader().getSortOrder());
		// reset the orders
		header.setGroupOrder(groupOrder);
		header.setSortOrder(sortOrder);
		SAMFileWriter out = OUT_IS_SAM ? writerFac.makeSAMWriter(header, false, new File(outFile)) : writerFac.makeBAMWriter(header, false, new File(outFile));

		// write SAMHeader
		String prevID = null;
		SAMRecord prevRecord = null;
		List<SAMRecord> recordList = new ArrayList<SAMRecord>();
		// check each alignment
		SAMRecordIterator results = in.iterator();
		if(verbose > 0) {
			System.err.println("Filtering alignments ...");
			statusTask.reset();
			statusTask.setInfo("alignments processed");
		}
		while(results.hasNext()) {
			SAMRecord record = results.next();
			if(verbose > 0)
				statusTask.updateStatus();
			String ID = record.getReadName();

			// fix read and quality string for this read, if is a secondary hit from multiple hits, used for BWA alignment
			if(ID.equals(prevID) && record.getReadLength() == 0)
				SAMAlignFixer.fixSAMRecordRead(record, prevRecord);
			if(chrFilter != null && !chrFilter.contains(record.getReferenceName())) {
				prevID = ID;
				prevRecord = record;
				continue;
			}
			
			// fix alignment, ignore if failed (unmapped or empty)
			if(!SAMAlignFixer.fixSAMRecord(record, knownVCF, DO_1DP)) {
				prevID = ID;
				prevRecord = record;
				continue;
			}

			if(!ID.equals(prevID) && prevID != null || !results.hasNext()) { // a non-first new ID meet, or end of alignments
				// calculate Bayesian based posterior probabilities
				calcHitPostP(recordList);
				// filter hits
				filterHits(recordList, MIN_INSERT, MAX_SEED_MIS, MAX_SEED_INDEL, MAX_ALL_MIS, MAX_ALL_INDEL, MIN_IDENTITY, MIN_MAPQ);
				// sort the list using the mapQ, using DESCREASING order
				Collections.sort(recordList, Collections.reverseOrder(recordComp));
				if(MAX_BEST > 0 && recordList.size() > MAX_BEST) // too much best hits, ignore this read
					recordList.clear();
				// report remaining alignments, up-to MAX_REPORT
				for(int i = 0; i < recordList.size() && (MAX_REPORT == 0 || i < MAX_REPORT); i++) {
					if(doUpdateBit)
						recordList.get(i).setNotPrimaryAlignmentFlag(i != 0);
					out.addAlignment(recordList.get(i));
				}
				// reset list
				recordList.clear();
			}
			// update
			if(!ID.equals(prevID)) {
				prevID = ID;
				prevRecord = record;
			}
			recordList.add(record);
		}

		// close files
		try {
			in.close();
			out.close();
		}
		catch(IOException e) {
			System.err.println(e.getMessage());
		}
		
		// Terminate the monitor task and monitor
		if(verbose > 0) {
			statusTask.cancel();
			statusTask.finish();
			processMonitor.cancel();
		}
	}

	// a nested class for sorting SAMRecord using align score
	static class SAMRecordMapQComparator implements Comparator<SAMRecord> {
		public int compare(SAMRecord r1, SAMRecord r2) {
			return r1.getMappingQuality() - r2.getMappingQuality();
		}
	}

	private static void printUsage() { 
		System.err.println("Usage:   java -jar " + progFile + " run filterSE " +
				"<-in SAM|BAM-INPUT> <-out SAM|BAM-OUTPUT> [options]" + newLine +
				"Options:    --min-insert INT  minimum insert length (excluding adapters) of a read for unamabiguous alignment [15]" + newLine +
				"            --seed-len INT  seed length for Burrows-Wheeler algorithm dependent aligners [25]" + newLine +
				"            --seed-mis DOUBLE  %mismatches allowed in seed region (by --seed-len) [4]" + newLine +
				"            --all-mis  DOUBLE  %mismatches allowed in the entire insert region (excluding clipped/intron regions) [6]" + newLine +
				"            --all-indel DOUBLE  %in-dels allowed in the entire insert region [2]" + newLine +
				"            -i/--identity DOUBLE  mimimum %identity allowd for the alignment as 100 - (%mismatches+%in-dels) [0]" + newLine +
				"            --1DP FLAG  enable 1-dimentional dymamic programming insert re-assesment, useful for non-local aligners, i.e. bowtie" + newLine +
				"            --match-score INT  match score for 1DP and calculating mapQ [1]" + newLine +
				"            --mis-score INT  mismatch score for 1DP and calculating mapQ [-2]" + newLine +
				"            --gap-open-penalty INT  gap open penalty for 1DP and calculating mapQ [4]" + newLine +
				"            --gap-ext-penalty INT  gap extension penalty for 1DP and calculating mapQ [1]" + newLine +
				"            --clip-penalty INT  additional penalty for soft or hard clipped bases for calculating mapQ [0]" + newLine +
				"            --known-SNP-penalty INT  known SNP penalty for calculating mapQ [1]" + newLine +
				"            --known-INDEL-penalty INT  known IN-DEL penalty for calculating mapQ [2]" + newLine +
				"            --known-MULTISUBSTITUTION-penalty INT  known large/multi-substitution penalty for calculating mapQ [3]" + newLine +
				"            --out-SAM FLAG  write SAM text output instead of BAM binary output" + newLine +
				"            --silent FLAG  ignore certain SAM format errors such as empty reads" + newLine +
				"            --min-mapQ INT  min mapQ calculated with Bayesian method [0]" + newLine +
				"            --max-best INT  max allowed best-stratum hits to report for a given read, 0 for no limit [0]" + newLine +
				"            --max-report INT  max report hits for all valid best stratum hits determined by --min-mapQ and --max-best, 0 for no limit [0]" + newLine +
				"            --no-update-bit FLAG  do not update the secondary alignment bit flag (0x100) after filtering" + newLine +
				"            --best-only FLAG  only report unique best hit, equivelant to --max-best 1 --max-report 1" + newLine +
				"            --best FLAG  report the best hit, ignore any secondary hit, equivelant to --max-best 0 --max-report 1" + newLine +
				"            --sort-method STRING  sorting method for output SAM/BAM file, must be \"none\", \"name\" or \"coordinate\" [none]" + newLine +
				"            --chrom-list FILE  pre-filtering file containing one chromosome name per-line" + newLine +
				"            --known-SNP FILE  known SNP file in vcf/gvcf format (v4.0+), will be used for calculating mapQ" + newLine +
				"            -v FLAG  show verbose information"
				);
	}

	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++)
			if(args[i].equals("-in"))
				inFile = args[++i];
			else if(args[i].equals("-out"))
				outFile = args[++i];
			else if(args[i].equals("--min-insert"))
				MIN_INSERT = Integer.parseInt(args[++i]);
			else if(args[i].equals("--seed-len"))
				SAMAlignFixer.setSEED_LEN(Integer.parseInt(args[++i]));
			else if(args[i].equals("--seed-mis"))
				MAX_SEED_MIS = Double.parseDouble(args[++i]);
			else if(args[i].equals("--all-mis"))
				MAX_ALL_MIS = Double.parseDouble(args[++i]);
			else if(args[i].equals("--all-indel"))
				MAX_ALL_INDEL = Double.parseDouble(args[++i]);
			else if(args[i].equals("-i") || args[i].equals("--identity"))
				MIN_IDENTITY = Double.parseDouble(args[++i]);
			else if(args[i].equals("--1DP"))
				DO_1DP = true;
			else if(args[i].equals("--match-score"))
				SAMAlignFixer.setMATCH_SCORE(Integer.parseInt(args[++i]));
			else if(args[i].equals("--mis-score"))
				SAMAlignFixer.setMIS_SCORE(Integer.parseInt(args[++i]));
			else if(args[i].equals("--gap-open-penalty"))
				SAMAlignFixer.setGAP_OPEN_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--gap-ext-penalty"))
				SAMAlignFixer.setGAP_EXT_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--clip-penalty"))
				SAMAlignFixer.setCLIP_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--known-SNP-penalty"))
				SAMAlignFixer.setKNOWN_SNP_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--known-INDEL-penalty"))
				SAMAlignFixer.setKNOWN_INDEL_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--known-MULTISUBSTITUTION-penalty"))
				SAMAlignFixer.setKNOWN_MULTISUBSTITUTION_PENALT(Integer.parseInt(args[++i]));
			else if(args[i].equals("--out-SAM"))
				OUT_IS_SAM = true;
			else if(args[i].equals("--silent"))
				isSilent = true;
			else if(args[i].equals("--min-mapQ"))
				MIN_MAPQ = Integer.parseInt(args[++i]);
			else if(args[i].equals("--max-best"))
				MAX_BEST = Integer.parseInt(args[++i]);
			else if(args[i].equals("--max-report"))
				MAX_REPORT = Integer.parseInt(args[++i]);
			else if(args[i].equals("--no-update-bit"))
				doUpdateBit = false;
			else if(args[i].equals("--best-only")) {
				MAX_BEST = 1;
				MAX_REPORT = 1;
			}
			else if(args[i].equals("--best")) {
				MAX_BEST = 0;
				MAX_REPORT = 1;
			}
			else if(args[i].equals("--sort-method")) {
				switch(args[++i]) {
				case "none":
					groupOrder = GroupOrder.none;
					sortOrder = SortOrder.unsorted;
					break;
				case "name":
					groupOrder = GroupOrder.query;
					sortOrder = SortOrder.queryname;
					break;
				case "coordinate":
					groupOrder = GroupOrder.reference;
					sortOrder = SortOrder.coordinate;
					break;
				default:
					throw new IllegalArgumentException("--sort-method must be \"none\", \"name\" or \"coordinate\"");
				}
			}
			else if(args[i].equals("--chrom-list"))
				chrFile = args[++i];
			else if(args[i].equals("--known-SNP"))
				knownSnpFile = args[++i];
			else if(args[i].equals("-v"))
				verbose++;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'");
		// Check required options
		if(inFile == null)
			throw new IllegalArgumentException("-in must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-out must be specified");
		// Check other options
		if(MAX_SEED_MIS < 0 || MAX_SEED_MIS > 100)
			throw new IllegalArgumentException("--seed-mis must be between 0 to 100");
		if(MAX_ALL_MIS < 0 || MAX_ALL_MIS > 100)
			throw new IllegalArgumentException("--all-mis must be between 0 to 100");
		if(MAX_ALL_INDEL < 0 || MAX_ALL_INDEL > 100)
			throw new IllegalArgumentException("--all-indel must be between 0 to 100");
		if(!(MIN_IDENTITY >= 0 && MIN_IDENTITY <= 100))
			throw new IllegalArgumentException("-i/--identity must be between 0 to 100");
		else
			MIN_IDENTITY /= 100.0; // use absolute identity
		if(OUT_IS_SAM && outFile.endsWith(".bam"))
			System.err.println("Warning: output file '" + outFile + "' might not be SAM format");
		if(MIN_MAPQ < 0)
			throw new IllegalArgumentException("--min-mapQ must be non negative integer");
		if(MAX_BEST < 0)
			throw new IllegalArgumentException("--max-best must be non negative integer");
		if(MAX_REPORT < 0)
			throw new IllegalArgumentException("--max-report must be non negative integer");
	}

	/** get align length from AlignerBoost internal tag
	 * @return the identity if tag "XA" exists
	 * throws {@RuntimeException} if tag "XA" not exists
	 */
	static int getSAMRecordAlignLen(SAMRecord record) throws RuntimeException {
		return record.getIntegerAttribute("XA");
	}
	
	/** get align insert length from AlignerBoost internal tag
	 * @return the identity if tag "XL" exists
	 * throws {@RuntimeException} if tag "XL" not exists
	 */
	static int getSAMRecordInsertLen(SAMRecord record) throws RuntimeException {
		return record.getIntegerAttribute("XL");
	}

	/** get align identity from AlignerBoost internal tag
	 * @return the identity if tag "XI" exists
	 * throws {@RuntimeException} if tag "XI" not exists
	 */
	static float getSAMRecordIdentity(SAMRecord record) throws RuntimeException {
		return record.getFloatAttribute("XI");
	}

	/** get align %seed mismatch from AlignerBoost internal tag
	 * @return the %seed mismatch if tags "YX" and "YL" exist
	 * throws {@RuntimeException} if tags "YX" and "YL" not exist
	 */
	static float getSAMRecordPercentSeedMis(SAMRecord record) throws RuntimeException {
		return 100f * record.getIntegerAttribute("YX") / record.getIntegerAttribute("YL");
	}

	/** get align %seed indel from AlignerBoost internal tag
	 * @return the %seed indel if tags "YG" and "YL" exist
	 * throws {@RuntimeException} if tags "YG" and "YL" not exist
	 */
	static float getSAMRecordPercentSeedIndel(SAMRecord record) throws RuntimeException {
		return 100f * record.getIntegerAttribute("YG") / record.getIntegerAttribute("YL");
	}

	/** get align %seed mismatch from AlignerBoost internal tag
	 * @return the %seed mismatch if tags "ZX" and "XL" exist
	 * throws {@RuntimeException} if tags "ZX" and "XL" not exist
	 */
	static float getSAMRecordPercentAllMis(SAMRecord record) throws RuntimeException {
		return 100f * record.getIntegerAttribute("ZX") / record.getIntegerAttribute("XL");
	}

	/** get align %seed indel from AlignerBoost internal tag
	 * @return the %seed indel if tags "ZG" and "XL" exist
	 * throws {@RuntimeException} if tags "ZG" and "XL" not exist
	 */
	static float getSAMRecordPercentAllIndel(SAMRecord record) throws RuntimeException {
		return 100f * record.getIntegerAttribute("ZG") / record.getIntegerAttribute("XL");
	}
	
	/**
	 * get internal align score
	 * @param record  SAMRecord to look at
	 * @return  align score
	 */
	static int getSAMRecordAlignScore(SAMRecord record) {
		return record.getIntegerAttribute("XQ");
	}

	/**
	 * Filter hits by removing hits not satisfying the user-specified criteria
	 * @param recordList  an array of ALignRecord
	 * @param minInsert  min insert length
	 * @param maxSeedMis  max %seed-mismatches
	 * @param maxSeedIndel  max %seed-indels
	 * @param maxAllMis  max %all-mismatches
	 * @param maxAllIndel  max %all-indels
	 * @param minIdentity  min identity 
	 * @param minQ 
	 */
	private static int filterHits(List<SAMRecord> recordList, int minInsert,
			double maxSeedMis, double maxSeedIndel, double maxAllMis, double maxAllIndel, double minIdentity, int minQ) {
		int n = recordList.size();
		int removed = 0;
		for(int i = n - 1; i >= 0; i--) { // search backward for maximum performance
			SAMRecord record = recordList.get(i);
			if(!(getSAMRecordInsertLen(record) >= minInsert
					&& getSAMRecordPercentSeedMis(record) <= maxSeedMis && getSAMRecordPercentSeedIndel(record) <= maxSeedIndel
					&& getSAMRecordPercentAllMis(record) <= maxAllMis && getSAMRecordPercentAllIndel(record) <= maxAllIndel
					&& getSAMRecordIdentity(record) >= minIdentity && record.getMappingQuality() >= minQ)) {
				recordList.remove(i);
				removed++;
			}
		}
		return removed;
	}

	/**
	 * Calculate the posterior probability mapQ value (in phred scale) using the Bayesian method
	 * @param recordList
	 * 
	 */
	static double[] calcHitPostP(List<SAMRecord> recordList) {
		if(recordList == null) // return null for null list
			return null;
		if(recordList.isEmpty())
			return new double[0]; // return empty array for empty list
		
		int nHits = recordList.size();
		// get un-normalized posterior probs
		double[] postP = new double[nHits];
		for(int i = 0; i < nHits; i++)
			// get postP as priorP * likelihood, with prior proportional to the alignLength
			postP[i] = getSAMRecordAlignLen(recordList.get(i)) * Math.pow(10.0,  getSAMRecordAlignLikelihood(recordList.get(i)));
		// normalize postP
		FilterSAMAlignSE.normalizePostP(postP);
		// reset the mapQ values
		for(int i = 0; i < nHits; i++) {
			//recordList.get(i).setAttribute("XP", Double.toString(postP[i]));
			// add the "XP" flag showing the mapQ value
			double mapQ = Math.round(Stats.phredP2Q(1 - postP[i]));
			if(Double.isNaN(mapQ) || Double.isInfinite(mapQ)) // is NaN or isInfinite
				recordList.get(i).setMappingQuality(INVALID_MAPQ);
			else {
				if(mapQ > MAX_MAPQ)
					mapQ = MAX_MAPQ;
				recordList.get(i).setMappingQuality((int) mapQ);
				//System.err.println(recordList.get(i).getSAMString());
			}
		}
		return postP;
	}

	/**
	 * get align likelihood
	 * @param record  SAMRecord to look at
	 * @return  log10 likelihood
	 */
	static double getSAMRecordAlignLikelihood(SAMRecord record) {
		return Double.parseDouble(record.getStringAttribute("XH"));
	}
	
/**
	 * Normalize posterior probabilities values
	 * @param postP
	 * @return  normalization constant pi
	 */
	static double normalizePostP(double[] postP) {
		double pi = 0; // normalization constant
		for(double p : postP)
			if(p >= 0)
				pi += p;
		for(int i = 0; i < postP.length; i++)
			postP[i] /= pi;
		return pi;
	}

/*	*//**
	 * get align postP
	 * @param record  SAMRecord to look at
	 * @return  posterior probability of this alignment
	 *//*
	static double getSAMRecordAlignPostP(SAMRecord record) {
		return Double.parseDouble(record.getStringAttribute("XP"));
	}*/
	
	static final int INVALID_MAPQ = 255;
	static final double MAX_MAPQ = 250; // MAX meaniful mapQ value, if not 255
	private static String inFile;
	private static String outFile;
	private static String chrFile;
	private static String knownSnpFile;
	// filter options
	private static int MIN_INSERT = 15;
	private static double MAX_SEED_MIS = 4; // max % seed mismatch
	private static final double MAX_SEED_INDEL = 0; // seed indel is always not allowed
	private static double MAX_ALL_MIS = 6; // max % all mismatch
	private static double MAX_ALL_INDEL = 0; // max % all indel
	private static double MIN_IDENTITY = 0; // min identity
	private static boolean DO_1DP;
	private static boolean isSilent; // ignore SAM warnings?
	// best stratum options
	private static int MIN_MAPQ = 0; // min map postP
	private static int MAX_BEST = 0;
	private static int MAX_REPORT = 0;
	private static boolean doUpdateBit = true;
	private static int verbose; // verbose level
	private static SAMRecordMapQComparator recordComp = new SAMRecordMapQComparator();
	private static Set<String> chrFilter;
	private static VCFFileReader knownVCF;
	// general options
	private static GroupOrder groupOrder = GroupOrder.none;
	private static SortOrder sortOrder = SortOrder.unsorted;
	private static boolean OUT_IS_SAM; // outFile is SAM format?
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static final int statusFreq = 30000;
}
