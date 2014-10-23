package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.*;

import java.io.*;
import java.util.*;

import htsjdk.samtools.*;

/** Filter SAM/BAM single-end (SE) alignments as well as do best-stratum selection to remove too divergent hits
 * @author Qi Zheng
 * @version 1.1
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

		SamReaderFactory readerFac = SamReaderFactory.makeDefault();
		SAMFileWriterFactory writerFac = new SAMFileWriterFactory();
		if(!isSilent)
			readerFac.validationStringency(ValidationStringency.LENIENT); // use LENIENT stringency
		else
			readerFac.validationStringency(ValidationStringency.SILENT); // use SILENT stringency

		SamReader in = readerFac.open(new File(inFile));
		SAMFileWriter out = OUT_IS_SAM ? writerFac.makeSAMWriter(in.getFileHeader(), false, new File(outFile)) : writerFac.makeBAMWriter(in.getFileHeader(), false, new File(outFile));

		// write SAMHeader
		String prevID = null;
		List<SAMRecord> recordList = new ArrayList<SAMRecord>();
		// check each alignment
		SAMRecordIterator results = in.iterator();
		while(results.hasNext()) {
			SAMRecord record = results.next();
			// fix alignment, ignore if failed (unmapped or empty)
			if(!SAMAlignFixer.fixSAMRecord(record, DO_1DP))
				continue;

			String ID = record.getReadName();
			if(!ID.equals(prevID) && prevID != null || !results.hasNext()) { // a non-first new ID meet, or end of alignments
				filterHits(recordList, MIN_INSERT, MAX_SEED_MIS, MAX_SEED_INDEL, MAX_ALL_MIS, MAX_ALL_INDEL);
				// sort the list first with an anonymous class of comparator, with DESCREASING order
				Collections.sort(recordList, Collections.reverseOrder(recordComp));
				int bestScore = !recordList.isEmpty() ? getSAMRecordAlignScore(recordList.get(0)) : 0;
				// remove non-best hits
				removeSecondaryHits(recordList, bestScore, MAX_DIV);
				if(MAX_BEST > 0 && recordList.size() > MAX_BEST) // too much best hits, ignore this read
					recordList.clear();
				// report remaining alignments, up-to MAX_REPORT
				for(int i = 0; i < recordList.size() && i < MAX_REPORT; i++)
					out.addAlignment(recordList.get(i));
				// reset list
				recordList.clear();
			}
			// update
			prevID = ID;
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
	}

	// a nested class for sorting SAMRecord using align score
	static class SAMRecordIdentityComparator implements Comparator<SAMRecord> {
		public int compare(SAMRecord r1, SAMRecord r2) {
			return getSAMRecordAlignScore(r1) - getSAMRecordAlignScore(r2);
		}
	}

	private static void printUsage() { 
		System.err.println("Usage:   java -jar " + progFile + " run filterSE " +
				"<-in SAM|BAM-INPUT> <-out SAM|BAM-OUTPUT> [options]" + newLine +
				"Options:    --min-insert  minimum insert length (excluding adapters) of a read to allow amabiguous alignment, default 15" + newLine +
				"            --seed-len seed length for Burrows-Wheeler algorithm dependent aligners, default 25" + newLine +
				"            --seed-mis %mismatches allowed in seed region, default 4" + newLine +
				"            --all-mis %mismatches allowed in the entire insert region (excluding masked regions and Ns), default 6" + newLine +
				"            --all-indel %in-dels allowed in the entire insert region, default 2" + newLine +
				"            --1DP enable 1-dimentional dymamic programming (1DP) re-aligning, useful for non-local aligners, i.e. bowtie" + newLine +
				"            --match-score match score for 1DP, default 1" + newLine +
				"            --mis-score mismatch score for 1DP, default -2" + newLine +
				"            --gap-open-penalty gap open penalty for 1DP, default 4" + newLine +
				"            --gap-ext-penalty gap extension penalty, default 1" + newLine +
				"            --out-SAM write SAM text output instead of BAM binary output" + newLine +
				"            --silent ignore certain SAM format errors such as empty reads" + newLine +
				"            --max-div max %divergent allowed for best stratum hits comparing to the top hit as for the identity%, default 4" + newLine +
				"            --max-best max allowed best-stratum hits to report for a given read, set to 0 for no limit, default 0" + newLine +
				"            --max-report max report hits for all valid best stratum hits determined by --max-div and --max-best, set to 0 for no limit, default 10, default 6"
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
			else if(args[i].equals("--out-SAM"))
				OUT_IS_SAM = true;
			else if(args[i].equals("--silent"))
				isSilent = true;
			else if(args[i].equals("--max-div"))
				MAX_DIV = Float.parseFloat(args[++i]);
			else if(args[i].equals("--max-best"))
				MAX_BEST = Integer.parseInt(args[++i]);
			else if(args[i].equals("--max-report")) {
				MAX_REPORT = Integer.parseInt(args[++i]);
				if(MAX_REPORT == 0)
					MAX_REPORT = Integer.MAX_VALUE;
			}
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
		if(OUT_IS_SAM && outFile.endsWith(".bam"))
			System.err.println("Warning: output file '" + outFile + "' might not be SAM format");
		if(MAX_DIV < 0)
			throw new IllegalArgumentException("--max-div must be non negative");
		if(MAX_BEST < 0)
			throw new IllegalArgumentException("--max-best must be non negative integer");
		if(MAX_REPORT < 0)
			throw new IllegalArgumentException("--max-report must be non negative integer");
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
	 * @param recordList
	 * @param minInsert
	 * @param maxSeedMis
	 * @param maxSeedIndel
	 * @param maxAllMis
	 * @param maxAllIndel
	 */
	private static int filterHits(List<SAMRecord> recordList, int minInsert,
			double maxSeedMis, double maxSeedIndel, double maxAllMis, double maxAllIndel) {
		int n = recordList.size();
		int removed = 0;
		for(int i = n - 1; i >= 0; i--) {// search backward
			SAMRecord record = recordList.get(i);
			if(!(getSAMRecordIdentity(record) >= MIN_IDENTITY && getSAMRecordInsertLen(record) >= minInsert
					&& getSAMRecordPercentSeedMis(record) <= maxSeedMis && getSAMRecordPercentSeedIndel(record) <= maxSeedIndel
					&& getSAMRecordPercentAllMis(record) <= maxAllMis && getSAMRecordPercentAllIndel(record) <= maxAllIndel)) {
				recordList.remove(i);
				removed++;
			}
		}
		return removed;
	}

	/** Remove secondary hits from sorted List of records to allow only best-stratum hits
	 */
	private static int removeSecondaryHits(List<SAMRecord> recordList, int bestScore, float maxDiv) {
		int n = recordList.size();
		int removed = 0;
		for(int i = n - 1; i >= 0; i--) { // search backward
			//System.err.printf("n:%d removed:%d bestIden:%f iden:%f maxDiv:%f%n", n, removed, bestIden, getSAMRecordIdentity(recordList.get(i)), maxDiv);
			if((float) getSAMRecordAlignScore(recordList.get(i)) / bestScore < 1 - maxDiv / 100f) {
				recordList.remove(i);
				removed++;
			}
			else
				break;
		}
		return removed;
	}

	private static String inFile;
	private static String outFile;
	// filter options
	private static int MIN_INSERT = 15;
	private static double MAX_SEED_MIS = 4; // max % seed mismatch
	private static final double MAX_SEED_INDEL = 0; // seed indel is always not allowed
	private static final float MIN_IDENTITY = 0;
	private static double MAX_ALL_MIS = 6; // max % all mismatch
	private static double MAX_ALL_INDEL = 0; // max % all indel
	private static boolean DO_1DP;
	private static boolean isSilent; // ignore SAM warnings?
	// best stratum options
	private static float MAX_DIV = 1; // max divergent%
	private static int MAX_BEST = 0; // no limits
	private static int MAX_REPORT = 10;
	private static SAMRecordIdentityComparator recordComp = new SAMRecordIdentityComparator();
	// general options
	private static boolean OUT_IS_SAM; // outFile is SAM format?
}
