package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.*;

import java.io.*;
import java.util.*;

import htsjdk.samtools.*;

/** Filter SAM/BAM single-end (SE) alignments as well as do best-stratum selection to remove too divergent hits
 * @author Qi Zheng
 * @version 1.2
 * @since 1.1
 */
public class FilterSAMAlignPE {
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
		SAMRecord prevRecord = null;
		List<SAMRecord> alnList = new ArrayList<SAMRecord>();
		// check each alignment
		SAMRecordIterator results = in.iterator();
		while(results.hasNext()) {
			SAMRecord record = results.next();
			String ID = record.getReadName();
			// fix read and quality string for this read, if is a secondary hit from multiple hits, used for BWA alignment
			if(ID.equals(prevID) && record.getReadLength() == 0)
				SAMAlignFixer.fixSAMRecordRead(record, prevRecord);
			
			// fix alignment, ignore if failed (unmapped or empty)
			if(!SAMAlignFixer.fixSAMRecord(record, DO_1DP)) {
				prevID = ID;
				prevRecord = record;
				continue;
			}
			if(!record.getReadPairedFlag()) {
				System.err.println("Error: alignment is not from a paired-end read at\n" + record.getSAMString());
				out.close();
				return;
			}

			if(!ID.equals(prevID) && prevID != null || !results.hasNext()) { // a non-first new ID meet, or end of alignments
				// create alnPEList from filtered alnList
				List<SAMRecordPair> alnPEList = createAlnPEListFromAlnList(alnList);
				// calculate posterior mapQ for each pair
				calcPEHitPostP(alnPEList);
				// filter PEhits
				filterPEHits(alnPEList, MIN_INSERT, MAX_SEED_MIS, MAX_SEED_INDEL, MAX_ALL_MIS, MAX_ALL_INDEL, MIN_MAPQ);
				// sort the list first with an anonymous class of comparator, with DESCREASING order
				Collections.sort(alnPEList, Collections.reverseOrder());
				if(MAX_BEST > 0 && alnList.size() > MAX_BEST) // too much best hits, ignore this read
					alnList.clear();
				// report remaining secondary alignments, up-to MAX_REPORT
				for(int i = 0; i < alnPEList.size() && (MAX_REPORT == 0 || i < MAX_REPORT); i++) {
					if(alnPEList.get(i).fwdRecord != null)
						out.addAlignment(alnPEList.get(i).fwdRecord);
					if(alnPEList.get(i).revRecord != null)
						out.addAlignment(alnPEList.get(i).revRecord);
				}
				// reset list
				alnList.clear();
			}
			// update
			if(!ID.equals(prevID)) {
				prevID = ID;
				prevRecord = record;
			}
			alnList.add(record);
		} // end while
		try {
			in.close();
			out.close();
		}
		catch(IOException e) {
			System.err.println(e.getMessage());
		}
	}

	// a nested class for keeping a pair of SAMRecord for PE alignment
	static class SAMRecordPair implements Comparable<SAMRecordPair> {
		public SAMRecordPair(SAMRecord fwdRecord, SAMRecord revRecord) throws IllegalArgumentException {
			if(fwdRecord == null && revRecord == null)
				throw new IllegalArgumentException("forward and reverse SAMRecord cannot be both null");
			this.fwdRecord = fwdRecord;
			this.revRecord = revRecord;
		}

		/** Get whether this SAMRecordPair is actually paired
		 * @return  true if both forward and reverse record are not null
		 */
		public boolean isPaired() {
			return fwdRecord != null && revRecord != null;
		}
		
		/** Get paired-end identity for a AlignmentPair
		 * @return overall identity of the pair
		 */
		public float getPEIdentity() {
			int PEInsertLen = 0;
			int PEMis = 0;
			int PEIndel = 0;
			if(fwdRecord != null) {
				PEInsertLen += getSAMRecordInsertLen(fwdRecord);
				PEMis += getSAMRecordPercentAllMis(fwdRecord);
				PEIndel += getSAMRecordPercentAllIndel(fwdRecord);
			}
			if(revRecord != null) {
				PEInsertLen += getSAMRecordInsertLen(revRecord);
				PEMis += getSAMRecordPercentAllMis(revRecord);
				PEIndel += getSAMRecordPercentAllIndel(revRecord);
			}
			return 1 - (PEMis + PEIndel) / (float) PEInsertLen;
		}
		
		/** Get paired-end align score for a AlignmentPair
		 * @return sum of the align score of the pair
		 */
		public int getPEAlignScore() {
			int PEScore = 0;
			if(fwdRecord != null)
				PEScore += getSAMRecordAlignScore(fwdRecord);
			if(revRecord != null)
				PEScore += getSAMRecordAlignScore(revRecord);
			return PEScore;
		}

		/** Get PE insert length
		 * @return PE insert length as the inferred insert size if paired, or the insertLen of not paired
		 */
		public int getPEInsertLen() {
			int len = 0;
			if(fwdRecord != null)
				len += getSAMRecordInsertLen(fwdRecord);
			if(revRecord != null)
				len += getSAMRecordInsertLen(revRecord);
			return len;
		}
		
		/** Get PE log-likelihood
		 * @return PE align log-likelihood
		 */
		public double getPEAlignLik() {
			if(isPaired())
				return FilterSAMAlignSE.getSAMRecordAlignLikelihood(fwdRecord) + FilterSAMAlignSE.getSAMRecordAlignLikelihood(revRecord);
			double log10Lik = fwdRecord != null ? FilterSAMAlignSE.getSAMRecordAlignLikelihood(fwdRecord) :
				FilterSAMAlignSE.getSAMRecordAlignLikelihood(revRecord); 
			byte[] qual = fwdRecord != null ? fwdRecord.getBaseQualities() : revRecord.getBaseQualities();
			// treat the missing mate as all SOFT-CLIPPED with same quality
			for(int q : qual)
				log10Lik += q / - SAMAlignFixer.PHRED_SCALE * SAMAlignFixer.CLIP_PENALTY;
			return log10Lik;
		}

		/**
		 * Get PE mapQ, which is the same for both forward and reverse read
		 * @return  mapQ either from fwdRecord or revRecord, which one is not null
		 */
		public int getPEMapQ() {
			return fwdRecord != null ? fwdRecord.getMappingQuality() : revRecord.getMappingQuality(); 
		}
		
		/**
		 * Set mapQ to an AlignRecordPair
		 * @param mapQ  mapQ to be set to both pairs
		 */
		public void setPEMapQ(int mapQ) {
			if(fwdRecord != null)
				fwdRecord.setMappingQuality(mapQ);
			if(revRecord != null)
				revRecord.setMappingQuality(mapQ);
		}
		
		/**
		 * Get SAMString for this pair
		 * @return  SAMString for non-null mate
		 */
		public String getSAMString() {
			StringBuilder sam = new StringBuilder();
			if(fwdRecord != null)
				sam.append(fwdRecord.getSAMString());
			if(revRecord != null)
				sam.append(revRecord.getSAMString());
			return sam.toString();
		}
		
		/** implements the Comparable method
		 * @return  the difference between the mapQ value
		 */
		public int compareTo(SAMRecordPair that) {
			return getPEMapQ() - that.getPEMapQ();
		}

		private SAMRecord fwdRecord;
		private SAMRecord revRecord;

	}
	
	public static List<SAMRecordPair> createAlnPEListFromAlnList(List<SAMRecord> alnList) {
		if(alnList == null)
			return null;
		List<SAMRecordPair> alnPEList = new ArrayList<SAMRecordPair>();
		for(int i = 0; i < alnList.size(); i++) {
			SAMRecord currAln = alnList.get(i);
			SAMRecord nextAln = i + 1 < alnList.size() ? alnList.get(i+1) : null;
			boolean currIsFirst = currAln.getFirstOfPairFlag();
			boolean nextIsFirst = nextAln != null && nextAln.getFirstOfPairFlag();
			int currTLen = alnList.get(i).getInferredInsertSize();
			int nextTLen = i + 1 < alnList.size() ? alnList.get(i+1).getInferredInsertSize() : 0;
			if(currIsFirst && !nextIsFirst // this is forward, next is not null and is reverse
					&& Math.abs(currTLen) > 0 && Math.abs(currTLen) == Math.abs(nextTLen)) { // is a align pair on the same Chromosome
				alnPEList.add(new SAMRecordPair(currAln, nextAln));
				i++; // advance through nextAln
			}
			else 
				if(!noMix)// not a pair, deal with current one only
					alnPEList.add(currIsFirst ? new SAMRecordPair(currAln, null) : new SAMRecordPair(null, currAln));
		}
		return alnPEList;
	}
	
	private static void printUsage() { 
		System.err.println("Usage:   java -jar " + progFile + " run filterPE " +
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
				"            --no-mix suppress unpaired alignments for paired reads, by default unpaired alignments are treated separately" + newLine +
				"            --min-mapQ min mapQ calculated with Bayesian method, default 0 (no limit)" + newLine +
				"            --max-best max allowed best-stratum pairs to report for a given read, default 0 (no limit)" + newLine +
				"            --max-report max report pairs for all valid best stratum pairs determined by --min-mapQ and --max-best, default 0 (no limit)"
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
			else if(args[i].equals("--no-mix"))
				noMix = true;
			else if(args[i].equals("--min-mapQ"))
				MIN_MAPQ = Integer.parseInt(args[++i]);
			else if(args[i].equals("--max-best"))
				MAX_BEST = Integer.parseInt(args[++i]);
			else if(args[i].equals("--max-report"))
				MAX_REPORT = Integer.parseInt(args[++i]);
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
		if(MIN_MAPQ < 0)
			throw new IllegalArgumentException("--max-div must be non negative");
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
	 * Calculate the posterior probability mapQ value (in phred scale) using the Bayesian method
	 * @param recordList
	 * 
	 */
	private static double[] calcPEHitPostP(List<SAMRecordPair> alnPEList) {
		if(alnPEList == null) // return null for null list
			return null;
		if(alnPEList.isEmpty())
			return new double[0]; // return empty array for empty list
		
		int nPairs = alnPEList.size();
		// get un-normalized posterior probs
		double[] postP = new double[nPairs];
		for(int i = 0; i < nPairs; i++)
			// get postP as priorP * likelihood, with prior proportional to the alignLength
			postP[i] = alnPEList.get(i).getPEInsertLen() * Math.pow(10.0,  alnPEList.get(i).getPEAlignLik());
		// normalize postP
		FilterSAMAlignSE.normalizePostP(postP);
		// reset the mapQ values
		for(int i = 0; i < nPairs; i++) {
			//recordList.get(i).setAttribute("XP", Double.toString(postP[i]));
			double mapQ = Math.round(SAMAlignFixer.phredP2Q(1 - postP[i]));
			if(Double.isNaN(mapQ) || Double.isInfinite(mapQ)) // is NaN or isInfinite
				alnPEList.get(i).setPEMapQ(INVALID_MAPQ);
			else
				alnPEList.get(i).setPEMapQ(mapQ > MAX_MAPQ ? MAX_MAPQ : (int) mapQ);
		}
		return postP;
	}

	private static int filterPEHits(List<SAMRecordPair> alnPEList, int minInsert,
			double maxSeedMis, double maxSeedIndel, double maxAllMis, double maxAllIndel, int minQ) {
		int n = alnPEList.size();
		int removed = 0;
		for(int i = n - 1; i >= 0; i--) { // search backward for maximum performance
			SAMRecordPair pair = alnPEList.get(i);
			if(!(  (pair.fwdRecord == null || getSAMRecordInsertLen(pair.fwdRecord) >= minInsert
					&& getSAMRecordPercentSeedMis(pair.fwdRecord) <= maxSeedMis
					&& getSAMRecordPercentSeedIndel(pair.fwdRecord) <= maxSeedIndel
					&& getSAMRecordPercentAllMis(pair.fwdRecord) <= maxAllMis
					&& getSAMRecordPercentAllIndel(pair.fwdRecord) <= maxAllIndel)
				&& (pair.revRecord == null || getSAMRecordInsertLen(pair.revRecord) >= minInsert
					&& getSAMRecordPercentSeedMis(pair.revRecord) <= maxSeedMis
					&& getSAMRecordPercentSeedIndel(pair.revRecord) <= maxSeedIndel
					&& getSAMRecordPercentAllMis(pair.revRecord) <= maxAllMis
					&& getSAMRecordPercentAllIndel(pair.revRecord) <= maxAllIndel) 
				&& pair.getPEMapQ() >= minQ) ) {
				alnPEList.remove(i);
/*				System.err.println("Removing pair:\n" + pair.getSAMString());
				if(pair.fwdRecord != null)
					System.err.printf("fwd: seedMis:%f seedIndel:%f allMis:%f allIndel:%f mapQ:%d%n", 
							getSAMRecordPercentSeedMis(pair.fwdRecord), 
							getSAMRecordPercentSeedIndel(pair.fwdRecord), 
							getSAMRecordPercentAllIndel(pair.fwdRecord), 
							getSAMRecordPercentAllIndel(pair.fwdRecord),
							pair.getPEMapQ());
				if(pair.revRecord != null)
					System.err.printf("rev: seedMis:%f seedIndel:%f allMis:%f allIndel:%f mapQ:%d%n", 
							getSAMRecordPercentSeedMis(pair.revRecord), 
							getSAMRecordPercentSeedIndel(pair.revRecord), 
							getSAMRecordPercentAllIndel(pair.revRecord), 
							getSAMRecordPercentAllIndel(pair.revRecord),
							pair.getPEMapQ());
				System.err.printf("minInsert:%d maxSeedMis:%f maxSeedIndel:%f maxAllMis:%f maxAllIndel:%f minQ:%d%n",
						minInsert, maxSeedMis, maxSeedIndel, maxAllMis, maxAllIndel, minQ);*/
				removed++;
			}
		}
		return removed;
	}

	private static final int INVALID_MAPQ = 255;
	private static final int MAX_MAPQ = 250; // MAX meaniful mapQ value, if not 255
	private static String inFile;
	private static String outFile;
	// filter options
	private static int MIN_INSERT = 15;
	private static double MAX_SEED_MIS = 4; // max % seed mismatch
	private static final double MAX_SEED_INDEL = 0; // seed indel is always not allowed
	private static double MAX_ALL_MIS = 6; // max % all mismatch
	private static double MAX_ALL_INDEL = 0; // max % all indel
	private static boolean DO_1DP;
	private static boolean isSilent; // ignore SAM warnings?
	private static boolean noMix; // do not allow unpaired alignments for paired reads?
	// best stratum options
	private static int MIN_MAPQ = 0; // max divergent
	private static int MAX_BEST = 0; // no limits
	private static int MAX_REPORT = 0;
	// general options
	private static boolean OUT_IS_SAM; // outFile is SAM format?
}
