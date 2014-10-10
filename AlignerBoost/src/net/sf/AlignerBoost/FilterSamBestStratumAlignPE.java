package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.progFile;

import java.io.*;
import java.util.*;

import htsjdk.samtools.*;

/** Parse SAM/BAM paired-end (PE) alignments for best-stratum final alignments to remove too devergent hits
 * Picard Java API is required in classpath
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class FilterSamBestStratumAlignPE {
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

		SamReader in = readerFac.open(new File(inFile));
		SAMFileWriter out = OUT_IS_SAM ? writerFac.makeSAMWriter(in.getFileHeader(), false, new File(outFile)) : writerFac.makeBAMWriter(in.getFileHeader(), false, new File(outFile));

		// write SAMHeader
		String prevID = null;
		List<SAMRecord> alnList = new ArrayList<SAMRecord>();
		// check each alignment
		for(SAMRecord record : in) {
			int readLen = record.getReadLength();
			if(record.getReferenceIndex() == -1 || readLen == 0) // non mapped read or 0-length read
				continue;
			if(!record.getReadPairedFlag()) {
				System.err.println("Error! alignment is not from a paired-end read at\n" + record.getSAMString());
				System.exit(-1);
			}

			// calculate identity for alignment
			String ID = record.getReadName();
			if(!ID.equals(prevID) && prevID != null) { // a new ID meet and not first one
				List<SAMRecordPair> alnPEList = createAlnPEListFromAlnList(alnList); // create alnPEList from alnList
				// sort the PEList
				Collections.sort(alnPEList, Collections.reverseOrder()); // sort in reverse order
				float bestIdentity = alnPEList.get(0).getPEIdentity();
				// remove non-best hits
				removeSecondaryPEHits(alnPEList, bestIdentity, MAX_DIV);
				if(MAX_BEST > 0 && alnPEList.size() > MAX_BEST) // too much best PE-hit, discard all alignments for this read
					alnPEList.clear();
				// report remaining secondary alignments
				for(int i = 0; i < alnPEList.size() && i < MAX_REPORT; i++) {
					if(alnPEList.get(i).fwdRecord != null)
						out.addAlignment(alnPEList.get(i).fwdRecord);
					if(alnPEList.get(i).revRecord != null)
						out.addAlignment(alnPEList.get(i).revRecord);
				}
				// reset alnList
				alnList.clear();
			}

			// update
			prevID = ID;
			alnList.add(record);
		} // end of while

		// process last ID
		List<SAMRecordPair> alnPEList = createAlnPEListFromAlnList(alnList); // create alnPEList from alnList
		// sort the PEList
		Collections.sort(alnPEList, Collections.reverseOrder()); // sort in reverse order
		float bestIdentity = alnPEList.get(0).getPEIdentity();
		// remove non-best hits
		removeSecondaryPEHits(alnPEList, bestIdentity, MAX_DIV);
		if(MAX_BEST > 0 && alnPEList.size() > MAX_BEST) // too much best PE-hit, discard all alignments for this read
			alnPEList.clear();
		// filter secondary alignments
		for(int i = 0; i < alnPEList.size() && i < MAX_REPORT; i++) {
			if(bestIdentity - alnPEList.get(i).getPEIdentity() <= MAX_DIV) {
				if(alnPEList.get(i).fwdRecord != null)
					out.addAlignment(alnPEList.get(i).fwdRecord);
				if(alnPEList.get(i).revRecord != null)
					out.addAlignment(alnPEList.get(i).revRecord);
			}
		}
		// reset alnList
		alnList.clear();

		// close files
		try {
			in.close();
			out.close();
		}
		catch(IOException e) {
			System.err.println(e.getMessage());
		}
	}

	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++)
			if(args[i].equals("-in"))
				inFile = args[++i];
			else if(args[i].equals("-out"))
				outFile = args[++i];
			else if(args[i].equals("--out-SAM"))
				OUT_IS_SAM = true;
			else if(args[i].equals("--max-div"))
				MAX_DIV = Float.parseFloat(args[++i]) / 100;
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
		if(OUT_IS_SAM && outFile.endsWith(".bam"))
			System.err.println("Warning: output file '" + outFile + "' might not be SAM format");
		if(MAX_DIV < 0)
			throw new IllegalArgumentException("--max-div must be non negative");
		if(MAX_BEST < 0)
			throw new IllegalArgumentException("--max-best must be non negative integer");
		if(MAX_REPORT < 0)
			throw new IllegalArgumentException("--max-report must be non negative integer");
	}

	// a nested class for keeping a pair of SAMRecord for PE alignment
	static class SAMRecordPair implements Comparable<SAMRecordPair> {
		public SAMRecordPair(SAMRecord fwdRecord, SAMRecord revRecord) throws IllegalArgumentException {
			if(fwdRecord == null && revRecord == null)
				throw new IllegalArgumentException("forward and reverse SAMRecord cannot be both null");
			this.fwdRecord = fwdRecord;
			this.revRecord = revRecord;
		}

		/** Get paired-end identity for a AlignmentPair
		 * @return overal identity of the pair
		 */
		public float getPEIdentity() throws RuntimeException {
			int PEInsertLen = 0;
			int PEMis = 0;
			int PEIndel = 0;
			if(fwdRecord != null) {
				PEInsertLen += fwdRecord.getIntegerAttribute("XL");
				PEMis += fwdRecord.getIntegerAttribute("XA");
				PEIndel += fwdRecord.getIntegerAttribute("XG");
			}
			if(revRecord != null) {
				PEInsertLen += revRecord.getIntegerAttribute("XL");
				PEMis += revRecord.getIntegerAttribute("XA");
				PEIndel += revRecord.getIntegerAttribute("XG");
			}
			return 1 - (PEMis + PEIndel) / (float) PEInsertLen;
		}

		/** Get PE insert length
		 * @return PE insert length as the sum of the pairs
		 */
		public int getPEInsertLen() throws RuntimeException {
			int PEInsertLen = 0;
			if(fwdRecord != null)
				PEInsertLen += fwdRecord.getIntegerAttribute("XL");
			if(revRecord != null)
				PEInsertLen += revRecord.getIntegerAttribute("XL");
			return PEInsertLen;
		}

		/** implements the Comparable method
		 * return -1 if PE identity is lower, 1 if higher, ties are broken by PEInsertLen
		 */
		public int compareTo(SAMRecordPair that) {
			float iden = getPEIdentity();
			float idenThat = that.getPEIdentity();
			if(iden < idenThat)
				return -1;
			else if(iden > idenThat)
				return 1;
			else
				return getPEInsertLen() - that.getPEInsertLen();
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
			if(currIsFirst && nextAln != null && !nextIsFirst && Math.abs(currTLen) == Math.abs(nextTLen)) { // is a align pair
				alnPEList.add(new SAMRecordPair(currAln, nextAln));
				i++;
			}
			else // not a pair, deal with current one only
				alnPEList.add(currIsFirst ? new SAMRecordPair(currAln, null) : new SAMRecordPair(null, currAln));
		}
		return alnPEList;
	}

	private static String inFile;
	private static String outFile;
	private static boolean OUT_IS_SAM; // outFile is SAM format?
	private static float MAX_DIV = 0.04f;
	private static int MAX_BEST = 0; // no limits
	private static int MAX_REPORT = 10;

	private static void printUsage() { 
		System.err.println("Usage:   java -jar " + progFile + " run stratumPE " +
				"<-in SAM|BAM-INPUT> <-out SAM|BAM-OUTPUT> [options]" +
				"Options:    --out-SAM write SAM text output intead of BAM binary output" +
				"            --max-div max %divergent allowed for best stratum hits comparing to the top hit as for the identity%, default 4" +
				"            --max-best max allowed best-stratum hits to report for a given read, set to 0 for no limit, default 0" +
				"            --max-report max report hits for all valid best stratum hits determined by --max-div and --max-best, set to 0 for no limit, default 10, default 6"
		);
	}
	
	/** Calculate read align identity using extra tags created in previous filter and fix step
	 * @return true identity if exists all the tags, or throws RuntimeException (none is checked)
	 */
	static float calcSAMRecordIdentity(SAMRecord record) throws RuntimeException {
		return record.getFloatAttribute("XI");
	}

	/** Remove secondary hits from sorted List of records to allow only best-stratum hits
	 */
	private static void removeSecondaryPEHits(List<SAMRecordPair> alnPEList, float bestIden, float maxDiv) {
		int n = alnPEList.size();
		for(int i = n - 1; i >= 0; i--) // search from back
			if(alnPEList.get(i).getPEIdentity() < bestIden - maxDiv)
				alnPEList.remove(i);
			else
				return;
	}
}
