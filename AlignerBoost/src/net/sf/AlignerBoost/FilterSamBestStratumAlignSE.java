package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.progFile;

import java.io.*;
import java.util.*;

import htsjdk.samtools.*;

/** Parse SAM/BAM single-end (SE) alignments for best-stratum final alignments to remove too divergent hits
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class FilterSamBestStratumAlignSE {
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
		List<SAMRecord> recordList = new ArrayList<SAMRecord>();
		// check each alignment
		for(SAMRecord record : in) {
			int readLen = record.getReadLength();
			if(record.getReferenceIndex() == -1 || readLen == 0) // non mapped read or 0-length read
				continue;
			// calculate identity for alignment
			String ID = record.getReadName();
			if(!ID.equals(prevID) && prevID != null) { // a new ID meet and not first one
				// sort the list first with an anonymous class of comparator
				Collections.sort(recordList, recordComp);
				float bestIdentity = calcSAMRecordIdentity(recordList.get(0));
				// remove non-best hits
				removeSecondaryHits(recordList, bestIdentity, MAX_DIV);
				if(MAX_BEST > 0 && recordList.size() > MAX_BEST) // too much best hit, ignore this read
					recordList.clear();
				// report remaining alignments
				for(int i = 0; i < recordList.size() && i < MAX_REPORT; i++)
					out.addAlignment(recordList.get(i));
				// reset list
				recordList.clear();
			}
			// update
			prevID = ID;
			recordList.add(record);
		}
		// process last ID
		Collections.sort(recordList, Collections.reverseOrder(recordComp)); // sort the recordList in reverse order
		float bestIdentity = calcSAMRecordIdentity(recordList.get(0));
		removeSecondaryHits(recordList, bestIdentity, MAX_DIV);
		if(MAX_BEST > 0 && recordList.size() > MAX_BEST) // too much best hit, ignore this read
			recordList.clear();
		// filter secondary alignments
		for(int i = 0; i < recordList.size() && i < MAX_REPORT; i++)
			if(bestIdentity - calcSAMRecordIdentity(recordList.get(i)) <= MAX_DIV)
				out.addAlignment(recordList.get(i));
		// reset list
		recordList.clear();

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

	// a nested class for sorting SAMRecord using identity, ties are broken by the insert length
	static class SAMRecordIdentityComparator implements Comparator<SAMRecord> {
		public int compare(SAMRecord r1, SAMRecord r2) {
			float iden1 = calcSAMRecordIdentity(r1);
			float iden2 = calcSAMRecordIdentity(r2);
			if(iden1 < iden2)
				return -1;
			else if(iden1 > iden2)
				return 1;
			else
				return r1.getIntegerAttribute("XL") - r2.getIntegerAttribute("XL"); // higher identity record rank higher
		}
	}

	private static String inFile;
	private static String outFile;
	private static boolean OUT_IS_SAM; // outFile is SAM format?
	private static float MAX_DIV = 0.04f;
	private static int MAX_BEST = 0; // no limits
	private static int MAX_REPORT = 10;
	private static SAMRecordIdentityComparator recordComp = new SAMRecordIdentityComparator();

	private static void printUsage() { 
		System.err.println("Usage:   java -jar " + progFile + " run stratumSE " +
				"<-in SAM|BAM-INPUT> <-out SAM|BAM-OUTPUT> [options]" +
				"Options:    --out-SAM write SAM text output intead of BAM binary output" +
				"            --max-div max %divergent allowed for best stratum hits comparing to the top hit as for the identity%, default 4" +
				"            --max-best max allowed best-stratum hits to report for a given read, set to 0 for no limit, default 0" +
				"            --max-report max report hits for all valid best stratum hits determined by --max-div and --max-best, set to 0 for no limit, default 10, default 6"
		);
	}

	/** Calculate read align identity using extra tags created in previous filter and fix step
	 * @return true identity if exists all the tags, or throws RuntimeException (non checked)
	 */
	static float calcSAMRecordIdentity(SAMRecord record) throws RuntimeException {
		return record.getFloatAttribute("XI");
	}

	/** Remove secondary hits from sorted List of records to allow only best-stratum hits
	 */
	private static void removeSecondaryHits(List<SAMRecord> recordList, float bestIden, float maxDiv) {
		int n = recordList.size();
		for(int i = n - 1; i >= 0; i--) // search from back
			if(calcSAMRecordIdentity(recordList.get(i)) < bestIden - maxDiv)
				recordList.remove(i);
			else
				return;
	}
}
