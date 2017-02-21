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
 * A utility class to format SAM/BAM files to customized tab-delimited cover file
 */
package edu.upenn.egricelab.AlignerBoost.utils;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.*;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.*;

/** Generate base-wise coverage summaries from SAM/BAM file
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class SamToCoverSumm {
	public static void main(String[] args) {
		if(args.length == 0) {
			printUsage();
			return;
		}
		// Parse options
		try {
			parseOptions(args);
		}
		catch(IllegalArgumentException e) {
			System.err.println("Error: " + e.getMessage());
			printUsage();
			return;
		}
		
		/* determine cover breaks */

		chrIdx = new HashMap<String, int[]>();
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader samIn = null;
		BufferedWriter out = null;
		BufferedReader bedIn = null;
		try {
			samIn = factory.open(new File(samInFile));
			out = new BufferedWriter(new FileWriter(outFile));

			SAMRecordIterator results = null;
			Map<String, List<QueryInterval>> chrSeen = new HashMap<String, List<QueryInterval>>(); // chromosomes seen so far
			if(bedFile == null) // no -R specified
				results = samIn.iterator();
			else {
				bedIn = new BufferedReader(new FileReader(bedFile));
				bedRegions = new ArrayList<QueryInterval>();
				String line = null;
				while((line = bedIn.readLine()) != null) {
					String[] fields = line.split("\t");
					if(fields.length < 3) // ignore header lines
						continue;
					String chr = fields[0];
					int chrI = samIn.getFileHeader().getSequenceIndex(chr);
					int start = Integer.parseInt(fields[1]) + 1; // bed start is 0-based
					int end = Integer.parseInt(fields[2]);
					if(chrI != -1) { // this Region is in the aligned chromosomes
						QueryInterval interval = new QueryInterval(chrI, start, end);
						bedRegions.add(interval); // add to global search regions
						if(!chrSeen.containsKey(chr))
							chrSeen.put(chr, new ArrayList<QueryInterval>());
						chrSeen.get(chr).add(interval);
					}
				}
				if(verbose > 0)
					System.err.println("Read in " + bedRegions.size() + " regions from BED file");
				bedIn.close();
				QueryInterval[] intervals = new QueryInterval[bedRegions.size()];
				intervals = bedRegions.toArray(intervals); // dump List to array[]
				intervals = QueryInterval.optimizeIntervals(intervals); // optimize and sort the query intervals
				results = samIn.query(intervals, false);
			}
			
			// Initialize chrom-index
			if(verbose > 0)
				System.err.println("Initialize chrom-index ...");
			for(SAMSequenceRecord headSeq : samIn.getFileHeader().getSequenceDictionary().getSequences()) {
				String chr = headSeq.getSequenceName();
				if(bedFile != null && !chrSeen.containsKey(chr)) // bed file specified and not in the regions
					continue;
				int len = headSeq.getSequenceLength();
				chrIdx.put(chr, new int[len + 1]);  // Position 0 is dummy
				if(bedFile == null) { // no -R specified
					chrSeen.put(chr, new ArrayList<QueryInterval>());
					chrSeen.get(chr).add(new QueryInterval(headSeq.getSequenceIndex(), 1, headSeq.getSequenceLength()));
				}
				if(verbose > 0)
					System.err.println("  " + chr + ": " + len);
			}

			if(verbose > 0) {
				// Start the processMonitor to monitor the process
				processMonitor = new Timer();
				// Start the ProcessStatusTask
				statusTask = new ProcessStatusTask("alignment(s) processed");

				// Schedule to show the status every 1 second
				processMonitor.scheduleAtFixedRate(statusTask, 0, statusFreq);
			}
			
			// Scan SAM/BAM file
			if(verbose > 0)
				System.err.println("Scan SAM/BAM file ...");
			while(results.hasNext()) {
				SAMRecord record = results.next();
				if(verbose > 0)
					statusTask.updateStatus(); // Update status
				int readLen = record.getReadLength();
				if(record.getReferenceIndex() == -1 || readLen == 0) // non mapped read or 0-length read
					continue;
				String chr = record.getReferenceName();
				// check strand
				int strand = record.getReadNegativeStrandFlag() ? 2 : 1;
				if((strand & myStrand) == 0)
					continue;
				if(record.getMappingQuality() < minMapQ)
					continue;
				int clone = 1;
				if(doNR) {
					Matcher match = nrPat.matcher(record.getReadName()); // whether match interval nrID pattern
					clone = match.find() ? Integer.parseInt(match.group(1)) : 1;
				}
				
				int[] idx = chrIdx.get(chr);
				int start = record.getUnclippedStart();
				Cigar cigar = record.getCigar();
				int pos = 0; // relative pos to start
				for(CigarElement cigEle : cigar.getCigarElements()) {
					int cigLen = cigEle.getLength();
					CigarOperator cigOp = cigEle.getOperator();
					switch(cigOp) {
					case M: case EQ: case X: case D:
						for(int i = 0; i < cigLen; i++) {
							idx[start + pos] += clone;
							pos++;
						}
						break;
					case S: // soft clip included by default
						if(countSoft) {
							for(int i = 0; i < cigLen; i++) {
								idx[start + pos] += clone;
								pos++;
							}
						}
						else
							pos += cigLen; // ignore soft-clip
						break;
					case N: case H: // ignored bases
						pos += cigLen;
						break;
						//case I: case P: // not present in reference at all
						//  break; // do nothing
					default: // case I or case P
						break;
					}
				} // end each cigEle
			} // end each record
			if(verbose > 0)
				statusTask.finish();

			/* determine cover max */
			if(verbose > 0) {
				System.err.println("Determine coverage range ...");
				statusTask.setInfo("chromosome scanned");
				statusTask.reset();
			}
			for(Map.Entry<String, int[]> entry : chrIdx.entrySet()) {
				int[] idx = entry.getValue();
				for(int val : idx)
					if(maxCover < val)
						maxCover = val;
				if(verbose > 0)
					statusTask.updateStatus();
			}
			// Terminate the monitor task and monitor
			if(verbose > 0) {
				statusTask.cancel();
				statusTask.finish();
				processMonitor.cancel();
			}
			
			int minBreak = Integer.MAX_VALUE;
			int maxBreak = 0;
			breaks = new ArrayList<Integer>();
			for(String br : breakStr.split(",")) {
				int b = Integer.parseInt(br);
				breaks.add(b);
				if(b < minBreak)
					minBreak = b;
				if(b > maxBreak)
					maxBreak = b;
			}
			if(minBreak > minCover) // additional break need at begin
				breaks.add(0, 0);
			if(maxBreak < maxCover) // additional break need at end
				breaks.add(Integer.MAX_VALUE);
			binCoverSumm = new long[breaks.size() - 1];
			
			// summary coverages
			if(verbose > 0)
				System.err.println("Checking basewise coverages ...");
			for(Map.Entry<String, List<QueryInterval>> entry : chrSeen.entrySet()) {
				String chr = entry.getKey();
				QueryInterval[] intervals = new QueryInterval[entry.getValue().size()]; // array to be dumped
				intervals = chrSeen.get(chr).toArray(intervals);
				if(bedFile != null) // optimization required
					intervals = QueryInterval.optimizeIntervals(intervals);
				int[] idx = chrIdx.get(chr);
				if(idx == null)
					continue; // ignore region without index
				
				for(QueryInterval interval : intervals) {
					for(int i = interval.start; i <= interval.end && i < idx.length; i++ /* always do coverage in single bp */) {
						totalCover++;
						int cover = idx[i];
						if(cover == 0)
							continue;
						int k = whichBin(cover, breaks);
						assert(k < breaks.size() - 1);
						binCoverSumm[k]++;
					}
				}
			}
			/* output coverage summaries */
			if(verbose > 0)
				System.err.println("Output summaries ...");
			out.write("bin_name\tbin_min\tbin_max\tcover_length\n");
			for(int k = 0; k < breaks.size() - 1; k++) {
				int bin_min = breaks.get(k);
				int bin_max = breaks.get(k + 1);
				String min = Integer.toString(bin_min);
				String max = bin_max != Integer.MAX_VALUE ? Integer.toString(bin_max) : INF_STR;
				String name = "(" + min + "," + max + "]";
				out.write(name + "\t" + min + "\t" + max + "\t" + binCoverSumm[k] + "\n");
			}
			if(doTotal)
				out.write("total\t0\tInf\t" + totalCover + "\n");
		}
		catch(IOException e) {
			System.err.println(e.getMessage());
		}
		catch(IllegalArgumentException e) {
			System.err.println(e.getMessage());
		}
		finally {
			try {
				if(samIn != null)
					samIn.close();
				if(out != null)
					out.close();
				if(bedIn != null)
					bedIn.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void printUsage() {
		System.err.println("java -jar " + progFile + " utils samToAbsCover " +
				"<-i SAM|BAM-INFILE> <-o OUTFILE> [options]" + newLine +
				"Options:    -s INT  genome strand(s) to look at, 1: plus, 2: minus, 3: both [" + myStrand + "]" + newLine +
				"            --count-soft FLAG  including soft-masked regions as covered region" + newLine +
				"            --nr FLAG  treat read as non-redundant tags, in which their clone information are embedded" + newLine +
				"            -Q/--min-mapQ  INT minimum mapQ cutoff [" + minMapQ + "]" + newLine +
				"            -R FILE  genome regions to search provided as a BED file; if provided the -i file must be a sorted BAM file with pre-built index" + newLine +
				"            -b/--breaks INT1,...,INTN breaks when reporting the basewise coverage [" + DEFAULT_BREAKS + "]" + newLine +
				"            -n/--no-total FLAG  do not report total coverage at the last line" + newLine +
				"            -v FLAG  show verbose information"
				);
	}
	
	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("-i"))
				samInFile = args[++i];
			else if(args[i].equals("-o"))
				outFile = args[++i];
			else if(args[i].equals("-s"))
				myStrand = Integer.parseInt(args[++i]);
			else if(args[i].equals("-R"))
				bedFile = args[++i];
			else if(args[i].equals("--count-soft"))
				countSoft = true;
			else if(args[i].equals("--nr")) {
				doNR = true;
				nrPat = Pattern.compile("^(?:tr|un:nr)\\d+:(\\d+):\\d+");
			}
			else if(args[i].equals("-Q") || args[i].equals("--min-mapQ"))
				minMapQ = Integer.parseInt(args[++i]);
			else if(args[i].equals("-b") || args[i].equals("--breaks")) // customized breaks provided
				breakStr = args[++i];
			else if(args[i].equals("-n") || args[i].equals("--no-total"))
				doTotal = false;
			else if(args[i].equals("-v"))
				verbose++;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'.");
		}
		// Check required options
		if(samInFile == null)
			throw new IllegalArgumentException("-i must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-o must be specified");
		// Reformat myStrand
		if(!(myStrand >= 1 && myStrand <= 3))
			throw new IllegalArgumentException("Unknown -s option, must be 1, 2 or 3");
	}
	
	/**
	 * get first index in which a given value belongs to
	 * @param x  value to check
	 * @param bins  list of breaks
	 * @return  first index for which x in (bins[k], bins[k+1]], or bins.length - 1
	 */
	public static int whichBin(int x, List<Integer> bins) {
		int k;
		for(k = 0; k < bins.size() - 1; k++)
			if(x > bins.get(k) && x <= bins.get(k + 1))
				break;
		return k;
	}

	public static final String DEFAULT_BREAKS = "0,5,10,20,30,100";

	private static String samInFile;
	private static String outFile;
	private static String bedFile;
	private static int myStrand = 3;
	private static boolean countSoft; // whether to count soft-clipped bases
	private static boolean doNR; // whether treat read as NR-tag
	private static int minMapQ;
	private static List<QueryInterval> bedRegions; // bed file regions as the query intervals
	private static String breakStr = DEFAULT_BREAKS;
	private static List<Integer> breaks; // break points, with ( break[i], break[i+1] ] represent the current bin range
	private static boolean doTotal = true;
	private static int verbose;

	private static int minCover = Integer.MAX_VALUE;
	private static int maxCover;
	private static Map<String, int[]> chrIdx;
	private static long[] binCoverSumm;
	private static long totalCover;

	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static Pattern nrPat;
//	private static Pattern nrPat = Pattern.compile("^(?:tr|un:nr)\\d+:(\\d+):\\d+");

	private static final int statusFreq = 10000;
	public static final String INF_STR = "Inf";
}
