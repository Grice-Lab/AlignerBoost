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

/** Format SAM/BAM file to simple tsv cover file
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class SamToAbsCover {
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

		chrIdx = new HashMap<String, int[]>();
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader samIn = null;
		BufferedWriter out = null;
		BufferedReader bedIn = null;
		try {
			samIn = factory.open(new File(samInFile));
			out = new BufferedWriter(new FileWriter(outFile));

			// get total alignments, if -norm is set
			if(normRPM) {
				if(verbose > 0)
					System.err.print("Determining total number of alignments ...");
				SAMRecordIterator allResults = samIn.iterator();
				while(allResults.hasNext()) {
					totalNum++;
					allResults.next();
				}
				allResults.close();
				if(verbose > 0)
					System.err.println(totalNum);
			}
			
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

			// Start the processMonitor to monitor the process
			processMonitor = new Timer();
			// Start the ProcessStatusTask
			statusTask = new ProcessStatusTask("alignment(s) processed");

			// Schedule to show the status every 1 second
			processMonitor.scheduleAtFixedRate(statusTask, 0, statusFreq);
			
			// Scan SAM/BAM file
			if(verbose > 0)
				System.err.println("Scan SAM/BAM file ...");
			while(results.hasNext()) {
				SAMRecord record = results.next();
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

			// Terminate the monitor task and monitor
			statusTask.cancel();
			statusTask.finish();
			processMonitor.cancel();

			// Output
			if(verbose > 0)
				System.err.println("Output ...");
			out.write("chrom\tstart\tend\tcover\n");
			for(Map.Entry<String, List<QueryInterval>> entry : chrSeen.entrySet()) {
				String chr = entry.getKey();
				QueryInterval[] intervals = new QueryInterval[entry.getValue().size()]; // array to be dumped
				intervals = chrSeen.get(chr).toArray(intervals);
				if(bedFile != null) // optimization required
					intervals = QueryInterval.optimizeIntervals(intervals);
				int[] idx = chrIdx.get(chr);
				if(idx == null)
					continue; // ignore region without index
				
				if(step == 1) { // no mean required
					for(QueryInterval interval : intervals) {
						for(int i = interval.start; i <= interval.end && i < idx.length; i++) {
							int val = idx[i];
							if(val >= minCover) {
								out.write(chr + "\t" + i + "\t" + (i + 1) + "\t" +
										(!normRPM ? val : (float) (1e6 * val / totalNum)) + "\n");
							}
						}
					}
				}
				else {
					for(QueryInterval interval : intervals) {
						for(int start = interval.start; start <= interval.end && start < idx.length; start += step) {
							int end = Math.min(start + step, idx.length);
							double val = Stats.mean(idx, start, end);
							if(val >= minCover) {
								if(normRPM)
									val /= totalNum / 1e6;
								out.write(chr + "\t" + start + "\t" + (end - 1) + "\t" + (float) val + "\n");
							}
						}
					}
				}
			}
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
				"Options:    -i  FILE                 SAM/BAM input, required" + newLine +
				"            -o  FILE                 TSV coverage output, required" + newLine +
				"            -s  INT                  genome strand(s) to look at, 1: plus, 2: minus, 3: both [" + myStrand + "]" + newLine +
				"            --norm-rpm  FLAG         normalize the coverage to RPM by total mapped read number" + newLine +
				"            --count-soft  FLAG       including soft-masked regions as covered region" + newLine +
				"            --nr  FLAG               treat read as non-redundant tags, in which their clone information are embedded" + newLine +
				"            -Q/--min-mapQ  INT       minimum mapQ cutoff [" + minMapQ + "]" + newLine +
				"            -R  FILE                 genome regions to search provided as a BED file; if provided the -i file must be a sorted BAM file with pre-built index" + newLine +
				"            -step  INT               step width for calculating the coverage or average coverages [" + step + "]" + newLine +
				"            --min-cover  INT         mimimum cover value to report [" + minCover + "]" + newLine +
//				"            -k/--keep-uncover  FLAG  keep 0-covered regions" + newLine +
				"            -v  FLAG                 show verbose information"
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
			else if(args[i].equals("--norm-rpm"))
				normRPM = true;
			else if(args[i].equals("--count-soft"))
				countSoft = true;
			else if(args[i].equals("--nr")) {
				doNR = true;
				nrPat = Pattern.compile("^(?:tr|un:nr)\\d+:(\\d+):\\d+");
			}
			else if(args[i].equals("-Q") || args[i].equals("--min-mapQ"))
				minMapQ = Integer.parseInt(args[++i]);
			else if(args[i].equals("-step"))
				step = Integer.parseInt(args[++i]);
			else if(args[i].equals("--min-cover"))
				minCover = Integer.parseInt(args[++i]);
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
		if(minCover >= 0)
			throw new IllegalArgumentException("--min-cover must be non-negative");
	}

	private static String samInFile;
	private static String outFile;
	private static String bedFile;
	private static int myStrand = 3;
	private static boolean normRPM;
	private static boolean countSoft; // whether to count soft-clipped bases
	private static boolean doNR; // whether treat read as NR-tag
	private static int minMapQ;
	private static List<QueryInterval> bedRegions; // bed file regions as the query intervals
	private static int step = 1;
	private static int minCover = 1;
//	private static boolean keep0;
	private static int verbose;

	private static long totalNum;
	private static Map<String, int[]> chrIdx;

	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static Pattern nrPat;
	private static final int statusFreq = 10000;
}
