/**
 * a util class to format SAM/BAM files to costomized tab-delimited cover file
 */
package net.sf.AlignerBoost.utils;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.*;
import static net.sf.AlignerBoost.EnvConstants.*;

/** Format SAM/BAM file to simple tsv cover file
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class SamToCover {
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
						bedRegions.add(interval);
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
				if(verbose > 0)
					System.err.println("  " + chr + ": " + len);
				chrIdx.put(chr, new int[len + 1]);  // Position 0 is dummy
			}

			// Start the processMonitor to monitor the process
			processMonitor = new Timer();
			// Start the ProcessStatusTask
			statusTask = new ProcessStatusTask("alignment(s) processed");

			// Schedule to show the status every 1 second
			processMonitor.scheduleAtFixedRate(statusTask, 0, 1000);
			
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
				Matcher match = nrPat.matcher(record.getReadName()); // whether match interval nrID pattern
				int clone = match.find() ? Integer.parseInt(match.group(1)) : 1;
				
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
				totalNum++;
			} // end each record

			// Terminate the monitor task and monitor
			statusTask.cancel();
			statusTask.finish();
			processMonitor.cancel();

			// Output
			if(verbose > 0)
				System.err.println("Output ...");
			out.write("chrom\tloc\tcover\n");
			for(Map.Entry<String, int[]> entry : chrIdx.entrySet()) {
				String chr = entry.getKey();
				int[] idx = entry.getValue();
				if(bedFile == null) { // bedFile not specified
					// Find all non-zero positions and output
					for(int i = 1; i < idx.length; i++) {  // Position 0 is dummy
						if(idx[i] != 0) { // position is covered
							int val = idx[i];
							if(!normRPM) // raw count
								out.write(chr + "\t" + i + "\t" + val + "\n");
							else // do RPM normalization
								out.write(chr + "\t" + i + "\t" + (float) val / totalNum * 1e6f + "\n");
						}
					}
				}
				else { // bedFile provided, only output restricted regions
					if(!chrSeen.containsKey(chr))
						continue; // ignore this chr if not in bedFile regions
					QueryInterval[] intervals = new QueryInterval[chrSeen.get(chr).size()];
					intervals = chrSeen.get(chr).toArray(intervals); // dump List to array of QueryIntervals
					intervals = QueryInterval.optimizeIntervals(intervals); // optimize intervals
					// output coverage for each interval separately
					for(QueryInterval interval : intervals) {
						for(int i = interval.start; i <= interval.end && i < idx.length; i++) {
							if(idx[i] != 0) { // position is covered
								int val = idx[i];
								if(!normRPM) // raw count
									out.write(chr + "\t" + i + "\t" + val + "\n");
								else // do RPM normalization
									out.write(chr + "\t" + i + "\t" + (float) val / totalNum * 1e6f + "\n");
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
		System.err.println("java -jar " + progFile + " utils samToCover " +
				"<-i SAM|BAM-INFILE> <-o OUTFILE> [options]" + newLine +
				"Options:    -s strand to look at, 1 for plus, 2 for minus and 3 for both, default: 3" + newLine +
				"            --norm-rpm normalize the coverage with RPM values of total mapped read number" + newLine +
				"            --count-soft including soft-masked regions as covered region, excluded by default" + newLine +
				"            -R genome regions to search provided as a BED file; if provided the -i file must be a sorted BAM file with index pre-built" + newLine +
				"            -v show verbose information"
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

	private static String samInFile;
	private static String outFile;
	private static String bedFile;
	private static int myStrand = 3;
	private static boolean normRPM;
	private static boolean countSoft; // whether to count soft-clipped bases
	private static List<QueryInterval> bedRegions; // bed file regions as the query intervals
	private static int verbose;

	private static long totalNum;
	private static Map<String, int[]> chrIdx;

	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static Pattern nrPat = Pattern.compile("^(?:tr|un:nr)\\d+:(\\d+):\\d+");
}
