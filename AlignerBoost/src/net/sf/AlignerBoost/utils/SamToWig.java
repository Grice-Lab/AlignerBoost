package net.sf.AlignerBoost.utils;

import static net.sf.AlignerBoost.EnvConstants.newLine;
import static net.sf.AlignerBoost.EnvConstants.progFile;

import java.io.*;
import java.util.*;

import htsjdk.samtools.*;

/** Format SAM/BAM file to UCSC Wiggle (wig) file fixed step format
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class SamToWig {
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

			// Scan SAM/BAM file
			System.err.println("Scan SAM/BAM file ...");
			SAMRecordIterator results = null;
			Set<Integer> chrSeen = new HashSet<Integer>(); // chromosome seen so far
			if(bedFile == null) // no -R specified
				results = samIn.iterator();
			else {
				bedIn = new BufferedReader(new FileReader(bedFile));
				String line = null;
				while((line = bedIn.readLine()) != null) {
					String[] fields = line.split("\t");
					if(fields.length < 3)
						continue;
					int chrI = samIn.getFileHeader().getSequenceIndex(fields[0]);
					int start = Integer.parseInt(fields[1]);
					int end = Integer.parseInt(fields[2]);
					if(chrI != -1) // this Region is in the aligned chromosomes
						bedRegions.add(new QueryInterval(chrI, start, end));
				}
				System.err.println("Read in " + bedRegions.size() + " regions from BED file");
				bedIn.close();
				QueryInterval[] intervals = new QueryInterval[bedRegions.size()];
				results = samIn.query(bedRegions.toArray(intervals), false);
			}
			
			// Initialize chrom-index
			System.err.println("Initialize chrom-index ...");
			for(SAMSequenceRecord headSeq : samIn.getFileHeader().getSequenceDictionary().getSequences()) {
				String chr = headSeq.getSequenceName();
				if(bedFile != null && !chrSeen.contains(headSeq.getSequenceIndex())) // bed file specified and not in the regions
					continue;
				int len = headSeq.getSequenceLength();
				System.err.println("  " + chr + ": " + len);
				chrIdx.put(chr, new int[len + 1]);  // Position 0 is dummy
			}

			// Start the processMonitor to monitor the process
			processMonitor = new Timer();
			// Start the ProcessStatusTask
			ProcessStatusTask statusTask = new ProcessStatusTask();

			// Schedual to show the status every 1 second
			processMonitor.scheduleAtFixedRate(statusTask, 0, 1000);
			
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
							idx[start + pos]++;
							pos++;
						}
						break;
					case S: // soft clip included by default
						if(countSoft) {
							for(int i = 0; i < cigLen; i++) {
								idx[start + pos]++;
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
			processMonitor.cancel();

			// Output
			System.err.println("Output ...");
			if(includeTrack) // output track line
				out.write(trackHeader + "\n");
			for(Map.Entry<String, int[]> entry : chrIdx.entrySet()) {
				String chr = entry.getKey();
				int[] idx = entry.getValue();
				// Find all non-zero positions
				for(int i = 1; i < idx.length; i++) {  // Position 0 is dummy
					if(idx[i] != 0) {  // Set position
						if(idx[i-1] == 0)  // Dummy position 0 is helpful to find the first occurence
							out.write("fixedStep chrom=" + chr + " start=" + i + " step=1\n");
						int val = idx[i];
						if(normRPM) // do RPM normalization, if specified
							out.write((float) val / totalNum * 1e6f + "\n");
						else
							out.write(val + "\n");
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
		System.err.println("java -jar " + progFile + " utils samToWig " +
				"<-i SAM|BAM-INFILE> <-o OUTFILE> [options]" + newLine +
				"Options:    -s strand to look at, 1 for plus, 2 for minus, 3 for both, default: 3" + newLine +
				"            --norm-rpm normalize the coverage with RPM values of total mapped read number" + newLine +
				"            --count-soft including soft-masked regions as covered region, excluding by default" + newLine +
				"            --no-track do not include the 'track-line' as the first line of the Wiggle file as the UCSC required" + newLine + 
				"            -name the track name used to display in UCSC Genome Browser, default is to use the OUTFILE name" + newLine +
				"            -desc the description of the track used to display in UCSC Genome Browser, default to use the track name" + newLine +
				"            -R genome regions to search provided as a BED file, the -i file must be a sorted BAM file with index pre-built"
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
			else if(args[i].equals("--no-track"))
				includeTrack = false;
			else if(args[i].equals("-name"))
				trackName = args[++i];
			else if(args[i].equals("-desc"))
				trackDesc = args[++i];
			else if(args[i].equals("-count-soft"))
				countSoft = true;
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
		// Set track name and desc
		if(trackName == null)
			trackName = outFile;
		if(trackDesc == null)
			trackDesc = trackName;
		// Finalize track header
		trackHeader += " name=" + trackName + " description=" + trackDesc;
	}

	private static String samInFile;
	private static String outFile;
	private static String bedFile;
	private static int myStrand = 3;
	private static boolean normRPM;
//	private static boolean isLog;
	private static boolean includeTrack = true; // include track line by default
	private static String trackName;
	private static String trackDesc;
	private static String trackHeader = "track type=wiggle_0";
	private static boolean countSoft; // whether to count soft-clipped bases
	private static List<QueryInterval> bedRegions; // bed file regions as the query intervals

	private static long totalNum;
	private static Map<String, int[]> chrIdx;

	private static Timer processMonitor;

}
