/**
 * a util class to format SAM/BAM files to costomized tab-delimited cover file
 */
package net.sf.AlignerBoost.utils;
import java.io.*;
import java.util.*;

import htsjdk.samtools.*;
import static net.sf.AlignerBoost.EnvConstants.*;

/** Format SAM/BAM file to simple tsv cover file
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class QuickSamClassify {
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
		BufferedReader gffIn = null;
		BufferedWriter out = null;
		BufferedReader bedIn = null;
		try {
			// open all required files
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
			if(verbose > 0) {
				processMonitor = new Timer();
				// Start the ProcessStatusTask
				statusTask = new ProcessStatusTask("GFF features read");

				// Schedule to show the status every 1 second
				processMonitor.scheduleAtFixedRate(statusTask, 0, statusFreq);
			}
			
			// Read and index GFF files
			int initBit = 01;
			typeMask = new TreeMap<String, Integer>();
			for(String gffFile : gffFiles) {
				gffIn = new BufferedReader(new FileReader(gffFile));
				String line = null;
				while((line = gffIn.readLine()) != null) {
					if(line.startsWith("#")) // comment line
						continue;
					String[] fields = line.split("\t");
					String chr = fields[0];
					String type = fields[2];
					int start = Integer.parseInt(fields[3]);
					int end = Integer.parseInt(fields[4]);
					// set bitMask
					int bitMask = 0;
					if(!typeMask.containsKey(type)) { // a new type encountered
						if(typeMask.size() >= Integer.SIZE) { // no more free bits available
							System.err.println("QuickSamClassify doesn't support more than " + (Integer.SIZE - 1) + " types");
							gffIn.close();
							out.close();
							return;
						}
						bitMask = initBit;
						typeMask.put(type, bitMask);
						initBit <<= 1; // use a higher bit for the next new type
					}
					else
						bitMask = typeMask.get(type);
					// mask this genome region
					maskRegion(chrIdx.get(chr), start, end, bitMask);
					if(verbose > 0)
						statusTask.updateStatus();
				}
				gffIn.close();
			}
			// reset statusTask
			if(verbose > 0) {
				statusTask.finish();
				statusTask.reset();
			}
			
			// Scan SAM/BAM file and output
			if(verbose > 0) {
				System.err.println("Scanning SAM/BAM file ...");
				statusTask.setInfo("alignments scanned");
			}
			out.write("name\tchrom\tstrand\tstart\tend\ttype" + newLine);
			while(results.hasNext()) {
				SAMRecord record = results.next();
				if(verbose > 0)
					statusTask.updateStatus(); // Update status
				int readLen = record.getReadLength();
				if(record.getReadUnmappedFlag() || record.getReferenceIndex() == -1 || readLen == 0) // non mapped read or 0-length read
					continue;
				int typeBit = 0;
				String id = record.getReadName();
				String chr = record.getReferenceName();
				String strand = record.getReadNegativeStrandFlag() ? "-" : "+";
				int start = record.getAlignmentStart();
				int end = record.getAlignmentEnd();
				int[] idx = chrIdx.get(chr);
				// check each alignment block
				for(AlignmentBlock block : record.getAlignmentBlocks()) {
					int blockStart = block.getReferenceStart();
					int blockLen = block.getLength();
					// pad this region
					for(int j = blockStart; j < blockStart + blockLen; j++)
						typeBit |= idx[j];
				}
				// output
				out.write(id + "\t" + chr + "\t" + strand + "\t" + start + "\t" + end + "\t" + unmask(typeBit) + newLine);
			} // end each record
			out.close();
			// Terminate the monitor task and monitor
			if(verbose > 0) {
				statusTask.cancel();
				statusTask.finish();
				processMonitor.cancel();
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
				if(gffIn != null)
					gffIn.close();
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
		System.err.println("java -jar " + progFile + " utils quickClassify " +
				"<-i SAM|BAM-INFILE> <-gff GFF-FILE> [-gff GFF-FILE2 -gff ...] <-o OUT-FILE> [options]" + newLine +
				"Options:    -R FILE  genome regions to search provided as a BED file; if provided the -i file must be a sorted BAM file with pre-built index" + newLine +
				"            -v FLAG  show verbose information"
				);
	}
	
	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("-i"))
				samInFile = args[++i];
			else if(args[i].equals("-o"))
				outFile = args[++i];
			else if(args[i].equals("-gff"))
				gffFiles.add(args[++i]);
			else if(args[i].equals("-R"))
				bedFile = args[++i];
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
		if(gffFiles.isEmpty())
			throw new IllegalArgumentException("-gff must be specified");
	}

	private static void maskRegion(int[] idx, int start, int end, int bitMask) {
		if(idx == null) // the region is outside of the alignment
			return;
		for(int i = start; i <= end; i++)
			idx[i] |= bitMask;
	}

	private static String unmask(int bit, String unkName) {
		if(bit == 0)
			return unkName;
	    StringBuilder seqType = new StringBuilder();
	    // Test each class
	    for(String type : typeMask.keySet())
	      if((bit & typeMask.get(type)) != 0)
	        seqType.append(seqType.length() == 0 ? type : "," + type);
	    return seqType.toString();
	}
	
	private static String unmask(int bit) {
		return unmask(bit, "intergenic");
	}

	private static String samInFile;
	private static String outFile;
	private static List<String> gffFiles = new ArrayList<String>();
	private static String bedFile;
	private static List<QueryInterval> bedRegions; // bed file regions as the query intervals
	private static int verbose;

	private static Map<String, int[]> chrIdx;
	private static Map<String, Integer> typeMask;

	private static final int statusFreq = 30000;
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
}
