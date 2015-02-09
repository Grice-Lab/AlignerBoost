/**
 * a util class to format SAM/BAM files to costomized tab-delimited cover file
 */
package net.sf.AlignerBoost.utils;
import java.io.*;
import java.util.*;

import net.sf.AlignerBoost.utils.StringUtils;
import static net.sf.AlignerBoost.EnvConstants.*;

/** Format SAM/BAM file to simple tsv cover file
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class ClassifyBED {
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
		BufferedReader bedIn = null;
		BufferedReader gffIn = null;
		BufferedWriter out = null;
		BufferedReader chrIn = null;
		try {
			// open all required files
			bedIn = new BufferedReader(new FileReader(bedInFile));
			chrIn = new BufferedReader(new FileReader(chrLenFile));
			out = new BufferedWriter(new FileWriter(outFile));

			// Read chrLenFile and initialize chrom-index
			if(verbose > 0)
				System.err.println("Initialize chrom-index ...");
			String line = null;
			while((line = chrIn.readLine()) != null) {
				String[] fields = line.split("\t");
				String chr = fields[0];
				int len = Integer.parseInt(fields[1]);
				if(verbose > 0)
					System.err.println("  " + chr + ": " + len);
				chrIdx.put(chr, new int[len + 1]);  // Position 0 is dummy
			}

			if(verbose > 0) {
				// Start the processMonitor to monitor the process
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
							System.err.println("ClassifySAM doesn't support more than " + (Integer.SIZE - 1) + " types");
							gffIn.close();
							out.close();
							chrIn.close();
							return;
						}
						bitMask = initBit;
						typeMask.put(type, bitMask);
						initBit <<= 1; // use a higher bit for the next bit
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
			statusTask.finish();
			statusTask.reset();
			
			// Scan BED file and write BED-DETAIL file
			if(verbose > 0)
				System.err.println("Scanning BED file ...");
			statusTask.setInfo(" regions read");
			boolean hasTrackLine = false;
			boolean isHeader = true;
			long id = 0;
			while((line = bedIn.readLine()) != null) {
				String[] fields = line.split("\t");
				if(fields.length < MIN_N_FIELDS) { // non-record lines
					if(line.startsWith("track")) { // track line exists
						line = line.replaceFirst("name=(?:\\\"[^\"=]+\\\"|\\S+)", "name=\"" + trackName + "\"");
						line = line.replaceFirst("description=(?:\\\"[^\"=]+\\\"|\\S+)", "description=\"" + trackDesc + "\"");
						line = line.replaceFirst("type=\\w+", "type=" + trackType);
						hasTrackLine = true;
					}
					out.write(line + "\n");
					continue;
				}
				if(isHeader && !hasTrackLine) { // no track-line found
					out.write("track name=\"" + trackName + "\" type=" + trackType + " description=\"" + trackDesc + "\"" + "\n");
					isHeader = false;
				}
				
				String chr = fields[0];
				int start = Integer.parseInt(fields[1]) + 1; // start is 0-based
				int end = Integer.parseInt(fields[2]);

				// get annotation bit
				int[] idx = chrIdx.get(chr);
				if(idx == null) {
					if(verbose > 0)
						System.err.println("Ignoring out-of dictionary region: " + line);
					continue;
				}
				if(fix) {
					if(start < 1)
						start = 1;
					if(end > idx.length)
						end = idx.length;
				}
				if(!detail)
					out.write(line + "\t" + (++id) + "\tType=" + unmask(idx, start, end) + "\n");
				else
					out.write(line + "\t" + (++id) + "\tOverlap=" + unmaskDetail(idx, start, end) + "\n");
				if(verbose > 0)
					statusTask.updateStatus(); // Update status
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
				if(bedIn != null)
					bedIn.close();
				if(chrIn != null)
					chrIn.close();
				if(gffIn != null)
					gffIn.close();
				if(out != null)
					out.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void printUsage() {
		System.err.println("java -jar " + progFile + " utils classifyBED " +
				"<-i BED-INFILE> <-g CHR-SIZE-FILE> <-gff GFF-FILE> [-gff GFF-FILE2 -gff ...] <-o BED-DETAIL-OUTFILE> [options]" + newLine +
				"Options:    -name STRING name attribute of the track-line, will override the original value [outfile name]" + newLine +
				"            -desc STRING description attribute of the track-line, will override the original value [-name]" + newLine +
				"            -detail FLAG write detail type overlapping information [not-enable]" + newLine +
				"            -v FLAG  show verbose information" + newLine +
				"            -fix try to fix BED coordinates instead of aborting execution"
				);
	}
	
	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("-i"))
				bedInFile = args[++i];
			else if(args[i].equals("-o"))
				outFile = args[++i];
			else if(args[i].equals("-g"))
				chrLenFile = args[++i];
			else if(args[i].equals("-gff"))
				gffFiles.add(args[++i]);
			else if(args[i].equals("-name"))
				trackName = args[++i];	
			else if(args[i].equals("-desc"))
				trackDesc = args[++i];			
			else if(args[i].equals("-detail"))
				detail = true;	
			else if(args[i].equals("-v"))
				verbose++;
			else if(args[i].equals("-fix"))
				fix = true;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'.");
		}
		// Check required options
		if(bedInFile == null)
			throw new IllegalArgumentException("-i must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-o must be specified");
		if(chrLenFile == null)
			throw new IllegalArgumentException("-g must be specified");
		if(gffFiles.isEmpty())
			throw new IllegalArgumentException("-gff must be specified");
		// set default values
		if(trackName == null)
			trackName = outFile.replaceFirst("\\.bed$", "");
		if(trackDesc == null)
			trackDesc = trackName;
	}

	private static void maskRegion(int[] idx, int start, int end, int bitMask) {
		if(idx == null) // the region is outside of the alignment
			return;
		for(int i = start; i <= end; i++)
			idx[i] |= bitMask;
	}

	private static String unmask(int[] idx, int start, int end, String unkName) {
		assert start >= 0 && start < idx.length && end >= start && end < idx.length;
		int bit = 0;
		for(int i = start; i <= end; i++)
			bit |= idx[i];
		if(bit == 0)
			return unkName;
	    StringBuilder seqType = new StringBuilder();
	    // Test each class
	    for(String type : typeMask.keySet())
	      if((bit & typeMask.get(type)) != 0)
	        seqType.append(seqType.length() == 0 ? type : "," + type);
	    return seqType.toString();
	}
	
	private static String unmask(int[] idx, int start, int end) {
		return unmask(idx, start, end, "intergenic");
	}

	private static String unmaskDetail(int[] idx, int start, int end, String unkName) {
		assert start >= 0 && start < idx.length && end >= start && end < idx.length;
		Map<String, Integer> overlapDetail = new TreeMap<String, Integer>();
		for(String type : typeMask.keySet()) {
			// Initiate
			int mask = typeMask.get(type);
			for(int i = start; i <= end; i++)
				if((idx[i] & mask) != 0) {
					if(!overlapDetail.containsKey(type))
						overlapDetail.put(type, 0);
					overlapDetail.put(type, overlapDetail.get(type) + 1);
				}
		}
	    return overlapDetail.isEmpty() ? unkName + ":" + (end - start + 1) : StringUtils.join(overlapDetail);
	}

	private static String unmaskDetail(int[] idx, int start, int end) {
	    return unmaskDetail(idx, start, end, "intergenic");
	}
	
	private static final int MIN_N_FIELDS = 3; // MIN # of fields for a record line
	private static String chrLenFile;
	private static String bedInFile;
	private static String outFile;
	private static List<String> gffFiles = new ArrayList<String>();
	private static boolean detail;
	private static int verbose;
	private static boolean fix;

	private static Map<String, int[]> chrIdx;
	private static Map<String, Integer> typeMask;

	private static final int statusFreq = 30000;
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static final String trackType = "bedDetail"; // required track line info
	private static String trackName;
	private static String trackDesc;
}
