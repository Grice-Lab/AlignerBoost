/**
 * a util class to format SAM/BAM files to costomized tab-delimited cover file
 */
package net.sf.AlignerBoost.utils;
import java.io.*;
import java.util.*;

import net.sf.AlignerBoost.*;
import static net.sf.AlignerBoost.EnvConstants.*;

/** Format SAM/BAM file to simple tsv cover file
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class ClassifyVCF {
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
		BufferedReader vcfIn = null;
		BufferedReader gffIn = null;
		BufferedWriter out = null;
		BufferedReader chrIn = null;
		try {
			// open all required files
			vcfIn = new BufferedReader(new FileReader(vcfInFile));
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
							System.err.println("QuickSamClassify doesn't support more than " + (Integer.SIZE - 1) + " types");
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
					maskRegion(chrIdx.get(chr), chr, start, end, bitMask);
					if(verbose > 0)
						statusTask.updateStatus();
				}
				gffIn.close();
			}
			// reset statusTask
			statusTask.finish();
			statusTask.reset();
			
			// Scan VCF file and output
			if(verbose > 0)
				System.err.println("Scanning VCF file ...");
			statusTask.setInfo("variants scanned");
			boolean infoWritten = false;
			String prevLine = "";
			while((line = vcfIn.readLine()) != null) {
				if(line.startsWith("##")) // comment line
					out.write(line + "\n");
				else if(prevLine.startsWith("##INFO") && !line.startsWith("##INFO") || // INFO field end
						!infoWritten && line.startsWith("#CHROM")) { /// no INFO field found, header line found
					out.write("##" + classInfo + "\n");
					out.write(line + "\n");
					infoWritten = true;
				}
				else {
					String[] fields = line.split("\t");
					// contruct a SNP from this record
					String chr = fields[0];
					int loc = Integer.parseInt(fields[1]);
					SNP snp = new SNP(chr, loc, fields[3], fields[4]);
					if(doSimplify)
						snp.simplify();
					// get annotation bit
					int typeBit = 0;
					int[] idx = chrIdx.get(chr);
					for(int i = snp.getLoc(); i < snp.getLoc() + snp.chrLen(); i++)
						typeBit |= idx[i];
					// Update old info and output
					if(!fields[7].equals(".")) // a non-empty INFO field
						fields[7] += ";" + classID + "=" + unmask(typeBit);
					else
						fields[7] = classID + "=" + unmask(typeBit);
					// output
					out.write(StringUtils.join("\t", fields) + "\n");
					if(verbose > 0)
						statusTask.updateStatus(); // Update status
				}
				prevLine = line;
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
				if(vcfIn != null)
					vcfIn.close();
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
		System.err.println("java -jar " + progFile + " utils classifyVCF " +
				"<-i SAM|BAM-INFILE> <-g CHR-SIZE-FILE> <-gff GFF-FILE,[GFF-FILE2,...]> <-o OUT-FILE> [options]" + newLine +
				"Options:    -v FLAG  show verbose information" + newLine +
				"            --no-simplify FLAG  do not try to simpify insertion/deletion/multi-substitution type of variations for accurate positioning"
				);
	}
	
	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("-i"))
				vcfInFile = args[++i];
			else if(args[i].equals("-o"))
				outFile = args[++i];
			else if(args[i].equals("-g"))
				chrLenFile = args[++i];
			else if(args[i].equals("-gff"))
				gffFiles = args[++i].split(",");
			else if(args[i].equals("-v"))
				verbose++;
			else if(args[i].equals("--no-simplify"))
				doSimplify = false;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'.");
		}
		// Check required options
		if(vcfInFile == null)
			throw new IllegalArgumentException("-i must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-o must be specified");
		if(chrLenFile == null)
			throw new IllegalArgumentException("-g must be specified");
		if(gffFiles == null)
			throw new IllegalArgumentException("-gff must be specified");
	}

	private static void maskRegion(int[] idx, String chr, int start, int end, int bitMask) {
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

	private static String chrLenFile;
	private static String vcfInFile;
	private static String outFile;
	private static String[] gffFiles;
	private static int verbose;
	private static boolean doSimplify = true;

	private static Map<String, int[]> chrIdx;
	private static Map<String, Integer> typeMask;

	private static final int statusFreq = 30000;
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static final String classID = "TP";
	private static final String classInfo = "INFO=<ID=TP,Number=.,Type=String,Description=\"Element Type\">";
}
