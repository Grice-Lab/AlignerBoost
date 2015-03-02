/**
 * a util class to format SAM/BAM files to costomized tab-delimited cover file
 */
package net.sf.AlignerBoost.utils;
import java.io.*;
import java.util.*;

import static net.sf.AlignerBoost.EnvConstants.*;
import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;


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
		VCFFileReader vcfIn = new VCFFileReader(new File(vcfInFile), false);
		BufferedReader gffIn = null;
		VariantContextWriterBuilder outBuilder = (new VariantContextWriterBuilder()).setOutputFile(outFile);
		if(buildIndex)
			outBuilder.setOption(Options.INDEX_ON_THE_FLY);
		else
			outBuilder.unsetOption(Options.INDEX_ON_THE_FLY);
		VariantContextWriter out = outBuilder.build();
		BufferedReader chrIn = null;
		try {
			// open all required files
			chrIn = new BufferedReader(new FileReader(chrLenFile));

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
			
			// Scan VCF file and output
			if(verbose > 0)
				System.err.println("Scanning VCF file ...");
			statusTask.setInfo("variants scanned");
			VCFHeader outHeader = new VCFHeader(vcfIn.getFileHeader()); // use a deep copy of the vcfInFile header
			// Add a new Info field
			outHeader.addMetaDataLine(new VCFInfoHeaderLine("GTYPE", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Genetic Type"));
			out.writeHeader(outHeader);
			// process each VariantContext record
			for(VariantContext var : vcfIn) {
				VariantContextBuilder varBuilder = new VariantContextBuilder(var);  // builder derived from old record
				// get GTYPE
				int typeBit = 0;
				int[] idx = chrIdx.get(var.getChr());
				for(int i = var.getStart(); i <= var.getEnd(); i++)
					typeBit |= idx[i];
				String gType = unmask(typeBit);
				varBuilder.attribute("GTYPE", gType); // add the new attribute to the INFO field
				out.add(varBuilder.make());
				if(verbose > 0)
					statusTask.updateStatus();
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
				"<-i VCF-INFILE> <-g CHR-SIZE-FILE> <-gff GFF-FILE> [-gff GFF-FILE2 -gff ...] <-o OUT-FILE> [options]" + newLine +
				"Options:    -v          FLAG  show verbose information" + newLine +
				"            --no-index  FLAG  do not build VCF index on-the-fly"
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
				gffFiles.add(args[++i]);
			else if(args[i].equals("-v"))
				verbose++;
			else if(args[i].equals("--no-index"))
				buildIndex = false;
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

	private static String chrLenFile;
	private static String vcfInFile;
	private static String outFile;
	private static List<String> gffFiles = new ArrayList<String>();
	private static int verbose;
	private static boolean buildIndex = true;

	private static Map<String, int[]> chrIdx;
	private static Map<String, Integer> typeMask;

	private static final int statusFreq = 10000;
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
}
