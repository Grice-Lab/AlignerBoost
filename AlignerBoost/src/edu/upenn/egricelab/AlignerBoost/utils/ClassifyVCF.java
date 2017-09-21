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
 * a util class to format SAM/BAM files to costomized tab-delimited cover file
 */
package edu.upenn.egricelab.AlignerBoost.utils;
import java.io.*;
import java.util.*;

import edu.upenn.egricelab.ucsc.GFF;
import edu.upenn.egricelab.ucsc.GFF3;
import edu.upenn.egricelab.ucsc.GTF;

import static edu.upenn.egricelab.AlignerBoost.EnvConstants.*;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;


/** Format SAM/BAM file to simple tsv cover file
 * @author Qi Zheng
 * @version 1.2
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

		gtypeIdx = new GTypeIndex();
		VCFFileReader vcfIn = new VCFFileReader(new File(vcfInFile), false); // do not require index
		BufferedReader gffIn = null;
		VariantContextWriterBuilder outBuilder = (new VariantContextWriterBuilder()).setOutputFile(outFile);
		if(buildIndex)
			outBuilder.setReferenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(dictFile)).setOption(Options.INDEX_ON_THE_FLY);
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
//				if(verbose > 0)
//					System.err.println("  " + chr + ": " + len);
				gtypeIdx.addChr(chr, len);
			}

			if(verbose > 0) {
				System.err.println("Reading GFF annotation files");
				// Start the processMonitor to monitor the process
				processMonitor = new Timer();
				// Start the ProcessStatusTask
				statusTask = new ProcessStatusTask("GFF features read");

				// Schedule to show the status every 1 second
				processMonitor.scheduleAtFixedRate(statusTask, 0, statusFreq);
			}	
			// Read and index GFF files
			for(String gffFile : gffFiles) {
				/* guess GFF specification */
				int gffSpecs = 0;
				if(gffFile.endsWith(".gtf"))
					gffSpecs = 2;
				else if(gffFile.endsWith(".gff") || gffFile.endsWith(".gff3"))
					gffSpecs = 3;
				else
					throw new IOException("Unrecognized GFF file extension" + gffFile);
				gffIn = new BufferedReader(new FileReader(gffFile));
				while((line = gffIn.readLine()) != null) {
					if(line.startsWith("#")) // comment line
						continue;
					GFF record = gffSpecs == 2 ? new GTF(line) : new GFF3(line);
					String chr = record.getSeqname();
					if(!gtypeIdx.hasChr(chr)) /* this chromosome doesn't exist */
						continue;
					String type = record.getType();
					int start = record.getStart(); /* GFF start is 1-based */
					int end = record.getEnd();   /* GFF end is 1-based */
					if(tagName != null && !tagName.isEmpty() && record.hasAttr(tagName))
						type = record.getAttr(tagName);
					// mask index
					gtypeIdx.maskRegion(chr, start - 1, end, type);
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
				// get GTYPE of this variant region
				Map<String, Integer> typeSum = gtypeIdx.unmaskSum(var.getContig(), var.getStart() - 1, var.getEnd());
				String typeStr;
				if(!showSumm)
					typeStr = typeSum.isEmpty() ? unType : StringUtils.join(",", typeSum.keySet());
				else
					typeStr = typeSum.isEmpty() ? unType + ":" + (var.getEnd() - var.getStart() + 1) : StringUtils.join(",", typeSum);
				varBuilder.attribute("GTYPE", typeStr); // add the new attribute to the INFO field
				
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
				"Options:    -g  FILE                chrom size file with 1st column the chromosome names and 2nd column their sizes, required" + newLine +
				"            -i  FILE                VCF/GVCF input file, required" + newLine +
				"            -gff  FILE              GTF/GFF3 annotation file(s) used for classification, required" + newLine +
				"            -o  FILE                VCF/GVCF output file with added classification information in the INFO field, required" + newLine +
				"            -v  FLAG                show verbose information" + newLine +
				"            --sum  FLAG             show summary of mapped feature types" + newLine +
				"            --unclassified  STRING  name for unclassified alignments [" + DEFAULT_UNCLASSIFIED_GTYPE + "]" + newLine +
				"            --tag  STRING           use value of given tag in the attrubute field (9th) instead of type field (3rd) as the genetic type, if available" + newLine +
				"            --no-index  FLAG        do not build VCF index on-the-fly" + newLine +
				"            -d  FILE                SAM reference dictionary file, required unless --no-index"
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
			else if(args[i].equals("--sum"))
				showSumm = true;
			else if(args[i].equals("--unclassified"))
				unType = args[++i];
			else if(args[i].equals("--tag"))
				tagName = args[++i];
			else if(args[i].equals("-v"))
				verbose++;
			else if(args[i].equals("--no-index"))
				buildIndex = false;
			else if(args[i].equals("-d"))
				dictFile = new File(args[++i]);
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
		if(buildIndex && dictFile == null)
			throw new IllegalArgumentException("-d must be specified unless use --no-index");
	}

	private static final String DEFAULT_UNCLASSIFIED_GTYPE = "intergenic";

	private static String chrLenFile;
	private static String vcfInFile;
	private static String outFile;
	private static List<String> gffFiles = new ArrayList<String>();
	private static int verbose;
	private static boolean buildIndex = true;
	private static File dictFile;
	private static boolean showSumm;
	private static String unType = DEFAULT_UNCLASSIFIED_GTYPE;
	private static String tagName;

	private static GTypeIndex gtypeIdx;

	private static final int statusFreq = 10000;
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
}
