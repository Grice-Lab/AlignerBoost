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

import htsjdk.samtools.*;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.*;
import edu.upenn.egricelab.ucsc.*;

/** Classify SAM/BAM files using fast per-bp indices
 * generate new SAM/BAM outputs with added aux tags
 * New Tag introduced to SAM/BAM file
 * XT  Z  GTYPE(s), separated by ','
 * @author Qi Zheng
 * @version 1.2
 * @since 1.1
 */
public class ClassifySAM {
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
		SamReaderFactory inFactory = SamReaderFactory.makeDefault();
		SAMFileWriterFactory outFactory = new SAMFileWriterFactory();
		SamReader samIn = null;
		BufferedReader gffIn = null;
		SAMFileWriter samOut = null;
		BufferedReader bedIn = null;
		try {
			// open input
			samIn = inFactory.open(new File(samInFile));
			// clone and modify the header
			SAMFileHeader header = samIn.getFileHeader().clone(); // copy the inFile header as outFile header
			// Add new programHeader
			SAMProgramRecord progRec = new SAMProgramRecord(progName + " utils classifySAM");
			progRec.setProgramName(progName + " utils classifySAM");
			progRec.setProgramVersion(progVer);
			progRec.setCommandLine(StringUtils.join(" ", args));
			header.addProgramRecord(progRec);
			
			samOut = outFactory.makeSAMOrBAMWriter(header, true, new File(outFile));

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
				if(verbose > 0)
					System.err.println("Optimizing customized searching intervals");
				QueryInterval[] intervals = new QueryInterval[bedRegions.size()];
				intervals = bedRegions.toArray(intervals); // dump List to array[]
				intervals = QueryInterval.optimizeIntervals(intervals); // optimize and sort the query intervals
				results = samIn.query(intervals, false);
			}
			
			// Initialize gtype-index
			if(verbose > 0)
				System.err.println("Initialize chrom-index ...");
			for(SAMSequenceRecord headSeq : samIn.getFileHeader().getSequenceDictionary().getSequences()) {
				String chr = headSeq.getSequenceName();
				if(bedFile != null && !chrSeen.containsKey(chr)) // bed file specified and not in the regions
					continue;
				int len = headSeq.getSequenceLength();
//				if(verbose > 0)
//					System.err.println("  " + chr + ": " + len);
				gtypeIdx.addChr(chr, len);
			}

			// Start the processMonitor to monitor the process
			if(verbose > 0) {
				System.err.println("Reading GFF annotation files");
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
				
				String line = null;
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
			if(verbose > 0) {
				statusTask.finish();
				statusTask.reset();
			}
			
			// Scan SAM/BAM file and output
			if(verbose > 0) {
				System.err.println("Scanning SAM/BAM file ...");
				statusTask.setInfo("alignments scanned");
			}
			while(results.hasNext()) {
				SAMRecord record = results.next();
				if(verbose > 0)
					statusTask.updateStatus(); // Update status
				int readLen = record.getReadLength();
				if(record.getReadUnmappedFlag() || record.getReferenceIndex() == -1 || readLen == 0) // non mapped read or 0-length read
					continue;
				if(record.getMappingQuality() < minMapQ)
					continue;
				String chr = record.getReferenceName();
				Map<String, Integer> typeSum = new HashMap<String, Integer>();
				int alignLen = 0;
				// check each alignment block
				for(AlignmentBlock block : record.getAlignmentBlocks()) {
					int blockStart = block.getReferenceStart();
					int blockLen = block.getLength(); /* SAM start is 1-based */
					alignLen += blockLen;
					// add into typeSum
					for(Map.Entry<String, Integer> pair : gtypeIdx.unmaskSum(chr, blockStart - 1, blockStart + blockLen).entrySet())
						typeSum.put(pair.getKey(), typeSum.getOrDefault(pair.getKey(), 0) + pair.getValue());
				}
				// output
				String typeStr;
				if(!showSumm)
					typeStr = typeSum.isEmpty() ? unType : StringUtils.join(",", typeSum.keySet());
				else
					typeStr = typeSum.isEmpty() ? unType + ":" + alignLen : StringUtils.join(",", typeSum);
				record.setAttribute(GTYPE_TAG, typeStr);
				samOut.addAlignment(record);
			} // end each record
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
				if(samOut != null)
					samOut.close();
				if(bedIn != null)
					bedIn.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void printUsage() {
		System.err.println("java -jar " + progFile + " utils ClassifySAM " +
				"<-i SAM|BAM-INFILE> <-gff GFF-FILE> [-gff GFF-FILE2 -gff ...] <-o OUT-FILE> [options]" + newLine +
				"Options:    -i  FILE                SAM/BAM input file, required" + newLine +
				"            -gff  FILE              GTF/GFF3 annotation file(s) used for classification, required" + newLine +
				"            -o  FILE                SAM/BAM output file, required" + newLine +
				"            -R  FILE  genome        regions to search provided as a BED file; if provided the -i file must be a sorted BAM file with pre-built index" + newLine +
				"            -Q/--min-mapQ  INT      minimum mapQ cutoff" + newLine +
				"            --sum  FLAG             show summary of mapped feature types" + newLine +
				"            --unclassified  STRING  name for unclassified alignments [" + DEFAULT_UNCLASSIFIED_GTYPE + "]" + newLine +
				"            --tag  STRING           use value of given tag in the attrubute field (9th) instead of type field (3rd) as the genetic type, if available" + newLine +
				"            -v  FLAG                show verbose information"
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
			else if(args[i].equals("-Q") || args[i].equals("--min-mapQ"))
				minMapQ = Integer.parseInt(args[++i]);
			else if(args[i].equals("--sum"))
				showSumm = true;
			else if(args[i].equals("--unclassified"))
				unType = args[++i];
			else if(args[i].equals("--tag"))
				tagName = args[++i];
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

	private static final String DEFAULT_UNCLASSIFIED_GTYPE = "intergenic";
	private static final String GTYPE_TAG = "XT";
	
	private static String samInFile;
	private static String outFile;
	private static List<String> gffFiles = new ArrayList<String>();
	private static String bedFile;
	private static List<QueryInterval> bedRegions; // bed file regions as the query intervals
	private static int minMapQ;
	private static boolean showSumm;
	private static String unType = DEFAULT_UNCLASSIFIED_GTYPE;
	private static String tagName;
	private static int verbose;

	private static GTypeIndex gtypeIdx;

	private static final int statusFreq = 10000;
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
}
