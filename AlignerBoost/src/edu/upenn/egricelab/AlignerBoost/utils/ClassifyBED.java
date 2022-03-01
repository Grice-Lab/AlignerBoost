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

package edu.upenn.egricelab.AlignerBoost.utils;
import java.io.*;
import java.util.*;

import edu.upenn.egricelab.AlignerBoost.utils.StringUtils;
import edu.upenn.egricelab.ucsc.GFF;
import edu.upenn.egricelab.ucsc.GFF3;
import edu.upenn.egricelab.ucsc.GTF;

import static edu.upenn.egricelab.AlignerBoost.EnvConstants.*;

/** Classify UCSC BED files using fast per-bp indices
 * generate new BED outputs with additional 'GTYPE' field at the last field of each line
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

		gtypeIdx = new GTypeIndex();
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
//				if(verbose > 0)
//					System.err.println("  " + chr + ": " + len);
				gtypeIdx.addChr(chr, len);
			}
			chrIn.close();

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
				else {
					out.close();
					throw new IOException("Unrecognized GFF file extension" + gffFile);
				}
				gffIn = new BufferedReader(new FileReader(gffFile));
				while((line = gffIn.readLine()) != null) {
					if(line.startsWith("#")) // comment line
						continue;
					GFF record = gffSpecs == 2 ? new GTF(line) : new GFF3(line);
					// mask this GFF region
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
			
			// Scan BED file and write BED-DETAIL file
			if(verbose > 0) {
				System.err.println("Scanning BED file ...");
				statusTask.setInfo(" regions read");
			}
			
			boolean isHeader = true;
			String trackLine;
			while((line = bedIn.readLine()) != null) {
				if(isHeader) { // first line
					if(line.startsWith("track")) { // is track line
						trackLine = line;
						trackLine = trackLine.replaceFirst("name=(?:\\\"[^\"=]+\\\"|\\S+)", "name=\"" + trackName + "\"");
						trackLine = line.replaceFirst("description=(?:\\\"[^\"=]+\\\"|\\S+)", "description=\"" + trackDesc + "\"");
						trackLine = line.replaceFirst("type=\\w+", "type=" + trackType);
					}
					else {
						trackLine = "track name=\"" + trackName + "\" type=" + trackType + " description=\"" + trackDesc + "\"" ;
					}
					if(keepTrack)
						out.write(trackLine + "\n");
					isHeader = false;
				}
				else { // record line
					String[] fields = line.split("\t");
					assert(fields.length >= MIN_N_FIELDS);

					String chr = fields[0];
					int start = Integer.parseInt(fields[1]); // BED start is 0-based
					int end = Integer.parseInt(fields[2]);   // BED end is 1-based

					if(fix) {
						if(start < 0)
							start = 0;
						if(end > gtypeIdx.getChrLen(chr))
							end = gtypeIdx.getChrLen(chr);
					}
					// get typeSum
					Map<String, Integer> typeSum = gtypeIdx.unmaskSum(chr, start, end);
					String typeStr;
					if(!detail)
						typeStr = typeSum.isEmpty() ? unType : StringUtils.join(",", typeSum.keySet());
					else
						typeStr = typeSum.isEmpty() ? unType + ":" + (end - start) : StringUtils.join(",", typeSum);

					out.write(line + "\tGType=" + typeStr + "\n");
				} // end record line
			} // while
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
				"Options:    -g  FILE                chrom size file with 1st column the chromosome names and 2nd column their sizes, required" + newLine +
				"            -i  FILE                BED6 input file, required" + newLine +
				"            -gff  FILE              GTF/GFF3 annotation file(s) used for classification, required" + newLine +
				"            -o  FILE                BED output file with added field of genetic-type summary, required" + newLine +
				"            --keep-track  FLAG      keep track-line (and update) if exists, otherwise no track-line will be kept" + newLine +
				"            --name  STRING          name attribute of the track-line, will override the original value [outfile name]" + newLine +
				"            --desc  STRING          description attribute of the track-line, will override the original value [-name]" + newLine +
				"            --detail  FLAG          show summary of mapped feature types" + newLine +
				"            --unclassified  STRING  name for unclassified alignments [" + DEFAULT_UNCLASSIFIED_GTYPE + "]" + newLine +
				"            --tag  STRING           use value of given tag in the attrubute field (9th) instead of type field (3rd) as the genetic type, if available" + newLine +
				"            -v  FLAG                show verbose information" + newLine +
				"            -fix  FLAG              try to fix BED coordinates instead of aborting execution"
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
			else if(args[i].equals("--keep-track"))
				keepTrack = true;
			else if(args[i].equals("--name"))
				trackName = args[++i];	
			else if(args[i].equals("--desc"))
				trackDesc = args[++i];			
			else if(args[i].equals("--detail"))
				detail = true;
			else if(args[i].equals("--unclassified"))
				unType = args[++i];
			else if(args[i].equals("--tag"))
				tagName = args[++i];
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
	
	private static final String DEFAULT_UNCLASSIFIED_GTYPE = "intergenic";
	
	private static final int MIN_N_FIELDS = 3; // MIN # of fields for a record line
	private static String chrLenFile;
	private static String bedInFile;
	private static String outFile;
	private static List<String> gffFiles = new ArrayList<String>();
	private static boolean detail;
	private static String unType = DEFAULT_UNCLASSIFIED_GTYPE;
	private static String tagName;
	private static int verbose;
	private static boolean fix;

	private static GTypeIndex gtypeIdx;

	private static final int statusFreq = 10000;
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static final String trackType = "bedDetail"; // required track line info
	private static boolean keepTrack;
	private static String trackName;
	private static String trackDesc;
}
