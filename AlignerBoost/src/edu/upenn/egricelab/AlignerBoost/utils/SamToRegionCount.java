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
 * A util class to count reads from SAM/BAM files in given BED regions 
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
 * @version 1.2
 * @since 1.2
 */
public class SamToRegionCount {
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

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader samIn = null;
		BufferedWriter out = null;
		BufferedReader bed6In = null;
		try {
			samIn = factory.open(new File(samInFile));
			out = new BufferedWriter(new FileWriter(outFile));
			SAMSequenceDictionary samDict = samIn.getFileHeader().getSequenceDictionary();

			// check each region and output
			bed6In = new BufferedReader(new FileReader(bed6File));
			if(verbose > 0) {
				// Start the processMonitor to monitor the process
				processMonitor = new Timer();
				// Start the ProcessStatusTask
				statusTask = new ProcessStatusTask("regions scanned");
				// Schedule to show the status every 1 second
				processMonitor.scheduleAtFixedRate(statusTask, 0, statusFreq);
				System.err.println("Scanning BED6 regions and output ...");
			}
			
			out.write("chrom\tstart\tend\tname\tcount\tstrand\n");
			String line = null;
			while((line = bed6In.readLine()) != null) {
				String[] fields = line.split("\t");
				if(fields.length < 6) // ignore header lines
					continue;
				if(verbose > 0)
					statusTask.updateStatus(); // Update status
				String chr = fields[0];
				int chrI = samIn.getFileHeader().getSequenceIndex(chr);
				if(chrI == -1) // this Region is not in the aligned chromosomes
					continue;
				int regionStart = Integer.parseInt(fields[1]) + 1; // BED start is 0-based
				int regionEnd = Integer.parseInt(fields[2]);
				int regionLen = regionEnd - regionStart + 1;
				if(regionLen <= 0)
					continue;
				int chrLen = samDict.getSequence(chrI).getSequenceLength();
				String name = fields[3];
				String regionStrand = fields[5];

				// Initialize scan-region
				int scanStart = regionStart - maxFlank >= 1 ? regionStart - maxFlank : 1;
				int scanEnd = regionEnd + maxFlank <= chrLen ? regionEnd + maxFlank : chrLen;
				// Query SAM file on the fly
				int count = 0;
				SAMRecordIterator results = samIn.query(chr, scanStart, scanEnd, false);
				while(results.hasNext()) {
					SAMRecord record = results.next();
					int readLen = record.getReadLength();
					if(record.getReferenceIndex() == -1 || readLen == 0) // non mapped read or 0-length read
						continue;
					// check relative strand
					String strand = record.getReadNegativeStrandFlag() ? "-" : "+";
					int relStrand = regionStrand.equals(".") ? 3 /* unknown */ : strand.equals(regionStrand) ? 1 /* sense */ : 2 /* antisense */;
					if((relStrand & myStrand) == 0) // unmatched strands
						continue;
					if(record.getMappingQuality() < minMapQ)
						continue;
					Matcher match = nrPat.matcher(record.getReadName()); // whether match interval nrID pattern
					int clone = match.find() ? Integer.parseInt(match.group(1)) : 1;
					
					int start = record.getAlignmentStart();
					int end = record.getAlignmentEnd();
					if(overlapRate(scanStart, scanEnd, start, end) >= minRate)
						count += clone;
				} // end each record
				results.close();
				// output
				out.write(chr + "\t" + regionStart + "\t" + regionEnd + "\t" + name + "\t" + count + "\t" + regionStrand + "\n");
			} // end each region
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
				if(out != null)
					out.close();
				if(bed6In != null)
					bed6In.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static double overlapRate(int start1, int end1, int start2, int end2) {
		// calculate overlap rate between two regions, relative to region1
		// both regions are in 1-based coordinates
		assert start1 <= end1 && start2 <= end2;
		if(!(start1 <= end2 && end1 >= start2)) // no overlap
			return 0;
		int overStart = start1 > start2 ? start1 : start2;
		int overEnd = end1 < end2 ? end1 : end2;
		return (double) (overEnd - overStart + 1) / (end1 - start1 + 1); 
	}

	private static void printUsage() {
		System.err.println("java -jar " + progFile + " utils samToRelCover " +
				"<-i SAM|BAM-INFILE> <-R BED6-FILE> <-o OUTFILE> [options]" + newLine +
				"Options:    -s INT  relative strand(s) to look at, must be 1: sense, 2: antisense or 3: [3]" + newLine +
				"            --count-soft FLAG  including soft-masked regions as covered region" + newLine +
				"            -f FLOAT minimum proportion of overlap to the region [1e-9]" + newLine +
				"            -flank INT max upsteam/downsteam positions to look at [0]" + newLine +
				"            -Q/--min-mapQ  INT minimum mapQ cutoff" + newLine +
				"            -v FLAG  show verbose information"
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
				bed6File = args[++i];
			else if(args[i].equals("-f"))
				minRate = Double.parseDouble(args[++i]);
			else if(args[i].equals("-flank"))
				maxFlank = Integer.parseInt(args[++i]);
			else if(args[i].equals("-Q") || args[i].equals("--min-mapQ"))
				minMapQ = Integer.parseInt(args[++i]);
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
		if(bed6File == null)
			throw new IllegalArgumentException("-R must be specified");
		// Reformat myStrand
		if(!(myStrand >= 1 && myStrand <= 3))
			throw new IllegalArgumentException("Unknown -s option, must be 1, 2 or 3");
	}

	private static String samInFile;
	private static String outFile;
	private static String bed6File;
	private static int myStrand = 3;
	private static double minRate = 1e-9;
	private static int maxFlank;
	private static int minMapQ;
	private static int verbose;

	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static Pattern nrPat = Pattern.compile("^(?:tr|un:nr)\\d+:(\\d+):\\d+");
	private static final int statusFreq = 10000;
}
