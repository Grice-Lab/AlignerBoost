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

import static edu.upenn.egricelab.AlignerBoost.EnvConstants.newLine;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.progFile;

import java.io.*;
import java.util.*;

/** Format SAM/BAM file to UCSC Wiggle (wig) file fixed step format
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class BedToWig {
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
		BufferedReader chrIn = null;
		BufferedReader bedIn = null;
		BufferedWriter out = null;
		try {
			chrIn = new BufferedReader(new FileReader(chrInFile));
			bedIn = new BufferedReader(new FileReader(bedInFile));
			out = new BufferedWriter(new FileWriter(outFile));
			// Initialize chrom-index
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
			
			// Start the processMonitor to monitor the process
			if(verbose > 0) {
				processMonitor = new Timer();
				// Start the ProcessStatusTask
				statusTask = new ProcessStatusTask("BED records processed");

				// Schedual to show the status every 1 second
				processMonitor.scheduleAtFixedRate(statusTask, 0, statusFreq);
				System.err.println("Scan BED6 file ...");
			}
			
			while((line = bedIn.readLine()) != null) {
				String[] fields = line.split("\t");
				if(fields.length < MIN_N_FIELDS) // non-record lines
					continue;
				String chr = fields[0];
				int start = Integer.parseInt(fields[1]) + 1; // start is 0-based
				int end = Integer.parseInt(fields[2]);
				int clone = useCol5Val ? Integer.parseInt(fields[4]) : 1; 
				
				if(verbose > 0)
					statusTask.updateStatus(); // Update status
				if(!chrIdx.containsKey(chr))
					continue;
				
				// check strand
				int strand = fields[5].equals("+") ? 1 : fields[5].equals("-") ? 2 : 3;
				if((strand & myStrand) == 0)
					continue;
				int[] idx = chrIdx.get(chr);
				for(int i = start + 1; i <= end; i++)
					idx[i] += clone;
				totalNum += clone;
			} // end each record

			// Terminate the monitor task and monitor
			statusTask.cancel();
			processMonitor.cancel();
			statusTask.finish();

			// Output
			if(verbose > 0)
				System.err.println("Output ...");
			if(includeTrack) // output track line
				out.write(trackHeader + "\n");
			
			for(Map.Entry<String, int[]> entry : chrIdx.entrySet()) {
				String chr = entry.getKey();
				int[] idx = entry.getValue();
				int prevStart = 0;
				// output coverage
				for(int start = 1; start < idx.length; start += step) {
					int end = start + step <= idx.length ? start + step : idx.length;
					double val = Stats.mean(idx, start, end);
					if(normRPM)
						val /= totalNum * 1e6;
					if(keep0 || val > 0) {
						if(prevStart == 0 || start - prevStart != step) // write not consecutive loc			
							out.write("fixedStep chrom=" + chr + " start=" + start + " step=" + step + "\n");
						out.write((float) val + "\n");
						prevStart = start;
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
				if(bedIn != null)
					bedIn.close();
				if(out != null)
					out.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void printUsage() {
		System.err.println("java -jar " + progFile + " utils samToWig " +
				"<-g CHR-SIZE-FILE> <-i BED6-INFILE> <-o OUTFILE> [options]" + newLine +
				"Options:    -s INT  genome strand(s) to look at, 1: plus, 2: minus, 3: both [3]" + newLine +
				"            -c/--clone-value FLAG  use BED file column 5 value as read clone" + newLine +
				"            --norm-rpm FLAG  normalize the coverage to RPM by total read number" + newLine +
				"            --no-track FLAG  do not include the 'track-line' as the first line of the Wiggle file as the UCSC required [false]" + newLine + 
				"            -name STRING  the track name used to display in UCSC Genome Browser [OUTFILE]" + newLine +
				"            -desc STRING  the description of the track used to display in UCSC Genome Browser [track name]" + newLine +
				"            -step INT step width for calculating the coverage or average coverages [1]" + newLine +
				"            -k/--keep-uncover FLAG keep 0-covered regions in wigFile [false]" + newLine +
				"            -v FLAG  show verbose information"
				);
	}
	
	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("-g"))
				chrInFile = args[++i];
			else if(args[i].equals("-i"))
				bedInFile = args[++i];
			else if(args[i].equals("-o"))
				outFile = args[++i];
			else if(args[i].equals("-s"))
				myStrand = Integer.parseInt(args[++i]);
			else if(args[i].equals("--norm-rpm"))
				normRPM = true;
			else if(args[i].equals("--no-track"))
				includeTrack = false;
			else if(args[i].equals("-name"))
				trackName = args[++i];
			else if(args[i].equals("-desc"))
				trackDesc = args[++i];
			else if(args[i].equals("-step"))
				step = Integer.parseInt(args[++i]);
			else if(args[i].equals("-k") || args[i].equals("--keep-uncover"))
				keep0 = true;
			else if(args[i].equals("c") || args[i].equals("--clone-value"))
					useCol5Val = true;
			else if(args[i].equals("-v"))
				verbose++;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'.");
		}
		// Check required options
		if(bedInFile == null)
			throw new IllegalArgumentException("-i must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-o must be specified");
		// Reformat myStrand
		if(!(myStrand >= 1 && myStrand <= 3))
			throw new IllegalArgumentException("Unknown -s option, must be 1, 2 or 3");
		// Set track name and desc
		if(trackName == null)
			trackName = outFile.replaceFirst("\\.wig$", "");
		if(trackDesc == null)
			trackDesc = trackName;
		// Finalize track header
		trackHeader += " name=" + trackName + " description=" + trackDesc;
	}

	private static String chrInFile;
	private static String bedInFile;
	private static String outFile;
	private static int myStrand = 3;
	private static boolean normRPM;
	private static int step = 1; // fixedWig step
	private static boolean useCol5Val;
	private static boolean keep0;
//	private static boolean isLog;
	private static boolean includeTrack = true; // include track line by default
	private static String trackName;
	private static String trackDesc;
	private static String trackHeader = "track type=wiggle_0";
	private static int verbose;

	private static long totalNum;
	private static Map<String, int[]> chrIdx;

	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static final int statusFreq = 10000;
	private static final int MIN_N_FIELDS = 6; // MIN # of fields in a BED6 file
}
