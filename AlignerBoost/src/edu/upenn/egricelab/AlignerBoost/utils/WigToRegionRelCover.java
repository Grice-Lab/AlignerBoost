/*******************************************************************************
 * This file is part of AlignerBoost, a generalized software toolkit to boost
 * the NextGen sequencing (NGS) aligner precision and sensitivity.
 * Copyright (C) 2016  Qi Zheng
 *
 * AlignerBoost is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * AlignerBoost is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
/**
 * A utility class to filter UCSC Wiggle Fixed format file
 */
package edu.upenn.egricelab.AlignerBoost.utils;

import static edu.upenn.egricelab.AlignerBoost.EnvConstants.newLine;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.progFile;

import java.io.*;
import java.util.*;
import java.util.regex.*;

/**
 * Convert UCSC Wiggle fixed format file to relative region cover file in tab-delimited format
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class WigToRegionRelCover {
	/**
	 * @param args
	 */
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

		chrIdx = new HashMap<String, float[]>();
		BufferedReader chrIn = null;
		BufferedReader wigIn = null;
		BufferedReader regionIn = null;
		BufferedWriter out = null;
		try {
			chrIn = new BufferedReader(new FileReader(chrInFile));
			wigIn = new BufferedReader(new FileReader(wigInFile));
			regionIn = new BufferedReader(new FileReader(regionInFile));
			out = new BufferedWriter(new FileWriter(outFile));
			
			/* read given regions */
			String line = null;
			chrSeen = new HashSet<String>();
			while((line = regionIn.readLine()) != null) {
				String[] fields = line.split("\t");
				if(fields.length < MIN_BED4_FIELDS) // ignore header lines
					continue;
				String chr = fields[0];
				chrSeen.add(chr);
			}
			if(verbose > 0)
				System.err.println("User specified regions read from BED file");
			// reopen the file
			regionIn.close();
			regionIn = new BufferedReader(new FileReader(regionInFile));
			
			// Initialize chrom-index
			if(verbose > 0)
				System.err.println("Initialize chrom-index ...");
			
			while((line = chrIn.readLine()) != null) {
				String[] fields = line.split("\t");
				String chr = fields[0];
				int len = Integer.parseInt(fields[1]);
				if(!chrSeen.contains(chr))
					continue;
				chrIdx.put(chr, new float[len + 1]);  // Position 0 is dummy
				if(verbose > 0)
					System.err.println("  " + chr + ": " + len);
			}
			
			// Start the processMonitor to monitor the process
			if(verbose > 0) {
				processMonitor = new Timer();
				// Start the ProcessStatusTask
				statusTask = new ProcessStatusTask("WIG values processed");

				// Schedule to show the status every 1 second
				processMonitor.scheduleAtFixedRate(statusTask, 0, statusFreq);
				System.err.println("Scanning WIG file ...");
			}
			
			String chr = "";
			int loc = 0;
			float[] idx = null;
			int span = 1;
			
			while((line = wigIn.readLine()) != null) {
				if(line.startsWith("#")) // ignore comments
					continue;
				else if(line.startsWith("track")) // ignore track lines
					continue;
				else if(line.startsWith("variableStep"))
					throw new RuntimeException("Variable WIG format file found, expecting Fixed WIG format");
				else if(line.startsWith("fixedStep")) {
					Matcher match = wigFixHead.matcher(line);
					match.find();
					/* update current position and index */
					chr = match.group(1);
					loc = Integer.parseInt(match.group(2));
					span = Integer.parseInt(match.group(3));
					if(chrIdx.containsKey(chr))
						idx = chrIdx.get(chr);
					else
						throw new RuntimeException(chr + " found but not exists in the region file specified by -G");
					continue;
				}
				else { // record line
					idx[loc] = Float.parseFloat(line);
					loc += span;
				
					if(verbose > 0)
						statusTask.updateStatus(); // Update status
				}
			}

			if(verbose > 0) {
//				statusTask.cancel();
//				processMonitor.cancel();
				statusTask.finish();
				statusTask.reset();
				statusTask.setInfo("coverage records written");
				System.err.println("Output ...");
			}

			// Output
			out.write(header + "\n"); // always use Unix newline
			
			while((line = regionIn.readLine()) != null) {
				String[] fields = line.split("\t");
				if(fields.length < MIN_BED4_FIELDS)
					continue;
				chr = fields[0];
				idx = chrIdx.get(chr);
				if(idx == null)
					continue;  // no index on this chrom
				
				int rgStart = Integer.parseInt(fields[1]) + 1; // BED start is 0-based
				int rgEnd = Integer.parseInt(fields[2]);
				String rgName = fields[3];

				int start = rgStart - flank;
				int end = rgEnd + flank - 1;
				if(start < 1)
					start = 1;
				if(end >= idx.length)
					end = idx.length - 1;
				// output coverage
				for(int span_start = start; span_start <= end; span_start += step) {
					int span_end = span_start + step;
					int span_mid = (span_start + span_end) / 2;
					float val = (float) Stats.mean(idx, span_start, span_end);
					int start_dist = span_mid -rgStart;
					int end_dist = span_mid - rgEnd;
					// output
					out.write(rgName + "\t" + chr + "\t" + span_start + "\t" + (span_end - 1) + "\t" + span_mid + "\t" +
							+ step + "\t" + start_dist + "\t" + end_dist + "\t" + val + "\n");
					if(verbose > 0)
						statusTask.updateStatus();
				}
			}
			// Terminate the monitor task and monitor
			if(verbose > 0) {
				statusTask.cancel();
				processMonitor.cancel();
				statusTask.finish();
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
				if(chrIn != null)
					chrIn.close();
				if(wigIn != null)
					wigIn.close();
				if(out != null)
					out.close();
				if(regionIn != null)
					regionIn.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	private static void printUsage() {
		System.err.println("java -jar " + progFile + " utils filterWigFix " +
				"<-g CHR-SIZE-FILE> <-i WIG-INFILE> <-R REGION-BEDFILE> <-o OUTFILE> [options]" + newLine +
				"Options:    -step INT  step used for output the WigFix coverages [" + step + "]" + newLine +
				"            -flank INT  up/down stream flanking size for searching [" + flank + "]" + newLine +
				"            -v FLAG  show verbose information"
				);
	}
	
	/* static methods */
	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("-g"))
				chrInFile = args[++i];
			else if(args[i].equals("-i"))
				wigInFile = args[++i];
			else if(args[i].equals("-R"))
				regionInFile = args[++i];
			else if(args[i].equals("-o"))
				outFile = args[++i];
			else if(args[i].equals("-step"))
				step = Integer.parseInt(args[++i]);
			else if(args[i].equals("-flank"))
				flank = Integer.parseInt(args[++i]);
			else if(args[i].equals("-v"))
				verbose++;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'.");
		}
		// Check required options
		if(chrInFile == null)
			throw new IllegalArgumentException("-g must be specified");
		if(wigInFile == null)
			throw new IllegalArgumentException("-i must be specified");
		if(regionInFile == null)
			throw new IllegalArgumentException("-R must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-o must be specified");
	}
	
	private static String chrInFile;
	private static String wigInFile;
	private static String regionInFile;
	private static String outFile;
	private static int step = 1; // output step
	private static int flank = 0; // up/down stream flanking
	private static int verbose;

	private static Set<String> chrSeen;
	private static Map<String, float[]> chrIdx;

	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static final int statusFreq = 10000; // status update frequency in millisecond
	private static Pattern wigFixHead = Pattern.compile("chrom=(\\w+) start=(\\d+) step=(\\d+)");
	private static String header = "name\tchrom\tstart\tend\tmid\twidth\tstart_dist\tend_dist\tcover"; // optional track line message
	private static int MIN_BED4_FIELDS = 4;
}
