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
 * Filter UCSC Wiggle fixed format file
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class FilterWigFix {

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
			
			/* read given regions, if specified */
			Map<String, List<GenomeInterval>> chrRegions = new HashMap<String, List<GenomeInterval>>(); // chromosomes seen so far
			String line = null;
			while((line = regionIn.readLine()) != null) {
				String[] fields = line.split("\t");
				if(fields.length < 3) // ignore header lines
					continue;
				String chr = fields[0];
				int start = Integer.parseInt(fields[1]) + 1; // BED start is 0-based
				int end = Integer.parseInt(fields[2]);
				if(!chrRegions.containsKey(chr))
					chrRegions.put(chr, new ArrayList<GenomeInterval>());
				chrRegions.get(chr).add(new GenomeInterval(chr, start, end));
			}
			if(verbose > 0)
				System.err.println("User specified regions read from BED file");
			
			// Initialize chrom-index
			if(verbose > 0)
				System.err.println("Initialize chrom-index ...");
			
			while((line = chrIn.readLine()) != null) {
				String[] fields = line.split("\t");
				String chr = fields[0];
				int len = Integer.parseInt(fields[1]);
				if(!chrRegions.containsKey(chr))
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
			int step = 1;
			
			while((line = wigIn.readLine()) != null) {
				if(line.startsWith("#")) // ignore comments
					continue;
				else if(line.startsWith("track"))
					trackLine = line; // use the original trackLine
				else if(line.startsWith("variableStep"))
					throw new RuntimeException("Variable WIG format file found, expecting Fixed WIG format");
				else if(line.startsWith("fixedStep")) {
					Matcher match = wigFixHead.matcher(line);
					match.find();
					/* update current position and index */
					chr = match.group(1);
					loc = Integer.parseInt(match.group(2));
					step = Integer.parseInt(match.group(3));
					if(chrIdx.containsKey(chr))
						idx = chrIdx.get(chr);
					else
						throw new RuntimeException(chr + " found but not exists in the region file specified by -G");
					continue;
				}
				else { // record line
					idx[loc] = Float.parseFloat(line);
					loc += step;
				
					if(verbose > 0)
						statusTask.updateStatus(); // Update status
				}
			}

			if(verbose > 0) {
//				statusTask.cancel();
//				processMonitor.cancel();
				statusTask.finish();
				statusTask.reset();
				statusTask.setInfo("WIG records written");
				System.err.println("Output ...");
			}

			// Output
			if(trackLine != null)
				out.write(trackLine + "\n"); // always use Unix newline
			
			for(Map.Entry<String, List<GenomeInterval>> entry : chrRegions.entrySet()) {
				chr = entry.getKey();
				List<GenomeInterval> intervals = GenomeInterval.optimizeIntervals(entry.getValue());
				idx = chrIdx.get(chr);
				if(idx == null)
					continue;  // no index on this chrom

				for(GenomeInterval interval : intervals) {
					float prevVal = 0;
					// output header
					if(keep0)
						out.write("fixedStep chrom=" + chr + " start=" + interval.start + " step=" + oStep + "\n");

					// output coverage
					for(int start = interval.start; start < interval.end && start < idx.length; start += oStep) {
						int end = start + oStep <= idx.length ? start + oStep : idx.length;
						float val = (float) Stats.mean(idx, start, end);
						// output header
						if(!keep0 && val != 0 && prevVal == 0)
						out.write("fixedStep chrom=" + chr + " start=" + start + " step=" + oStep + "\n");
						if(keep0 || val != 0) {
							out.write(val + "\n");
							if(verbose > 0)
								statusTask.updateStatus();
						}
						prevVal = val;
					}
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
				"Options:    -step INT  step used for output the WigFix coverages [1]" + newLine +
				"            -k/--keep-uncover FLAG  keep 0-covered regions in WigFix file [false]" + newLine +
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
				oStep = Integer.parseInt(args[++i]);
			else if(args[i].equals("-k") || args[i].equals("--keep-uncover"))
				keep0 = true;
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
	private static int oStep = 1; // fixedWig step
	private static boolean keep0; // do not ignore 0 values
	private static int verbose;

	private static Map<String, float[]> chrIdx;

	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static final int statusFreq = 10000; // status update frequency in millisecond
	private static Pattern wigFixHead = Pattern.compile("chrom=(\\w+) start=(\\d+) step=(\\d+)");
	private static String trackLine; // optional track line message
}