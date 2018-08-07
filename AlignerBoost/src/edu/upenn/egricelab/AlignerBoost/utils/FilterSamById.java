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
 * a util class to format SAM/BAM files to customized tab-delimited cover file
 */
package edu.upenn.egricelab.AlignerBoost.utils;
import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

import htsjdk.samtools.*;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.*;

/** Filter SAM/BAM file with a given list
 * @author Qi Zheng
 * @version 1.2
 * @since 1.2
 */
public class FilterSamById {
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

		SamReaderFactory inFactory = SamReaderFactory.makeDefault();
		SAMFileWriterFactory outFactory = new SAMFileWriterFactory();
		SamReader samIn = null;
		SAMFileWriter samOut = null;
		BufferedReader idIn = null;
		idFilter = new HashSet<String>();
		try {
			samIn = inFactory.open(new File(inFile));
			// clone and modify the header
			SAMFileHeader header = samIn.getFileHeader().clone(); // copy the inFile header as outFile header
			// Add new programHeader
			SAMProgramRecord progRec = new SAMProgramRecord(progName + " utils filterSamById");
			progRec.setProgramName(progName + " utils filterSamById");
			progRec.setProgramVersion(progVer);
			progRec.setCommandLine(StringUtils.join(" ", args));
			header.addProgramRecord(progRec);
			
			samOut = outFactory.makeSAMOrBAMWriter(header, true, new File(outFile));
			idIn = new BufferedReader(new FileReader(idFile));
			
			// Read in ID list
			System.err.println("Reading in ID-list");
			String id = null;
			while((id = idIn.readLine()) != null)
				idFilter.add(id);
			
			// Filter SAM/BAM file
			System.err.println("Filtering SAM/BAM file ...");
			
			for(SAMRecord record : samIn) {
				String readName = record.getReadName();
				if(noDesc)
					readName = namePat.matcher(readName).group(); // replace readName with matched pattern
				boolean flag = idFilter.contains(readName);
				if(inverse)
					flag = !flag;
				if(flag)
					samOut.addAlignment(record);
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
				if(samOut != null)
					samOut.close();
				if(idIn != null)
					idIn.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void printUsage() {
		System.err.println("java -jar " + progFile + " utils filterSamById " +
				"<-i SAM|BAM-INFILE> <-l ID-LISTFILE> <-o SAM|BAM-OUTFILE> [options]" + newLine +
				"Options:    -i  FILE         SAM/BAM input file, required" + newLine +
				"            -l  FILE         ID list, required" + newLine +
				"            -o  FILE         filtered SAM/BAM file, required" + newLine +
				"            -v  FLAG         inverse the filter, only show SAMRecords that is NOT in the ID-list" + newLine +
				"            --no-desc  FLAG  remove description from SAMRecord readname (anything after first white space) before filtering"
				);
	}
	
	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("-i"))
				inFile = args[++i];
			else if(args[i].equals("-o"))
				outFile = args[++i];
			else if(args[i].equals("-l"))
				idFile = args[++i];
			else if(args[i].equals("-v"))
				inverse = true;
			else if(args[i].equals("--no-desc"))
				noDesc = true;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'.");
		}
		// Check required options
		if(inFile == null)
			throw new IllegalArgumentException("-i must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-o must be specified");
		if(idFile == null)
			throw new IllegalArgumentException("-l must be specified");
	}

	private static String inFile;
	private static String outFile;
	private static String idFile;
	private static boolean inverse;
	private static boolean noDesc;
	private static Set<String> idFilter; // bed file regions as the query intervals

	private static Pattern namePat = Pattern.compile("^\\S+");
}
