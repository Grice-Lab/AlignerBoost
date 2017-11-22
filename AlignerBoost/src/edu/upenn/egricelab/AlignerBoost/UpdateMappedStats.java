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
 * 
 */
package edu.upenn.egricelab.AlignerBoost;

import static edu.upenn.egricelab.AlignerBoost.EnvConstants.newLine;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.progFile;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.GroupOrder;

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 *
 */
public class UpdateMappedStats {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// parse options
		if(args.length != 2) {
			printUsage();
			return;
		}
		inFile = args[0];
		outFile = args[1];
		
		BufferedReader in = null;
		BufferedWriter out = null;
		SamReader samIn = null;
		String line = null;
		String header = null;
		Map<String, String> libInfo = new HashMap<String, String>();
		Map<String, Boolean> isProcessed = new HashMap<String, Boolean>();
		try {
			// Get configs
			configs = NGSExpDesign.createNGSExpDesignsFromFile(inFile);
			// test outfile first
			if(!(new File(outFile)).exists())
				throw new IllegalArgumentException("OUTFILE must has been exists from 'stats total' step");
			in = new BufferedReader(new FileReader(outFile));
			header = in.readLine();
			while((line = in.readLine()) != null) {
				String[] fields = line.split("\t");
				isProcessed.put(fields[0], fields.length >= 6);
				libInfo.put(fields[0],  line);
			}
			in.close();
			
			// update output
			out = new BufferedWriter(new FileWriter(outFile));
			if(header.contains("total_mapped"))
				out.write(header + newLine); // keep old header
			else
				out.write(header + "\ttotal_mapped" + newLine);
			
			for(NGSExpDesign conf : configs) {
				System.out.println(conf.libName);
				if(isProcessed.get(conf.libName)) {// processed lib, keep old info
					out.write(libInfo.get(conf.libName) + newLine);
					continue;
				}
				if(conf.refGenome.equals("NA")){
					out.write(libInfo.get(conf.libName) + "\tNA" + newLine);
					continue;
				}
				// check alignments
				long totalMapped = 0;
				int maxReport = conf.getMaxReport(); // max reported hits possible in the alignment file
				
				String prevID = "";
				Map<String, Integer> readHit = null;
				SamReaderFactory readerFac = SamReaderFactory.makeDefault();
				readerFac.validationStringency(ValidationStringency.SILENT); // set validation level to silent
				samIn = readerFac.open(new File(conf.getAlignFilteredFileName()));
				GroupOrder inOrder = samIn.getFileHeader().getGroupOrder();
				boolean isUniq = maxReport == 1 && !conf.isPaired; /* Single-end, uniquely reported alignments */
				
				if(!isUniq && inOrder != GroupOrder.query)
					readHit = new HashMap<String, Integer>();
				for(SAMRecord record : samIn) {
					String id = record.getReadName();
					int clone = 1; // default is for read
					if(conf.doNR) { // we are looking at NR tags
						Matcher match = nrPat.matcher(id);
						if(match.find())
							clone = Integer.parseInt(match.group(1));
					}
					if(isUniq || inOrder == GroupOrder.query && !id.equals(prevID)) // query-ordered, new ID found
						totalMapped += clone;
					else { // non-unique and not query-ordered, need recording found reads
						// first make sure it is AlignerBoost filtered reads
						if(record.getStringAttribute("XH") == null)
							throw new RuntimeException("Cannot get mapped read summary for non-query sorted non-AlignerBoost processed BAM file");
						if(!readHit.containsKey(id)) { // this id hasn't been seen yet
							totalMapped += clone;
							readHit.put(id, 0);
						}
						readHit.put(id, readHit.get(id) + 1); // increment hit
						if(!record.getReadPairedFlag() && readHit.get(id) >= maxReport /* this read is SE and won't be found */
								|| record.getReadPairedFlag() && readHit.get(id) >= 2 * maxReport) /* this read is PE won't be found */
							readHit.remove(id);
					}
					prevID = id;
				}
				samIn.close();
				out.write(libInfo.get(conf.libName) + "\t" + totalMapped + newLine);
			}
		}
		catch(IOException e) {
			System.err.println("Error: " + e.getMessage());
		}
		catch(IllegalArgumentException e) {
			System.err.println("Error: " + e.getMessage());
		}
		finally {
			try {
				if(in != null)
					in.close();
				if(out != null)
					out.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	private static void printUsage() {
		System.err.println("Usage:    java -jar " + progFile + " stats mapped <EXPERIMENT-CONFIG-INFILE> <OUTFILE>");
	}

	private static String inFile;
	private static String outFile;
	private static List<NGSExpDesign> configs;
	private static Pattern nrPat = Pattern.compile("^(?:tr|un:nr)\\d+:(\\d+):\\d+");

}
