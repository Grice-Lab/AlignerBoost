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

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 *
 */
public class UpdateNRStats {

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
				isProcessed.put(fields[0], fields.length >= 5);
				libInfo.put(fields[0],  line);
			}
			in.close();
			
			// update output
			out = new BufferedWriter(new FileWriter(outFile));
			if(header.contains("total_NR"))
				out.write(header + newLine); // keep old header
			else
				out.write(header + "\ttotal_NR" + newLine);
			
			for(NGSExpDesign conf : configs) {
				System.out.println(conf.libName);
				if(isProcessed.get(conf.libName)) {// processed lib, keep old info
					out.write(libInfo.get(conf.libName) + newLine);
					continue;
				}
				if(!conf.doNR){
					out.write(libInfo.get(conf.libName) + "\tNA" + newLine);
					continue;
				}
				// check trimmed FASTQ file
				int totalNR = 0;
				in = new BufferedReader(new FileReader(conf.getNRReadFileName()));
				while((line = in.readLine()) != null) {
					if(line.startsWith(">"))
						totalNR++;
				}
				in.close();
				out.write(libInfo.get(conf.libName) + "\t" + totalNR + newLine);
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
		System.err.println("Usage:    java -jar " + progFile + " stats NR <EXPERIMENT-CONFIG-INFILE> <OUTFILE>");
	}

	private static String inFile;
	private static String outFile;
	private static List<NGSExpDesign> configs;

}
