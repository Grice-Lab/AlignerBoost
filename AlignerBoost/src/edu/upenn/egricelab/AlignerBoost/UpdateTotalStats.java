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
 * a class to initialize or update library total read information
 */
package edu.upenn.egricelab.AlignerBoost;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.*;

import java.io.*;
import java.util.*;
import java.util.zip.*;

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 *
 */
public class UpdateTotalStats {

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
		try {
			// Get configs
			configs = NGSExpDesign.createNGSExpDesignsFromFile(inFile);
			// test outfile first
			if((new File(outFile)).exists()) {
				System.err.println("Output file '" + outFile + "' already exists, keep old information");
				in = new BufferedReader(new FileReader(outFile));
				header = in.readLine();
				while((line = in.readLine()) != null) {
					String[] fields = line.split("\t");
					if(fields.length >= 3)
						libInfo.put(fields[0],  line);
				}
				in.close();
			} // end if
			else
				header = "lib_name\tread_len\ttotal_num";
			
			// update output
			out = new BufferedWriter(new FileWriter(outFile));
			if(header.contains("total_num"))
				out.write(header + newLine); // keep old header
			else
				out.write("lib_name\tread_len\ttotal_num" + newLine);
			for(NGSExpDesign conf : configs) {
				System.out.println(conf.libName);
				if(libInfo.containsKey(conf.libName)) {// processed lib, keep old info
					out.write(libInfo.get(conf.libName) + newLine);
					continue;
				}
				// check FASTQ readFile
				int totalNum = 0;
				String readFile = conf.getReadFile();
				in = !readFile.endsWith(".gz") ? new BufferedReader(new FileReader(conf.getReadFile())) :
					new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readFile))));
				while((line = in.readLine()) != null) {
					if(line.startsWith("@")) {
						totalNum++;
						int readLen = in.readLine().length();
						if(readLen > conf.readLen) {
							System.err.println("Warning: Found a read for lib \"" + conf.libName + "\" of length "
									+ readLen + " that is longer than the read_len option in the config file, need to be fixed before following steps");
						}
						in.readLine();
						in.readLine();
					}
				}
				in.close();
				out.write(conf.libName + "\t" + conf.readLen + "\t" + totalNum + newLine);
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
		System.err.println("Usage:    java -jar " + progFile + " stats total <EXPERIMENT-CONFIG-INFILE> <OUTFILE>");
	}

	private static String inFile;
	private static String outFile;
	private static List<NGSExpDesign> configs;

}
