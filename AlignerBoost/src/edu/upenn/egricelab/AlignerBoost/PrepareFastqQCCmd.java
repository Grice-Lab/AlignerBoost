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
 * A helper class to prepare sh/bash FASTQ command file from a formated experimental design config file 
 */
package edu.upenn.egricelab.AlignerBoost;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.newLine;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.progFile;
import static edu.upenn.egricelab.AlignerBoost.NGSExpDesign.PROJECT_DIR;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 *
 */
public class PrepareFastqQCCmd {

	/**
	 * @param args  command-line opts passed to this method  (directly or via the AlignerBoost class)
	 */
	public static void main(String[] args) {
		// parse options
		if(args.length != 2) {
			printUsage();
			return;
		}
		inFile = args[0];
		outFile = args[1];
		
		BufferedWriter out = null;
		try {
			// Get configs
			configs = NGSExpDesign.createNGSExpDesignsFromFile(inFile);
			shPath = NGSExpDesign.SH_PATH;
			
			out = new BufferedWriter(new FileWriter(outFile));
			out.write("#!" + shPath + newLine + newLine);
			// process each lib conf
			for(NGSExpDesign conf : configs) {
				String in = "-in " + conf.readFile;
				if(conf.isPaired)
					in += " -mate " + conf.mateFile;
				String mateLen = conf.isPaired ? " -mateLen " + conf.mateLen : " ";
				String outFn = conf.libName + "_QC.txt";
				if(!PROJECT_DIR.equals("."))
					outFn = PROJECT_DIR + "/" + outFn;
				String cmd = "java -jar " + progFile + " run fastqQC " + in + " -out " + outFn +
						" -readLen " + conf.readLen + mateLen + newLine;
		
				if(!(new File(outFn)).exists())
					out.write(cmd);
				else {
					System.err.println("QC output file already exists, won't override");
					out.write("#" + cmd);
				}
			} // end each config
			// chmod of the outFn
			(new File(outFile)).setExecutable(true);
		}
		catch(IOException e) {
			System.err.println("Error: " + e.getMessage());
		}
		catch(IllegalArgumentException e) {
			System.err.println("Error: " + e.getMessage());
		}
		finally {
			try {
				if(out != null)
					out.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void printUsage() {
		System.err.println("Usage:    java -jar " + progFile + " prepare readQC <EXPERIMENT-CONFIG-INFILE> <BASH-OUTFILE>");
	}

	private static String shPath;
	private static String inFile;
	private static String outFile;
	private static List<NGSExpDesign> configs;
}