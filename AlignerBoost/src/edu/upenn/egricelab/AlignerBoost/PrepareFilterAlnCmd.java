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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * @author Qi Zheng
 *
 */
public class PrepareFilterAlnCmd {

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
		
		BufferedWriter out = null;
		try {
			// Get configs
			configs = NGSExpDesign.createNGSExpDesignsFromFile(inFile);
			shPath = NGSExpDesign.SH_PATH;
			
			out = new BufferedWriter(new FileWriter(outFile));
			out.write("#!" + shPath + newLine + newLine);
			// process each lib conf
			for(NGSExpDesign conf : configs) {
				if(conf.refGenome.equals("NA"))
					continue;
				String inFn = conf.getAlignRawFileName();
				String outFn = conf.getAlignFilteredFileName();
				String dp = conf.aligner.equals("bowtie") ? " --1DP " : " "; // always enable 1DP for non-SW supported aligners
				String silent = conf.isPaired ? " --silent " : " ";
				String prog = !conf.isPaired ? "filterSE" : "filterPE";
				float minRate = conf.minAlignRate;
				String maxHit = NGSExpDesign.supportMaxHit(conf.aligner) ? " -N " + conf.maxHit : " "; 
				String knownSnp = conf.knownSnpFile != null ? " --known-SNP " + conf.knownSnpFile + " " : " ";
				String fixMD = conf.aligner.equals("seqalto") ? " --fix-MD " : " ";
				String ignoreClip = conf.hasSpliced() && !NGSExpDesign.isRNAAligner(conf.aligner) ? " --ignore-clip-penalty " : " ";
				String fragLen = conf.isPaired && !conf.hasSpliced ? " --min-frag-len " + conf.minFragLen + " --max-frag-len " + conf.maxFragLen + " " : " ";
				String est = conf.isPaired && !conf.hasSpliced ? " " : " --no-estimate ";
/*				ClipHandlingMode clipHandle = !conf.hasSpliced || NGSExpDesign.isRNAAligner(conf.aligner)
						? ClipHandlingMode.USE : ClipHandlingMode.IGNORE;*/
	/*			String cmd = "java -jar " + progFile + " run " + prog + " -r " + minRate +
						" --seed-len " + conf.seedLen + " --seed-mis " + conf.seedMis + " --seed-indel " + conf.seedIndel +
						" --all-mis " + conf.allMis + " --all-indel " + conf.allIndel + dp + silent + maxHit + ignoreClip +
						" --min-mapQ " + conf.minMapQ + " --max-best " + conf.maxBest + " --max-report " + conf.maxReport +
						" --sort-method " + conf.sortMethod + " " + knownSnp + fixMD + conf.otherFilterOpts + " -in " + inFn + " -out " + outFn;*/
				String cmd = "java -jar " + progFile + " run " + prog + " -r " + minRate +
						" --seed-len " + conf.seedLen + " --max-sensitivity " + 
						dp + silent + maxHit + ignoreClip + fragLen + est +
						" --min-mapQ " + conf.minMapQ + " --max-best " + conf.maxBest + " --max-report " + conf.maxReport +
						" --sort-method " + conf.sortMethod + " " + knownSnp + fixMD + conf.otherFilterOpts + " -in " + inFn + " -out " + outFn;
				if(!(new File(outFn)).exists())
					out.write(cmd + newLine);
				else {
					System.err.println("Filtered SAM/BAM file already exists, won't override");
					cmd = cmd.replace("\n", "\n#");
					out.write("#" + cmd + newLine);
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
		System.err.println("Usage:    java -jar " + progFile + " prepare filter <EXPERIMENT-CONFIG-INFILE> <BASH-OUTFILE>");
	}

	private static String shPath;
	private static String inFile;
	private static String outFile;
	private static List<NGSExpDesign> configs;
}
