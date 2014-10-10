/**
 * 
 */
package net.sf.AlignerBoost;

import static net.sf.AlignerBoost.EnvConstants.newLine;
import static net.sf.AlignerBoost.EnvConstants.progFile;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * @author Qi Zheng
 *
 */
public class PrepareBestStratumCmd {

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
				String inFn = conf.getAlignFilteredFileName();
				String outFn = conf.getAlignBestStratumFileName();
				String prog = !conf.isPaired ? "stratumSE" : "stratumPE";
				String cmd = "java -jar " + progFile + " run " + prog +
						" --max-div " + conf.maxDiv * 100 + " --max-best " + conf.maxBest + " --max-report " + conf.maxReport +
						" -in " + inFn + " -out " + outFn;

				if(!(new File(outFn)).exists())
					out.write(cmd + newLine);
				else {
					System.err.println("Trimmed output file already exists, won't override");
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
		System.err.println("Usage:    java -jar " + progFile + " prepare stratum <EXPERIMENT-CONFIG-INFILE> <BASH-OUTFILE>");
	}

	private static String shPath;
	private static String inFile;
	private static String outFile;
	private static List<NGSExpDesign> configs;

}
