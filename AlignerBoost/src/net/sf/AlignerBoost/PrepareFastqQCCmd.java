/**
 * A helper class to prepare sh/bash FASTQ command file from a formated experimental design config file 
 */
package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.newLine;
import static net.sf.AlignerBoost.EnvConstants.progFile;
import static net.sf.AlignerBoost.NGSExpDesign.PROJECT_DIR;

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
					in += "-mate" + conf.mateFile;
				String outFn = conf.libName + "_QC.txt";
				if(!PROJECT_DIR.equals("."))
					outFn = PROJECT_DIR + "/" + outFn;
				String cmd = "java -jar " + progFile + " prepare readQC " + in + " -out " + outFn + " -readLen " + conf.readLen + newLine;
		
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