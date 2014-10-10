/**
 * a class to initialize or update library total read information
 */
package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.*;

import java.io.*;
import java.util.*;

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
				in = new BufferedReader(new FileReader(conf.getReadFile()));
				while((line = in.readLine()) != null) {
					if(line.startsWith("@")) {
						totalNum++;
						in.readLine();
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
