/**
 * 
 */
package net.sf.AlignerBoost;

import static net.sf.AlignerBoost.EnvConstants.newLine;
import static net.sf.AlignerBoost.EnvConstants.progFile;

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
				// check trimmed FASTQ file
				int totalMapped = 0;
				SamReaderFactory readerFac = SamReaderFactory.makeDefault();
				samIn = readerFac.open(new File(conf.getAlignFilteredFileName()));
				for(SAMRecord record : samIn) {
					Matcher match = nrPat.matcher(record.getReadName());
					if(match.find()) // is a NR id
						totalMapped += Integer.parseInt(match.group(1));
					else
						totalMapped++;
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
