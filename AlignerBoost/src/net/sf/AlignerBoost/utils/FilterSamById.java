/**
 * a util class to format SAM/BAM files to customized tab-delimited cover file
 */
package net.sf.AlignerBoost.utils;
import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

import htsjdk.samtools.*;
import static net.sf.AlignerBoost.EnvConstants.*;

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
			samOut = outFactory.makeSAMOrBAMWriter(samIn.getFileHeader().clone(), true, new File(outFile));
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
				"Options:    -v FLAG  inverse the filter, only show SAMRecords that is NOT in the ID-list" + newLine +
				"            --no-desc FLAG  remove description from SAMRecord readname (anything after first white space) before filtering"
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
