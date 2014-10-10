/**
 * a class to get NR-tags from FASTQ reads/pairs
 */
package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.newLine;
import static net.sf.AlignerBoost.EnvConstants.progFile;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 *
 */
public class Fastq2NR {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Parse command-line options
		try {
			parseOptions(args);
		}
		catch(IllegalArgumentException e) {
			System.err.println("Error: " + e.getMessage());
			printUsage();
			return;
		}
		
		try {
			// get seqCounts
			Map<String, Integer> seqCounts = !isPaired ? extractNRFromRead(readInFile) :
				extractNRFromRead(readInFile, mateInFile);
			// output seqCounts
			if(!isPaired)
				writeNR(seqCounts, readOutFile);
			else
				writeNR(seqCounts, readOutFile, mateOutFile);
		}
		catch(IOException e) {
			System.err.println("IO error: " + e.getMessage());
		}
		finally {
			try {
				if(readIn != null)
					readIn.close();
				if(mateIn != null)
					mateIn.close();
				if(readOut != null)
					readOut.close();
				if(mateOut != null)
					mateOut.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	public static void writeNR(Map<String, Integer> seqCounts, String readOutFile)
			throws IOException {
		readOut = new BufferedWriter(new FileWriter(readOutFile));
		int nr = 0;
		int tr = 0;
		int un = 0;
		for(Map.Entry<String, Integer> entry : seqCounts.entrySet()) {
			String seq = entry.getKey();
			int count = entry.getValue();
			int len = seq.length();
			String nrID = readLen == 0 ? "nr" + (++nr) : len < readLen ? "tr" + (++tr) : "un" + (++un);
			readOut.write(">" + nrID + ":" + count + ":" + len + "\n" + seq + "\n");
		}
		readOut.close();
	}

	public static void writeNR(Map<String, Integer> seqCounts, String readOutFile, String mateOutFile)
			throws IOException {
		readOut = new BufferedWriter(new FileWriter(readOutFile));
		mateOut = new BufferedWriter(new FileWriter(mateOutFile));
		int nr = 0;
		int tr = 0;
		int un = 0;
		for(Map.Entry<String, Integer> entry : seqCounts.entrySet()) {
			String seq = entry.getKey();
			int count = entry.getValue();
			int len = seq.length() - 1;
			String[] seqs = seq.split(":");
			String nrID = readLen == 0 ? "nr" + (++nr) : len < readLen * 2 ? "tr" + (++tr) : "un" + (++un);
			readOut.write(">" + nrID + ":" + count + ":" + len + "/1\n" + seqs[0] + "\n");
			mateOut.write(">" + nrID + ":" + count + ":" + len + "/2\n" + seqs[1] + "\n");
		}
		readOut.close();
		mateOut.close();
	}

	/**
	 * method to extract NR from reads
	 * @param readInFile
	 * @return  map of seqCounts
	 * @throws IOException
	 */
	public static Map<String, Integer> extractNRFromRead(String readInFile)
			throws IOException {
		readIn = new BufferedReader(new FileReader(readInFile));
		Map<String, Integer> seqCounts = new HashMap<String, Integer>();
		String line = null;
		while((line = readIn.readLine()) != null) {
			if(line.startsWith("@")) { // def line
				String seq = readIn.readLine(); // next line
				if(!seqCounts.containsKey(seq))
					seqCounts.put(seq, 0); // init with 0 if not exists
				seqCounts.put(seq, seqCounts.get(seq) + 1);
				readIn.readLine(); // ignore 2 lines
				readIn.readLine();
			}
		}
		readIn.close();
		return seqCounts;
	}

	/**
	 * Overloaded method to extract NR from read pairs
	 * @param readInFile
	 * @param mateInFile
	 * @return  map of seqCounts
	 * @throws IOException
	 */
	public static Map<String, Integer> extractNRFromRead(String readInFile, String mateInFile)
			throws IOException {
		readIn = new BufferedReader(new FileReader(readInFile));
		mateIn = new BufferedReader(new FileReader(mateInFile));
		Map<String, Integer> seqCounts = new HashMap<String, Integer>();
		String line1 = null;
		String line2 = null;
		while((line1 = readIn.readLine()) != null && (line2 = mateIn.readLine()) != null) {
			if(line1.startsWith("@") && line2.startsWith("@")) { // def lines
				String seq = readIn.readLine() + ":" + mateIn.readLine(); // next lines
				if(!seqCounts.containsKey(seq))
					seqCounts.put(seq, 0); // init with 0 if not exists
				seqCounts.put(seq, seqCounts.get(seq) + 1);
				readIn.readLine(); // ignore 2 lines
				readIn.readLine();
				mateIn.readLine();
				mateIn.readLine();
			}
		}
		readIn.close();
		mateIn.close();
		return seqCounts;		
	}

	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("-in"))
				readInFile = args[++i];
			else if(args[i].equals("--mate-in")) {
				mateInFile = args[++i];
				isPaired = true;
			}
			else if(args[i].equals("-out"))
				readOutFile = args[++i];
			else if(args[i].equals("--mate-out"))
				mateOutFile = args[++i];
			else if(args[i].equals("-readLen"))
				readLen = Integer.parseInt(args[++i]);
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'");
		}
		// Check options
		if(readInFile == null)
			throw new IllegalArgumentException("-in must be specified");
		if(readOutFile == null)
			throw new IllegalArgumentException("-out must be specified");
		if(readLen < 0)
			throw new IllegalArgumentException("-readLen must be non-negative integer");
		if(mateInFile == null ^ mateOutFile == null)
			throw new IllegalArgumentException("--mate-in and --mate-out must be specified at the same time");
	}

	private static void printUsage() {
		System.err.println(
				"Usage:    java -jar " + progFile + " run NR <-in FASTQ-INFILE> <-out FAS-OUTFILE>" +
						" [--mate-in <MATE-INFILE> <--mate-out MATE-FAS-OUTFILE> [-readLen <int>]" + newLine +
						"Options:    -in FASTQ file for (forward) reads" + newLine +
						"            --mate-in FASTQ file for reverse reads, optional" + newLine +
						"            -out FASTA OUTPUT file for (forward) reads" + newLine + 
						"            --mate-out FASTA OUTPUT file for reverse reads, optional" + newLine +
						"            -readLen int value of read length, optional"
				);
	}

	private static String readInFile;
	private static String mateInFile;
	private static String readOutFile;
	private static String mateOutFile;
	private static boolean isPaired;
	private static int readLen;
	
	private static BufferedReader readIn;
	private static BufferedReader mateIn;
	private static BufferedWriter readOut;
	private static BufferedWriter mateOut;
}
