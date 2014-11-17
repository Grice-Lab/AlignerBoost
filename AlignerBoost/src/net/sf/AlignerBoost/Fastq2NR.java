/**
 * a class to get NR-tags from FASTQ reads/pairs
 */
package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.newLine;
import static net.sf.AlignerBoost.EnvConstants.progFile;
import static net.sf.AlignerBoost.utils.Stats.*;

import java.io.*;
import java.util.zip.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

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
			if(!isPaired) {
				Map<String, NRTag> seq2NR = extractNRFromRead(readInFile);
				// output seqCounts
				writeNR(seq2NR, readOutFile);
			}
			else {
				Map<String, NRPair> seq2NRPair = extractNRFromRead(readInFile, mateInFile);
				// output seqCounts
				writeNR(seq2NRPair, readOutFile, mateOutFile);				
			}
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

	public static void writeNR(Map<String, NRTag> seq2NR, String readOutFile)
			throws IOException {
		readOut = new BufferedWriter(new FileWriter(readOutFile));
		int nr = 0;
		int tr = 0;
		int un = 0;
		for(Entry<String, NRTag> entry : seq2NR.entrySet()) {
			NRTag tag = entry.getValue();
			int len = tag.seq.length();
			String nrID = readLen == 0 ? "nr" + (++nr) : len < readLen ? "tr" + (++tr) : "un" + (++un);
			readOut.write("@" + nrID + ":" + tag.clone + ":" + len + "\n" + tag.seq + "\n" +
					"+\n" + tag.getNRQual() + "\n");
		}
		readOut.close();
	}

	public static void writeNR(Map<String, NRPair> seq2NRPair, String readOutFile, String mateOutFile)
			throws IOException {
		readOut = new BufferedWriter(new FileWriter(readOutFile));
		mateOut = new BufferedWriter(new FileWriter(mateOutFile));
		int nr = 0;
		int tr = 0;
		int un = 0;
		for(Entry<String, NRPair> entry : seq2NRPair.entrySet()) {
			NRPair pair = entry.getValue();
			int len = pair.getPairLength();
			String nrID = readLen == 0 ? "nr" + (++nr) : len < readLen * 2 ? "tr" + (++tr) : "un" + (++un);
			readOut.write("@" + nrID + ":" + pair.clone + ":" + len + "/1\n" + pair.readSeq + "\n" +
					"+\n" + pair.getNRReadQual() + "\n");
			mateOut.write("@" + nrID + ":" + pair.clone + ":" + len + "/2\n" + pair.mateSeq + "\n" +
					"+\n" + pair.getNRMateQual() + "\n");
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
	@SuppressWarnings("resource")
	public static Map<String, NRTag> extractNRFromRead(String readInFile)
			throws IOException {
		readIn = !readInFile.endsWith(".gz") ? new BufferedReader(new FileReader(readInFile)) :
			new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readInFile))));
		Map<String, NRTag> seq2NR = new HashMap<String, NRTag>();
		String line = null;
		while((line = readIn.readLine()) != null) {
			if(line.startsWith("@")) { // def line
				String seq = readIn.readLine(); // next line
				readIn.readLine(); // ignore 1 line
				String qual = readIn.readLine(); // next line
				if(!seq2NR.containsKey(seq))
					seq2NR.put(seq, new NRTag(seq, qual)); // construct a NRTag if not exists
				else
					seq2NR.get(seq).addRead(seq,  qual);
			}
		}
		readIn.close();
		return seq2NR;
	}

	/**
	 * Overloaded method to extract NR from read pairs
	 * @param readInFile
	 * @param mateInFile
	 * @return  map of seqCounts
	 * @throws IOException
	 */
	@SuppressWarnings("resource")
	public static Map<String, NRPair> extractNRFromRead(String readInFile, String mateInFile)
			throws IOException {
		readIn = !readInFile.endsWith(".gz") ? new BufferedReader(new FileReader(readInFile)) :
			new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readInFile))));
		mateIn = !mateInFile.endsWith(".gz") ? new BufferedReader(new FileReader(mateInFile)) :
			new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(mateInFile))));
		Map<String, NRPair> seq2NRPair = new HashMap<String, NRPair>();
		String line1 = null;
		String line2 = null;
		while((line1 = readIn.readLine()) != null && (line2 = mateIn.readLine()) != null) {
			if(line1.startsWith("@") && line2.startsWith("@")) { // def lines
				String seq1 = readIn.readLine();
				String seq2 = mateIn.readLine();
				readIn.readLine(); // Ignore sep line
				mateIn.readLine(); // Ignore sep line
				String qual1 = readIn.readLine();
				String qual2 = mateIn.readLine();
				String seq = seq1 + ":" + seq2;
				if(!seq2NRPair.containsKey(seq))
					seq2NRPair.put(seq, new NRPair(seq1, qual1, seq2, qual2)); // construct a NRPair if if not exists
				else
					seq2NRPair.get(seq).addPair(seq1, qual1, seq2, qual2);
			}
		}
		readIn.close();
		mateIn.close();
		return seq2NRPair;		
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
			else if(args[i].equals("--ascii-offset"))
				asciiOffset = Integer.parseInt(args[++i]);
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
						" [--mate-in <MATE-INFILE> <--mate-out MATE-FAS-OUTFILE> [-readLen <int>] [--ascii-offset <int>]" + newLine +
						"Options:    -in FASTQ file for (forward) reads (support .gz compressed file)" + newLine +
						"            --mate-in FASTQ file for reverse reads, optional" + newLine +
						"            -out FASTQ OUTPUT file for (forward) reads (support .gz compressed file)" + newLine + 
						"            --mate-out FASTQ OUTPUT file for reverse reads, optional" + newLine +
						"            -readLen int value of read length, optional" + newLine +
						"            --ascii-offset int value of the ascii-offset of this read, default 33"
				);
	}

	/**
	 * A Nested class member for a NR-tag
	 * @author Qi Zheng
	 *
	 */
	static class NRTag {
		/**
		 * Construct a NR-tag from a FASTQ read
		 * @param read  read seq
		 * @param qual  read quality
		 */
		public NRTag(String read, String qual) {
			assert read.length() == qual.length();
			seq = read;
			clone = 1;
			nrQual = new int[read.length()];
			for(int i = 0; i < qual.length(); i++)
				nrQual[i] += qual.charAt(i) - asciiOffset;
		}

		/**
		 * Add a read to this NRTag
		 * @param read  another read, must with identical seq
		 * @param qual  quality of the read
		 */
		public boolean addRead(String read, String qual) {
			if(!read.equals(seq))
				return false; // not added
			assert read.length() == qual.length();
			clone++;
			for(int i = 0; i < qual.length(); i++)
				nrQual[i] += qual.charAt(i) - asciiOffset;
			return true;
		}
		
		/**
		 * Compare whether this NRTag equals other object
		 * @return  true if the two NRTag has identical strings 
		 */
		@Override
		public boolean equals(Object otherObj) {
			if(!(otherObj != null && otherObj instanceof NRTag))
				return false;
			NRTag other = (NRTag) otherObj;
			return seq.equals(other.seq); // two NR-tag is differentiated only by the sequence
		}
		
		
		/**
		 * Get hashCode for this NRTag
		 * @return  hashCode from the internal seq
		 */
		@Override
		public int hashCode() {
			return seq.hashCode();
		}
		
		/**
		 * Get the ascii quality string for this NRTag
		 * @return  quality string of this NRTag
		 */
		public String getNRQual() {
			StringBuilder str = new StringBuilder();
			for(int q : nrQual) {
				int b = q + ASCII_OFFSET;
				if(b > Byte.MAX_VALUE)
					b = Byte.MAX_VALUE;
				str.append((char) b);
			}
			return str.toString();
		}
		
		private String seq;
		private int clone;
		private int[] nrQual; // collapsed log-error rate
	}
	
	/**
	 * A Nested class member for a NR-pair
	 * @author Qi Zheng
	 *
	 */
	static class NRPair {
		/**
		 * Construct a NR-pair from a FASTQ read pair
		 * @param read  read seq
		 * @param qual  read quality
		 */
		public NRPair(String readSeq, String readQual, String mateSeq, String mateQual) {
			assert readSeq.length() == readQual.length() && mateSeq.length() == mateQual.length();
			this.readSeq = readSeq;
			this.mateSeq = mateSeq;
			clone = 1;
			readNRQual = new int[readSeq.length()];
			mateNRQual = new int[mateSeq.length()];
			for(int i = 0; i < readQual.length(); i++)
				readNRQual[i] += readQual.charAt(i) - asciiOffset;
			for(int i = 0; i < mateQual.length(); i++)
				mateNRQual[i] += mateQual.charAt(i) - asciiOffset;
		}
	
		/**
		 * Add a read to this NRTag
		 * @param read  another read, must with identical seq
		 * @param qual  quality of the read
		 */
		public boolean addPair(String readSeq, String readQual, String mateSeq, String mateQual) {
			if(!(this.readSeq.equals(readSeq) && this.mateSeq.equals(mateSeq)))
				return false;
			assert readSeq.length() == readQual.length() && mateSeq.length() == mateQual.length();
			clone++;
			for(int i = 0; i < readQual.length(); i++)
				readNRQual[i] += readQual.charAt(i) - asciiOffset;
			for(int i = 0; i < mateQual.length(); i++)
				mateNRQual[i] += mateQual.charAt(i) - asciiOffset;
			return true;
		}
		
		/**
		 * Get total length of this NRPair
		 * @return  length of both read and mate seqs
		 */
		public int getPairLength() {
			return readSeq.length() + mateSeq.length();
		}

		/**
		 * Compare whether this NRPair equals other object
		 * @return  true if the two NRPair has identical strings 
		 */
		@Override
		public boolean equals(Object otherObj) {
			if(!(otherObj != null && otherObj instanceof NRPair))
				return false;
			NRPair other = (NRPair) otherObj;
			return readSeq.equals(other.readSeq) && mateSeq.equals(other.mateSeq); // two NR-tag is differentiated only by the sequence
		}
		
		
		/**
		 * Get hashCode for this NRTag
		 * @return  hashCode from the internal concatenation of read and mate seq
		 */
		@Override
		public int hashCode() {
			return (readSeq + ":" + mateSeq).hashCode();
		}
		
		/**
		 * Get the ascii quality string for the forward read
		 * @return  forward read quality string of this NRPair
		 */
		public String getNRReadQual() {
			StringBuilder str = new StringBuilder();
			for(int q : readNRQual) {
				int b = q + ASCII_OFFSET;
				if(b > Byte.MAX_VALUE)
					b = Byte.MAX_VALUE;
				str.append((char) b);
			}
			return str.toString();
		}
		
		/**
		 * Get the ascii quality string for the reverse mate
		 * @return  reverse mate quality string of this NRPair
		 */
		public String getNRMateQual() {
			StringBuilder str = new StringBuilder();
			for(int q : mateNRQual) {
				int b = q + ASCII_OFFSET;
				if(b > Byte.MAX_VALUE)
					b = Byte.MAX_VALUE;
				str.append((char) b);
			}
			return str.toString();
		}
		
		private String readSeq;
		private String mateSeq;
		private int clone;
		private int[] readNRQual; // collapsed log-error rate
		private int[] mateNRQual; // collapsed log-error rate
	}

	private static String readInFile;
	private static String mateInFile;
	private static String readOutFile;
	private static String mateOutFile;
	private static boolean isPaired;
	private static int readLen;
	private static int asciiOffset = 33;
	
	private static BufferedReader readIn;
	private static BufferedReader mateIn;
	private static BufferedWriter readOut;
	private static BufferedWriter mateOut;
}
