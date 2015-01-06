package net.sf.AlignerBoost;
import static net.sf.AlignerBoost.EnvConstants.newLine;
import static net.sf.AlignerBoost.EnvConstants.progFile;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.*;

/** A util class to get base-to-base QC report of Fastq read files
 * The variance will be calculated as one-pass algorithm VAR(X) = E(X^2) - E(X)^2
 * @Author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class FastqReadQC {
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

		// Read in all read quality
		System.err.println("Reading read quaility ...");
		BufferedReader readIn = null;
		BufferedReader mateIn = null;
		BufferedWriter out = null;
		long[] N1 = null;
		long[] N2 = null;
		long[] qualS_1 = null;
		long[] qualS_2 = null;
		long[] qualSS_1 = null;
		long[] qualSS_2 = null;
		int[] qualMin_1 = null;
		int[] qualMin_2 = null;
		int[] qualMax_1 = null;
		int[] qualMax_2 = null;
		try {
			out = new BufferedWriter(new FileWriter(outFile));
			// processing readFiles
			boolean isFirst = true;
			for(String readFile: readFiles) {
				System.err.println("Processing readFile: " + readFile);
				readIn = !readFile.endsWith(".gz") ? new BufferedReader(new FileReader(readFile)) :
					new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readFile))));
				String line = null;
				while((line = readIn.readLine()) != null) {
					if(line.startsWith("+")) { // will always meet the + first, even a qual line starts with a + addidentally
						String qual = readIn.readLine(); // Get next qual line
						if(isFirst) {
							isFirst = false;
							if(readLen == 0) // readLen not specified
								readLen = qual.length();
							N1 = new long[readLen];
							qualS_1 = new long[readLen];
							qualSS_1 = new long[readLen];
							qualMin_1 = new int[readLen];
							qualMax_1 = new int[readLen];
							// Init the min and max -1
							for(int i = 0; i < readLen; i++) {
								qualMin_1[i] = -1;
								qualMax_1[i] = -1;
							}
						}
						for(int i = 0; i < qual.length(); i++) {
							N1[i]++;
							int Q = qual.charAt(i);
							qualS_1[i] += Q;
							qualSS_1[i] += Q * Q;
							if(qualMin_1[i] == -1 || Q < qualMin_1[i]) // update min
								qualMin_1[i] = Q;
							if(qualMax_1[i] == -1 || Q > qualMax_1[i]) // update max
								qualMax_1[i] = Q;
						}
					}
				} // end of each line
				readIn.close();
			} // end of each readFile

			// processing mateFiles
			isFirst = true;
			for(String mateFile : mateFiles) {
				System.err.println("Processing mateFile: " + mateFile);
				mateIn = !mateFile.endsWith(".gz") ? new BufferedReader(new FileReader(mateFile)) :
					new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(mateFile))));
				String line = null;
				while((line = mateIn.readLine()) != null) {
					if(line.startsWith("+")) {
						String qual = mateIn.readLine(); // Get next qual line
						if(isFirst) {
							isFirst = false;
							if(mateLen == 0) // first qual line
								mateLen = qual.length();
							if(mateLen != readLen) {
								mateIn.close();
								out.close();
								throw new IllegalArgumentException("Mate length '" + mateLen + "' is different to the read length '" + readLen + "' in mateFile: '" + mateFile + "'");
							}
							N2 = new long[mateLen];
							qualS_2 = new long[mateLen];
							qualSS_2 = new long[mateLen];
							qualMin_2 = new int[mateLen];
							qualMax_2 = new int[mateLen];
							// Init the min and max with -1
							for(int i = 0; i < mateLen; i++) {
								qualMin_2[i] = -1;
								qualMax_2[i] = -1;
							}
						}
						for(int i = 0; i < qual.length(); i++) {
							N2[i]++;
							int Q = qual.charAt(i);
							qualS_2[i] += Q;
							qualSS_2[i] += Q * Q;
							if(qualMin_2[i] == -1 || Q < qualMin_2[i]) // update min
								qualMin_2[i] = Q;
							if(qualMax_2[i] == -1 || Q > qualMax_2[i]) // update max
								qualMax_2[i] = Q;
						}
					}
				} // end of each mate line
				mateIn.close();
			} // end of each mateFile

			// Auto detect qBase
			if(doAutoDetect) {
				System.err.print("Auto detecting qBase ... ");
				int minQ = Integer.MAX_VALUE;
				for(int q : qualMin_1)
					if(q != -1 && q < minQ)
						minQ = q; // Update minQ
				if(minQ > 64) // Must be in Illumina qBase
					qBase = 64;
				else
					qBase = 33;
				System.err.println(qBase);
			}

			// Output
			System.err.println("Output ...");
			// Output the header line
			out.write("pos\tN\tmin\tmax\tmean\tsd" + newLine);
			// Output each base pos, note that qBase will be subtracted from mean, min and max, but not sd because sd(X - c) = sd(X)
			for(int i = 0; i < readLen; i++) {
				int pos = i + 1;
				double mean = (double) qualS_1[i] / N1[i];
				double sd = Math.sqrt(((double) qualSS_1[i] / N1[i] - mean * mean) * N1[i] / (N1[i] - 1));  // use sample var as VAR(X) = N / (N-1) * (E(X^2) - E(X)^2)
				int min1 = qualMin_1[i] != -1 ? qualMin_1[i] - qBase : -1;
				int max1 = qualMax_1[i] != -1 ? qualMax_1[i] - qBase : -1;
				out.write(pos + "\t" + N1[i] + "\t" + min1 + "\t" + max1 + "\t" + (mean - qBase) + "\t" + sd + newLine);
			}
			// Output each mate pos, if exists
			for(int i = mateLen - 1; !mateFiles.isEmpty() && i >= 0; i--) {
				int pos = - i - 1;
				double mean = (double) qualS_2[i] / N2[i];
				double sd = Math.sqrt(((double) qualSS_2[i] / N2[i] - mean * mean) * N2[i] / (N2[i] - 1));
				int min2 = qualMin_2[i] != -1 ? qualMin_2[i] - qBase : -1;
				int max2 = qualMax_2[i] != -1 ? qualMax_2[i] - qBase : -1;
				out.write(pos + "\t" + N2[i] + "\t" + min2 + "\t" + max2 + "\t" + (mean - qBase) + "\t" + sd + newLine);
			}
			System.err.println("Done!");
		}
		catch(IOException e) {
			System.err.println(e.getMessage());
		}
		finally {
			try {
				if(readIn != null)
					readIn.close();
				if(mateIn != null)
					mateIn.close();
				if(out != null)
					out.close();
			}
			catch(IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void parseOptions(String[] args) throws IllegalArgumentException {
		readFiles = new ArrayList<String>();
		mateFiles = new ArrayList<String>();
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("-in"))
				while(i + 1 < args.length && !args[i+1].startsWith("-"))
					readFiles.add(args[++i]);
			else if(args[i].equals("-mate")) {
				while(i + 1 < args.length && !args[i+1].startsWith("-"))
					mateFiles.add(args[++i]);
			}
			else if(args[i].equals("-out"))
				outFile = args[++i];
			else if(args[i].equals("-Sanger")) {
				qBase = 33;
				doAutoDetect = false;
			}
			else if(args[i].equals("-Illumina")) {
				qBase = 64;
				doAutoDetect = false;
			}
			else if(args[i].equals("-qBase")) {
				qBase = Integer.parseInt(args[++i]);
				doAutoDetect = false;
			}
			else if(args[i].equals("-readLen")) {
				readLen = mateLen = Integer.parseInt(args[++i]);
			}
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'");
		}
		// Check options
		if(readFiles.isEmpty())
			throw new IllegalArgumentException("-in must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-out must be specified");
	}

	private static void printUsage() {
		System.err.println(
				"Usage:    java -jar " + progFile + " run fastqQC <-in FASTQ-INFILE [FASTQ-INFILE2 ...]>" +
						" [-mate <MATE-INFILE> [MATE-INFILE2 ...]] -<-out OUTFILE>" +
						" [-Sanger] [-Illumina] [-qBase <int>] [-readLen <int>]" + newLine +
						"Options:    -in REQUIRED FILE  FASTQ files for single-end or forward paired-end reads, multiple files should be separated by space (support .gz compressed files)" + newLine +
						"            -mate FASTQ FILE  files for reverse paired-end reads, multiple files should be separated by space (support .gz compressed files)" + newLine +
						"            -out REQUIRED FILE  OUTPUT file" + newLine + 
						"            -Sanger FLAG  use Sanger ascii offset, equivilent to set qBase=33; default is to auto-detect" + newLine +
						"            -Illumina FLAG  use Illumina 1.5x ascii offset, equivilent to set qBase=64; default is to auto-detect" + newLine +
						"            -qBase INT  ascii offset, override -Sanger or -Illumina [auto-detect]" + newLine +
						"            -readLen INT  read length, REQUIRED if reads are of variable length (i.e. from Illumina MiSeq or PacBio)"
				);
	}

	private static String outFile;
	private static List<String> readFiles;
	private static List<String> mateFiles;
	private static boolean doAutoDetect = true;
	private static int qBase = 33; // default w/ Sanger qBase
	private static int readLen;
	private static int mateLen;
}
