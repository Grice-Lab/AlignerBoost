/**
 * Main class for AlignerBooster
 */
package net.sf.AlignerBoost;

import java.util.Arrays;

/**
 * @author Qi Zheng
 * @version 1.1.1
 * @since 1.1
 *
 */
public class AlignerBoost {

	/**
	 * @param args arguments passed to main method
	 */
	public static void main(String[] args) {
		// Parse args list
		if(args.length == 0) {
			printProgHeader();
			printUsage();
			return;
		}
		String cmdGrp = args[0];
		String[] cmds = Arrays.copyOfRange(args, 1, args.length);
		if(args.length == 1 && !cmdGrp.equals("help")) { // no command provided and not help
			printProgHeader();
			System.err.println("Error! No \"command\" provided");
			printUsage();
			return;
		}
		else {
			// dispatch the main method calls
			String cmd = cmds[0];
			String[] opts = Arrays.copyOfRange(cmds, 1, cmds.length);
			switch(cmdGrp) {
			case "prepare":
				dispatchPrepare(cmd, opts);
				break;
			case "stats":
				dispatchStats(cmd, opts);
				break;
/*			case "utils":
				dispatchUtils(cmd, opts);
				break;*/
			case "run":
				dispatchRun(cmd, opts);
				break;
			case "help":
				printProgHeader();
				printUsage();
				break;
			default:
				break;
			}
		}
	}
	
	/**
	 * Dispatch prepare group commands
	 * @param cmd 
	 * @param cmd command + options from main args
	 */
	private static void dispatchPrepare(String cmd, String[] opts) {
		// Dispatch prepare commands by forwarding to corresponding class main method calls
		switch(cmd) {
		case "readQC":
			PrepareFastqQCCmd.main(opts);
			break;
		case "trim":
			PrepareTrimCmd.main(opts);
			break;
		case "NR":
			PrepareNRCmd.main(opts);
			break;
		case "align":
			PrepareMapCmd.main(opts);
			break;
		case "filter":
			PrepareFilterCmd.main(opts);
			break;
		case "stratum":
			PrepareBestStratumCmd.main(opts);
			break;
		default:
			System.err.println("Unknown prepare command '" + cmd + "'");
			printUsage();
			return;
		}
	}

	/**
	 * Dispatch stats group commands
	 * @param cmd 
	 * @param cmd command + options from main args
	 */
	private static void dispatchStats(String cmd, String[] opts) {
		// Dispatch stats commands by forwarding to corresponding class main method calls
		switch(cmd) {
		case "total":
			UpdateTotalStats.main(opts);
			break;
		case "trimmed":
			UpdateTrimmedStats.main(opts);
			break;
		case "NR":
			UpdateNRStats.main(opts);
			break;
		case "mapped":
			PrepareMappedStats.main(opts);
			break;
		default:
			System.err.println("Unknown stats command '" + cmd + "'");
			printUsage();
			return;
		}
	}
	
	/**
	 * Dispatch utils group commands
	 * @param cmd 
	 * @param cmd command + options from main args
	 */
/*	private static void dispatchUtils(String cmd, String[] opts) {
		// Dispatch prepare commands by forwarding to corresponding class main method calls
		switch(cmd) {
		case "sam2ABLoc":
			Sam2Loc.main(opts);
			break;
		case "sam2ABCover":
			Sam2Cover.main(opts);
			break;
		case "sam2Wig":
			Sam2Wig.main(opts);
			break;
		case "bed2ABLoc":
			Bed2Loc.main(opts);
			break;
		case "bed2ABCover":
			Bed2Cover.main(opts);
			break;
		case "bed2Wig":
			Bed2Wig.main(opts);
			break;
		default:
			System.err.println("Unknown utils command '" + cmd + "'");
			printUsage();
			return;
		}
	}*/
	
	/**
	 * Dispatch run group commands, whose main method should be generally called by BASH scripts indirectly.
	 * @param cmd 
	 * @param cmd command + options from main args
	 */
	private static void dispatchRun(String cmd, String[] opts) {
		// Dispatch prepare commands by forwarding to corresponding class main method calls
		switch(cmd) {
		case "readQC":
			PrepareFastqQCCmd.main(opts);
			break;
		case "trim":
			PrepareTrimCmd.main(opts);
			break;
		case "NR":
			PrepareNRCmd.main(opts);
			break;
		case "align":
			PrepareMapCmd.main(opts);
			break;
		case "filter":
			PrepareFilterCmd.main(opts);
			break;
		case "stratum":
			PrepareBestStratumCmd.main(opts);
			break;
		default:
			System.err.println("Unknown prepare command '" + cmd + "'");
			printUsage();
			return;
		}
	}
	
	/**
	 * Print the program header
	 */
	public static void printProgHeader() {
		String newLine = System.getProperty("line.separator", "\n");
		System.err.println("Program: " + progName + " (" + progDesc + ")");
		System.err.println("Version: " + progVer + newLine);
	}
	
	/**
	 * Print the program usage message
	 */
	public static void printUsage() {
		String newLine = System.getProperty("line.separator", "\n");
		System.err.println(
				"Usage:    java -jar " + progFile + " <command-group> [<command> [options]]" + newLine + 
				"Command-group:    prepare  prepare bash scripts for various " + progName + " procedures" + newLine +
				"                  stats    create/update per-project statistic summary after each " + progName + " step" + newLine +
				"                  utils    run varials utilities/apps come with the " + progName + " boundle" + newLine +
				"                  run      internal run commands recommended to be called withing bash scripts generated by \"prepare\" commands" + newLine +
				"                  help     print this help message and quit" + newLine + newLine +		
				"Command:          prepare  readQC       prepare NGS read QC commands" + newLine +
				"                  prepare  trim         prepare NGS read trimming commands" + newLine +
				"                  prepare  NR           prepare NGS read Non-redundant collapsing commands" + newLine +
				"                  prepare  align        prepare NGS read aligning (mapping) commands" + newLine +
				"                  prepare  filter       prepare NGS alignment filtering commands" + newLine +
				"                  prepare  stratum      prepare NGS alignment \"best-stratum\" filtering commands" + newLine +
				"                  stats    total        init total read stats of a project" + newLine +
				"                  stats    trimmed      add/update trimmed read stats" + newLine +
				"                  stats    NR           add/update NR tag stats" + newLine +
				"                  stats    mapped       add/update final mapped unique reads that passed all the filters" + newLine +
				"                  utils    sam2ABLoc    convert SAM/BAM files to AlignBoost tab-delimited alignment files" + newLine +
				"                  utils    sam2ABCover  convert SAM/BAM files to AlignBoost tab-delimited read-coverage files" + newLine +
				"                  utils    sam2Wig      convert SAM/BAM files to UCSC Wiggle files" +
				"                  utils    bed2ABLoc    convert UCSC bed files to AlignBoost tab-delimited alignment files" + newLine +
				"                  utils    bed2ABCover  convert UCSC bed files to AlignBoost tab-delimited read-coverage files" + newLine +
				"                  utils    bed2Wig      convert UCSC bed files to UCSC Wiggle files" + newLine +
                "                  run      fastqQC      get NGS read QC from FASTQ files" + newLine +
				"                  run      NR           collapse NGS reads to NR-tags" + newLine +
				"                  run      filter       re-calibarate and filter NGS alignments with various criterial" + newLine +
				"                  run      stratum      boost NGS alignments accuracy by picking only \"best-stratum\" hits"
				);
	}
	
	static final String progName = "AlignerBooster";
	static final String progVer = "v1.1";
	static final String progDesc = "A tool for boosting the accuracy of NextGen-seq aligner";
	static final String progFile = "AlignerBooster-v1.1.jar";
}
