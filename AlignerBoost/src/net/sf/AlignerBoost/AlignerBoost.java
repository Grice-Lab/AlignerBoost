/**
 * Main class for AlignerBooster
 */
package net.sf.AlignerBoost;

import net.sf.AlignerBoost.utils.*;
import static net.sf.AlignerBoost.EnvConstants.newLine;
import static net.sf.AlignerBoost.EnvConstants.progDesc;
import static net.sf.AlignerBoost.EnvConstants.progFile;
import static net.sf.AlignerBoost.EnvConstants.progName;
import static net.sf.AlignerBoost.EnvConstants.progVer;

import java.util.Arrays;

/**
 * @author Qi Zheng
 * @version 1.3
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
		if(cmdGrp.equals("help")) {
			printProgHeader();
			printUsage();
			return;
		}
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
			case "utils":
				dispatchUtils(cmd, opts);
				break;
			case "run":
				dispatchRun(cmd, opts);
				break;
			default:
				printProgHeader();
				System.err.println("Error! Unknown command-group '" + cmd + "'");
				printUsage();
				return;
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
			PrepareFilterAlnCmd.main(opts);
			break;
/*		case "stratum":
			PrepareBestStratumCmd.main(opts);
			break;*/
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
			UpdateMappedStats.main(opts);
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
	private static void dispatchUtils(String cmd, String[] opts) {
		// Dispatch prepare commands by forwarding to corresponding class main method calls
		switch(cmd) {
/*		case "sam2Loc":
			Sam2Loc.main(opts);
			break;*/
		case "sam2AbsCover":
			SamToAbsCover.main(opts);
			break;
		case "sam2RelCover":
			SamToRelCover.main(opts);
			break;
		case "sam2BinCover":
			SamToBinCover.main(opts);
			break;
		case "sam2RegCount":
			SamToRegionCount.main(opts);
			break;
		case "sam2Wig":
			SamToWig.main(opts);
			break;
		case "filterSamById":
			FilterSamById.main(opts);
			break;
/*		case "bed2Loc":
			Bed2Loc.main(opts);
			break;
		case "bed2Cover":
			Bed2Cover.main(opts);
			break;
		case "bed2Wig":
			Bed2Wig.main(opts);
			break;*/
		case "classifySAM":
			ClassifySAM.main(opts);
			break;
		case "classifyVCF":
			ClassifyVCF.main(opts);
			break;
		case "classifyBED":
			ClassifyBED.main(opts);
			break;
		default:
			System.err.println("Unknown utils command '" + cmd + "'");
			printUsage();
			return;
		}
	}
	
	/**
	 * Dispatch run group commands, whose main method should be generally called by BASH scripts indirectly.
	 * @param cmd 
	 * @param cmd command + options from main args
	 */
	private static void dispatchRun(String cmd, String[] opts) {
		// Dispatch prepare commands by forwarding to corresponding class main method calls
		switch(cmd) {
		case "fastqQC":
			FastqReadQC.main(opts);
			break;
		case "NR":
			Fastq2NR.main(opts);
			break;
		case "filterSE":
			FilterSAMAlignSE.main(opts);
			break;
		case "filterPE":
			FilterSAMAlignPE.main(opts);
			break;
		default:
			System.err.println("Unknown run command '" + cmd + "'");
			printUsage();
			return;
		}
	}
	
	/**
	 * Print the program header
	 */
	public static void printProgHeader() {
		System.err.println("Program: " + progName + " (" + progDesc + ")");
		System.err.println("Version: " + progVer + newLine);
	}
	
	/**
	 * Print the program usage message
	 */
	private static void printUsage() {
		System.err.println(
				"Usage:    java -jar " + progFile + " <command-group> [<command> [options]]" + newLine + newLine +
				"Command-group:    prepare  prepare bash scripts for various " + progName + " procedures" + newLine +
				"                  stats    create/update per-project statistic summary after each " + progName + " step" + newLine +
				"                  utils    run varials utilities/apps come with the " + progName + " boundle" + newLine +
				"                  run      internal run commands recommended to be called withing bash scripts generated by \"prepare\" commands" + newLine +
				"                  help     print this help message and quit" + newLine + newLine +		
				"Command:          prepare  readQC        prepare NGS read QC commands" + newLine +
				"                  prepare  trim          prepare NGS read trimming commands" + newLine +
				"                  prepare  NR            prepare NGS read Non-redundant collapsing commands" + newLine +
				"                  prepare  align         prepare NGS read aligning (mapping) commands" + newLine +
				"                  prepare  filter        prepare NGS alignment filtering commands" + newLine +
				"                  stats    total         init total read stats of a project" + newLine +
				"                  stats    trimmed       add/update trimmed read stats" + newLine +
				"                  stats    NR            add/update NR tag stats" + newLine +
				"                  stats    mapped        add/update final mapped unique reads that passed all the filters" + newLine +
//				"                  utils    sam2Loc       convert SAM/BAM files to AlignBoost tab-delimited alignment files" + newLine +
				"                  utils    sam2AbsCover  convert SAM/BAM file to tab-delimited coverage file w/ absolute loc" + newLine +
				"                  utils    sam2RelCover  convert SAM/BAM file to tab-delimited coverage file w/ relative pos" + newLine +
				"                  utils    sam2BinCover  convert SAM/BAM file to tab-delimited coverage file w/ binned regions" + newLine +
				"                  utils    sam2RegCount  count reads from SAM/BAM file in given regions in BED file" + newLine +
				"                  utils    sam2Wig       convert SAM/BAM files to UCSC Wiggle files" + newLine +
				"                  utils    filterSamById filter SAM/BAM files with a given ID list" + newLine +
//				"                  utils    bed2Loc       convert UCSC bed files to AlignBoost tab-delimited alignment files" + newLine +
//				"                  utils    bed2Cover     convert UCSC bed files to AlignBoost tab-delimited read-coverage files" + newLine +
//				"                  utils    bed2Wig       convert UCSC bed files to UCSC Wiggle files" + newLine +
				"                  utils    classifySAM quick classify SAM/BAM file alignments given genomic annotations" + newLine +
				"                  utils    classifyVCF   quick classify VCF variation file given genomic annotations" + newLine +
				"                  utils    classifyBED   quick classify BED file to BED detail file given genomic annotations" + newLine +
                "                  run      fastqQC       get NGS read QC from FASTQ files" + newLine +
				"                  run      NR            collapse NGS FASTQ reads/pairs to NR-tags/NR-pairs" + newLine +
				"                  run      filterSE      boost single-end alignment accuracy using a Bayesian-based mapQ estimation method" + newLine +
				"                  run      filterPE      boost paired-end alignment accuracy using a Bayesian-based mapQ estimation method"
				);
	}
}
