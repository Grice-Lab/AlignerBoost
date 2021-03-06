/*******************************************************************************
 *     This file is part of AlignerBoost, a generalized software toolkit to boost
 *     the NextGen sequencing (NGS) aligner precision and sensitivity.
 *     Copyright (C) 2015  Qi Zheng
 *
 *     AlignerBoost is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     AlignerBoost is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
/**
 * Main class for AlignerBooster
 */
package edu.upenn.egricelab.AlignerBoost;

import static edu.upenn.egricelab.AlignerBoost.EnvConstants.newLine;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.progDesc;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.progFile;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.progName;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.progVer;

import java.util.Arrays;

import edu.upenn.egricelab.AlignerBoost.utils.*;

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
		else if(cmds.length == 0) { // command group only
			printProgHeader();
			printUsage(cmdGrp);
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
		case "sam2CoverSumm":
			SamToCoverSumm.main(opts);
			break;
		case "sam2Wig":
			SamToWig.main(opts);
			break;
		case "bed2Wig":
			BedToWig.main(opts);
			break;
		case "bed2AbsCover":
			BedToAbsCover.main(opts);
			break;
		case "filterSamById":
			FilterSamById.main(opts);
			break;
		case "classifySAM":
			ClassifySAM.main(opts);
			break;
		case "classifyVCF":
			ClassifyVCF.main(opts);
			break;
		case "classifyBED":
			ClassifyBED.main(opts);
			break;
		case "classSummSAM":
			ClassSummSAM.main(opts);
			break;
		case "filterWigFix":
			FilterWigFix.main(opts);
			break;
		case "filterWigVar":
			FilterWigVar.main(opts);
			break;
		case "wigFix2RelCover":
			WigFixToRegionRelCover.main(opts);
			break;
		case "wigVar2RelCover":
			WigVarToRegionRelCover.main(opts);
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
	 * Print the global program usage message
	 */
	private static void printUsage() {
		System.err.println(
				usageMsg + newLine + newLine +
				cmdGrpMsg + newLine + newLine +
				prepareCmdMsg + newLine +
				statsCmdMsg + newLine +
				coreCmdMsg + newLine +
				utilsCmdMsg + newLine);
	}

	/**
	 * Print the group usage message
	 * @param group  command group
	 */
	private static void printUsage(String group) {
		switch(group) {
		case "prepare":
			System.err.println(
					usageMsg + newLine + newLine + 
					prepareCmdMsg + newLine);
			break;
		case "stats":
			System.err.println(
					usageMsg + newLine + newLine +
					statsCmdMsg + newLine);
			break;
		case "run":
			System.err.println(
					usageMsg + newLine + newLine +
					coreCmdMsg + newLine);
			break;
		case "utils":
			System.err.println(
					usageMsg + newLine + newLine +
					utilsCmdMsg + newLine);
			break;
		default:
			System.err.println("Error! Unknown Command-Group '" + group + "' invoked");
			printUsage();
			break;
		}
	}
	
	private static final String usageMsg =
			"Usage:    java -jar " + progFile + " <command-group> [<command> [options]]";
	private static final String cmdGrpMsg =
			"Command-Groups:   prepare  prepare bash scripts for various " + progName + " procedures" + newLine +
			"                  stats    create/update per-project statistic summary after each " + progName + " step" + newLine +
			"                  utils    run varials utilities/apps come with the " + progName + " boundle" + newLine +
			"                  run      core commands recommended to be called withing bash scripts generated by the \"prepare\" commands" + newLine +
			"                  help     print this help message and quit";
	private static final String prepareCmdMsg =
			"Prepare Commands: prepare  readQC        prepare NGS read QC commands" + newLine +
			"                  prepare  trim          prepare NGS read trimming commands" + newLine +
			"                  prepare  NR            prepare NGS read Non-redundant collapsing commands" + newLine +
			"                  prepare  align         prepare NGS read aligning (mapping) commands" + newLine +
			"                  prepare  filter        prepare NGS alignment filtering commands";
	private static final String statsCmdMsg =
			"Stats Commands:   stats    total         init total read stats of a project" + newLine +
			"                  stats    trimmed       add/update trimmed read stats" + newLine +
			"                  stats    NR            add/update NR tag stats" + newLine +
			"                  stats    mapped        add/update final mapped unique reads that passed all the filters";
	private static final String coreCmdMsg =
            "Core Commands:    run      fastqQC       get NGS read QC from FASTQ files" + newLine +
			"                  run      NR            collapse NGS FASTQ reads/pairs to NR-tags/NR-pairs" + newLine +
			"                  run      filterSE      boost single-end mapping precision & sensitivity by calculating and filtering alignment mapQ" + newLine +
			"                  run      filterPE      boost paired-end mapping precision & sensitivity by calculating and filtering alignment mapQ";
	private static final String utilsCmdMsg =
			"Utils Commands:   utils    sam2AbsCover  convert SAM/BAM file to tab-delimited coverage file w/ absolute loc" + newLine +
			"                  utils    sam2RelCover  convert SAM/BAM file to tab-delimited coverage file w/ relative pos" + newLine +
			"                  utils    sam2BinCover  convert SAM/BAM file to tab-delimited coverage file w/ binned regions" + newLine +
			"                  utils    sam2RegCount  count reads from SAM/BAM file in given regions in BED file" + newLine +
			"                  utils    sam2Wig       convert SAM/BAM files to UCSC Wiggle files" + newLine +
			"                  utils    sam2CoverSumm generate basewise coverage depth summaries for SAM/BAM file" + newLine +
			"                  utils    bed2Wig       convert BED6 file to UCSC Wiggle file"+ newLine +
			"                  utils    bed2AbsCover  convert BED6 file to UCSC tab-delimited coverage file w/ absolute loc" + newLine +
			"                  utils    filterSamById filter SAM/BAM files with a given ID list" + newLine +
			"                  utils    classifySAM   quick classify SAM/BAM file given genomic annotations in GFF file(s)" + newLine +
			"                  utils    classifyVCF   quick classify VCF variation file given genomic annotations in GFF file(s)" + newLine +
			"                  utils    classifyBED   quick classify BED file given genomic annotations in GFF file(s)" + newLine +
			"                  utils    classSummSAM  quick classify and summarize SAM/BAM file given genomic annotations in GFF file(s)" + newLine +
			"                  utils    filterWigFix  filter UCSC Wiggle Fixed format file(s) with given regions in BED file" + newLine +
			"                  utils    filterWigVar  filter UCSC Wiggle Variable format file(s) with given regions in BED file" + newLine +
			"                  utils    wigFix2RelCover  convert UCSC Wiggle Fixed format file to tax-delimited coverage file in given regions" + newLine +
			"                  utils    wigVar2RelCover  convert UCSC Wiggle Variable format file to tax-delimited coverage file in given regions";
; 
}
