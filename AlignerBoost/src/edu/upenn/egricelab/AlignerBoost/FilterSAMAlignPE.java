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
package edu.upenn.egricelab.AlignerBoost;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.*;

import java.io.*;
import java.util.*;

import edu.upenn.egricelab.AlignerBoost.utils.ProcessStatusTask;
import edu.upenn.egricelab.AlignerBoost.utils.Stats;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.GroupOrder;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.variant.vcf.VCFFileReader;

/** Filter SAM/BAM single-end (SE) alignments as well as do best-stratum selection to remove too divergent hits
 * @author Qi Zheng
 * @version 1.2
 * @since 1.1
 */
public class FilterSAMAlignPE {
	public static void main(String[] args) {
		if(args.length == 0) {
			printUsage();
			return;
		}
		try {
			parseOptions(args);
		}
		catch(IllegalArgumentException e) {
			System.err.println("Error: " + e.getMessage());
			printUsage();
			return;
		}

		// Read in chrList, if specified
		if(chrFile != null) {
			chrFilter = new HashSet<String>();
			try {
				BufferedReader chrFilterIn = new BufferedReader(new FileReader(chrFile));
				String chr = null;
				while((chr = chrFilterIn.readLine()) != null)
					chrFilter.add(chr);
				chrFilterIn.close();
				if(verbose > 0)
					System.err.println("Only looking at alignments on " + chrFilter.size() + " specified chromosomes");
			}
			catch(IOException e) {
				System.err.println("Error: " + e.getMessage());
				return;
			}
		}

		if(verbose > 0) {
			// Start the processMonitor
			processMonitor = new Timer();
			// Start the ProcessStatusTask
			statusTask = new ProcessStatusTask();
			// Schedule to show the status every 1 second
			processMonitor.scheduleAtFixedRate(statusTask, 0, statusFreq);
		}
		
		// Read in known SNP file, if specified
		if(knownSnpFile != null) {
			if(verbose > 0)
				System.err.println("Checking known SNPs from user specified VCF file");
			knownVCF = new VCFFileReader(new File(knownSnpFile));
		}
		
		SamReaderFactory readerFac = SamReaderFactory.makeDefault();
		SAMFileWriterFactory writerFac = new SAMFileWriterFactory();
		if(!isSilent)
			readerFac.validationStringency(ValidationStringency.LENIENT); // use LENIENT stringency
		else
			readerFac.validationStringency(ValidationStringency.SILENT); // use SILENT stringency

		SamReader in = readerFac.open(new File(inFile));
		SAMFileHeader inHeader = in.getFileHeader();
		if(inHeader.getGroupOrder() == GroupOrder.reference && inHeader.getSortOrder() == SortOrder.coordinate)
			System.err.println("Warning: Input file '" + inFile + "' might be sorted by coordinate and cannot be correctly processed!");

		SAMFileHeader header = inHeader.clone(); // copy the inFile header as outFile header
		//System.err.println(inFile + " groupOrder: " + in.getFileHeader().getGroupOrder() + " sortOrder: " + in.getFileHeader().getSortOrder());
		// reset the orders
		header.setGroupOrder(groupOrder);
		header.setSortOrder(sortOrder);

		SAMFileWriter out = OUT_IS_SAM ? writerFac.makeSAMWriter(header, false, new File(outFile)) : writerFac.makeBAMWriter(header, false, new File(outFile));

		// write SAMHeader
		String prevID = null;
		SAMRecord prevRecord = null;
		List<SAMRecord> alnList = new ArrayList<SAMRecord>();
		List<SAMRecordPair> alnPEList = null;
		// check each alignment
		SAMRecordIterator results = in.iterator();
		if(verbose > 0) {
			System.err.println("Filtering alignments ...");
			statusTask.reset();
			statusTask.setInfo("alignments processed");
		}
		while(results.hasNext()) {
			SAMRecord record = results.next();
			if(verbose > 0)
				statusTask.updateStatus();
			String ID = record.getReadName();
			// fix read and quality string for this read, if is a secondary hit from multiple hits, used for BWA alignment
			if(ID.equals(prevID) && record.getReadLength() == 0)
				SAMAlignFixer.fixSAMRecordRead(record, prevRecord);
			if(chrFilter != null && !chrFilter.contains(record.getReferenceName())) {
				prevID = ID;
				prevRecord = record;
				continue;
			}
			// fix alignment, ignore if failed (unmapped or empty)
			if(!SAMAlignFixer.fixSAMRecord(record, knownVCF, DO_1DP)) {
				prevID = ID;
				prevRecord = record;
				continue;
			}
			if(!record.getReadPairedFlag()) {
				System.err.println("Error: alignment is not from a paired-end read at\n" + record.getSAMString());
				out.close();
				return;
			}

			if(!ID.equals(prevID) && prevID != null || !results.hasNext()) { // a non-first new ID meet, or end of alignments
				// create alnPEList from filtered alnList
				alnPEList = createAlnPEListFromAlnList(alnList);
				//System.err.printf("%d alignments for %s transformed to %d alnPairs%n", alnList.size(), prevID, alnPEList.size());
				// calculate posterior mapQ for each pair
				calcPEHitPostP(alnPEList);
				// filter PEhits
				filterPEHits(alnPEList, MIN_INSERT, MAX_SEED_MIS, MAX_SEED_INDEL, MAX_ALL_MIS, MAX_ALL_INDEL, MIN_IDENTITY, MIN_MAPQ);
				// sort the list first with an anonymous class of comparator, with DESCREASING order
				Collections.sort(alnPEList, Collections.reverseOrder());
				if(MAX_BEST > 0 && alnPEList.size() > MAX_BEST) {// too much best hits, ignore this read
					alnList.clear();
					alnPEList.clear();
				}
				// report remaining secondary alignments, up-to MAX_REPORT
				for(int i = 0; i < alnPEList.size() && (MAX_REPORT == 0 || i < MAX_REPORT); i++) {
					if(doUpdateBit)
						alnPEList.get(i).setNotPrimaryAlignmentFlags(i != 0);
					if(alnPEList.get(i).fwdRecord != null)
						out.addAlignment(alnPEList.get(i).fwdRecord);
					if(alnPEList.get(i).revRecord != null)
						out.addAlignment(alnPEList.get(i).revRecord);
				}
				// reset list
				alnList.clear();
				alnPEList.clear();
			}
			// update
			if(!ID.equals(prevID)) {
				prevID = ID;
				prevRecord = record;
			}
			alnList.add(record);
		} // end while
		try {
			in.close();
			out.close();
		}
		catch(IOException e) {
			System.err.println(e.getMessage());
		}
		// Terminate the monitor task and monitor
		if(verbose > 0) {
			statusTask.cancel();
			statusTask.finish();
			processMonitor.cancel();
		}
	}

	// a nested class for keeping a pair of SAMRecord for PE alignment
	static class SAMRecordPair implements Comparable<SAMRecordPair> {
		public SAMRecordPair(SAMRecord fwdRecord, SAMRecord revRecord) throws IllegalArgumentException {
			if(fwdRecord == null && revRecord == null)
				throw new IllegalArgumentException("forward and reverse SAMRecord cannot be both null");
			this.fwdRecord = fwdRecord;
			this.revRecord = revRecord;
		}

		/** Get whether this SAMRecordPair is actually paired
		 * @return  true if both forward and reverse record are not null
		 */
		public boolean isPaired() {
			return fwdRecord != null && revRecord != null;
		}
		
		/** Set the secondary flag for both SAMRecord of this Pair
		 * @param flag  flag to set
		 */
		public void setNotPrimaryAlignmentFlags(boolean flag) {
			if(fwdRecord != null)
				fwdRecord.setNotPrimaryAlignmentFlag(flag);
			if(revRecord != null)
				revRecord.setNotPrimaryAlignmentFlag(flag);
		}

		/** Get paired-end identity for a AlignmentPair
		 * @return overall identity of the pair
		 */
		public float getPEIdentity() {
			int PEInsertLen = 0;
			int PEMis = 0;
			int PEIndel = 0;
			if(fwdRecord != null) {
				PEInsertLen += getSAMRecordInsertLen(fwdRecord);
				PEMis += getSAMRecordPercentAllMis(fwdRecord);
				PEIndel += getSAMRecordPercentAllIndel(fwdRecord);
			}
			if(revRecord != null) {
				PEInsertLen += getSAMRecordInsertLen(revRecord);
				PEMis += getSAMRecordPercentAllMis(revRecord);
				PEIndel += getSAMRecordPercentAllIndel(revRecord);
			}
			return 1 - (PEMis + PEIndel) / (float) PEInsertLen;
		}
		
		/** Get paired-end align score for a AlignmentPair
		 * @return sum of the align score of the pair
		 */
		public int getPEAlignScore() {
			int PEScore = 0;
			if(fwdRecord != null)
				PEScore += getSAMRecordAlignScore(fwdRecord);
			if(revRecord != null)
				PEScore += getSAMRecordAlignScore(revRecord);
			return PEScore;
		}

		/** Get PE insert length
		 * @return PE insert length as the inferred insert size if paired, or the insertLen of not paired
		 */
		public int getPEInsertLen() {
			int len = 0;
			if(fwdRecord != null)
				len += getSAMRecordInsertLen(fwdRecord);
			if(revRecord != null)
				len += getSAMRecordInsertLen(revRecord);
			return len;
		}
		
		/** Get PE log-likelihood
		 * @return PE align log-likelihood
		 */
		public double getPEAlignLik() {
			if(isPaired())
				return FilterSAMAlignSE.getSAMRecordAlignLikelihood(fwdRecord) + FilterSAMAlignSE.getSAMRecordAlignLikelihood(revRecord);
			double log10Lik = fwdRecord != null ? FilterSAMAlignSE.getSAMRecordAlignLikelihood(fwdRecord) :
				FilterSAMAlignSE.getSAMRecordAlignLikelihood(revRecord); 
			byte[] qual = fwdRecord != null ? fwdRecord.getBaseQualities() : revRecord.getBaseQualities();
			// treat the missing mate as all SOFT-CLIPPED with same quality
			for(int q : qual)
				log10Lik += q / - Stats.PHRED_SCALE * SAMAlignFixer.CLIP_PENALTY;
			return log10Lik;
		}

		/**
		 * Get PE mapQ, which is the same for both forward and reverse read
		 * @return  mapQ either from fwdRecord or revRecord, which one is not null
		 */
		public int getPEMapQ() {
			return fwdRecord != null ? fwdRecord.getMappingQuality() : revRecord.getMappingQuality(); 
		}
		
		/**
		 * Set mapQ to an AlignRecordPair
		 * @param mapQ  mapQ to be set to both pairs
		 */
		public void setPEMapQ(int mapQ) {
			if(fwdRecord != null)
				fwdRecord.setMappingQuality(mapQ);
			if(revRecord != null)
				revRecord.setMappingQuality(mapQ);
		}
		
		/**
		 * Get SAMString for this pair
		 * @return  SAMString for non-null mate
		 */
		public String getSAMString() {
			StringBuilder sam = new StringBuilder();
			if(fwdRecord != null)
				sam.append(fwdRecord.getSAMString());
			if(revRecord != null)
				sam.append(revRecord.getSAMString());
			return sam.toString();
		}
		
		/** implements the Comparable method
		 * @return  the difference between the mapQ value
		 */
		@Override
		public int compareTo(SAMRecordPair that) {
			return getPEMapQ() - that.getPEMapQ();
		}

		/** Override the hashCode method
		 * @return  the mapQ as the hashCode
		 */
		@Override
		public int hashCode() {
			return getPEMapQ();
		}
		
		private SAMRecord fwdRecord;
		private SAMRecord revRecord;

	}
	
	public static List<SAMRecordPair> createAlnPEListFromAlnList(List<SAMRecord> alnList) {
		if(alnList == null)
			return null;
		List<SAMRecordPair> alnPEList = new ArrayList<SAMRecordPair>();
		for(int i = 0; i < alnList.size(); i++) {
			SAMRecord currAln = alnList.get(i);
			SAMRecord nextAln = i + 1 < alnList.size() ? alnList.get(i+1) : null;
			boolean currIsFirst = currAln.getFirstOfPairFlag();
			boolean nextIsFirst = nextAln != null && nextAln.getFirstOfPairFlag();
			int currTLen = alnList.get(i).getInferredInsertSize();
			int nextTLen = i + 1 < alnList.size() ? alnList.get(i+1).getInferredInsertSize() : 0;
			if(currIsFirst ^ nextIsFirst // this is forward, next is reverse or vice versa
					&& Math.abs(currTLen) > 0 && Math.abs(currTLen) == Math.abs(nextTLen)) { // is a align pair on the same Chromosome
				alnPEList.add(currIsFirst ? new SAMRecordPair(currAln, nextAln) : new SAMRecordPair(nextAln, currAln));
				i++; // advance through nextAln
			}
			else 
				if(!noMix)// not a pair, deal with current one only
					alnPEList.add(currIsFirst ? new SAMRecordPair(currAln, null) : new SAMRecordPair(null, currAln));
		}
		return alnPEList;
	}
	
	private static void printUsage() { 
		System.err.println("Usage:   java -jar " + progFile + " run filterSE " +
				"<-in SAM|BAM-INPUT> <-out SAM|BAM-OUTPUT> [options]" + newLine +
				"Options:    --min-insert INT  minimum insert length (excluding adapters) of a read for unamabiguous alignment [15]" + newLine +
				"            --seed-len INT  seed length for Burrows-Wheeler algorithm dependent aligners [25]" + newLine +
				"            --seed-mis DOUBLE  %mismatches allowed in seed region (by --seed-len) [4]" + newLine +
				"            --all-mis  DOUBLE  %mismatches allowed in the entire insert region (excluding clipped/intron regions) [6]" + newLine +
				"            --all-indel DOUBLE  %in-dels allowed in the entire insert region [2]" + newLine +
				"            -i/--identity DOUBLE  mimimum %identity allowd for the alignment as 100 - (%mismatches+%in-dels) [0]" + newLine +
				"            --1DP FLAG  enable 1-dimentional dymamic programming insert re-assesment, useful for non-local aligners, i.e. bowtie" + newLine +
				"            --match-score INT  match score for 1DP and calculating mapQ [1]" + newLine +
				"            --mis-score INT  mismatch score for 1DP and calculating mapQ [-2]" + newLine +
				"            --gap-open-penalty INT  gap open penalty for 1DP and calculating mapQ [4]" + newLine +
				"            --gap-ext-penalty INT  gap extension penalty for 1DP and calculating mapQ [1]" + newLine +
				"            --clip-penalty INT  additional penalty for soft or hard clipped bases for calculating mapQ [0]" + newLine +
				"            --known-SNP-penalty INT  known SNP penalty for calculating mapQ [1]" + newLine +
				"            --known-INDEL-penalty INT  known IN-DEL penalty for calculating mapQ [2]" + newLine +
				"            --known-MULTISUBSTITUTION-penalty INT  known large/multi-substitution penalty for calculating mapQ [3]" + newLine +
				"            --out-SAM FLAG  write SAM text output instead of BAM binary output" + newLine +
				"            --silent FLAG  ignore certain SAM format errors such as empty reads" + newLine +
				"            --min-mapQ INT  min mapQ calculated with Bayesian method [0]" + newLine +
				"            --max-best INT  max allowed best-stratum hits to report for a given read, 0 for no limit [0]" + newLine +
				"            --max-report INT  max report valid hits determined by --min-mapQ and --max-best, 0 for no limit [0]" + newLine +
				"            --no-update-bit FLAG  do not update the secondary alignment bit flag (0x100) after filtering" + newLine +
				"            --best-only FLAG  only report unique best hit, equivelant to --max-best 1 --max-report 1" + newLine +
				"            --best FLAG  report the best hit, ignore any secondary hit, equivelant to --max-best 0 --max-report 1" + newLine +
				"            --sort-method STRING  sorting method for output SAM/BAM file, must be \"none\", \"name\" or \"coordinate\" [none]" + newLine +
				"            --chrom-list FILE  pre-filtering file containing one chromosome name per-line" + newLine +
				"            --known-SNP FILE  known SNP file in vcf/gvcf format (v4.0+, .gz supported), used for calculating mapQ" + newLine +
				"            --AF-tag STRING  Allele Frequency Tag in VCF file to check/use for determining penaltyScores for known SNPs, use NULL to disable [AF]" + newLine +
				"            -v FLAG  show verbose information"
				);
	}

	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++)
			if(args[i].equals("-in"))
				inFile = args[++i];
			else if(args[i].equals("-out"))
				outFile = args[++i];
			else if(args[i].equals("--min-insert"))
				MIN_INSERT = Integer.parseInt(args[++i]);
			else if(args[i].equals("--seed-len"))
				SAMAlignFixer.setSEED_LEN(Integer.parseInt(args[++i]));
			else if(args[i].equals("--seed-mis"))
				MAX_SEED_MIS = Double.parseDouble(args[++i]);
			else if(args[i].equals("--all-mis"))
				MAX_ALL_MIS = Double.parseDouble(args[++i]);
			else if(args[i].equals("--all-indel"))
				MAX_ALL_INDEL = Double.parseDouble(args[++i]);
			else if(args[i].equals("-i") || args[i].equals("--identity"))
				MIN_IDENTITY = Double.parseDouble(args[++i]);
			else if(args[i].equals("--1DP"))
				DO_1DP = true;
			else if(args[i].equals("--match-score"))
				SAMAlignFixer.setMATCH_SCORE(Integer.parseInt(args[++i]));
			else if(args[i].equals("--mis-score"))
				SAMAlignFixer.setMIS_SCORE(Integer.parseInt(args[++i]));
			else if(args[i].equals("--gap-open-penalty"))
				SAMAlignFixer.setGAP_OPEN_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--gap-ext-penalty"))
				SAMAlignFixer.setGAP_EXT_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--clip-penalty"))
				SAMAlignFixer.setCLIP_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--known-SNP-penalty"))
				SAMAlignFixer.setKNOWN_SNP_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--known-INDEL-penalty"))
				SAMAlignFixer.setKNOWN_INDEL_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--out-SAM"))
				OUT_IS_SAM = true;
			else if(args[i].equals("--silent"))
				isSilent = true;
			else if(args[i].equals("--no-mix"))
				noMix = true;
			else if(args[i].equals("--min-mapQ"))
				MIN_MAPQ = Integer.parseInt(args[++i]);
			else if(args[i].equals("--max-best"))
				MAX_BEST = Integer.parseInt(args[++i]);
			else if(args[i].equals("--max-report"))
				MAX_REPORT = Integer.parseInt(args[++i]);
			else if(args[i].equals("--no-update-bit"))
				doUpdateBit = false;
			else if(args[i].equals("--best-only")) {
				MAX_BEST = 1;
				MAX_REPORT = 1;
			}
			else if(args[i].equals("--best")) {
				MAX_BEST = 0;
				MAX_REPORT = 1;
			}
			else if(args[i].equals("--sort-method")) {
				switch(args[++i]) {
				case "none":
					groupOrder = GroupOrder.none;
					sortOrder = SortOrder.unsorted;
					break;
				case "name":
					groupOrder = GroupOrder.query;
					sortOrder = SortOrder.queryname;
					break;
				case "coordinate":
					groupOrder = GroupOrder.reference;
					sortOrder = SortOrder.coordinate;
					break;
				default:
					throw new IllegalArgumentException("--sort-method must be \"none\", \"name\" or \"coordinate\"");
				}
			}
			else if(args[i].equals("--chrom-list"))
				chrFile = args[++i];
			else if(args[i].equals("--known-SNP"))
				knownSnpFile = args[++i];
			else if(args[i].equals("--AF-tag"))
				SAMAlignFixer.setAFTag(args[++i]);
			else if(args[i].equals("-v"))
				verbose++;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'");
		// Check required options
		if(inFile == null)
			throw new IllegalArgumentException("-in must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-out must be specified");
		// Check other options
		if(MAX_SEED_MIS < 0 || MAX_SEED_MIS > 100)
			throw new IllegalArgumentException("--seed-mis must be between 0 to 100");
		if(MAX_ALL_MIS < 0 || MAX_ALL_MIS > 100)
			throw new IllegalArgumentException("--all-mis must be between 0 to 100");
		if(MAX_ALL_INDEL < 0 || MAX_ALL_INDEL > 100)
			throw new IllegalArgumentException("--all-indel must be between 0 to 100");
		if(!(MIN_IDENTITY >= 0 && MIN_IDENTITY <= 100))
			throw new IllegalArgumentException("-i/--identity must be between 0 to 100");
		else
			MIN_IDENTITY /= 100.0; // use absolute identity
		if(OUT_IS_SAM && outFile.endsWith(".bam"))
			System.err.println("Warning: output file '" + outFile + "' might not be SAM format");
		if(MIN_MAPQ < 0)
			throw new IllegalArgumentException("--max-div must be non negative");
		if(MAX_BEST < 0)
			throw new IllegalArgumentException("--max-best must be non negative integer");
		if(MAX_REPORT < 0)
			throw new IllegalArgumentException("--max-report must be non negative integer");
	}

	/** get align length from AlignerBoost internal tag
	 * @return the identity if tag "XA" exists
	 * throws {@RuntimeException} if tag "XA" not exists
	 */
	static int getSAMRecordAlignLen(SAMRecord record) throws RuntimeException {
		return record.getIntegerAttribute("XA");
	}
	
	/** get align insert length from AlignerBoost internal tag
	 * @return the identity if tag "XL" exists
	 * throws {@RuntimeException} if tag "XL" not exists
	 */
	static int getSAMRecordInsertLen(SAMRecord record) throws RuntimeException {
		return record.getIntegerAttribute("XL");
	}

	/** get align identity from AlignerBoost internal tag
	 * @return the identity if tag "XI" exists
	 * throws {@RuntimeException} if tag "XI" not exists
	 */
	static float getSAMRecordIdentity(SAMRecord record) throws RuntimeException {
		return record.getFloatAttribute("XI");
	}

	/** get align %seed mismatch from AlignerBoost internal tag
	 * @return the %seed mismatch if tags "YX" and "YL" exist
	 * throws {@RuntimeException} if tags "YX" and "YL" not exist
	 */
	static float getSAMRecordPercentSeedMis(SAMRecord record) throws RuntimeException {
		return 100f * record.getIntegerAttribute("YX") / record.getIntegerAttribute("YL");
	}

	/** get align %seed indel from AlignerBoost internal tag
	 * @return the %seed indel if tags "YG" and "YL" exist
	 * throws {@RuntimeException} if tags "YG" and "YL" not exist
	 */
	static float getSAMRecordPercentSeedIndel(SAMRecord record) throws RuntimeException {
		return 100f * record.getIntegerAttribute("YG") / record.getIntegerAttribute("YL");
	}

	/** get align %seed mismatch from AlignerBoost internal tag
	 * @return the %seed mismatch if tags "ZX" and "XL" exist
	 * throws {@RuntimeException} if tags "ZX" and "XL" not exist
	 */
	static float getSAMRecordPercentAllMis(SAMRecord record) throws RuntimeException {
		return 100f * record.getIntegerAttribute("ZX") / record.getIntegerAttribute("XL");
	}

	/** get align %seed indel from AlignerBoost internal tag
	 * @return the %seed indel if tags "ZG" and "XL" exist
	 * throws {@RuntimeException} if tags "ZG" and "XL" not exist
	 */
	static float getSAMRecordPercentAllIndel(SAMRecord record) throws RuntimeException {
		return 100f * record.getIntegerAttribute("ZG") / record.getIntegerAttribute("XL");
	}
	
	/**
	 * get internal align score
	 * @param record  SAMRecord to look at
	 * @return  align score
	 */
	static int getSAMRecordAlignScore(SAMRecord record) {
		return record.getIntegerAttribute("XQ");
	}

	/**
	 * Calculate the posterior probability mapQ value (in phred scale) using the Bayesian method
	 * @param recordList
	 * 
	 */
	private static double[] calcPEHitPostP(List<SAMRecordPair> alnPEList) {
		if(alnPEList == null) // return null for null list
			return null;
		if(alnPEList.isEmpty())
			return new double[0]; // return empty array for empty list
		
		int nPairs = alnPEList.size();
		// get un-normalized posterior probs
		double[] postP = new double[nPairs];
		for(int i = 0; i < nPairs; i++)
			// get postP as priorP * likelihood, with prior proportional to the alignLength
			postP[i] = alnPEList.get(i).getPEInsertLen() * Math.pow(10.0,  alnPEList.get(i).getPEAlignLik());
		// normalize postP
		FilterSAMAlignSE.normalizePostP(postP);
		// reset the mapQ values
		for(int i = 0; i < nPairs; i++) {
			//recordList.get(i).setAttribute("XP", Double.toString(postP[i]));
			double mapQ = Math.round(Stats.phredP2Q(1 - postP[i]));
			if(Double.isNaN(mapQ) || Double.isInfinite(mapQ)) // is NaN or isInfinite
				alnPEList.get(i).setPEMapQ(INVALID_MAPQ);
			else
				alnPEList.get(i).setPEMapQ(mapQ > MAX_MAPQ ? MAX_MAPQ : (int) mapQ);
		}
		return postP;
	}

	private static int filterPEHits(List<SAMRecordPair> alnPEList, int minInsert,
			double maxSeedMis, double maxSeedIndel, double maxAllMis, double maxAllIndel, double minIdentity, int minQ) {
		int n = alnPEList.size();
		int removed = 0;
		for(int i = n - 1; i >= 0; i--) { // search backward for maximum performance
			SAMRecordPair pair = alnPEList.get(i);
			if(!(  (pair.fwdRecord == null || getSAMRecordInsertLen(pair.fwdRecord) >= minInsert
					&& getSAMRecordPercentSeedMis(pair.fwdRecord) <= maxSeedMis
					&& getSAMRecordPercentSeedIndel(pair.fwdRecord) <= maxSeedIndel
					&& getSAMRecordPercentAllMis(pair.fwdRecord) <= maxAllMis
					&& getSAMRecordPercentAllIndel(pair.fwdRecord) <= maxAllIndel
					&& getSAMRecordIdentity(pair.fwdRecord) >= minIdentity)
				&& (pair.revRecord == null || getSAMRecordInsertLen(pair.revRecord) >= minInsert
					&& getSAMRecordPercentSeedMis(pair.revRecord) <= maxSeedMis
					&& getSAMRecordPercentSeedIndel(pair.revRecord) <= maxSeedIndel
					&& getSAMRecordPercentAllMis(pair.revRecord) <= maxAllMis
					&& getSAMRecordPercentAllIndel(pair.revRecord) <= maxAllIndel
					&& getSAMRecordIdentity(pair.revRecord)>= minIdentity)
				&& pair.getPEMapQ() >= minQ) ) {
				alnPEList.remove(i);
/*				System.err.println("Removing pair:\n" + pair.getSAMString());
				if(pair.fwdRecord != null)
					System.err.printf("fwd: seedMis:%f seedIndel:%f allMis:%f allIndel:%f mapQ:%d%n", 
							getSAMRecordPercentSeedMis(pair.fwdRecord), 
							getSAMRecordPercentSeedIndel(pair.fwdRecord), 
							getSAMRecordPercentAllIndel(pair.fwdRecord), 
							getSAMRecordPercentAllIndel(pair.fwdRecord),
							pair.getPEMapQ());
				if(pair.revRecord != null)
					System.err.printf("rev: seedMis:%f seedIndel:%f allMis:%f allIndel:%f mapQ:%d%n", 
							getSAMRecordPercentSeedMis(pair.revRecord), 
							getSAMRecordPercentSeedIndel(pair.revRecord), 
							getSAMRecordPercentAllIndel(pair.revRecord), 
							getSAMRecordPercentAllIndel(pair.revRecord),
							pair.getPEMapQ());
				System.err.printf("minInsert:%d maxSeedMis:%f maxSeedIndel:%f maxAllMis:%f maxAllIndel:%f minQ:%d%n",
						minInsert, maxSeedMis, maxSeedIndel, maxAllMis, maxAllIndel, minQ);*/
				removed++;
			}
		}
		return removed;
	}

	private static final int INVALID_MAPQ = 255;
	private static final int MAX_MAPQ = 250; // MAX meaniful mapQ value, if not 255
	private static String inFile;
	private static String outFile;
	private static String chrFile;
	private static String knownSnpFile;
	// filter options
	private static int MIN_INSERT = 15;
	private static double MAX_SEED_MIS = 4; // max % seed mismatch
	private static final double MAX_SEED_INDEL = 0; // seed indel is always not allowed
	private static double MAX_ALL_MIS = 6; // max % all mismatch
	private static double MAX_ALL_INDEL = 0; // max % all indel
	private static double MIN_IDENTITY = 0; // min identity
	private static boolean DO_1DP;
	private static boolean isSilent; // ignore SAM warnings?
	private static boolean noMix; // do not allow unpaired alignments for paired reads?
	// best stratum options
	private static int MIN_MAPQ = 0; // max divergent
	private static int MAX_BEST = 0; // no limits
	private static int MAX_REPORT = 0;
	private static boolean doUpdateBit = true;
	private static int verbose; // verbose level
	private static Set<String> chrFilter;
	private static VCFFileReader knownVCF;
	// general options
	private static GroupOrder groupOrder = GroupOrder.none;
	private static SortOrder sortOrder = SortOrder.unsorted;
	private static boolean OUT_IS_SAM; // outFile is SAM format?
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static final int statusFreq = 10000;
}