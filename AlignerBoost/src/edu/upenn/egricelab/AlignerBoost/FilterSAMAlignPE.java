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
import org.apache.commons.math3.distribution.*;

//import org.apache.commons.math3.linear.IllConditionedOperatorException;

import edu.upenn.egricelab.AlignerBoost.utils.ProcessStatusTask;
import edu.upenn.egricelab.AlignerBoost.utils.Stats;
import edu.upenn.egricelab.AlignerBoost.utils.StringUtils;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.GroupOrder;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.variant.vcf.VCFFileReader;

/** Filter SAM/BAM single-end (SE) alignments as well as do best-stratum selection to remove too divergent hits
 * By filtering SE read it will modifiy the NH tag and add the XN tag
 * Tag  Type  Description
 * NH   i     Number of reported alignments
 * XN   i     Number of total alignments satisfying the user-specified criteria except for the mapQ limitation
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
		// Add new programHeader
		SAMProgramRecord progRec = new SAMProgramRecord(progName);
		progRec.setProgramName(progName);
		progRec.setProgramVersion(progVer);
		progRec.setCommandLine(StringUtils.join(" ", args));
		header.addProgramRecord(progRec);
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
		
		// Estimate fragment length distribution by scan one-pass through the alignments
		SAMRecordIterator results = in.iterator();
		if(!NO_ESTIMATE) {
			if(verbose > 0) {
				System.err.println("Estimating insert fragment size distribution ...");
				statusTask.reset();
				statusTask.setInfo("alignments scanned");
			}
			long N = 0;
			double fragL_S = 0; // fragLen sum
			double fragL_SS = 0; // fragLen^2 sum
			while(results.hasNext()) {
				SAMRecord record = results.next();
				if(verbose > 0)
					statusTask.updateStatus();
				if(record.getFirstOfPairFlag() && !record.isSecondaryOrSupplementary()) {
					double fragLen = Math.abs(record.getInferredInsertSize());
					if(fragLen != 0 && fragLen >= MIN_FRAG_LEN && fragLen <= MAX_FRAG_LEN) { // only consider certain alignments
						N++;
						fragL_S += fragLen;
						fragL_SS += fragLen * fragLen;
					}
					// stop estimate if already enough
					if(MAX_ESTIMATE_SCAN > 0 && N >= MAX_ESTIMATE_SCAN)
						break;
				}
			}
			if(verbose > 0)
				statusTask.finish();
			// estimate fragment size
			if(N >= MIN_ESTIMATE_BASE) { // override command line values
				MEAN_FRAG_LEN = fragL_S / N;
				SD_FRAG_LEN = Math.sqrt( (N * fragL_SS - fragL_S * fragL_S) / (N * (N - 1)) );
				if(verbose > 0)
					System.err.printf("Estimated fragment size distribution: N(%.1f, %.1f)%n", MEAN_FRAG_LEN, SD_FRAG_LEN);
			}
			else {
				System.err.println("Unable to estimate the fragment size distribution due to too few observed alignments");
				System.err.println("You have to specify the '--mean-frag-len' and '--sd-frag-len' on the command line and re-run this step");
				statusTask.cancel();
				processMonitor.cancel();
				out.close();
				return;
			}
			// Initiate the normal model
			normModel = new NormalDistribution(MEAN_FRAG_LEN, SD_FRAG_LEN);
			// reset the iterator, if necessary
			if(in.type() == SamReader.Type.SAM_TYPE) {
				try {
					in.close();
				}
				catch(IOException e) {
					System.err.println(e.getMessage());
				}
				in = readerFac.open(new File(inFile));
			}
			results.close();
			results = in.iterator();
		} // end of NO_ESTIMATE

		// check each alignment again
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
			
			// fix MD:Z string for certain aligners with invalid format (i.e. seqAlto)
			if(fixMD)
				SAMAlignFixer.fixMisStr(record);
			
			// fix alignment, ignore if failed (unmapped or empty)
			if(!SAMAlignFixer.fixSAMRecord(record, knownVCF, DO_1DP)) {
				prevID = ID;
				prevRecord = record;
				continue;
			}
			if(!record.getReadPairedFlag()) {
				System.err.println("Error: alignment is not from a paired-end read at\n" + record.getSAMString());
				out.close();
				statusTask.cancel();
				processMonitor.cancel();
				return;
			}

			if(!ID.equals(prevID) && prevID != null || !results.hasNext()) { // a non-first new ID meet, or end of alignments
				// create alnPEList from filtered alnList
				alnPEList = createAlnPEListFromAlnList(alnList);
				//System.err.printf("%d alignments for %s transformed to %d alnPairs%n", alnList.size(), prevID, alnPEList.size());
				int totalPair = alnPEList.size();
				// filter highly unlikely PEhits
				filterPEHits(alnPEList, MIN_ALIGN_RATE, MIN_IDENTITY);
				// calculate posterior mapQ for each pair
				calcPEHitPostP(alnPEList, totalPair, MAX_HIT);
				// filter hits by mapQ
				if(MIN_MAPQ > 0)
					filterPEHits(alnPEList, MIN_MAPQ);
				
				// sort the list first with an anonymous class of comparator, with DESCREASING order
				Collections.sort(alnPEList, Collections.reverseOrder());				
				// control max-best
				if(MAX_BEST != 0 && alnPEList.size() > MAX_BEST) { // potential too much best hits
					int nBestStratum = 0;
					int bestMapQ = alnPEList.get(0).getPEMapQ(); // best mapQ from first PE
					for(SAMRecordPair pr : alnPEList)
						if(pr.getPEMapQ() == bestMapQ)
							nBestStratum++;
						else
							break; // stop searching for sorted list
					if(nBestStratum > MAX_BEST)
						alnPEList.clear();
				}
				// filter alignments with auxiliary filters
				if(!MAX_SENSITIVITY)
					filterPEHits(alnPEList, MAX_SEED_MIS, MAX_SEED_INDEL, MAX_ALL_MIS, MAX_ALL_INDEL);

				// report remaining secondary alignments, up-to MAX_REPORT
				for(int i = 0; i < alnPEList.size() && (MAX_REPORT == 0 || i < MAX_REPORT); i++) {
					SAMRecordPair repPair = alnPEList.get(i);
					if(doUpdateBit)
						repPair.setNotPrimaryAlignmentFlags(i != 0);
					int nReport = MAX_REPORT == 0 ? Math.min(alnPEList.size(), MAX_REPORT) : alnPEList.size();
					int nFiltered = alnPEList.size();
					if(repPair.fwdRecord != null) {
						repPair.fwdRecord.setAttribute("NH", nReport);
						repPair.fwdRecord.setAttribute("XN", nFiltered);
						out.addAlignment(repPair.fwdRecord);
					}
					if(repPair.revRecord != null) {
						repPair.revRecord.setAttribute("NH", nReport);
						repPair.revRecord.setAttribute("XN", nFiltered);
						out.addAlignment(repPair.revRecord);
					}
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
		
		/** Get PE inferredInsertSize (template length)
		 * @return PE inferred insert size as the distance between leftmost start and rightmost end
		 */
		public int getPETemplateLen() {
			int tLen = fwdRecord != null ? fwdRecord.getInferredInsertSize() : revRecord.getInferredInsertSize();
			return Math.abs(tLen);
		}
		
		/**
		 * Get PE log-likelihood
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
		 * Get the PE paring probability according to learned Gaussion distribution
		 * @return probability density that this pair is in given inferred template length
		 */
		public double getPEPairPr() {
			return normModel.density(getPETemplateLen());
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
		 * Get PE postP, which is the same for both forward and reverse read
		 * @return  postP either from fwdRecord or revRecord, which one is not null
		 */
		public double getPEPostP() {
			return fwdRecord != null ? Double.parseDouble(fwdRecord.getStringAttribute("XP")) : Double.parseDouble(revRecord.getStringAttribute("XP")); 
		}
		
		/**
		 * Set postP to an AlignRecordPair
		 * @param postP  postP to be set to both pair components
		 */
		public void setPEPostP(double postP) {
			if(fwdRecord != null)
				fwdRecord.setAttribute("XP", Double.toString(postP));
			if(revRecord != null)
				revRecord.setAttribute("XP", Double.toString(postP));
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
		 * @return  the difference between the mapQ value, ties are broken by PEinsert length
		 */
		@Override
		public int compareTo(SAMRecordPair that) {
			return Double.compare(getPEPostP(), that.getPEPostP());
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
				"Options:    -r/--min-align-rate DOUBLE  minimum fraction of align length relative to the read length [" + MIN_ALIGN_RATE + "]" + newLine +
				"            --seed-len INT  seed length for Burrows-Wheeler algorithm dependent aligners [" + SAMAlignFixer.SEED_LEN + "]" + newLine +
				"            --seed-mis DOUBLE  %mismatches allowed in seed region [" + MAX_SEED_MIS + "]" + newLine +
				"            --seed-indel DOUBLE  %indels allowed in seed region [" + MAX_SEED_INDEL + "]" + newLine +
				"            --all-mis  DOUBLE  %mismatches allowed in the entire insert region (excluding clipped/intron regions) [" + MAX_ALL_MIS + "]" + newLine +
				"            --all-indel DOUBLE  %in-dels allowed in the entire insert region [" + MAX_ALL_INDEL + "]" + newLine +
				"            -i/--identity DOUBLE  mimimum %identity allowd for the alignment as 100 - (%mismatches+%in-dels) [" + MIN_IDENTITY + "]" + newLine +
				"            --no-estimate FLAG  do not try to estimate the fragment size distribution; paring probability is ignored" + newLine +
				"            --min-frag-len DOUBLE  estimated minimum fragment (insert) length [" + MIN_FRAG_LEN + "]" + newLine +
				"            --max-frag-len DOUBLE  estimated maximum fragment (insert) length [" + MAX_FRAG_LEN + "]" + newLine +
				"            --max-estimate-scan INT  maximum alignment records to use for estimate fragment distribution [" + MAX_ESTIMATE_SCAN + "]" + newLine +
				"            --clip-handle STRING  how to treat soft/hard-clipped bases as mismathes, USE for use all, IGNORE for ignore, END5 for only use 5' clipped, END3 for only use 3' clipped [" + SAMAlignFixer.CLIP_MODE + "]" + newLine +
				"            --1DP FLAG  enable 1-dimentional dymamic programming insert re-assesment, useful for non-local aligners, i.e. bowtie" + newLine +
				"            --1DP-gap-open-penalty INT  gap open penalty for 1DP [" + SAMAlignFixer.GAP_OPEN_PENALTY_1DP + "]" + newLine +
				"            --1DP-gap-ext-penalty INT  gap extension penalty for 1DP [" + SAMAlignFixer.GAP_EXT_PENALTY_1DP + "]" + newLine +
				"            --match-score INT  match score for 1DP and calculating mapQ [" + SAMAlignFixer.MATCH_SCORE + "]" + newLine +
				"            --mis-score INT  mismatch score for 1DP and calculating mapQ [" + SAMAlignFixer.MIS_SCORE + "]" + newLine +
				"            --gap-open-penalty INT  gap open penalty for calculating mapQ [" + SAMAlignFixer.GAP_OPEN_PENALTY + "]" + newLine +
				"            --gap-ext-penalty INT  gap extension penalty for calculating mapQ [" + SAMAlignFixer.GAP_EXT_PENALTY + "]" + newLine +
				"            -rindel/--relative-indel-penalty FLAG  use relative indel penalty instead of absolute penalty" + newLine +
				"            --clip-penalty INT  additional penalty for soft or hard clipped bases for calculating mapQ [" + SAMAlignFixer.CLIP_PENALTY + "]" + newLine +
				"            --ignore-clip-penalty FLAG  ignore clip penalties completley, good for RNA-seq alignment with DNA-seq aligners" + newLine +
				"            --known-SNP-penalty INT  known SNP penalty for calculating mapQ [" + SAMAlignFixer.KNOWN_SNP_PENALTY + "]" + newLine +
				"            --known-INDEL-penalty INT  known IN-DEL penalty for calculating mapQ [" + SAMAlignFixer.KNOWN_INDEL_PENALTY + "]" + newLine +
				"            --known-MULTISUBSTITUTION-penalty INT  known large/multi-substitution penalty for calculating mapQ [" + SAMAlignFixer.KNOWN_MULTISUBSTITUTION_PENALTY + "]" + newLine +
				"            --out-SAM FLAG  write SAM text output instead of BAM binary output" + newLine +
				"            --silent FLAG  ignore certain SAM format errors such as empty reads" + newLine +
				"            -N/--max-hit INT  max-hit value used during the mapping step, 0 for no limit [" + MAX_HIT + "]" + newLine +
				"            --min-mapQ INT  min mapQ calculated with Bayesian method [" + MIN_MAPQ + "]" + newLine +
				"            --max-best INT  max allowed best-stratum hits to report for a given read, 0 for no limit [" + MAX_BEST + "]" + newLine +
				"            --max-report INT  max report valid hits determined by --min-mapQ and --max-best, 0 for no limit [" + MAX_REPORT + "]" + newLine +
				"            --no-update-bit FLAG  do not update the secondary alignment bit flag (0x100) after filtering" + newLine +
				"            --best-only FLAG  only report unique best hit, equivelant to --max-best 1 --max-report 1" + newLine +
				"            --best FLAG  report the best hit, ignore any secondary hit, equivelant to --max-best 0 --max-report 1" + newLine +
				"            --max-sensitivity FLAG  maximaze sensitivity by ignoring the mismatch, indel options" + newLine +
				"            --sort-method STRING  sorting method for output SAM/BAM file, must be \"none\", \"name\" or \"coordinate\" [none]" + newLine +
				"            --chrom-list FILE  pre-filtering file containing one chromosome name per-line" + newLine +
				"            --known-SNP FILE  known SNP file in vcf/gvcf format (v4.0+, .gz supported), used for calculating mapQ" + newLine +
				"            --AF-tag STRING  Allele Frequency Tag in VCF file to check/use for determining penaltyScores for known SNPs, use NULL to disable [AF]" + newLine +
				"            --fix-MD FLAG  try to fix the MD:Z string format for certain NGS aligners that generate invalid tags" + newLine +
				"            -v FLAG  show verbose information"
				);
	}

	private static void parseOptions(String[] args) throws IllegalArgumentException {
		for(int i = 0; i < args.length; i++)
			if(args[i].equals("-in"))
				inFile = args[++i];
			else if(args[i].equals("-out"))
				outFile = args[++i];
			else if(args[i].equals("-r") || args[i].equals("--min-align-rate"))
				MIN_ALIGN_RATE = Double.parseDouble(args[++i]);
			else if(args[i].equals("--seed-len"))
				SAMAlignFixer.setSEED_LEN(Integer.parseInt(args[++i]));
			else if(args[i].equals("--seed-mis"))
				MAX_SEED_MIS = Double.parseDouble(args[++i]);
			else if(args[i].equals("--seed-indel"))
				MAX_SEED_INDEL = Double.parseDouble(args[++i]);
			else if(args[i].equals("--all-mis"))
				MAX_ALL_MIS = Double.parseDouble(args[++i]);
			else if(args[i].equals("--all-indel"))
				MAX_ALL_INDEL = Double.parseDouble(args[++i]);
			else if(args[i].equals("-i") || args[i].equals("--identity"))
				MIN_IDENTITY = Double.parseDouble(args[++i]);
			else if(args[i].equals("--no-estimate"))
				NO_ESTIMATE = true;
			else if(args[i].equals("--min-frag-len")) {
				MIN_FRAG_LEN = Double.parseDouble(args[++i]);
				if(!(MIN_FRAG_LEN >= 0))
					throw new IllegalArgumentException("--min-frag-len must be non-negative");
			}
			else if(args[i].equals("--max-frag-len")) {
				MAX_FRAG_LEN = Double.parseDouble(args[++i]);
				if(!(MAX_FRAG_LEN > 0))
					throw new IllegalArgumentException("--max-frag-len must be positive");
			}
			else if(args[i].equals("--max-estimate-scan")) {
				MAX_ESTIMATE_SCAN = Long.parseLong(args[++i]);
			}
			else if(args[i].equals("--clip-handle"))
				SAMAlignFixer.CLIP_MODE = SAMAlignFixer.ClipHandlingMode.valueOf(args[++i]);
			else if(args[i].equals("--1DP"))
				DO_1DP = true;
			else if(args[i].equals("--1DP-gap-open-penalty"))
				SAMAlignFixer.setGAP_OPEN_PENALTY_1DP(Integer.parseInt(args[++i]));
			else if(args[i].equals("--1DP-gap-ext-penalty"))
				SAMAlignFixer.setGAP_EXT_PENALTY_1DP(Integer.parseInt(args[++i]));
			else if(args[i].equals("--match-score"))
				SAMAlignFixer.setMATCH_SCORE(Integer.parseInt(args[++i]));
			else if(args[i].equals("--mis-score"))
				SAMAlignFixer.setMIS_SCORE(Integer.parseInt(args[++i]));
			else if(args[i].equals("--gap-open-penalty"))
				SAMAlignFixer.setGAP_OPEN_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--gap-ext-penalty"))
				SAMAlignFixer.setGAP_EXT_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("-rindel") || args[i].equals("--relative-indel-penalty"))
				SAMAlignFixer.INDEL_MODE = SAMAlignFixer.IndelPenaltyMode.RELATIVE;
			else if(args[i].equals("--clip-penalty"))
				SAMAlignFixer.setCLIP_PENALTY(Integer.parseInt(args[++i]));
			else if(args[i].equals("--ignore-clip-penalty"))
				SAMAlignFixer.setIGNORE_CLIP_PENALTY(true);
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
			else if(args[i].equals("-N") || args[i].equals("--max-hit"))
				MAX_HIT = Integer.parseInt(args[++i]);
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
			else if(args[i].equals("--max-sensitivity")) {
				MAX_SENSITIVITY = true;
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
			else if(args[i].equals("--fix-MD"))
				fixMD = true;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'");
		// Check required options
		if(inFile == null)
			throw new IllegalArgumentException("-in must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-out must be specified");
		// Check other options
		if(MIN_ALIGN_RATE < 0 || MIN_ALIGN_RATE > 1)
			throw new IllegalArgumentException("-r/--min-align-rate must be between 0 and 1");
		if(MAX_SEED_MIS < 0 || MAX_SEED_MIS > 100)
			throw new IllegalArgumentException("--seed-mis must be between 0 and 100");
		if(MAX_ALL_MIS < 0 || MAX_ALL_MIS > 100)
			throw new IllegalArgumentException("--all-mis must be between 0 and 100");
		if(MAX_SEED_INDEL < 0 || MAX_SEED_INDEL > 100)
			throw new IllegalArgumentException("--seed-indel must be between 0 and 100");
		if(MAX_ALL_INDEL < 0 || MAX_ALL_INDEL > 100)
			throw new IllegalArgumentException("--all-indel must be between 0 and 100");
		if(!(MIN_IDENTITY >= 0 && MIN_IDENTITY <= 100))
			throw new IllegalArgumentException("-i/--identity must be between 0 and 100");
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

	/** get align insert rate relative to the read length from AlignerBoost internal tag
	 * @return the identity if tag "XL" exists
	 * throws {@RuntimeException} if tag "XL" not exists
	 */
	static double getSAMRecordInsertRate(SAMRecord record) throws RuntimeException {
		return (double) record.getIntegerAttribute("XL") / record.getReadLength();
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
	private static double[] calcPEHitPostP(List<SAMRecordPair> alnPEList, int totalPair, int maxPair) {
		if(alnPEList == null) // return null for null list
			return null;
		if(alnPEList.isEmpty())
			return new double[0]; // return empty array for empty list
		
		int nPairs = alnPEList.size();
		// get un-normalized posterior probs
		double[] postP = new double[nPairs];
		if(nPairs == 1) {
			postP[0] = UNIQ_MAPQ;
			alnPEList.get(0).setPEMapQ(UNIQ_MAPQ);
			return postP;
		}
		
		for(int i = 0; i < nPairs; i++) {
			// get postP as priorP * likelihood, with prior proportional to the alignLength
			postP[i] = alnPEList.get(i).getPEInsertLen() * Math.pow(10.0,  alnPEList.get(i).getPEAlignLik());
			if(!NO_ESTIMATE) // pairing probability needs to be considered
				postP[i] *= alnPEList.get(i).getPEPairPr();
		}
		
		// normalize postP
		Stats.normalizePostP(postP, maxPair == 0 || totalPair < maxPair ? 0 : Math.sqrt(maxPair));
		// reset the mapQ values
		for(int i = 0; i < nPairs; i++) {
			alnPEList.get(i).setPEPostP(postP[i]);
			double mapQ = Stats.phredP2Q(1 - postP[i]);
			if(Double.isNaN(mapQ)) // is NaN
				alnPEList.get(i).setPEMapQ(INVALID_MAPQ);
			else
				alnPEList.get(i).setPEMapQ(mapQ > MAX_MAPQ ? MAX_MAPQ : (int) Math.round(mapQ));
		}
		return postP;
	}

	private static int filterPEHits(List<SAMRecordPair> alnPEList,
			double maxSeedMis, double maxSeedIndel, double maxAllMis, double maxAllIndel) {
		int n = alnPEList.size();
		int removed = 0;
		for(int i = n - 1; i >= 0; i--) { // search backward for maximum performance
			SAMRecordPair pair = alnPEList.get(i);
			if(!(  (pair.fwdRecord == null || getSAMRecordPercentSeedMis(pair.fwdRecord) <= maxSeedMis
					&& getSAMRecordPercentSeedIndel(pair.fwdRecord) <= maxSeedIndel
					&& getSAMRecordPercentAllMis(pair.fwdRecord) <= maxAllMis
					&& getSAMRecordPercentAllIndel(pair.fwdRecord) <= maxAllIndel)
				&& (pair.revRecord == null || getSAMRecordPercentSeedMis(pair.revRecord) <= maxSeedMis
					&& getSAMRecordPercentSeedIndel(pair.revRecord) <= maxSeedIndel
					&& getSAMRecordPercentAllMis(pair.revRecord) <= maxAllMis
					&& getSAMRecordPercentAllIndel(pair.revRecord) <= maxAllIndel)
				)) {
				alnPEList.remove(i);
				removed++;
			}
		}
		return removed;
	}
	
	private static int filterPEHits(List<SAMRecordPair> alnPEList, double minAlignRate, double minIdentity) {
		int n = alnPEList.size();
		int removed = 0;
		for(int i = n - 1; i >= 0; i--) { // search backward for maximum performance
			SAMRecordPair pair = alnPEList.get(i);
			if(!(  (pair.fwdRecord == null || getSAMRecordInsertRate(pair.fwdRecord)>= minAlignRate 
					&& getSAMRecordIdentity(pair.fwdRecord) >= minIdentity)
				&& (pair.revRecord == null || getSAMRecordInsertRate(pair.revRecord)>= minAlignRate 
					&& getSAMRecordIdentity(pair.revRecord)>= minIdentity)
				)) {
				alnPEList.remove(i);
				removed++;
			}
		}
		return removed;
	}

	private static int filterPEHits(List<SAMRecordPair> alnPEList, int minMapQ) {
		int n = alnPEList.size();
		int removed = 0;
		for(int i = n - 1; i >= 0; i--) { // search backward for maximum performance
			if(alnPEList.get(i).getPEMapQ() < minMapQ) {
				alnPEList.remove(i);
				removed++;
			}
		}
		return removed;
	}
	
	private static final int INVALID_MAPQ = 255;
	private static final int MAX_MAPQ = 200; // MAX meaniful mapQ value, if not 255
	private static final int UNIQ_MAPQ = 250;
	private static String inFile;
	private static String outFile;
	private static String chrFile;
	private static String knownSnpFile;
	// filter options
	private static double MIN_ALIGN_RATE = 0.9;
	private static double MAX_SEED_MIS = 4; // max % seed mismatch
	private static double MAX_SEED_INDEL = 0; // max % seed indel
	private static double MAX_ALL_MIS = 6; // max % all mismatch
	private static double MAX_ALL_INDEL = 0; // max % all indel
	private static double MIN_IDENTITY = 0; // min identity
	private static boolean MAX_SENSITIVITY; // enable max-sensitivity?
	private static boolean DO_1DP;
	private static boolean isSilent; // ignore SAM warnings?
	private static boolean noMix; // do not allow unpaired alignments for paired reads?
	// mate pair options
	private static boolean NO_ESTIMATE;
	private static double MIN_FRAG_LEN = 50;
	private static double MAX_FRAG_LEN = 750;
	private static double MEAN_FRAG_LEN;
	private static double SD_FRAG_LEN;
	private static final long MIN_ESTIMATE_BASE = 1000; // MIN alignment number to make an accurate estimate
	private static long MAX_ESTIMATE_SCAN;
	// best stratum options
	private static int MAX_HIT = 10; // MAX_HIT used during the mapping step
	private static int MIN_MAPQ = 10;
	private static int MAX_BEST = 1;
	private static int MAX_REPORT = 1;
	private static boolean doUpdateBit = true;
	private static int verbose; // verbose level
	private static boolean fixMD = false;
	private static Set<String> chrFilter;
	private static VCFFileReader knownVCF;
	// general options
	private static GroupOrder groupOrder = GroupOrder.none;
	private static SortOrder sortOrder = SortOrder.unsorted;
	private static boolean OUT_IS_SAM; // outFile is SAM format?
	private static Timer processMonitor;
	private static ProcessStatusTask statusTask;
	private static final int statusFreq = 10000;
	// Math related static members
	private static NormalDistribution normModel;
}
