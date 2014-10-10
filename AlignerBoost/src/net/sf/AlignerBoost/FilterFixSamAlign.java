package net.sf.AlignerBoost;
import java.io.*;
import java.util.regex.*;

import htsjdk.samtools.*;
import static net.sf.AlignerBoost.EnvConstants.progFile;

/** Parse SAM/BAM alignment file with Picard java packages in CLASSPATH
 * A 1-dimensional DP algorithm (1DP) is implemented to re-estimate the insert length of the alignment
 * The calculated metrix are stored in application-specified SAM tags as below:
 * Tag  Type  Description
 * XL   i     # alignment (insert) length, including M,=,X,I,D but not S,H,P,N, determined by Cigar or 1DP
 * XF   i     # actual insert from (start) relative to read
 * XS   i     # of seed mismatches
 * XA   i     # of all mismatches
 * XG   i     # of all indels, Ns are not included
 * XI   f     # alignment identity 1 - (XA + XG) / XL
 * YS   f     # alignment score from AlignerBoost
 * YQ   f     # quality weighted alignment score from AlignerBoost
 * @author Qi Zheng
 * @version 1.1
 */
public class FilterFixSamAlign {
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

		SamReaderFactory readerFac = SamReaderFactory.makeDefault();
		SAMFileWriterFactory writerFac = new SAMFileWriterFactory();
		if(!isSilent)
			readerFac.validationStringency(ValidationStringency.LENIENT); // use LENIENT strigency
		else
			readerFac.validationStringency(ValidationStringency.SILENT); // use SILENT strigency

		SamReader in = readerFac.open(new File(inFile));
		SAMFileWriter out = OUT_IS_SAM ? writerFac.makeSAMWriter(in.getFileHeader(), false, new File(outFile)) : writerFac.makeBAMWriter(in.getFileHeader(), false, new File(outFile));

		// check each alignment
		for(SAMRecord record : in) {
			int readLen = record.getReadLength();
			if(record.getReferenceIndex() == -1 || readLen == 0) // non mapped read or 0-length read
				continue;
			int cigarLen = record.getCigarLength();
			Cigar cigar = record.getCigar();
			// valid cigar first
			assert isValidCigarLength(cigar, readLen);

			// calcuate alignment length
			int alnLen = calcAlnLenByCigar(cigar);
			// get align status index
			char[] status = getAlnStatusBySAMRecord(record, alnLen);
			// get insert from relative to read, ignore soft-clipped region
			int insertFrom = calcInsertFromByCigar(cigarLen, cigar, record.getReadNegativeStrandFlag()); // 0-based
			// calculate insert length, either as align length or by 1DP
			int insertLen = 0;
			if(!DO_1DP)
				insertLen = calcInsertLenByCigar(cigar);
			else {
				insertLen = calcInsertLenBy1DP(status, insertFrom);
				if(!NO_FIX) // fx SAM is requested
					fixSAMRecordCigarTag(record, alnLen, status, insertFrom, insertLen);
			}

			// calculate nmismatches and indels
			int nSeedMis = 0;
			int nSeedIndel = 0;
			int nAllMis = 0;
			int nAllIndel = 0;
			for(int i = 0; i < insertLen; i++) {
				if(status[i] == 'X') { // mismatch
					if(i < SEED_LEN)
						nSeedMis++;
					nAllMis++;
				}
				else if(status[i] == 'I' || status[i] == 'D') { // indel 
					if(i < SEED_LEN)
						nSeedIndel++;
					nAllIndel++;
				}
				else
					continue;
			}
			// calculate identity
			// output
			if(insertLen >= MIN_INSERT &&
					nSeedMis <= MAX_SEED_MIS * SEED_LEN / 100 && nSeedIndel <= MAX_SEED_INDEL * SEED_LEN / 100 &&
					nAllMis <= MAX_ALL_MIS * insertLen / 100 && nAllIndel <= MAX_ALL_INDEL * insertLen / 100) {
				// add additional tags
				float identity = 1 - (nAllMis + nAllIndel) / ((float) insertLen); 
				record.setAttribute("XL", insertLen);
				record.setAttribute("XF", insertFrom);
				record.setAttribute("XS", nSeedMis);
				record.setAttribute("XA", nAllMis);
				record.setAttribute("XG", nAllIndel);
				record.setAttribute("XI", identity);
				//record.setAttribute("XS", calcAlnScore(status));
				//record.setAttribute("XQ", calcAlnScore(status, record.getBaseQualities()));
				//System.out.println("strand: " + (record.getReadNegativeStrandFlag() ? "-" : "+") + " " + record.getCigarString() + "  " + record.getStringAttribute("MD") + " insertFrom: " + insertFrom + " insertLen: " + insertLen);
				out.addAlignment(record);
			}
		} // end each record
		try {
			in.close();
			out.close();
		}
		catch(IOException e) {
			System.err.println(e.getMessage());
		}
	}

	private static String inFile;
	private static String outFile;
	private static int MIN_INSERT = 15;
	private static int SEED_LEN = 25;
	private static double MAX_SEED_MIS = 4;
	private static final double MAX_SEED_INDEL = 0; // seed indel is always not allowed
	private static double MAX_ALL_MIS = 6;
	private static double MAX_ALL_INDEL = 0;
	private static boolean DO_1DP;
	private static int MATCH_SCORE = 1;
	private static int MIS_SCORE = -2;
	private static int GAP_OPEN_PENALTY = 4;
	private static int GAP_EXT_PENALTY = 1;
	private static int REF_QUAL = 40; // reference quality for deletions
	private static boolean NO_FIX; // do not fix the SAM file if 1DP requested ?
	private static boolean OUT_IS_SAM; // outFile is SAM format?
	private static boolean isSilent; // ignore SAM warnings?

	private static Pattern misPat1 = Pattern.compile("(\\d+)(.*)");
	private static Pattern misPat2 = Pattern.compile("([A-Z]|\\^[A-Z]+)(\\d+)");
	private static Pattern misPat3 = Pattern.compile("\\d+|[A-Z]|\\^[A-Z]+");

	private static void printUsage() { 
		System.err.println("Usage:   java -jar " + progFile + " run filter " +
				"<-in SAM|BAM-INPUT> <-out SAM|BAM-OUTPUT> [options]" +
				"Options:    --min-insert  minimum insert length (excluding adapters) of a read to allow amabiguous alignment, default 15" +
				"            --seed-len seed length for Burrows-Wheeler algorithm dependent aligners, default 25" +
				"            --seed-mis %mismatches allowed in seed region, default 4" +
				"            --all-mis %mismatches allowed in the entire insert region (excluding masked regions and Ns), default 6" +
				"            --all-indel %in-dels allowed in the entire insert region, default 2" +
				"            --1DP enable 1-dimentional dymamic programming (1DP) re-aligning, useful for non-local aligners, i.e. bowtie" +
				"            --match-score match score for 1DP, default 1" +
				"            --mis-score mismatch score for 1DP, default -2" +
				"            --gap-open-penalty gap open penalty for 1DP, default 4" +
				"            --gap-ext-penalty gap extension penalty, default 1" +
				"            --no-fix AlignerBoost by default tries to fix the MD:z tag to include correct mismatch information, set it to disable this behaviour"+
				"            --out-SAM write SAM text format output instead of BAM binary output" +
				"            --silent ignore certain SAM format errors such as empty reads"
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
				SEED_LEN = Integer.parseInt(args[++i]);
			else if(args[i].equals("--seed-mis"))
				MAX_SEED_MIS = Double.parseDouble(args[++i]);
			else if(args[i].equals("--all-mis"))
				MAX_ALL_MIS = Double.parseDouble(args[++i]);
			else if(args[i].equals("--all-indel"))
				MAX_ALL_INDEL = Double.parseDouble(args[++i]);
			else if(args[i].equals("--1DP"))
				DO_1DP = true;
			else if(args[i].equals("--match-score"))
				MATCH_SCORE = Integer.parseInt(args[++i]);
			else if(args[i].equals("--mis-score"))
				MIS_SCORE = Integer.parseInt(args[++i]);
			else if(args[i].equals("--gap-open-penalty"))
				GAP_OPEN_PENALTY = Integer.parseInt(args[++i]);
			else if(args[i].equals("--gap-ext-penalty"))
				GAP_EXT_PENALTY = Integer.parseInt(args[++i]);
			else if(args[i].equals("--no-fix"))
				NO_FIX = true;
			else if(args[i].equals("--out-SAM"))
				OUT_IS_SAM = true;
			else if(args[i].equals("--silent"))
				isSilent = true;
			else
				throw new IllegalArgumentException("Unknown option '" + args[i] + "'");
		// Check required options
		if(inFile == null)
			throw new IllegalArgumentException("-in must be specified");
		if(outFile == null)
			throw new IllegalArgumentException("-out must be specified");;
		// Check other options
		if(SEED_LEN <= 0)
			throw new IllegalArgumentException("--seed-len must be positive");
		if(MAX_SEED_MIS < 0 || MAX_SEED_MIS > 100)
			throw new IllegalArgumentException("--seed-mis must be between 0 to 100");
		if(MAX_ALL_MIS < 0 || MAX_ALL_MIS > 100)
			throw new IllegalArgumentException("--all-mis must be between 0 to 100");
		if(MAX_ALL_INDEL < 0 || MAX_ALL_INDEL > 100)
			throw new IllegalArgumentException("--all-indel must be between 0 to 100");
		if(MATCH_SCORE <= 0)
			throw new IllegalArgumentException("--match-score must be positive");
		if(MIS_SCORE > 0)
			throw new IllegalArgumentException("--mis-score must be non-positive");
		if(GAP_OPEN_PENALTY < 0)
			throw new IllegalArgumentException("--gap-open-penalty must be non negative");
		if(GAP_EXT_PENALTY <= 0)
			throw new IllegalArgumentException("--gap-ext-penalty must be positive");
		if(OUT_IS_SAM && outFile.endsWith(".bam"))
			System.err.println("Warning: output file '" + outFile + "' might not be SAM format");
	}

	/** Valid if a cigar length match the read length
	 * @return true if cigar is valid for its length
	 */
	static boolean isValidCigarLength(Cigar cigar, int readLen) {
		int length = 0;
		for(CigarElement cigEle : cigar.getCigarElements()) {
			CigarOperator cigOp = cigEle.getOperator();
			if(cigOp == CigarOperator.M || cigOp == CigarOperator.I || cigOp == CigarOperator.EQ ||
					cigOp == CigarOperator.X || cigOp == CigarOperator.S)
				length += cigEle.getLength();
		}
		return length == readLen;
	}

	/** calculate alignment length using Cigar
	 * @return actual alignment length including M,=,X,I,D,S but not H,P and N
	 */
	static int calcAlnLenByCigar(Cigar cigar) {
		int alnLen = 0;
		for(CigarElement cigEle : cigar.getCigarElements()) {
			switch(cigEle.getOperator()) {
			case M: case EQ: case X: case I: case D: case S:
				alnLen += cigEle.getLength();
				break;
			default: // H,P or N
				break; // do nothing
			}
		}
		return alnLen;
	}

	/** calculate insert length by alnLen and Cigar, ignore soft-clipped bases
	 * @return insert length without soft-clip
	 */
	static int calcInsertLenByCigar(Cigar cigar) {
		if(cigar == null || cigar.numCigarElements() == 0)
			return 0;
		int insertLen = 0;
		for(CigarElement cigEle : cigar.getCigarElements()) {
			switch(cigEle.getOperator()) {
			case M: case EQ: case X: case I: case D:
				insertLen += cigEle.getLength();
				break;
			default: // S,H,P or N
				break; // do nothing
			}
		}
		return insertLen;
	}

	/** calculate relative read from position by Cigar
	 * @return 0-based from excluding soft-clipped bases
	 */
	static int calcInsertFromByCigar(int cigarLen, Cigar cigar, boolean isMinus) {
		if(cigar == null || cigar.numCigarElements() == 0)
			return 0;
		CigarElement firstCig = cigar.getCigarElement(0);
		CigarElement lastCig = cigar.getCigarElement(cigarLen - 1);
		if(!isMinus && firstCig.getOperator() == CigarOperator.S)
			return firstCig.getLength();
		else if(isMinus && lastCig.getOperator() == CigarOperator.S)
			return lastCig.getLength();
		else
			return 0;
	}

	/** get align status ('M', '=', 'X', 'I', 'S') given cigar and mismatch tag, if exists
	 * @return a char array index with length = alnLen, with 0 as dummy position
	 */
	static char[] getAlnStatusBySAMRecord(SAMRecord record, int alnLen) {
		Cigar cigar = record.getCigar();
		String misStr = record.getStringAttribute("MD");
		boolean isMinus = record.getReadNegativeStrandFlag();

		// init align status index
		char[] status = new char[alnLen];
		int shift = 0;
		// build status from cigar first
		for(CigarElement cigEle : cigar.getCigarElements()) {
			int cigLen = cigEle.getLength();
			CigarOperator cigOp = cigEle.getOperator();
			switch(cigOp) {
			case M: case EQ: case X: case I: case D: case S:
				for(int i = 0; i < cigLen; i++)
					status[shift++] = (char) CigarOperator.enumToCharacter(cigOp);
				break;
			default: // H, N or P, do nothing
				break;
			}
		}
		// check misStr, if exists
		if(misStr != null) {
			// set the pos to the first non 'S'
			int pos = 0;
			while(status[pos] == 'S')
				pos++;
			Matcher match1 = misPat1.matcher(misStr);
			match1.find();
			pos = advanceReferencePos(status, pos, Integer.parseInt(match1.group(1)));
			String other = match1.group(2);
			Matcher match2 = misPat2.matcher(other); // parse other part
			while(match2.find()) {
				String misSeq = match2.group(1);
				int followLen = Integer.parseInt(match2.group(2));
				if(misSeq.length() == 1) { // a mismatch
					/*if(!(status[pos] == '=' || status[pos] == 'M')) {
						System.out.println(pos + ":" + status[pos] + " " + match2.group());
						System.out.println(record.getSAMString());
						System.out.println(inFile);
						System.exit(-1);
					}*/
					assert status[pos] == '=' || status[pos] == 'M';
					status[pos] = 'X';
					pos = advanceReferencePos(status, pos);
				}
				else // a deletion
					pos += misSeq.length() - 1;
				// advance the pos by followLen relative to read
				pos = advanceReferencePos(status, pos, followLen);
				//pos += followLen;
			}
			assert pos <= alnLen;
		}
		// reverse status index, if on minus strand, so it is always in read orientation
		if(isMinus)
			reverse(status);
		return status;
	}

	
	/** reverse an char array
	 */
	private static void reverse(char[] arr) {
		int len = arr.length;
		for(int i = 0; i < len / 2; i++) {
			char tmp = arr[i];
			arr[i] = arr[len - 1 - i];
			arr[len - 1 - i] = tmp;
		}
	}

	/** private method to advance reference pos to a given step, ignore I (insertion) that is not present on reference
	 */
	private static int advanceReferencePos(char[] status, int oldPos, int step) {
		int shift = 0;
		int pos = oldPos;
		while(pos < status.length) {
			if(shift == step && status[pos] != 'I')
				break;
			if(status[pos++] != 'I') // insertion not present on reference
				shift++;
		}
		return pos;
	}

	private static int advanceReferencePos(char[] status, int oldPos) {
		return advanceReferencePos(status, oldPos, 1);
	}

	/** calculate actual insert size by 1-D dynamic programming (1DP)
	 * @param status align status index always in read orientation
	 * @return actual insert length
	 */
	static int calcInsertLenBy1DP(char[] status, int insertFrom) {
		// do 1DP to determine insertLen
		int alnLen = status.length;
		int[] dpScore = new int[alnLen + 1]; // Position 0 is dummy
		int maxScore = Integer.MIN_VALUE;
		int insertLen = 0;
		// dynamic-programming to get dp-score, highest score and insert-length simultaneously
		for(int i = insertFrom; i < alnLen; i++) {
			int score = 0;
			switch(status[i]) {
			case 'M': case '=':
				score = MATCH_SCORE;
				break;
			case 'X':
				score = MIS_SCORE;
				break;
			case 'I': case 'D': // indel
				score = i == 0 || status[i - 1] != 'I' && status[i - 1] != 'D' ? GAP_OPEN_PENALTY : GAP_EXT_PENALTY;
				break;
			default: // do nothing
				break;
			}
			dpScore[i + 1] = dpScore[i] + score; // dpScore is 1-based
			if(dpScore[i + 1] > maxScore) {
				insertLen = i + 1;   // Update insertLen
				maxScore = dpScore[i + 1];  // Update maxScore
			}
		}
		return insertLen - insertFrom;
	}

	/** Calc refAlnLength
	 * @param insertLen 
	 * @param from 
	 */
	static int calcReferenceAlnLen(char[] status, int from, int alnLen) {
		int refAlnLen = 0;
		for(int i = from; i < status.length; i++)
			if(status[i] != 'S' && status[i] == 'I')
				refAlnLen++;
		return refAlnLen;
	}

	/** Calc readInsertLength by subtracting deletions from the insertLen
	 */
	static int calcReferenceInsertLen(char[] status, int insertFrom, int insertLen) {
		assert insertFrom + insertLen <= status.length;
		// status is in read orientation
		int refInsertLen = insertLen;
		for(int i = insertFrom; i < insertFrom + insertLen; i++)
			if(status[i] == 'I')
				refInsertLen--;
		return refInsertLen;
	}

	/** Fix SAMRecord Cigar and mismatch tag (MD:Z), given the insert from and insert length
	 *
	 */
	static void fixSAMRecordCigarTag(SAMRecord record, int alnLen, char[] status, int from, int insertLen) {
		if(from + insertLen == alnLen) // no fix needed
			return;
		boolean isMinus = record.getReadNegativeStrandFlag();
		int readLen = record.getReadLength();
		int refAlnLen = calcReferenceAlnLen(status, from, alnLen);
		int refInsertLen = calcReferenceInsertLen(status, from, insertLen);
		
		int to = from + insertLen; // 1-based
		int refClipLen = refAlnLen - refInsertLen; // clipped refInsertLen by 1DP
		int clipLen = alnLen - to; // clipped insertLen by 1DP
		// calculate readClipLen, ignore deletion
		int readClipLen = clipLen;
		for(int i = to; i < status.length; i++)
			if(status[i] == 'D')
				readClipLen--;

		// fix Cigar
		Cigar oldCig = record.getCigar();
		Cigar newCig = new Cigar();
		int pos = 0; // relative pos to reference
		// note for Cigar soft-clip always happens at the end of the read
		if(isMinus) // - strand, add clip at the begining
			newCig.add(new CigarElement(readClipLen, CigarOperator.S));
		for(CigarElement cigEle : oldCig.getCigarElements()) {
			CigarOperator cigOp = cigEle.getOperator();
			int cigLen = cigEle.getLength();
			if(isMinus) { // - strand, clip at the beginning
				switch(cigOp) {
				case M: case EQ: case X: case I: case D:
					if(pos + cigLen <= clipLen) // competely in soft-clipped region
						;
					else if(pos < clipLen && pos + cigLen > clipLen) // partially inside
						newCig.add(new CigarElement(pos + cigLen - clipLen, cigOp));
					else // not in clipped region
						newCig.add(new CigarElement(cigLen, cigOp)); // make a copy
					pos += cigLen;
					break;
				case S:
					if(pos + cigLen > clipLen) // not in clipped region
						newCig.add(new CigarElement(cigLen, cigOp)); // make a copy
					pos += cigLen; // count but doesn't record
					break;
				default: // ('N', 'H', 'P')
					if(pos >= clipLen) // not clipped
						newCig.add(new CigarElement(cigLen, cigOp));
					break;
				}
			} // end - strand
			else { // + strand, clip at the end
				switch(cigOp) {
				case M: case EQ: case X: case I: case D:
					if(pos + cigLen <= to) // not in clipped region
						newCig.add(new CigarElement(cigLen, cigOp)); // make a copy
					else if(pos < to && pos + cigLen > to) // partially inside
						newCig.add(new CigarElement(to - pos, cigOp));
					else // competely in soft-clipped region
						;
					pos += cigLen;
					break;
				case S:
					if(pos + cigLen <= to) // not in clipped region
						newCig.add(new CigarElement(cigLen, cigOp)); // make a copy
					pos += cigLen; // count but doesn't record
					break;
				default: // 'N', 'H', 'P'
					if(pos < to) // not in clipped region
						newCig.add(new CigarElement(cigLen, cigOp));
					break;
				}
			} // end + strand
		} // end each cigEle
		if(!isMinus) // + strand, add clip at the end
			newCig.add(new CigarElement(readClipLen, CigarOperator.S));
		/*if(!isValidCigarLength(newCig, record.getReadLength())) {
			System.err.println(oldCig);
			System.err.println(newCig);
			System.err.println(isMinus + ":refAlnLen" + refAlnLen + " refInsertLen:" + refInsertLen + " refClipLen:" + refClipLen + " clipLen:" + clipLen + " readClipLen:" + readClipLen + " to:" + to + " readLen:" + record.getReadLength());
			System.err.println(record.getStringAttribute("MD"));
			System.exit(-1);
		}*/
		assert isValidCigarLength(newCig, readLen);

		// fix mismatch string, if exists
		String oldMisStr = record.getStringAttribute("MD");
		// fix refClipLen by removing S
		/*CigarElement firstCig = newCig.getCigarElement(0);
		CigarElement lastCig = newCig.getCigarElement(newCig.numCigarElements() - 1);
		if(isMinus && firstCig.getOperator() == CigarOperator.S)
			refClipLen -= firstCig.getLength();
		else if(!isMinus && lastCig.getOperator() == CigarOperator.S)
			refClipLen -= lastCig.getLength();
		else
		;*/

		if(oldMisStr != null && refClipLen > 0) {
			StringBuilder newMisStr = new StringBuilder(); // use StringBuilder for performance
			Matcher match = misPat3.matcher(oldMisStr);
			pos = 0;
			while(match.find()) {
				String s = match.group();
				if(isMinus) { // - strand, clip at the beginning
					if(isDecimal(s)) { // a span met
						int span = Integer.parseInt(s);
						if(pos >= refClipLen) // not in clipped region
							newMisStr.append(span); // make a copy
						else if(pos < refClipLen && pos + span > refClipLen) // partially inside
							newMisStr.append(pos + span - refClipLen);
						else // competely in soft-clipped region
							; // do nothing
						pos += span;
					}
					else { // a mismatch or deletion
						if(pos >= refClipLen) // not in clip region, note deletion won't partially in clipped region
							newMisStr.append(s);
						if(s.length() == 1) // a mismatch
							pos++;
						else // a deletion
							pos += s.length() - 1;
					}
				}
				else { // + strand, clip at the end
					if(isDecimal(s)) { // a span met
						int span = Integer.parseInt(s);
						if(pos + span <= refInsertLen) // not in clipped region
							newMisStr.append(span); // make a copy
						else if(pos < refInsertLen && pos + span > refInsertLen) // partially inside
							newMisStr.append(pos + span - refInsertLen);
						else // competely in clipped region
							; // do nothing
						pos += span;
					}
					else { // a mismatch or deletion
						if(pos < refInsertLen) // inside
							newMisStr.append(s);
						if(s.length() == 1) // a mismatch
							pos++;
						else // a deletion
							pos += s.length() - 1;
					}
				}
			} // end while find
			/*if(!isMatchedCigarMisStr(newCig, newMisStr.toString())) {
				System.err.println(isMinus + " refClipLen:" + refClipLen);
				System.err.println(" refAlnLen:" + refAlnLen + " refInsertLen:" + refInsertLen + " insertLen:" + insertLen);
				System.err.println("newCig:" + newCig);
				System.err.println("newMisStr:" + newMisStr);
				System.err.println(record.getSAMString());
				System.err.println("inFile:" + inFile);
				System.exit(-1);
			}
			 */
			assert isMatchedCigarMisStr(newCig, newMisStr.toString());
			record.setCigar(newCig); // update the cigar
			record.setAttribute("MD", newMisStr.toString()); // update misStr

		} // end if
	}

	/** validate whether cigar and MD:Z mismatch string matches
	 * @return true if their length matches
	 */
	static boolean isMatchedCigarMisStr(Cigar cigar, String misStr) {
		int cigRefAlnLen = 0;
		int misRefAlnLen = 0;
		for(CigarElement cigEle : cigar.getCigarElements()) {
			CigarOperator cigOp = cigEle.getOperator();
			int cigLen = cigEle.getLength();
			switch(cigOp) {
			case M: case EQ: case X: case D:
				cigRefAlnLen += cigLen;
				break;
			default:
				break;
			}
		}

		Matcher match1 = misPat1.matcher(misStr);
		boolean found = match1.find();
		if(!found && cigar.numCigarElements() != 0)
			return false;
		misRefAlnLen = Integer.parseInt(match1.group(1)); // starting pos in the MD:Z string, can be 0
		String other = match1.group(2);
		Matcher match2 = misPat2.matcher(other); // parse other part
		while(match2.find()) {
			String misSeq = match2.group(1);
			int followLen = Integer.parseInt(match2.group(2));
			if(misSeq.length() == 1) // mismatch
				misRefAlnLen++;
			else // deletion
				misRefAlnLen += misSeq.length() - 1;
			misRefAlnLen += followLen;
		}
		return cigRefAlnLen == misRefAlnLen;
	}

	/** determine whether a string is a decimal integer
	 * @return true if is a valid decimal integer
	 */
	public static boolean isDecimal(String s) {
		if(s.isEmpty())
			return false;
		for(int i = 0; i < s.length(); i++) {
			if(i == 0 && s.charAt(i) == '-') {
				if(s.length() == 1)
					return false;
				else
					continue;
			}
			if(!Character.isDigit(s.charAt(i)))
				return false;
		}
		return true;
	}

	/** calculate un-weighted alignment score
	 * @param status  status index
	 * @return  alignment score
	 */
	private static int calcAlnScore(char[] status) {
		if(status == null)
			return 0;
		int alnScore = 0;
		for(int i = 0; i < status.length; i++) {
			int score = 0;
			switch(status[i]) {
			case 'M': case '=':
				score = MATCH_SCORE;
				break;
			case 'X':
				score += MIS_SCORE;
				break;
			case 'S': case 'H': case 'P': case 'N':
				break;
			case 'I': case 'D':
				score = (status[i-1] != 'I' && status[i-1] != 'D' ? GAP_OPEN_PENALTY : GAP_EXT_PENALTY); 
			default: // S,H,P or N
				break; // do nothing
			}
			alnScore += score;
		}
		return alnScore;
	}

	/** calculate quality weighted alignment score
	 * @param status  status index
	 * @param qual  quality scores in Phred scale
	 * @return  alignment score
	 */
	private static int calcAlnScore(char[] status, byte[] qual) {
		if(status == null)
			return 0;
		assert status.length >= qual.length;
		int alnScore = 0;
		int pos = 0; // relative pos on read
		for(int i = 0; i < status.length; i++) {
			int score = 0;
			int weight = 0;
			switch(status[i]) {
			case 'M': case '=':
				score = MATCH_SCORE;
				weight = qual[pos++];
				break;
			case 'X':
				score += MIS_SCORE;
				weight = qual[pos++];
				break;
			case 'S': case 'H': case 'P': case 'N':
				break;
			case 'I':
				score = (status[i-1] != 'I' ? GAP_OPEN_PENALTY : GAP_EXT_PENALTY);
				weight += qual[pos++];
			case 'D':
				score = (status[i-1] != 'D' ? GAP_OPEN_PENALTY : GAP_EXT_PENALTY);
				weight += REF_QUAL;
			default: // S,H,P or N
				break; // do nothing
			}
			alnScore += score * weight;
		}
		return alnScore;
	}
}
