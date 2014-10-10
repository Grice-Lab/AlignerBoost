/** 
 * a class to provide static method to filter and fix SAM/BAM alignment
 */
package net.sf.AlignerBoost;
import java.util.regex.*;

import htsjdk.samtools.*;

/** Parse SAM/BAM alignment file with Picard java packages in CLASSPATH
 * A 1-dimensional DP algorithm (1DP) is implemented to re-estimate the insert length of the alignment
 * The calculated metrix are stored in application-specified SAM tags,
 * X?: alignment tags, Y? seed tags, Z? total tags
 * Tag  Type  Description
 * XL   i     insert length, including M,=,X,I,D but not S,H,P,N, determined by Cigar or 1DP
 * XF   i     actual insert from (start) relative to read
 * XI   f     alignment identity 1 - (XA + XG) / XL
 * XS   f     alignment score from AlignerBoost
 * XQ   f     quality weighted alignment score from AlignerBoost
 * YL	i     seed length for calculating seed mismatches and indels
 * YX   i     # of seed mismatches
 * YG   i     # of seed indels
 * ZX   i     # of all mismatches
 * ZG   i     # of all indels
 * @author Qi Zheng
 * @version 1.1
 */
public class SAMAlignFixer {
	/**
	 * Fix a SAMRecord by adding the AlignerBoost-specific tags above, and optionally fix the alignment by 1DP 
	 * @param record  SAMRecord alignment to be fixed
	 * @param do1DP  whether do additional 1DP fixing?
	 * @return  false if this SAMRecord does not need to be fixed because it is not-mapped or empty.
	 */
	public static boolean fixSAMRecord(SAMRecord record, boolean do1DP) {
		int readLen = record.getReadLength();
		if(record.getReferenceIndex() == -1 || readLen == 0) // non mapped read or 0-length read
			return false;
		int cigarLen = record.getCigarLength();
		Cigar cigar = record.getCigar();
		// valid cigar first
		assert isValidCigarLength(cigar, readLen);

		// Calculate alignment length
		int alnLen = calcAlnLenByCigar(cigar);
		// get align status index
		char[] status = getAlnStatusBySAMRecord(record, alnLen);
		// get insert from relative to read, ignore soft-clipped region
		int insertFrom = calcInsertFromByCigar(cigarLen, cigar, record.getReadNegativeStrandFlag()); // 0-based
		// calculate insert length, either as align length or by 1DP
		int insertLen = 0;
		if(!do1DP)
			insertLen = calcInsertLenByCigar(cigar);
		else {
			insertLen = calcInsertLenBy1DP(status, insertFrom);
			// alignment need to be fixed after 1DP
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
		// add customized tags
		float identity = 1 - (nAllMis + nAllIndel) / ((float) insertLen);
		// alignment tags
		record.setAttribute("XL", insertLen);
		record.setAttribute("XF", insertFrom);
		record.setAttribute("XI", identity);
		//record.setAttribute("XS", calcAlnScore(status));
		//record.setAttribute("XQ", calcAlnScore(status, record.getBaseQualities()));
		// seed tags
		record.setAttribute("YL", SEED_LEN);
		record.setAttribute("YX", nSeedMis);
		record.setAttribute("YG", nSeedIndel);
		// all tags
		record.setAttribute("ZX", nAllMis);
		record.setAttribute("ZG", nAllIndel);
		return true;
	}
	
	/**
	 * fix SAMRecord without doing 1DP
	 * @param record  SAMRecord to be fixed
	 * @return false if the record does not need to be fixed due to unmapped or empty
	 */
	public static boolean fixSAMRecord(SAMRecord record) {
		return fixSAMRecord(record, false);
	}

	/** test if a cigar length match the read length
	 * @return true if cigar is valid for its length
	 */
	private static boolean isValidCigarLength(Cigar cigar, int readLen) {
		int length = 0;
		for(CigarElement cigEle : cigar.getCigarElements()) {
			CigarOperator cigOp = cigEle.getOperator();
			switch(cigOp) {
			case M: case I: case EQ: case X: case S:
				length += cigEle.getLength();
				break;
			default: // H, P, D or N
				break;
			}
		}
		return length == readLen;
	}

	/** calculate alignment length using Cigar
	 * @return actual alignment length including M,=,X,I,D,S but not H,P and N
	 */
	private static int calcAlnLenByCigar(Cigar cigar) {
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
	private static int calcInsertLenByCigar(Cigar cigar) {
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
	private static int calcInsertFromByCigar(int cigarLen, Cigar cigar, boolean isMinus) {
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
	private static char[] getAlnStatusBySAMRecord(SAMRecord record, int alnLen) {
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
	private static int calcInsertLenBy1DP(char[] status, int insertFrom) {
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
	private static int calcReferenceAlnLen(char[] status, int from, int alnLen) {
		int refAlnLen = 0;
		for(int i = from; i < status.length; i++)
			if(status[i] != 'S' && status[i] != 'I') // not S or I
				refAlnLen++;
		return refAlnLen;
	}

	/** Calc readInsertLength by subtracting deletions from the insertLen
	 */
	private static int calcReferenceInsertLen(char[] status, int insertFrom, int insertLen) {
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
	private static void fixSAMRecordCigarTag(SAMRecord record, int alnLen, char[] status, int from, int insertLen) {
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
/*			System.err.println(isMinus + " refClipLen:" + refClipLen);
			System.err.println(" refAlnLen:" + refAlnLen + " refInsertLen:" + refInsertLen + " insertLen:" + insertLen);
			System.err.println("newCig:" + newCig);
			System.err.println("newMisStr:" + newMisStr);
			System.err.println(record.getSAMString());
			System.exit(-1);*/
			 
			assert isMatchedCigarMisStr(newCig, newMisStr.toString());
			record.setCigar(newCig); // update the cigar
			record.setAttribute("MD", newMisStr.toString()); // update misStr
		} // end if
	}

	/** validate whether cigar and MD:Z mismatch string matches
	 * @return true if their length matches
	 */
	private static boolean isMatchedCigarMisStr(Cigar cigar, String misStr) {
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

	/**
	 * @return the sEED_LEN
	 */
	public static int getSEED_LEN() {
		return SEED_LEN;
	}

	/**
	 * @param seedLen the seedLen to set
	 */
	public static void setSEED_LEN(int seedLen) {
		if(seedLen <= 0)
			throw new IllegalArgumentException("--seed-len must be positive");
		SEED_LEN = seedLen;
	}

	/**
	 * @return the MATCH_SCORE
	 */
	public static int getMATCH_SCORE() {
		return MATCH_SCORE;
	}

	/**
	 * @param matchScore the matchScore to set
	 */
	public static void setMATCH_SCORE(int matchScore) {
		if(matchScore <= 0)
			throw new IllegalArgumentException("--match-score must be positive");
		MATCH_SCORE = matchScore;
	}

	/**
	 * @return the MIS_SCORE
	 */
	public static int getMIS_SCORE() {
		return MIS_SCORE;
	}

	/**
	 * @param misScore the misScore to set
	 */
	public static void setMIS_SCORE(int misScore) {
		if(misScore > 0)
			throw new IllegalArgumentException("--mis-score must be non-positive");
		MIS_SCORE = misScore;
	}

	/**
	 * @return the GAP_OPEN_PENALTY
	 */
	public static int getGAP_OPEN_PENALTY() {
		return GAP_OPEN_PENALTY;
	}

	/**
	 * @param gapOpenPenalty the gapOpenPenalty to set
	 */
	public static void setGAP_OPEN_PENALTY(int gapOpenPenalty) {
		if(gapOpenPenalty < 0)
			throw new IllegalArgumentException("--gap-open-penalty must be non negative");
		GAP_OPEN_PENALTY = gapOpenPenalty;
	}

	/**
	 * @return the GAP_EXT_PENALTY
	 */
	public static int getGAP_EXT_PENALTY() {
		return GAP_EXT_PENALTY;
	}

	/**
	 * @param gapExtPenalty the gapExtPenalty to set
	 */
	public static void setGAP_EXT_PENALTY(int gapExtPenalty) {
		if(gapExtPenalty <= 0)
			throw new IllegalArgumentException("--gap-ext-penalty must be positive");
		GAP_EXT_PENALTY = gapExtPenalty;
	}

	// default 1DP parameters
	private static int SEED_LEN = 25;
	private static int MATCH_SCORE = 1;
	private static int MIS_SCORE = -2;
	private static int GAP_OPEN_PENALTY = 4;
	private static int GAP_EXT_PENALTY = 1;
	private static final int REF_QUAL = 40; // reference quality for deletions
	// mismatch string patterns
	private static final Pattern misPat1 = Pattern.compile("(\\d+)(.*)");
	private static final Pattern misPat2 = Pattern.compile("([A-Z]|\\^[A-Z]+)(\\d+)");
	private static final Pattern misPat3 = Pattern.compile("\\d+|[A-Z]|\\^[A-Z]+");
}
