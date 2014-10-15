/** 
 * a class to provide static method to filter and fix SAM/BAM alignment
 */
package net.sf.AlignerBoost;
import java.util.List;
import java.util.regex.*;

import htsjdk.samtools.*;

/** Parse SAM/BAM alignment file with Picard java packages in CLASSPATH
 * A 1-dimensional DP algorithm (1DP) is implemented to re-estimate the insert length of the alignment
 * The calculated metrix are stored in application-specified SAM tags,
 * X?: alignment tags, Y? seed tags, Z? total tags
 * Tag  Type  Description
 * XL   i     insert length, including M,=,X,I,D but not S,H,P,N, determined by Cigar or 1DP
 * XF   i     actual insert from (start) relative to reference
 * XI   f     alignment identity 1 - (YX + YG) / XL
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
		Cigar cigar = record.getCigar();
		// valid cigar first
		assert isValidCigarLength(cigar, readLen);

		// Calculate alignment length
		int alnLen = calcAlnLenByCigar(cigar);
		// get align status index
		char[] status = getAlnStatusBySAMRecord(record, alnLen);
		// get InsertRegion from either Cigar or 1DP
		InsertRegion oldInsReg = calcInsertRegionByAlnStatus(status);
		int insertFrom = oldInsReg.from; // 0-based
		int insertTo = oldInsReg.to; // 1-based
		// alignment need to be fixed after 1DP
		if(do1DP) {
			InsertRegion newInsReg = calcInsertRegionBy1DP(status);
			if(!(oldInsReg.from == newInsReg.from && oldInsReg.to == newInsReg.to)) { // if insertRegion changed
				insertFrom = newInsReg.from;
				insertTo = newInsReg.to;
				fixSAMRecordCigarMisStr(record, alnLen, status, newInsReg, oldInsReg);
			}
		}
		int insertLen = insertTo - insertFrom;
		// calculate nmismatches and indels
		int nSeedMis = 0;
		int nSeedIndel = 0;
		int nAllMis = 0;
		int nAllIndel = 0;
		for(int i = insertFrom; i < insertTo; i++) {
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

	/** get align status ('M', '=', 'X', 'I', 'S') given cigar and mismatch tag, if exists
	 * @return a char array index with length = alnLen, and always in the reference orientation
	 */
	private static char[] getAlnStatusBySAMRecord(SAMRecord record, int alnLen) {
		Cigar cigar = record.getCigar();
		String misStr = record.getStringAttribute("MD");

		// Initiate align status index
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
/*		// reverse status index, if on minus strand, so it is always in read orientation
		if(isMinus)
			reverse(status);*/
		return status;
	}
/*

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

	/**
	 * calculate insert region with given alnStatus
	 * @param alnStatus  array of char index of insert status
	 * @return insertRegion
	 */
	private static InsertRegion calcInsertRegionByAlnStatus(char[] status) {
		if(status == null || status.length == 0)
			return null;
		int alnLen = status.length;
		int from = 0;
		int to = alnLen;
		// determine from by searching alnStatus 5'
		for(from = 0; from < alnLen && status[from] == 'S'; from++)
			continue;
		// determine to by searching alnStatus 3'
		for(to = alnLen; to > 0 && status[to - 1] == 'S'; to--)
			continue;
		return new InsertRegion(from, to);
	}
	
	/**
	 * calculate insert region with given alnStatus, using original status or by 1DP
	 * @param alnStatus  array of char index of insert status
	 * @param do1DP  whether to do1DP
	 * @return insertRegion
	 */
	private static InsertRegion calcInsertRegionBy1DP(char[] status) {
		if(status == null || status.length == 0)
			return null;
		int alnLen = status.length;
		int from = 0;
		int to = alnLen;
		// do 1DP to determine insert region
		int[] dpScore = new int[alnLen + 1]; // Position 0 is dummy
		int maxScore = Integer.MIN_VALUE;
		// dynamic-programming to get dp-score, highest score and insert region simultaneously
		for(int i = 0; i < alnLen; i++) {
			int score;
			switch(status[i]) {
			case 'M': case '=':
				score = MATCH_SCORE;
				break;
			case 'X':
				score = MIS_SCORE;
				break;
			case 'I': case 'D': // indel
				score = i == 0 || status[i - 1] != 'I' && status[i - 1] != 'D' ? -GAP_OPEN_PENALTY : -GAP_EXT_PENALTY;
				break;
			case 'S':
				score = 0;
				break;
			default: // do nothing
				score = 0;
				break;
			}
			dpScore[i + 1] = dpScore[i] + score >= 0 ? dpScore[i] + score : 0; // dpScore is 1-based
			if(dpScore[i + 1] > maxScore) {
				to = i + 1;   // update to 1-based
				maxScore = dpScore[i + 1];  // Update maxScore
			}
		} // end each status[]
		// track back to get from
		for(from = to; from >= 0 && dpScore[from] > 0; from--)
			continue;
		return new InsertRegion(from, to);
	}

	/** Calculate refFrom from from
	 * @param status  alignment status index array
	 * @param from  align from
	 * @return refFrom relative to the reference
	 */
	private static int calcReferenceFrom(char[] status, int from) {
		int refFrom = from;
		for(int i = 0; i < from; i++)
			if(status[i] == 'I') // I
				refFrom--;
		return refFrom;
	}

	/** Calculate refTo from to
	 * @param status  alignment status index array
	 * @param to  align to
	 * @return refTo relative to the reference
	 */
	private static int calcReferenceTo(char[] status, int to) {
		int refTo = to;
		for(int i = 0; i < to; i++)
			if(status[i] == 'I') // I
				refTo--;
		return refTo;
	}

	/** Fix SAMRecord Cigar and mismatch tag (MD:Z), given the insert from and insert length
	 *
	 */
	private static void fixSAMRecordCigarMisStr(SAMRecord record, int alnLen, char[] status,
			InsertRegion newInsReg, InsertRegion oldInsReg) {
		int from = newInsReg.from; // 0-based
		int to = newInsReg.to; // 1-based
		if(from == 0 && to == alnLen) // no fix needed
			return;
		
		int readLen = record.getReadLength();
		
		// calculate readFrom and readTo
/*		int readFrom = calcReadInsertFrom(status, from, alnLen);
		int readTo = calcReadInsertTo(status, to, alnLen);*/
		
		// calculate refFrom and refTo
		int refFrom = calcReferenceFrom(status, from);
		int refTo = calcReferenceTo(status, to);
		// fix Cigar string of the alignment
		Cigar oldCig = record.getCigar();
		Cigar newCig = new Cigar();
		if(from > 0) // soft-clip exists at 5'
			newCig.add(new CigarElement(from, CigarOperator.S));
		int pos = 0; // relative pos to alignment
		for(CigarElement cigEle : oldCig.getCigarElements()) {
			CigarOperator cigOp = cigEle.getOperator();
			int cigLen = cigEle.getLength();
			switch(cigOp) {
			case M: case EQ: case X: case I: case D:
				if(pos + cigLen <= from || pos >= to) // competely in soft-clipped region
					; // do nothing
				else if(pos < to && pos + cigLen > from) { // partially clipped
					int clipFrom = pos > from ? pos : from;
					int clipTo = pos + cigLen < to ? pos + cigLen : to;
					newCig.add(new CigarElement(clipTo - clipFrom, cigOp));
				}
				else // not in clipped region
					newCig.add(new CigarElement(cigLen, cigOp)); // make a copy
				pos += cigLen;
				break;
			case S:
				pos += cigLen; // count length but do nothing else
				break;
			default: // ('N', 'H', 'P')
				if(pos >= from && pos + cigLen <= to) // not clipped
					newCig.add(new CigarElement(cigLen, cigOp));
				break;
			}
		} // end each cigEle
		if(to < alnLen) // soft-clip exists at 3'
			newCig.add(new CigarElement(alnLen - to, CigarOperator.S));
		//assert isValidCigarLength(newCig, readLen);
/*		if(calcAlnLenByCigar(newCig) != calcAlnLenByCigar(oldCig)) {
			System.err.println("Cigar length doesn't match at:\n" + oldCig + " <-> " + newCig);
			System.err.printf("from:%d to:%d%n", from, to);
			System.err.println(record.getSAMString());
			System.exit(-1);
		}*/
		assert calcAlnLenByCigar(newCig) == calcAlnLenByCigar(oldCig);
		
		// fix mismatch string, if exists
		String oldMisStr = record.getStringAttribute("MD");
		if(oldMisStr != null) {
			StringBuilder newMisStr = new StringBuilder(); // use StringBuilder for performance
			Matcher match = misPat3.matcher(oldMisStr);
			int refPos = oldInsReg.from; // relative pos to the reference, excluding 'I'
			while(match.find()) {
				String s = match.group();
				if(isDecimal(s)) { // a span met
					int span = Integer.parseInt(s);
					if(refPos + span <= refFrom || refPos >= refTo) // completely in soft-clipped region
						; // do nothing
					else if(refPos < refTo && refPos + span > refFrom) { // partially clipped
						int clipFrom = refPos > refFrom ? refPos : refFrom;
						int clipTo = refPos + span < refTo ? refPos + span : refTo;
						newMisStr.append(clipTo - clipFrom);
					}
					else // not in clipped region
						newMisStr.append(span);
					refPos += span;
				}
				else { // a mismatch or deletion tag
					if(refPos >= refFrom && refPos < refTo) // not in clip region, note deletion won't partially in clipped region
						newMisStr.append(s);
					if(s.length() == 1) // a mismatch
						refPos++;
					else // a deletion
						refPos += s.length() - 1;
				}
			} // end while find
			if(!Character.isDigit(newMisStr.charAt(newMisStr.length() - 1))) // if the new misStr deosn't end with number
					newMisStr.append("0");
/*			if(!isMatchedCigarMisStr(newCig, newMisStr.toString())) {
				System.err.printf("oldCig:%s oldMisStr:%s%nnewCig:%s newMisStr:%s%n", oldCig, oldMisStr, newCig, newMisStr);
				System.err.printf("from:%d to:%d%n", from, to);
				System.err.printf("refFrom:%d refTo:%d%n", refFrom, refTo);
				System.exit(-1);
			}*/
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
/*		if(cigRefAlnLen != misRefAlnLen)
			System.err.printf("cigRefAlnLen:%d misRefAlnLen:%d%n", cigRefAlnLen, misRefAlnLen);*/
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
				score = (status[i-1] != 'I' && status[i-1] != 'D' ? -GAP_OPEN_PENALTY : -GAP_EXT_PENALTY); 
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
				score = (status[i-1] != 'I' ? -GAP_OPEN_PENALTY : -GAP_EXT_PENALTY);
				weight += qual[pos++];
			case 'D':
				score = (status[i-1] != 'D' ? -GAP_OPEN_PENALTY : -GAP_EXT_PENALTY);
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

	/**
	 * a nested static class member holding POD of a insertRegion
	 * @author Qi Zheng
	 */
	public static class InsertRegion {
		/**
		 * construct a InsertRegion
		 * @param from
		 * @param to
		 */
		public InsertRegion(int from, int to) {
			this.from = from;
			this.to = to;
		}
		
		int from;
		int to;
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
