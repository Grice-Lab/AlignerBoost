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
 * a class to provide static method to filter and fix SAM/BAM alignment
 */
package edu.upenn.egricelab.AlignerBoost;
import java.util.Arrays;
import java.util.regex.*;

import edu.upenn.egricelab.AlignerBoost.utils.Stats;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import static edu.upenn.egricelab.AlignerBoost.EnvConstants.*;
import static edu.upenn.egricelab.AlignerBoost.utils.Stats.PHRED_SCALE;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;

/** Parse SAM/BAM alignment file with Picard java packages in CLASSPATH
 * A 1-dimensional DP algorithm (1DP) is implemented to re-estimate the insert length of the alignment
 * The calculated metrix are stored in application-specified SAM tags,
 * X?: alignment tags, Y? seed tags, Z? total tags
 * Tag  Type  Description
 * XA   i     alignment length, including M,=,X,I,D,S but not H,P,N
 * XL   i     insert length, including M,=,X,I,D but not S,H,P,N, determined by Cigar or 1DP
 * XF   i     actual insert from (start) relative to reference
 * XI   f     alignment identity as 1 - (YX + YG) / XL
 * XH   Z     alignment likelihood given this mapping loc and quality, in string format to preserve double precision
 * XV   i     known SNVs found in this alignment
 * YL	i     seed length for calculating seed mismatches and indels
 * YX   i     # of seed mismatches
 * YG   i     # of seed indels
 * ZX   i     # of all mismatches
 * ZG   i     # of all indels
 * @author Qi Zheng
 * @version 1.2
 */
public class SAMAlignFixer {
	/* Nested types and enums */
	public static enum ClipHandlingMode {
		USE, // count clipped bases as mismatches
		IGNORE, // ignore clipped bases
		END5, // only count 5' clipped bases
		END3; // only count 3' clipped bases
	}
	
	/**
	 * Fix a SAMRecord read sequence and quality, according to its previous record (with read and quality provided)
	 * @param record  the record to be fixed
	 * @param prevRecord  previos record in the SAM file that contains read and quality information
	 * @return  true if the record can be fixed
	 */
	public static boolean fixSAMRecordRead(SAMRecord record, SAMRecord prevRecord) {
		if(record.getReadLength() != 0) // no fix neccessary
			return false;
		if(!(prevRecord != null && prevRecord.getReadLength() != 0 && prevRecord.getBaseQualities() != null &&
				record.getReadName().equals(prevRecord.getReadName()))) // cannot be fixed
			return false;
		String read = prevRecord.getReadString(); // a copy of readString, can be modified
		String qual = prevRecord.getBaseQualityString(); // a copy of readString, can be modified
		int clippedStart = 0;
		int clippedEnd = read.length();
		// check cigar to see whether need hard-clip
		Cigar cigar = record.getCigar();
		int cigNum = record.getCigarLength();
		for(int i = 0; i < cigNum; i++) {
			CigarElement cigEle = cigar.getCigarElement(i);
			if(cigEle.getOperator() == CigarOperator.H) { // a hard-clip, read needs to be clipped
				if(i == 0) // 5' clip
					clippedStart += cigEle.getLength();
				else if(i == cigNum - 1)
					clippedEnd -= cigEle.getLength();
				else
					throw new RuntimeException("Invalid cigar operator at SAMRecord:" + newLine + record.getSAMString());
			}
		}
		if(clippedEnd - clippedStart == read.length()) { // clip not required
			record.setReadString(read);
			record.setBaseQualityString(qual);
		}
		else { // clip required
			record.setReadString(read.substring(clippedStart, clippedEnd));
			record.setBaseQualityString(qual.substring(clippedStart, clippedEnd));
		}
		return true;
	}

	/**
	 * Fix a SAMRecord by adding the AlignerBoost-specific tags above, and optionally fix the alignment by 1DP 
	 * @param record  SAMRecord alignment to be fixed
	 * @param knownVCF  knownSnp to look at, ignored if null
	 * @param do1DP  whether do additional 1DP fixing?
	 * @return  false if this SAMRecord does not need to be fixed because it is not-mapped or empty.
	 */
	public static boolean fixSAMRecord(SAMRecord record, VCFFileReader knownVCF, boolean do1DP) {
		int readLen = record.getReadLength();
		if(record.getReadUnmappedFlag() || record.getReferenceIndex() == -1 || readLen == 0) // non mapped read or 0-length read
			return false;
		Cigar cigar = record.getCigar();

		// Calculate alignment length
		int alnLen = calcAlnLenByCigar(cigar);
		if(alnLen == 0)
			return false;
		// valid cigar first
		assert isValidCigarLength(cigar, readLen);
		
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
				// fix Cigar and misStr
				fixSAMRecordCigarMisStr(record, alnLen, status, newInsReg, oldInsReg);
				// fix alignStart
				record.setAlignmentStart(record.getAlignmentStart() + insertFrom - oldInsReg.from);
			}
		}
		int insertLen = insertTo - insertFrom;


		// set log-likelihood tag
		double log10lik = calcAlignLik(status, record.getBaseQualities(), cigar); // likelihood w/o knownSnp
		char[] bestStatus = status; // record the best likelihood status
		if(knownVCF != null) { // try every alternative alignment w/ known varation
			String qSeq = record.getReadString();
			
			CloseableIterator<VariantContext> vars = knownVCF.query(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd());
			while(vars.hasNext()) {
				VariantContext var = vars.next();
				char[] updatedStatus = status.clone(); // make a new copy of the old status
				int[] penaltyScore = new int[status.length];
				Arrays.fill(penaltyScore, -1); // indicating default value
				if(updateKnownSnv(updatedStatus, penaltyScore, record.getReferenceName(), record.getAlignmentStart(), alnLen, qSeq, var) > 0) { // is a ready known SNP
					double updatedLog10lik = calcAlignLik(updatedStatus, record.getBaseQualities(), cigar, penaltyScore);
					if(updatedLog10lik > log10lik) { // a better likelihood
						record.setAttribute("XV", var.getChr() + ":" + var.getStart() + "-" + var.getEnd() + ":" + var.getID());
						log10lik = updatedLog10lik; // update
						bestStatus = updatedStatus; // update
					}
				}
			} // end each var
			vars.close();
		}
		record.setAttribute("XH", Double.toString(log10lik)); // use the best likelihood

		// calculate mismatches and indels
		boolean isMinus = record.getReadNegativeStrandFlag();
		int nSeedMis = 0;
		int nSeedIndel = 0;
		int nAllMis = 0;
		int nAllIndel = 0;
		int nCigars = cigar.numCigarElements();
		for(int i = insertFrom; i < insertTo; i++) {
			if(bestStatus[i] == 'X') { // mismatch
				if(!isMinus && i < SEED_LEN || isMinus && i >= alnLen - SEED_LEN)
					nSeedMis++;
				nAllMis++;
			}
			else if(bestStatus[i] == 'I' || bestStatus[i] == 'D') { // indel 
				if(!isMinus && i < SEED_LEN || isMinus && i >= alnLen - SEED_LEN)
					nSeedIndel++;
				nAllIndel++;
			}
			else
				continue;
		}
		// check soft-clips in a second pass
		if(clipMode != ClipHandlingMode.IGNORE) {
			for(int i = 0; i < alnLen && bestStatus[i] == 'S'; i++) {
				boolean clipFlag = clipMode == ClipHandlingMode.USE ||
						!isMinus && clipMode == ClipHandlingMode.END5 ||
						isMinus && clipMode == ClipHandlingMode.END3;
				if(clipFlag) {
					if(!isMinus && i < SEED_LEN || isMinus && i >= alnLen - SEED_LEN)
						nSeedMis++;
					nAllMis++;
				}
			}
			for(int i = alnLen - 1; i >= 0 && bestStatus[i] == 'S'; i--) {
				boolean clipFlag = clipMode == ClipHandlingMode.USE ||
						!isMinus && clipMode == ClipHandlingMode.END3 ||
						isMinus && clipMode == ClipHandlingMode.END5;
				if(clipFlag) {
					if(!isMinus && i < SEED_LEN || isMinus && i >= alnLen - SEED_LEN)
						nSeedMis++;
					nAllMis++;
				}
			}
			// add hard-clipped bases if neccessary
			CigarElement firstCigar = cigar.getCigarElement(0);
			CigarElement lastCigar = cigar.getCigarElement(nCigars - 1);
			if(firstCigar.getOperator() == CigarOperator.H && (clipMode == ClipHandlingMode.USE ||
					!isMinus && clipMode == ClipHandlingMode.END5 ||
					isMinus && clipMode == ClipHandlingMode.END3)) {
				nAllMis += firstCigar.getLength();
				nSeedMis += Math.min(firstCigar.getLength(), SEED_LEN);
			}
			if(lastCigar.getOperator() == CigarOperator.H && (clipMode == ClipHandlingMode.USE ||
					!isMinus && clipMode == ClipHandlingMode.END3 ||
					isMinus && clipMode == ClipHandlingMode.END5)) {
				nAllMis += lastCigar.getLength();
				nSeedMis += Math.min(lastCigar.getLength(), SEED_LEN);
			}
		}
		
		// add additional customized tags
		float identity = 1 - (nAllMis + nAllIndel) / ((float) insertLen);
		// alignment tags
		record.setAttribute("XA", alnLen);
		record.setAttribute("XL", insertLen);
		record.setAttribute("XF", insertFrom);
		record.setAttribute("XI", identity);
		//record.setAttribute("XQ", calcAlignScore(status, record.getBaseQualities()));

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
	 * @param knownVCF  knownSnp to look at, ignored if null
	 * @return false if the record does not need to be fixed due to unmapped or empty
	 */
	public static boolean fixSAMRecord(SAMRecord record, VCFFileReader knownVCF) {
		return fixSAMRecord(record, knownVCF, false);
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
		if(cigar == null)
			return 0;
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

	/** calculate alignment hard-clipped length by Cigar
	 * @param cigar  cigar to use to calculate hard-clipped length
	 * @return hard-clipped length (H)
	 */
/*	private static int calcSAMRecordHardClippedLenByCigar(Cigar cigar) {
		int hClipLen = 0;
		for(CigarElement cigEle : cigar.getCigarElements()) {
			switch(cigEle.getOperator()) {
			case H:
				hClipLen += cigEle.getLength();
				break;
			default: // H,P or N
				break; // do nothing
			}
		}
		return hClipLen;
	}*/

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

	/** Fix the MD:Z string, if it is not in the correct format
	 * @param record  a SAMRecord with MD:Z string to be fixed
	 * @return true  if the format of the MD:Z string is wrong but fixable, or false if not
	 */
	static boolean fixMisStr(SAMRecord record) {
		String misStr = record.getStringAttribute("MD");
		if(misStr == null)
			return false;
		Matcher match4 = misPat4.matcher(misStr); // parse the entire misStr
		StringBuilder newMisStr = new StringBuilder();
		boolean fixedFlag = false;
		while(match4.find()) {
			String word = match4.group();
			if(word.length() > 1 && !Character.isDigit(word.charAt(0)) && word.charAt(0) != '^') { // a consecutive mismatch string, not a valid format
				for(int i = 0; i < word.length(); i++) {
					newMisStr.append(word.charAt(i));
					if(i != word.length() - 1)
						newMisStr.append(0);
				}
				fixedFlag = true;
			}
			else
				newMisStr.append(word);
		}
		/* fixed the misStr start and end */
		if(Character.isLetter(newMisStr.charAt(0))) {
			newMisStr.insert(0, '0');
			fixedFlag = true;
		}
		if(Character.isLetter(newMisStr.charAt(newMisStr.length() - 1))) {
			newMisStr.append(0);
			fixedFlag = true;
		}
		if(fixedFlag) {
			//System.err.println("oldMisStr:" + misStr + " newMisStr:" + newMisStr);
			record.setAttribute("MD", newMisStr.toString());
		}
		return fixedFlag;
	}
	
	/** Fix SAMRecord Cigar and mismatch tag (MD:Z), given the insert from and insert length
	 *
	 */
	private static void fixSAMRecordCigarMisStr(SAMRecord record, int alnLen, char[] status,
			InsertRegion newInsReg, InsertRegion oldInsReg) {
		if(newInsReg.from == oldInsReg.from && newInsReg.to == oldInsReg.to) // no fix needed if the insert region doesn't change
			return;
		int from = newInsReg.from; // 0-based
		int to = newInsReg.to; // 1-based
		//int readLen = record.getReadLength();
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
				if(Stats.isDecimal(s)) { // a span met
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

	/** calculate Alignment score
	 * @param oldStatus  status index
	 * @param chr 
	 * @param qual  quality scores in Phred scale
	 * @return  MAPQ style score in phred scale
	 */
/*	private static int calcAlignScore(char[] status, byte[] qual) {
		if(status == null)
			return 0;
		assert status.length >= qual.length;
		int alnScore = 0; // current map prob
		int pos = 0; // relative pos on read
		for(int i = 0; i < status.length; i++) {
			float score = 0;
			switch(status[i]) {
			case 'M': case '=': // treat as match
				score = MATCH_SCORE * qual[pos++]; // treat match prob == 1
				break;
			case 'X': // mismatch
				score = MIS_SCORE * qual[pos++];
				break;
			case 'S': case 'H': case 'P': case 'N':
				break;
			case 'I': // insert consumes read
				score = i == 0 || status[i-1] != 'I' && status[i-1] != 'D' ? -GAP_OPEN_PENALTY * REF_QUAL : -GAP_EXT_PENALTY * REF_QUAL;
				pos++;
				break;
			case 'D':
				score = i == 0 || status[i-1] != 'D' && status[i-1] != 'D' ? -GAP_OPEN_PENALTY * REF_QUAL : -GAP_EXT_PENALTY * REF_QUAL;
				break;
			default:
				break; // do nothing
			}
			alnScore += score;
		}
		return alnScore;
	}*/

	private static int updateKnownSnv(char[] status, int penaltyScore[], String chr, int alnStart,
			int alnLen, String qSeq, VariantContext var) {
		assert status.length == penaltyScore.length;
		assert status.length == alnLen;

		char[] oldStatus = status.clone(); // make a copy of old status
		if(var.isFiltered() || !(var.getChr().equals(chr) && var.getStart() < alnStart + alnLen && var.getEnd() >= alnStart))
			return 0; // no status updated
		int nUpdated = 0;
		int varStart = 0;
		int varEnd = 0;
		int varLen = 0;
		String varAllele = null;
		for(int i = 0, loc = alnStart, pos = 0; i < alnLen; i++) {
			// loc is on reference
			// pos is on read
			switch(oldStatus[i]) {
			case 'N':
				loc++;
				break;
			case 'H': case 'P':
				break;
			case 'M': case '=': case 'S':
				loc++;
				pos++;
				break;
			case 'X': // an SNV/SNP
				if(i == 0 || oldStatus[i-1] != 'X')
					varStart = loc;
				if(i == oldStatus.length - 1 || oldStatus[i+1] != 'X') {
					varEnd = loc;
					varLen = varEnd - varStart + 1;
					varAllele = qSeq.substring(pos - varLen + 1, pos + 1);
					// check this variant
					if(var.getStart() == varStart && var.getEnd() == varEnd &&
							(varLen == 1 && var.isSNP() || varLen > 1 && var.isMNP()) // matched SNV or MNP found
							&& var.hasAlternateAllele(Allele.create(varAllele))) {
						nUpdated += varLen;
						Arrays.fill(status, i - varLen + 1, i + 1, (varLen == 1 ? 'V' : 'B')); // a known SNP/MNP
						Arrays.fill(penaltyScore, i - varLen + 1, i + 1, getPenaltyScoreFromInfo(var, varAllele, AFTag));
					}
				}
				loc++;
				pos++;
				break;
			case 'I': // an insertion
				if(i == 0 || oldStatus[i-1] != 'I') // an insert start
					varStart = loc;
				if(i == oldStatus.length - 1 || oldStatus[i+1] != 'I') { // an insert end
					varEnd = loc;
					varLen = varEnd - varStart + 1;
					varAllele = qSeq.substring(pos - varLen + 1, pos + 1);
					// check each variant here
					if(var.getStart() == varStart && var.getEnd() == varEnd &&
							var.isIndel() && var.hasAlternateAllele(Allele.create(varAllele))) { // known insertion
						nUpdated += varLen;
						Arrays.fill(status, i - varLen + 1, i + 1, '-');
						Arrays.fill(penaltyScore, i - varLen + 1, i + 1, getPenaltyScoreFromInfo(var, varAllele, AFTag));
					}
				}
				pos++; // insertion takes pos on read
				break;
			case 'D': // a deletion
				if(i == 0 || oldStatus[i-1] != 'D') // a deletion start
					varStart = loc;
				if(i == oldStatus.length - 1 || oldStatus[i+1] != 'D') { // a deletion end
					varEnd = loc;
					varLen = varEnd - varStart + 1;
					varAllele = qSeq.substring(pos - varLen + 1, pos + 1);
					// check each variant here
					if(var.getStart() == varStart && var.getEnd() == varEnd &&
							var.isIndel() && var.getReference().length() == varLen) { // a known deletion
						nUpdated += varLen;
						Arrays.fill(status, i - varLen + 1, i + 1, '-');
						Arrays.fill(penaltyScore, i - varLen + 1, i + 1, getPenaltyScoreFromInfo(var, varAllele, AFTag));
					}
				}
				loc++; // deletion takes loc on reference
				break;
			default:
				break;
			} // end switch
		} // end for
		return nUpdated;
	}

	/** calculate Alignment log-likelihood given the alignment status and quality
	 * @param status  status index
	 * @param qual  quality scores in Phred scale
	 * @param hClipLen  hard-clipped length
	 * @param penaltyScore  customized penaltyScores for known SNPs
	 * @return  log-likelihood of this alignment
	 */
	private static double calcAlignLik(char[] status, byte[] qual, Cigar cigar, int[] penaltyScore) {
		if(status == null)
			return Double.NaN;
		assert qual == null || status.length >= qual.length && (penaltyScore == null || status.length == penaltyScore.length);
		// make local copy of qual, or dummy array if no qual available
		byte[] baseQ;
		if(qual != null) {
			baseQ = new byte[qual.length];
			for(int i = 0; i < qual.length; i++)
				baseQ[i] = qual[i] >= MIN_PHRED_QUAL ? qual[i] : MIN_PHRED_QUAL;
		}
		else {
			// create a dummy qual array with REF_QUAL
			baseQ = new byte[status.length];
			Arrays.fill(baseQ, REF_QUAL);
		}
		
		double log10Lik = 0;
		int pos = 0; // relative pos on read
		for(int i = 0; i < status.length; i++) {
			switch(status[i]) {
			case 'M': case '=': // treat as match
				if(pos >= baseQ.length){
					System.err.println("pos:" + pos + " baseQ:" + baseQ.length + " cigar:" + cigar.toString() + " status:" + String.valueOf(status));
				}
				log10Lik += Stats.phredP2Q(1 - Stats.phredQ2P(baseQ[pos++]), -1); // use non-error prob
				break;
			case 'X': // mismatch
				log10Lik += baseQ[pos++] / -PHRED_SCALE; // use error prob directly
				break;
			case 'S': // soft-clipped
				log10Lik += baseQ[pos++] / -PHRED_SCALE - CLIP_PENALTY;
				break;
			case 'H': case 'P': case 'N': // not possible
				break;
			case 'I': // insert consumes read
				if(i == 0 || status[i-1] != 'I') // gap-open
					log10Lik -= GAP_OPEN_PENALTY;
				log10Lik -= GAP_EXT_PENALTY;
				pos++;
				break;
			case 'D':
				if(i == 0 || status[i-1] != 'D')
					log10Lik -= GAP_OPEN_PENALTY;
				log10Lik -= GAP_EXT_PENALTY;
				break;
			case 'V': // known SNP/SNV position, treat similar as match
				int snpPenalty = penaltyScore != null && penaltyScore[i] >= 0 ? penaltyScore[i] : KNOWN_SNP_PENALTY;  
				log10Lik += Stats.phredP2Q(1 - Stats.phredQ2P(baseQ[pos++]), -1) - snpPenalty;
				break;
			case '-': // known in-del position
				int indelPenalty = penaltyScore != null && penaltyScore[i] >= 0 ? penaltyScore[i] : KNOWN_INDEL_PENALTY;
				log10Lik += Stats.phredP2Q(1 - Stats.phredQ2P(baseQ[pos++]), -1) - indelPenalty;
				break;
			case 'B': // known multi-substitution
				int mnpPenalty = penaltyScore != null && penaltyScore[i] >= 0 ? penaltyScore[i] : KNOWN_MULTISUBSTITUTION_PENALTY;
				log10Lik += Stats.phredP2Q(1 - Stats.phredQ2P(baseQ[pos++]), -1) - mnpPenalty;
			default:
				break; // do nothing
			}
		}
		if(cigar.getCigarElement(0).getOperator() == CigarOperator.H) // hard-clips exist at the beginning
			log10Lik += cigar.getCigarElement(0).getLength() * (Stats.mean(qual, 0, HCLIP_SAMPLE_LEN) / -PHRED_SCALE - CLIP_PENALTY); // use average quality at 5'
		if(cigar.getCigarElement(cigar.numCigarElements()-1).getOperator() == CigarOperator.H) // hard-clips exist at the end
			log10Lik += cigar.getCigarElement(cigar.numCigarElements()-1).getLength() * (Stats.mean(qual, qual.length - HCLIP_SAMPLE_LEN, qual.length) / -PHRED_SCALE - CLIP_PENALTY); // use average quality at 5'
		return log10Lik;
	}

	/** calculate Alignment log-likelihood given the alignment status and quality, using default penalties
	 * @param status  status index
	 * @param qual  quality scores in Phred scale
	 * @param hClipLen  hard-clipped length
	 * @return  log-likelihood of this alignment
	 */
	private static double calcAlignLik(char[] status, byte[] qual, Cigar cigar) {
		return calcAlignLik(status, qual, cigar, null);
	}
	
/*	*//** calculate Alignment log-likelihood given only alignment status with no quality (from a FASTA alignment)
	 * @param status  status index
	 * @return  log-likelihood of this alignment
	 *//*
	private static double calcAlignLik(char[] status, int hClipLen) {
		if(status == null)
			return 0;
		// make local copy of qual
		double log10Lik = 0;
		for(int i = 0; i < status.length; i++) {
			switch(status[i]) {
			case 'M': case '=': // treat as match
				log10Lik += Stats.phredP2Q(1 - Stats.phredQ2P(REF_QUAL), -1); // use non-error prob
				break;
			case 'X': // mismatch
				log10Lik += REF_QUAL / -PHRED_SCALE; // use error prob directly
				break;
			case 'S': // soft-clipped
				log10Lik += REF_QUAL / -PHRED_SCALE * CLIP_PENALTY;
				break;
			case 'H': case 'P': case 'N': // not possible
				break;
			case 'I': case 'D': // gap
				log10Lik += i == 0 || status[i-1] != 'I' && status[i-1] != 'D' ?
						GAP_OPEN_PENALTY * REF_QUAL / -PHRED_SCALE : GAP_EXT_PENALTY * REF_QUAL / -PHRED_SCALE;
				break;
			default:
				break; // do nothing
			}
		}
		if(hClipLen > 0) // hard-clips exist
			log10Lik += hClipLen * REF_QUAL / -PHRED_SCALE * CLIP_PENALTY;
		return log10Lik;
	}
*/
	
	private static int getPenaltyScoreFromInfo(VariantContext var, String varAllele, String tag) {
		if(!var.hasAttribute(tag))
			return -1; // use -1 indicator
		String[] afVals = var.getAttributeAsString(tag, "").split(",");
		assert var.getNAlleles() == afVals.length + 1; // afVals only contains altAlleles
		// check each allele of this variations
		int alleleIdx = var.getAlleleIndex(Allele.create(varAllele));
		if(alleleIdx != -1) {
			String val = afVals[alleleIdx - 1];
			if(val.startsWith("["))
				val = val.substring(1);
			if(val.endsWith("]"))
				val = val.substring(0, val.length() - 1);
			double logLik = - Math.log10(Double.parseDouble(val));
			if(!Double.isInfinite(logLik) && logLik >= 0)
				return (int) Math.round(logLik);
		}
		return -1;
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
			throw new IllegalArgumentException("SEED_LEN must be positive");
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
			throw new IllegalArgumentException("MATCH_SCORE must be positive");
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
			throw new IllegalArgumentException("MIS_SCORE must be non-positive");
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
			throw new IllegalArgumentException("GAP_OPEN_PENALTY must be non negative");
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
			throw new IllegalArgumentException("GAP_EXT_PENALTY must be positive");
		GAP_EXT_PENALTY = gapExtPenalty;
	}

	/**
	 * @return the CLIP_PENALTY
	 */
	public static int getCLIP_PENALTY() {
		return CLIP_PENALTY;
	}

	/**
	 * @param clipPenalty the clipPenalty to set
	 */
	public static void setCLIP_PENALTY(int clipPenalty) {
		if(clipPenalty < 0)
			throw new IllegalArgumentException("CLIP_PENALTY must be non-negative");
		CLIP_PENALTY = clipPenalty;
	}

	/**
	 * @return the KNOWN_SNP_PENALTY
	 */
	public static int getKNOWN_SNP_PENALTY() {
		return KNOWN_SNP_PENALTY;
	}

	/**
	 * @param knownSnpPenalty the knownSnpPenalty to set
	 */
	public static void setKNOWN_SNP_PENALTY(int knownSnpPenalty) {
		if(knownSnpPenalty < 0)
			throw new IllegalArgumentException("KNOWN_SNP_PENALTY must be non-negative");
		KNOWN_SNP_PENALTY = knownSnpPenalty;
	}

	/**
	 * @return the KNOWN_INDEL_PENALTY
	 */
	public static int getKNOWN_INDEL_PENALTY() {
		return KNOWN_INDEL_PENALTY;
	}

	/**
	 * @param knownIndelPenalty the knownIndelPenalty to set
	 */
	public static void setKNOWN_INDEL_PENALTY(int knownIndelPenalty) {
		if(knownIndelPenalty < 0)
			throw new IllegalArgumentException("KNOWN_INDEL_PENALTY must be non-negative");
		KNOWN_INDEL_PENALTY = knownIndelPenalty;
	}

	/**
	 * @return the KNOWN_MULTISUBSTITUTION_PENALTY
	 */
	public static int getKNOWN_MULTISUBSTITUTION_PENALT() {
		return KNOWN_MULTISUBSTITUTION_PENALTY;
	}

	/**
	 * @param knownMultiSubstitutionPenalty the knownMultiSubstitutionPenalty to set
	 */
	public static void setKNOWN_MULTISUBSTITUTION_PENALT(int knownMultiSubstitutionPenalty) {
		KNOWN_MULTISUBSTITUTION_PENALTY = knownMultiSubstitutionPenalty;
	}

	/**
	 * @return the AFTag
	 */
	public static String getAFTag() {
		return AFTag;
	}

	/**
	 * @param AFTag the AFTag to set
	 */
	public static void setAFTag(String tag) {
		AFTag = tag;
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
	
	static int SEED_LEN = 25;
	static int MATCH_SCORE = 1;
	static int MIS_SCORE = -2;
	static int GAP_OPEN_PENALTY = 4;
	static int GAP_EXT_PENALTY = 1;
	static int CLIP_PENALTY = 0; // additional CLIP_PENALTY except for the mismatch penalty
	static ClipHandlingMode clipMode = ClipHandlingMode.END5;
	static int KNOWN_SNP_PENALTY = 0;
	static int KNOWN_INDEL_PENALTY = 2;
	static int KNOWN_MULTISUBSTITUTION_PENALTY = 2;
	static String AFTag = "AF"; // alleleFrequency tag

	private static final byte REF_QUAL = 40; // reference quality for deletions
//	private static final byte AVG_READ_QUAL = 25;
	private static final byte MIN_PHRED_QUAL = 1; // min phred qual to avoid -Inf
	private static final int HCLIP_SAMPLE_LEN = 5; // sampling length for estimating the average quality of hard-clipped regions
	// mismatch string patterns
	private static final Pattern misPat1 = Pattern.compile("(^\\d+)(.*)");
	private static final Pattern misPat2 = Pattern.compile("([A-Z]|\\^[A-Z]+)(\\d+)");
	private static final Pattern misPat3 = Pattern.compile("\\d+|[A-Z]|\\^[A-Z]+");
	private static final Pattern misPat4 = Pattern.compile("\\d+|[A-Z]+|\\^[A-Z]+");
}
