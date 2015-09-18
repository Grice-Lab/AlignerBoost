/**
 * A Utility class for Methods for manipulating IUPAC Degenerative nucleotides
 */
package edu.upenn.egricelab.AlignerBoost.utils;

/**
 * @author zhengqi
 */
public class IUPACNucl {
	public static String revcom(String seq)
	throws IllegalArgumentException {
		StringBuilder rcSeq = new StringBuilder(seq.length());
		for(int i = seq.length() - 1; i >= 0; i--) {
			char b = '\0';
			switch(seq.charAt(i)) {
			case 'A':
				b = 'T';
				break;
			case 'T':
				b = 'A';
				break;
			case 'U':
				b = 'A';
				break;
			case 'G':
				b = 'C';
				break;
			case 'C':
				b = 'G';
				break;
			case 'Y':
				b = 'R';
				break;
			case 'R':
				b = 'Y';
				break;
			case 'S':
				b = 'S';
				break;
			case 'W':
				b = 'W';
				break;
			case 'K':
				b = 'M';
				break;
			case 'M':
				b = 'K';
				break;
			case 'B':
				b = 'V';
				break;
			case 'D':
				b = 'H';
				break;
			case 'H':
				b = 'D';
				break;
			case 'V':
				b = 'B';
				break;
			case 'N':
				b = 'N';
				break;
			case 'a':
				b = 't';
				break;
			case 't':
				b = 'a';
				break;
			case 'u':
				b = 'a';
				break;
			case 'g':
				b = 'c';
				break;
			case 'c':
				b = 'g';
				break;
			case 'y':
				b = 'r';
				break;
			case 'r':
				b = 'y';
				break;
			case 's':
				b = 's';
				break;
			case 'w':
				b = 'w';
				break;
			case 'k':
				b = 'm';
				break;
			case 'm':
				b = 'k';
				break;
			case 'b':
				b = 'v';
				break;
			case 'd':
				b = 'h';
				break;
			case 'h':
				b = 'd';
				break;
			case 'v':
				b = 'b';
				break;
			case 'n':
				b = 'n';
				break;
			default:
				throw new IllegalArgumentException("Invalid IUPAC Nucliotide character found in " + seq);
			} /* end of switch */
			rcSeq.append(b);
		} /* end of for */
		return rcSeq.toString();
	}
}
