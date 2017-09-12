/**
 * 
 */
package edu.upenn.egricelab.ucsc;

import java.util.*;

/**
 * @author zhengqi
 * @version v1.1
 */
public class BED implements Comparable<BED> {
	/* constructors */
	/** default constructor */
	public BED() { 	}
	
	/**
	 * construct a BED record from given basic info
	 * @throws IllegalArgumentException if start or end coordinate is invalid
	 */
	public BED(String chrom, int start, int end) throws IllegalArgumentException {
		init(chrom, start, end);
	}
	
	/**
	 * constrct a BED record from a given string
	 * @throws ArrayIndexOutOfBoundsException if line has fewer than 3 records
	 * @throws IllegalArgumentException if start or end coordinate is invalid
	 */
	public BED(String line)
			throws IllegalArgumentException, ArrayIndexOutOfBoundsException
	{
		String[] fields = line.split(sep);
		init(fields[chromIndex], Integer.parseInt(fields[startIndex]), Integer.parseInt(fields[endIndex]));
	}

	/**
	 * Test whether this BED object equals an given object using their basic 3 fields
	 */
	@Override
	public boolean equals(Object o) {
		if(!(o != null && o instanceof BED))
			return false;
		BED other = (BED) o;
		return chrom.equals(other.chrom) && start == other.start && end == other.end; 
	}
	
	/**
	 * Get hashCode using the basic 3 fields of the BED object
	 * guarantee to generate same hashCode if two BED objects are equal
	 */
	@Override
	public int hashCode() {
		return Arrays.hashCode(new Object[] {chrom, start, end});
	}
	
	/**
	 * compare a BED object to another using chrom, than start, than end
	 */
	@Override
	public int compareTo(BED other) {
		if(!chrom.equals(other.chrom))
			return chrom.compareTo(chrom);
		else if(start != other.start)
			return start - other.start;
		else
			return end - other.end;
	}
	
	/**
	 * get the string representation of this BED record
	 * @override the Object version
	 */
	@Override
	public String toString() {
		return chrom + sep + start + sep + end;
	}
	
	/* helper methods */
	/**
	 * initiate a bed record using given basic information
	 * @throws IllegalArgumentException if start or end coordinate is invalid
	 */
	protected void init(String chrom, int start, int end)
	throws IllegalArgumentException
	{
		if(!(0 <= start && start <= end))
			throw new IllegalArgumentException("BED format coordinates must be non-nagative and left-closed, right-open");
		this.chrom = chrom;
		this.start = start;
		this.end = end;
	}

	/* member fields */
	private String chrom;
	private int start; /* 0-based chrom start */
	private int end;   /* 1-based chrom end */
	
	/* class constants */
	public static final String sep = "\t";
	public static final int chromIndex = 0;
	public static final int startIndex = 1;
	public static final int endIndex = 2;
}
