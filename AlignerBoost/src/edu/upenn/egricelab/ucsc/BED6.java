/**
 * 
 */
package edu.upenn.egricelab.ucsc;

/**
 * @author zhengqi
 *
 */
public class BED6 extends BED4 {
	/* constructors */
	/** default constructor, do nothing */
	public BED6() { }
	
	public BED6(String line)
			throws IllegalArgumentException, ArrayIndexOutOfBoundsException
	{
		String[] fields = line.split(sep);
		init(fields[chromIndex], Integer.parseInt(fields[startIndex]), Integer.parseInt(fields[endIndex]),
				fields[nameIndex], Integer.parseInt(fields[scoreIndex]), fields[strandIndex]);
	}
	
	
	/* helper methods */
	protected void init(String chrom, int start, int end, String name, int score, String strand)
			throws IllegalArgumentException
	{
		super.init(chrom, start, end, name);
		if(score < 0)
			throw new IllegalArgumentException("score must be non-negative");
		if(!(strand.equals(".") || strand.equals("+") || strand.equals("-")))
			throw new IllegalArgumentException("strand must be one of '.', '+' or '-'");
		this.score = score;
		this.strand = strand;
	}
	
	/**
	 * get the string representation of this BED6 record
	 * @override the BED version
	 */
	@Override
	public String toString() {
		return super.toString() + sep + score + sep + strand;
	}
	
	/* member fields */
	private int score;
	private String strand = ".";
	
	/* class constants */
	public static final int scoreIndex = 4;
	public static final int strandIndex = 5;
}
