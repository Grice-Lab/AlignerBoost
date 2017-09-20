/**
 * 
 */
package edu.upenn.egricelab.ucsc;

/**
 * @author zhengqi
 *
 */
public class BED4 extends BED {
	/* constructors */
	/** default constructor, do nothing */
	public BED4() { }
	
	public BED4(String line)
			throws NumberFormatException, ArrayIndexOutOfBoundsException
	{
		String[] fields = line.split(sep);
		init(fields[chromIndex], Integer.parseInt(fields[startIndex]), Integer.parseInt(fields[endIndex]),
				fields[nameIndex]);
	}
	
	
	/* helper methods */
	protected void init(String chrom, int start, int end, String name)
			throws IllegalArgumentException
	{
		super.init(chrom, start, end);
		this.name = name;
	}
	
	/**
	 * get the string representation of this BED6 record
	 * @override the BED version
	 */
	@Override
	public String toString() {
		return super.toString() + sep + name;
	}
	
	/* member fields */
	private String name;
	
	/* class constants */
	public static final int nameIndex = 3;
}
