/**
 * 
 */
package edu.upenn.egricelab.ucsc;

import edu.upenn.egricelab.AlignerBoost.utils.StringUtils;
import edu.upenn.egricelab.ucsc.ColorRGB;

/**
 * @author zhengqi
 *
 */
public class BED12 extends BED6 {
	/* constructors */
	/** default constructor, do nothing */
	public BED12() { }
	
	public BED12(String line)
			throws IllegalArgumentException, ArrayIndexOutOfBoundsException
	{
		String[] fields = line.split(sep);
		init(fields[chromIndex], Integer.parseInt(fields[startIndex]), Integer.parseInt(fields[endIndex]),
				fields[nameIndex], Integer.parseInt(fields[scoreIndex]), fields[strandIndex],
				Integer.parseInt(fields[thickStartIndex]), Integer.parseInt(fields[thickEndIndex]), fields[rgbIndex],
				Integer.parseInt(fields[blockCountIndex]),
				fields[blockSizesIndex].split(valueSep), fields[blockStartsIndex].split(valueSep));
	}
	
	/* helper methods */
	protected void init(String chrom, int start, int end, String name, int score, String strand,
			int thickStart, int thickEnd, String rgb,
			int blockCount, String[] blockSizeStrs, String[] blockStartStrs)
			throws IllegalArgumentException, ArrayIndexOutOfBoundsException
	{
		super.init(chrom, start, end, name, score, strand);
		if(!(0 <= thickStart && thickStart <= thickEnd))
			throw new IllegalArgumentException("thickStart and thickEnd must be non-nagative and left-closed, right-open");
		itemRGB = new ColorRGB(rgb);
		if(blockCount < 0)
			throw new IllegalArgumentException("blockCount must be non-negative");
		this.blockCount = blockCount;
		blockSizes = new int[blockCount];
		blockStarts = new int[blockCount];
		blockEnds = new int[blockCount];
		for(int i = 0; i < blockCount; i++) {
			blockSizes[i] = Integer.parseInt(blockSizeStrs[i]);
			blockStarts[i] = Integer.parseInt(blockStartStrs[i]);
			blockEnds[i] = blockStarts[i] + blockSizes[i];
		}
	}
	
	/**
	 * get the string representation of this BED6 record
	 * @override the BED version
	 */
	@Override
	public String toString() {
		return super.toString() + sep + thickStart + sep + thickEnd + sep + itemRGB + sep +
				blockCount + sep +
				StringUtils.join(valueSep, blockSizes) + sep + StringUtils.join(valueSep, blockStarts);
	}
	
	/* member fields */
	private int thickStart;
	private int thickEnd;
	private ColorRGB itemRGB;
	private int blockCount;
	private int[] blockSizes;
	private int[] blockStarts;
	private int[] blockEnds;
	
	/* class constants */
	public static final String valueSep = ",";
	public static final int thickStartIndex = 6;
	public static final int thickEndIndex = 7;
	public static final int rgbIndex = 8;
	public static final int blockCountIndex = 9;
	public static final int blockSizesIndex = 10;
	public static final int blockStartsIndex = 11;
}
