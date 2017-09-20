/**
 * 
 */
package edu.upenn.egricelab.ucsc;

import java.util.*;
/**
 * A General Feature Format (GFF) format represents a genetic feature or annotation
 * @author zhengqi
 * @version v1.1
 *
 */
public abstract class GFF implements Comparable<GFF> {
	/* constructors */
	/** default constructor, do nothing */
	public GFF() {  }
	
	/** construct a GFF record from String */
	public GFF(String line) throws IllegalArgumentException, ArrayIndexOutOfBoundsException
	{
		String[] fields = line.split(sep);
		init(fields[seqnameIndex], fields[sourceIndex], fields[typeIndex],
				Integer.parseInt(fields[startIndex]), Integer.parseInt(fields[endIndex]),
				fields[scoreIndex].equals(".") ? Double.NaN : Double.parseDouble(fields[scoreIndex]), /* score */
				fields[strandIndex], /* strand */
				fields[frameIndex].equals(".") ? -1 : Integer.parseInt(fields[frameIndex]), /* frame */
				fields[attrsIndex]);
	}
	
	/**
	 * @return the seqname
	 */
	public String getSeqname() {
		return seqname;
	}

	/**
	 * @param seqname the seqname to set
	 */
	public void setSeqname(String seqname) {
		this.seqname = seqname;
	}

	/**
	 * @return the source
	 */
	public String getSource() {
		return source;
	}

	/**
	 * @param source the source to set
	 */
	public void setSource(String source) {
		this.source = source;
	}

	/**
	 * @return the type
	 */
	public String getType() {
		return type;
	}

	/**
	 * @param type the type to set
	 */
	public void setType(String type) {
		this.type = type;
	}

	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}

	/**
	 * @param start the start to set
	 */
	public void setStart(int start) {
		this.start = start;
	}

	/**
	 * @return the end
	 */
	public int getEnd() {
		return end;
	}

	/**
	 * @param end the end to set
	 */
	public void setEnd(int end) {
		this.end = end;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	public void setScore(double score) {
		this.score = score;
	}

	/**
	 * @return the strand
	 */
	public String getStrand() {
		return strand;
	}

	/**
	 * @param strand the strand to set
	 */
	public void setStrand(String strand) {
		this.strand = strand;
	}

	/**
	 * @return the frame
	 */
	public int getFrame() {
		return frame;
	}

	/**
	 * @param frame the frame to set
	 */
	public void setFrame(int frame) {
		this.frame = frame;
	}

	/**
	 * @return the attrNames
	 */
	public List<String> getAttrNames() {
		return attrNames;
	}
	
	/* member methods */
	public boolean hasAttr(String name) {
		return attrMap.containsKey(name);
	}
	
	public String getAttr(String name) {
		return attrMap.get(name);
	}
	
	public void setAttr(String name, String value) {
		if(!attrNames.contains(name))
			attrNames.add(name);
		attrMap.put(name, value);
	}

	@Override
	public boolean equals(Object o) {
		if(!(o != null && o instanceof GFF))
			return false;
		GFF other = (GFF) o;
		return seqname.equals(other.seqname) && start == other.start && end == other.end
				&& source.equals(other.source) && type.equals(other.type);
	}
	
	@Override
	public int compareTo(GFF other) {
		if(!seqname.equals(other.seqname))
			return seqname.compareTo(other.seqname);
		else if(start != other.start)
			return start - other.start;
		else if(end != other.end)
			return end - other.end;
		else if(!source.equals(other.source))
			return source.compareTo(other.source);
		else
			return type.compareTo(other.type);
	}
	
	@Override
	public int hashCode() {
		return Arrays.hashCode(new Object[] {
				seqname,
				start,
				end,
				source,
				type
		});
	}
	
	@Override
	public String toString() {
		return seqname + sep + source + sep + type + sep
				+ start + sep + end + sep
				+ (Double.isNaN(score) ? "." : Double.toString(score)) + sep
				+ strand + sep
				+ (frame >= 0 ? Integer.toString(frame) : ".") + sep
				+ writeAttributes();
	}
	
	/** read attributes into a map from a String */
	public abstract void readAttributes(String attrStr);
	
	/** write attributes in a given order */
	public abstract String writeAttributes();

	/* helper methods */
	protected void init(String seqname, String source, String type,
			int start, int end, double score,
			String strand, int frame, String attrStr) throws IllegalArgumentException
	{
		if(!(0 <= start && start <= end))
			throw new IllegalArgumentException("GFF coordinates must be non-negative");
		if(!(strand.equals(".") || strand.equals("+") || strand.equals("-")))
			throw new IllegalArgumentException("GFF strand must be one of '.', '+' or '-'");
		this.seqname = seqname;
		this.source = source;
		this.type = type;
		this.start = start;
		this.end = end;
		this.score = score;
		this.strand = strand;
		this.frame = frame;
		readAttributes(attrStr);
	}
	
	/* member fields */
	private String seqname;
	private String source;
	private String type;
	private int start;
	private int end;
	private double score;
	private String strand;
	private int frame;
	protected List<String> attrNames; /* attribute names */
	protected Map<String, String> attrMap; /* attribute name->value map */
	
	/* class constants */
	public static final String sep = "\t";
	public static final int seqnameIndex = 0;
	public static final int sourceIndex = 1;
	public static final int typeIndex = 2;
	public static final int startIndex = 3;
	public static final int endIndex = 4;
	public static final int scoreIndex = 5;
	public static final int strandIndex = 6;
	public static final int frameIndex = 7;
	public static final int attrsIndex = 8;
}
