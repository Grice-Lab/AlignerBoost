/**
 * A class representing genomic intervals that mimic the QueryIntervals coming from HTSJDK
 * All coordiantes are 1-based inclusing
 */
package edu.upenn.egricelab.AlignerBoost.utils;

import java.util.*;

/**
 * @author zhengqi
 *
 */
public class GenomeInterval
implements Comparable<GenomeInterval>, Cloneable {
	/* constructors */
	/**
	 * Construct a GenomeInterval using given region
	 * @param chr
	 * @param start
	 * @param end
	 * @throws IllegalArgumentException
	 */
	public GenomeInterval(String chr, int start, int end)
	throws IllegalArgumentException {
		if(start > end)
			throw new IllegalArgumentException("start must not greater than end");
		this.chr = chr;
		this.start = start;
		this.end = end;
	}
	
	/**
	 * compare this GenomeInterval to another interva,
	 * first based on chr, than start, then end
	 */
	public int compareTo(GenomeInterval other) {
		return !chr.equals(other.chr) ? chr.compareTo(other.chr): 
			start != other.start ? start - other.start : end - other.end;
	}
	
	/**
	 * Test whether this GenomeInterval is equals to another interval
	 * @return  true if all fields are equal
	 */
	@Override
	public boolean equals(Object other) {
		if(other == null)
			return false;
		if(!(other instanceof GenomeInterval))
			return false;
		GenomeInterval oInt = (GenomeInterval) other;
		return chr.equals(oInt.chr) && start == oInt.start && end == oInt.end; 
	}
	
	/**
	 * Return the hashCode of this interval.
	 * Two equal intervals are guarenteed to return the same hashCode
	 */
	@Override
	public int hashCode() {
		return Arrays.hashCode(new Object[] {
				chr,
				start, // auto-boxed
				end    // auto-boxed
		});
	}
	
	/**
	 * Clone this interval
	 * @return  a new interval copy with all fields set the same
	 */
	@Override
	public GenomeInterval clone() {
		return new GenomeInterval(chr, start, end);
	}
	
	/**
	 * Check whether this interval is overlap with another interval
	 * @param other  interval to check overlapping
	 * @return  true if both intervals share the chromosome and with overlapping coordinates
	 */
	public boolean isOverlap(GenomeInterval other) {
		return chr.equals(other.chr) && start <= other.end && end >= other.start; 
	}
	
	/**
	 * Merge this interval to another interval,
	 * modify the content of this interval if merged successfully
	 * @param other  interval to be merged
	 * @return  this interval
	 */
	public GenomeInterval mergeWith(GenomeInterval other) {
		if(!isOverlap(other))
			return this;
		/* update the coordinates */
		start = start < other.start ? start : other.start;
		end = end > other.end ? end : other.end;
		return this;
	}
	
	/* Static methods */
	/**
	 * Merge two given intervals, if possible
	 * @param intervalA  first interval
	 * @param intervalB  second interval
	 * @return  return the merged interval, if overlaps, or a copy of the first interval
	 */
	public static GenomeInterval merge(GenomeInterval intervalA, GenomeInterval intervalB) {
		GenomeInterval merged = intervalA.clone(); // make a copy of the first interval
		return merged.mergeWith(intervalB);
	}
	
	/**
	 * Optimize a list of intervals
	 * @param intervals  a list of intervals
	 * @return  a new list of intervals with each element not overlapping to any other element
	 */
	public static List<GenomeInterval> optimizeIntervals(List<GenomeInterval> intervals) {
		Stack<GenomeInterval> s = new Stack<GenomeInterval>();
		if(intervals.isEmpty())
			return s;
		// sort the intervals
		Collections.sort(intervals);
		s.push(intervals.get(0));
		for(int i = 1; i < intervals.size(); i++) {
			GenomeInterval top = s.peek();
			if(!top.isOverlap(intervals.get(i))) // not overlap
					s.push(intervals.get(i));
			else // overlap
				top.mergeWith(intervals.get(i));
		}
		return s;
	}
	
	String chr;
	int start;
	int end;
}
