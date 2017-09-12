/**
 * 
 */
package edu.upenn.egricelab.AlignerBoost.utils;

import java.util.*;

/**
 * A GeneticType Index class using a BitVector as underlying storage
 * @author zhengqi
 * @version v1.1
 *
 */
public class GTypeIndex {
	/* constructors */
	/** Default constructor, initiate index with default size */
	public GTypeIndex() {
		bitShift = new ArrayList<String>(DEFAULT_NBIT);
		chrIdx = new HashMap<String, BitSet>();
	}
	
	/* member methods */
	/** get number of gtypes */
	public int getNumGTypes() {
		return bitShift.size();
	}
	
	/** get shift for this gtype */
	public int getShift(String gtype) {
		return bitShift.indexOf(gtype);
	}
	
	public int getChrLen(String chr) {
		return chrIdx.get(chr).size() / maxBits;
	}
	
	/**
	 * add a new Gtype with increasing shift, if not exists
	 * @gtype  gtype to be added
	 */
	void addGType(String gtype) {
		if(!bitShift.contains(gtype)) {
			bitShift.add(gtype);
			updateIndex();
		}
	}

	/** add a new chrom with given length */
	public void addChr(String chr, int len) {
		chrIdx.put(chr, new BitSet(len * maxBits));
	}
	
	/** remove a chrom from current index */
	public BitSet removeChr(String chr) {
		return chrIdx.remove(chr);
	}
	
	/** update chrIdx if current number of gtypes exceeds the maxBits used in the index */
	void updateIndex() {
		if(bitShift.size() <= maxBits) /* no update necessary */
			return;
		int oldMaxBits = maxBits;
		maxBits += (int) Math.ceil(maxBits * LOADING_FACTOR);
		for(Map.Entry<String, BitSet> pair : chrIdx.entrySet()) {
			String chr = pair.getKey();
			BitSet oldIdx = pair.getValue();
			int chrLen = getChrLen(chr);
			/* create a new index */
			BitSet idx = new BitSet(chrLen * maxBits);
			/* copy the old index bits */
			for(int i = 0; i < chrLen; i++)
				for(int j = 0; j < oldMaxBits; j++)
					idx.set(i *  maxBits + j, oldIdx.get(i * oldMaxBits + j));
			pair.setValue(idx); /* update index */
		}
	}
	
	/**
	 * mask a chrom region as a given type
	 * @param chr  chrom
	 * @param start  start (inclusive)
	 * @param end  end (inclusive)
	 * @param gtype  genetic type
	 */
	public void maskRegion(String chr, int start, int end, String gtype) {
		addGType(gtype);
		maskRegion(chr, start, end, getShift(gtype));
	}

	/**
	 * mask a given chrom region with given shift as checked
	 * @param chr  chrom
	 * @param start  0-based start
	 * @param end  0-based end
	 * @param shift  bit shift
	 */
	public void maskRegion(String chr, int start, int end, int shift) {
		BitSet idx = chrIdx.get(chr);
		for(int i = start; i <= end; i++)
			maskLoc(idx, i, shift);
	}
	
	/**
	 * mask the index at given genome loc
	 * @param idx  BitSet index
	 * @param loc  0-based location
	 * @param shift  bit shift
	 */
	public void maskLoc(BitSet idx, int loc, int shift) {
		idx.set(loc * maxBits + shift);
	}
	
	/**
	 * get bit from given index at given genome location
	 * @param idx
	 * @param loc
	 * @param shift
	 * @return
	 */
	public boolean getLocBit(BitSet idx, int loc, int shift) {
		return idx.get(loc * maxBits + shift);
	}
	
	/**
	 * get BitSet from a given genome region
	 * @param chr  chrom
	 * @param start  0-based start
	 * @param end  0-based end
	 * @return  a new BitSet representing this region
	 */
	public BitSet getRegionBitSet(String chr, int start, int end) {
		int numBits = getNumGTypes();
		BitSet idx = chrIdx.get(chr);
		BitSet bits = new BitSet(numBits);
		for(int i = start; i <= end; i++)
			for(int j = 0; j < numBits; j++)
				if(getLocBit(idx, i, j))
					bits.set(j);
		return bits;
	}
	
	/**
	 * unmask a given genome region
	 * @param chr  chrom
	 * @param start  0-based start
	 * @param end  0-based end
	 * @return  a set of matched gtypes
	 */
	public Set<String> unmask(String chr, int start, int end) {
		int numBits = bitShift.size();
		Set<String> hits = new HashSet<String>(numBits);
		BitSet idx = chrIdx.get(chr);
		for(int i = start; i <= end; i++)
			for(int j = 0; j < numBits; j++)
				if(getLocBit(idx, i, j))
					hits.add(bitShift.get(j));
		return hits;
	}
	
	/**
	 * unmask a given genome region and get the total matched length of each gtype
	 * @param chr  chrom
	 * @param start  0-based start
	 * @param end  0-based end
	 * @return  a map of matched lengths of each gtype
	 */
	public Map<String, Integer> unmaskSum(String chr, int start, int end) {
		int numBits = bitShift.size();
		Map<String, Integer> hitSum = new HashMap<String, Integer>(numBits);
		BitSet idx = chrIdx.get(chr);
		for(int i = start; i <= end; i++)
			for(int j = 0; j < numBits; j++)
				if(getLocBit(idx, i, j))
					hitSum.put(bitShift.get(j), hitSum.getOrDefault(bitShift.get(j), 0) + 1);
		return hitSum;
	}
	
	/* member fields */
	private Map<String, BitSet> chrIdx; /* per-chromosome indices */
	private List<String> bitShift;
	private int maxBits = DEFAULT_NBIT;
	
	/* class constants */
	public static final int DEFAULT_NBIT = 16;
	public static final double LOADING_FACTOR = 0.75;
	
	/* static utility methods */
}
