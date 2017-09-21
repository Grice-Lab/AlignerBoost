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
	/* embedded types and enums */
	public static class BitMask {
		/* constructors */
		public BitMask(int length) {
			this.length = length;
			gtypeIdx = new HashMap<String, BitSet>(DEFAULT_NTYPES);
		}
		
		/* member methods */
		/** get length of this BitMask */
		public int getLength() {
			return length;
		}
		
		/** get number of gtypes */
		public int getNumGType() {
			return gtypeIdx.size();
		}
		
		/** add a new gtype into this BitMask, if not exists */
		public void add(String gtype) {
			if(!gtypeIdx.containsKey(gtype))
				gtypeIdx.put(gtype, new BitSet(length));
		}
		
		/** remove a gtype from this BitMask */
		public void remove(String gtype) {
			gtypeIdx.remove(gtype);
		}
		
		/** get bit of given gtype at given location */
		public boolean getBit(String gtype, int loc) {
			return gtypeIdx.get(gtype).get(loc);
		}
		
		/** set bit of given gtype at given location */
		public void setBit(String gtype, int loc) {
			gtypeIdx.get(gtype).set(loc);
		}
		
		/**
		 * set all bits of a given gtype in given genome region
		 * @param gtype  gtype to be set
		 * @param start  0-based start
		 * @param end    1-based end
		 */
		public void setBit(String gtype, int start, int end) {
			BitSet idx = gtypeIdx.get(gtype);
			idx.set(start, end);
		}
		
		/** unmask the gtypes matched at a given location */
		public Set<String> unmask(int loc) {
			Set<String> hits = new HashSet<String>(getNumGType());
			for(Map.Entry<String, BitSet> pair : gtypeIdx.entrySet()) {
				String gtype = pair.getKey();
				BitSet idx = pair.getValue();
				if(idx.get(loc))
					hits.add(gtype);
			}
			return hits;
		}
		
		/**
		 * unmask the gtypes matched at a given genomic region
		 * @param start  0-based start
		 * @param end  1-based end
		 * @return  a Set of matched gtypes
		 */
		public Set<String> unmask(int start, int end) {
			Set<String> hits = new HashSet<String>(getNumGType());
			for(Map.Entry<String, BitSet> pair : gtypeIdx.entrySet()) {
				String gtype = pair.getKey();
				BitSet idx = pair.getValue();
				for(int i = start; i < end; i++)
					if(idx.get(i)) {
						hits.add(gtype); /* no search needed */
						break;
					}
			}
			return hits;
		}
		
		/**
		 * unmask the gtypes matched summary at a given genomic region
		 * @param start  0-based start
		 * @param end  1-based end
		 * @return  a map of matching summary of gtypes
		 */
		public Map<String, Integer> unmaskSum(int start, int end) {
			Map<String, Integer> hits = new HashMap<String, Integer>(getNumGType());
			for(Map.Entry<String, BitSet> pair : gtypeIdx.entrySet()) {
				int sum = pair.getValue().get(start, end).cardinality();
				if(sum > 0)
					hits.put(pair.getKey(), sum); /* only summarize significant gtypes */
			}
			return hits;
		}
		
		/* member fields */
		private int length;  /* length of this BitMask */
		private Map<String, BitSet> gtypeIdx; /* per-gtype indices */
		
		/* class constants */
		public static final int DEFAULT_NTYPES = 16;
	}
	
	/* constructors */
	/** Default constructor, initiate index with default size */
	public GTypeIndex() {
		chrIdx = new HashMap<String, BitMask>();
	}
	
	/* member methods */
	public boolean hasChr(String chr) {
		return chrIdx.containsKey(chr);
	}
	
	public int getChrLen(String chr) {
		return chrIdx.get(chr).getLength();
	}
	
	/** add a new chrom with given length */
	public void addChr(String chr, int len) {
		chrIdx.put(chr, new BitMask(len));
	}
	
	/** remove a chrom from current index */
	public BitMask removeChr(String chr) {
		return chrIdx.remove(chr);
	}
	
	/**
	 * mask a chrom region as a given type
	 * @param chr  chrom
	 * @param start  0-based start
	 * @param end  1-based end
	 * @param gtype  genetic type
	 */
	public void maskRegion(String chr, int start, int end, String gtype) {
		BitMask bitMask = chrIdx.get(chr);
		if(bitMask != null && 0 <= start && start <= end && end <= bitMask.getLength()) { /* a valid region */
			bitMask.add(gtype);
			bitMask.setBit(gtype, start, end);
		}
	}

	/**
	 * unmask a given genome region
	 * @param chr  chrom
	 * @param start  0-based start
	 * @param end  1-based end
	 * @return  a set of matched gtypes, or null if this chrom does not exist
	 */
	public Set<String> unmask(String chr, int start, int end) {
		return chrIdx.containsKey(chr) ? chrIdx.get(chr).unmask(start, end) : null;
	}
	
	/**
	 * unmask a given genome region and get the total matched length of each gtype
	 * @param chr  chrom
	 * @param start  0-based start
	 * @param end  1-based end
	 * @return  a map of matched lengths of each gtype, or null if this chrom does not exist
	 */
	public Map<String, Integer> unmaskSum(String chr, int start, int end) {
		return chrIdx.containsKey(chr) ? chrIdx.get(chr).unmaskSum(start,  end) : null;
	}
	
	/* member fields */
	private Map<String, BitMask> chrIdx; /* per-chromosome BitMask */
}
