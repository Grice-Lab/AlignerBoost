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
 * A class representing a Single-nucleotide polymorphism (SNP), also called as Single-nucleotide Variation (SNV)
 * Note that this SNP class also contains short in-dels 
 */
package edu.upenn.egricelab.AlignerBoost;

/**
 * @author Qi Zheng
 * @version 1.2
 * @since 1.2
 *
 */
public class SNP {
	
	/**
	 * A nested Enum member class for possible SNPTypes
	 * @author Qi Zheng
	 */
	static enum SNPType { SIMPLE, INSERTION, DELETION, MULTISUBSTITUTION };
	
	/**
	 * Construct a SNP with given info and determine the type automatically
	 * @param chr
	 * @param loc
	 * @param refAllele
	 * @param altAllele
	 */
	public SNP(String chr, int loc, String refAllele, String altAllele) {
		this.chr = chr;
		this.loc = loc;
		this.refAllele = refAllele.toUpperCase();
		this.altAllele = altAllele.toUpperCase();
		if(refAllele.length() == 1 && altAllele.length() == 1)
			type = SNPType.SIMPLE;
		else if(refAllele.length() < altAllele.length())
			type = SNPType.INSERTION;
		else if(refAllele.length() > altAllele.length())
			type = SNPType.DELETION;
		else
			type = SNPType.MULTISUBSTITUTION;
	}
	
	/**
	 * @return the chr
	 */
	public String getChr() {
		return chr;
	}

	/**
	 * @return the loc
	 */
	public int getLoc() {
		return loc;
	}

	/**
	 * @return the refAllele
	 */
	public String getRefAllele() {
		return refAllele;
	}

	/**
	 * @return the altAllele
	 */
	public String getAltAllele() {
		return altAllele;
	}

	/**
	 * @return the type
	 */
	public SNPType getType() {
		return type;
	}

	/** Simplify this SNP by reducing to shorted ref and alt alleles
	 * @return  this simplified SNP
	 */
	public SNP simplify() {
		if(isSimplified()) // already simplified
			return this;
		int refLen = refAllele.length();
		int altLen = altAllele.length();
		switch(type) {
		case INSERTION:
			int insLen = altLen - refLen;
			loc += refLen;
			refAllele = "";
			altAllele = altAllele.substring(refLen, refLen + insLen);
			break;
		case DELETION:
			int delLen = refLen - altLen;
			loc += altLen;
			refAllele = refAllele.substring(altLen, altLen + delLen);
			altAllele = "";
			break;
		case MULTISUBSTITUTION:
			loc++;
			refAllele = refAllele.substring(1);
			altAllele = altAllele.substring(1);
			break;
		default:
			break;
		}
		return this;
	}
	
	/**
	 * test whether this SNP has been simplified
	 * @return
	 */
	public boolean isSimplified() {
		return refAllele.length() == altAllele.length() || refAllele.isEmpty() || altAllele.isEmpty();
	}
	
	/**
	 * Get the MapLoc of this SNP
	 * @return  MapLoc for this SNP
	 */
	public String mapLoc() {
		return chr + ":" + loc;
	}
	
	/**
	 * Get the size of this SNP
	 * @return  1 if is SIMPLE, actual size if INSERTION or DELETION or MULTISUBSTITUTION
	 */
	public int size() {
		switch(type) {
		case SIMPLE:
			return 1;
		case INSERTION: case DELETION:
			return Math.abs(refAllele.length() - altAllele.length());
		case MULTISUBSTITUTION:
			return altAllele.length();
		default:
			return 0;
		}
	}

	/**
	 * Get the chrom covered lenth of this SNP
	 * @return  1 if is SIMPLE, actual size if DELETION or MULTISUBSTITUTION, 0 if INSERTION
	 */
	public int chrLen() {
		switch(type) {
		case SIMPLE:
			return 1;
		case INSERTION:
			return 0;
		case DELETION:
			return refAllele.length() - altAllele.length();
		case MULTISUBSTITUTION:
			return altAllele.length();
		default:
			return 0;
		}
	}
	
	/**
	 * Return the String representation of this SNP
	 */
	@Override
	public String toString() {
		return chr + ":" + loc + ":" + type + ":" + refAllele + "/" + altAllele;
	}

	/**
	 * test whether this SNP equals a given object
	 * @return  true only if thatObj is also a SNP and they have all same fields
	 */
	@Override
	public boolean equals(Object thatObj) {
		if(!(thatObj instanceof SNP))
			return false;
		SNP that = (SNP) thatObj;
		return chr.equals(that.chr) && loc == that.loc
				&& refAllele.equals(that.refAllele) && altAllele.equals(that.altAllele)
				&& type == that.type;
	}
	
	/**
	 * get the hashCode of this SNP
	 * @return  hashCode based on chr, loc, refAllele, altAllele and type. Equal SNPs will garantee to have equal hashCodes
	 */
	@Override
	public int hashCode() {
		return chr.hashCode() + loc + refAllele.hashCode() + altAllele.hashCode() + type.hashCode();
	}
	
	String chr;
	int loc;
	String refAllele = "";
	String altAllele = "";
	SNPType type;
}
