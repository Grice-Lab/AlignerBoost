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
 * A utility class with simple Statistic methods
 */
package edu.upenn.egricelab.AlignerBoost.utils;

/**
 * @author Qi Zheng
 * @version 1.2
 * @since 1.2
 */
public class Stats {

	/** determine whether a string is a decimal integer
	 * @return true if is a valid decimal integer
	 */
	public static boolean isDecimal(String s) {
		if(s.isEmpty())
			return false;
		for(int i = 0; i < s.length(); i++) {
			if(i == 0 && s.charAt(i) == '-') {
				if(s.length() == 1)
					return false;
				else
					continue;
			}
			if(!Character.isDigit(s.charAt(i)))
				return false;
		}
		return true;
	}

	/**
	 * Transfer phred q-value to p-value
	 * @param q  q-value in log10-scale
	 * @param scale  phred scale
	 * @return  p-value
	 */
	public static double phredQ2P(double q, double scale) {
		if(q < 0)
			q = 0;
		return Math.pow(10.0, q / -scale);
	}

	/**
	 * Transfer phred q-value to p-value in -10 scale
	 * @param q  q-value in log10-scale
	 * @return  p-value
	 */
	public static double phredQ2P(double q) {
		return phredQ2P(q, PHRED_SCALE);
	}

	/**
	 * Transfer phred p-value to q-value
	 * @param p  p-value
	 * @param scale  phred scale
	 * @return  q-value in log10-scale
	 */
	public static double phredP2Q(double p, double scale) {
		return -scale * Math.log10(p);
	}

	/**
	 * Transfer phred p-value to q-value
	 * @param p  p-value
	 * @return  q-value in log10-scale
	 */
	public static double phredP2Q(double p) {
		return phredP2Q(p, PHRED_SCALE);
	}
	
	/**
	 * Get the mean value of a byte array of given region
	 * @param array  byte array
	 * @param start  0-based start index
	 * @param end  1-based end index
	 * @return  mean value in double
	 */
	public static double mean(byte[] array, int start, int end) {
		assert(start <= end);
		if(array == null)
			return Double.NaN;
		double x = 0;
		for(int i = start; i < end; i++)
			x += array[i];
		return x / (end - start);
	}
	
	/**
	 * Get the mean value of a byte array
	 * @param array  byte array
	 * @return  mean value in double
	 */
	public static double mean(byte[] array) {
		return mean(array, 0, array.length);
	}

	//	public static final int MAX_QUAL = 255; // max mapQ
	/**
	 * Get the mean value of a byte array
	 * @param array  byte array
	 * @param start  0-based start
	 * @param end  1-based end
	 * @return  mean value in double
	 */
	public static double mean(int[] array, int start, int end) {
		if(array == null || array.length == 0)
			return Double.NaN;
		assert end > start;
		double x = 0;
		for(int i = start; i < end; i++)
			x += array[i];
		return x / (end - start);
	}

	//	public static final int MAX_QUAL = 255; // max mapQ
	/**
	 * Get the mean value of a byte array
	 * @param array  byte array
	 * @param start  0-based start
	 * @param end  1-based end
	 * @return  mean value in double
	 */
	public static double mean(float[] array, int start, int end) {
		if(array == null || array.length == 0)
			return Double.NaN;
		assert end > start;
		double x = 0;
		for(int i = start; i < end; i++)
			x += array[i];
		return x / (end - start);
	}

	public static final double PHRED_SCALE = 10; // scaling factor for phred scores
	public static final int ASCII_OFFSET = 33;
}
