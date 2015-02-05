/**
 * A utility class with simple Statistic methods
 */
package net.sf.AlignerBoost.utils;

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
	 * Get the mean value of a byte array
	 * @param array  byte array
	 * @return  mean value in double
	 */
	public static double mean(byte[] array) {
		if(array == null || array.length == 0)
			return Double.NaN;
		double x = 0;
		for(byte b : array)
			x += b;
		return x / array.length;
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
