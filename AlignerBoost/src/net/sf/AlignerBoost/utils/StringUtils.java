/**
 * 
 */
package net.sf.AlignerBoost.utils;

/**
 * @author Qi Zheng
 * @version  v1.2
 * @since  v1.2
 *
 */
public class StringUtils {

	/**
	 * Join an array of Strings using given separator
	 * @param sep  separator
	 * @param arr  array to be join
	 * @return  joined String, or null if arr is null
	 */
	public static String join(String sep, String[] arr) {
		if(arr == null)
			return null;
		StringBuilder str = new StringBuilder();
		for(int i = 0; i < arr.length; i++)
			str.append(i == 0 ? arr[i] : sep + arr[i]);
		return str.toString();
	}
	
	/**
	 * Join an array of Strings using space
	 * @param arr  array to be join
	 * @return  joined String, or null if arr is null
	 */
	public static String join(String[] arr) {
		return join(" ", arr);
	}

}
