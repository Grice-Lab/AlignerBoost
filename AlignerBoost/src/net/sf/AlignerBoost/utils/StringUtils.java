/**
 * 
 */
package net.sf.AlignerBoost.utils;

import java.util.Map;

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

	/**
	 * A generic method to join a Map of any objects using : and given separator
	 * @param  map any map to be joined to string
	 * @param  sep String separator
	 * @return  joined string as "key:value[<sep>key:value]", or null if map is null
	 */
	public static <K, V> String join(Map<K, V> map, String sep) {
		if(map == null)
			return null;
		StringBuilder str = new StringBuilder();
		for(Map.Entry<K, V> entry : map.entrySet())
			str.append(str.length() == 0 ? entry.getKey() + ":" + entry.getValue() : sep + entry.getKey() + ":" + entry.getValue());
		return str.toString();
	}

	/**
	 * A generic method to join a Map of any objects using : and default separator ','
	 * @param  map any map to be joined to string
	 * @return  joined string as "key:value[,key:value]", or null if map is null
	 */
	public static <K, V> String join(Map<K, V> map) {
		return join(map, ",");
	}
}
