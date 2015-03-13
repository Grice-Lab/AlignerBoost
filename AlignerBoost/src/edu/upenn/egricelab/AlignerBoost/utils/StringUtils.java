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
 * 
 */
package edu.upenn.egricelab.AlignerBoost.utils;

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
