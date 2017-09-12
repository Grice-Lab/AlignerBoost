/**
 * 
 */
package edu.upenn.egricelab.ucsc;

import java.util.*;
import java.util.regex.*;

/**
 * @author zhengqi
 *
 */
public class GTF extends GFF {
	/* constructors */
	/** default constructor, do nothing */
	public GTF() {  }
	
	/** construct a GTF record from String, using its super class constructor */
	public GTF(String line) throws IllegalArgumentException, ArrayIndexOutOfBoundsException
	{
		super(line);
	}

	/** read attributes into a map from a String */
	@Override
	public void readAttributes(String attrStr) {
		attrNames = new ArrayList<String>();
		attrMap = new HashMap<String, String>();
		for(String field : attrStr.split(sep)) {
			Matcher match = attrPat.matcher(field);
			while(match.find())
				setAttr(match.group(1), match.group(2));
		}
	}

	/** write attributes in a given order */
	@Override
	public String writeAttributes() {
		if(attrNames == null || attrMap == null)
			return "";
		StringBuilder attrStr = new StringBuilder();
		for(int i = 0; i < attrNames.size(); i++) {
			if(i > 0)
				attrStr.append(sep);
			attrStr.append(attrNames.get(i) + " " + quote + getAttr(attrNames.get(i)) + quote);
		}
		return attrStr.toString();
	}
	
	/* class constants */
	public static final String quote = "\"";
	public static final String sep = "; ";
	public static final Pattern attrPat = Pattern.compile("(\\w+) \"([^\"]+)\"");
}
