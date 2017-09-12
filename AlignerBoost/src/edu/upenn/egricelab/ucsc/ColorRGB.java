/**
 * 
 */
package edu.upenn.egricelab.ucsc;

/**
 * A RGB color class used by many UCSC formats
 * @author zhengqi
 * @version v1.1
 */
public class ColorRGB {
	/* constructor */
	/** default constructor, do nothing */
	public ColorRGB() {  }

	/**
	 * construct RGB from given values
	 * @throws IllegalArgumentException if R/G/B values are not in [0~255]
	 */
	public ColorRGB(int red, int green, int blue) throws IllegalArgumentException
	{
		init(red, green, blue);
	}
	
	public ColorRGB(String rgb) throws IllegalArgumentException, ArrayIndexOutOfBoundsException
	{
		String[] rgbFields = rgb.split(sep, numFields);
		init(Integer.parseInt(rgbFields[0]), Integer.parseInt(rgbFields[1]), Integer.parseInt(rgbFields[2]));
	}

	public int red;
	public int green;
	public int blue;

	@Override
	public boolean equals(Object o) {
		if(!(o != null && o instanceof ColorRGB))
			return false;
		ColorRGB other = (ColorRGB) o;
		return red == other.red && green == other.green && blue == other.blue;
	}
	
	@Override
	public String toString() {
		return red + sep + green + sep + blue;
	}
	
	/* helper methods */
	protected void init(int red, int green, int blue) throws IllegalArgumentException
	{
		if(!(0 <= red && red <= MAX_VALUE && 0 <= green && green <= MAX_VALUE && 0 <= blue && blue <= MAX_VALUE))
			throw new IllegalArgumentException("R/G/B values must be in [0, " + MAX_VALUE + "]");
		this.red = red;
		this.green = green;
		this.blue = blue;
	}
	
	/* class constants */
	public static final int MAX_VALUE = 255;
	public static final String sep = ",";
	public static final int numFields = 3;
}
