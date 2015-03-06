/**
 * A class for Enviroment variables, versions and other apps constants
 */
package net.sf.AlignerBoost;

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 *
 */
public class EnvConstants {
	// prog versions
	public static final String progName = "AlignerBoost";
	public static final String progVer = "v1.3.2";
	public static final String progDesc = "A tool for boosting the accuracy of NextGen-seq aligner";
	public static final String progFile = "AlignerBoost.jar";
	// environment variables
	public static final String newLine = System.getProperty("line.separator", "\n");
	public static final String osName = System.getProperty("os.name");
}
