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
 * A class for Enviroment variables, versions and other apps constants
 */
package edu.upenn.egricelab.AlignerBoost;

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 *
 */
public class EnvConstants {
	// prog versions
	public static final String progName = "AlignerBoost";
	public static final String progVer = "v1.3.4";
	public static final String progDesc = "A tool for boosting the precision and sensitivity of NextGen-seq aligners";
	public static final String progFile = progName + ".jar";
	// environment variables
	public static final String newLine = System.getProperty("line.separator", "\n");
	public static final String osName = System.getProperty("os.name");
}
