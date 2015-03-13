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
 * a TimerTask idea for use of displaying a process status in stderr showing the # of processed items,
 * it should be run inside a Timer or separate thread
 */
package edu.upenn.egricelab.AlignerBoost.utils;

import java.util.TimerTask;

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class ProcessStatusTask extends TimerTask {

	/**
	 * Overloaded construct a ProcessStatusTask with given info String
	 * @param info  information String
	 */
	public ProcessStatusTask(String info) {
		this.info = info;
	}
	
	/**
	 * Overloaded default constructor to set the initial status to 0
	 */
	public ProcessStatusTask() { // Default constructor
		this(DEFAULT_INFO);
	}

	/**
	 * get status
	 * @return  status
	 */
	public long getStatus() {
		return status;
	}

	/**
	 * update the status by increment the status count
	 */
	public void updateStatus() {
		status++;
	}

	/**
	 * @return the info
	 */
	public String getInfo() {
		return info;
	}

	/**
	 * @param info  the display info to set
	 */
	public void setInfo(String info) {
		this.info = info;
	}
	
	/**
	 * reset this ProcessStatusTask to initial status
	 */
	public void reset() {
		status = 0;
		info = DEFAULT_INFO;
		isFinished = false;
	}

	/**
	 * re-implementation of the run method of a TimerTask
	 * @see java.util.TimerTask#run()
	 */
	@Override
	public void run() {
		if(!isFinished && status > 0 && status != prevStatus) {
			System.err.println(status + " " + info);
			prevStatus = status;
		}
	}

	/**
	 * Show terminal status up-on finish this task
	 */
	public void finish() {
		isFinished = true;
		System.err.println("Total " + status + " " + info);
	}

	private static final String DEFAULT_INFO = "processed"; // default display info
	private volatile long status;
	private volatile long prevStatus;
	private volatile boolean isFinished;
	private String info;
}
