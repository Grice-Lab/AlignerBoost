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
 * a TimerTask idea for use of displaying a process status in stderr showing the % of process of a task,
 * it should be run inside a Timer or separate thread
 */
package net.sf.AlignerBoost.utils;

import java.util.TimerTask;

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class PercentProcessStatusTask extends TimerTask {
	/**
	 * Overloaded construct a ProcessStatusTask with a given total and info String
	 * @param total  total status
	 * @param info  information String
	 * @throws  an {@link IllegalArgumentException} if the total status is not positive integers 
	 */
	public PercentProcessStatusTask(long total, String info)
			throws IllegalArgumentException {
		if(total < 0)
			throw new IllegalArgumentException("total must be positive");
		this.total = total;
		this.info = info;
	}

	/**
	 * Overloaded construct a ProcessStatusTask with a given total
	 * @param total  total status
	 * @throws  an {@link IllegalArgumentException} if the total status is not positive integers 
	 */
	public PercentProcessStatusTask(long total)
			throws IllegalArgumentException {
		this(total, DEFAULT_INFO);
	}
	
	/**
	 * get number of processed tasks
	 * @return  processed
	 */
	public long getProcessed() {
		return processed;
	}

	/**
	 * update the status by increment the status count
	 */
	public void update() {
		processed++;
	}

	/**
	 * @return the info
	 */
	public String getInfo() {
		return info;
	}

	/**
	 * @param info the info to set
	 */
	public void setInfo(String info) {
		this.info = info;
	}

	/**
	 * Reset this PercentProcessTask to initial status
	 */
	public void reset() {
		processed = 0;
		info = DEFAULT_INFO;
		isFinished = false;
	}
	
	/**
	 * re-implementation of the run method of a TimerTask
	 * @see java.util.TimerTask#run()
	 */
	@Override
	public void run() {
		if(!isFinished && processed > 0 && processed != prevProcessed) {
			System.err.printf("%.1f%% %s%n", 100.0 * processed / total, info);
			prevProcessed = processed;
		}
	}

	/**
	 * Show terminal status up-on finish this task
	 */
	public void finish() {
		isFinished = true;
		System.err.printf("Total %.1f%% %s%n", 100.0 * processed / total, info);
	}
	
	private static final String DEFAULT_INFO = "processed";
	private volatile long processed;
	private volatile long prevProcessed;
	private volatile boolean isFinished;
	private long total;
	private String info;
}
