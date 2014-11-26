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
	 * construct a ProcessStatusTask with given initial status, total and info String
	 * @param processed   initial processed
	 * @param total  total status
	 * @param info  information String 
	 * @throws  an {@link IllegalArgumentException} if the total status is not positive integers 
	 */
	public PercentProcessStatusTask(long processed, long total, String info)
			throws IllegalArgumentException {
		if(total <= 0)
			throw new IllegalArgumentException("total must be a positive integer");
		this.processed = processed;
		this.total = total;
		this.info = info;
	}
	
	/**
	 * Overloaded construct a ProcessStatusTask with a given total and info String
	 * @param total  total status
	 * @param info  information String
	 * @throws  an {@link IllegalArgumentException} if the total status is not positive integers 
	 */
	public PercentProcessStatusTask(long total, String info)
			throws IllegalArgumentException {
		this(0, total, info);
	}

	/**
	 * Overloaded construct a ProcessStatusTask with a given total
	 * @param total  total status
	 * @throws  an {@link IllegalArgumentException} if the total status is not positive integers 
	 */
	public PercentProcessStatusTask(long total)
			throws IllegalArgumentException {
		this(total, "processed");
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
	 * re-implementation of the run method of a TimerTask
	 * @see java.util.TimerTask#run()
	 */
	@Override
	public void run() {
		if(processed > 0 && processed != prevProcessed) {
			System.err.printf("%.1f%% %s%n", 100.0 * processed / total, info);
			prevProcessed = processed;
		}
	}

	/**
	 * Show terminal status up-on finish this task
	 */
	public void finish() {
		System.err.printf("%.1f%% %s%n", 100.0 * processed / total, info);
	}
	
	private volatile long processed;
	private volatile long prevProcessed;
	private long total;
	private String info;
}
