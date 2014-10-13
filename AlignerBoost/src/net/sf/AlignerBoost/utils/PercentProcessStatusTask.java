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
	 * construct a ProcessStatusTask with given initial status
	 * @param status  initial status
	 * @throws  an {@link IllegalArgumentException} if the total status is not positive integers 
	 */
	public PercentProcessStatusTask(long total) throws IllegalArgumentException {
		if(total <= 0)
			throw new IllegalArgumentException("total must be a positive integer");
		this.total = total;
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
	 * re-implementation of the run method of a TimerTask
	 * @see java.util.TimerTask#run()
	 */
	@Override
	public void run() {
		System.err.printf("%.1f%% processed%n", processed / (double) total * 100);
	}

	/**
	 * Show terminal status up-on finish this task
	 */
	public void finish() {
		System.err.printf("%.1f%% processed%n", processed / (double) total * 100);
	}
	
	private volatile long processed;
	private long total;

}
