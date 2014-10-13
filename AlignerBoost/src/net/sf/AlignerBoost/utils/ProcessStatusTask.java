/**
 * a TimerTask idea for use of displaying a process status in stderr showing the # of processed items,
 * it should be run inside a Timer or separate thread
 */
package net.sf.AlignerBoost.utils;

import java.util.TimerTask;

/**
 * @author Qi Zheng
 * @version 1.1
 * @since 1.1
 */
public class ProcessStatusTask extends TimerTask {

	/**
	 * construct a ProcessStatusTask with given initial status
	 * @param status  initial status
	 */
	public ProcessStatusTask(long status) {
		this.status = status;
	}

	/**
	 * default constructor to set the initial status to 0
	 */
	public ProcessStatusTask() { } // Default constructor, do nothing

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
	 * re-implementation of the run method of a TimerTask
	 * @see java.util.TimerTask#run()
	 */
	@Override
	public void run() {
		System.err.println(status + " processed");
	}

	private volatile long status;

	/**
	 * Show terminal status up-on finish this task
	 */
	public void finish() {
		System.err.println("Total " + status + " processed successfully");
	}

}
