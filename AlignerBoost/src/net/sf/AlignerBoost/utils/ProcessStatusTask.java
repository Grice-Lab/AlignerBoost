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
	}

	/**
	 * re-implementation of the run method of a TimerTask
	 * @see java.util.TimerTask#run()
	 */
	@Override
	public void run() {
		if(status > 0 && status != prevStatus) {
			System.err.println(status + " " + info);
			prevStatus = status;
		}
	}

	/**
	 * Show terminal status up-on finish this task
	 */
	public void finish() {
		System.err.println("Total " + status + " " + info);
	}

	private static final String DEFAULT_INFO = "processed"; // default display info
	private volatile long status;
	private volatile long prevStatus;
	private String info;
}
