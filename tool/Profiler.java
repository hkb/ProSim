package tool;

import java.util.HashMap;

/**
 * Simple manual profiler that is meant to measure the usage of selected areas of the code.
 *
 * @author hkb
 */
public class Profiler {
	
	private HashMap<String,Long> time;
	private HashMap<String,Long> timer;
	private HashMap<String,Integer> counter;
	
	/**
	 * Create new profiler.
	 */
	public Profiler() {
		this.time = new HashMap<String,Long>();
		this.timer = new HashMap<String,Long>();
		this.counter = new HashMap<String,Integer>();
	}
	
	/**
	 * To be called when the execution enters a selected area.
	 * 
	 * @param area The name of the area.
	 * @require The execution must not already be in an area with the same name.
	 */
	public void enter(String area) {
		if(this.timer.containsKey(area)) {
			throw new IllegalArgumentException("Already timing " + area);
		}
		
		if(this.counter.containsKey(area)) {
			this.counter.put(area, this.counter.get(area)+1);
		} else {
			this.counter.put(area, 1);
		}
		
		this.timer.put(area, System.nanoTime());
	}

	/**
	 * To be called when the execution exits a selected area.
	 * 
	 * @param area The name of the area.
	 * @require The execution must be in an area with the same name.
	 */
	public void exit(String area) {
		long time = System.nanoTime();
		
		if(!this.timer.containsKey(area)) {
			throw new IllegalArgumentException("Isn't timing " + area);
		}
		
		time -= this.timer.get(area);
		
		if(this.time.containsKey(area)) {
			this.time.put(area, this.time.get(area) + time);
		} else {
			this.time.put(area, time);
		}
		
		this.timer.remove(area);
	}
	
	/**
	 * Clears the profiler.
	 */
	public void clear() {
		this.counter.clear();
		this.time.clear();
		this.timer.clear();
	}

	/**
	 * The current execution statistics.
	 * 
	 * @return A string with the execution statistics.
	 */
	public String stats() {
		StringBuilder output = new StringBuilder();
		
		for (String area : this.counter.keySet()) {
			output.append(area + ": ");
			output.append("count: " + this.counter.get(area));
			output.append(" total time: " + this.time.get(area)/1000000.0);
			output.append("ms relative time: " + (this.time.get(area)/1000000.0)/this.counter.get(area));
			output.append("ms\n");
		}
		
		return output.toString();
	}
}
