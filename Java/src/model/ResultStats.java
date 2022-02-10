package model;

import tools.Utils;

public class ResultStats 
{
	public float min;
	public float q1;
	public float median;
	public float mean;
	public float q3;
	public float max;
	
	public void log2p()
	{
		this.min = (float)Utils.log2(1 + this.min);
		this.q1 = (float)Utils.log2(1 + this.q1);
		this.median = (float)Utils.log2(1 + this.median);
		this.mean = (float)Utils.log2(1 + this.mean);
		this.q3 = (float)Utils.log2(1 + this.q3);
		this.max = (float)Utils.log2(1 + this.max);
	}
}
