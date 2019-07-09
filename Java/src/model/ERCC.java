package model;

import java.util.HashMap;

public class ERCC 
{
	private float[][] values = null; // Counts
	private String[] names = null; // Names of ERCCs
	private boolean flag = false;
	
	private HashMap<Long, Integer> ercc_indexes = new HashMap<Long, Integer>(); // Keep tracks of which indexes are ERCCs
	private int localIndex;
	
	public boolean isERCC(long index)
	{
		if(ercc_indexes.get(index) != null) return true;
		return false;
	}
	
	public void callERCC(long index, String name)
	{
		this.ercc_indexes.put(index, this.localIndex);
		this.names[this.localIndex] = name;
		this.localIndex++;
	}
	
	public ERCC(int nbCells) 
	{
		this.values = new float[92][nbCells]; // Create array of 92 at the beginning (should be enough in most cases)
		this.names = new String[92];
		this.localIndex = 0;
	}
	
	public void set(long erccIndex, long cellIndex, float value) 
	{
		int localI = this.ercc_indexes.get(erccIndex);
		if(localI > values.length - 1) // If current array size is not enough anymore...
		{
			float[][] valuesTmp = new float[values.length * 2][this.values[0].length];
			String[] namesTmp = new String[values.length * 2];
			for(int k = 0; k < values.length; k++) 
			{
				namesTmp[k] = names[k];
				for(int l = 0; l < values[k].length; l++) valuesTmp[k][l] = values[k][l];
			}
			this.values = valuesTmp;
			this.names = namesTmp;
		}
		this.values[localI][(int)cellIndex] = value; // TODO Handle big arrays
	}
	
	public void set(int i, int j, float value, String name) 
	{
		if(i > values.length - 1) // If current array size is not enough anymore...
		{
			float[][] valuesTmp = new float[values.length * 2][this.values[0].length];
			String[] namesTmp = new String[values.length * 2];
			for(int k = 0; k < values.length; k++) 
			{
				namesTmp[k] = names[k];
				for(int l = 0; l < values[k].length; l++) valuesTmp[k][l] = values[k][l];
			}
			this.values = valuesTmp;
			this.names = namesTmp;
		}
		this.values[i][j] = value;
		this.names[i] = name;
	}
	
	private void prepareWriting() // Remove extra rows
	{
		flag = true;
		int maxLength = 0;
		
		for(int k = 0; k < names.length; k++) 
		{
			if(names[k] == null) break;
			maxLength = k+1;
		}
		
		float[][] valuesTmp = new float[maxLength][values[0].length];
		String[] namesTmp = new String[maxLength];
		
		for(int k = 0; k < maxLength; k++) 
		{
			namesTmp[k] = names[k];
			for(int l = 0; l < values[k].length; l++) valuesTmp[k][l] = values[k][l];
		}
		this.values = valuesTmp;
		this.names = namesTmp;
	}
	
	public float[][] getArray()
	{
		if(!flag) prepareWriting();
		return values;
	}
	
	public String[] getNames()
	{
		if(!flag) prepareWriting();
		return names;
	}
}
