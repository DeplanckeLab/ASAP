package model;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

import bigarrays.StringArray64;
import tools.Utils;

public class Metadata 
{
	public String path;
	public Metatype type;
	public MetaOn on;
	public HashSet<String> categories;
	public HashMap<String, Long> categoriesMap;
	public StringArray64 values = null;
	public String[][] matrixValues = null; // TODO handle this if big?
	public long nbcol = -1;
	public long nbrow = -1;
	public long size = -1;
	
	public Metadata()
	{
		
	}
	
	public Metadata(String name, Metatype type, MetaOn on, long nbCol, long nbRow, StringArray64 values) 
	{
		this(name, type, on, nbCol, nbRow);
		this.values = values;
	}
	
	public Metadata(String name, Metatype type, MetaOn on, long nbCol, long nbRow, HashMap<String, Long> categoriesMap) 
	{
		this(name, type, on, nbCol, nbRow);
		this.categoriesMap = categoriesMap;
		this.categories = new HashSet<>();
		for(String cat:this.categoriesMap.keySet()) this.categories.add(cat);
	}
	
	public Metadata(String name, Metatype type, MetaOn on, long nbCol, long nbRow) 
	{
		this(name, type, on);
		this.nbcol = nbCol;
		this.nbrow = nbRow;
	}
	
	public Metadata(String name, Metatype type, MetaOn on) 
	{
		this(name);
		this.type = type;
		this.on = on;
	}
	
	public Metadata(String name) 
	{
		this();
		this.path = name;
	}
	
	public String toString()
	{
		return "\"name\":\""+path+"\",\"on\":\""+on+"\",\"type\":\""+type+"\"";
	}
	
	public boolean isCategorical() // Hard to find some kind of threshold for this
	{
		return isCategorical(this.values.size());
	}
	
	public boolean isCategorical(long length) // Hard to find some kind of threshold for this
	{
		if(this.categories == null) return false;
		return this.categories.size() <= length * 0.10;
	}
	
	public void inferType()
	{
		this.categories = new HashSet<>();
		this.categoriesMap = new HashMap<String, Long>();
		if(this.values == null) return;
		boolean isNumeric = true;
		for(long i = 0; i < this.values.size(); i++)
		{
			String v = this.values.get(i);
			if(isNumeric) // Check if numeric
			{
				try
				{
					Float.parseFloat(v);
				}
				catch(NumberFormatException nfe)
				{
					isNumeric = false;
				}
			}
			boolean added = this.categories.add(v);
			if(added) this.categoriesMap.put(v, 1L);
			else this.categoriesMap.put(v, this.categoriesMap.get(v) + 1);
			if(!isCategorical()) break; // To alleviate the calculations/ Map size
		}
		if(isCategorical()) this.type = Metatype.DISCRETE;
		else
		{
			this.categories = null;
			if(isNumeric) this.type = Metatype.NUMERIC;
			else this.type = Metatype.STRING;
		}
	}
	
	public static StringBuilder toString(Collection<Metadata> meta)
	{
		StringBuilder sb = new StringBuilder();
		String prefixT = "";
		sb.append("\"metadata\":[");
    	for(Metadata m:meta) 
    	{
    		sb.append(prefixT);
    		m.addMeta(sb, true, null, 10);
    		prefixT = ",";
    	}
    	sb.append("]");
		return sb;
	}
	
	public HashMap<String, Long> fillMap()
	{
		if(this.categoriesMap != null && !this.categoriesMap.isEmpty()) return this.categoriesMap;
		this.inferType();
		return this.categoriesMap;
	}
	
	public void addMeta(StringBuilder sb, boolean writeValues, StringArray64 cellNames, long max_nb_values)
	{
		sb.append("{");
		if(path != null) sb.append("\"name\":\"").append(this.path).append("\"");
		if(this.on != null) sb.append(",\"on\":\"").append(this.on).append("\"");
		if(this.type != null) sb.append(",\"type\":\"").append(this.type).append("\"");
		
		if(this.nbcol == -1) // Try to save the day?
		{
			if(this.values != null) // Vector
			{
				if(this.on == MetaOn.CELL) this.nbcol = this.values.size();
				else if(this.on == MetaOn.GENE) this.nbcol = 1;
			}
			else if(this.matrixValues != null)// Matrix
			{
				if(this.on == MetaOn.CELL) this.nbcol = this.matrixValues.length;
				else if(this.on == MetaOn.GENE) this.nbcol = this.matrixValues[0].length;
			}
		}
		
		if(this.nbrow == -1)// Try to save the day?
		{
			if(this.values != null) // Vector
			{
				if(this.on == MetaOn.GENE) this.nbrow = this.values.size();
				else if(this.on == MetaOn.CELL) this.nbrow = 1;
			}
			else if(this.matrixValues != null)// Matrix
			{
				if(this.on == MetaOn.GENE) this.nbrow = this.matrixValues.length;
				else if(this.on == MetaOn.CELL) this.nbrow = this.matrixValues[0].length;
			}
		}
		
		if(this.nbcol != -1) sb.append(",\"nber_cols\":").append(this.nbcol);
		if(this.nbrow != -1) sb.append(",\"nber_rows\":").append(this.nbrow);	
		if(this.size != -1) sb.append(",\"dataset_size\":").append(this.size);
		
		if(this.type == Metatype.DISCRETE)
		{	
			if(this.categoriesMap != null && !this.categoriesMap.isEmpty())
			{
				sb.append(",\"categories\":{");
				String prefix = "\"";
				for(String cat:this.categoriesMap.keySet()) 
				{
					sb.append(prefix);
					prefix = ",\"";
					sb.append(cat).append("\":").append(this.categoriesMap.get(cat));
				}
				sb.append("}");
			}
			else if(this.categories != null)
			{
				sb.append(",\"categories\":{");
				String prefix = "\"";
				for(String cat:this.categories) 
				{
					sb.append(prefix);
					prefix = ",\"";
					sb.append(cat).append("\":").append(-1);
				}
				sb.append("}");
			}
		}
		
		if((this.values != null || this.matrixValues != null) && writeValues)
		{
			sb.append(",\"values\":[");
			if(this.values != null) // If vector
			{
				String prefix = "\"";
				if(this.type == Metatype.NUMERIC) prefix = "";
				long cnt = 0; // To count the number of values to store
				
				for(long i = 0; i < this.values.size(); i++)
				{
					String val = this.values.get(i);
					if(this.type == Metatype.NUMERIC) val = Utils.format(val);
					sb.append(prefix).append(val);
					if(this.type != Metatype.NUMERIC) prefix = "\",\"";
					else prefix = ",";
					cnt++;
					if(cnt >= max_nb_values) break;
				}
				
				if(this.type != Metatype.NUMERIC) sb.append("\"");
			}
			else	// If matrix
			{
				String prefix2 = "";
				long cnt_y = 0; // To count the number of values to store by cols
				
				for(int j = 0; j < this.matrixValues[0].length; j++)
				{
	    			String prefix = "\"";
	    			if(this.type == Metatype.NUMERIC) prefix = "";
	    			long cnt_x = 0; // To count the number of values to store by rows
	    			
					sb.append(prefix2).append("[");
	    			for(int i = 0; i < this.matrixValues.length; i++)
	    			{
	    				String val = this.matrixValues[i][j];
	    				if(this.type == Metatype.NUMERIC) val = Utils.format(val);
	    				sb.append(prefix).append(val);
	    				if(this.type != Metatype.NUMERIC) prefix = "\",\"";
	    				else prefix = ",";
	    				cnt_x++;
	    				if(cnt_x >= max_nb_values) break;
					}
	    			if(this.type == Metatype.NUMERIC) sb.append("]");
	    			else sb.append("\"]");
	    			
	    			prefix2 = ",";
	    			cnt_y++;
	    			if(cnt_y >= max_nb_values) break;
				}
			}
			sb.append("]");
	
			addCellNames(sb, ",", cellNames, max_nb_values);
		}
		sb.append("}");
	}
	
	public static void addCellNames(StringBuilder sb, String pre, StringArray64 cellNames, long max_nb_values)
	{
		if(cellNames != null)
		{
			sb.append(pre).append("\"cells\":[");
			String prefix = "\"";
			long cnt = 0; // To count the number of values to store
			for(String val:cellNames) 
			{
				sb.append(prefix);
				prefix = "\",\"";
				sb.append(val);
				cnt++;
				if(cnt >= max_nb_values) break;
			}
			sb.append("\"]");
		}
	}
	
	public static void addGeneNames(StringBuilder sb, String pre, StringArray64 geneNames, long max_nb_values)
	{
		if(geneNames != null)
		{
			sb.append(pre).append("\"genes\":[");
			String prefix = "\"";
			long cnt = 0; // To count the number of values to store
			for(String val:geneNames) 
			{
				sb.append(prefix);
				prefix = "\",\"";
				sb.append(val);
				cnt++;
				if(cnt >= max_nb_values) break;
			}
			sb.append("\"]");
		}
	}
	
}