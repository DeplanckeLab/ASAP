package json;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import hdf5.loom.LoomData;
import model.Parameters;

public class FilterJSON
{
	public long nber_filtered_genes = 0;
	public ArrayList<DataPlots> list_plots = new ArrayList<DataPlots>();
	public ArrayList<DataWarnings> warnings = new ArrayList<DataWarnings>();
	public String info = "";
	
    public void writeJSON(LoomData data)
    {
    	if(data.nber_genes == nber_filtered_genes) new ErrorJSON("The output matrix has no more genes. You should put less stringent thresholds.");
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
        	bw.write("{\"nber_rows\":" + (data.nber_genes - nber_filtered_genes) + ",");
        	bw.write("\"nber_cols\":" + data.nber_cells + ",");
        	bw.write("\"nber_filtered_genes\":" + nber_filtered_genes + ",");
        	bw.write("\"nber_zeros\":" + data.nber_zeros + ",");
        	if(!info.equals("")) bw.write("\"info\":\"" + info + "\",");
        	bw.write("\"list_plots\":" + DataPlots.toString(list_plots) + ",");
        	bw.write("\"warnings\":" + DataWarnings.toString(warnings));
        	bw.write("}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
    
    public void addPlot(String name, String description)
    {
    	DataPlots d = new DataPlots();
    	d.description = description;
    	d.name = name;
    	list_plots.add(d);
    }
    
    public void addWarning(String message)
    {
    	DataWarnings w = new DataWarnings();
    	w.message = message;
    	warnings.add(w);
    }
}

class DataPlots
{
	String name;
	String description;
	
	public static String toString(ArrayList<DataPlots> plots)
	{
		String res = "[";
		for(DataPlots p:plots)
		{
			res += "{\"name\":\"" + p.name + "\",\"description\":\"" + p.description + "\"},"; 
		}
		if(!res.equals("[")) res = res.substring(0, res.length() - 1);
		return res + "]";
	}
}

class DataWarnings
{
	String message;
	
	public static String toString(ArrayList<DataWarnings> warnings)
	{
		String res = "[";
		for(DataWarnings w:warnings)
		{
			res += "{\"message\":\"" + w.message + "\"},"; 
		}
		if(!res.equals("[")) res = res.substring(0, res.length() - 1);
		return res + "]";
	}
}