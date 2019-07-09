package json;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.stream.JsonReader;

import hdf5.loom.LoomData;
import model.Metadata;
import model.Parameters;

public class FilterCellsJSON 
{
	// I know I have 3 different things that could be called the same (or put in 3 different classes), but Meh...
	public long[] selected_cells;
	public int[] kept_genes;
	public String[] filtered_cells;
	public long[] discarded_cols;
	
	public static FilterCellsJSON parseJSON(String jsonFile)
	{
		FilterCellsJSON res = null;
		try
		{
			Gson gson = new GsonBuilder().create();
			JsonReader reader = new JsonReader(new FileReader(jsonFile));
			res = gson.fromJson(reader, FilterCellsJSON.class);
			reader.close();
		}
		catch(FileNotFoundException nfe)
		{
			new ErrorJSON("The JSON cell list was not found at the given path: " + jsonFile);
		}
		catch(Exception e)
		{
			new ErrorJSON("Problem detected when reading the JSON cell list: "+e);
		}
		return res;
	}
	
	public static void writeOutputJSON(LoomData data)
    {
    	try
    	{
    		//data.meta.forEach((m)->m.fillMap());
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
    		bw.write("{\"detected_format\":\"LOOM\",");
    		bw.write("\"nber_rows\":" + data.nber_genes + ",");
        	bw.write("\"nber_cols\":" + data.nber_cells + ",");
        	bw.write("\"nber_zeros\":" + data.nber_zeros + ",");
        	bw.write("\"nber_ercc\":" + data.nber_ercc + ",");
         	bw.write("\"is_count_table\":" + (data.is_count_table?1:0) + ",");
         	bw.write(Metadata.toString(data.meta) + "}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		new ErrorJSON(ioe.getMessage());
    	}
    }
}
