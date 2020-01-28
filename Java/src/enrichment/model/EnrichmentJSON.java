package enrichment.model;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.stream.JsonReader;

import hdf5.loom.LoomData;
import json.ErrorJSON;
import model.Metadata;
import model.Parameters;

public class EnrichmentJSON
{
	// I know I have 3 different things that could be called the same (or put in 3 different classes), but Meh...
	public long[] up;
	public long[] down;
	public long[] indexes_match;
	
	public static EnrichmentJSON parseJSON(String jsonFile)
	{
		EnrichmentJSON res = null;
		try
		{
			Gson gson = new GsonBuilder().create();
			JsonReader reader = new JsonReader(new FileReader(jsonFile));
			res = gson.fromJson(reader, EnrichmentJSON.class);
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
		// Prepare output JSON
		StringBuffer sb = new StringBuffer();
		sb.append("{\"detected_format\":\"LOOM\",");	
		sb.append("\"loom_version\":\"").append(Parameters.loomVersion).append("\",");
		sb.append("\"nber_rows\":").append(data.nber_genes).append(",");
		sb.append("\"nber_cols\":").append(data.nber_cells).append(",");
		sb.append("\"nber_zeros\":").append(data.nber_zeros).append(",");
		sb.append("\"is_count_table\":").append(data.is_count_table?1:0).append(",");
		sb.append(Metadata.toString(data.meta)).append("}");
		
		// Write output JSON
        try
        {
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
    		bw.write(sb.toString());
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		new ErrorJSON(ioe.getMessage());
    	}
    }
}
