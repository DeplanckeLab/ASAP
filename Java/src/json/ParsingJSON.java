package json;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.TypeAdapter;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.google.gson.stream.JsonWriter;

import hdf5.loom.LoomData;
import hdf5.loom.LoomFile;
import model.Metadata;
import model.Parameters;
import parsing.model.FileType;

public class ParsingJSON
{
	public LoomFile loom;
	public LoomData data;
	public String message;
	
	public ParsingJSON() 
	{
		loom = new LoomFile("w", Parameters.outputFolder + "output.loom");
		data = new LoomData(Parameters.nGenes, Parameters.nCells);
		message = null;
	}
	
	public void writeOutputJSON()
	{
		writeOutputJSON(true);
	}
	
    public void writeOutputJSON(boolean doSize)
    {
    	int empty_columns = 0;
    	if(data.is_count_table)
    	{
	    	for(int i = 0; i < data.depth.size(); i++) 
	    	{
	    		if(data.depth.get(i) == 0) empty_columns++;
	    	}
    	}

    	if(doSize) data.meta.forEach((m)->m.size = loom.getSizeInBytes(m.path));
    	
    	// Prepare String
    	StringBuilder sb = new StringBuilder();
	    sb.append("{\"detected_format\":\"").append(Parameters.fileType).append("\",");
	    if(Parameters.fileType == FileType.LOOM) sb.append("\"loom_version\":\"").append(Parameters.loomVersion).append("\",");
    	if(message != null) sb.append("\"message\":\"").append(message).append("\",");
    	sb.append("\"nber_rows\":").append(data.nber_genes).append(",");
    	sb.append("\"nber_cols\":").append(data.nber_cells).append(",");
    	sb.append("\"nber_not_found_genes\":").append(data.nber_not_found_genes).append(",");
    	sb.append("\"nber_zeros\":").append(data.nber_zeros).append(",");
    	if(data.is_count_table && empty_columns != 0) sb.append("\"empty_columns\":").append(empty_columns).append(",");
    	sb.append("\"is_count_table\":").append(data.is_count_table?1:0);	
    	
    	// Handle Metadata generated
    	if(data.meta != null) sb.append(",").append(Metadata.toString(data.meta));
        
        // Handle eventual additional metadata
        if(data.existing_meta != null) sb.append(",").append(Metadata.toString(data.existing_meta, "existing_metadata"));
        
        // End
        sb.append("}");
        	
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
    
    public static ParsingJSON loadJSON(String jsonFile)
    {
    	ParsingJSON res = null;
		try
		{
			Gson gson = new GsonBuilder().registerTypeAdapter(Boolean.class, booleanAsIntAdapter).registerTypeAdapter(boolean.class, booleanAsIntAdapter).create();
			JsonReader reader = new JsonReader(new FileReader(jsonFile));
			res = gson.fromJson(reader, ParsingJSON.class); // contains the whole infos
			reader.close();
		}
		catch(FileNotFoundException nfe)
		{
			System.err.println("The JSON gene list was not found at the given path: " + Parameters.outputFolder + "output.json" + "\nStopping program...");
			System.exit(-1);
		}
		catch(Exception e)
		{
			System.out.println(e);
			System.err.println("Problem detected when reading the JSON gene list. Stopping program...");
			System.exit(-1);
		}
		return res;
    }
    
    private static final TypeAdapter<Boolean> booleanAsIntAdapter = new TypeAdapter<Boolean>() // This is only used to convert JSON int into booleans
    {
    	public void write(JsonWriter out, Boolean value) throws IOException 
    	{
    		if (value == null) out.nullValue();
    		else out.value(value);
    	}
    	  
    	public Boolean read(JsonReader in) throws IOException 
    	{
    		JsonToken peek = in.peek();
    		switch (peek) 
    		{
    			case BOOLEAN: return in.nextBoolean();
    			case NULL: in.nextNull(); return null;
    			case NUMBER: return in.nextInt() != 0;
    			case STRING: return Boolean.parseBoolean(in.nextString());
    			default: throw new IllegalStateException("Expected BOOLEAN or NUMBER but was " + peek);
    		}
    	}
	};
}
