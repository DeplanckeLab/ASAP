package model;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.TypeAdapter;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.google.gson.stream.JsonWriter;

class Batch
{
	public int nber_lines_parsed;
	public int nber_lines_unparsed;
	public int nber_groups;
	
	@Override
	public String toString() {
		return "{\"nber_lines_parsed\":" + nber_lines_parsed + ",\"nber_lines_unparsed\":" + nber_lines_unparsed + ",\"nber_groups\":" + nber_groups + "}";
	}
}

public class OutputJSON
{
	public int nber_cells = 0;
	public int nber_genes = 0;
	public int nber_not_found_genes = 0;
	public int nber_ercc = 0;
	public boolean is_count_table = true;
	public long nber_zeros = 0;
	public int nber_duplicated_genes = 0;
	public int nber_unique_genes = 0;
	public int nber_all_duplicated_genes = 0;
	public Batch batch_file = null;
	public long nber_total_biotypes = 0;
	public long bio_protein_coding = 0;
	public long bio_rRNA = 0;
	public long total_chrs = 0;
	public long chr_MT = 0;
	public HashMap<String, Integer> biotypes = new HashMap<>();
	public HashMap<String, Integer> chrs = new HashMap<>();
	
	public void setBiotypes()
	{
		for(String bio:biotypes.keySet()) nber_total_biotypes += biotypes.get(bio);
		bio_protein_coding = (biotypes.get("protein_coding") == null)?0:biotypes.get("protein_coding");
		bio_rRNA = ((biotypes.get("rRNA") == null)?0:biotypes.get("rRNA")) + ((biotypes.get("Mt_rRNA") == null)?0:biotypes.get("Mt_rRNA")) + ((biotypes.get("rRNA_pseudogene") == null)?0:biotypes.get("rRNA_pseudogene"));
		for(String chr:chrs.keySet()) total_chrs += chrs.get(chr);
		chr_MT = ((chrs.get("MT") == null)?0:chrs.get("MT"));
	}
	
    public void writeOutputJSON()
    {
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
        	bw.write("{\"nber_genes\":" + nber_genes + ",");
        	bw.write("\"nber_cells\":" + nber_cells + ",");
        	bw.write("\"nber_not_found_genes\":" + nber_not_found_genes + ",");
        	bw.write("\"nber_duplicated_genes\":" + nber_duplicated_genes + ",");
        	bw.write("\"nber_all_duplicated_genes\":" + nber_all_duplicated_genes + ",");
        	bw.write("\"nber_zeros\":" + nber_zeros + ",");
        	bw.write("\"nber_ercc\":" + nber_ercc + ",");
        	bw.write("\"nber_unique_genes\":" + nber_unique_genes + ",");
        	if(is_count_table)
        	{
	        	bw.write("\"nber_total_biotypes\":" + nber_total_biotypes + ",");
	        	bw.write("\"nber_protein_coding\":" + bio_protein_coding + ",");
	        	bw.write("\"nber_rRNA\":" + bio_rRNA + ",");
	        	bw.write("\"nber_total_chr\":" + total_chrs + ",");
	        	bw.write("\"nber_MT\":" + chr_MT + ",");
        	}
        	bw.write("\"is_count_table\":" + (is_count_table?1:0) + ",");
        	bw.write("\"batch_file\":" + batch_file + "}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
    
    public static OutputJSON loadOutputJSON()
    {
    	OutputJSON res = null;
		try
		{
			Gson gson = new GsonBuilder().registerTypeAdapter(Boolean.class, booleanAsIntAdapter).registerTypeAdapter(boolean.class, booleanAsIntAdapter).create();
			JsonReader reader = new JsonReader(new FileReader(Parameters.outputFolder + "output.json"));
			res = gson.fromJson(reader, OutputJSON.class); // contains the whole infos
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
