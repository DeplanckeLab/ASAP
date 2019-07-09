package json;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import model.Parameters;

public class DimReducJSON 
{
	public ArrayList<DataWarnings> warnings = new ArrayList<DataWarnings>();
	
    public void writeJSON()
    {
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
        	bw.write("{\"text\":" + " TODO " + ",");
        	bw.write("{\"PC1\":" + " TODO " + ",");
        	bw.write("{\"PC2\":" + " TODO " + ",");
        	bw.write("{\"PC3\":" + " TODO " + ",");
        	bw.write("{\"PC4\":" + " TODO " + ",");
        	bw.write("{\"PC5\":" + " TODO " + ",");
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
        
    public void addWarning(String message)
    {
    	DataWarnings w = new DataWarnings();
    	w.message = message;
    	warnings.add(w);
    }
}
