package json;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import model.Parameters;

public class ErrorJSON 
{
	public String displayed_error = "Error";
	
	public ErrorJSON(String errorMessage, String path) 
	{
		this.displayed_error = errorMessage;
		writeJSON(path);
		System.err.println(this.displayed_error);
		System.exit(-1);
	}
	
	public ErrorJSON(String errorMessage) 
	{
		this(errorMessage, Parameters.outputFolder + "output.json");
	}
    
    private void writeJSON(String path)
    {
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(path));
        	bw.write("{\"displayed_error\":\"" + displayed_error + "\"}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
}
