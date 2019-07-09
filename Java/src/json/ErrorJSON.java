package json;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import model.Parameters;

public class ErrorJSON 
{
	public String displayed_error = "Error";
	
	public ErrorJSON(String errorMessage) 
	{
		this.displayed_error = errorMessage;
		writeJSON();
		System.err.println(this.displayed_error);
		System.exit(-1);
	}
	
	public ErrorJSON(String errorMessage, String additionalOptions) 
	{
		this.displayed_error = errorMessage;
		writeJSON(additionalOptions);
		System.err.println(this.displayed_error);
		System.exit(-1);
	}
	
    private void writeJSON()
    {
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
        	bw.write("{\"displayed_error\":\"" + displayed_error + "\"}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
    
    private void writeJSON(String additionalOptions)
    {
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
        	bw.write("{\"displayed_error\":\"" + displayed_error + "\","+ additionalOptions +"}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
}
