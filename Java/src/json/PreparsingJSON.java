package json;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import model.Parameters;
import parsing.model.GroupPreparse;
import parsing.model.MetaPreparse;

public class PreparsingJSON
{
    public static void writeOutputJSON(MetaPreparse meta)
    {
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
        	bw.write("{\"detected_format\":\"" + Parameters.fileType + "\"," + meta + "}");   	
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
    
	
    public static void writeOutputJSON(ArrayList<GroupPreparse> groups)
    {
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
        	bw.write("{\"detected_format\":\"" + Parameters.fileType + "\",\"list_groups\":" + toString(groups) + "}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
    
    public static void writeListingJSON(ArrayList<String> files)
    {
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
        	bw.write("{\"detected_format\":\"" + Parameters.fileType + "\",\"list_files\":[");
        	boolean putComma = false;
       		for(String f:files)
        	{
       			if(putComma) bw.write(",");
       			bw.write("{\"filename\":\"" + f + "\"}");
       			putComma = true;
        	}
       		bw.write("]}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
     	
   	public static String toString(ArrayList<GroupPreparse> stuff)
   	{
   		StringBuffer sb = new StringBuffer("[");
   		String prefix0 = "";
   		for(GroupPreparse g:stuff)
    	{
    		sb.append(prefix0).append("{\"group\":\"").append(g.name).append("\",\"nber_cols\":").append(g.nbCells).append(",\"nber_rows\":").append(g.nbGenes).append(",\"is_count\":").append(g.isCount?1:0); 
    		if(g.geneNames != null && g.geneNames.length != 0) 
    		{ 
    			sb.append(",\"genes\":[");
    			String prefix = "";
    			for(String gene:g.geneNames) { sb.append(prefix).append("\"").append(gene).append("\""); prefix = ","; }
    			sb.append("]");
    		}
    		if(g.cellNames != null && g.cellNames.length != 0) 
    		{ 
    			sb.append(",\"cells\":[");
    			String prefix = "";
    			for(String cell:g.cellNames) { sb.append(prefix).append("\"").append(cell).append("\""); prefix = ","; }
    			sb.append("]");
    		}
    		if(g.matrix != null)
    		{
    			sb.append(",\"matrix\":[");
    			String prefix1 = "";
    			for(int i = 0; i< g.matrix.length; i++)
    			{
    				sb.append(prefix1).append("[");
    				String prefix2 = "";
       				if(g.matrix[i].length != 0) for(int j = 0; j < g.matrix[i].length; j++) { sb.append(prefix2).append(g.matrix[i][j]); prefix2 = ","; }
     				sb.append("]");
    				prefix1 = ",";
    			}
    			sb.append("]");
    		}
    		if(g.additionalMetadataPath != null && g.additionalMetadataPath.size() > 0)
    		{
    			sb.append(",\"existing_metadata\":[");
    			String prefix = "";
    			for(String meta:g.additionalMetadataPath) { sb.append(prefix).append("\"").append(meta).append("\""); prefix = ","; }
    			sb.append("]");
    		}
    		sb.append("}");
    		prefix0 = ",";
    	}
   		sb.append("]");
    	return sb.toString();
    }
}
