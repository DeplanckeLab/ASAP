package json;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import model.Metadata;
import model.Parameters;
import parsing.model.FileType;
import parsing.model.GroupPreparse;
import tools.Utils;

public class PreparsingJSON
{
    /** 
     * Main function creating the JSON of the Preparse steps
     * @param groups
     */
    public static void writeOutputJSON(ArrayList<GroupPreparse> groups)
    {
    	// Prepare String
    	StringBuilder sb = new StringBuilder();
	    sb.append("{\"detected_format\":\"").append(Parameters.fileType).append("\",");
	    if(Parameters.fileType == FileType.LOOM) sb.append("\"loom_version\":\"").append(Parameters.loomVersion).append("\",");
    	sb.append("\"list_groups\":").append(toString(groups)).append("}");
    	
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
    
    /** 
     * Creating the JSON of the Preparsing step for H5AD files
     * @param groups Matrices found in file (X, raw.X, raw/X)
     * @param geneMetadata Gene Metadata found in file (/var and /varm)
     * @param cellMetadata Cell Metadata found in file (/obs and /obsm)
     * @param ignoredMetadata Gene & Cell Metadata found in file but ignored for parsing (/varp /obsp and others)
     */
    public static void writeH5ADOutputJSON(ArrayList<GroupPreparse> groups, List<Metadata> okMetadata, List<Metadata> otherMetadata)
    {
    	// Prepare String
    	StringBuilder sb = new StringBuilder();
	    sb.append("{\"detected_format\":\"").append(Parameters.fileType).append("\",");
    	sb.append("\"list_groups\":").append(toString(groups)).append(",").append(Metadata.toString(okMetadata, "metadata")).append(",").append(Metadata.toString(otherMetadata, "other_metadata")).append("}");
    	
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
   
    /**
     * If archive containing multiple files
     * First: Just list the files
     * 
     * @param files
     */
    public static void writeListingJSON(ArrayList<String> files)
    {
    	// Prepare String
    	StringBuilder sb = new StringBuilder();
    	sb.append("{\"detected_format\":\"").append(Parameters.fileType).append("\",");
	    if(Parameters.fileType == FileType.LOOM) sb.append("\"loom_version\":\"").append(Parameters.loomVersion).append("\",");
    	sb.append("\"list_files\":[");
    	String prefix = "";
   		for(String f:files)
    	{
   			sb.append(prefix).append("{\"filename\":\"").append(f).append("\"}");
   			prefix = ",";
    	}
   		sb.append("]}");
    	
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
     	
   	public static String toString(ArrayList<GroupPreparse> stuff)
   	{
   		StringBuilder sb = new StringBuilder("[");
   		String prefix0 = "";
   		for(GroupPreparse g:stuff)
    	{
    		sb.append(prefix0).append("{\"group\":\"").append(g.name).append("\",\"nber_cols\":").append(g.nbCells).append(",\"nber_rows\":").append(g.nbGenes).append(",\"is_count\":").append(g.isCount?1:0); 
    		if(g.geneNames != null && g.geneNames.length != 0) 
    		{ 
    			sb.append(",\"genes\":[");
    			String prefix = "";
    			for(String gene:g.geneNames) { sb.append(prefix).append("\"").append(Utils.handleSpecialCharacters(gene)).append("\""); prefix = ","; }
    			sb.append("]");
    		}
    		if(g.cellNames != null && g.cellNames.length != 0) 
    		{ 
    			sb.append(",\"cells\":[");
    			String prefix = "";
    			for(String cell:g.cellNames) { sb.append(prefix).append("\"").append(Utils.handleSpecialCharacters(cell)).append("\""); prefix = ","; }
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
    			for(String meta:g.additionalMetadataPath) { sb.append(prefix).append("\"").append(Utils.handleSpecialCharacters(meta)).append("\""); prefix = ","; }
    			sb.append("]");
    		}
    		sb.append("}");
    		prefix0 = ",";
    	}
   		sb.append("]");
    	return sb.toString();
    }
}
