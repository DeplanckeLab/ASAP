package json;

import java.util.ArrayList;

public class WarningJSON 
{
	private static ArrayList<String> warnings = new ArrayList<String>();
	
	public static void addWarning(String warningMessage) 
	{
		warnings.add(warningMessage);
	}
	
	public static boolean isAnyWarning()
	{
		return !warnings.isEmpty();
	}
	
    public static String getJSON()
    {
    	StringBuilder sb = new StringBuilder();
    	sb.append("\"warnings\":[");
    	String prefix = "";
    	for(String warning:warnings) 
    	{
    		sb.append(prefix).append("\"").append(warning).append("\"");
    		prefix = ",";
    	}
    	sb.append("]");
    	return sb.toString();
    }
}
