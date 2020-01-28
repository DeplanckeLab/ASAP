package json;

import java.io.FileNotFoundException;
import java.io.FileReader;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.stream.JsonReader;

public class CopyMetaJSON 
{
	private String[] meta;
	
	public static String[] parseJSON(String jsonFile)
	{
		try
		{
			Gson gson = new GsonBuilder().create();
			JsonReader reader = new JsonReader(new FileReader(jsonFile));
			CopyMetaJSON obj = gson.fromJson(reader, CopyMetaJSON.class);
			reader.close();
			return obj.meta;
		}
		catch(FileNotFoundException nfe)
		{
			new ErrorJSON("The JSON metadata list was not found at the given path: " + jsonFile);
		}
		catch(Exception e)
		{
			new ErrorJSON("Problem detected when reading the JSON metadata list: " + e);
		}
		return null;
	}
}
