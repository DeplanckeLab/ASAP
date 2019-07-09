package config;

import java.io.InputStream;
import java.util.Properties;

import json.ErrorJSON;
import model.Parameters;

public class Config 
{
	private Properties configFile;
	private static Config instance;
	
	private Config()
	{
		configFile = new java.util.Properties();
		try 
		{
			InputStream is = Config.class.getClassLoader().getResourceAsStream(Parameters.configFile);
			if(is == null) is = Config.class.getClassLoader().getResourceAsStream("config/" + Parameters.configFile);
			configFile.load(is);
		}
		catch(Exception eta)
		{
		    new ErrorJSON(eta.getMessage());
		}
	}

	private String getValue(String key) 
	{
		return configFile.getProperty(key);
	}
	
	public static String getProperty(String key)
	{
		if (instance == null) instance = new Config();
		return instance.getValue(key);
	}
}
