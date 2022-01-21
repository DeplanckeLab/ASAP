package config;

import java.io.InputStream;
import java.util.Properties;

import json.ErrorJSON;

public class Config 
{
	private Properties configFile;
	private String dbHost = "postgres:5434"; // Local access from the machine
	public static final String driver = "org.postgresql.Driver";
	
	private Config(String config_file, String dbHost)
	{
		this(config_file);
		this.dbHost = dbHost;
	}
	
	private Config(String config_file)
	{
		this.configFile = new java.util.Properties();
		try 
		{
			InputStream is = Config.class.getClassLoader().getResourceAsStream(config_file);
			if(is == null) is = Config.class.getClassLoader().getResourceAsStream("config/" + config_file);
			this.configFile.load(is);
		}
		catch(Exception eta)
		{
		    new ErrorJSON(eta.getMessage());
		}
	}
	
	public static Config ConfigMAIN()
	{
		return new Config("asap.conf");
	}
	
	public static Config ConfigDEV()
	{
		return new Config("asap.conf", "asap.epfl.ch:5433"); // Access via intranet - No internet access
	}
		
	private String getValue(String key) 
	{
		return configFile.getProperty(key);
	}
	
	public String getProperty(String key)
	{
		return this.getValue(key);
	}
	
	public String getURL(String dbName)
	{
		return "jdbc:postgresql://" + this.dbHost + "/" + dbName + "?user=" + this.getProperty("mDbUser") + "&password=" + this.getProperty("mDbPwds");
	}
	
	public String getURLFromHost(String dbHostName)
	{
		return "jdbc:postgresql://" + dbHostName + "?user=" + this.getProperty("mDbUser") + "&password=" + this.getProperty("mDbPwds");
	}
}
