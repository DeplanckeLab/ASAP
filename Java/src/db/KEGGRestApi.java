package db;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;

import json.ErrorJSON;
import model.Parameters;
import tools.Utils;

public class KEGGRestApi 
{
	public static HashMap<String, String> keggDescription = null;
	public static HashMap<String, ArrayList<String>> keggAnnot = null;
	public static HashMap<String, String> keggIDtoEnsembl = null;
	public static HashMap<String, String> keggIDtoNCBI = null; // Entrez
	
	public static void generateKEGGDB() 
	{
		
		HashMap<Integer, String> keggDB = fetchKEGGOrganisms();
		System.out.println(keggDB.size() + " organisms are found in KEGG DB");
	
		DBManager.connect();
		HashMap<Integer, Integer> taxonsInDB = DBManager.listAllOrganisms();
		System.out.println("There are " + taxonsInDB.size() + " organisms in ASAP DB.");
		DBManager.disconnect();
		
		Object[] keys = keggDB.keySet().toArray();
		for(Object taxon:keys) if(!taxonsInDB.containsKey(taxon)) keggDB.remove(taxon);
		System.out.println(taxonsInDB.size() + " taxons are shared between KEGG and ASAP DBs");
		
		// Processing all taxons
		for(int taxon:keggDB.keySet())
		{
			long t = System.currentTimeMillis();
			int ASAPId = taxonsInDB.get(taxon);
			String keggShortName = keggDB.get(taxon);
			System.out.println("\nProcessing " + keggShortName + " : " + taxon + " (" + ASAPId + ") ...");
			new File(Parameters.outputFolder + "kegg_id/").mkdirs();
			new File(Parameters.outputFolder + "ensembl_id/").mkdirs();
			new File(Parameters.outputFolder + "ncbi_id/").mkdirs();
			generateKEGGDB(keggShortName, taxon);
			System.out.println("Processed in " + Utils.toReadableTime(System.currentTimeMillis() - t));
		}	
	}
	
	public static void generateKEGGDB(String keggShortName, int taxon)
	{
		keggDescription = new HashMap<String, String>();
		keggAnnot = new HashMap<String, ArrayList<String>> ();
		keggIDtoNCBI = new HashMap<String, String>();
		keggIDtoEnsembl = new HashMap<String, String>();
		fetchKEGGPathways(keggShortName);
		System.out.println(keggDescription.size() + " KEGG pathways were found.");
		for(String keggID:keggDescription.keySet()) fetchGenesInPathway(keggShortName, keggID);
		
		// Create GMT with KEGG ids
		try
		{
			BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "kegg_id/kegg."+taxon + ".gmt"));
			for(String geneset:keggAnnot.keySet())
			{
				ArrayList<String> genes = keggAnnot.get(geneset);
				bw.write(geneset+"\t"+keggDescription.get(geneset)+"\thttp://www.genome.jp/dbget-bin/www_bget?"+geneset);
				for (String gene:genes) bw.write("\t"+gene); // List gene
				bw.write("\n");
			}
			bw.close();
		}
		catch(IOException ioe)
		{
			new ErrorJSON(ioe.getMessage());
		}
		
		// Fetch Gene names from KEGG ids
		geneNameConversionToNCBI(keggShortName);

		// Create GMT from REST data
		try
		{
			BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "ncbi_id/kegg."+taxon + ".gmt"));
			for(String geneset:keggAnnot.keySet())
			{
				ArrayList<String> genes = keggAnnot.get(geneset);
				bw.write(geneset+"\t"+keggDescription.get(geneset)+"\thttp://www.genome.jp/dbget-bin/www_bget?"+geneset);
				for (String gene:genes) bw.write("\t"+keggIDtoNCBI.get(gene)); // List gene
				bw.write("\n");
			}
			bw.close();
		}
		catch(IOException ioe)
		{
			new ErrorJSON(ioe.getMessage());
		}
		
		// Fetch Gene names from KEGG ids
		/*geneNameConversionToEnsembl(keggShortName);

		// Create GMT from REST data
		try
		{
			BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "ensembl_id/kegg."+taxon + ".gmt"));
			for(String geneset:keggAnnot.keySet())
			{
				ArrayList<String> genes = keggAnnot.get(geneset);
				bw.write(geneset+"\t"+keggDescription.get(geneset)+"\thttp://www.genome.jp/dbget-bin/www_bget?"+geneset);
				for (String gene:genes) bw.write("\t"+keggIDtoEnsembl.get(gene)); // List gene
				bw.write("\n");
			}
			bw.close();
		}
		catch(IOException ioe)
		{
			new ErrorJSON(ioe.getMessage());
		}*/
	}
	
    public static void geneNameConversionToNCBI(String keggShortName)
    {
    	long t = System.currentTimeMillis();
    	System.out.println(keggIDtoNCBI.size() + " KEGG Genes to transform to Entrez IDs...");
    	
    	try
    	{
		    URL url = new URL("http://rest.genome.jp/link/ncbi-geneid/"+keggShortName);
		    URLConnection urlc = url.openConnection();
		    urlc.setDoOutput(true);
		    urlc.setAllowUserInteraction(false);
	
		    //get result
		    BufferedReader br = new BufferedReader(new InputStreamReader(urlc.getInputStream()));
		    String l = null;
		    while((l=br.readLine())!=null) 
		   	{
		    	String[] results = l.split("\t");
		    	if(results.length >= 2)
		    	{
		    		keggIDtoNCBI.put(results[0], results[1].substring("ncbi-geneid:".length()));
		    	}
		 	}
		  	br.close();
    	}
    	catch(IOException ioe)
    	{
    		new ErrorJSON(ioe.getMessage());
    	}
    	
    	System.out.println("Conversion took " + Utils.toReadableTime(System.currentTimeMillis() - t));
    }
    
    public static void geneNameConversionToEnsembl(String keggShortName)
    {
    	long t = System.currentTimeMillis();
    	
    	System.out.println(keggIDtoEnsembl.size() + " KEGG Genes to transform to Ensembl...");
    	try
    	{
    		for(String keggID:keggIDtoEnsembl.keySet())
    		{
			    URL url = new URL("https://www.kegg.jp/entry/"+keggID);
			    URLConnection urlc = url.openConnection();
			    urlc.setDoOutput(true);
			    urlc.setAllowUserInteraction(false);
		
			    //get result
			    BufferedReader br = new BufferedReader(new InputStreamReader(urlc.getInputStream()));
			    String l = null;
			    while((l=br.readLine())!=null) 
			   	{
			    	int index = l.indexOf("Ensembl");
			    	if(index != -1)
			    	{
			    		System.out.println(index);
			    		//keggIDtoEnsembl.put(results[0], results[1].substring("ncbi-geneid:".length()));
			    	}
			 	}
			  	br.close();
    		}
    	}
    	catch(IOException ioe)
    	{
    		new ErrorJSON(ioe.getMessage());
    	}
    	
    	System.out.println("Conversion took " + Utils.toReadableTime(System.currentTimeMillis() - t));
    }

    public static void fetchGenesInPathway(String organism, String pathway)
    {
    	try
    	{
	        URL url = new URL("http://rest.kegg.jp/link/"+organism+"/" + pathway);
	        URLConnection urlc = url.openConnection();
	        urlc.setDoOutput(true);
	        urlc.setAllowUserInteraction(false);
	
	        //get result
	        BufferedReader br = new BufferedReader(new InputStreamReader(urlc.getInputStream()));
	        String l = null;
	        ArrayList<String> genes = new ArrayList<String>();
	        while ((l=br.readLine())!=null) 
	        {
	        	String name = l.split("\t")[1];
	        	keggIDtoNCBI.put(name, null);
	        	keggIDtoEnsembl.put(name, null);
	        	genes.add(name);
	        }
	        keggAnnot.put(pathway, genes);
	        br.close();
    	}
    	catch(IOException ioe)
    	{
    		new ErrorJSON(ioe.getMessage());
    	}
    }
    
    public static void fetchKEGGPathways(String organism)
    {
    	try
    	{
	        URL url = new URL("http://rest.kegg.jp/list/pathway/" + organism);
	        URLConnection urlc = url.openConnection();
	        urlc.setDoOutput(true);
	        urlc.setAllowUserInteraction(false);
	
	        //get result
	        BufferedReader br = new BufferedReader(new InputStreamReader(urlc.getInputStream()));
	        String l = null;
	        while ((l=br.readLine())!=null) 
	        {
	        	String[] results = l.split("\t");
	        	keggDescription.put(results[0].substring(5), results[1].substring(0, results[1].lastIndexOf(" -")));
	        }
	        br.close();
    	}
    	catch(IOException ioe)
    	{
    		new ErrorJSON(ioe.getMessage());
    	}
    }
    
    public static HashMap<Integer, String> fetchKEGGOrganisms()
    {
    	HashMap<Integer, String> map = new HashMap<>();
    	
    	try
    	{
	        URL url = new URL("http://rest.kegg.jp/list/genome");
	        URLConnection urlc = url.openConnection();
	        urlc.setDoOutput(true);
	        urlc.setAllowUserInteraction(false);
	
	        //get result
	        BufferedReader br = new BufferedReader(new InputStreamReader(urlc.getInputStream()));
	        String l = null;
	        while ((l=br.readLine())!=null) 
	        {
	        	String[] results = l.split("\t")[1].split(",|;");
	        	if(results.length > 1)
	        	{
		        	Integer taxon = null;
		        	try
		        	{
		        		taxon = Integer.parseInt(results[1].trim());
		        	}
		        	catch(NumberFormatException nfe)
		        	{
		        		taxon = Integer.parseInt(results[2].trim());
		        	}
		        	map.put(taxon, results[0]);
	        	}
	        }
	        br.close();
    	}
    	catch(IOException ioe)
    	{
    		new ErrorJSON(ioe.getMessage());
    	}
        
        return map;
    }
}
