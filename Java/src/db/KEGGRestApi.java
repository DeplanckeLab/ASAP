package db;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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
	public static HashMap<String, String> keggIDtoGene = null;
	
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
			generateKEGGDB(keggShortName, Parameters.outputFolder + "kegg."+taxon+".gmt");
	    	System.out.println("Processed in " + Utils.toReadableTime(System.currentTimeMillis() - t));
		}
	}
	
	public static void generateKEGGDB(String keggShortName, String outputGMTFile)
	{
		keggDescription = new HashMap<String, String>();
		keggAnnot = new HashMap<String, ArrayList<String>> ();
		keggIDtoGene = new HashMap<String, String>();
		fetchKEGGPathways(keggShortName);
		System.out.println(keggDescription.size() + " KEGG pathways were found.");
		for(String keggID:keggDescription.keySet()) fetchGenesInPathway(keggShortName, keggID);
		
		// Fetch Gene names from KEGG ids
		geneNameConversion();
		
		// Create GMT from REST data
		try
		{
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputGMTFile));
			for(String geneset:keggAnnot.keySet())
			{
				ArrayList<String> genes = keggAnnot.get(geneset);
				bw.write(geneset+"\t"+keggDescription.get(geneset)+"\thttp://www.genome.jp/dbget-bin/www_bget?"+geneset);
				for (String gene:genes) bw.write("\t"+keggIDtoGene.get(gene)); // List gene
				bw.write("\n");
			}
			bw.close();
		}
		catch(IOException ioe)
		{
			new ErrorJSON(ioe.getMessage());
		}
	}
	
    public static void geneNameConversion()
    {
    	System.out.println(keggIDtoGene.size() + " KEGG Genes to transform to Entrez IDs...");
    	
    	String[] keggIDs = keggIDtoGene.keySet().toArray(new String[] {});
    	int i = 0;

    	// Process 100 by 100
    	while(i < keggIDs.length)
    	{
	    	StringBuilder sb = new StringBuilder("http://rest.kegg.jp/conv/ncbi-geneid/");
	    	String prefix = "";
	    	do
	    	{ 
	    		sb.append(prefix).append(keggIDs[i]); prefix = "+"; 
	    		i++;
	    	} while(i % 100 != 0 && i < keggIDs.length);
	    	
	    	try
	    	{
			    URL url = new URL(sb.toString());
			    URLConnection urlc = url.openConnection();
			    urlc.setDoOutput(true);
			    urlc.setAllowUserInteraction(false);
		
			    //get result
			    BufferedReader br = new BufferedReader(new InputStreamReader(urlc.getInputStream()));
			    String l = null;
			    while((l=br.readLine())!=null) 
			   	{
			    	String[] results = l.split("\t");
			    	keggIDtoGene.put(results[0], results[1].substring("ncbi-geneid:".length()));
			 	}
			  	br.close();
	    	}
	    	catch(IOException ioe)
	    	{
	    		new ErrorJSON(ioe.getMessage());
	    	}
    	}
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
	        	keggIDtoGene.put(name, null);
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
