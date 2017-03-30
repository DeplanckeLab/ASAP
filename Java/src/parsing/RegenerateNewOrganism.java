package parsing;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import com.google.gson.Gson;
import com.google.gson.stream.JsonReader;

import db.DBManager;
import model.OutputJSON;
import model.Parameters;
import parsing.model.Gene;
import tools.Utils;

public class RegenerateNewOrganism 
{
	public static OutputJSON json = null;
	public static HashMap<String, ArrayList<Gene>> dbGenes = null;
	public static BufferedWriter bw_NF = null;
	public static ArrayList<String> genes = null;
	public static HashMap<String, Integer> geneDups = null;
	
	public static void regenerateJSON()
	{
		long t = System.currentTimeMillis();
		dbGenes = DBManager.getGenesInDB();
		System.out.println("Accessing DB time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
		t = System.currentTimeMillis();
		if(RegenerateNewOrganism.dbGenes.size() > 0)
		{
			try
			{
				geneDups = new HashMap<String, Integer>();
				json = OutputJSON.loadOutputJSON();
				System.out.println("There was " + json.nber_not_found_genes + " not found genes in output.json.");
				json.nber_not_found_genes = 0;
				System.out.println("There was " + json.nber_duplicated_genes + " duplicated genes [unique] in output.json.");
				json.nber_duplicated_genes = 0;
				System.out.println("There was " + json.nber_all_duplicated_genes + " duplicated genes in output.json.");
				json.nber_all_duplicated_genes = 0;
				genes = getGenesFromJSON();
				System.out.println("Find " + genes.size() + " genes in JSON.");
				bw_NF = new BufferedWriter(new FileWriter(Parameters.outputFolder + "not_found_genes.txt"));
				for(int i = 0; i < genes.size(); i++) genes.set(i, generateName(genes.get(i))); // Replace the element according to what is found in DB
				bw_NF.close();
				System.out.println("Now, there is " + json.nber_not_found_genes + " not found genes in output.json.");
				System.out.println("Now, there is " + json.nber_duplicated_genes + " duplicated genes [unique] in output.json.");
				System.out.println("Now, there is " + json.nber_all_duplicated_genes + " duplicated genes in output.json.");
				json.nber_unique_genes = geneDups.size();
				for(String gene:geneDups.keySet()) if(geneDups.get(gene) != 1) json.nber_duplicated_genes++;
				json.writeOutputJSON();
				writeGeneNamesJSON();
				writeDuplicatedGenes(geneDups);
			}
			catch(IOException ioe)
			{
				ioe.printStackTrace();
			}
		}
		else System.err.println("No gene found in database for this organism. Aborted.");
	}
	
    public static void writeDuplicatedGenes(HashMap<String, Integer> geneDups)
    {
    	try
    	{
        	BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "duplicated_genes.txt"));
        	for(String geneKey:geneDups.keySet())
        	{
        		Integer nb = geneDups.get(geneKey);
        		if(nb > 1) bw.write(geneKey + "\t" + nb + "\n");
        	}
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
	
	public static String generateName(String gene) throws IOException
	{
		ArrayList<Gene> dbHit = dbGenes.get(gene.toUpperCase());
		String ensIdList = "";
		String geneIdList = "";
		if(dbHit == null) 
		{
			bw_NF.write(gene + "\n"); // I cannot regenerate this perfectly since I don't know the line it was found in
			json.nber_not_found_genes++;
		}
		else
		{
			for(Gene gHit:dbHit) // There can be multiple time the same name
			{
				ensIdList += gHit.ensembl_id + ",";
				geneIdList += gHit.name + ",";
			}
			ensIdList = ensIdList.substring(0, ensIdList.length()-1); // remove last "," => not pretty but Meh
			geneIdList = geneIdList.substring(0, geneIdList.length()-1); // remove last "," => not pretty but Meh
		}
		String res = "[\"" + gene + "\",\"" + ensIdList + "\",\"" + geneIdList + "\"]";
		Integer count = geneDups.get(res);
		if(count == null) geneDups.put(res, 1);// First time I see this gene 
		else
		{
			count++;
			json.nber_all_duplicated_genes++;
			geneDups.put(res, count);
			res = "[\"" + gene + "," + count + "\",\"" + ensIdList + "\",\"" + geneIdList + "\"]"; // Make the name unique
		}
		return res;
	}
	
    public static void writeGeneNamesJSON()
    {
    	try
    	{
        	BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "gene_names.json"));
        	if(genes.size() == 0) bw.write("{}");
        	else
        	{
	        	bw.write("[" + genes.get(0));
	        	for(int i = 1; i < genes.size(); i++) bw.write("," + genes.get(i)); // Already in JSON format
	        	bw.write("]");
        	}
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    		System.exit(-1);
    	}
    }
    
	public static ArrayList<String> getGenesFromJSON()
	{
		ArrayList<String> res = new ArrayList<>();
		try
		{
			Gson gson = new Gson();
			JsonReader reader = new JsonReader(new FileReader(Parameters.JSONFileName));
			String[][] listGenes = gson.fromJson(reader, String[][].class); // contains the whole genes lists
			for (int i = 0; i < listGenes.length; i++) 
			{
				String originalName = listGenes[i][0];
				int index = originalName.indexOf(",");
				if(index != -1) originalName = originalName.substring(0, index);
				res.add(originalName);
			}
			reader.close();
		}
		catch(FileNotFoundException nfe)
		{
			System.err.println("The JSON gene list was not found at the given path: " + Parameters.JSONFileName + "\nStopping program...");
			System.exit(-1);
		}
		catch(Exception e)
		{
			System.out.println(e);
			System.err.println("Problem detected when reading the JSON gene list. Stopping program...");
			System.exit(-1);
		}
		return res;
	}
}
