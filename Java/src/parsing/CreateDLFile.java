package parsing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import com.google.gson.Gson;
import com.google.gson.stream.JsonReader;

public class CreateDLFile 
{
	public static ArrayList<String> gene_names;
	
	public static void create(String inputMatrix, String outputFile) // JSON file
	{
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(inputMatrix));
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
			String line = br.readLine(); // header
			bw.write(line + "\n");
			line = br.readLine();
			while(line != null)
			{
				int indexTab = line.indexOf("\t");
				int id = Integer.parseInt(line.substring(0, indexTab)); // Should start at 0
				bw.write(gene_names.get(id) + "\t" + line.substring(indexTab+1) + "\n");
				line = br.readLine();
			}
			br.close();
			bw.close();
		}
		catch(IOException ioe)
		{
			ioe.printStackTrace();
		}
	}
	
	public static void getGenesFromJSON(String JSONFile) // JSON file
	{
		gene_names = new ArrayList<>();
		try
		{
			Gson gson = new Gson();
			JsonReader reader = new JsonReader(new FileReader(JSONFile));
			String[][] listGenes = gson.fromJson(reader, String[][].class); // contains the whole genes lists
			for (int i = 0; i < listGenes.length; i++) 
			{
				//String originalName = listGenes[i][0];
				//originalName = originalName.substring(0, originalName.indexOf(",")); // can be -1
				gene_names.add(listGenes[i][0] + "|" + listGenes[i][1] + "|" + listGenes[i][2]);
			}
			reader.close();
		}
		catch(FileNotFoundException nfe)
		{
			System.err.println("The JSON gene list was not found at the given path: " + JSONFile + "\nStopping program...");
			System.exit(-1);
		}
		catch(Exception e)
		{
			System.out.println(e);
			System.err.println("Problem detected when reading the JSON gene list. Stopping program...");
			System.exit(-1);
		}
	}
}

