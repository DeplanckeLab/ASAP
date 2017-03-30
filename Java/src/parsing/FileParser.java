package parsing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import model.OutputJSON;
import model.Parameters;
import parsing.model.Gene;

public class FileParser
{
	public static OutputJSON json = null;
	public static String header = "";
	public static HashMap<String, ArrayList<Gene>> dbGenes = null;
	
    public static void parse()
    {
    	json = new OutputJSON();
    	HashSet<String> ERCCs = new HashSet<String>();
    	ArrayList<String> genes = new ArrayList<String>();
    	HashMap<String, Integer> geneDups = new HashMap<String, Integer>();
    	ArrayList<String> cell_names = new ArrayList<>();
        System.out.println("Parsing file : " + Parameters.fileName);
        try
        {
        	int current_line = 0;
        	int num_columns = -1;
        	int num_columns_header = -1;
        	boolean headerWritten = false;
            BufferedReader br = null;
            if(Parameters.fileName.endsWith(".gz")) br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(Parameters.fileName)))); 
            else if(Parameters.fileName.endsWith(".zip")) 
            {
            	ZipInputStream zis = new ZipInputStream(new FileInputStream(Parameters.fileName));
            	zis.getNextEntry();
            	br = new BufferedReader(new InputStreamReader(zis));
            }
            else br = new BufferedReader(new FileReader(Parameters.fileName));
        	BufferedWriter bw_tab = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.tab"));
        	BufferedWriter bw_dl_tab = new BufferedWriter(new FileWriter(Parameters.outputFolder + "dl_output.tab"));
        	BufferedWriter bw_ercc = new BufferedWriter(new FileWriter(Parameters.outputFolder + "ercc.tab"));
        	BufferedWriter bw_NF = new BufferedWriter(new FileWriter(Parameters.outputFolder + "not_found_genes.txt"));
        	BufferedWriter bw_exonLength = new BufferedWriter(new FileWriter(Parameters.outputFolder + "sum_exon_length.txt"));
        	String line = br.readLine();
        	
        	// Skip the first lines if needed
        	while(current_line < Parameters.skip_line) { line = br.readLine(); current_line++; } // Skipping the first "skip_line" lines
        	
        	// Read header if exist
        	if(Parameters.has_header) // Parsing Header
        	{
        		String[] tokens = line.split(Parameters.delimiter);
        		for(String t:tokens)
        		{
        			String name = t.trim();
        			if(name.startsWith("\"") || name.startsWith("'")) name = name.substring(1);
        			if(name.endsWith("\"") || name.endsWith("'")) name = name.substring(0, name.length() - 1);
        			if(!name.equals("")) cell_names.add(name);
        		}
        		if(cell_names.size() == 0) Parameters.has_header = false; // Header is empty => not taken into account
        		else num_columns_header = cell_names.size();
        		line = br.readLine();
        	}
        	
        	// Read the rest of the file
            while(line != null)
            {
            	current_line++; // To shift the count from [0, n-1] TO [1, n]
            	String[] tokens = line.split(Parameters.delimiter);         	
            	ArrayList<String> rowValues = new ArrayList<String>();
            	int col = 0;
            	
            	// First store the values as Strings into the array (remove empty values and format)
        		for(String t:tokens)
        		{
        			String value = t.trim();
        			if(value.startsWith("\"") || value.startsWith("'")) value = value.substring(1);
        			if(value.endsWith("\"") || value.endsWith("'")) value = value.substring(0, value.length() - 1);
        			if(value.equals("NA")) error("Value at line " + current_line+ " col " + col + " is NA, but this is not allowed in the current version of ASAP.");
        			if(!value.equals("")) rowValues.add(value);
        			col++;
        		}
        		int l = rowValues.size();
        		
        		// Then parse the values stored and check consistency with number of columns
        		if(rowValues.size() != 0) // if the line is not empty
        		{
        			if(num_columns == -1) num_columns = l;
        			if(num_columns != l) error("Row " + current_line+ " contains a different number of values (" + rowValues.size() + ") than the other rows (" + num_columns + ")");
        			if(num_columns_header != -1 && num_columns != num_columns_header && num_columns != (num_columns_header + 1)) error("Row " + current_line+ " contains a different number of values (" + rowValues.size() + ") than the header (" + num_columns_header + ")");
        			if(!headerWritten) { createHeader(cell_names, num_columns, num_columns_header); bw_tab.write(header + "\n"); bw_dl_tab.write(header + "\n"); bw_ercc.write(header + "\n"); headerWritten = true; } // Only called once, after first line is parsed
        			
        			// If I am here, it means that the number of values is correct => I grab the gene name if exists
        			int start = 0;
        			boolean isERCC = false;
        			int end = num_columns;
        			String gene = null;
					String biotype = null;
					String chr = null;
    				switch(Parameters.name_column)
    				{
    					case FIRST:
    						gene = rowValues.get(0);
    						start = 1;
    						break;
    					case LAST:
    						gene = rowValues.get(end - 1);
    						end = num_columns - 1;
    						break;
    					case NONE:
    						gene = "Gene_" + (json.nber_genes + 1); // Starts at 1 is more pretty
    				}
    				if(gene.startsWith("ERCC-")) // This is recognized as an ERCC
    				{
    					if(ERCCs.contains(gene)) error("Duplicated ERCCs are not allowed (at line " + current_line + ")");
    					ERCCs.add(gene);
    					bw_ercc.write(gene);
    					isERCC = true;
    					json.nber_ercc++;
    				}
    				else
    				{
    					ArrayList<Gene> dbHit = dbGenes.get(gene.toUpperCase());
    					String ensIdList = "";
    					String geneIdList = "";
    					long sumExonLength = -1;
    					if(dbHit == null) 
    					{
    						bw_NF.write(gene + "\t" + current_line + "\n");
    						json.nber_not_found_genes++;
    					}
    					else
    					{
	    					for(Gene gHit:dbHit) // There can be multiple time the same name
	    					{
	    						ensIdList += gHit.ensembl_id + ",";
	    						geneIdList += gHit.name + ",";
	    						biotype = gHit.biotype; // if there are multiple biotypes across annotations, I keep one at random...
	    						chr = gHit.chr; // if there are multiple CHR across annotations, I keep one at random...
	    						sumExonLength = gHit.sum_exon_length; // if there are multiple CHR across annotations, I keep one at random...
	    					}
	    					ensIdList = ensIdList.substring(0, ensIdList.length()-1); // remove last "," => not pretty but Meh
	    					geneIdList = geneIdList.substring(0, geneIdList.length()-1); // remove last "," => not pretty but Meh
    					}
    					String res = gene + "|" + ensIdList + "|" + geneIdList;
    					Integer count = geneDups.get(res);
    					if(count == null) geneDups.put(res, 1);// First time I see this gene 
    					else
    					{
    						count++;
    						json.nber_all_duplicated_genes++;
    						geneDups.put(res, count);
    						res = gene + "," + count + "|" + ensIdList + "|" + geneIdList; // Make the name unique
    					}
    					genes.add(res); // Add unique name
    					bw_exonLength.write(res+"\t"+sumExonLength+"\n");
    					bw_dl_tab.write(res);
    					bw_tab.write(""+json.nber_genes); // This is the ID number
    					json.nber_genes++;
    				}
    				
    				// Then I look at the values and write them
        			for(int i = start; i < end; i++)
        			{
        				String value = rowValues.get(i);
	    				try
	    				{
	        				double v = Double.parseDouble(value);
	        				double rv = Math.round(v);
	        				if(Math.abs(rv - v) < 1E-5) // if the diff between the value and the rounded is lower than 1E-5 I round it (to avoid Java bugs with float representation)
	        				{
	        					v = rv;
	        					if(!isERCC)
	        					{
	        						if(biotype != null)
	        						{
		        						Integer c = json.biotypes.get(biotype);
		        						if(c == null) c = 0;
		        						json.biotypes.put(biotype, c + (int)v);
	        						}
	        						if(chr != null)
	        						{
		        						Integer c = json.chrs.get(chr);
		        						if(c == null) c = 0;
		        						json.chrs.put(chr, c + (int)v);
	        						}
	        					}
	        				}
	        				else json.is_count_table = false; // if this is not the case, then the table contains float numbers
	        				if(v == 0) json.nber_zeros++;
	        				if(isERCC) 
	        				{
	        					if(json.is_count_table) bw_ercc.write("\t" + (int)v);
	        					else bw_ercc.write("\t" + v);
	        				}
	        				else 
	        				{
	        					if(json.is_count_table)
	        					{
		        					bw_dl_tab.write("\t" + (int)v);
		        					bw_tab.write("\t" + (int)v);
	        					}
	        					else
	        					{
		        					bw_dl_tab.write("\t" + v);
		        					bw_tab.write("\t" + v);
	        					}
	        				}
	    				}
	    				catch(NumberFormatException nfe) // in case we don't have a number
	    				{
	    					error("Value at line " + current_line+ " col " + col + " is not a number.");
	    				}
        			}
        			if(isERCC) bw_ercc.write("\n");
        			else 
        			{
        				bw_dl_tab.write("\n");
        				bw_tab.write("\n");
        			}
        		}
                line = br.readLine();
            }
            bw_tab.close();
            bw_ercc.close();
            bw_NF.close();
            bw_dl_tab.close();
            bw_exonLength.close();
            br.close();
        } 
        catch(IOException ioe)
        {
        	System.err.println(ioe.getMessage());
        	try
        	{
            	BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
            	bw.write("{\"displayed_error\":\"" + ioe.getMessage() + "\"}");
            	bw.close();
        	}
        	catch(IOException ioe2)
        	{
        		System.err.println(ioe2.getMessage());
        	}
        	System.exit(-1);
        }
        writeDuplicatedGenes(geneDups);
        writeGeneNamesJSON(genes);
        json.nber_unique_genes = geneDups.size();
        for(String gene:geneDups.keySet()) if(geneDups.get(gene) != 1) json.nber_duplicated_genes++;
        if(json.nber_ercc == 0) new File(Parameters.outputFolder + "ercc.tab").delete(); // If no ERCC, I delete the ERCC file
        json.setBiotypes();
        json.writeOutputJSON();
    }
    
    public static void createHeader(ArrayList<String> cell_names, int num_columns, int num_columns_header) throws IOException // Create and store the header in static variable 'header'
    {
    	header = "Genes";  // This should always be written in the beginning
		if(Parameters.has_header)
		{
			if(num_columns == num_columns_header) // Header has no missing value
			{
				switch(Parameters.name_column)
				{
					case FIRST:
						for(int i = 1; i < cell_names.size(); i++) header += "\t" + cell_names.get(i); // Ignore first value
						json.nber_cells = num_columns - 1;
						break;
					case LAST:
						for(int i = 0; i < cell_names.size() - 1; i++) header += "\t" + cell_names.get(i); // Ignore last value
						json.nber_cells = num_columns - 1;
						break;
					case NONE:
						for(int i = 0; i < cell_names.size(); i++) header += "\t" + cell_names.get(i); // Take all
						json.nber_cells = num_columns;
				}
			}
			else // num_columns == num_columns_header - 1 => One missing column in the header
			{
				switch(Parameters.name_column)
				{
					case FIRST:
					case LAST:
						for(int i = 0; i < cell_names.size(); i++) header += "\t" + cell_names.get(i); // Take all
						json.nber_cells = num_columns - 1;
						break;
					case NONE:
						error("You stated that there are no Gene Names, but header is missing one column ( " + num_columns_header + " ) as compared to first line ( " + num_columns + " ).");
				}
			}
		}
		else  // Create my own header
		{
			switch(Parameters.name_column)
			{
				case FIRST:
				case LAST:
					json.nber_cells = num_columns - 1;
					break;
				case NONE:
					json.nber_cells = num_columns;
			}
			for(int i = 1; i <= json.nber_cells; i++) header += "\tCell_" + i;
		}
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
    
    public static void writeGeneNamesJSON(ArrayList<String> genes)
    {
    	try
    	{
        	BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "gene_names.json"));
        	if(genes.size() == 0) bw.write("{}");
        	else
        	{
	        	bw.write("[" + geneNameToJSON(genes.get(0)));
	        	for(int i = 1; i < genes.size(); i++) bw.write("," + geneNameToJSON(genes.get(i)));
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
    
    public static String geneNameToJSON(String name)
    {
    	int i1 = name.indexOf("|");
    	int i2 = name.indexOf("|", i1+1);
    	String originalGene = name.substring(0, i1);
    	String ensGene = name.substring(i1+1, i2);
    	String otherGenes = name.substring(i2+1, name.length());
    	return "[\"" + originalGene + "\",\"" + ensGene + "\",\"" + otherGenes + "\"]";
    }
    
    public static void error(String error)
    {
    	try
    	{
        	BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
        	bw.write("{\"displayed_error\":\"" + error + "\"}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		System.err.println(ioe.getMessage());
    	}
    	System.exit(-1);
    }
    

}
