package parsing;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import bigarrays.StringArray64;
import compression.CompressionHandler;
import hdf5.HDF5Tools;
import hdf5.h510x.H510xHandler;
import hdf5.loom.LoomFile;
import hdf5.loom.LoomHandler;
import json.ErrorJSON;
import json.ParsingJSON;
import json.PreparsingJSON;
import model.MapGene;
import model.MetaOn;
import model.Metadata;
import model.Parameters;
import parsing.model.ColumnName;
import parsing.model.FileType;
import parsing.model.Gene;
import parsing.model.GroupPreparse;
import parsing.model.MetaPreparse;
import tools.Utils;

public class FileParser
{
    public static void parse()
    {
    	System.out.println("Parsing file : " + Parameters.fileName + " as a " + Parameters.fileType + " file.");
    	switch(Parameters.fileType) // We assume here that the file type is correct (because we preparsed the file)
    	{
    		case RAW_TEXT:
				try
				{
					parseText(new BufferedReader(new FileReader(Parameters.fileName)));
				}
				catch(NumberFormatException nfe) { new ErrorJSON("Your data matrix contains non-numeric values." + (Parameters.has_header?"":"No header is specified in the options, but maybe there is one?") + (Parameters.name_column == ColumnName.NONE?"":"You did not specify a column with gene_name, but maybe there is one?") + (Parameters.name_column == ColumnName.FIRST?"":"You specified first column as the one containing gene_name, but maybe it's another one?") + (Parameters.name_column == ColumnName.LAST?"":"You specified last column as the one containing gene_name, but maybe it's another one?")); }
				catch(IOException ioe) { new ErrorJSON("Tried to read your file as a PLAIN TEXT file, but it failed."); }
    			break;
    		case ARCHIVE:
    		case ARCHIVE_COMPRESSED:
    		case COMPRESSED:
				try
				{
					parseText(new BufferedReader(new InputStreamReader(CompressionHandler.getReader())));
				}
				catch(NumberFormatException nfe) { new ErrorJSON("Your data matrix contains non-numeric values. Did you forget there is a header?\",\"detected_format\":\"" + Parameters.fileType); }
				catch(IOException ioe) { new ErrorJSON("Tried to read your file as a ARCHIVE/COMPRESSED file, but it failed."); }
				break;
    		case H5_10x:
    			H510xHandler.parse();
    			break;
    		case LOOM:
    			LoomHandler.parse();
    			break;
    		default:
    			new ErrorJSON("This file format is not handled.");
    	}
    }
    
    public static void preparse()
    {
    	System.out.println("Preparsing file : " + Parameters.fileName);
    	HDF5Tools.checkIfHDF5orLoom(); // Check if it is an HDF5 archive
    	if(Parameters.fileType == FileType.H5_10x) H510xHandler.preparse();
    	else if(Parameters.fileType == FileType.LOOM) LoomHandler.preparse();
    	else // It's a text file, or an archive
    	{
	    	long nRows = -1;
			if(Parameters.selection == null)
			{
				ArrayList<String> files = CompressionHandler.getList();
				if(files == null) 
				{
					System.out.println("File format not detected. Assumed to be plain text.");
					Parameters.fileType = FileType.RAW_TEXT;
					try
					{
						nRows = Utils.countNonEmptyLines(new BufferedInputStream(new FileInputStream(Parameters.fileName)));
						preparseText(nRows, new BufferedReader(new FileReader(Parameters.fileName)));
					}
					catch(NumberFormatException nfe) { new ErrorJSON("Your data matrix contains non-numeric values. Did you forget there is a header?\",\"detected_format\":\"" + Parameters.fileType); }
					catch(IOException ioe) { new ErrorJSON("Tried to read your file as a PLAIN TEXT file, but it failed.\",\"detected_format\":\"" + Parameters.fileType); }
				}
				else
				{
					if(files.size() == 0) new ErrorJSON("Your archive/compressed file is empty or corrupted.\",\"detected_format\":\"" + Parameters.fileType);
					else if(files.size() == 1)
					{
						System.out.println("A single file was detected in your archive. Processing automatically the preparsing...");
						Parameters.selection = files.get(0);
	    				try
	    				{
	    					nRows = Utils.countNonEmptyLines(new BufferedInputStream(CompressionHandler.getReader()));
	    					preparseText(nRows, new BufferedReader(new InputStreamReader(CompressionHandler.getReader())));
	    				}
	    				catch(NumberFormatException nfe) { new ErrorJSON("Your data matrix contains non-numeric values. Did you forget there is a header?\",\"detected_format\":\"" + Parameters.fileType); }
	    				catch(IOException ioe) { new ErrorJSON("Tried to read your file as a COMPRESSED/ARCHIVED file, but it failed.\",\"detected_format\":\"" + Parameters.fileType); }
					}
					else
					{
						System.out.println("Several files are detected in your archive, but no selection was made using -sel option. Therefore creating a JSON listing the archive content.");
						PreparsingJSON.writeListingJSON(files);
					}
				}
			}
			else
			{
				try
				{
					nRows = Utils.countNonEmptyLines(new BufferedInputStream(CompressionHandler.getReader()));
					preparseText(nRows, new BufferedReader(new InputStreamReader(CompressionHandler.getReader())));
				}
				catch(NumberFormatException nfe) { new ErrorJSON("Your data matrix contains non-numeric values. Did you forget there is a header?\",\"detected_format\":\"" + Parameters.fileType); }
				catch(IOException ioe) { new ErrorJSON("Tried to read your file as a COMPRESSED/ARCHIVED file, but it failed.\",\"detected_format\":\"" + Parameters.fileType); }
			}
    	}
    }
    
    public static void parseMetadata()
    {
    	
    }
    
    public static void preparseMetadata()
    {
    	System.out.println("Preparsing metadata file : " + Parameters.fileName);
    	HDF5Tools.checkIfHDF5orLoom(); // Check if it is an HDF5 archive
    	if(Parameters.fileType == FileType.H5_10x || Parameters.fileType == FileType.LOOM) new ErrorJSON("This file is a HDF5 file. ASAP cannot extract metadata from these files atm.");
    	
		if(Parameters.selection == null)
		{
			ArrayList<String> files = CompressionHandler.getList();
			if(files == null) 
			{
				System.out.println("File format not detected. Assumed to be plain text.");
				Parameters.fileType = FileType.RAW_TEXT;
				try
				{
					preparseMeta(new BufferedReader(new FileReader(Parameters.fileName)));
				}
				catch(NumberFormatException nfe) { new ErrorJSON("Your metadata matrix contains non-numeric values. Did you forget there is a header?\",\"detected_format\":\"" + Parameters.fileType); }
				catch(IOException ioe) { new ErrorJSON("Tried to read your file as a PLAIN TEXT file, but it failed.\",\"detected_format\":\"" + Parameters.fileType); }
			}
			else
			{
				if(files.size() == 0) new ErrorJSON("Your archive/compressed file is empty or corrupted.\",\"detected_format\":\"" + Parameters.fileType);
				else if(files.size() == 1)
				{
					System.out.println("A single file was detected in your archive. Processing automatically the preparsing...");
					Parameters.selection = files.get(0);
    				try
    				{
    					preparseMeta(new BufferedReader(new InputStreamReader(CompressionHandler.getReader())));
    				}
    				catch(NumberFormatException nfe) { new ErrorJSON("Your metadata matrix contains non-numeric values. Did you forget there is a header?\",\"detected_format\":\"" + Parameters.fileType); }
    				catch(IOException ioe) { new ErrorJSON("Tried to read your file as a COMPRESSED/ARCHIVED file, but it failed.\",\"detected_format\":\"" + Parameters.fileType); }
				}
				else
				{
					System.out.println("Several files are detected in your archive, but no selection was made using -sel option. Therefore creating a JSON listing the archive content.");
					PreparsingJSON.writeListingJSON(files);
				}
			}
		}
		else
		{
			try
			{
				preparseMeta(new BufferedReader(new InputStreamReader(CompressionHandler.getReader())));
			}
			catch(NumberFormatException nfe) { new ErrorJSON("Your metadata matrix contains non-numeric values. Did you forget there is a header?\",\"detected_format\":\"" + Parameters.fileType); }
			catch(IOException ioe) { new ErrorJSON("Tried to read your file as a COMPRESSED/ARCHIVED file, but it failed.\",\"detected_format\":\"" + Parameters.fileType); }
		}
    	
    }
    
    public static void preparseMeta(BufferedReader br) throws IOException
    {
    	// JSON object
    	MetaPreparse g;
    	if(Parameters.selection == null) g = new MetaPreparse(Parameters.fileName.substring(Parameters.fileName.lastIndexOf("/") + 1));
    	else g = new MetaPreparse(Parameters.selection);
    	
    	// Initialize objects that will be output in JSON
    	LoomFile loom = new LoomFile("r", Parameters.loomFile);
    	g.not_found = new HashSet<>();
    	g.ambiguous = new HashSet<>();
    	if(Parameters.which == MetaOn.CELL)
    	{
    		g.cellNames = loom.getCellNames();
    		g.cellMap = new HashMap<>();
    		for(int i = 0; i < g.cellNames.size(); i++) g.cellMap.put(g.cellNames.get(i), i); // For faster comparison
    	}
    	else
    	{
    		g.geneNames = loom.getGeneNames();
    		g.geneMap = new HashMap<>();
    		for(int i = 0; i < g.geneNames.length; i++) 
    		{
    			Gene gene = g.geneNames[i];
				String[] tokens;
				if(gene.name != null)
				{	
					tokens = gene.name.split(",");
					for(String t:tokens) 
					{
						List<Integer> list = g.geneMap.get(t);
						if(list == null) list = new ArrayList<Integer>();
						list.add(i);
						g.geneMap.put(t, list); // For faster comparison
					}
				}
				if(gene.ensembl_id != null)
				{
					tokens = gene.ensembl_id.split(",");
					for(String t:tokens) 
					{
						List<Integer> list = g.geneMap.get(t);
						if(list == null) list = new ArrayList<Integer>();
						list.add(i);
						g.geneMap.put(t, list); // For faster comparison
					}
				}
				// We do not store alt names in the Loom file anyway
    		}
    	}
    	g.metaIndex = new HashMap<>();
    	loom.close();
    	
		// How to read each line?
		int startI = 0;
		int endI = 0;
		int nameI = 0;
    	
    	// Start to read file
    	String line = br.readLine();
    	String[] tokens = line.split(Parameters.delimiter);

    	if(Parameters.has_header) 
    	{
    		String[] header = tokens;
    		line = br.readLine();
    		tokens = line.split(Parameters.delimiter);
    		if(header.length != tokens.length && header.length - 1 != tokens.length) new ErrorJSON("Header size (" + header.length +  ") does not match first line size (" + tokens.length + ").");
    		
    		// First, need to read header for metadata names. Thus, search for header names
    		endI = header.length;
    		switch(Parameters.name_column)
    		{
    			case FIRST:
    				if(header.length == tokens.length) startI = 1;
    				break;
    			case LAST:
    				if(header.length == tokens.length) endI = header.length - 1;
    				break;
    			case NONE:
    				if(header.length != tokens.length) new ErrorJSON("Header size (" + header.length +  ") does not match first line size (" + tokens.length + ").");
    		}
    		
    		// Create a Metadata object for each header, except for the Gene/Cell Name
    		for (int i = startI; i < endI; i++) 
    		{
    			Metadata m = new Metadata();
    			m.path = header[i].trim();
    			m.on = Parameters.which;
    			if(Parameters.which == MetaOn.CELL) m.values = new StringArray64(g.cellNames.size());
    			else m.values = new StringArray64(g.geneNames.length);
    			g.metaIndex.put(i, m);
			}
    		
    		// Now, we will be reading the Metadata values  		
    		switch(Parameters.name_column)
    		{
    			case FIRST:
    				startI = 1;
    				nameI = 0;
    				endI = header.length;
    				break;
    			case LAST:
      				startI = 0;
    				nameI = header.length - 1;
    				endI = header.length - 1;
    				break;
    			case NONE:
      				startI = 0;
    				nameI = -1;
    				endI = header.length;
    		}
    	}
    	else
    	{  		
    		endI = tokens.length;
    		switch(Parameters.name_column)
    		{
    			case FIRST:
    				startI = 1;
    				break;
    			case LAST:
    				nameI = tokens.length - 1;
    				endI = tokens.length - 1;
    				break;
    			case NONE:
    				nameI = -1;
    		}
    		
    		// Create a Metadata object for each value, except for the Gene/Cell Name
    		for (int i = startI; i < endI; i++) 
    		{
    			Metadata m = new Metadata();
    			m.path = "Metadata_"+(i+1); // Generated name
    			m.on = Parameters.which;
    			if(Parameters.which == MetaOn.CELL) m.values = new StringArray64(g.cellNames.size());
    			else m.values = new StringArray64(g.geneNames.length);
    			g.metaIndex.put(i, m);
			}
    	}
   
    	// Read remaining of the file
    	long l = 0;
  		while(line != null)
		{
    		for (int i = startI; i < endI; i++) 
    		{
    			Metadata m = g.metaIndex.get(i);
    			if(nameI == -1) m.values.set(l, tokens[i]); // Assumed to be ordered the same way	
    			else
    			{
    				String name = tokens[nameI];
    				if(Parameters.which == MetaOn.CELL)
    				{
		    			Integer index = g.cellMap.get(name);
		    			if(index == null) g.not_found.add(name);
		    			else m.values.set(index, tokens[i]);
    				}
    				else
    				{
		    			List<Integer> indexes = g.geneMap.get(name);
		    			if(indexes == null) g.not_found.add(name);
		    			else if(indexes.size() == 1) m.values.set(indexes.get(0), tokens[i]);
		    			else g.ambiguous.add(name);
    				}
    			}
			}
    		l++;
    		line = br.readLine();
    		if(line != null) tokens = line.split(Parameters.delimiter);
		}
		
		// Inferring metadata type
		for(Integer index:g.metaIndex.keySet()) g.metaIndex.get(index).inferType(); 
    		
    	PreparsingJSON.writeOutputJSON(g);
    }
    
    private static void preparseText(long nRows, BufferedReader br) throws IOException, NumberFormatException
    {
    	String line = br.readLine();
    	ArrayList<GroupPreparse> res = new ArrayList<>(); // Will just contain one element here
    	GroupPreparse g;
    	if(Parameters.selection == null) g = new GroupPreparse(Parameters.fileName.substring(Parameters.fileName.lastIndexOf("/") + 1));
    	else g = new GroupPreparse(Parameters.selection);

    	res.add(g);
    	String[] tokens = line.split(Parameters.delimiter);
    	if(Parameters.has_header) 
    	{
    		String[] header = tokens;
    		line = br.readLine();
    		if(line == null) new ErrorJSON("There is no data in your file. Possible reasons: Not a matrix or weird encoding (is it UTF-8?)\",\"detected_format\":\"" + Parameters.fileType);
    		tokens = line.split(Parameters.delimiter);
    		g.nbGenes = nRows - 1;
    		switch(Parameters.name_column)
    		{
    			case FIRST:
    				int startI = 1;
    				if(header.length == tokens.length) g.nbCells = header.length - 1;
    				else if(header.length == tokens.length - 1) { g.nbCells = header.length; startI = 0;}
    				else new ErrorJSON("Header has " + header.length + " elements, while l.2 has " + tokens.length + " elements. Please correct file or parsing parameters.\",\"detected_format\":\"" + Parameters.fileType);
    				g.cellNames = new String[(int)Math.min(10, g.nbCells)]; // TODO if too many cells (more than 2,147,483,647), it will not work. Use the String64 array instead
    				for(int i = startI;i < g.cellNames.length + startI; i++) g.cellNames[i-startI] = header[i];
    				g.geneNames = new String[(int)Math.min(10, g.nbGenes)]; // Should not be > 2billions
    				g.matrix = new float[g.geneNames.length][g.cellNames.length];
    				for(int i=0; i < g.matrix.length; i++)
    				{
    					if(tokens.length != g.nbCells + 1) new ErrorJSON("There should be " + g.nbCells + " cells, but l."+ (i+2) + " has " + tokens.length + " elements. Please correct file or parsing parameters.\",\"detected_format\":\"" + Parameters.fileType);
        				g.geneNames[i] = tokens[0];
        				for(int j = 1; j < g.matrix[i].length + 1; j++)
        				{
        					g.matrix[i][j - 1] = Float.parseFloat(tokens[j]);
        					if(g.matrix[i][j - 1] != (int)g.matrix[i][j - 1]) g.isCount = false;
        				}
        				tokens = br.readLine().split(Parameters.delimiter);
    				}
    				break;
    			case LAST:
    				if(header.length == tokens.length) g.nbCells = header.length - 1;
    				else if(header.length == tokens.length - 1) g.nbCells = header.length;
    				else new ErrorJSON("Header has " + header.length + " elements, while l.2 has " + tokens.length + " elements. Please correct file or parsing parameters.\",\"detected_format\":\"" + Parameters.fileType);
    				g.cellNames = new String[(int)Math.min(10, g.nbCells)]; // TODO if too many cells (more than 2,147,483,647), it will not work. Use the String64 array instead
    				for(int i = 0;i < g.cellNames.length; i++) g.cellNames[i] = header[i];
    				g.geneNames = new String[(int)Math.min(10, g.nbGenes)]; // Should not be > 2billions
    				g.matrix = new float[g.geneNames.length][g.cellNames.length];
    				for(int i=0; i < g.matrix.length; i++)
    				{
    					if(tokens.length != g.nbCells + 1) new ErrorJSON("There should be " + g.nbCells + " cells, but l."+ (i+2) + " has " + tokens.length + " elements. Please correct file or parsing parameters.\",\"detected_format\":\"" + Parameters.fileType);
        				g.geneNames[i] = tokens[(int)g.nbCells];  // TODO if too many cells (more than 2,147,483,647), it will not work. Use the String64 array instead
        				for(int j = 0; j < g.matrix[i].length; j++)
        				{
        					g.matrix[i][j] = Float.parseFloat(tokens[j]);
        					if(g.matrix[i][j] != (int)g.matrix[i][j]) g.isCount = false;
        				}
        				tokens = br.readLine().split(Parameters.delimiter);
    				}
    				break;
    			case NONE:
    				if(header.length == tokens.length) g.nbCells = header.length;
    				else new ErrorJSON("Header has " + header.length + " elements, while l.2 has " + tokens.length + " elements. Please correct file or parsing parameters.\",\"detected_format\":\"" + Parameters.fileType);
    				g.cellNames = new String[(int)Math.min(10, g.nbCells)]; // TODO if too many cells (more than 2,147,483,647), it will not work. Use the String64 array instead
    				for(int i = 0;i < g.cellNames.length; i++) g.cellNames[i] = header[i];
    				g.matrix = new float[(int)Math.min(10, g.nbGenes)][g.cellNames.length];
    				for(int i = 0; i < g.matrix.length; i++)
    				{
    					if(tokens.length != g.nbCells) new ErrorJSON("There should be " + g.nbCells + " cells, but l."+ (i+2) + " has " + tokens.length + " elements. Please correct file or parsing parameters.\",\"detected_format\":\"" + Parameters.fileType);
        				for(int j = 0; j < g.matrix[i].length; j++)
        				{
        					g.matrix[i][j] = Float.parseFloat(tokens[j]);
        					if(g.matrix[i][j] != (int)g.matrix[i][j]) g.isCount = false;
        				}
        				tokens = br.readLine().split(Parameters.delimiter);
    				}
    		}
    	}
    	else
    	{
    		g.nbGenes = nRows;
    		switch(Parameters.name_column)
    		{
    			case FIRST:
      				g.nbCells = tokens.length - 1;
    				g.geneNames = new String[(int)Math.min(10, g.nbGenes)]; // Should not be > 2billions
    				g.matrix = new float[g.geneNames.length][(int)Math.min(10, g.nbCells)];
    				for(int i=0; i < g.matrix.length; i++)
    				{
    					if(tokens.length != g.nbCells + 1) new ErrorJSON("There should be " + g.nbCells + " cells, but l."+ (i+1) + " has " + tokens.length + " elements. Please correct file or parsing parameters.\",\"detected_format\":\"" + Parameters.fileType);
        				g.geneNames[i] = tokens[0];
        				for(int j = 1; j < g.matrix[i].length + 1; j++)
        				{
         					g.matrix[i][j - 1] = Float.parseFloat(tokens[j]);
        					if(g.matrix[i][j - 1] != (int)g.matrix[i][j - 1]) g.isCount = false;
        				}
        				tokens = br.readLine().split(Parameters.delimiter);
    				}
    				break;
    			case LAST:
    				g.nbCells = tokens.length - 1;
    				g.geneNames = new String[(int)Math.min(10, g.nbGenes)]; // Should not be > 2billions
    				g.matrix = new float[g.geneNames.length][(int)Math.min(10, g.nbCells)];
    				for(int i=0; i < g.matrix.length; i++)
    				{
    					if(tokens.length != g.nbCells + 1) new ErrorJSON("There should be " + g.nbCells + " cells, but l."+ (i+1) + " has " + tokens.length + " elements. Please correct file or parsing parameters.\",\"detected_format\":\"" + Parameters.fileType);
        				g.geneNames[i] = tokens[(int)g.nbCells];  // TODO if too many cells (more than 2,147,483,647), it will not work. Use the String64 array instead
        				for(int j = 0; j < g.matrix[i].length; j++)
        				{
        					g.matrix[i][j] = Float.parseFloat(tokens[j]);
        					if(g.matrix[i][j] != (int)g.matrix[i][j]) g.isCount = false;
        				}
        				tokens = br.readLine().split(Parameters.delimiter);
    				}
    				break;
    			case NONE:
    				g.nbCells = tokens.length;
    				g.matrix = new float[(int)Math.min(10, g.nbGenes)][(int)Math.min(10, g.nbCells)];
    				for(int i = 0; i < g.matrix.length; i++)
    				{
    					if(tokens.length != g.nbCells) new ErrorJSON("There should be " + g.nbCells + " cells, but l."+ (i+1) + " has " + tokens.length + " elements. Please correct file or parsing parameters.\",\"detected_format\":\"" + Parameters.fileType);
        				for(int j = 0; j < g.matrix[i].length; j++)
        				{
        					g.matrix[i][j] = Float.parseFloat(tokens[j]);
        					if(g.matrix[i][j] != (int)g.matrix[i][j]) g.isCount = false;
        				}
        				tokens = br.readLine().split(Parameters.delimiter);
    				}
    		}
    	}
    	if(g.geneNames != null) for (int i = 0; i < g.geneNames.length; i++) g.geneNames[i] = g.geneNames[i].replaceAll("\"", "");
    	if(g.cellNames != null) for (int i = 0; i < g.cellNames.length; i++) g.cellNames[i] = g.cellNames[i].replaceAll("\"", "");
    	PreparsingJSON.writeOutputJSON(res);
    }
    
    private static void parseText(BufferedReader br) throws IOException
    {	
    	ParsingJSON json = new ParsingJSON();
    	String line = br.readLine();
    	
	    // Prepare the writing in the file
    	json.loom.createEmptyFloat32MatrixDataset("/matrix", Parameters.nGenes, Parameters.nCells, Parameters.defaultChunkX, Parameters.defaultChunkY); // "/matrix", 
    	json.data.ens_names = new StringArray64(json.data.nber_genes); // After checking the database, this is filled appropriately
    	json.data.gene_names = new StringArray64(json.data.nber_genes); // After checking the database, this is filled appropriately
    	
    	// Header
    	if(Parameters.has_header) // Parsing Header if exist
    	{
    		String[] header = line.split(Parameters.delimiter); // TODO does not work if too many cells
    		if(header.length == Parameters.nCells) json.data.cell_names = new StringArray64(header);
    		else if(header.length == Parameters.nCells + 1)
    		{
	    		switch(Parameters.name_column)
	    		{
	    			case FIRST: json.data.cell_names = new StringArray64(Arrays.copyOfRange(header, 1, header.length)); break; // TODO does not work if too many cells
	    			case LAST: json.data.cell_names = new StringArray64(Arrays.copyOfRange(header, 0, header.length - 1)); break; // TODO does not work if too many cells
	    			case NONE: new ErrorJSON("Header should contain " + Parameters.nCells + " elements, it has " + header.length); break;
	    		}
    		} 
    		else new ErrorJSON("Header should contain " + Parameters.nCells + " elements, it has " + header.length);
    		for (int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, json.data.cell_names.get(i).trim().replaceAll("'|\"", ""));
    		line = br.readLine();
    	}
    	else for(int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, "Cell_"+(i+1)); // Create fake names if none exist
    	
    	// How many columns to expect per line
    	long ncols = Parameters.nCells + 1;
    	long start = 0;
    	long end = ncols;
    	long current_line = 0;
		switch(Parameters.name_column)
		{
			case FIRST: start++; break;
			case LAST: end--; break;
			case NONE: end--; ncols--; break;
		}

		// Prepare the blocks that will contain all columns (because we write nbChunk Genes x All Cells)
		int nbBlocks = 0;
		int indexRow = 0;
		float[][] subMatrix = new float[Parameters.defaultChunkY][(int)Parameters.nCells];
        while(line != null)
        {
        	current_line++; // To shift the count from [0, n-1] TO [1, n]
        	
        	String[] tokens = line.split(Parameters.delimiter);    
        	if(tokens.length != ncols) new ErrorJSON("Row " + current_line+ " contains a different number of values (" + tokens.length + ") than the other rows (" + ncols + ")");
        	
 			// If I am here, it means that the number of values is correct => I grab the gene name if exists
			String gene = null;
			switch(Parameters.name_column)
			{
				case FIRST: gene = tokens[0].trim().replaceAll("'|\"", ""); break;
				case LAST: gene = tokens[tokens.length - 1].trim().replaceAll("'|\"", ""); break;
				case NONE: gene = "Gene_" + current_line; // Starts at 1 is more pretty
			}
        			
			// Reading the gene name, and checking in DB
			MapGene.addGene(gene, json.data, current_line - 1);
			MapGene.parseGene(json.data, current_line - 1);
	
			// Handle biotypes/Mito
			boolean isProteinCoding = false;
			boolean isMito = false;
			boolean isRibo = false;
			String biotype = json.data.biotypes.get(current_line - 1);
			String chr = json.data.chromosomes.get(current_line - 1);
			if(biotype.equals("protein_coding")) isProteinCoding = true;
			else if(biotype.equals("rRNA")) isRibo = true;
			if(chr.equals("MT")) isMito = true;
			
			// Reading the rest of the values
			for(long j = start; j < end; j++)
			{
				// Get the next value
				String stringValue = tokens[(int)j].trim().replaceAll(",", ".").replaceAll("'|\"", ""); // If French representation of float numbers // TODO fails if too big array
				
				// Check it is not NA or empty
    			if(stringValue.equals("")) new ErrorJSON("A value at line " + current_line+ " is empty, but this is not allowed in the current version of ASAP.");
    			if(stringValue.equals("NA")) new ErrorJSON("A value at line " + current_line + " is NA, but this is not allowed in the current version of ASAP.");
    			
    			// Parse to a float and fill the appropriate arrays
				try
				{
    				float v = Float.parseFloat(stringValue);
    				int rv = Math.round(v);
    				if(Math.abs(rv - v) < 1E-5) v = rv;// if the diff between the value and the rounded is lower than 1E-5 I round it (to avoid Java bugs with float representation)
    				else json.data.is_count_table = false;

    				subMatrix[indexRow][(int)j - (int)start] = v;
    				
    				// Handle biotype/chromosome count per cell
    				if(isProteinCoding) json.data.proteinCodingContent.set(j - start, json.data.proteinCodingContent.get(j - start) + v); 
    				if(isRibo) json.data.ribosomalContent.set(j - start, json.data.ribosomalContent.get(j - start) + v); 
    				if(isMito) json.data.mitochondrialContent.set(j - start, json.data.mitochondrialContent.get(j - start) + v);
    				
    				// Generate the sums
    				json.data.depth.set(j - start, json.data.depth.get(j - start) + v);
    				json.data.sum.set(current_line - 1, json.data.sum.get(current_line - 1) + v);
    				
    				// Number of zeroes / detected genes      				
    				if(v == 0) json.data.nber_zeros++;
    				else json.data.detected_genes.set(j - start, json.data.detected_genes.get(j - start) + 1);
				}
				catch(NumberFormatException nfe) // in case we don't have a number
				{
					new ErrorJSON("Value '"+stringValue+"' at line " + current_line+ " col " + (j+1) + " is not a number.");
				}
			}
			
			// One gene is finished
			indexRow++;
			line = br.readLine();
			
			// Write content of buffer if full
			if(indexRow >= subMatrix.length || line == null) // If buffer is full or if no more line to read
			{
				json.loom.writeFloatBlockDataset("/matrix", subMatrix, nbBlocks, 0);
				subMatrix = new float[Parameters.defaultChunkY][(int)Parameters.nCells];
				indexRow = 0;
				nbBlocks++;
			}
        }
        json.loom.resizeDataset("/matrix", Parameters.nGenes, Parameters.nCells); // Cause writing fixed-size blocks can extend the matrix size with 0
        LoomFile.fillLoomFile(json.loom, json.data);
        json.writeOutputJSON();
        json.loom.close();
    }
}
