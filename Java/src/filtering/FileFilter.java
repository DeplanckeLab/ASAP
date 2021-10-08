package filtering;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import bigarrays.DoubleArray64;
import bigarrays.StringArray64;
import hdf5.loom.LoomData;
import hdf5.loom.LoomFile;
import json.ErrorJSON;
import json.FilterCellsJSON;
import json.FilterJSON;
import model.MetaOn;
import model.Metadata;
import model.Metatype;
import model.Parameters;
import tools.Utils;

public class FileFilter 
{
	public static FilterJSON filtJSON = null;
	public static int[] expressedGenesPerSample = null;
	public static double[] rowVar = null; // var
	public static double[] rowCoeffOfVar = null; // coeffofvar
	public static double[] sizeFactors = null; // scLVM
	public static double[] loggeomeans = null;
	public static String yCol = "Nb Expressed Genes [count > 0]";
	
	public static void filterCells()
    {
		System.out.println("Filtering Cells of file : " + Parameters.loomFile);
		LoomFile loom = new LoomFile("r", Parameters.loomFile);
	   	// First: check empty cells
    	if(!loom.exists("/col_attrs/_Depth")) new ErrorJSON("Error in the Loom file. Path /col_attrs/_Depth should exist!");
    	DoubleArray64 depth = loom.readDoubleArray("/col_attrs/_Depth");
    	HashSet<Long> cellIndexesToFilter = new HashSet<Long>();
    	for(long i = 0; i < depth.size(); i++) if(depth.get(i) == 0) cellIndexesToFilter.add(i);
    	System.out.println(cellIndexesToFilter.size() +  " empty cells to filter.");
		System.out.print("Reading JSON file...");
		if(!Parameters.isIndex)
		{
			String[] toFilter = FilterCellsJSON.parseJSON(Parameters.JSONFileName).filtered_cells;
			HashSet<String> toFilterMap = new HashSet<>();
			for(String s:toFilter) toFilterMap.add(s);
			System.out.println(" OK! " + toFilterMap.size() + " cell names in the JSON to filter.");
			System.out.print("Checking if the cell names are indeed in the metadata... ");
			StringArray64 cellNames = loom.getCellNames();
			for (int i = 0; i < cellNames.size(); i++) {
				if(toFilterMap.contains(cellNames.get(i))) cellIndexesToFilter.add((long)i);
			}
		}
		else
		{
			long[] toFilter = FilterCellsJSON.parseJSON(Parameters.JSONFileName).discarded_cols;
			if(toFilter == null) new ErrorJSON("Field 'discarded_cols' is missing from JSON");
			for(long l:toFilter) cellIndexesToFilter.add(l);
		}

		System.out.println(cellIndexesToFilter.size() + " cells were found and will be filtered.");
		
		// Create new Loom filtered
		long[] dim = loom.getDimensions();
		System.out.println("Current /matrix size = " + dim[0] + " x " + dim[1]);
		if(dim[1] - cellIndexesToFilter.size() == 0) new ErrorJSON("No more cells after filtering...");
		System.out.println("Creating New Loom file with /matrix size = " + dim[0] + " x " + (dim[1] - cellIndexesToFilter.size()));
		LoomData data = new LoomData(dim[0], dim[1] - cellIndexesToFilter.size());
		data.is_count_table = true;
		LoomFile loomNew = new LoomFile("w", Parameters.outputFolder + "output.loom");
		
		// Process Main Matrix
		LoomFile.copyFloatMatrixByGene("/matrix", loom, loomNew, cellIndexesToFilter, data); // data should be updated

		// Process metadata
		System.out.print("Now processing the Metadata...");
		List<Metadata> meta = loom.listMetadata();
		System.out.println(meta.size() + " metadata found");
		for(Metadata m:meta) 
		{
			if(!m.path.equals("/row_attrs/_Sum")) 
			{
				LoomFile.copyMetadata(m, null, cellIndexesToFilter, loom, loomNew);
				data.meta.add(m);
			}
		}
		loomNew.writeDoubleArray("/row_attrs/_Sum", data.sum);
		Metadata m = new Metadata("/row_attrs/_Sum", Metatype.NUMERIC, MetaOn.GENE, 1, data.sum.size());
		m.size = loomNew.getSizeInBytes("/row_attrs/_Sum");
		data.meta.add(m);
		loom.close();
		loomNew.close();
		FilterCellsJSON.writeOutputJSON(data);
    }
	
	private static void filterGenes(LoomFile loom, HashSet<Long> geneIndexesToFilter)
    {
    	// Create new Loom filtered
		long[] dim = loom.getDimensions();
		System.out.println("Current /matrix size = " + dim[0] + " x " + dim[1]);
		if(dim[0] - geneIndexesToFilter.size() == 0) new ErrorJSON("No more genes after filtering...");
		System.out.println("Creating New Loom file with /matrix size = " + (dim[0] - geneIndexesToFilter.size()) + " x " + dim[1]);
		LoomData data = new LoomData(dim[0] - geneIndexesToFilter.size(), dim[1]);
		data.is_count_table = true;
		LoomFile loomNew = new LoomFile("w", Parameters.outputFolder + "output.loom");
		
		// Process Main Matrix
		LoomFile.copyFloatMatrixByCell("/matrix", loom, loomNew, geneIndexesToFilter, data); // data should be updated

		// Process metadata
		System.out.print("Now processing the Metadata...");
		List<Metadata> meta = loom.listMetadata();
		System.out.println(meta.size() + " metadata found");
		for(Metadata m:meta) 
		{
			if(!m.path.equals("/col_attrs/_Depth") && !m.path.equals("/col_attrs/_Detected_Genes"))
			{
				LoomFile.copyMetadata(m, geneIndexesToFilter, null, loom, loomNew);
				data.meta.add(m);
			}
		}
		loomNew.writeDoubleArray("/col_attrs/_Depth", data.depth);
		Metadata m = new Metadata("/col_attrs/_Depth", Metatype.NUMERIC, MetaOn.CELL, data.depth.size(), 1);
		m.size = loomNew.getSizeInBytes("/col_attrs/_Depth");
		data.meta.add(m);
		loomNew.writeIntArray("/col_attrs/_Detected_Genes", data.detected_genes);
		m = new Metadata("/col_attrs/_Detected_Genes", Metatype.NUMERIC, MetaOn.CELL, data.detected_genes.size(), 1);
		m.size = loomNew.getSizeInBytes("/col_attrs/_Detected_Genes");
		data.meta.add(m);
		loom.close();
		loomNew.close();
		FilterCellsJSON.writeOutputJSON(data);
    }
	
    public static void filterGenes()
    {
        try
        {
        	switch(Parameters.filtModel)
        	{
    			case BASIC:
    				filterBASIC();
    				break;
        		case CPM:
        			filterCPM();
        			break;
        		case KEEP:
        			filterKEEP();
        			break;
        		default:
        			System.err.println("Not implemented yet");
        	}
        } 
        catch(Exception ioe) // IOException ioe)
        {
        	new ErrorJSON(ioe.getMessage());
        }
    }
    
    public static void filterDEMetadata()
    {
    	LoomFile loom = new LoomFile("r", Parameters.loomFile);
    	// Check empty columns / empty rows
    	if(!Parameters.iAnnot.startsWith("/row_attrs")) new ErrorJSON("Your DE metadata should starts with /row_attrs (gene metadata)");
    	if(!loom.exists(Parameters.iAnnot)) new ErrorJSON("Error in the Loom file. Path "+Parameters.iAnnot+" should exist!");
    	float[][] deResults = loom.readFloatMatrix(Parameters.iAnnot);
    	if(deResults[0].length != 5) new ErrorJSON("Your DE metadata should contain 5 columns, it contains " + deResults[0].length);
    	loom.close();
    	
    	HashMap<Integer, Float> geneUpIndexesToKeep = new HashMap<Integer, Float>();
    	HashMap<Integer, Float> geneDownIndexesToKeep = new HashMap<Integer, Float>();
    	for(int i = 0; i < deResults.length; i++)
    	{
    		if(Math.pow(2, Math.abs(deResults[i][0])) >= Parameters.fcThreshold && deResults[i][1] <= Parameters.pThreshold && deResults[i][2] <= Parameters.fdrThreshold)
    		{
    			// Passing all thresholds
    			if(deResults[i][0] >= 0) geneUpIndexesToKeep.put(i, deResults[i][0]);
    			else geneDownIndexesToKeep.put(i, deResults[i][0]);
    		}
    	}
    	
    	// Keep only the top X values, if requested
    	int[] upIndexes = new int[Math.min(Parameters.topThreshold , geneUpIndexesToKeep.size())];
    	int[] sortedkeys = Utils.sortF(geneUpIndexesToKeep, true);
    	for(int i=0; i < upIndexes.length; i++) upIndexes[i] = sortedkeys[i];

    	// Keep only the top X values, if requested
    	int[] downIndexes = new int[Math.min(Parameters.topThreshold , geneDownIndexesToKeep.size())];
    	sortedkeys = Utils.sortF(geneDownIndexesToKeep, false);
    	for(int i=0; i < downIndexes.length; i++) downIndexes[i] = sortedkeys[i];
    	
    	// Create output file as a String
    	StringBuilder sb = new StringBuilder();
    	if(Parameters.displayValues)
    	{
			sb.append("{\"up\":[");
			String prefix = "";
			for(int i=0; i < upIndexes.length; i++) 
			{
				sb.append(prefix).append(upIndexes[i]);
				prefix = ",";
			}
			sb.append("],\"down\":[");
			prefix = "";
			for(int i=0; i < downIndexes.length; i++) 
			{
				sb.append(prefix).append(downIndexes[i]);
				prefix = ",";
			}
			sb.append("]}");
    	}
    	else
    	{
    		sb.append("{\"up\":").append(upIndexes.length).append(",\"down\":").append(downIndexes.length).append("}");
    	}
		
		// Writing results
		Utils.writeJSON(sb, Parameters.JSONFileName);
    }
    
    public static void filterKEEP() // Only remove empty genes
    {
    	LoomFile loom = new LoomFile("r", Parameters.loomFile);
    	// Check empty columns / empty rows
    	if(!loom.exists("/row_attrs/_Sum")) new ErrorJSON("Error in the Loom file. Path /row_attrs/_Sum should exist!");
    	DoubleArray64 sum = loom.readDoubleArray("/row_attrs/_Sum");
    	HashSet<Long> geneIndexesToFilter = new HashSet<Long>();
    	for(long i = 0; i < sum.size(); i++) if(sum.get(i) == 0) geneIndexesToFilter.add(i);
    	System.out.println(geneIndexesToFilter.size() +  " empty genes to filter.");
    	
    	int[] toKeep = FilterCellsJSON.parseJSON(Parameters.JSONFileName).kept_genes;
    	if(toKeep == null) new ErrorJSON("JSON file does not have 'kept_genes' field");
    	HashSet<Integer> toKeepMap = new HashSet<>();
    	for(int s:toKeep) toKeepMap.add(s);
    	
    	for(long i = 0; i < sum.size(); i++) if(sum.get(i) == 0 || !toKeepMap.contains((int)i)) geneIndexesToFilter.add(i);
    	System.out.println(geneIndexesToFilter.size() +  " genes to filter.");
    	
    	filterGenes(loom, geneIndexesToFilter);
    }
    
    public static void filterBASIC() // Only remove empty genes
    {
    	LoomFile loom = new LoomFile("r", Parameters.loomFile);
    	// Check empty columns / empty rows
    	if(!loom.exists("/row_attrs/_Sum")) new ErrorJSON("Error in the Loom file. Path /row_attrs/_Sum should exist!");
    	DoubleArray64 sum = loom.readDoubleArray("/row_attrs/_Sum");
    	HashSet<Long> geneIndexesToFilter = new HashSet<Long>();
    	for(long i = 0; i < sum.size(); i++) if(sum.get(i) == 0) geneIndexesToFilter.add(i);
    	System.out.println(geneIndexesToFilter.size() +  " empty genes to filter.");
    	filterGenes(loom, geneIndexesToFilter);
    }
    
    public static void filterCPM()
    {
    	LoomFile loom = new LoomFile("r", Parameters.loomFile);
    	// Check empty columns / empty rows
    	if(!loom.exists("/row_attrs/_Sum")) new ErrorJSON("Error in the Loom file. Path /row_attrs/_Sum should exist!");
    	DoubleArray64 sum = loom.readDoubleArray("/row_attrs/_Sum");
    	HashSet<Long> geneIndexesToFilter = new HashSet<Long>();
    	for(long i = 0; i < sum.size(); i++) if(sum.get(i) == 0) geneIndexesToFilter.add(i);
    	System.out.println(geneIndexesToFilter.size() +  " empty genes to filter.");
     	
    	// Read the file block per block according to natural storage
		long[] dim = loom.getDimensions();
		//float[][] cpm = new float[(int)dim[0]][(int)dim[1]];
		int[] blockSize = loom.getChunkSizes();
		int nbTotalBlocks = (int)Math.ceil((double)dim[0] / blockSize[0]);
		System.out.print("Reading " + nbTotalBlocks + " independent blocks...");
		DoubleArray64 depth = loom.readDoubleArray("/col_attrs/_Depth");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{		
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = loom.readFloatBlock("/matrix", blockSize[0], (int)dim[1], nbBlocks, 0l);
			
			// Parsing Data and generating summary annotations
			
			for(int x = 0; x < subMatrix.length; x++)
			{
				long i = x + nbBlocks * blockSize[0]; // Original gene index
				if(!geneIndexesToFilter.contains(i))
				{
					int detected = 0;
					for(int j = 0; j < subMatrix[0].length; j++)
					{
						float value = subMatrix[x][j];
							
						//cpm[(int)i][(int)j] = (value / (float)depth.get(j)) * 1000000;
						//if(cpm[(int)i][(int)j] > Parameters.nbCountsPerCell) detected++;
						
						float cpm = (value / (float)depth.get(j)) * 1000000;
						if(cpm > Parameters.nbCountsPerCell) detected++;
					}
					if(detected < Parameters.nbCellsDetected) geneIndexesToFilter.add(i);
				}
			}
		}
		//loom.writeLayer("_CPM", cpm);
		filterGenes(loom, geneIndexesToFilter);
    }

}
