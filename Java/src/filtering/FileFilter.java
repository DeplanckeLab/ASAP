package filtering;

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
			for(long l:toFilter) cellIndexesToFilter.add(l);
		}

		System.out.println(cellIndexesToFilter.size() + " cells were found and will be filtered.");
		
		// Create new Loom filtered
		System.out.println("Filtering the LOOM file...");
		long[] dim = loom.getDimensions();
		System.out.println("Current /matrix size = " + dim[0] + " x " + dim[1]);
		System.out.println("Creating New Loom file with /matrix size = " + dim[0] + " x " + (dim[1] - cellIndexesToFilter.size()));
		
		LoomData data = new LoomData(dim[0], dim[1] - cellIndexesToFilter.size());
		LoomFile loomNew = new LoomFile("w", Parameters.outputFolder + "output.loom");
		
		// Process Main Matrix
		loomNew.createEmptyMatrix(data.nber_genes, data.nber_cells);
    	    	
    	// Read the file block per block according to natural storage
		int[] blockSize = loom.getChunkSizes();
		int nbTotalBlocks = (int)Math.ceil((double)data.nber_genes / blockSize[0]);
		System.out.print("Writing " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{		
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = loom.readBlock(blockSize[0], (int)dim[1], nbBlocks);
			
			// Parsing Data and generating summary annotations
			for(int x = 0; x < subMatrix.length; x++)
			{
				for(int j = 0; j < subMatrix[0].length; j++)
				{
					if(!cellIndexesToFilter.contains((long)j)) // Because I don't count the filtered cells
					{
						int i = x + nbBlocks * blockSize[0]; // Original index
						
						float value = subMatrix[x][j];
						
						if(data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) data.is_count_table = false;
	        				
	        			// Re-Generate the sums
						data.sum.set(i, data.sum.get(i) + value);
						
						// Number of zeroes / detected genes 
						if(value == 0) data.nber_zeros++;
					}
				}
			}
			
			// Writing this merged block to output
			loomNew.writeBlockInMatrix(subMatrix, data.nber_cells, nbBlocks * blockSize[0], cellIndexesToFilter);
		}
		System.out.println(" OK!");
		
		// Process metadata
		System.out.print("Now processing the Metadata...");
		List<Metadata> meta = loom.listMetadata();
		System.out.println(meta.size() + " metadata found");
		for(Metadata m:meta) 
		{
			if(!m.path.equals("/row_attrs/_Sum")) 
			{
				LoomFile.copyMetadata(m, cellIndexesToFilter, true, loom, loomNew);
				data.meta.add(m);
			}
		}
		loomNew.writeDoubleArray("/row_attrs/_Sum", data.sum);
		Metadata m = new Metadata("/row_attrs/_Sum", Metatype.NUMERIC, MetaOn.GENE, 1, data.sum.size());
		m.size = loomNew.getSizeInBytes("/row_attrs/_Sum");
		data.meta.add(m);
		data.nber_ercc = loomNew.getNbERCCs(); // For final JSON
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
    	HashSet<Integer> toKeepMap = new HashSet<>();
    	for(int s:toKeep) toKeepMap.add(s);
    	
    	for(long i = 0; i < sum.size(); i++) if(sum.get(i) == 0 || !toKeepMap.contains((int)i)) geneIndexesToFilter.add(i);
    	System.out.println(geneIndexesToFilter.size() +  " genes to filter.");
    	
    	filterGenes(loom, geneIndexesToFilter);
    	loom.close();
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
    	loom.close();
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
			float[][] subMatrix = loom.readBlock(blockSize[0], (int)dim[1], nbBlocks);
			
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
		loom.close();
    }
    
    private static void filterGenes(LoomFile loom, HashSet<Long> geneIndexesToFilter)
    {
    	// Create new Loom filtered
		System.out.println("Filtering the LOOM file...");
		long[] dim = loom.getDimensions();
		System.out.println("Current /matrix size = " + dim[0] + " x " + dim[1]);
		System.out.println("Creating New Loom file with /matrix size = " + (dim[0] - geneIndexesToFilter.size()) + " x " + dim[1]);
		
		LoomData data = new LoomData(dim[0] - geneIndexesToFilter.size(), dim[1]);
		LoomFile loomNew = new LoomFile("w", Parameters.outputFolder + "output.loom");
		
		// Process Main Matrix
		loomNew.createEmptyMatrix(data.nber_genes, data.nber_cells);
    	    	
    	// Read the file block per block according to natural storage
		long indexRowToWrite = 0;
		int[] blockSize = loom.getChunkSizes();
		int nbTotalBlocks = (int)Math.ceil((double)dim[0] / blockSize[0]);
		System.out.print("Reading " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{		
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = loom.readBlock(blockSize[0], (int)dim[1], nbBlocks);
			
			// Parsing Data and generating summary annotations
			for(int x = 0; x < subMatrix.length; x++)
			{
				long i = x + nbBlocks * blockSize[0]; // Original index
				if(!geneIndexesToFilter.contains(i))
				{
					for(int j = 0; j < subMatrix[0].length; j++)
					{
						float value = subMatrix[x][j];
							
						if(data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) data.is_count_table = false;
		        				
		        		// Re-Generate the sums
						data.depth.set(j, data.depth.get(j) + value);
						
						// Number of zeroes / detected genes 
						if(value == 0) data.nber_zeros++;
						else data.detected_genes.set(j, data.detected_genes.get(j) + 1);
					}
					// Writing this row to output
					loomNew.writeFloatRow("/matrix", subMatrix[x], indexRowToWrite);
					indexRowToWrite++;
				}
			}
		}
		System.out.println(" OK!");
		
		// Process metadata
		System.out.print("Now processing the Metadata...");
		List<Metadata> meta = loom.listMetadata();
		System.out.println(meta.size() + " metadata found");
		for(Metadata m:meta) 
		{
			if(!m.path.equals("/col_attrs/_Depth") && !m.path.equals("/col_attrs/_Detected_Genes"))
			{
				LoomFile.copyMetadata(m, geneIndexesToFilter, false, loom, loomNew);
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
		data.nber_ercc = loomNew.getNbERCCs(); // For final JSON
		loomNew.close();
		FilterCellsJSON.writeOutputJSON(data);
    }
    
    /*public static void filterEXPRESSED() throws IOException
    {
    	BufferedWriter bw_filt = new BufferedWriter(new FileWriter(Parameters.outputFolder + "filtered.genes.txt"));
    	filtJSON.nber_filtered_genes = 0;
    	filtJSON.nber_zeros = 0;
    	// Header + Stats (First Read of file)
    	getStatsOnFile(false);
    	// Analysis
    	expressedGenesPerSample = new int[filtJSON.nber_cells];
    	double threshold = Utils.quartile(rowSum, Parameters.pcKept);
    	BufferedReader br = new BufferedReader(new FileReader(Parameters.fileName));  // Read twice (because it's 50% faster, writing the file is much faster)
    	String line = br.readLine(); // Skip Header
    	line = br.readLine();
    	int nbGenes = 0;
    	while(line != null)
    	{
    		int[] arrayDetected = new int[filtJSON.nber_cells];
    		String[] tokens = line.split("\t");
    		double sum = 0;
    		int nbZeros = 0;
    		for(int i = 1; i < tokens.length; i++)
    		{
    			double val = Double.parseDouble(tokens[i]);
    			if(val == 0) nbZeros++;
    			else arrayDetected[i-1]=1;
    			sum += val;
    		}
    		if(sum > threshold) 
    		{
    			if(parsJSON.loom.is_count_table) for(int i = 0; i < (tokens.length - 1); i++) expressedGenesPerSample[i] += arrayDetected[i];
    			else for(int i = 0; i < (tokens.length - 1); i++) expressedGenesPerSample[i]++;
    			filtJSON.nber_zeros += nbZeros;
    			bw.write(line + "\n");
    		}
       		else 
       		{
       			bw_filt.write(geneNames[nbGenes] + "\n");
       			filtJSON.nber_filtered_genes++;
       		}
    		nbGenes++;
            line = br.readLine();
        }
       	yCol = "Nb Expressed Genes [reads > 0]";
       	filtJSON.info = filtJSON.nber_filtered_genes + " genes where sum of expression <= " + Utils.format(threshold) + " across " + filtJSON.nber_cells + " cells [mean <= "+ Utils.format(threshold / filtJSON.nber_cells) +"] were filtered out.";
       	bw_filt.close();
       	br.close();
    }
    
    public static void filterVAR() throws IOException
    {
    	BufferedWriter bw_filt = new BufferedWriter(new FileWriter(Parameters.outputFolder + "filtered.genes.txt"));
    	filtJSON.nber_filtered_genes = 0;
    	filtJSON.nber_zeros = 0;
    	// Header + Stats (First Read of file)
    	getStatsOnFile(false);
    	// Analysis
    	expressedGenesPerSample = new int[filtJSON.nber_cells];
    	double threshold = Utils.quartile(rowVar, Parameters.pcKept);
    	BufferedReader br = new BufferedReader(new FileReader(Parameters.fileName));  // Read twice (because it's 50% faster, writing the file is much faster)
    	String line = br.readLine(); // Skip Header
    	line = br.readLine();
    	int nbGenes = 0;
    	while(line != null)
    	{
    		int[] arrayDetected = new int[filtJSON.nber_cells];
    		String[] tokens = line.split("\t");
    		int nbZeros = 0;
    		double mean = 0;
    	    double M2 = 0;
    		for(int i = 1; i < tokens.length; i++)
    		{
    			double val = Double.parseDouble(tokens[i]);
    			if(val == 0) nbZeros++;
    			else arrayDetected[i-1]=1;
    			double delta = val - mean;
    			mean = mean + delta/i;
    			M2 = M2 + delta*(val - mean);
    		}
    		double var = (M2 / (tokens.length - 2)); // M2 / tokens.length == Variance
    		if(var > threshold) 
    		{
    			if(parsJSON.loom.is_count_table) for(int i = 0; i < (tokens.length - 1); i++) expressedGenesPerSample[i] += arrayDetected[i];
    			else for(int i = 0; i < (tokens.length - 1); i++) expressedGenesPerSample[i]++;
    			filtJSON.nber_zeros += nbZeros;
    			bw.write(line + "\n");
    		}
       		else 
       		{
       			bw_filt.write(geneNames[nbGenes] + "\n");
       			filtJSON.nber_filtered_genes++;
       		}
    		nbGenes++;
            line = br.readLine();
        }
       	yCol = "Nb Expressed Genes [reads > 0]";
       	filtJSON.info = filtJSON.nber_filtered_genes + " genes where variance <= " + Utils.format(threshold) + " across " + filtJSON.nber_cells + " cells were filtered out.";
       	bw_filt.close();
       	br.close();
    }
    
    public static void filterCOEFFOFVAR() throws IOException
    {
    	BufferedWriter bw_filt = new BufferedWriter(new FileWriter(Parameters.outputFolder + "filtered.genes.txt"));
    	filtJSON.nber_filtered_genes = 0;
    	filtJSON.nber_zeros = 0;
    	// Header + Stats (First Read of file)
    	getStatsOnFile(false);
    	// Analysis
    	expressedGenesPerSample = new int[filtJSON.nber_cells];
     	double threshold = Utils.quartile(rowCoeffOfVar, Parameters.pcKept);
    	BufferedReader br = new BufferedReader(new FileReader(Parameters.fileName));  // Read twice (because it's 50% faster, writing the file is much faster)
    	String line = br.readLine(); // Skip Header
    	line = br.readLine();
    	int nbGenes = 0;
    	while(line != null)
    	{
    		int[] arrayDetected = new int[filtJSON.nber_cells];
    		String[] tokens = line.split("\t");
    		int nbZeros = 0;
    		double mean = 0;
    	    double M2 = 0;
    		for(int i = 1; i < tokens.length; i++)
    		{
    			double val = Double.parseDouble(tokens[i]);
    			if(val == 0) nbZeros++;
    			else arrayDetected[i-1]=1;
    			double delta = val - mean;
    			mean = mean + delta/i;
    			M2 = M2 + delta*(val - mean);
    		}
    		double ecartType = Math.sqrt(M2 / (tokens.length - 2)); // M2 / (tokens.length - 2) == Variance
    		double coeffOfVar = ecartType / mean;
    		if(mean == 0) coeffOfVar = 0;
    		if(coeffOfVar > threshold)
    		{
    			if(parsJSON.loom.is_count_table) for(int i = 0; i < (tokens.length - 1); i++) expressedGenesPerSample[i] += arrayDetected[i];
    			else for(int i = 0; i < (tokens.length - 1); i++) expressedGenesPerSample[i]++;
    			filtJSON.nber_zeros += nbZeros;
    			bw.write(line + "\n");
    		}
       		else 
       		{
       			bw_filt.write(geneNames[nbGenes] + "\n");
       			filtJSON.nber_filtered_genes++;
       		}
    		nbGenes++;
            line = br.readLine();
        }
       	yCol = "Nb Expressed Genes [reads > 0]";
       	filtJSON.info = filtJSON.nber_filtered_genes + " genes where coefficient of variation <= " + Utils.format(threshold) + " across " + filtJSON.nber_cells + " cells were filtered out.";
       	bw_filt.close();
       	br.close();
    }
    
    public static void filterSCLVM() throws IOException
    {
    	BufferedWriter bw_filt = new BufferedWriter(new FileWriter(Parameters.outputFolder + "filtered.genes.txt"));
    	filtJSON.nber_filtered_genes = 0;
    	filtJSON.nber_zeros = 0;
    	// Header + Stats (First Read of file)
    	double[][] dataset = getStatsOnFile(true);
    	// Generate sizeFactors & normalize
    	double[] sizeFactors = new double[filtJSON.nber_cells];
    	for(int i = 0; i < filtJSON.nber_cells; i++)
    	{
    		// TODO "every gene contains at least one zero"
    		ArrayList<Double> sub = new ArrayList<Double>();
    		for(int j = 0; j < parsJSON.loom.nber_genes; j++) if(Double.isFinite(loggeomeans[j]) && dataset[i][j] > 0) sub.add(Math.log(dataset[i][j]) - loggeomeans[j]);
    		sizeFactors[i] = Math.exp(Utils.median(sub));
    		for(int j = 0; j < parsJSON.loom.nber_genes; j++) dataset[i][j] /= sizeFactors[i]; // Normalize original dataset
    	}
    	System.out.println();
    	
    	Stats.startRHandle(true);
    	Stats.scLVM(dataset);
    	Stats.stopRHandle();

    	/*expressedGenesPerSample = new int[filtJSON.nber_cells];
     	double threshold = Utils.quartile(rowCoeffOfVar, Parameters.pcKept);
    	BufferedReader br = new BufferedReader(new FileReader(Parameters.fileName));  // Read twice (because it's 50% faster, writing the file is much faster)
    	String line = br.readLine(); // Skip Header
    	line = br.readLine();
    	int nbGenes = 0;
    	while(line != null)
    	{
    		int[] arrayDetected = new int[filtJSON.nber_cells];
    		String[] tokens = line.split("\t");
    		int nbZeros = 0;
    		double mean = 0;
    	    double M2 = 0;
    		for(int i = 1; i < tokens.length; i++)
    		{
    			double val = Double.parseDouble(tokens[i]);
    			if(val == 0) nbZeros++;
    			else arrayDetected[i-1]=1;
    			double delta = val - mean;
    			mean = mean + delta/i;
    			M2 = M2 + delta*(val - mean);
    		}
    		double ecartType = Math.sqrt(M2 / (tokens.length - 2)); // M2 / (tokens.length - 2) == Variance
    		double coeffOfVar = ecartType / mean;
    		if(mean == 0) coeffOfVar = 0;
    		if(coeffOfVar > threshold)
    		{
    			if(parsJSON.is_count_table) for(int i = 0; i < (tokens.length - 1); i++) expressedGenesPerSample[i] += arrayDetected[i];
    			else for(int i = 0; i < (tokens.length - 1); i++) expressedGenesPerSample[i]++;
    			filtJSON.nber_zeros += nbZeros;
    			bw.write(line + "\n");
    		}
       		else 
       		{
       			bw_filt.write(geneNames[nbGenes] + "\n");
       			filtJSON.nber_filtered_genes++;
       		}
    		nbGenes++;
            line = br.readLine();
        }
       	yCol = "Nb Expressed Genes [reads > 0]";
       //	filtJSON.info = filtJSON.nber_filtered_genes + " genes where coefficient of variation <= " + Utils.format(threshold) + " across " + filtJSON.nber_cells + " cells were filtered out.";
       	bw_filt.close();
       	//br.close();
    }*/
    
   /* public static double[][] getStatsOnFile(boolean getdataset) throws IOException
    {
    	BufferedReader br = new BufferedReader(new FileReader(Parameters.fileName));
    	String line = br.readLine(); // header
    	readHeader(line);
    	double [][] dataset = null;
    	if(Parameters.nbCellsDetected > filtJSON.nber_cells) new ErrorJSON("'Min Detected' should be smaller than the total number of cells/samples i.e. <=" + filtJSON.nber_cells);
    	if(getdataset) dataset = new double[filtJSON.nber_cells][(int)parsJSON.loom.nber_genes];// TODO Error if data too big
    	geneNames = new String[(int)parsJSON.loom.nber_genes];// TODO Error if data too big
    	colSum = new double[(int)filtJSON.nber_cells]; // TODO Error if data too big
    	rowSum = new double[(int)parsJSON.loom.nber_genes];// TODO Error if data too big
    	rowVar = new double[(int)parsJSON.loom.nber_genes];// TODO Error if data too big
    	rowCoeffOfVar = new double[(int)parsJSON.loom.nber_genes];// TODO Error if data too big
    	loggeomeans = new double[(int)parsJSON.loom.nber_genes]; // scLVM// TODO Error if data too big
    	sizeFactors = new double[filtJSON.nber_cells]; // scLVM
    	int nbGenes = 0;
    	line = br.readLine();
    	while(line != null)
    	{
    		String[] tokens = line.split("\t");
    		double mean = 0;
    	    double M2 = 0;
    		geneNames[nbGenes] = tokens[0];
     		for(int i = 1; i < tokens.length; i++)
    		{
    			double val = Double.parseDouble(tokens[i]);
    			double delta = val - mean;
    			mean = mean + delta/i;
    			M2 = M2 + delta*(val - mean);
    			colSum[i-1] += val;
    			rowSum[nbGenes] += val;
    			loggeomeans[nbGenes] += Math.log(val);
    			if(getdataset) dataset[i-1][nbGenes] = val;
    		}
    		loggeomeans[nbGenes] /= (tokens.length - 1);
    		rowVar[nbGenes] = M2 / (tokens.length - 2);
    		if(mean == 0) rowCoeffOfVar[nbGenes] = 0;
    		else rowCoeffOfVar[nbGenes] = Math.sqrt(rowVar[nbGenes]) / mean;
    		nbGenes++;
            line = br.readLine();
        }
    	br.close();
    	filtJSON.nber_genes = nbGenes;
    	if(filtJSON.nber_genes != parsJSON.loom.nber_genes) new ErrorJSON("Detected different number of genes between parsingJSON("+parsJSON.loom.nber_genes+") and Data Matrix("+filtJSON.nber_genes+")");
    	return dataset;
    }*/
}
