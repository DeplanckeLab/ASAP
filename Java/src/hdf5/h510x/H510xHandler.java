package hdf5.h510x;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import bigarrays.StringArray64;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import hdf.hdf5lib.exceptions.HDF5JavaException;
import hdf5.ProgressiveReader;
import hdf5.loom.LoomFile;
import json.ErrorJSON;
import json.ParsingJSON;
import json.PreparsingJSON;
import model.MapGene;
import model.Parameters;
import parsing.model.GroupPreparse;

public class H510xHandler 
{
	private static HashSet<String> validGroups = null;
	
	public static boolean is10xFormatOK()
	{
		validGroups = new HashSet<>();
		IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName);
		List<String> groups = reader.object().getAllGroupMembers("/");
		for(String g:groups) 
		{
			boolean valid = true;
			if(!reader.object().isDataSet("/" + g + "/barcodes")) valid = false;
			if(!reader.object().isDataSet("/" + g + "/data")) valid = false;
			if(!reader.object().isDataSet("/" + g + "/indices")) valid = false;
			if(!reader.object().isDataSet("/" + g + "/indptr")) valid = false;
			if(!reader.object().isDataSet("/" + g + "/shape")) valid = false;
			if(!reader.object().isDataSet("/" + g + "/gene_names") && !reader.object().isDataSet("/" + g + "/genes")) 
			{
				if(!reader.object().isDataSet("/" + g + "/features/id")) valid = false; //sc-ATAC-seq
			}
			if(valid) validGroups.add(g);
		}
		reader.close();
		if(validGroups.isEmpty()) return false;
		return true;
	}
	
	private static StringArray64 readStringArray(IHDF5Reader reader, String path) // Because it is not a Loom file
	{
		if(reader == null) new ErrorJSON("Please open the HDF5 file first");
		long length = reader.getDataSetInformation(path).getDimensions()[0];
		StringArray64 res = new StringArray64(length); // Know the actual size of the array, and create it sufficiently big
		int nbChunks = (int)(length / StringArray64.chunkSize()) + 1;
		for(int i = 0; i < nbChunks; i++) res.set(i, reader.string().readArrayBlock(path, StringArray64.chunkSize(), i));
		return res;
	}
	
	public static void parse()
	{
		try
		{
			System.out.println("Parsing the 10x file...");
			if(Parameters.selection == null) new ErrorJSON("No selection was made for parsing H5_10x file. It is mandatory.");
			ParsingJSON json = new ParsingJSON(); // Do NOT init the matrix in RAM
			IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName); // Here, it is not a Loom file
			
			// First handle the gene annotation
			if(reader.object().isDataSet("/" + Parameters.selection + "/gene_names")) // scRNA-seq
			{
				json.data.gene_names = readStringArray(reader, "/" + Parameters.selection + "/gene_names");
				json.data.ens_names = readStringArray(reader, "/" + Parameters.selection + "/genes");
				json.data.original_gene_names = json.data.gene_names.copy();
			}
			else // scATAC-seq
			{
				json.data.gene_names = readStringArray(reader, "/" + Parameters.selection + "/features/id");
				json.data.ens_names = new StringArray64(json.data.gene_names.size());
				json.data.original_gene_names = json.data.gene_names.copy();
				// /features/feature_type = DISCRETE [peaks, ...?]
				// /features/id = peak name
				// /features/name = peak name
			}
			MapGene.parseGenes(json.data); // Add Info from DB
			
			// Then handle the cell names
			json.data.cell_names = readStringArray(reader, "/" + Parameters.selection + "/barcodes");
			for (int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, json.data.cell_names.get(i).trim().replaceAll("'|\"", ""));
			json.data.is_count_table = true; // It is always count data from 10x datasets
			
			// Prepare the writing in the file
			// Since this format is not compatible with our chunking specs, we first create a chunked matrix
			json.loom.createEmptyFloat32MatrixDataset("/matrix", json.data.nber_genes, json.data.nber_cells, Parameters.defaultChunkX, Parameters.defaultChunkY);
			
			// Now read the file block per block (1000 cells at a time, all genes)
			ProgressiveReader pr = new ProgressiveReader(reader, json.data.nber_genes, json.data.nber_cells, Parameters.defaultChunkX, Parameters.defaultChunkY);
			
			// Use indptr to compute the number of detected genes
			for(int i = 0; i < pr.indptr.size() - 1; i++) json.data.detected_genes.set(i, (int)(pr.indptr.get(i + 1) - pr.indptr.get(i)));
			
			for(int nbBlocks = 0; nbBlocks < pr.nbTotalBlocks; nbBlocks++)
			{
				// Retrieve the blocks that will contain all columns (because we write per blocks of 1000 cells x all genes)
				//System.out.println("Block " + (nbBlocks + 1));
				
				// First find out the index of the columns
				long k = (nbBlocks + 1) * Parameters.defaultChunkX;
				if( k > json.data.nber_cells) k = json.data.nber_cells;
				
				// Then extract the necessary data (& convert to dense matrix in the process)
				float[][] subMatrix = pr.readSubMatrix(nbBlocks * Parameters.defaultChunkX, k, json, true);
						
				// Writing this merged block to output
				json.loom.writeFloatBlockDataset("/matrix", subMatrix, 0, nbBlocks);
			}
			
			json.loom.resizeDataset("/matrix", json.data.nber_genes, json.data.nber_cells); // Cause writing fixed-size blocks can extend the matrix size with 0
			
			long nbNonZeroValues = reader.getDataSetInformation("/" + Parameters.selection + "/data").getDimensions()[0]; // Know the actual size of the array

			reader.close();
			
			json.data.nber_zeros = json.data.nber_cells * json.data.nber_genes - nbNonZeroValues;
			LoomFile.fillLoomFile(json.loom, json.data);
			json.writeOutputJSON();
			json.loom.close();
		}
		catch(HDF5JavaException e)
		{
			new ErrorJSON(e.getMessage());
		}
	}
	
    public static void preparse()
    {
    	System.out.println("10x file was detected, and format is OK. [" + validGroups.size() + " valid group" + (validGroups.size() == 1?"":"s") + "]");
    	ArrayList<GroupPreparse> foundDatasets = new ArrayList<GroupPreparse>();
    	for(String group:validGroups)
    	{
    		GroupPreparse g = new GroupPreparse(group);
        	IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName);
        	long[] dim = reader.int64().readArray("/" + group + "/shape");
        	g.nbGenes = dim[0];
        	g.nbCells = dim[1];
        	// Computing the dense submatrix
        	int nGenes = (int)Math.min(10, g.nbGenes);
        	int nCells = (int)Math.min(10, g.nbCells);
        	g.matrix = new float[nGenes][nCells]; // Dense matrix to be computed from sparse format
	        // Read 3 datasets required for recreating the dense matrix     
	        long[] indptr = reader.int64().readArrayBlock("/" + group + "/indptr", nCells + 1, 0); // submatrix [0, ncol+1]
	        long[] indices = reader.int64().readArrayBlock("/" + group + "/indices", (int)indptr[indptr.length - 1], 0); // submatrix [0, ncol] // TODO if the index is too big, this fails
	        int[] data = reader.int32().readArrayBlock("/" + group + "/data", (int)indptr[indptr.length - 1], 0); // TODO if the index is too big, this fails
	        // Fill the dense matrix with values != 0
	        for(int i = 0; i < nCells; i++)
			{
				for(long j = indptr[i]; j < indptr[i+1]; j++)
				{
					if(indices[(int)j] < nGenes) g.matrix[(int)indices[(int)j]][i] = data[(int)j]; // TODO if the index is too big, this fails
				}
			}
        	g.isCount = true; // It is always count data from 10x datasets
        	g.cellNames = reader.string().readArrayBlock("/" + group + "/barcodes", 10, 0);
        	for (int i = 0; i < g.cellNames.length; i++) g.cellNames[i] = g.cellNames[i].replaceAll("\"", "");
        	
			if(reader.object().isDataSet("/" + group + "/gene_names")) g.geneNames = reader.string().readArrayBlock("/" + group + "/gene_names", 10, 0); // scRNA-seq
			else g.geneNames = reader.string().readArrayBlock("/" + group + "/features/id", 10, 0); // scATAC-seq
			for (int i = 0; i < g.geneNames.length; i++) g.geneNames[i] = g.geneNames[i].replaceAll("\"", "");
        	reader.close();
        	foundDatasets.add(g);
    	}
    	PreparsingJSON.writeOutputJSON(foundDatasets);
    }
}
