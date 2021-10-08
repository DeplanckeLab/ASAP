package hdf5.h5ad;

import java.util.ArrayList;

import ch.systemsx.cisd.hdf5.HDF5CompoundDataMap;
import ch.systemsx.cisd.hdf5.HDF5CompoundMemberInformation;
import ch.systemsx.cisd.hdf5.HDF5DataSetInformation;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import hdf.hdf5lib.exceptions.HDF5JavaException;
import json.ErrorJSON;
import json.PreparsingJSON;
import model.Parameters;
import parsing.model.GroupPreparse;

public class H5ADHandler 
{
	public static boolean isH5ADFormatOK() // I only check the main sparse matrix here
	{
		IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName);
		if(!reader.exists("/X"))
		{
			if(!reader.exists("/raw.X")) return false;
			if(reader.object().isGroup("/raw.X")) 
			{
				if(!reader.object().isDataSet("/raw.X/indices")) return false;
				if(!reader.object().isDataSet("/raw.X/indptr")) return false;
				if(!reader.object().isDataSet("/raw.X/data")) return false;
			}
		}
		else if(reader.object().isGroup("/X")) 
		{
			if(!reader.object().isDataSet("/X/indices")) return false;
			if(!reader.object().isDataSet("/X/indptr")) return false;
			if(!reader.object().isDataSet("/X/data")) return false;
		}
		reader.close();
		return true;
	}
	
	/*private static StringArray64 readStringArray(IHDF5Reader reader, String path) // Because it is not a Loom file
	{
		if(reader == null) new ErrorJSON("Please open the HDF5 file first");
		long length = reader.getDataSetInformation(path).getDimensions()[0];
		StringArray64 res = new StringArray64(length); // Know the actual size of the array, and create it sufficiently big
		int nbChunks = (int)(length / StringArray64.chunkSize()) + 1;
		for(int i = 0; i < nbChunks; i++) res.set(i, reader.string().readArrayBlock(path, StringArray64.chunkSize(), i));
		return res;
	}*/
	
	public static void parse()
	{
		/*try
		{
			System.out.println("Parsing the h5ad file...");
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
				float[][] subMatrix = pr.readSubMatrix(nbBlocks * Parameters.defaultChunkX, k, json);
						
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
		}*/
	}
	
    public static void preparse()
    {
    	try
    	{
    		IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName);
	    	ArrayList<GroupPreparse> foundDatasets = new ArrayList<GroupPreparse>();
	    	
	    	// Read /raw.X and /X matrices (handle dense / sparse csr csc formats)
	    	GroupPreparse g = readGroup("/raw.X", reader);
	    	if(g != null) foundDatasets.add(g);
	    	g = readGroup("/X", reader);
	    	if(g != null) foundDatasets.add(g);
	    	
	    	
	    	// List Metadata
			/*List<String> m = reader.getGroupMembers("/row_attrs");
			for(String mem:m) if(!reader.isGroup("/row_attrs/" + mem)) g.additionalMetadataPath.add("/row_attrs/" + mem);
			m = reader.getGroupMembers("/col_attrs");
			for(String mem:m) if(!reader.isGroup("/col_attrs/" + mem)) g.additionalMetadataPath.add("/col_attrs/" + mem);
			if(reader.exists("/attrs"))
			{
				m = reader.getGroupMembers("/attrs");
				for(String mem:m) if(!reader.isGroup("/attrs/" + mem)) g.additionalMetadataPath.add("/attrs/" + mem);
			}
			if(reader.exists("/layers"))
			{
				m = reader.getGroupMembers("/layers");
				for(String mem:m) if(!reader.isGroup("/layers/" + mem)) g.additionalMetadataPath.add("/layers/" + mem);
			}
        	

        	
			if(reader.object().isDataSet("/" + group + "/gene_names")) g.geneNames = reader.string().readArrayBlock("/" + group + "/gene_names", 10, 0); // scRNA-seq
			else g.geneNames = reader.string().readArrayBlock("/" + group + "/features/id", 10, 0); // scATAC-seq
			for (int i = 0; i < g.geneNames.length; i++) g.geneNames[i] = g.geneNames[i].replaceAll("\"", "");
        	*/
        	
			// Close handle
			reader.close();
        	
        	// Write output
        	foundDatasets.add(g);
        	PreparsingJSON.writeOutputJSON(foundDatasets);
    	}
		catch(HDF5JavaException e)
		{
			new ErrorJSON(e.getMessage());
		}
    }
    
    private static GroupPreparse readGroup(String path, IHDF5Reader reader)
    {
    	GroupPreparse g = new GroupPreparse(path); // We may create one raw.X and one X dataset, depending on what is present in the .h5ad file

    	// Getting the size of the dataset from the attributes (is it always here? dunno)
    	String sparse_format = "csr";
    	String type = "DENSE";
		if(!reader.exists(path)) return null;

		if(reader.object().isGroup(path)) 
		{
			type = "SPARSE";
	    	if(reader.object().hasAttribute(path, "h5sparse_format")) sparse_format = reader.string().getAttr(path, "h5sparse_format");
	    	if(reader.object().hasAttribute(path, "h5sparse_shape"))
	    	{
	    		int[] shape = reader.int32().getArrayAttr(path, "h5sparse_shape");
	    		g.nbCells = shape[0];
	    		g.nbGenes = shape[1];
	    	}
		}
		else // DENSE
		{
	    	HDF5DataSetInformation info = reader.getDataSetInformation(path); // maxDim / Chunks / Dimensions / type float32... etc..
	    	g.nbGenes = info.getDimensions()[1];
	    	g.nbCells = info.getDimensions()[0];
		}	
		
		// TODO if here we still don't know the size of the matrix, we need to infer it from "var" or "obs" datasets
    	if(g.nbCells == 0) new ErrorJSON("Could not read the H5ad file (0 cells found. Probably because a feature is not fully implemented yet. Please contact support.)");
    	if(g.nbGenes == 0) new ErrorJSON("Could not read the H5ad file (0 genes found. Probably because a feature is not fully implemented yet. Please contact support.)");
    	if(!sparse_format.equals("csr") && !sparse_format.equals("csc")) new ErrorJSON("Could not read the H5ad file (Sparse Matrix NOT in CSR or CSC format. Reading type not implemented yet. Please contact support.)");

    	// Computing the dense submatrix
    	g.isCount = true; // Is it always count data from 10x datasets?
    	if(type.equals("DENSE")) 
    	{
	    	// Get submatrix
    		g.matrix = reader.float32().readMatrixBlock(path, 10, 10, 0, 0); // First is nb row, Second is nb col, two others are offsets I suppose?      
	    	// Is it count matrix?
    		for (int i = 0; i < g.matrix.length; i++) {
				for (int j = 0; j < g.matrix[i].length; j++) {
					if(g.matrix[i][j] != (int)g.matrix[i][j]) { g.isCount = false; break;}
				}
			}
    	}
    	else if(sparse_format.equals("csr")) // CSR
    	{
        	int nGenes = (int)Math.min(10, g.nbGenes);
        	int nCells = (int)Math.min(10, g.nbCells);
        	g.matrix = new float[nGenes][nCells]; // Dense matrix to be computed from sparse format
	        // Read 3 datasets required for recreating the dense matrix     
	        long[] indptr = reader.int64().readArrayBlock(path + "/indptr", nCells + 1, 0); // submatrix [0, ncol+1]
	        long[] indices = reader.int64().readArrayBlock(path + "/indices", (int)indptr[indptr.length - 1], 0); // submatrix [0, ncol] // TODO if the index is too big, this fails
	        float[] data = reader.float32().readArrayBlock(path + "/data", (int)indptr[indptr.length - 1], 0); // TODO if the index is too big, this fails
	        // Fill the dense matrix with values != 0
	        for(int i = 0; i < nCells; i++)
			{
				for(long j = indptr[i]; j < indptr[i+1]; j++)
				{
					if(indices[(int)j] < nGenes) 
					{
						if(data[(int)j] != (int)data[(int)j]) g.isCount = false; 
						g.matrix[(int)indices[(int)j]][i] = data[(int)j]; // TODO if the index is too big, this fails
					}
				}
			}
    	}
    	else if(sparse_format.equals("csc")) // CSC
    	{
        	int nGenes = (int)Math.min(10, g.nbGenes);
        	int nCells = (int)Math.min(10, g.nbCells);
        	g.matrix = new float[nGenes][nCells]; // Dense matrix to be computed from sparse format
	        // Read 3 datasets required for recreating the dense matrix     
	        long[] indptr = reader.int64().readArrayBlock(path + "/indptr", nCells + 1, 0); // submatrix [0, ncol+1]
	        long[] indices = reader.int64().readArrayBlock(path + "/indices", (int)indptr[indptr.length - 1], 0); // submatrix [0, ncol] // TODO if the index is too big, this fails
	        float[] data = reader.float32().readArrayBlock(path + "/data", (int)indptr[indptr.length - 1], 0); // TODO if the index is too big, this fails
	        // Fill the dense matrix with values != 0
	        for(int i = 0; i < nGenes; i++)
			{
				for(long j = indptr[i]; j < indptr[i+1]; j++)
				{
					if(indices[(int)j] < nCells) 
					{
						if(data[(int)j] != (int)data[(int)j]) g.isCount = false; 
						g.matrix[i][(int)indices[(int)j]] = data[(int)j]; // TODO if the index is too big, this fails
					}
				}
			}
    	}
		
    	// Get gene names and cell names
    	if(path.equals("/X"))
    	{
    		HDF5CompoundDataMap mapping = reader.compound().read("/obs", HDF5CompoundDataMap.class);
    		System.out.println(mapping.size());
    		for(String key:mapping.keySet()) 
    		{
    			System.out.println(key);
    			if(key.equals("index"))
    			{
    				Object o = mapping.get(key);
    				System.out.println(o);
    			}
    		}
    		/*HDF5CompoundMemberInformation[] tt = reader.compound().getDataSetInfo("/obs");
    		for(HDF5CompoundMemberInformation info:tt)
    		{
    			System.out.println(info.getType());
    			System.out.println(info.getName());
    			
    			//reader.compound().getDataSetType("/obs", , arg2)
    			System.out.println();
    			//g.cellNames = reader.compound().readArrayBlock("/obs", info.getType(), 10, 0);
    		}*/
    		//CompoundElement ds;
    		//g.cellNames = reader.compound().readArrayBlock("/obs", tt[0].getType(), 10, 0);
    		System.out.println();
    	}
    	// TODO obs.index and var.index for gene names and cell names (or raw.var)
    	/*
    	g.cellNames = reader.string().readArrayBlock("/" + group + "/barcodes", 10, 0);
    	for (int i = 0; i < g.cellNames.length; i++) g.cellNames[i] = g.cellNames[i].replaceAll("\"", "");
    	
    	// Cell names (in obs)
    	if(reader.exists("/col_attrs/CellID")) 
    	{
    		g.cellNames = reader.string().readArrayBlock("/col_attrs/CellID", 10, 0);
    		for (int i = 0; i < g.cellNames.length; i++) g.cellNames[i] = g.cellNames[i].trim().replaceAll("\"", "");
    	}
    	// Gene names (in var or row.var)
    	if(reader.exists("/row_attrs/Gene")) 
    	{
    		g.geneNames = reader.string().readArrayBlock("/row_attrs/Gene", 10, 0);
    		for (int i = 0; i < g.geneNames.length; i++) g.geneNames[i] = g.geneNames[i].trim().replaceAll("\"", "");
    	}*/
    	
		return g;
    }
}
