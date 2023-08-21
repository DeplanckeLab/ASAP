package hdf5.h5ad;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import bigarrays.StringArray64;
import ch.systemsx.cisd.hdf5.HDF5DataClass;
import ch.systemsx.cisd.hdf5.HDF5DataSetInformation;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import hdf5.ProgressiveReader;
import hdf5.loom.LoomData;
import hdf5.loom.LoomFile;
import json.ErrorJSON;
import json.ParsingJSON;
import json.PreparsingJSON;
import model.MapGene;
import model.MetaOn;
import model.Metadata;
import model.Metatype;
import model.Parameters;
import parsing.model.GroupPreparse;
import tools.Utils;

public class H5ADHandler 
{
	public static boolean isH5ADFormatOK() // I only check the main sparse matrix here
	{
		IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName);
		if(!reader.exists("/X"))
		{
			if(!reader.exists("/raw.X") && !reader.exists("/raw/X")) return false; // Backward compatible with old raw/X
			if(reader.object().isGroup("/raw.X")) 
			{
				if(!reader.object().isDataSet("/raw.X/indices")) return false;
				if(!reader.object().isDataSet("/raw.X/indptr")) return false;
				if(!reader.object().isDataSet("/raw.X/data")) return false;
			}
			else if(reader.object().isGroup("/raw/X")) 
			{
				if(!reader.object().isDataSet("/raw/X/indices")) return false;
				if(!reader.object().isDataSet("/raw/X/indptr")) return false;
				if(!reader.object().isDataSet("/raw/X/data")) return false;
			} else return false;
		}
		else if(reader.object().isGroup("/X")) 
		{
			if(!reader.object().isDataSet("/X/indices")) return false;
			if(!reader.object().isDataSet("/X/indptr")) return false;
			if(!reader.object().isDataSet("/X/data")) return false;
		}
		if(!reader.exists("/var")) return false;
		if(!reader.exists("/obs")) return false;
		reader.close();
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
		System.out.println("Parsing the h5ad file...");
		if(Parameters.selection.endsWith("/")) Parameters.selection = Parameters.selection.substring(0, Parameters.selection.length() - 1);
		
		ParsingJSON json = new ParsingJSON(); // Do NOT init the matrix in RAM
		IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName); // Here, it is not a Loom file
		
		// 1. First handle the gene annotation
		if(Parameters.rowNames != null)
		{
			if(!reader.exists(Parameters.rowNames)) new ErrorJSON(Parameters.rowNames + " does not exist in the input H5ad file.");
			json.data.original_gene_names = readStringArray(reader, Parameters.rowNames);
			if(json.data.original_gene_names.size() != json.data.nber_genes) new ErrorJSON(Parameters.rowNames + "(" + json.data.original_gene_names.size() + ") is not of the expected size (" + json.data.nber_genes + ")");
			json.data.ens_names = json.data.original_gene_names.copy();
			json.data.gene_names = json.data.original_gene_names.copy();
		}
		else
		{
			// Simply creates the arrays. They will be filled with 'parseGenes' function
			json.data.original_gene_names = new StringArray64(json.data.nber_genes);
			json.data.ens_names = new StringArray64(json.data.nber_genes);
			json.data.gene_names = new StringArray64(json.data.nber_genes);
		}
		MapGene.parseGenes(json.data); // Add Info from DB
		
		// 2. Then handle the cell names
		if(Parameters.colNames != null)
		{
			if(!reader.exists(Parameters.colNames)) new ErrorJSON(Parameters.colNames + " does not exist in the input H5ad file.");
			json.data.cell_names = readStringArray(reader, Parameters.colNames);
			if(json.data.cell_names.size() != json.data.nber_cells) new ErrorJSON(Parameters.colNames + "(" + json.data.cell_names.size() + ") is not of the expected size (" + json.data.nber_cells + ")");
			for (int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, json.data.cell_names.get(i).trim().replaceAll("'|\"", ""));
		}
    	else
    	{
    		json.data.cell_names = new StringArray64(json.data.nber_cells);
    		for (int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, "Cell_"+(i+1)); // Create fake names if none exist
    	}
		
		// 3. Cataloging metadata with correct size in /var /varm /varp /raw/var (genes) and /obs /obsm /obsp /raw/obs (cells)
		List<Metadata> okMetadata = new ArrayList<Metadata>();
		List<Metadata> otherMetadata = new ArrayList<Metadata>();
    	listMetadata(reader, okMetadata, otherMetadata, json.data.nber_genes, json.data.nber_cells);
		System.out.println(okMetadata.size() + " metadata imported " + Metadata.toString(okMetadata));
		System.out.println(otherMetadata.size() + " metadata NOT imported " + Metadata.toString(otherMetadata));
    	
    	// 4. Looking for the main matrix
    	if(Parameters.selection == null) new ErrorJSON("'-sel' parameter is mandatory for H5ad files");    	
    	GroupPreparse mainG = checkMatrix(Parameters.selection, reader);
    	if(mainG == null) new ErrorJSON("'-sel' parameter points to a non-existing or wrongly formatted dataset: " + Parameters.selection);
    	if(mainG.nbCells != json.data.nber_cells) new ErrorJSON("'-sel "+Parameters.selection+"' has " + mainG.nbCells + " cells, which does not match the '-ncells " + json.data.nber_cells + "' parameter");
    	if(mainG.nbGenes != json.data.nber_genes) new ErrorJSON("'-sel "+Parameters.selection+"' has " + mainG.nbGenes + " genes, which does not match the '-ngenes " + json.data.nber_genes + "' parameter");
    	
    	// Write it in /matrix
    	boolean success = writeMatrix(json, reader, mainG.name, "/matrix", true);
    	if(!success) new ErrorJSON("Could not write the main Matrix (/matrix) from " + mainG.name);
    	
    	// 5. Looking for other potential matrices and writing them in layers
    	if(reader.exists("/X") && !Parameters.selection.equals("/X")) 
    	{
    		success = writeMatrix(json, reader, "/X", "/layers/X", false);
    		// If not same size as main matrix, success = false
    		if(success) okMetadata.add(new Metadata("/layers/X", Metatype.NUMERIC, MetaOn.EXPRESSION_MATRIX, Parameters.nCells, Parameters.nGenes));
    	}
    	if(reader.exists("/raw.X") && !Parameters.selection.equals("/raw.X")) 
    	{
    		success = writeMatrix(json, reader, "/raw.X", "/layers/raw.X", false);
    		// If not same size as main matrix, success = false
    		if(success) okMetadata.add(new Metadata("/layers/raw.X", Metatype.NUMERIC, MetaOn.EXPRESSION_MATRIX, Parameters.nCells, Parameters.nGenes));
    	}
    	if(reader.exists("/raw/X") && !Parameters.selection.equals("/raw/X"))
    	{
    		success = writeMatrix(json, reader, "/raw/X", "/layers/raw/X", false);
    		// If not same size as main matrix, success = false
    		if(success) okMetadata.add(new Metadata("/layers/raw/X", Metatype.NUMERIC, MetaOn.EXPRESSION_MATRIX, Parameters.nCells, Parameters.nGenes));
    	}
		if(reader.exists("/layers"))
		{
			List <String> m = reader.getGroupMembers("/layers");
			for(String mem:m) if(!Parameters.selection.equals("/layers/" + mem)) 
			{
				success = writeMatrix(json, reader, "/layers/" + mem, "/layers/" + mem, false);
				if(success) okMetadata.add(new Metadata("/layers/" + mem, Metatype.NUMERIC, MetaOn.EXPRESSION_MATRIX, Parameters.nCells, Parameters.nGenes));
			}
		}
    	
    	// Should I add /uns stuff? In the main /attrs global metadata maybe?
		/*if(reader.exists("/uns"))
		{
			m = reader.getGroupMembers("/uns");
			for(String mem:m) if(!reader.isGroup("/uns/" + mem)) g.additionalMetadataPath.add("/uns/" + mem);
		}*/
    	      	
		// 6. Create metadata array and fill standard metadata (stable_id, depth, gene_names, ...)
		LoomFile.fillLoomFile(json.loom, json.data);
		
		// 7. Add other metadata
    	json.data.existing_meta = new ArrayList<>();
		for(Metadata m:okMetadata)
		{
	    	success = copyMetadata(json, reader, m);
	    	if(success) 
	    	{
	    		m = json.loom.fillInfoMetadata(m.path,true);
	    		json.data.existing_meta.add(m);
	    	}
		}
		if(json.data.existing_meta.size() == 0) json.data.existing_meta = null;
		
		// Close handle
		reader.close();				
		
		// Write output
		json.writeOutputJSON(false);
		json.loom.close();
	}
	
	public static boolean copyMetadata(ParsingJSON json, IHDF5Reader reader, Metadata m)
	{
		// If layer
		if(m.on == MetaOn.EXPRESSION_MATRIX)
		{
			// Already copied. I just get the size of the object in the created Loom
			m.size = json.loom.getSizeInBytes(m.path);
			return true;
		}
		
		// Where should we put it
		String output_path = m.path;
		String sub_name = m.path.substring(m.path.lastIndexOf("/"));
		if(m.on == MetaOn.CELL) output_path = "/col_attrs" + sub_name;
		else if(m.on == MetaOn.GENE) output_path = "/row_attrs" + sub_name;
		else new ErrorJSON("What is that? EXPRESSION_MATRIX?"); // EXPRESSION_MATRIX should already been handled

		// If it already exists, it's a system metadata? I don't want to overwrite it
		if(json.loom.exists(output_path)) 
		{
			System.out.println(output_path + " already in Loom file. Is it a system path? [Discarded]");
			return false;
		}
		
		// Check if it's a categorical metadata
		int[] codes = null;
		if(m.nbCat > 0)
		{
			if(!reader.exists(m.path + "/codes")) return false; // This should not happen
			codes = reader.readIntArray(m.path + "/codes");
			m.path = m.path + "/categories";
		}
		
		// Get dataset infos
		long[] dim = reader.getDataSetInformation(m.path).getDimensions();
		boolean flag_copied = false;
			
		// Depending on the data type, we handle it differently
		HDF5DataClass h5_class =  reader.getDataSetInformation(m.path).getTypeInformation().getDataClass();
		switch (h5_class) 
		{
			case FLOAT:	
				if(dim.length == 0) // Not an array => single value
				{
					if(codes != null) { System.out.println("Not Handled: " + m.path); return false; }
					json.loom.writeFloat(output_path, reader.readFloat(m.path));
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					if(codes != null) { System.out.println("Not Handled: " + m.path); return false; }
					float[][] values = reader.readFloatMatrix(m.path);
					json.loom.writeFloatMatrix(output_path, values);
					m.size = json.loom.getSizeInBytes(output_path);
				}
				else // Vector
				{
					float[] values = reader.readFloatArray(m.path);
					if(codes != null)
					{
						// CATEGORICAL
						float[] modified_values = new float[codes.length];
						for(int i = 0; i < modified_values.length; i++) modified_values[i] = values[codes[i]];
						m.categories = new HashSet<String>();
						for(float v:values) m.categories.add(""+v);
						values = modified_values;
					}
					json.loom.writeFloatArray(output_path, values);
					m.size = json.loom.getSizeInBytes(output_path);
				}
				flag_copied = true;
				m.path = output_path;
				break;
			case INTEGER:
				if(dim.length == 0) // Not an array => single value
				{
					if(codes != null) { System.out.println("Not Handled: " + m.path); return false; }
					json.loom.writeInt(output_path, reader.readInt(m.path));
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					if(codes != null) { System.out.println("Not Handled: " + m.path); return false; }
					int[][] values = reader.readIntMatrix(m.path);
					json.loom.writeIntMatrix(output_path, values);
					m.size = json.loom.getSizeInBytes(output_path);			
				}
				else // Vector
				{
					int[] values = reader.readIntArray(m.path);
					if(codes != null)
					{
						// CATEGORICAL
						int[] modified_values = new int[codes.length];
						for(int i = 0; i < modified_values.length; i++) modified_values[i] = values[codes[i]];
						m.categories = new HashSet<String>();
						for(int v:values) m.categories.add(""+v);
						values = modified_values;
					}
					json.loom.writeIntArray(output_path, values);
					m.size = json.loom.getSizeInBytes(output_path);
				}
				flag_copied = true;
				m.path = output_path;
				break;
			case STRING:
				if(dim.length == 0) // Not an array => single value
				{
					if(codes != null) { System.out.println("Not Handled: " + m.path); return false; }
					json.loom.writeString(output_path, reader.readString(m.path));
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					new ErrorJSON("Cannot create a 2D array of STRING");
				}
				else // Vector
				{
					String[] values = reader.readStringArray(m.path);
					if(codes != null)
					{
						// CATEGORICAL
						String[] modified_values = new String[codes.length];
						for(int i = 0; i < modified_values.length; i++) modified_values[i] = values[codes[i]];
						m.categories = new HashSet<String>();
						for(String v:values) m.categories.add(v);
						values = modified_values;
					}
					json.loom.writeStringArray(output_path, values);
					m.size = json.loom.getSizeInBytes(output_path);
				}
				flag_copied = true;
				m.path = output_path;
				break;
			default:
				m.type = Metatype.NOT_HANDLED;
				break;
		}
		return flag_copied;
	}
	
	public static boolean checkSize(IHDF5Reader reader, String path, long nbGenes, long nbCells)
	{
    	// Check type of matrix
		if(reader.object().isGroup(path)) // SPARSE 
		{
			int[] shape = null;
	    	if(reader.object().hasAttribute(path, "shape")) shape = reader.int32().getArrayAttr(path, "shape");
	    	else new ErrorJSON("Dataset " + path + " is missing the 'shape' parameter"); // It is expected in a H5AD file
	    	if(shape[0] != nbCells || shape[1] != nbGenes)
	    	{
	    		System.out.println(path + ": Not same size [Discarded]");
	    		return false;
	    	}
		}
		else // DENSE
		{
			long[] shape = reader.getDataSetInformation(path).getDimensions();
			if(shape[0] != nbCells || shape[1] != nbGenes)
			{
				System.out.println(path + ": Not same size [Discarded]");
	    		return false;
			}
		}
		return true;
	}
	
	public static boolean writeMatrix(ParsingJSON json, IHDF5Reader reader, String path, String out, boolean isMain)
	{
		if(path.endsWith("/")) path = path.substring(0, path.length() - 1);
		
		// Getting the size of the dataset from the attributes (is it always here? dunno)
    	String sparse_format = "csr"; // default
		if(!reader.exists(path)) return false;

		// Check dimensions
		boolean check = checkSize(reader, path, json.data.nber_genes, json.data.nber_cells);
		if(!check) 
		{
			if(isMain) new ErrorJSON("Dataset " + path + " is not of the expected size: " + json.data.nber_genes + " x " + json.data.nber_cells);
			return false;
		}
		
    	// Check type of matrix
		if(reader.object().isGroup(path)) // SPARSE 
		{	
			if(!reader.object().isDataSet(path + "/indices")) return false;
			if(!reader.object().isDataSet(path + "/indptr")) return false;
			if(!reader.object().isDataSet(path + "/data")) return false;
			// Properly formatted H5AD sparse matrix
	    	if(reader.object().hasAttribute(path, "h5sparse_format")) sparse_format = reader.string().getAttr(path, "h5sparse_format");
	    	else if(reader.object().hasAttribute(path, "encoding-type")) sparse_format = reader.string().getAttr(path, "encoding-type");
	    	if(sparse_format != null && sparse_format.equals("csc_matrix")) sparse_format = "csc";
	    	if(sparse_format != null && sparse_format.equals("csr_matrix")) sparse_format = "csr";
	    	if(sparse_format == null || (!sparse_format.equals("csr") && !sparse_format.equals("csc"))) new ErrorJSON("Could not read the H5ad file (Sparse Matrix NOT in CSR or CSC format. Reading type not implemented yet. Please contact support.)");
	    	// Read matrix
	    	if(sparse_format.equals("csr")) // CSR
	    	{
	    		// Create empty matrix
	    		json.loom.createEmptyFloat32MatrixDataset(out, json.data.nber_genes, json.data.nber_cells, Parameters.defaultChunkX, Parameters.defaultChunkY);
	    		
	    		// Now read the file block per block (1000 cells at a time, all genes)
				ProgressiveReader pr = new ProgressiveReader(reader, path, json.data.nber_genes, json.data.nber_cells, Parameters.defaultChunkX, Parameters.defaultChunkY);
				
				// Use indptr to compute the number of detected genes
				if(isMain) for(int i = 0; i < pr.indptr.size() - 1; i++) json.data.detected_genes.set(i, (int)(pr.indptr.get(i + 1) - pr.indptr.get(i)));
				
				// Write Matrix
				for(int nbBlocks = 0; nbBlocks < pr.nbTotalBlocks; nbBlocks++)
				{			
					// First find out the index of the columns
					long k = (nbBlocks + 1) * Parameters.defaultChunkX;
					if( k > json.data.nber_cells) k = json.data.nber_cells;
					
					// Then extract the necessary data (& convert to dense matrix in the process)
					float[][] subMatrix = pr.readSubMatrix(nbBlocks * Parameters.defaultChunkX, k, json, isMain);
					
					// Writing this merged block to output
					json.loom.writeFloatBlockDataset(out, subMatrix, 0, nbBlocks);
				}
				json.loom.resizeDataset(out, json.data.nber_genes, json.data.nber_cells); // Cause writing fixed-size blocks can extend the matrix size with 0
	    		
				if(isMain)
				{
					// Since it's sparse it's easy to know the number of 0 values
					long nbNonZeroValues = reader.getDataSetInformation(path + "/data").getDimensions()[0]; // Know the actual size of the array
					json.data.nber_zeros = json.data.nber_cells * json.data.nber_genes - nbNonZeroValues;
				}
	    	}
	    	else if(sparse_format.equals("csc")) // CSC
	    	{
	    		// TODO Handle this...
	    		if(isMain) new ErrorJSON("[FORMAT] This file contains a dataset ("+path+") in the 'CSC' sparse format, which is not handled. Please use a newer version of AnnData to submit your H5AD file.");
	    		System.out.println(path+": Dataset in the 'CSC' sparse format, which is not handled [Discarded]");
	    		return false;
	    	}
		}
		else // DENSE
		{
			LoomData data = null;
			if(isMain) data = json.data;
	    	copyFloatMatrix(path, reader, out, json.loom, data, json.data.nber_genes, json.data.nber_cells);
		}			
		return true;
	}
	
	public static void copyFloatMatrix(String from_path, IHDF5Reader from, String to_path, LoomFile to, LoomData data, long nberGenes, long nberCells)
	{	   	
		// Create empty matrix
		to.createEmptyFloat32MatrixDataset(to_path, nberGenes, nberCells, Parameters.defaultChunkX, Parameters.defaultChunkY);
		
		// Read all thing is RAM (LAZYYYYYYYYY)
		float[][] matrix = Utils.t(from.float32().readMatrix(from_path)); // Transposed
		
		// Compute stats if needed (for parsing mainly)
		if(data != null)
		{
			// Parsing Data and generating summary annotations
			for(int i = 0; i < matrix.length; i++) // Original index of gene
			{
				for(int j = 0; j < matrix[0].length; j++) // cells
				{
					float value = matrix[i][j];
					
					data.is_count_table = data.is_count_table && Utils.isInteger(value);
					
	    			// Handle biotype count per cell
	    			String biotype = data.biotypes.get(i);
	    			if(biotype.equals("protein_coding")) data.proteinCodingContent.set(j, data.proteinCodingContent.get(j) + value);
		    		else if(biotype.equals("rRNA")) data.ribosomalContent.set(j, data.ribosomalContent.get(j) + value);
	    			if(data.chromosomes.get(i).equals("MT")) data.mitochondrialContent.set(j, data.mitochondrialContent.get(j) + value);
	        			
	        		// Generate the sums
	    			data.depth.set(j, data.depth.get(j) + value);
	    			data.sum.set(i, data.sum.get(i) + value);
						
					// Number of zeroes / detected genes 
					if(value == 0) data.nber_zeros++;
					else data.detected_genes.set(j, data.detected_genes.get(j) + 1);
				}
			}
		}
		
		// Write Matrix
		to.writeFloatMatrixChunked(to_path, matrix);
	}
	
    public static void preparse()
    {
		IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName);
    	ArrayList<GroupPreparse> foundDatasets = new ArrayList<GroupPreparse>();
    	
    	// Metadata are shared for all matrices (not stored in GroupPreparse object)
		List<Metadata> okMetadata = new ArrayList<Metadata>();
		List<Metadata> otherMetadata = new ArrayList<Metadata>();
    	   		
    	// 1. First step: Looking for the main matrices
		GroupPreparse mainG = null;
    	GroupPreparse g = checkMatrix("/X", reader);
    	if(g != null) foundDatasets.add(g);
    	g = checkMatrix("/raw.X", reader);
    	if(g != null) foundDatasets.add(g);
    	g = checkMatrix("/raw/X", reader); // Backward compatible
    	if(g != null) foundDatasets.add(g);

    	// 2. Find which matrix to load as Main (/matrix)
    	if(foundDatasets.size() == 0)
    	{
			reader.close();
        	new ErrorJSON("No Matrix found in h5ad file?");
    	}
    	else if(Parameters.selection == null && foundDatasets.size() > 1)
    	{
			// If multiple matrices without any selected one, finish the preparsing
			reader.close();
        	PreparsingJSON.writeH5ADOutputJSON(foundDatasets, okMetadata, otherMetadata);
        	return;
    	}
    	else if(Parameters.selection == null && foundDatasets.size() == 1)
    	{
    		// Only one matrix to load
    		mainG = foundDatasets.get(0);
    		Parameters.selection = mainG.name;
    	}
    	else
    	{
    		// Only remaining case: Parameters.selection != null
    		if(Parameters.selection.endsWith("/")) Parameters.selection = Parameters.selection.substring(0, Parameters.selection.length() - 1);
    		for(GroupPreparse gp:foundDatasets)
    		{
    			if(gp.name.equals(Parameters.selection)) mainG = gp;
    		}
    	}
    	if(mainG == null) new ErrorJSON("The selected path is no a viable matrix in the h5ad file: " + Parameters.selection);
    	System.out.println("Found dataset: " + mainG.name);
    	
    	// 3. Cataloging metadata
    	listMetadata(reader, okMetadata, otherMetadata, mainG.nbGenes, mainG.nbCells);
    	
		// 3.5. Intermediary step: Looking for cell/gene names
		String[] colNames = null;
		String[] rowNames = null;
		if(Parameters.colNames == null)
		{
			// Try to guess the actual colnames
			if(mainG.name.startsWith("/raw"))
			{
				if(reader.object().hasAttribute("/raw/obs", "_index")) Parameters.colNames = reader.string().getAttr("/raw/obs", "_index");
				if(Parameters.colNames == null && reader.object().hasAttribute("/raw/obs", "index")) Parameters.colNames = reader.string().getAttr("/raw/obs", "index");
			}
			if(Parameters.colNames != null)
			{
				long[] dim = reader.getDataSetInformation(Parameters.colNames).getDimensions();
				if(dim[0] != mainG.nbCells) Parameters.colNames = null;
			}
			if(Parameters.colNames == null && reader.object().hasAttribute("/obs", "_index")) Parameters.colNames = reader.string().getAttr("/obs", "_index");
			if(Parameters.colNames == null && reader.object().hasAttribute("/obs", "index")) Parameters.colNames = reader.string().getAttr("/obs", "index");
			if(Parameters.colNames != null)
			{
				long[] dim = reader.getDataSetInformation(Parameters.colNames).getDimensions();
				if(dim[0] != mainG.nbCells) Parameters.colNames = null;
			}
			if(Parameters.colNames != null) colNames = reader.string().readArrayBlock(Parameters.colNames, 10, 0);
		}
		else
		{
			if(!reader.exists(Parameters.colNames)) new ErrorJSON(Parameters.colNames + " does not exist. Please change the '--col-names' parameter");
			long[] dim = reader.getDataSetInformation(Parameters.colNames).getDimensions();
			if(dim[0] != mainG.nbCells) new ErrorJSON(Parameters.colNames + " is not of the same size as the main matrix. Please change the '--col-names' parameter");
			colNames = reader.string().readArrayBlock(Parameters.colNames, 10, 0);
		}
		if(Parameters.rowNames == null)
		{
			// Try to guess the actual rownames
			if(mainG.name.startsWith("/raw"))
			{
				if(reader.object().hasAttribute("/raw/var", "_index")) Parameters.rowNames = reader.string().getAttr("/raw/var", "_index");
				if(Parameters.rowNames == null && reader.object().hasAttribute("/raw/var", "index")) Parameters.rowNames = reader.string().getAttr("/raw/var", "index");
			}
			if(Parameters.rowNames != null)
			{
				long[] dim = reader.getDataSetInformation(Parameters.rowNames).getDimensions();
				if(dim[0] != mainG.nbGenes) Parameters.rowNames = null;
			}
			if(Parameters.rowNames == null && reader.object().hasAttribute("/var", "_index")) Parameters.rowNames = reader.string().getAttr("/var", "_index");
			if(Parameters.rowNames == null && reader.object().hasAttribute("/var", "index")) Parameters.rowNames = reader.string().getAttr("/var", "index");
			if(Parameters.rowNames != null)
			{
				long[] dim = reader.getDataSetInformation(Parameters.rowNames).getDimensions();
				if(dim[0] != mainG.nbGenes) Parameters.rowNames = null;
			}
			if(Parameters.rowNames != null) rowNames = reader.string().readArrayBlock(Parameters.rowNames, 10, 0);

		}
		else
		{
			if(!reader.exists(Parameters.rowNames)) new ErrorJSON(Parameters.rowNames + " does not exist. Please change the '--row-names' parameter");
			long[] dim = reader.getDataSetInformation(Parameters.rowNames).getDimensions();
			if(dim[0] != mainG.nbGenes) new ErrorJSON(Parameters.rowNames + " is not of the same size as the main matrix. Please change the '--row-names' parameter");
			rowNames = reader.string().readArrayBlock(Parameters.rowNames, 10, 0);
		}
		
		// 4. Reading the main matrix
		foundDatasets = new ArrayList<>();
		mainG = readMatrix(mainG.name, reader, mainG.nbGenes, mainG.nbCells); // Should be the main matrix
    	if(mainG != null) 
    	{
    		mainG.cellNames = colNames;
    		mainG.geneNames = rowNames;
    		foundDatasets.add(mainG);
    	}
      	
		// Close handle
		reader.close();
    	
    	// Write output
    	PreparsingJSON.writeH5ADOutputJSON(foundDatasets, okMetadata, otherMetadata);
    }
    
    private static boolean sizeOK(Metadata m, long nbGenes, long nbCells)
    {
    	if(m.on == MetaOn.GENE && m.nbrow == nbGenes) return true;
    	if(m.on == MetaOn.CELL && m.nbcol == nbCells) return true;
    	if(m.on == MetaOn.EXPRESSION_MATRIX && m.nbcol == nbCells && m.nbrow == nbGenes) return true;
    	return false;
    }
    
    private static void _listMetadataG(String path, MetaOn on, IHDF5Reader reader, List<Metadata> okMetadata, List<Metadata> otherMetadata, long nbGenes, long nbCells)
    {
    	if(on != MetaOn.GENE && on != MetaOn.CELL) new ErrorJSON("Not handled");
    	
    	if(!path.endsWith("/")) path = path + "/"; // It's a group
    	
		List<String> m = reader.getGroupMembers(path);
		if(m.contains("__categories"))
		{
    		new ErrorJSON("[FORMAT] This file contains a group '__categories' which is an old format, which is not handled. Please use a newer version of AnnData to submit your H5AD file.");
		}
		for(String mem:m)
		{
			if(reader.isGroup(path + mem))
			{
				// Subdirectories are taken into account at only 1 depth
				List<String> m_2 = reader.getGroupMembers(path + mem);
				if(m_2.contains("categories") && m_2.contains("codes"))
				{
					Metadata tmp_m = new Metadata(path + mem);
					tmp_m.on = on; // Set before this function call (depends on /var /obs ...)
					
					// Categorical variable
					if(!reader.isGroup(path + mem + "/categories") && !reader.isGroup(path + mem + "/codes"))
					{
						long[] dim = reader.getDataSetInformation(path + mem + "/codes").getDimensions();
						
						// Set nbcol/nbrow
						if(tmp_m.on == MetaOn.GENE) 
						{
							tmp_m.nbrow = dim[0];
							if(dim.length > 1) tmp_m.nbcol = dim[1];
							else tmp_m.nbcol = 1;
						}
						else if(tmp_m.on == MetaOn.CELL)
						{
							tmp_m.nbcol = dim[0];
							if(dim.length > 1) tmp_m.nbrow = dim[1];
							else tmp_m.nbrow = 1;
						}
						
						// Retrieve extra information
						dim = reader.getDataSetInformation(path + mem + "/categories").getDimensions();
						tmp_m.nbCat = (int)dim[0];
						if(sizeOK(tmp_m, nbGenes, nbCells)) okMetadata.add(tmp_m);
						else 
						{
							System.out.println(tmp_m.path + ": Not same size [Discarded]");
							otherMetadata.add(tmp_m);
						}
					}
					else
					{
						System.out.println(tmp_m.path + ": categories/codes found, but groups? [Discarded]"); 
						otherMetadata.add(tmp_m);
					}
				}
				else if(m_2.contains("data") && m_2.contains("indices") && m_2.contains("indptr"))
				{
					Metadata tmp_m = new Metadata(path + mem);
					tmp_m.on = on; // Set before this function call (depends on /var /obs ...)
					System.out.println(tmp_m.path +  ": EXPRESSION_MATRIX found? (Not Handled for now) [Discarded]"); 
					otherMetadata.add(tmp_m);
				}
				else _listMetadata(path + mem, on, reader, okMetadata, otherMetadata, nbGenes, nbCells); // Recursive search
			}
			else
			{
				Metadata tmp_m = new Metadata(path + mem);
				tmp_m.on = on;
				
				long[] dim = reader.getDataSetInformation(path + mem).getDimensions();
				
				// Set correct nbcol/nbrow
				if(tmp_m.on == MetaOn.GENE) 
				{
					tmp_m.nbrow = dim[0];
					if(dim.length > 1) tmp_m.nbcol = dim[1];
					else tmp_m.nbcol = 1;
				}
				else if(tmp_m.on == MetaOn.CELL)
				{
					tmp_m.nbcol = dim[0];
					if(dim.length > 1) tmp_m.nbrow = dim[1];
					else tmp_m.nbrow = 1;
				}
				
				// Check size
				if(sizeOK(tmp_m, nbGenes, nbCells)) okMetadata.add(tmp_m);
				else 
				{
					System.out.println(tmp_m.path + ": Not same size [Discarded]");
					otherMetadata.add(tmp_m);
				}
			}
		}
	}
    
    private static void _listMetadata(String path, MetaOn on, IHDF5Reader reader, List<Metadata> okMetadata, List<Metadata> otherMetadata, long nbGenes, long nbCells)
    {
    	if(reader.isGroup(path)) _listMetadataG(path, on, reader, okMetadata, otherMetadata, nbGenes, nbCells);
    	else
    	{
	    	HDF5DataSetInformation inf = reader.getDataSetInformation(path);
	    	HDF5DataClass c = inf.getTypeInformation().getDataClass();
	    	if(c == HDF5DataClass.COMPOUND)
	    	{
	    		new ErrorJSON("[FORMAT] This file contains a dataset ("+path+") in the 'COMPOUND' format, which is not handled. Please use a newer version of AnnData to submit your H5AD file.");
	    	}
	    	else
	    	{
	    		new ErrorJSON("[FORMAT] This file contains datasets in a format which is not handled:" + path);
	    	}
    	}
	}
    
    private static void listMetadata(IHDF5Reader reader, List<Metadata> okMetadata, List<Metadata> otherMetadata, long nbGenes, long nbCells)
    {  	
    	// If Parsing && selection matrix is in /raw/, then we should prioritize the /raw/ metadata
    	if(Parameters.selection != null && Parameters.selection.startsWith("/raw/"))
    	{
    		// Cataloging
    		if(reader.object().exists("/raw/var")) // Mandatory in this case. For Backward compatibility
    		{
    			_listMetadata("/raw/var", MetaOn.GENE, reader, okMetadata, otherMetadata, nbGenes, nbCells);
    		}
    	}
    	
    	// Infer dimensions (Rows = genes) from var
		if(reader.object().exists("/var")) // Should exist since we've run isH5ADFormatOK()
		{
			_listMetadata("/var", MetaOn.GENE, reader, okMetadata, otherMetadata, nbGenes, nbCells);
		}
		else new ErrorJSON("/var path is not found");
		
    	// Cataloging
		if(reader.object().exists("/varm")) // Not mandatory here. I don't know if it should be mandatory?
		{
			_listMetadata("/varm", MetaOn.GENE, reader, okMetadata, otherMetadata, nbGenes, nbCells);
		}
		
		// Cataloging
		if(reader.object().exists("/varp")) // Not mandatory here. I don't know if it should be mandatory?
		{
			_listMetadata("/varp", MetaOn.GENE, reader, okMetadata, otherMetadata, nbGenes, nbCells);
		}
		
		// Cataloging
    	if(Parameters.selection == null || !Parameters.selection.startsWith("/raw/"))
    	{
			if(reader.object().exists("/raw/var")) // Not mandatory. For Backward compatibility
			{
				_listMetadata("/raw/var", MetaOn.GENE, reader, okMetadata, otherMetadata, nbGenes, nbCells);
			}
    	}
    	
    	// If Parsing && selection matrix is in /raw/, then we should prioritize the /raw/ metadata
    	if(Parameters.selection != null && Parameters.selection.startsWith("/raw/"))
    	{
    		// Cataloging
    		if(reader.object().exists("/raw/obs")) // Not mandatory. For Backward compatibility
    		{
    			_listMetadata("/raw/obs", MetaOn.CELL, reader, okMetadata, otherMetadata, nbGenes, nbCells);
    		}
    	}
    	
    	// Infer dimensions (Columns = cells) from obs
		if(reader.object().exists("/obs")) // Should exist since we've run isH5ADFormatOK()
		{
			_listMetadata("/obs", MetaOn.CELL, reader, okMetadata, otherMetadata, nbGenes, nbCells);
		} 
		else new ErrorJSON("/obs path is not found");
		
    	// Cataloging for obsm
		if(reader.object().exists("/obsm")) // Not mandatory here. I don't know if it should be mandatory?
		{
			_listMetadata("/obsm", MetaOn.CELL, reader, okMetadata, otherMetadata, nbGenes, nbCells);
		} 
				
    	// Cataloging for obsp
		if(reader.object().exists("/obsp")) // Not mandatory here. I don't know if it should be mandatory?
		{
			_listMetadata("/obsp", MetaOn.CELL, reader, okMetadata, otherMetadata, nbGenes, nbCells);
		} 
		
    	// Cataloging for /raw/obs
    	if(Parameters.selection == null || !Parameters.selection.startsWith("/raw/"))
    	{
			if(reader.object().exists("/raw/obs")) // Not mandatory here. I don't know if it should be mandatory?
			{
				_listMetadata("/raw/obs", MetaOn.CELL, reader, okMetadata, otherMetadata, nbGenes, nbCells);
			} 
    	}
    }
    
    private static GroupPreparse checkMatrix(String path, IHDF5Reader reader)
    {
		if(path.endsWith("/")) path = path.substring(0, path.length() - 1);
    
    	// Getting the size of the dataset from the attributes (is it always here? dunno)
    	String sparse_format = null; // default
		if(!reader.exists(path)) return null;

    	// Create object to return
    	GroupPreparse g = new GroupPreparse(path);
    	
    	// Check type of matrix
		if(reader.object().isGroup(path)) 
		{
			if(!reader.object().isDataSet(path + "/indices")) return null;
			if(!reader.object().isDataSet(path + "/indptr")) return null;
			if(!reader.object().isDataSet(path + "/data")) return null;
			// Properly formatted H5AD sparse matrix
	    	if(reader.object().hasAttribute(path, "h5sparse_format")) sparse_format = reader.string().getAttr(path, "h5sparse_format");
	    	else if(reader.object().hasAttribute(path, "encoding-type")) sparse_format = reader.string().getAttr(path, "encoding-type");
	    	if(sparse_format != null && sparse_format.equals("csc_matrix")) sparse_format = "csc";
	    	if(sparse_format != null && sparse_format.equals("csr_matrix")) sparse_format = "csr";
	    	if(sparse_format == null || (!sparse_format.equals("csr") && !sparse_format.equals("csc"))) return null;
	    	if(reader.object().hasAttribute(path, "shape")) 
	    	{
	    		long[] shape = reader.int64().getMDArrayAttr(path, "shape").getAsFlatArray();
	    		g.nbCells = shape[0];
	    		g.nbGenes = shape[1];
	    	}
	    	else return null;
		}
		else // DENSE
		{
    		long[] dim = reader.getDataSetInformation(path).getDimensions();
    		g.nbCells = dim[0];
    		g.nbGenes = dim[1];
		}
		return g;
    }
    
    private static GroupPreparse readMatrix(String path, IHDF5Reader reader, long nbGenes, long nbCells)
    {  	
    	if(path.endsWith("/")) path = path.substring(0, path.length() - 1);
    	
    	// Getting the size of the dataset from the attributes (is it always here? dunno)
    	String sparse_format = "csr"; // default
		if(!reader.exists(path)) return null;

    	// Create object to return
    	GroupPreparse g = new GroupPreparse(path);
		g.nbCells = nbCells;
		g.nbGenes = nbGenes;
		g.isCount = true; // Is it always count data ?
    	
    	// Check type of matrix
		if(reader.object().isGroup(path)) 
		{
			if(!reader.object().isDataSet(path + "/indices")) return null;
			if(!reader.object().isDataSet(path + "/indptr")) return null;
			if(!reader.object().isDataSet(path + "/data")) return null;
			// Properly formatted H5AD sparse matrix
	    	if(reader.object().hasAttribute(path, "h5sparse_format")) sparse_format = reader.string().getAttr(path, "h5sparse_format");
	    	else if(reader.object().hasAttribute(path, "encoding-type")) sparse_format = reader.string().getAttr(path, "encoding-type");
	    	if(sparse_format != null && sparse_format.equals("csc_matrix")) sparse_format = "csc";
	    	if(sparse_format != null && sparse_format.equals("csr_matrix")) sparse_format = "csr";
	    	if(sparse_format == null || (!sparse_format.equals("csr") && !sparse_format.equals("csc"))) new ErrorJSON("Could not read the H5ad file (Sparse Matrix NOT in CSR or CSC format. Reading type not implemented yet. Please contact support.)");
	    	if(sparse_format.equals("csr")) // CSR
	    	{
	        	int nGenes = (int)Math.min(10, g.nbGenes);
	        	int nCells = (int)Math.min(10, g.nbCells);
	        	g.matrix = new float[nGenes][nCells]; // Dense matrix to be computed from sparse format
		        // Read 3 datasets required for recreating the dense matrix     
		        long[] indptr = reader.int64().readArrayBlock(path + "/indptr", nCells + 1, 0); // submatrix [0, ncol+1]
		        long[] indices = reader.int64().readArrayBlock(path + "/indices", (int)indptr[indptr.length - 1], 0); // submatrix [0, ncol]
		        float[] data = reader.float32().readArrayBlock(path + "/data", (int)indptr[indptr.length - 1], 0);
		        // Fill the dense matrix with values != 0
		        for(int i = 0; i < nCells; i++)
				{
					for(long j = indptr[i]; j < indptr[i+1]; j++)
					{
						if(indices[(int)j] < nGenes) 
						{
							g.isCount = g.isCount && Utils.isInteger(data[(int)j]);
							g.matrix[(int)indices[(int)j]][i] = data[(int)j];
						}
					}
				}
	    	}
	    	else if(sparse_format.equals("csc")) // CSC
	    	{
	    		new ErrorJSON("[FORMAT] This file contains a dataset ("+path+") in the 'CSC' sparse format, which is not handled. Please use a newer version of AnnData to submit your H5AD file.");
	        	/**
	        	 */
	    		/*
	    		int nGenes = (int)Math.min(10, g.nbGenes);
	        	int nCells = (int)Math.min(10, g.nbCells);
	        	g.matrix = new float[nGenes][nCells]; // Dense matrix to be computed from sparse format
		        // Read 3 datasets required for recreating the dense matrix     
		        long[] indptr = reader.int64().readArrayBlock(path + "/indptr", nCells + 1, 0); // submatrix [0, ncol+1]
		        long[] indices = reader.int64().readArrayBlock(path + "/indices", (int)indptr[indptr.length - 1], 0);
		        float[] data = reader.float32().readArrayBlock(path + "/data", (int)indptr[indptr.length - 1], 0);
		        for(int i = 0; i < nGenes; i++)
				{
					for(long j = indptr[i]; j < indptr[i+1]; j++)
					{
						if(indices[(int)j] < nCells) 
						{
							g.isCount = g.isCount && Utils.isInteger(data[(int)j]);
							g.matrix[i][(int)indices[(int)j]] = data[(int)j];
						}
					}
				}
		        */
	    	}
		}
		else // DENSE
		{
	    	// Get submatrix
    		g.matrix = reader.float32().readMatrixBlock(path, 10, 10, 0, 0); // First is nb row, Second is nb col, two others are offsets I suppose?      
	    	// Is it count matrix? (Only from 10 first rows/cols is maybe not super representative :D MEN!
    		for (int i = 0; i < g.matrix.length; i++) {
				for (int j = 0; j < g.matrix[i].length; j++) {
					g.isCount = g.isCount && Utils.isInteger(g.matrix[i][j]);
					if(!g.isCount) break;
				}
			}
		}			
		return g;
    }
}
