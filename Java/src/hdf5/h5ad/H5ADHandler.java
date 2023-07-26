package hdf5.h5ad;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import bigarrays.StringArray64;
import ch.systemsx.cisd.hdf5.HDF5DataClass;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import hdf.hdf5lib.exceptions.HDF5JavaException;
import hdf5.ProgressiveReader;
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
		try
		{
			System.out.println("Parsing the h5ad file...");
			ParsingJSON json = new ParsingJSON(); // Do NOT init the matrix in RAM
			IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName); // Here, it is not a Loom file
			
			// 1. First handle the gene annotation
			if(Parameters.rowNames != null)
			{
				json.data.original_gene_names = readStringArray(reader, Parameters.rowNames);
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
				json.data.cell_names = readStringArray(reader, Parameters.colNames);
				for (int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, json.data.cell_names.get(i).trim().replaceAll("'|\"", ""));
			}
	    	else
	    	{
	    		json.data.cell_names = new StringArray64(json.data.nber_cells);
	    		for (int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, "Cell_"+(i+1)); // Create fake names if none exist
	    	}
			
			// 3. Check dimensions and cataloging metadata in /var /varm /varp /raw/var (genes) and /obs /obsm /obsp /raw/obs (cells)
			Set<Metadata> okMetadata = new HashSet<Metadata>();
			Set<Metadata> otherMetadata = new HashSet<Metadata>();
	    	listMetadata(reader, okMetadata, otherMetadata, Parameters.nGenes, Parameters.nCells);
			System.out.println(okMetadata.size() + " metadata imported " + Metadata.toString(okMetadata));
			System.out.println(otherMetadata.size() + " metadata NOT imported " + Metadata.toString(otherMetadata));
	    	
	    	// 4. Second step: Looking for the main matrix
	    	if(Parameters.selection == null)
	    	{
	    		if(reader.exists("/X")) Parameters.selection = "/X";
	    		else if(reader.exists("/raw.X")) Parameters.selection = "/raw.X";
	    		else if(reader.exists("/raw/X")) Parameters.selection = "/raw/X";
	    		else new ErrorJSON("No main Matrix found. Please specify it using '-sel'");
	    	}
	    	// Write it in /matrix
	    	boolean success = writeMatrix(json, reader, Parameters.selection, "/matrix", true);
	    	if(!success) new ErrorJSON("Could not write the main Matrix (/matrix) from " + Parameters.selection);
	    	
	    	// 5. Looking for layers and writing all other potential matrices in layers
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
	    	
	    	// TODO /uns?
			/*if(reader.exists("/uns"))
			{
				m = reader.getGroupMembers("/layers");
				for(String mem:m) if(!reader.isGroup("/layers/" + mem)) g.additionalMetadataPath.add("/layers/" + mem);
			}*/
        	      	
			// 6. Create metadata array and fill standard metadata (stable_id, depth, gene_names, ...)
			LoomFile.fillLoomFile(json.loom, json.data);
			
			// 7. Add other metadata
	    	json.data.existing_meta = new ArrayList<>();
			for(Metadata m:okMetadata)
			{
		    	success = copyMetadata(json, reader, m);
		    	if(success) json.data.existing_meta.add(m);
			}
			if(json.data.existing_meta.size() == 0) json.data.existing_meta = null;
			
			// Close handle
			reader.close();				
			
			// Write output
			json.writeOutputJSON(false);
			json.loom.close();
		}
		catch(HDF5JavaException e)
		{
			new ErrorJSON(e.getMessage());
		}
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
		
		// Get dataset infos
		long[] dim = reader.getDataSetInformation(m.path).getDimensions();
		boolean flag_copied = false;
		
		// Where should we put it
		String output_path = m.path;
		String sub_name = m.path.substring(m.path.lastIndexOf("/"));
		if(m.on == MetaOn.CELL) output_path = "/col_attrs" + sub_name;
		else if(m.on == MetaOn.GENE) output_path = "/row_attrs" + sub_name;
		else new ErrorJSON("What is that? EXPRESSION_MATRIX?"); // EXPRESSION_MATRIX should already been handled

		// If it already exists, it's a system metadata! I don't want to overwrite it
		if(json.loom.exists(output_path)) 
		{
			System.out.println(output_path + " already in Loom file. Is it a system path? [Discarded]");
			return false;
		}
			
		// Depending on the data type, we handle it differently
		HDF5DataClass h5_class =  reader.getDataSetInformation(m.path).getTypeInformation().getDataClass();
		switch (h5_class) 
		{
			case FLOAT:	
				if(dim.length == 0) // Not an array => single value
				{
					json.loom.writeFloat(output_path, reader.readFloat(m.path));
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					float[][] values = reader.readFloatMatrix(m.path);
					json.loom.writeFloatMatrix(output_path, values);
					m.size = json.loom.getSizeInBytes(output_path);			
				}
				else // Vector
				{
					float[] values = reader.readFloatArray(m.path);
					json.loom.writeFloatArray(output_path, values);
					m.size = json.loom.getSizeInBytes(output_path);
				}
				flag_copied = true;
				m.path = output_path;
				break;
			case INTEGER:
				if(dim.length == 0) // Not an array => single value
				{
					json.loom.writeInt(output_path, reader.readInt(m.path));
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					int[][] values = reader.readIntMatrix(m.path);
					json.loom.writeIntMatrix(output_path, values);
					m.size = json.loom.getSizeInBytes(output_path);			
				}
				else // Vector
				{
					int[] values = reader.readIntArray(m.path);
					json.loom.writeIntArray(output_path, values);
					m.size = json.loom.getSizeInBytes(output_path);
				}
				flag_copied = true;
				m.path = output_path;
				break;
			case STRING:
				if(dim.length == 0) // Not an array => single value
				{
					json.loom.writeString(output_path, reader.readString(m.path));
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					new ErrorJSON("Cannot create a 2D array of STRING");
				}
				else // Vector
				{
					String[] values = reader.readStringArray(m.path);
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
	
	public static boolean checkSize(IHDF5Reader reader, String path)
	{
    	// Check type of matrix
		if(reader.object().isGroup(path)) // SPARSE 
		{
			int[] shape = null;
	    	if(reader.object().hasAttribute(path, "shape")) shape = reader.int32().getArrayAttr(path, "shape");
	    	else new ErrorJSON("Dataset " + path + " is missing the 'shape' parameter"); // It is expected in a H5AD file
	    	if(shape[0] != Parameters.nCells || shape[1] != Parameters.nGenes)
	    	{
	    		System.out.println("Dataset " + path + " is not of the correct size [Ignored]");
	    		return false;
	    	}
		}
		else // DENSE
		{
			long[] shape = reader.getDataSetInformation(path).getDimensions();
			if(shape[0] != Parameters.nGenes || shape[1] != Parameters.nCells)
			{
	    		System.out.println("Dataset " + path + " is not of the correct size [Ignored]");
	    		return false;
			}
		}
		return true;
	}
	
	public static boolean writeMatrix(ParsingJSON json, IHDF5Reader reader, String path, String out, boolean isMain)
	{
		// Getting the size of the dataset from the attributes (is it always here? dunno)
    	String sparse_format = "csr"; // default
		if(!reader.exists(path)) return false;

		// Check dimensions
		boolean check = checkSize(reader, path);
		if(!check) 
		{
			if(isMain) new ErrorJSON("Dataset " + path + " is not of the expected size: " + Parameters.nGenes + " x " + Parameters.nCells);
			return false;
		}
		
		// Prepare the writing in the file
		// Since this format is not compatible with our chunking specs, we first create a chunked matrix
		json.loom.createEmptyFloat32MatrixDataset(out, json.data.nber_genes, json.data.nber_cells, Parameters.defaultChunkX, Parameters.defaultChunkY);
    	
    	// Check type of matrix
		if(reader.object().isGroup(path)) // SPARSE 
		{	
			if(!reader.object().isDataSet(path + "/indices")) return false;
			if(!reader.object().isDataSet(path + "/indptr")) return false;
			if(!reader.object().isDataSet(path + "/data")) return false;
			// Properly formatted H5AD sparse matrix
	    	if(reader.object().hasAttribute(path, "h5sparse_format")) sparse_format = reader.string().getAttr(path, "h5sparse_format");
	    	if(!sparse_format.equals("csr") && !sparse_format.equals("csc")) new ErrorJSON("Could not read the H5ad file (Sparse Matrix NOT in CSR or CSC format. Reading type not implemented yet. Please contact support.)");
	    	if(sparse_format.equals("csr")) // CSR
	    	{
	    		// Now read the file block per block (1000 cells at a time, all genes)
				ProgressiveReader pr = new ProgressiveReader(reader, path, json.data.nber_genes, json.data.nber_cells, Parameters.defaultChunkX, Parameters.defaultChunkY);
				
				// Use indptr to compute the number of detected genes
				for(int i = 0; i < pr.indptr.size() - 1; i++) json.data.detected_genes.set(i, (int)(pr.indptr.get(i + 1) - pr.indptr.get(i)));
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
	    		
	        	/*int nGenes = (int)Math.min(10, g.nbGenes);
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
							//if(data[(int)j] != (int)data[(int)j]) g.isCount = false; 
							if(Utils.isCount())
							g.matrix[(int)indices[(int)j]][i] = data[(int)j]; // TODO if the index is too big, this fails
						}
					}
				}*/
	    	}
	    	else if(sparse_format.equals("csc")) // CSC
	    	{
	        	/*int nGenes = (int)Math.min(10, g.nbGenes);
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
							//if(data[(int)j] != (int)data[(int)j]) g.isCount = false;
							if(Utils.isCount())
							g.matrix[i][(int)indices[(int)j]] = data[(int)j]; // TODO if the index is too big, this fails
						}
					}
				}*/
	    	}
		}
		else // DENSE
		{
	    	// Original file block sizes is kept (easier to perform the copy)
	    	//int[] blockSize = loom.getChunkSizes();
	    	//if(blockSize[0] / blockSize[1] > 2) LoomFile.copyFloatMatrixByCell("/matrix", loom, json.loom, json.data); // Handling HCA uneven chunks which is uncompatible with our original gene by gene process
	    	//else LoomFile.copyFloatMatrixByGene("/matrix", loom, json.loom, json.data);
	    	//loom.close();
		}			
		return true;
	}
	
    public static void preparse()
    {
    	try
    	{
    		IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName);
	    	ArrayList<GroupPreparse> foundDatasets = new ArrayList<GroupPreparse>();
	    	
	    	// Metadata are shared for all matrices (not stored in GroupPreparse object)
			Set<Metadata> okMetadata = new HashSet<Metadata>();
			Set<Metadata> otherMetadata = new HashSet<Metadata>();
	    	
			/// 1. First step: infer dimensions and cataloging metadata
	    	long[] res = listMetadata(reader, okMetadata, otherMetadata, -1, -1);
	    	long nbGenes = res[0];
	    	long nbCells = res[1];
	    			
			/// 1.5 Intermediary step: Looking for cell/gene names
			String[] colNames = null;
			String[] rowNames = null;
    		if(Parameters.colNames != null) colNames = reader.string().readArrayBlock(Parameters.colNames, 10, 0);
    		if(Parameters.rowNames != null) rowNames = reader.string().readArrayBlock(Parameters.rowNames, 10, 0);
			
	    	/// 2. Second step: Looking for the main matrices
	    	GroupPreparse g = null;
	    	g = readMatrix("/X", reader, nbGenes, nbCells); // Should be the main matrix
	    	if(g != null) 
	    	{
	    		g.cellNames = colNames;
	    		g.geneNames = rowNames;
	    		foundDatasets.add(g);
	    	}
	    	g = readMatrix("/raw.X", reader, nbGenes, nbCells);
	    	if(g != null) 
	    	{
	    		g.cellNames = colNames;
	    		g.geneNames = rowNames;
	    		foundDatasets.add(g);
	    	}
	    	g = readMatrix("/raw/X", reader, nbGenes, nbCells); // Backward compatible
	    	if(g != null) 
	    	{
	    		g.cellNames = colNames;
	    		g.geneNames = rowNames;
	    		foundDatasets.add(g);
	    	}
	    	        	      	
			// Close handle
			reader.close();
        	
        	// Write output
        	PreparsingJSON.writeH5ADOutputJSON(foundDatasets, okMetadata, otherMetadata);
    	}
		catch(HDF5JavaException e)		{
			new ErrorJSON(e.getMessage());
		}
    }
    
    private static long[] listMetadata(IHDF5Reader reader, Set<Metadata> okMetadata, Set<Metadata> otherMetadata, long nbGenes, long nbCells)
    {  	
    	// Infer dimensions (Rows = genes) from var
		if(reader.object().exists("/var")) // Should exist since we've run isH5ADFormatOK()
		{
			List<String> m = reader.getGroupMembers("/var");
			for(String mem:m)
			{
				if(reader.isGroup("/var/" + mem))
				{
					// Subdirectories are not taken into account (no recursive DFS here, I only list its content)
					List<String> m_2 = reader.getGroupMembers("/var/" + mem);
					for(String mem_2:m_2)
					{
						String p_2 = "/var/" + mem + "/" + mem_2;
						Metadata tmp_gene_m = new Metadata(p_2);
						if(reader.isGroup(p_2)) // I don't list this directory content
						{
							tmp_gene_m.nbrow = -1;
							tmp_gene_m.nbcol = -1;
						}
						else
						{
							long[] dim = reader.getDataSetInformation(p_2).getDimensions();
							tmp_gene_m.nbrow = dim[0];
							if(dim.length > 1) tmp_gene_m.nbcol = dim[1];
							else tmp_gene_m.nbcol = 1;
						}
						tmp_gene_m.on = MetaOn.GENE;
						otherMetadata.add(tmp_gene_m);
					}
				}
				else
				{
					Metadata tmp_gene_m = new Metadata("/var/" + mem); // I take the opportunity to generate the list of metadata :)
					long[] dim = reader.getDataSetInformation("/var/" + mem).getDimensions();
					if(dim.length != 1) new ErrorJSON("Type of '/var/" + mem + "' is not a vector (1D)?");
					if(nbGenes == -1) nbGenes = dim[0];
					else if(nbGenes != dim[0]) new ErrorJSON("Dimension of '/var/" + mem + "' should be " + nbGenes + ", it's " + dim[0]);
					tmp_gene_m.on = MetaOn.GENE;
					tmp_gene_m.nbcol = 1; // It's a vector
					tmp_gene_m.nbrow = dim[0];
					okMetadata.add(tmp_gene_m); 
				}
			}
		} else new ErrorJSON("/var path is not found");
		if(nbGenes == -1) new ErrorJSON("'/var' is empty. Could not estimate number of features (e.g. genes).");
		
    	// Cataloging and checking dimensions (Columns = cells) for varm
		if(reader.object().exists("/varm")) // Not mandatory here. I don't know if it should be mandatory?
		{
			List<String> m = reader.getGroupMembers("/varm");
			for(String mem:m)
			{
				if(reader.isGroup("/varm/" + mem))
				{
					// Subdirectories are not taken into account (no recursive DFS here, I only list its content)
					List<String> m_2 = reader.getGroupMembers("/varm/" + mem);
					for(String mem_2:m_2)
					{
						String p_2 = "/varm/" + mem + "/" + mem_2;
						Metadata tmp_gene_m = new Metadata(p_2);
						if(reader.isGroup(p_2)) // I don't list this directory content
						{
							tmp_gene_m.nbrow = -1;
							tmp_gene_m.nbcol = -1;
						}
						else
						{
							long[] dim = reader.getDataSetInformation(p_2).getDimensions();
							tmp_gene_m.nbrow = dim[0];
							if(dim.length > 1) tmp_gene_m.nbcol = dim[1];
							else tmp_gene_m.nbcol = 1;
						}
						tmp_gene_m.on = MetaOn.GENE;
						otherMetadata.add(tmp_gene_m);
					}
				}
				else
				{
					Metadata tmp_gene_m = new Metadata("/varm/" + mem); // I generate the list of metadata
					long[] dim = reader.getDataSetInformation("/varm/" + mem).getDimensions();
					if(dim.length != 2) new ErrorJSON("Type of '/varm/" + mem + "' is not a matrix (2D)?");
					if(nbGenes != dim[0]) new ErrorJSON("Dimension of '/varm/" + mem + "' should be " + nbGenes + " x " + dim[1] + " but it's " + dim[0] + " x " + dim[1]);
					tmp_gene_m.on = MetaOn.GENE;
					tmp_gene_m.nbcol = dim[1];
					tmp_gene_m.nbrow = dim[0];
					okMetadata.add(tmp_gene_m); 
				}
			}
		}
		
		// Cataloging and ignoring metadata for varp
		if(reader.object().exists("/varp")) // Not mandatory here. I don't know if it should be mandatory?
		{
			List<String> m = reader.getGroupMembers("/varp");
			for(String mem:m) // Everything in /varp is ignored for parsing
			{
				if(reader.isGroup("/varp/" + mem))
				{
					// Subdirectories are not taken into account (no recursive DFS here, I only list its content)
					List<String> m_2 = reader.getGroupMembers("/varp/" + mem);
					for(String mem_2:m_2)
					{
						String p_2 = "/varp/" + mem + "/" + mem_2;
						Metadata tmp_gene_m = new Metadata(p_2);
						if(reader.isGroup(p_2)) // I don't list this directory content
						{
							tmp_gene_m.nbrow = -1;
							tmp_gene_m.nbcol = -1;
						}
						else
						{
							long[] dim = reader.getDataSetInformation(p_2).getDimensions();
							tmp_gene_m.nbrow = dim[0];
							if(dim.length > 1) tmp_gene_m.nbcol = dim[1];
							else tmp_gene_m.nbcol = 1;
						}
						tmp_gene_m.on = MetaOn.GENE;
						otherMetadata.add(tmp_gene_m);
					}
				}
				else
				{
					Metadata tmp_gene_m = new Metadata("/varp/" + mem); // I generate the list of metadata
					long[] dim = reader.getDataSetInformation("/varp/" + mem).getDimensions();
					if(dim.length != 2) new ErrorJSON("Type of '/varp/" + mem + "' is not a matrix (2D)?");
					tmp_gene_m.on = MetaOn.GENE;
					tmp_gene_m.nbcol = dim[1];
					tmp_gene_m.nbrow = dim[0];
					otherMetadata.add(tmp_gene_m); 
				}
			}
		}
		
		// Cataloging and ignoring metadata for /raw/var
		if(reader.object().exists("/raw/var/")) // Not mandatory. For Backward compatibility
		{
			List<String> m = reader.getGroupMembers("/raw/var/");
			for(String mem:m) // Everything in /raw/var is ignored for parsing
			{
				if(reader.isGroup("/raw/var/" + mem))
				{
					// Subdirectories are not taken into account (no recursive DFS here, I only list its content)
					List<String> m_2 = reader.getGroupMembers("/raw/var/" + mem);
					for(String mem_2:m_2)
					{
						String p_2 = "/raw/var/" + mem + "/" + mem_2;
						Metadata tmp_gene_m = new Metadata(p_2);
						if(reader.isGroup(p_2)) // I don't list this directory content
						{
							tmp_gene_m.nbrow = -1;
							tmp_gene_m.nbcol = -1;
						}
						else
						{
							long[] dim = reader.getDataSetInformation(p_2).getDimensions();
							tmp_gene_m.nbrow = dim[0];
							if(dim.length > 1) tmp_gene_m.nbcol = dim[1];
							else tmp_gene_m.nbcol = 1;
						}
						tmp_gene_m.on = MetaOn.GENE;
						otherMetadata.add(tmp_gene_m);
					}
				}
				else
				{
					Metadata tmp_gene_m = new Metadata("/raw/var/" + mem); // I generate the list of metadata
					long[] dim = reader.getDataSetInformation("/raw/var/" + mem).getDimensions();
					if(dim.length != 1) new ErrorJSON("Type of '/raw/var/" + mem + "' is not a vector (1D)?");
					tmp_gene_m.on = MetaOn.GENE;
					tmp_gene_m.nbcol = 1;
					tmp_gene_m.nbrow = dim[0];
					otherMetadata.add(tmp_gene_m); 
				}
			}
		}
		
    	// Infer dimensions (Columns = cells) from obs
		if(reader.object().exists("/obs")) // Should exist since we've run isH5ADFormatOK()
		{
			List<String> m = reader.getGroupMembers("/obs");
			for(String mem:m)
			{
				if(reader.isGroup("/obs/" + mem))
				{
					// Subdirectories are not taken into account (no recursive DFS here, I only list its content)
					List<String> m_2 = reader.getGroupMembers("/obs/" + mem);
					for(String mem_2:m_2)
					{
						String p_2 = "/obs/" + mem + "/" + mem_2;
						Metadata tmp_cell_m = new Metadata(p_2);
						if(reader.isGroup(p_2)) // I don't list this directory content
						{
							tmp_cell_m.nbrow = -1;
							tmp_cell_m.nbcol = -1;
						}
						else
						{
							long[] dim = reader.getDataSetInformation(p_2).getDimensions();
							tmp_cell_m.nbcol = dim[0];
							if(dim.length > 1) tmp_cell_m.nbrow = dim[1];
							else tmp_cell_m.nbrow = 1;
						}
						tmp_cell_m.on = MetaOn.CELL;
						otherMetadata.add(tmp_cell_m);
					}
				}
				else
				{
					Metadata tmp_cell_m = new Metadata("/obs/" + mem); // I take the opportunity to generate the list of metadata :)
					long[] dim = reader.getDataSetInformation("/obs/" + mem).getDimensions();
					if(dim.length != 1) new ErrorJSON("Type of '/obs/" + mem + "' is not a vector (1D)?");
					if(nbCells == -1) nbCells = dim[0];
					else if(nbCells != dim[0]) new ErrorJSON("Dimension of '/obs/" + mem + "' should be " + nbCells + ", it's " + dim[0]);
					tmp_cell_m.on = MetaOn.CELL;
					tmp_cell_m.nbcol = dim[0];
					tmp_cell_m.nbrow = 1; // I's a vector
					okMetadata.add(tmp_cell_m); // I take the opportunity to generate the list of metadata :)
				}
			}
		} else new ErrorJSON("/obs path is not found");
		if(nbCells == -1) new ErrorJSON("'/obs' is empty. Could not estimate number of cells.");
		
    	// Cataloging and checking dimensions (Columns = cells) for obsm
		if(reader.object().exists("/obsm")) // Not mandatory here. I don't know if it should be mandatory?
		{
			List<String> m = reader.getGroupMembers("/obsm");
			for(String mem:m)
			{
				if(reader.isGroup("/obsm/" + mem))
				{
					// Subdirectories are not taken into account (no recursive DFS here, I only list its content)
					List<String> m_2 = reader.getGroupMembers("/obsm/" + mem);
					for(String mem_2:m_2)
					{
						String p_2 = "/obsm/" + mem + "/" + mem_2;
						Metadata tmp_cell_m = new Metadata(p_2);
						if(reader.isGroup(p_2)) // I don't list this directory content
						{
							tmp_cell_m.nbrow = -1;
							tmp_cell_m.nbcol = -1;
						}
						else
						{
							long[] dim = reader.getDataSetInformation(p_2).getDimensions();
							tmp_cell_m.nbcol = dim[0];
							if(dim.length > 1) tmp_cell_m.nbrow = dim[1];
							else tmp_cell_m.nbrow = 1;
						}
						tmp_cell_m.on = MetaOn.CELL;
						otherMetadata.add(tmp_cell_m);
					}
				}
				else
				{
					Metadata tmp_cell_m = new Metadata("/obsm/" + mem);  // I generate the list of metadata
					long[] dim = reader.getDataSetInformation("/obsm/" + mem).getDimensions();
					if(dim.length != 2) new ErrorJSON("Type of '/obsm/" + mem + "' is not a matrix (2D)?");
					if(nbCells != dim[0]) new ErrorJSON("Dimension of '/obsm/" + mem + "' should be " + nbCells + " x " + dim[1] + " but it's " + dim[0] + " x " + dim[1]);
					tmp_cell_m.on = MetaOn.CELL;
					tmp_cell_m.nbcol = dim[0];
					tmp_cell_m.nbrow = dim[1];
					okMetadata.add(tmp_cell_m); // I take the opportunity to generate the list of metadata :)
				}
			}
		}
		
		// Cataloging and ignoring metadata for obsp
		if(reader.object().exists("/obsp")) // Not mandatory here. I don't know if it should be mandatory?
		{
			List<String> m = reader.getGroupMembers("/obsp");
			for(String mem:m) // Everything in /obsp is ignored for parsing
			{
				if(reader.isGroup("/obsp/" + mem))
				{
					// Subdirectories are not taken into account (no recursive DFS here, I only list its content)
					List<String> m_2 = reader.getGroupMembers("/obsp/" + mem);
					for(String mem_2:m_2)
					{
						String p_2 = "/obsp/" + mem + "/" + mem_2;
						Metadata tmp_cell_m = new Metadata(p_2);
						if(reader.isGroup(p_2)) // I don't list this directory content
						{
							tmp_cell_m.nbrow = -1;
							tmp_cell_m.nbcol = -1;
						}
						else
						{
							long[] dim = reader.getDataSetInformation(p_2).getDimensions();
							tmp_cell_m.nbcol = dim[0];
							if(dim.length > 1) tmp_cell_m.nbrow = dim[1];
							else tmp_cell_m.nbrow = 1;
						}
						tmp_cell_m.on = MetaOn.CELL;
						otherMetadata.add(tmp_cell_m);
					}
				}
				else
				{
					Metadata tmp_cell_m = new Metadata("/obsp/" + mem); // I generate the list of metadata
					long[] dim = reader.getDataSetInformation("/obsp/" + mem).getDimensions();
					if(dim.length != 2) new ErrorJSON("Type of '/obsp/" + mem + "' is not a matrix (2D)?");
					tmp_cell_m.on = MetaOn.CELL;
					tmp_cell_m.nbcol = dim[1];
					tmp_cell_m.nbrow = dim[0];
					otherMetadata.add(tmp_cell_m); 
				}
			}
		}
		
		// Cataloging and ignoring metadata for /raw/obs
		if(reader.object().exists("/raw/obs")) // Not mandatory. For Backward compatibility
		{
			List<String> m = reader.getGroupMembers("/raw/obs");
			for(String mem:m) // Everything in /raw/obs is ignored for parsing
			{
				if(reader.isGroup("/raw/obs/" + mem))
				{
					// Subdirectories are not taken into account (no recursive DFS here, I only list its content)
					List<String> m_2 = reader.getGroupMembers("/raw/obs/" + mem);
					for(String mem_2:m_2)
					{
						String p_2 = "/raw/obs/" + mem + "/" + mem_2;
						Metadata tmp_cell_m = new Metadata(p_2);
						if(reader.isGroup(p_2)) // I don't list this directory content
						{
							tmp_cell_m.nbrow = -1;
							tmp_cell_m.nbcol = -1;
						}
						else
						{
							long[] dim = reader.getDataSetInformation(p_2).getDimensions();
							tmp_cell_m.nbcol = dim[0];
							if(dim.length > 1) tmp_cell_m.nbrow = dim[1];
							else tmp_cell_m.nbrow = 1;
						}
						tmp_cell_m.on = MetaOn.CELL;
						otherMetadata.add(tmp_cell_m);
					}
				}
				else
				{
					Metadata tmp_cell_m = new Metadata("/raw/obs/" + mem); // I generate the list of metadata
					long[] dim = reader.getDataSetInformation("/raw/obs/" + mem).getDimensions();
					if(dim.length != 1) new ErrorJSON("Type of '/raw/obs/" + mem + "' is not a vector (1D)?");
					tmp_cell_m.on = MetaOn.CELL;
					tmp_cell_m.nbcol = dim[0];
					tmp_cell_m.nbrow = 1;
					otherMetadata.add(tmp_cell_m); 
				}
			}
		}
		return new long[] { nbGenes, nbCells };
    }
    
    private static GroupPreparse readMatrix(String path, IHDF5Reader reader, long nbGenes, long nbCells)
    {  	
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
	    	if(!sparse_format.equals("csr") && !sparse_format.equals("csc")) new ErrorJSON("Could not read the H5ad file (Sparse Matrix NOT in CSR or CSC format. Reading type not implemented yet. Please contact support.)");
	    	if(sparse_format.equals("csr")) // CSR
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
							g.isCount = g.isCount && Utils.isInteger(data[(int)j]);
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
							g.isCount = g.isCount && Utils.isInteger(data[(int)j]);
							g.matrix[i][(int)indices[(int)j]] = data[(int)j]; // TODO if the index is too big, this fails
						}
					}
				}
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
