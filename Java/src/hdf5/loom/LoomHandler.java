package hdf5.loom;

import java.util.ArrayList;
import java.util.List;

import bigarrays.StringArray64;
import ch.systemsx.cisd.hdf5.HDF5DataSetInformation;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import hdf.hdf5lib.exceptions.HDF5JavaException;
import json.ErrorJSON;
import json.ParsingJSON;
import json.PreparsingJSON;
import model.MapGene;
import model.Parameters;
import parsing.model.GroupPreparse;

public class LoomHandler 
{
	public static void preparse()
    {
		try
		{
			System.out.println("Loom file is detected. Preparsing the file.");
	    	ArrayList<GroupPreparse> foundDatasets = new ArrayList<GroupPreparse>();
	    	GroupPreparse g = new GroupPreparse(Parameters.fileName.substring(Parameters.fileName.lastIndexOf("/") + 1)); // There should be only one dataset in a loom file
	    	IHDF5Reader reader = HDF5Factory.openForReading(Parameters.fileName);
	    	HDF5DataSetInformation info = reader.getDataSetInformation("/matrix"); // maxDim / Chunks / Dimensions / type float32... etc..
	    	g.nbGenes = info.getDimensions()[0];
	    	g.nbCells = info.getDimensions()[1];
	    	g.matrix = reader.float32().readMatrixBlock("/matrix", 10, 10, 0, 0); // First is nb row, Second is nb col, two others are offsets I suppose?      
	    	for (int i = 0; i < g.matrix.length; i++) {
				for (int j = 0; j < g.matrix[i].length; j++) {
					if(g.matrix[i][j] != (int)g.matrix[i][j]) { g.isCount = false; break;}
				}
			}
	    	if(reader.exists("/col_attrs/CellID")) 
	    	{
	    		g.cellNames = reader.string().readArrayBlock("/col_attrs/CellID", 10, 0);
	    		for (int i = 0; i < g.cellNames.length; i++) g.cellNames[i] = g.cellNames[i].trim().replaceAll("\"", "");
	    	}
	    	if(reader.exists("/row_attrs/Gene")) 
	    	{
	    		g.geneNames = reader.string().readArrayBlock("/row_attrs/Gene", 10, 0);
	    		for (int i = 0; i < g.geneNames.length; i++) g.geneNames[i] = g.geneNames[i].trim().replaceAll("\"", "");
	    	}
	    	// List Metadata
			List<String> m = reader.getGroupMembers("/row_attrs");
			for(String mem:m) if(!reader.isGroup("/row_attrs/" + mem)) g.additionalMetadataPath.add("/row_attrs/" + mem);
			m = reader.getGroupMembers("/col_attrs");
			for(String mem:m) if(!reader.isGroup("/col_attrs/" + mem)) g.additionalMetadataPath.add("/col_attrs/" + mem);
	    	// Close Loom file
			reader.close();
	    	foundDatasets.add(g);
	    	// Write final JSON file
	    	PreparsingJSON.writeOutputJSON(foundDatasets);
		}
		catch(HDF5JavaException e)
		{
			new ErrorJSON(e.getMessage());
		}
    }
	
	public static void copyMatrixByGene(String path, LoomFile in, ParsingJSON out)
	{
		int[] blockSize = in.getChunkSizes();
		blockSize[1] = 64; // We take all columns anyways, so we can change to whatever size we want
		if(blockSize[0] < 64 && 64 % blockSize[0] == 0) blockSize[0] = 64; // If chunk < 64 and multiple of 64, then better to pass to 64. If not, then we keep the original stuff
		int nbTotalBlocks = (int)Math.ceil((double)out.data.nber_genes / blockSize[0]);
    	
    	// Initialize main matrix in Loom with 0s
    	out.loom.createEmptyFloat32MatrixDataset(path, nbTotalBlocks * blockSize[0], out.data.nber_cells, blockSize[0], blockSize[1]); // I create a bit more rows, that will be filtered later on
    	
    	// Read the original file blockSize x totalNbCells according to natural storage
		System.out.println("Writing & Parsing " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = in.readFloatBlock(path, blockSize[0], (int)out.data.nber_cells, nbBlocks, 0l);
			
			// Adapting the size to block size if last block
			if(subMatrix.length < blockSize[0])
			{
				float[][] tmpMatrix = new float[blockSize[0]][(int)out.data.nber_cells];
				for(int i = 0; i < subMatrix.length; i++) tmpMatrix[i] = subMatrix[i];
				subMatrix = tmpMatrix;
			}

			// Parsing Data and generating summary annotations
			for(int x = 0; x < subMatrix.length; x++)
			{
				long i = x + nbBlocks * blockSize[0]; // Original index of genes
    			if(i < out.data.nber_genes) // In case the block is bigger than the number of genes
    			{
					for(int j = 0; j < subMatrix[0].length; j++) // cells
					{
						float value = subMatrix[x][j];
						
						if(out.data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) out.data.is_count_table = false;
							
		    			// Handle biotype count per cell
		    			String biotype = out.data.biotypes.get(i);
		    			if(biotype.equals("protein_coding")) out.data.proteinCodingContent.set(j, out.data.proteinCodingContent.get(j) + value);
		    			else if(biotype.equals("rRNA")) out.data.ribosomalContent.set(j, out.data.ribosomalContent.get(j) + value);
		    			if(out.data.chromosomes.get(i).equals("MT")) out.data.mitochondrialContent.set(j, out.data.mitochondrialContent.get(j) + value);
		        			
		        		// Generate the sums
		    			out.data.depth.set(j, out.data.depth.get(j) + value);
		    			out.data.sum.set(i, out.data.sum.get(i) + value);
							
						// Number of zeroes / detected genes 
						if(value == 0) out.data.nber_zeros++;
						else out.data.detected_genes.set(j, out.data.detected_genes.get(j) + 1);
					}
    			}
			}
			
			// Writing this block to output
			out.loom.writeFloatBlockDataset(path, subMatrix, nbBlocks, 0);
		}
		in.close();
	}
	
	public static void copyMatrixByCell(String path, LoomFile in, ParsingJSON out)
	{
		int[] blockSize = in.getChunkSizes();
		blockSize[0] = 64; // We take all rows anyways, so we can change to whatever size we want
		if(blockSize[1] < 64 && 64 % blockSize[1] == 0) blockSize[1] = 64; // If chunk < 64 and multiple of 64, then better to pass to 64. If not, then we keep the original stuff
		int nbTotalBlocks = (int)Math.ceil((double)out.data.nber_cells / blockSize[1]);
    		
    	// Initialize main matrix in Loom with 0s
    	out.loom.createEmptyFloat32MatrixDataset(path, out.data.nber_genes, nbTotalBlocks * blockSize[1], blockSize[0], blockSize[1]); // I create a bit more cols, that will be filtered later on
    	
    	// Read the original file blockSize x totalNbGenes
		System.out.println("Writing & Parsing " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = in.readFloatBlock(path, (int)out.data.nber_genes, blockSize[1], 0l, nbBlocks);

			// Adapting the size to block size if last block
			if(subMatrix[0].length < blockSize[1])
			{
				float[][] tmpMatrix = new float[(int)out.data.nber_genes][blockSize[1]];
				for(int i = 0; i < subMatrix.length; i++) for(int j = 0; j < subMatrix[i].length; j++) tmpMatrix[i][j] = subMatrix[i][j];
				subMatrix = tmpMatrix;
			}

			// Parsing Data and generating summary annotations
			for(int i = 0; i < subMatrix.length; i++) // Original index of gene
			{
				for(int y = 0; y < subMatrix[0].length; y++) // cells
				{
					long j = y + nbBlocks * blockSize[1]; // Original index of cells
					
	    			if(j < out.data.nber_cells) // In case the block is bigger than the number of cells
	    			{
						float value = subMatrix[i][y];
						
						if(out.data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) out.data.is_count_table = false;
							
		    			// Handle biotype count per cell
		    			String biotype = out.data.biotypes.get(i);
		    			if(biotype.equals("protein_coding")) out.data.proteinCodingContent.set(j, out.data.proteinCodingContent.get(j) + value);
		    			else if(biotype.equals("rRNA")) out.data.ribosomalContent.set(j, out.data.ribosomalContent.get(j) + value);
		    			if(out.data.chromosomes.get(i).equals("MT")) out.data.mitochondrialContent.set(j, out.data.mitochondrialContent.get(j) + value);
		        			
		        		// Generate the sums
		    			out.data.depth.set(j, out.data.depth.get(j) + value);
		    			out.data.sum.set(i, out.data.sum.get(i) + value);
							
						// Number of zeroes / detected genes 
						if(value == 0) out.data.nber_zeros++;
						else out.data.detected_genes.set(j, out.data.detected_genes.get(j) + 1);
					}
    			}
			}
			
			// Writing this block to output
			out.loom.writeFloatBlockDataset(path, subMatrix, 0, nbBlocks);
		}
		in.close();
	}
	
	public static void parse()
    {
		try
		{
			ParsingJSON json = new ParsingJSON();
			LoomFile loom = new LoomFile("r", Parameters.fileName);
			
	    	// Process cells
	    	StringArray64 cellNames = loom.getCellNames();
	    	if(cellNames != null) 
	    	{
	    		if(json.data.cell_names.size() != cellNames.size()) new ErrorJSON("Nb cells in args does not match the Loom file");
	    		json.data.cell_names = cellNames;
	    		for (int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, json.data.cell_names.get(i).replaceAll("'|\"", ""));
	    	}
	    	else for(int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, "Cell_"+(i+1)); // Create fake names if none exist
	    	
	    	// Process Gene names
	    	if(loom.exists("/row_attrs/Accession")) json.data.ens_names = loom.readStringArray("/row_attrs/Accession"); // EnsemblIDs
	    	else json.data.ens_names = new StringArray64(json.data.nber_genes); // After checking the database, this is filled appropriately
	    	if(loom.exists("/row_attrs/Gene")) json.data.gene_names = loom.readStringArray("/row_attrs/Gene"); // HGNC names
	    	else json.data.gene_names = new StringArray64(json.data.nber_genes); // After checking the database, this is filled appropriately
	    	MapGene.parseGenes(json.data);
	
	    	// Original file block sizes is kept (easier to perform the copy)
	    	json.data.is_count_table = true;
	    	int[] blockSize = loom.getChunkSizes();
	    	if(blockSize[0] / blockSize[1] > 2) copyMatrixByCell("/matrix", loom, json); // Handling HCA uneven chunks which is uncompatible with our original gene by gene process
	    	else copyMatrixByGene("/matrix", loom, json);
	    				
	    	json.loom.resizeDataset("/matrix", (int)json.data.nber_genes, (int)json.data.nber_cells); // Cause writing fixed-size blocks can extend the matrix size with 0
			LoomFile.fillLoomFile(json.loom, json.data);
			json.writeOutputJSON();
			json.loom.close();
		}
		catch(HDF5JavaException e)
		{
			new ErrorJSON(e.getMessage());
		}
    }
	
	public static void parse2()
    {
		try
		{
			ParsingJSON json = new ParsingJSON();
			LoomFile loom = new LoomFile("r", Parameters.fileName);
			
	    	// Process cells
	    	StringArray64 cellNames = loom.getCellNames();
	    	if(cellNames != null) 
	    	{
	    		if(json.data.cell_names.size() != cellNames.size()) new ErrorJSON("Nb cells in args does not match the Loom file");
	    		json.data.cell_names = cellNames;
	    		for (int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, json.data.cell_names.get(i).replaceAll("'|\"", ""));
	    	}
	    	else for(int i = 0; i < json.data.cell_names.size(); i++) json.data.cell_names.set(i, "Cell_"+(i+1)); // Create fake names if none exist
	    	
	    	// Process Gene names
	    	if(loom.exists("/row_attrs/Accession")) json.data.ens_names = loom.readStringArray("/row_attrs/Accession"); // EnsemblIDs
	    	else json.data.ens_names = new StringArray64(json.data.nber_genes); // After checking the database, this is filled appropriately
	    	if(loom.exists("/row_attrs/Gene")) json.data.gene_names = loom.readStringArray("/row_attrs/Gene"); // HGNC names
	    	else json.data.gene_names = new StringArray64(json.data.nber_genes); // After checking the database, this is filled appropriately
	    	MapGene.parseGenes(json.data);
	
	    	// Original file block sizes is kept (easier to perform the copy)
	    	json.data.is_count_table = true;
	    	int[] blockSize = loom.getChunkSizes();
	    	
	    	// Initialize main matrix in Loom with 0s
	    	json.loom.createEmptyFloat32MatrixDataset("/matrix", json.data.nber_genes, json.data.nber_cells, blockSize[0], blockSize[1]);
	    	
	    	// Read the original file blockSize x totalNbCells according to natural storage
			int nbTotalBlocks = (int)Math.ceil((double)json.data.nber_genes / blockSize[0]);
			System.out.println("Writing & Parsing " + nbTotalBlocks + " independent blocks...");
			for(int nbBlocks = 0; nbBlocks < nbTotalBlocks - 1; nbBlocks++)
			{
				// Retrieve the blocks that will contain all columns (because we will write all cells at once)
				float[][] subMatrix = loom.readFloatBlock("/matrix", blockSize[0], (int)json.data.nber_cells, nbBlocks, 0l);

				// Parsing Data and generating summary annotations
				for(int x = 0; x < subMatrix.length; x++)
				{
					for(int j = 0; j < subMatrix[0].length; j++) // cells
					{
						long i = x + nbBlocks * blockSize[0]; // Original index of genes
						
						float value = subMatrix[x][j];
						
						if(json.data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) json.data.is_count_table = false;
							
		    			// Handle biotype count per cell
		    			String biotype = json.data.biotypes.get(i);
		    			if(biotype.equals("protein_coding")) json.data.proteinCodingContent.set(j, json.data.proteinCodingContent.get(j) + value);
		    			else if(biotype.equals("rRNA")) json.data.ribosomalContent.set(j, json.data.ribosomalContent.get(j) + value);
		    			if(json.data.chromosomes.get(i).equals("MT")) json.data.mitochondrialContent.set(j, json.data.mitochondrialContent.get(j) + value);
		        			
		        		// Generate the sums
						json.data.depth.set(j, json.data.depth.get(j) + value);
						json.data.sum.set(i, json.data.sum.get(i) + value);
							
						// Number of zeroes / detected genes 
						if(value == 0) json.data.nber_zeros++;
						else json.data.detected_genes.set(j, json.data.detected_genes.get(j) + 1);
					}
				}
				
				// Writing this block to output
				json.loom.writeFloatBlockDataset("/matrix", subMatrix, nbBlocks, 0);
			}
			
			// The last block may be incorrectly written if nbRows is different than chunk size. we treat it separately
			// Reading last block
			float[][] subMatrix = loom.readFloatBlock("/matrix", blockSize[0], (int)json.data.nber_cells, nbTotalBlocks - 1, 0l);
			if(blockSize[0] > subMatrix.length) // if sizeX different than chunk sizeX (genes)
			{
				float[][] dataResized = new float[blockSize[0]][(int)json.data.nber_cells];
				for(int i = 0; i < subMatrix.length; i++) for(int j = 0; j < subMatrix[i].length; j++) dataResized[i][j] = subMatrix[i][j];
				subMatrix = dataResized;
			}
			// Writing this block to output
			json.loom.writeFloatBlockDataset("/matrix", subMatrix, nbTotalBlocks - 1, 0);
			loom.close();
			
			json.loom.resizeDataset("/matrix", (int)json.data.nber_genes, (int)json.data.nber_cells);
			LoomFile.fillLoomFile(json.loom, json.data);
			json.writeOutputJSON();
			json.loom.close();
		}
		catch(HDF5JavaException e)
		{
			new ErrorJSON(e.getMessage());
		}
    }
}
