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
	
	public static void parse()
    {
		try
		{
			System.out.println("Parsing the LOOM file...");
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
	
	    	// Process Main Matrix
		    // Prepare the writing in the file
	    	json.loom.createEmptyMatrix(json.data.nber_genes, json.data.nber_cells);
	    	
			// Parsing Data and generating summary annotations
			json.data.is_count_table = true;
	    	
	    	// Read the file block per block according to natural storage
			int[] blockSize = loom.getChunkSizes();
			int nbTotalBlocks = (int)Math.ceil((double)json.data.nber_genes / blockSize[0]);
			System.out.println("Writing & Parsing " + nbTotalBlocks + " independent blocks...");
			for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
			{
				// Retrieve the blocks that will contain all columns (because we write gene by gene)
				float[][] subMatrix = loom.readBlock(blockSize[0], (int)json.data.nber_cells, nbBlocks);  // TODO does not work if too big array
				
				// Parsing Data and generating summary annotations
				for(int x = 0; x < subMatrix.length; x++)
				{
					for(int j = 0; j < subMatrix[0].length; j++) // cells
					{
						long i = x + nbBlocks * blockSize[0]; // Original index of genes
						
						float value = subMatrix[x][j];
						
						if(json.data.erccs != null && json.data.erccs.isERCC(i)) json.data.erccs.set(i, j, value);
						else
						{
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
				}
				
				// Writing this merged block to output
				json.loom.writeBlockInMatrix(subMatrix, nbBlocks * blockSize[0], json.data.removed);
			}
			loom.close();
			
			json.loom.resizeDataset("/matrix", Parameters.nGenes - json.data.removed.size(), Parameters.nCells); // Remove extra annotations/ERCCs that were deleted from main matrix
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
