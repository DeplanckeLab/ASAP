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
	    	if(blockSize[0] / blockSize[1] > 2) LoomFile.copyFloatMatrixByCell("/matrix", loom, json.loom, json.data); // Handling HCA uneven chunks which is uncompatible with our original gene by gene process
	    	else LoomFile.copyFloatMatrixByGene("/matrix", loom, json.loom, json.data);
	    	loom.close();
	    	
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
