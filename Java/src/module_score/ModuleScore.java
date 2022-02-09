package module_score;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import bigarrays.FloatArray64;
import bigarrays.StringArray64;
import db.DBManager;
import hdf5.loom.LoomFile;
import json.ErrorJSON;
import model.GeneSet;
import model.MetaOn;
import model.Metadata;
import model.Metatype;
import model.Parameters;
import tools.Utils;

public class ModuleScore 
{
	/**
	 * ModuleScore main function (computes the score for each cell/sample for a given geneset)
	 */	
	public static void runModuleScore()
	{
		ArrayList<Long> selection = new ArrayList<Long>(); // Indexes in Loom file
		if(Parameters.geneset_id != -1)
		{
			// Fetching gene set in DB
			DBManager.connect();
			GeneSet geneset = DBManager.getUniqueGeneSet(Parameters.geneset_id);
			long organism_id = DBManager.getOrganismFromGeneSets(geneset.gene_set_id);
			HashMap<String, Long> genes = DBManager.getGenesInDBByID(organism_id);
			DBManager.disconnect();

			// Open Loom file in read-only
			LoomFile loom = new LoomFile("r", Parameters.loomFile);
	    	StringArray64 ens_ids = loom.readStringArray("/row_attrs/Accession");
			if(ens_ids == null) new ErrorJSON("No Ensembl Ids in Loom file.", Parameters.JSONFileName);
	    	loom.close();

	    	// Check which genes of Loom match genes of Database
			for(long i = 0; i < ens_ids.size(); i++) 
			{
				String ens = ens_ids.get(i);
				Long id = genes.get(ens);
				if(id != null && geneset.content.contains(id)) selection.add(i);
			}
		}
		else
		{
			if(Parameters.metaName.contains("col_attrs")) new ErrorJSON("The metadata should be on the genes, not the cells");
			// Fetching gene set in metadata
			LoomFile loom = new LoomFile("r", Parameters.loomFile);
			selection = loom.getIndexesWhereValueIs(Parameters.metaName, Parameters.selection);
			loom.close();
		}
		
		// Run scoring
		float[] scores = null;
    	if(selection.size() != 0)
    	{
    		switch(Parameters.moduleScoreModel)
    		{
    		case PCA:
    			scores = scorePCA(selection);
    			break;
    		case Seurat:
    			scores = scoreSeurat(selection, Parameters.nBackgroundGenes, Parameters.nBins); // default is 100  24
    			break;
    		default:
    			new ErrorJSON("This method is not implemented");
    		}
    	}
		
		// Check whether I need to add metadata to Loom file
    	StringBuilder meta = new StringBuilder(); // What will be in the JSON
    	if(selection.size() != 0) 
    	{
    		if(Parameters.oAnnot == null) meta.append(",\"scores\":").append(Utils.toString(scores)); // To be added in the JSON
    		else // Added in the Loom file
    		{  			
    			Metadata m = new Metadata(); // Create metadata (just for the JSON)
    			m.on = MetaOn.CELL;
    			m.path = Parameters.oAnnot;
    			m.type = Metatype.NUMERIC;
    			m.nbcol = scores.length;
    			m.nbrow = 1;
    			
    			LoomFile loom = new LoomFile("r+", Parameters.loomFile);
    			loom.writeMetadata(Parameters.oAnnot, scores);
    			m.size = loom.getSizeInBytes(Parameters.oAnnot); // Compute the size
    			loom.close();
    			
    			meta.append(",\"metadata\":[");
    			m.addMeta(meta, false, null, -1); // Add to String
    			meta.append("]");
    		}
    	}
    	
    	// Create output String
		StringBuilder sb = new StringBuilder();
    	sb.append("{").append("\"time_idle\":").append(Parameters.idleTime).append(",\"overlap\":").append(selection.size()).append(meta).append("}");
		
    	// Write output.json
    	Utils.writeJSON(sb, Parameters.JSONFileName);
	}
	
	private static float[] scorePCA(ArrayList<Long> selection)
	{	
		// TODO
		return null;
	}
	
	// In Seurat package, this method is ran on the normalized data if it exists. If not, on count data. Not on scaled data. 
	private static float[] scoreSeurat(ArrayList<Long> selection, int ctrl, int nbin)
	{	
		// Recuperate subMatrix from Loom file
		LoomFile loom = new LoomFile("r", Parameters.loomFile);
		FloatArray64 sum = loom.readFloatArray("/row_attrs/_Sum");
		loom.close();
		
		float[] data_sum = sum.toArray();
		if(data_sum == null) new ErrorJSON("'Quantiles' function is not implemented yet for big arrays");

		// Find quantiles
		int[] quantiles = Utils.quantiles(data_sum, nbin);
		
		// Index them for faster access
		HashMap<Integer, List<Integer>> indexQuantiles = new HashMap<Integer, List<Integer>>(nbin);
		for (int i = 0; i < quantiles.length; i++) 
		{
			Integer key = quantiles[i];
			List<Integer> listIndexes = indexQuantiles.get(key);
			if(listIndexes == null) listIndexes = new ArrayList<Integer>();
			listIndexes.add(i);
			indexQuantiles.put(key, listIndexes);
		}
		
		// Now, build the ctrl by picking randomly "ctrl" genes in each corresponding bin
		HashSet<Integer> ctrl_use = new HashSet<Integer>();
		// Go through each selected index
		for(Long sel:selection)
		{
			int quant = quantiles[sel.intValue()];
			List<Integer> otherIndexesWithSameQuantile = indexQuantiles.get(quant);
			int[] randomSample = Utils.sample(otherIndexesWithSameQuantile, ctrl);
			for(int s:randomSample) ctrl_use.add(s);
		}
		
		// Formatting in correct array type for readRows function
		ArrayList<Long> ctrl_use_format = new ArrayList<Long>();
		for(Integer c:ctrl_use) ctrl_use_format.add(c.longValue());
		
		// Recuperate subMatrices from Loom file
		loom = new LoomFile("r", Parameters.loomFile);
		float[][] features_subset = loom.readRows(selection, Parameters.metaForComputation);
		float[][] ctrl_subset = loom.readRows(ctrl_use_format, Parameters.metaForComputation);	
		loom.close();
		
		// TODO ColMeans could be done directly in readRows (not storing the submatrices)
		float[] features_scores = Utils.colMeans(features_subset);
		float[] ctrl_scores = Utils.colMeans(ctrl_subset);
		
		float[] res = new float[features_scores.length];
		for (int i = 0; i < res.length; i++) res[i] = features_scores[i] - ctrl_scores[i];
		
		return res;
	}
}






