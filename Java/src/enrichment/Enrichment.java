package enrichment;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import bigarrays.StringArray64;
import db.DBManager;
import enrichment.model.EnrichmentJSON;
import enrichment.model.ResultSet;
import hdf5.loom.LoomFile;
import json.ErrorJSON;
import model.GeneSet;
import model.Parameters;
import tools.FisherExactTest_2X2;
import tools.Utils;
import tools.FisherExactTest_2X2.Alternative;

public class Enrichment 
{
	/**
	 * Enrichment main function (computes the Fisher's Exact Test p-values for a given list of genesets)
	 */	
	public static void runEnrichment()
	{
		// Fetching DB
		DBManager.connect();
		ArrayList<GeneSet> genesets = DBManager.getGeneSets(Parameters.geneset_id);
		long organism_id = DBManager.getOrganismFromGeneSets(Parameters.geneset_id);
		HashMap<String, Long> genes = DBManager.getGenesInDBByID(organism_id);
		//System.out.println(genes.size() + " genes were fetched in DB");
		DBManager.disconnect();
		//System.out.println(genesets.size() + " genesets were found in ASAP DB for id == " + Parameters.geneset_id);
		if(genesets.size() == 0) new ErrorJSON("No genesets in DB corresponding to this id", Parameters.JSONFileName);
		HashSet<Long> genesInGenesets = new HashSet<Long>();
		for(GeneSet g:genesets) for(long gene:g.content) genesInGenesets.add(gene);
		//System.out.println(genesInGenesets.size() + " unique genes in selected geneset");
		
    	// Open Loom file in read-only
    	LoomFile loom = new LoomFile("r", Parameters.loomFile);
    	StringArray64 ens_ids = loom.readStringArray("/row_attrs/Accession");
		if(ens_ids == null) new ErrorJSON("No Ensembl Ids in Loom file.", Parameters.JSONFileName);
    	//System.out.println(ens_ids.size() + " genes found in the Loom file");
    	loom.close();
    	
    	HashMap<String, Long> dbIDByEnsemblId = new HashMap<String, Long>();
		for(String ens:ens_ids) 
		{
			Long id = genes.get(ens);
			if(id != null && genesInGenesets.contains(id)) dbIDByEnsemblId.put(ens, id);
		}
    	
		//System.out.println(dbIDByEnsemblId.size() + " unique genes left after overlap with background + Loom file (by EnsemblIds)");
		if(dbIDByEnsemblId.size() == 0) new ErrorJSON("No gene in the Loom file matches Ensembl Ids in the database. Nothing to enrich.", Parameters.JSONFileName);
		
		// Getting JSON containing the list(s) of genes to enrich 
		EnrichmentJSON json = EnrichmentJSON.parseJSON(Parameters.fileName);
		
		// Running enrichment
    	HashMap<String,  ArrayList<ResultSet>> listRes = new HashMap<String,  ArrayList<ResultSet>>();
		if(json.down != null) listRes.put("down", enrich(loomIndexesToDBIds(json.down, dbIDByEnsemblId, ens_ids), genesets, dbIDByEnsemblId));
		if(json.up != null) listRes.put("up", enrich(loomIndexesToDBIds(json.up, dbIDByEnsemblId, ens_ids), genesets, dbIDByEnsemblId));
		if(json.indexes_match != null) listRes.put("indexes_match", enrich(loomIndexesToDBIds(json.indexes_match, dbIDByEnsemblId, ens_ids), genesets, dbIDByEnsemblId));
    	
		// Prepare output String
		StringBuilder sb = new StringBuilder();
    	sb.append("{").append("\"time_idle\":").append(Parameters.idleTime);
    	sb.append(",\"headers\":[\"name\",\"description\",\"p-value\",\"").append(Parameters.adjMethod).append("\",\"effect size\",\"size geneset\",\"overlap w/ genes\"]");
		
		// Go through the results and create the output String
		for(String key:listRes.keySet())
		{
			sb.append(",\"").append(key).append("\":[");
			ArrayList<ResultSet> res = listRes.get(key);
			
			// First loop to get all pvalues and adjust
			double[] pvalues = new double[res.size()];
			for(int i = 0; i < pvalues.length; i++) pvalues[i] = res.get(i).p_value;
			double[] adj_p_value = Utils.p_adjust(pvalues, Parameters.adjMethod);
			int[] sortedIndexes = Utils.order(pvalues, false);
			
			// Second loop to generate the JSON file
			String prefix = "";
			for(int i:sortedIndexes)
			{
				ResultSet enrichRes = res.get(i);
				sb.append(prefix).append("[\"").append(enrichRes.name).append("\"");
				sb.append(",\"").append(enrichRes.description).append("\"");
				sb.append(",").append(enrichRes.p_value);
				sb.append(",").append(adj_p_value[i]);
				sb.append(",").append((enrichRes.effect_size == Double.MAX_VALUE)?"\"Inf\"":enrichRes.effect_size);
				sb.append(",").append(enrichRes.size);
				sb.append(",").append(enrichRes.overlap);
				sb.append("]");
				prefix = ",";
			} 
			sb.append("]");
		}
		sb.append("}");
		
    	// Write output.json
    	Utils.writeJSON(sb, Parameters.JSONFileName);
	}
	
	public static void runMarkerEnrichment()
	{
		Parameters.minGenesInPathway = 5;
		Parameters.maxGenesInPathway = 500;
		
		// Get list of output files from FindMarker's step
		ArrayList<File> files = Utils.listMarkerFiles(Parameters.fileName);
		if(files.size() == 0) new ErrorJSON("No marker files found!");
		
		// Fetching DB
		ArrayList<GeneSet> genesets = new ArrayList<GeneSet>();
		long organism_id = -1;
		DBManager.connect();
		for(long geneset_id:Parameters.geneset_ids)
		{
			genesets.addAll(DBManager.getGeneSets(geneset_id));
			long id = DBManager.getOrganismFromGeneSets(geneset_id);
			if(organism_id == -1) organism_id = id;
			if(organism_id != id) new ErrorJSON("Cannot perform Enrichment across multiple species");
		}
		HashMap<String, Long> genes = DBManager.getGenesInDBByID(organism_id);
		DBManager.disconnect();
	
		// Create unique gene hash
		if(genesets.size() == 0) new ErrorJSON("No genesets found in DB corresponding to this/these id(s)");
		HashSet<Long> genesInGenesets = new HashSet<Long>();
		for(GeneSet g:genesets) for(long gene:g.content) genesInGenesets.add(gene);

		// Perform enrichment for each .tsv file
		HashMap<String,  ArrayList<ResultSet>> listRes = new HashMap<String,  ArrayList<ResultSet>>();
		try
		{
			for(File tsv:files)
			{
				// Prepare lists
				HashMap<String, Long> dbIDByEnsemblId = new HashMap<String, Long>(); // background genes
				HashSet<Long> upRegulatedGenes = new HashSet<Long>(); // up-regulated genes
				
		    	// Read file
				BufferedReader br = new BufferedReader(new FileReader(tsv));
				br.readLine(); // Header: stable_id	original_gene	gene	ensembl	logfc	pval	fdr	avg_1	avg_not_1
				String line = br.readLine();
				try
				{
					while(line != null)
					{
						String[] tokens = line.split("\t");
						float logfc = Float.parseFloat(tokens[4]);
						float fdr = Float.parseFloat(tokens[6]);
						String ensembl = tokens[3];
						Long id = genes.get(ensembl);
						if(id != null && genesInGenesets.contains(id)) 
						{
							dbIDByEnsemblId.put(ensembl, id);
							if(logfc >= 1 && fdr <= 0.05) upRegulatedGenes.add(id);// Abritrary cutoff for FET. Not optimal. Should use a ranked test (like GSEA). Here I keep only UP-regulated genes
						}
						line = br.readLine();
					}
				}
				catch(NumberFormatException nfe)
				{
					new ErrorJSON(nfe.getMessage());
				}
				br.close();

				// Running enrichment
		    	listRes.put(tsv.getAbsolutePath(), enrich(upRegulatedGenes, genesets, dbIDByEnsemblId));
			}
		}
		catch(IOException ioe)
		{
			new ErrorJSON(ioe.getMessage());
		}

		// Prepare output String {TSV MODE}
		for(String key:listRes.keySet())
		{
			StringBuilder sb = new StringBuilder();
	    	sb.append("name\tdescription\tp-value\t").append(Parameters.adjMethod).append("\teffect_size\tsize_geneset\toverlap_w_genes\n");
			ArrayList<ResultSet> res = listRes.get(key);
			
			// First loop to get all pvalues and adjust
			double[] pvalues = new double[res.size()];
			for(int i = 0; i < pvalues.length; i++) pvalues[i] = res.get(i).p_value;
			double[] adj_p_value = Utils.p_adjust(pvalues, Parameters.adjMethod);
			int[] sortedIndexes = Utils.order(pvalues, false);
			
			// Second loop to generate the TSV content
			for(int i:sortedIndexes)
			{
				ResultSet enrichRes = res.get(i);
				sb.append(enrichRes.name).append("\t");
				sb.append(enrichRes.description).append("\t");
				sb.append(enrichRes.p_value).append("\t");
				sb.append(adj_p_value[i]).append("\t");
				sb.append((enrichRes.effect_size == Double.MAX_VALUE)?"\"Inf\"":enrichRes.effect_size).append("\t");
				sb.append(enrichRes.size).append("\t");
				sb.append(enrichRes.overlap).append("\n");;
			}
			
			// Write in TSV file
			try
			{
				BufferedWriter bw = new BufferedWriter(new FileWriter(key.replaceAll(".tsv", ".enrichment.tsv")));
				bw.write(sb.toString());
				bw.close();
			}
			catch(IOException ioe)
			{
				new ErrorJSON(ioe.getMessage());
			}
		}
		
		// Prepare output String {JSON MODE}
		/*StringBuilder sb = new StringBuilder();
    	sb.append("{").append("\"time_idle\":").append(Parameters.idleTime);
    	sb.append(",\"headers\":[\"name\",\"description\",\"p-value\",\"").append(Parameters.adjMethod).append("\",\"effect size\",\"size geneset\",\"overlap w/ genes\"]");
		
		// Go through the results and create the output String
		for(String key:listRes.keySet())
		{
			sb.append(",\"").append(key).append("\":[");
			ArrayList<ResultSet> res = listRes.get(key);
			
			// First loop to get all pvalues and adjust
			double[] pvalues = new double[res.size()];
			for(int i = 0; i < pvalues.length; i++) pvalues[i] = res.get(i).p_value;
			double[] adj_p_value = Utils.p_adjust(pvalues, Parameters.adjMethod);
			int[] sortedIndexes = Utils.order(pvalues, false);
			
			// Second loop to generate the JSON file
			String prefix = "";
			for(int i:sortedIndexes)
			{
				ResultSet enrichRes = res.get(i);
				sb.append(prefix).append("[\"").append(enrichRes.name).append("\"");
				sb.append(",\"").append(enrichRes.description).append("\"");
				sb.append(",").append(enrichRes.p_value);
				sb.append(",").append(adj_p_value[i]);
				sb.append(",").append((enrichRes.effect_size == Double.MAX_VALUE)?"\"Inf\"":enrichRes.effect_size);
				sb.append(",").append(enrichRes.size);
				sb.append(",").append(enrichRes.overlap);
				sb.append("]");
				prefix = ",";
			} 
			sb.append("]");
		}
		sb.append("}");
    	
    	// Write output.json
    	Utils.writeJSON(sb);*/
	}

	private static HashSet<Long> loomIndexesToDBIds(long[] loomIdx, HashMap<String, Long> dbIDByEnsemblId, StringArray64 ens_ids)
	{
		HashSet<Long> res = new HashSet<Long>();
		for(long idx:loomIdx)
		{
			String ensemblId = ens_ids.get(idx);
			if(ensemblId != null)
			{
				Long id = dbIDByEnsemblId.get(ensemblId);
				if(id != null) res.add(id);
			}
		}
		return res;
	}
	
	private static ArrayList<ResultSet> enrich(HashSet<Long> indexesInLoom, ArrayList<GeneSet> genesets, HashMap<String, Long> backgroundMap)
	{
		HashSet<Long> background = new HashSet<Long>();
		for(String key:backgroundMap.keySet()) background.add(backgroundMap.get(key));
		int indexesInLoomAndBackground = 0;
		for(Long idx:indexesInLoom) if(background.contains(idx)) indexesInLoomAndBackground++;
		//System.out.println("Performing Functional Enrichment on " + indexesInLoom.size() + " genes");
		ArrayList<ResultSet> res = new ArrayList<ResultSet>();
		for(GeneSet g:genesets)
		{
			ResultSet gRes = new ResultSet();
			switch(Parameters.enrichModel) 
			{
				case FET:
					gRes.name = g.identifier;
					gRes.description = g.name;
					int overlap = 0;
					int genesetSize = 0;
					for(long gene:g.content) 
					{
						if(background.contains(gene))
						{
							genesetSize++;
							if(indexesInLoom.contains(gene)) overlap++;
						}
					}
					if(genesetSize >= Parameters.minGenesInPathway && genesetSize <= Parameters.maxGenesInPathway) // Process only these
					{
						double[] resFet = doFExactTest(overlap, genesetSize - overlap, indexesInLoomAndBackground - overlap, background.size() - genesetSize - indexesInLoomAndBackground + overlap);
						gRes.overlap = overlap;
						gRes.p_value = resFet[0];
						gRes.effect_size = resFet[1];
						gRes.size = genesetSize;
						res.add(gRes);
					}
					break;
				default:
					new ErrorJSON("This model is not yet implemented", Parameters.JSONFileName);
			}
			
		}
		return res;
	}
	
    private static double[] doFExactTest(int a, int b, int c, int d) // res[0] = pvalue, res[1] = OR
    {
        double p = FisherExactTest_2X2.getSingleton().fisher(a, b, c, d, Alternative.GREATER); // Homemade function (should be similar to R output)
        if(b == 0 || c == 0) return new double[]{p, Double.MAX_VALUE}; // Special case, OR = Inf (but p can be computed)
        if(d == 0) return new double[]{1, 0}; // Special case, p = 1 & OR = 0
        double OR = (double)(a * d) / (b * c);
        return new double[]{p, OR};
    }
	
	/*private static ResultSet applyGSEA(HashMap<String, ArrayList<String>> pathways, HashMap<String, Double> exprs_1, HashMap<String, Double> exprs_2, int minSize, int maxSize)
	{	
		ResultSet res = new ResultSet();
		// Compute the FC Vector
		float[] fc = new float[exprs_1.size()];
        String[] names = new String[exprs_1.size()];
        int i = 0;
		for(String gene:exprs_1.keySet()) 
        {
			names[i] = gene;
			fc[i] = (float)Math.abs(exprs_2.get(gene) - exprs_1.get(gene));
			i++;
		}
        RankedList rl = RankedListGenerators.createBySorting("RankedGenes", names, fc, SortMode.REAL, Order.DESCENDING);
        // Generate the pathway set
		ArrayList<GeneSet> gsets = new ArrayList<>();
        for(String pathway:pathways.keySet()) 
		{	
			ArrayList<String> genes = pathways.get(pathway);
			if(genes.size() >= minSize && genes.size() <= maxSize)
			{
				gsets.add(new FSet(pathway, genes.toArray(new String[0])));
			}
		}
        // Run ssGSEA
        EnrichmentDb edb = new KSTests(nullStream).executeGsea(rl, gsets.toArray(new GeneSet[0]), Parameters.nbRepeat, rsg, null, new GeneSetScoringTableReqdParam().createGeneSetCohortGenerator(true));
        // Update the ResultSet
        PValueCalculator pvc = new PValueCalculatorImpls.GseaImpl("meandiv");
        EnrichmentResult[] results = pvc.calcNPValuesAndFDR(edb.getResults());
        EnrichmentDb edb_with_fdr = edb.cloneDeep(results);
        for (int j = 0; j < edb_with_fdr.getNumResults(); j++) 
        {
            EnrichmentResult result = edb_with_fdr.getResult(j);
            String pathway = result.getGeneSetName();
            double pvalue = (double)result.getScore().getNP();
            res.p_value.put(pathway, pvalue);
    		res.ES.put(pathway, (double)result.getScore().getES());
    		res.NES.put(pathway, (double)result.getScore().getNES());
    		res.Z_score.put(pathway, Stats.getZScore(pvalue));
    		res.FDR.put(pathway, (double)result.getScore().getFDR());
    		res.FWER.put(pathway, (double)result.getScore().getFWER());
        }
		return res;
	}*/
}






