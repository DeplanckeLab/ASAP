package enrichment;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import bigarrays.StringArray64;
import db.DBManager;
import enrichment.model.EnrichmentJSON;
import enrichment.model.ResultSet;
import hdf5.loom.LoomFile;
import jsc.contingencytables.ContingencyTable2x2;
import jsc.contingencytables.FishersExactTest;
import json.ErrorJSON;
import model.GeneSet;
import model.Parameters;
import tools.Utils;

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
		int organism_id = DBManager.getOrganismFromGeneSets(Parameters.geneset_id);
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
    	if(a == 0) return new double[]{1, 0};
    	if(b == 0 || c == 0) return new double[]{0, Double.MAX_VALUE};
    	if(d == 0) return new double[]{1, 0}; // I do this after, in case c & d are both equal to 0
        ContingencyTable2x2 table = new ContingencyTable2x2(a, b, c, d);
        FishersExactTest test = new FishersExactTest(table);
        double twotailedP = test.getOppositeTailProb() + test.getOneTailedSP();
        double OR = (double)(a * d) / (b * c);
        return new double[]{twotailedP, OR};
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






