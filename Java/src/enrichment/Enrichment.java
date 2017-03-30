package enrichment;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.Writer;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.stream.JsonReader;

import jsc.contingencytables.ContingencyTable2x2;
import jsc.contingencytables.FishersExactTest;
import model.Parameters;
import tools.Utils;

public class Enrichment 
{
	// Pathway data
	private static HashMap<String, Pathway> data_pathways = new HashMap<String, Pathway>();
	
	// Annotation Data
	//private static HashMap<String, Gene> ensToGene = new HashMap<String, Gene>();
	private static HashMap<String, Boolean> backgroundGenes = new HashMap<String, Boolean>(); // ID = Gene Name
	
	// Gene list to enrich
	private static HashSet<String> genesToEnrich = new HashSet<String>();
	
	// Potential warning message
	private static String warningMess = null;
	
	public static void readFiles()
	{
		long t1 = System.currentTimeMillis();
		
		// Load Background Genes
		if(!Parameters.isSilent) System.out.println("Reading the Background file...");
		loadBackgroundGenes(Parameters.backgroundFile);
		if(!Parameters.isSilent) System.out.println("Background file provided: "+backgroundGenes.size()+" IDs are found.");
		
		// Load Gene-Pathway mapping
		if(!Parameters.isSilent) System.out.println("Reading the Gene-Pathway mapping file...");
		loadDataPathways(Parameters.pathwayFile);
		if(!Parameters.isSilent) System.out.println("Pathway Map provided: "+data_pathways.keySet().size()+" pathways.");
		String[] copy = backgroundGenes.keySet().toArray(new String[backgroundGenes.size()]);
		for(String gene:copy) // Remove extra genes and put back false to extra genes
		{
			if(backgroundGenes.get(gene) == false) backgroundGenes.remove(gene);
			else backgroundGenes.put(gene, false);
		}
		if(!Parameters.isSilent) System.out.println("Background file provided: "+backgroundGenes.size()+" IDs are remaining.");
		
		// Load Genes to Enrich
		if(!Parameters.isSilent) System.out.println("Reading the Gene list to enrich [JSON]...");
		loadGeneListJSON(Parameters.listGenesFile);
		if(!Parameters.isSilent) System.out.println("Gene list to enrich [JSON] provided: "+genesToEnrich.size()+" genes.");
		warningMess = "After filtering out genes not present in pathway file: "+genesToEnrich.size()+" genes remain to be enriched.";
		
		if(!Parameters.isSilent) System.out.println("Loading file time: "+Utils.toReadableTime(System.currentTimeMillis() - t1));
	}
	
	/**
	 * Enrichment main function (computes the p-values for a given list of pathways)
	 * @param toto This is a toto!
	 */
	public static void runEnrichment()
	{
		long t1 = System.currentTimeMillis();
		// Initialization
		File folder = new File(Parameters.outputFolder); // This should be a folder specified, not a single file
		if(folder.exists() && !folder.isDirectory())
		{
			System.err.println("The Output folder specified is not a folder.\nStopping program...");
			System.exit(-1);
		}
		if(!Parameters.isSilent) System.out.println("\n" + Parameters.model.toString() + " model is used.\n");
		
		// Running Enrichment for each pathway
		ResultSet res = new ResultSet();
		res.init(data_pathways.size());
		int k = 0;
		if(!Parameters.isSilent) System.out.println("Computing the scores for each pathway...");
		switch(Parameters.model)
		{
		case FET:
			for (String path:data_pathways.keySet()) 
			{
				Pathway p = data_pathways.get(path);
				res.pathways[k] = p.id;
				res.descriptions[k] = p.description;
				res.urls[k] = p.url;
				int overlap = 0;
				for(String gene:p.listGenes)
				{
					for(String g:genesToEnrich)
					{
						if(gene.equals(g)) overlap++;
					}
				}
				double[] resFet = doFExactTest(overlap, p.listGenes.size() - overlap, genesToEnrich.size() - overlap, backgroundGenes.size() - p.listGenes.size() - genesToEnrich.size() + overlap);
				res.p_value[k] = resFet[0];
				res.OR[k] = resFet[1];
				k++;
			}
			break;
		default:
			System.err.println("This model is not yet implemented.");
			System.exit(-1);
			break;
		}
		// Adjust P-values (all of them need to be computed to do so)
		res.adj_p_value = p_adjust(res.p_value, Parameters.adjMethod);
		res.warning = Enrichment.warningMess;
		// Filtering Results
		int filtered = 0;
		ResultSet res_filtered = new ResultSet();
		for(int i = 0; i < data_pathways.size(); i++) if(res.adj_p_value[i] <= Parameters.probaCutoff) filtered++;
		System.out.println(filtered + " values passed the " + (Parameters.probaCutoff * 100) + "% threshold after p-value adjustment");
		res_filtered.init(filtered);
		int index = 0;
		for(int i = 0; i < data_pathways.size(); i++)
		{
			if(res.adj_p_value[i] <= Parameters.probaCutoff)
			{
				res_filtered.clone(res, i, index);
				index ++;
			}
		}
		// Writing Results
		if(!Parameters.isSilent) System.out.println("Writing the result file...");
		try
		{
			Writer writer = new FileWriter(Parameters.outputFolder + "output.json");
			Gson gson = new GsonBuilder().create();
			gson.toJson(res_filtered, writer);
			writer.close();
		}
		catch(Exception e)
		{
			System.err.println("Problem detected when writing the JSON result file. Stopping program...");
			System.exit(-1);
		}
		if(!Parameters.isSilent) System.out.println("Computation Done!");
		
		if(!Parameters.isSilent) System.out.println("Score & Writing computation time: "+ Utils.toReadableTime(System.currentTimeMillis() - t1));
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
	
    private static double[] p_adjust(final double[] pvalues, String adjMethod)
    {
        if (!adjMethod.equals("BH") && !adjMethod.equals("fdr") && !adjMethod.equals("bonferroni") && !adjMethod.equals("none"))
        {
            System.out.println("This adjustment method is not implemented.");
            System.exit(0);
        }
        if (adjMethod.equals("fdr")) adjMethod = "BH";
        if (adjMethod.equals("none")) return pvalues;

        Integer[] idx = new Integer[pvalues.length];
        for (int i = 0; i < idx.length; i++) idx[i] = i;
        double[] adj_p_value = new double[pvalues.length];

        Arrays.sort(idx, new Comparator<Integer>()
        { // I sort the indexes with respect to the double values
            @Override
            public int compare(Integer o1, Integer o2)
            {
                return Double.compare(pvalues[o2], pvalues[o1]);
            }
        });
        int[] index = new int[idx.length];
        for (int i = 0; i < index.length; i++)
        {
            index[idx[i]] = i;
        }

        for (int i = 0; i < pvalues.length; i++)
        {
            if (adjMethod.equals("BH"))
            {
                adj_p_value[i] = pvalues.length / (pvalues.length - (double) index[i]) * pvalues[i];
            }
            if (adjMethod.equals("bonferroni"))
            {
                adj_p_value[i] = pvalues.length * pvalues[i];
            }
        }

        double min = Double.MAX_VALUE;
        for (int i = 0; i < index.length; i++) // cummin
        {
            double adjP = adj_p_value[idx[i]];
            if (adjP < min)
            {
                min = adjP;
            } else
            {
                adjP = min;
            }
            adj_p_value[idx[i]] = Math.min(1, adjP);
        }

        return adj_p_value;
    }
	
	private static void loadGeneListJSON(String fileName) // JSON file
	{
		try
		{
			Gson gson = new Gson();
			JsonReader reader = new JsonReader(new FileReader(fileName));
			String[][] listGenes = gson.fromJson(reader, String[][] .class); // contains the whole genes lists
			for (int i = 0; i < listGenes.length; i++) 
			{
				for (int j = 0; j <3; j++) // I know there are only three categories
				{
					String[] gene = listGenes[i][j].split(",");
					for(String g:gene)
					{
						String ge = g.toUpperCase();
						if(backgroundGenes.get(ge) != null)
						{
							backgroundGenes.put(ge, true);
							genesToEnrich.add(ge);
						}
					}
				}
			}
			System.out.println(genesToEnrich.size() + " genes were matching other lists to be enriched.");
			reader.close();
		}
		catch(FileNotFoundException nfe)
		{
			System.err.println("The JSON gene list was not found at the given path: " + fileName + "\nStopping program...");
			System.exit(-1);
		}
		catch(Exception e)
		{
			System.out.println(e);
			System.err.println("Problem detected when reading the JSON gene list. Stopping program...");
			System.exit(-1);
		}
	}
	
	private static void loadDataPathways(String fileName) // GMT file
	{
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			String line = br.readLine();
			while(line != null)
			{
				String[] tokens = line.split("\t");
				Pathway p = new Pathway();
				p.id = tokens[0];
				p.description = tokens[1];
				p.url = tokens[2];
				for(int i = 3; i < tokens.length; i++) 
				{
					String gene = tokens[i].toUpperCase();
					Boolean inMatrix = backgroundGenes.get(gene);
					if(inMatrix != null) // If this gene was in the background file
					{
						backgroundGenes.put(gene, true); // If this gene was in the background file
						p.listGenes.add(gene);
					}
				}
				if(p.listGenes.size() <= Parameters.maxGenesInPathway && p.listGenes.size() >= Parameters.minGenesInPathway) data_pathways.put(p.id, p);
				line = br.readLine();
			}
			br.close();
		}
		catch(FileNotFoundException nfe)
		{
			System.err.println("The Pathway file was not found at the given path: " + fileName + "\nStopping program...");
			System.exit(-1);
		}
		catch(Exception e)
		{
			System.err.println("Problem detected when reading the Pathway file. Stopping program...");
			System.exit(-1);
		}
	}
	
	private static void loadBackgroundGenes(String filename)
	{
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String line = br.readLine(); // Header is skipped
			line = br.readLine();
			while(line != null)
			{
				String[] vals = line.substring(0, line.indexOf('\t')).split("\\|");
				if(vals.length > 0) for(String gene:vals[0].split(",")) backgroundGenes.put(gene.toUpperCase(), false); // Load Ens
				if(vals.length > 1) for(String gene:vals[1].split(",")) backgroundGenes.put(gene.toUpperCase(), false); // Load Genes
				if(vals.length > 2) for(String gene:vals[2].split(",")) backgroundGenes.put(gene.toUpperCase(), false); // Load Others
				
				line = br.readLine();
			}
			br.close();
		}
		catch(FileNotFoundException nfe)
		{
			System.err.println("The Background file was not found at the given path: " + filename + "\nStopping program...");
			System.exit(-1);
		}
		catch(Exception e)
		{
			System.err.println("Problem detected when reading the Background file. Stopping program...");
			e.printStackTrace();
			System.exit(-1);
		}
	}
}

class JSONListGenes
{
	public List<List<String>> list_genes;
}

class Pathway
{
	public String id;
	public String url;
	public String description;
	public HashSet<String> listGenes = new HashSet<String>();
}

class ResultSet
{
	double[] p_value;
	double[] adj_p_value;
	double[] OR;
	String[] pathways;
	String[] descriptions;
	String[] urls;
	String warning;
	
	public void init(int length)
	{
		p_value = new double[length];
		adj_p_value = new double[length];
		descriptions = new String[length];
		pathways = new String[length];
		urls = new String[length];
		OR = new double[length];
		warning = null;
	}
	
	public void clone(ResultSet res, int res_index, int this_index)
	{
		p_value[this_index] = res.p_value[res_index];
		adj_p_value[this_index] = res.adj_p_value[res_index];
		descriptions[this_index] = res.descriptions[res_index];
		pathways[this_index] = res.pathways[res_index];
		urls[this_index] = res.urls[res_index];
		OR[this_index] = res.OR[res_index];
		this.warning = res.warning;
	}
}




