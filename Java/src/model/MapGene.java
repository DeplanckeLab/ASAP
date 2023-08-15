package model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import hdf5.loom.LoomData;
import json.ErrorJSON;
import parsing.model.Gene;

public class MapGene 
{
	public static HashMap<String, ArrayList<Gene>> ensembl_db;
	public static HashMap<String, ArrayList<Gene>> gene_db;
	public static HashMap<String, ArrayList<Gene>> alt_db;
	public static HashMap<String, ArrayList<Gene>> obsolete_db;
	
	public static void init()
	{
		ensembl_db = new HashMap<>();
		gene_db = new HashMap<>();
		alt_db = new HashMap<>();
		obsolete_db = new HashMap<>();
	}
	
	public static Gene retrieveLatest(String gene, String ens, ArrayList<Gene> genes) // If two are equal, check the font case, then retrieve the first
	{
		if(genes.size() == 0) return null;
		if(genes.size() == 1) return genes.get(0);
		// More than 1 hit
		
		// First check latest Ensembl release
		ArrayList<Gene> latest = new ArrayList<Gene>();
		int max = 0;
		for(Gene g:genes) if(g.latest_ensembl_release > max) max = g.latest_ensembl_release;
		for(Gene g:genes) if(g.latest_ensembl_release == max) latest.add(g);

		// In case there are still many, check the raw font case for any difference
		if(latest.size() > 1)
		{
			for(Gene g:genes) if(gene != null && g.name.equals(gene)) return g; // Take first hit with exact same font case
		}
			
		return latest.get(0); // Take first one
	}
	
	public static void parseGenes(LoomData data) // Parsing all genes
    {
		try
		{
			for(long i = 0; i < data.nber_genes; i++) parseGene(data, i);
		}
		catch(IndexOutOfBoundsException iooe)
		{
			new ErrorJSON(iooe.getMessage());
		}
    }
	
	/***
	 * 
	 * @param gene
	 * @return 1 if Ensembl, 0 if Gene, -1 if 10 first were not found in any DB
	 */
	public static int isEnsembl(String[] gene)
	{
		for(int i = 0; i < 10; i++) // Check the 10 first genes max
		{
			ArrayList<Gene> dbHit = MapGene.ensembl_db.get(gene[i].toUpperCase()); // First test the Ensembl DB
			if(dbHit != null) return 1; // It's an Ensembl ID. So I consider the whole thing to be Ensembl genes.
			// If not, I check the Gene DB
			dbHit = MapGene.gene_db.get(gene[i].toUpperCase());
			if(dbHit != null) return 0; // It's a Gene ID. So I consider the whole thing to be Gene ids
			// if not found in any DB, I try with the next one
		}
		// If I reach here it means that none of the 10 first element was found, in any DB
		return -1;
	}
	
	public static void parseGene(LoomData data, long index) // Only parsing the ith gene
    {
		String ens = data.ens_names.get(index);
		if(ens != null && ens.equals("")) ens = null;
		String gene = data.gene_names.get(index);
		if(gene != null && gene.equals("")) gene = null;
		
    	// Get infos from DB (Ensembl name, alt names, etc...)
    	ArrayList<Gene> dbHit = null;
    	Gene gHit = null;
    	if(ens != null)
    	{
    		dbHit = MapGene.ensembl_db.get(ens.toUpperCase()); // First test the Ensembl ID
    		if(dbHit == null) 
    		{
    			dbHit = MapGene.ensembl_db.get(ens.replaceFirst("\\.\\d+", "").toUpperCase()); // Test the ensembl on the geneDB
    			if(dbHit != null) ens = ens.replaceFirst("\\.\\d+", "");
    		}
    	}
    	if(dbHit == null && ens != null) dbHit = MapGene.gene_db.get(ens.toUpperCase()); // Test the ensembl on the geneDB
    	if(dbHit == null && gene != null) dbHit = MapGene.ensembl_db.get(gene.toUpperCase()); // Test the gene on the EnsemblDB
    	if(dbHit == null && gene != null) dbHit = MapGene.gene_db.get(gene.toUpperCase()); // Else check the HGNC DB
    	if(dbHit == null && ens != null) dbHit = MapGene.alt_db.get(ens.toUpperCase()); // Test the ensembl on the alt_gene DB (last hope)
    	if(dbHit == null && gene != null) dbHit = MapGene.alt_db.get(gene.toUpperCase()); // Test the gene on the alt_gene DB (ultimate last hope)
    	if(dbHit == null && ens != null) dbHit = MapGene.obsolete_db.get(ens.toUpperCase()); // Test the ensembl on the alt_gene DB (ultimate last hope v2)
    	if(dbHit == null && gene != null) dbHit = MapGene.obsolete_db.get(gene.toUpperCase()); // Test the gene on the alt_gene DB (desperate move)
    	if(dbHit == null) // Not Found at all
    	{
			gHit = new Gene();
			if(ens == null && gene == null)
			{
				gene = new StringBuffer("Gene_").append(index + 1).toString();
				ens = gene;
			}
			if(ens == null) ens = "";
			if(gene == null) gene = ens; // Ensure Gene is never Null
			gHit.ensembl_id = ens;
			gHit.name = gene;
			gHit.alt_names = new HashSet<>();
			gHit.biotype = "__unknown";
			gHit.chr = "__unknown";
			gHit.sum_exon_length = 0;
			data.nber_not_found_genes++;
		}

    	// Retrieve gene with latest release in Ensembl
    	if(dbHit != null) gHit = MapGene.retrieveLatest(gene, ens, dbHit);
		
		// Update the final entries for the Loom
    	if(gHit.name == null) new ErrorJSON("Gene " + data.original_gene_names.get(index) + " is NULL from DB???");
		data.ens_names.set(index, gHit.ensembl_id);
		data.gene_names.set(index, gHit.name);
		data.biotypes.set(index, gHit.biotype);
		data.chromosomes.set(index, gHit.chr);
		data.sumExonLength.set(index, gHit.sum_exon_length);
    }
	
}
