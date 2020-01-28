package model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import bigarrays.StringArray64;
import hdf5.loom.LoomData;
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
	
	public static Gene retrieveLatest(ArrayList<Gene> genes) // If two are equal, retrieve the first
	{
		Gene res = null;
		int max = 0;
		for(Gene g:genes) 
		{
			if(g.latest_ensembl_release > max)
			{
				max = g.latest_ensembl_release;
				res = g;
			}
		}
		return res;
	}
	
	public static void parseGenes(LoomData data) // Parsing all genes
    {
		for(long i = 0; i < data.nber_genes; i++) parseGene(data, i);
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
    	if(ens != null) dbHit = MapGene.ensembl_db.get(ens.toUpperCase()); // First test the Ensembl ID
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
			if(ens == null && gene == null) gene = new StringBuffer("Gene_").append(index + 1).toString();
			gHit.ensembl_id = (ens == null)?"":ens;
			gHit.name  = (gene == null)?"":gene;
			gHit.alt_names = new HashSet<>();
			gHit.biotype = "__unknown";
			gHit.chr = "__unknown";
			gHit.sum_exon_length = 0;
			data.nber_not_found_genes++;
		}

    	// Retrieve gene with latest release in Ensembl
    	if(dbHit != null) gHit = MapGene.retrieveLatest(dbHit);
		
		// Update the final entries for the Loom
		data.ens_names.set(index, gHit.ensembl_id);
		data.biotypes.set(index, gHit.biotype);
		data.chromosomes.set(index, gHit.chr);
		data.sumExonLength.set(index, gHit.sum_exon_length);
		data.gene_names.set(index, gHit.name);
    }
	
	public static void addGene(String gene, LoomData data, long index)
	{
		if(data.gene_names == null) data.gene_names = new StringArray64(data.nber_genes);
		data.gene_names.set(index, gene);
	}
	
	
}
