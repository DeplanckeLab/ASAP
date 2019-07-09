package model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import bigarrays.IntArray64;
import bigarrays.StringArray64;
import hdf5.loom.LoomData;
import parsing.model.Gene;

public class MapGene 
{
	public static HashMap<String, ArrayList<Gene>> ensembl_db;
	public static HashMap<String, ArrayList<Gene>> gene_db;
	public static HashMap<String, ArrayList<Gene>> alt_db;
	
	public static void init()
	{
		ensembl_db = new HashMap<>();
		gene_db = new HashMap<>();
		alt_db = new HashMap<>();
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
	
	public static void parseGene(LoomData data, int index) // Based on 10x Ensembl Genes // json.loom.ens_names[index] is current gene to parse (Ensembl Name)
	{	
		// Handle duplicated names
		String ens = data.ens_names.get(index);
		
		// Get infos from DB (Ensembl name, alt names, etc...)
		ArrayList<Gene> dbHit = MapGene.ensembl_db.get(ens);
		Gene gHit = null;
		if(dbHit != null) gHit = MapGene.retrieveLatest(dbHit);
		else // This is not an EnsemblID
		{
			dbHit = MapGene.gene_db.get(ens);
			if(dbHit != null) gHit = MapGene.retrieveLatest(dbHit);
			else // Last Hope in ALT genes
			{
				dbHit = MapGene.alt_db.get(ens);
				if(dbHit != null) gHit = MapGene.retrieveLatest(dbHit);
				else // Not Found at all
				{
					gHit = new Gene();
					gHit.ensembl_id = ens;
					gHit.name = data.gene_names.get(index);
					gHit.alt_names = new HashSet<>();
					gHit.biotype = "__unknown";
					gHit.chr = "__unknown";
					gHit.sum_exon_length = -1;
					data.nber_not_found_genes++;
				}
			}
		}
		
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
	
	public static void parseGenes(LoomData data) // Parsing all genes
    {
		for(long i = 0; i < data.nber_genes; i++) parseGene(data, i);
    }
	
	public static boolean parseGene(LoomData data, long index) // Only parsing the ith gene
    {
		String ens = data.ens_names.get(index);
		String gene = data.gene_names.get(index);

    	if(gene.startsWith("__")) // HTseq annotation
    	{
    		switch(gene)
    		{
    			case "__no_feature":
    				if(data.__no_feature == null) data.__no_feature = new IntArray64(data.nber_cells);
    				data.removed.add(index); // plan to remove this row from final matrix
    	    		return false;
    			case "__ambiguous":
    				if(data.__ambiguous == null) data.__ambiguous = new IntArray64(data.nber_cells);
    				data.removed.add(index); // plan to remove this row from final matrix
    	    		return false;
    			case "__too_low_aQual":
    				if(data.__too_low_aQual == null) data.__too_low_aQual = new IntArray64(data.nber_cells);
    				data.removed.add(index); // plan to remove this row from final matrix
    	    		return false;
    			case "__not_aligned":
    				if(data.__not_aligned == null) data.__not_aligned = new IntArray64(data.nber_cells);
    				data.removed.add(index); // plan to remove this row from final matrix
    	    		return false;
    			case "__alignment_not_unique":
    				if(data.__alignment_not_unique == null) data.__alignment_not_unique = new IntArray64(data.nber_cells);
    				data.removed.add(index); // plan to remove this row from final matrix
    	    		return false;
    	    	default:
    	    		// Do nothing
    		}
    	}
    	
    	// ERCCs
    	if(gene.startsWith("ERCC-")) // This is recognized as an ERCC
		{
			if(data.erccs == null) data.erccs = new ERCC((int)data.nber_cells);
			data.erccs.callERCC(index, gene); // Save this index as an ERCC
			data.nber_ercc++;
			data.removed.add(index); // plan to remove this row from final matrix
			return false;
		}
    	
    	// Get infos from DB (Ensembl name, alt names, etc...)
    	ArrayList<Gene> dbHit = null;
    	Gene gHit = null;
    	if(!ens.equals("")) dbHit = MapGene.ensembl_db.get(ens.toUpperCase()); // First test the Ensembl ID
    	if(dbHit == null && !gene.equals("")) dbHit = MapGene.ensembl_db.get(gene.toUpperCase()); // Test the gne on the EnsemblDB
    	if(dbHit == null && !gene.equals("")) dbHit = MapGene.gene_db.get(gene.toUpperCase()); // Else check the HGNC DB
    	if(dbHit == null) // Not found in HGNC nor Ensembl
    	{
    		if(!gene.equals("")) dbHit = MapGene.alt_db.get(gene.toUpperCase()); // Last Hope in ALT genes
    		if(dbHit == null) // Not Found at all
    		{
    			gHit = new Gene();
    			if(!ens.equals("")) gHit.ensembl_id = ens; // Just not in our DB for unknown reason
    			else gHit.ensembl_id = "";
    			if(!gene.equals("")) gHit.name = gene;
    			else gHit.name = new StringBuffer("Gene_").append(index + 1).toString();
    			gHit.alt_names = new HashSet<>();
    			gHit.biotype = "__unknown";
    			gHit.chr = "__unknown";
    			gHit.sum_exon_length = -1;
    			data.nber_not_found_genes++;
    		}
    	}

    	// Retrieve gene with latest release in Ensembl
    	if(dbHit != null) gHit = MapGene.retrieveLatest(dbHit);
		
		// Update the final entries for the Loom
		data.ens_names.set(index, gHit.ensembl_id);
		data.biotypes.set(index, gHit.biotype);
		data.chromosomes.set(index, gHit.chr);
		data.sumExonLength.set(index, gHit.sum_exon_length);
		data.gene_names.set(index, gHit.name);
	
    	return true; // Gene parsing success
    }
}
