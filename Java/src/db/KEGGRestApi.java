package db;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;

public class KEGGRestApi 
{
	public static HashMap<String, String> keggDescription = null;
	public static HashMap<String, ArrayList<String>> keggAnnot = null;
	public static HashMap<String, String> keggIDtoGene = null;
	
	public static void main(String[] args)  throws IOException 
	{
		generateKEGGDB("hsa", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/kegg.hsa.gmt"); // mmu, rno (rat), dme (droso mela)
		generateKEGGDB("mmu", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/kegg.mmu.gmt"); // mmu, rno (rat), dme (droso mela)
	}
	
	public static void generateKEGGDB(String organism, String outputGMTFile) throws IOException
	{
		keggDescription = new HashMap<String, String>();
		keggAnnot = new HashMap<String, ArrayList<String>> ();
		keggIDtoGene = new HashMap<String, String>();
		fetchKEGGPathways(organism);
		System.out.println(keggDescription.size() + " KEGG pathways were found.");
		for(String keggID:keggDescription.keySet()) fetchGenesInPathway(keggID);
		// Fetch Gene names from KEGG ids
		geneNameConversion(organism);
		// Create GMT from REST data
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputGMTFile));
		for(String geneset:keggAnnot.keySet())
		{
			ArrayList<String> genes = keggAnnot.get(geneset);
			bw.write(geneset+"\t"+keggDescription.get(geneset)+"\thttp://www.genome.jp/dbget-bin/www_bget?"+geneset);
			for (String gene:genes) bw.write("\t"+keggIDtoGene.get(gene)); // List gene
			bw.write("\n");
		}
		bw.close();
	}
	
    public static void geneNameConversion(String organism) throws IOException // http://rest.kegg.jp/get/hsa:10458 gives everything, including Ensembl
    {
        URL url = new URL("http://rest.kegg.jp/list/"+organism);
        URLConnection urlc = url.openConnection();
        urlc.setDoOutput(true);
        urlc.setAllowUserInteraction(false);

        //get result
        BufferedReader br = new BufferedReader(new InputStreamReader(urlc.getInputStream()));
        String l = null;
        while ((l=br.readLine())!=null) 
        {
        	String[] results = l.split("\t");
        	String gene = results[1].replaceAll("uncharacterized ", "").replaceAll("protein ", "").split(";")[0].split(",")[0].toUpperCase().trim();
        	// Handle problems in mapping in KEGG hsa
        	switch(results[0])
        	{
        		case "hsa:105369274": gene = "LOC105369274"; break;
            	case "hsa:102723532": gene = "AC171558.1"; break;
            	case "hsa:100288562": gene = "LOC100288562"; break;
            	case "hsa:101929601": gene = "LOC101929601"; break;
            	case "hsa:101930111": gene = "LOC101930111"; break;
            	case "hsa:101929627": gene = "LOC101929627"; break;
            	case "hsa:643802": gene = "LOC643802"; break;
            	case "hsa:81691":  gene = "AC004381.6"; break;
            	case "hsa:102724428": gene = "CH507-42P11.8"; break;
            	case "hsa:102725035": gene = "LOC102725035"; break;
            	case "hsa:102724788": gene = "LOC102724788"; break;
            	case "hsa:102800317": gene = "CSNK1E"; break;
            	default:
    	        	gene = gene.replaceAll("DNA-DIRECTED RNA POLYMERASE III SUBUNIT RPC5", "LOC101060521"); // Bad annotation
    	        	gene = gene.replaceAll("EUKARYOTIC ELONGATION FACTOR 2 KINASE", "LOC101930123");
    	        	gene = gene.replaceAll("TBC1D7-LOC100130357 READTHROUGH", "TBC1D7");
    	        	gene = gene.replaceAll("ADENYLATE KINASE ISOENZYME 1-LIKE", "LOC390877");
    	        	gene = gene.replaceAll("CYTOCHROME P450 2D6-LIKE", "LOC107987478");
    	        	gene = gene.replaceAll("PROLINE DEHYDROGENASE 1", "AC007325.2");
    	        	gene = gene.replaceAll("NEUROPEPTIDE Y RECEPTOR TYPE 4", "LOC105379861");
    	        	gene = gene.replaceAll("PROTHYMOSIN ALPHA", "LOC100506248");
    	        	gene = gene.replaceAll("LEUKOSIALIN", "LOC105369247");
    	        	gene = gene.replaceAll("ICOS LIGAND", "CH507-9B2.1");
    	        	gene = gene.replaceAll("CYTOCHROME P450 2D6-LIKE", "LOC107987478");
    	        	gene = gene.replaceAll("STRUCTURE-SPECIFIC ENDONUCLEASE SUBUNIT SLX1", "LOC105369236"); // 177 
    	        	gene = gene.replaceAll("MULTIDRUG RESISTANCE-ASSOCIATED 6", "LOC105369239");
    	        	gene = gene.replaceAll("40S RIBOSOMAL S26", "LOC101929876");
        	}
        	keggIDtoGene.put(results[0], gene);
        }
        br.close();
    }

    public static void fetchGenesInPathway(String pathway) throws IOException // see http://rest.kegg.jp/list/organism for list of organisms
    {
        URL url = new URL("http://rest.kegg.jp/link/genes/" + pathway);
        URLConnection urlc = url.openConnection();
        urlc.setDoOutput(true);
        urlc.setAllowUserInteraction(false);

        //get result
        BufferedReader br = new BufferedReader(new InputStreamReader(urlc.getInputStream()));
        String l = null;
        ArrayList<String> genes = new ArrayList<String>();
        while ((l=br.readLine())!=null) genes.add(l.split("\t")[1]);
        keggAnnot.put(pathway, genes);
        br.close();
    }
    
    public static void fetchKEGGPathways(String organism) throws IOException // see http://rest.kegg.jp/list/organism for list of organisms
    {
        URL url = new URL("http://rest.kegg.jp/list/pathway/" + organism);
        URLConnection urlc = url.openConnection();
        urlc.setDoOutput(true);
        urlc.setAllowUserInteraction(false);

        //get result
        BufferedReader br = new BufferedReader(new InputStreamReader(urlc.getInputStream()));
        String l = null;
        while ((l=br.readLine())!=null) 
        {
        	String[] results = l.split("\t");
        	keggDescription.put(results[0].substring(5), results[1].substring(0, results[1].lastIndexOf(" -")));
        }
        br.close();
    }
}
