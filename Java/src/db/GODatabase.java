package db;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.ProgressMonitor;

import tools.Utils;

class Specie
{
	String common_name;
	String ncbi_taxa_id;
	String genus;
	String species;
	String taxonomic_rank;
}

public class GODatabase 
{
   private static HashMap<String, GOTerm> termsBP = null;
   private static HashMap<String, GOTerm> termsMF = null;
   private static HashMap<String, GOTerm> termsCC = null;

   public static void generateGODB() 
   {
	   // Connect to ASAP DB
	   
	   DBManager.connect();
	   HashMap<Integer, Integer> taxonsInDB = DBManager.listAllOrganisms();
	   System.out.println("There are " + taxonsInDB.size() + " organisms in ASAP DB.");
	   DBManager.disconnect();
	   
	   // Connect to GO DB
	   DBManager.JDBC_DRIVER = "com.mysql.jdbc.Driver";
	   DBManager.URL = "jdbc:mysql://mysql-eg-publicsql.ebi.ac.uk:4157/go_latest?user=go_select";
	   DBManager.connect();
	   
	   // Get all GO terms (common for all taxons)
	   fetchTerms();
	   System.out.println(termsBP.size() + " GO BP terms were fetched.");
	   System.out.println(termsMF.size() + " GO MF terms were fetched.");
	   System.out.println(termsCC.size() + " GO CC terms were fetched.");
	   
	   /*for(String goterm:termsCC.keySet())
	   {
		   GOTerm go = termsCC.get(goterm);
		   go.genes = fetchGOGenesDirect(go.id, taxonId);
		   go.descendants = fetchDescendant(go.id);
		   p.increment();
	   }
		//GODatabase.generateGODB(Parameters.outputFolder, Parameters.taxon);
	   generateGODB(Parameters.outputFolder + "go." + taxon + ".gmt", int taxonId)*/
   }
   
   public static void generateGODB(String outputGMTFolder, int taxonId)
   {
	   try
	   {
		   System.out.println("Computing ("+ taxonId + ")...");
		   if(!outputGMTFolder.endsWith("/")) outputGMTFolder += "/";
		   fetchTerms();
		   System.out.println(termsBP.size() + " GO BP terms were fetched.");
		   System.out.println(termsMF.size() + " GO MF terms were fetched.");
		   System.out.println(termsCC.size() + " GO CC terms were fetched.");
		   
		   long t1 = System.currentTimeMillis();
		   ProgressBar p = new ProgressBar("GO Cellular Component", termsCC.size());
		   BufferedWriter bw = new BufferedWriter(new FileWriter(outputGMTFolder + "GO_CC_"+taxonId+".gmt"));
		   for(String goterm:termsCC.keySet())
		   {
			   GOTerm go = termsCC.get(goterm);
			   ArrayList<String> genes = fetchGOGenesRecursively(go.id, taxonId);
			   if(genes.size() > 0) // If the GO is empty then it's useless to print it
			   {
				   bw.write(go.id + "\t" + go.description + "\t" + "http://amigo.geneontology.org/amigo/term/"+go.id);
				   for(String gene:genes) bw.write("\t" + gene);
				   bw.write("\n");
			   }
			   p.increment();
		   }
		   p.close();
		   bw.close();
		   System.out.println("GO CC processing time: "+ Utils.toReadableTime(System.currentTimeMillis() - t1));
		   
		   t1 = System.currentTimeMillis();
		   p = new ProgressBar("GO Molecular Function", termsMF.size());
		   bw = new BufferedWriter(new FileWriter(outputGMTFolder + "GO_MF_"+taxonId+".gmt"));
		   for(String goterm:termsMF.keySet())
		   {
			   GOTerm go = termsMF.get(goterm);
			   ArrayList<String> genes = fetchGOGenesRecursively(go.id, taxonId);
			   if(genes.size() > 0) // If the GO is empty then it's useless to print it
			   {
				   bw.write(go.id + "\t" + go.description + "\t" + "http://amigo.geneontology.org/amigo/term/"+go.id);
				   for(String gene:genes) bw.write("\t" + gene);
				   bw.write("\n");
			   }
			   p.increment();
		   }
		   p.close();
		   bw.close();
		   System.out.println("GO MF processing time: "+ Utils.toReadableTime(System.currentTimeMillis() - t1));
		   
		   t1 = System.currentTimeMillis();
		   p = new ProgressBar("GO Biological Process", termsBP.size());
		   bw = new BufferedWriter(new FileWriter(outputGMTFolder + "GO_BP_"+taxonId+".gmt"));
		   for(String goterm:termsBP.keySet())
		   {
			   GOTerm go = termsBP.get(goterm);
			   ArrayList<String> genes = fetchGOGenesRecursively(go.id, taxonId);
			   if(genes.size() > 0) // If the GO is empty then it's useless to print it
			   {
				   bw.write(go.id + "\t" + go.description + "\t" + "http://amigo.geneontology.org/amigo/term/"+go.id);
				   for(String gene:genes) bw.write("\t" + gene);
				   bw.write("\n");
			   }
			   p.increment();
		   }
		   p.close();
		   bw.close();
		   System.out.println("GO BP processing time: "+ Utils.toReadableTime(System.currentTimeMillis() - t1));
	   }
	   catch(IOException ioe)
	   {
		   ioe.printStackTrace();
	   }
   }
   
   public static void fetchTerms()
   {
	   termsBP = new HashMap<String, GOTerm>();
	   termsMF = new HashMap<String, GOTerm>();
	   termsCC = new HashMap<String, GOTerm>();
	   Statement stmt = null;
	   try
	   {   
		   	stmt = DBManager.conn.createStatement();
	      
		   	String sql = "SELECT * FROM term";
			ResultSet rs = stmt.executeQuery(sql);
		    while(rs.next())
		    {
		    	GOTerm go = new GOTerm();
		    	if(rs.getInt("is_obsolete") == 0)
		    	{
		    		go.id = rs.getString("acc");
		    		go.description = rs.getString("name");
			    	String term_type = rs.getString("term_type");
			    	if(term_type.equals("biological_process")) termsBP.put(go.id, go);
			    	else if(term_type.equals("molecular_function")) termsMF.put(go.id, go);
			    	else if(term_type.equals("cellular_component")) termsCC.put(go.id, go);
		    	}
		    }
		    rs.close();
	   }
	   catch(Exception e)
	   {
		   e.printStackTrace();
	   }
	   finally
	   {
		   try
		   {
			   if(stmt!=null) stmt.close();
		   }
		   catch(SQLException se2){ }// nothing we can do
	   }
   }
   
   public static ArrayList<String> fetchGOGenesRecursively(String goId, int taxon)
   {
	   ArrayList<String> genes = new ArrayList<String>();
	   Statement stmt = null;
	   try
	   {   
		   	stmt = DBManager.conn.createStatement();
	      
		   	String sql = ""
		   			+ "SELECT DISTINCT gene_product.symbol FROM term"
		   			+ " INNER JOIN graph_path ON (term.id=graph_path.term1_id)" // For all genes from subtrees
		   			+ " INNER JOIN association ON (graph_path.term2_id=association.term_id)"
		   			+ " INNER JOIN gene_product ON (association.gene_product_id=gene_product.id)"
		   			+ " INNER JOIN species ON (gene_product.species_id=species.id)"
		   			+ " INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id)" // For restricting results to Uniprot and not PDB
		   			+ " WHERE term.acc = '"+ goId +"' AND  ncbi_taxa_id = "+ taxon; // dbxref.xref_dbname = 'UniProtKB' AND
		   			//+ " INNER JOIN evidence ON (association.id=evidence.association_id)"
		   			//+ " WHERE evidence.code='IEA' AND dbxref.xref_dbname = 'UniProtKB' AND term.acc = '"+ goId +"' AND  ncbi_taxa_id = "+ taxon;
			ResultSet rs = stmt.executeQuery(sql);
		    while(rs.next()) genes.add(rs.getString("symbol").toUpperCase().trim());
		    rs.close();
	   }
	   catch(Exception e)
	   {
		   e.printStackTrace();
	   }
	   finally
	   {
		   try
		   {
			   if(stmt!=null) stmt.close();
		   }
		   catch(SQLException se2){ }// nothing we can do
	   }
	   return genes;
   }
   
   public static ArrayList<String> fetchGOGenesDirect(String goId, int taxon)
   {
	   ArrayList<String> genes = new ArrayList<String>();
	   Statement stmt = null;
	   try
	   {   
		   	stmt = DBManager.conn.createStatement();
	      
		   	String sql = ""
		   			+ "SELECT DISTINCT gene_product.symbol FROM term"
		   			+ " INNER JOIN association ON (term.id=association.term_id)"
		   			+ " INNER JOIN gene_product ON (association.gene_product_id=gene_product.id)"
		   			+ " INNER JOIN species ON (gene_product.species_id=species.id)"
		   			+ " INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id)" // For restricting results to Uniprot and not PDB
		   			+ " WHERE dbxref.xref_dbname = 'UniProtKB' AND term.acc = '"+ goId +"' AND  ncbi_taxa_id = "+ taxon;
		   			//+ " INNER JOIN evidence ON (association.id=evidence.association_id)"
		   			//+ " WHERE evidence.code='IEA' AND dbxref.xref_dbname = 'UniProtKB' AND term.acc = '"+ goId +"' AND  ncbi_taxa_id = "+ taxon;
			ResultSet rs = stmt.executeQuery(sql);
		    while(rs.next()) genes.add(rs.getString("gene_product.symbol").toUpperCase().trim());
		    rs.close();
	   }
	   catch(Exception e)
	   {
		   e.printStackTrace();
	   }
	   finally
	   {
		   try
		   {
			   if(stmt!=null) stmt.close();
		   }
		   catch(SQLException se2){ }// nothing we can do
	   }
	   return genes;
   }
   
   public static ArrayList<String> fetchDescendant(String goId) // Get nodes, subnodes, subsubnodes, etc...
   {
	   // TODO: eventually get only direct descendants? distance = 1?
	   ArrayList<String> go = new ArrayList<String>();
	   Statement stmt = null;
	   try
	   {   
		   	stmt = DBManager.conn.createStatement();
		   	String sql = ""
		   			+ "SELECT DISTINCT descendant.acc FROM term"
		   			+ " INNER JOIN graph_path ON (term.id=graph_path.term1_id)" // For all genes from subtrees
		   			+ " INNER JOIN term AS descendant ON (descendant.id=graph_path.term2_id)"
		   			+ " WHERE term.acc = '"+ goId +"' AND distance <> 0";
			ResultSet rs = stmt.executeQuery(sql);
		    while(rs.next()) go.add(rs.getString("descendant.acc").trim());
		    rs.close();
	   }
	   catch(Exception e)
	   {
		   e.printStackTrace();
	   }
	   finally
	   {
		   try
		   {
			   if(stmt!=null) stmt.close();
		   }
		   catch(SQLException se2){ }// nothing we can do
	   }
	   return go;
   }
}

class GOTerm
{
	String id;
	String description;
	ArrayList<String> genes;
	ArrayList<String> descendants;
}

class ProgressBar
{
	public ProgressMonitor progressBar;
	public int max;
	public int step;
	
	public ProgressBar(String dbname, int max)
	{
		this.max = max;
		this.step = 0;
	    progressBar = new ProgressMonitor(null, "Fetching database "+dbname, null, 0, max);
	}
	
	public void increment()
	{
		step++;
	    progressBar.setProgress(step);
	    progressBar.setNote(step + " GOCC Terms fetched over "+max);
	    if(progressBar.isCanceled()) close();
	}
	
	public void close()
	{
		progressBar.close();
	}
}