package db;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import json.ErrorJSON;
import model.MapGene;
import parsing.model.Gene;

public class DBManager 
{
	public static String JDBC_DRIVER = null;
	public static String URL = null;
	public static Connection conn = null;
	
	public static void generateEnrichmentDB()
	{
			/*DrugBankParser.parse("hsa", "C:/Users/gardeux/Dropbox/GeneSets/DrugBank/drugbank.xml", "C:/Users/gardeux/Dropbox/ASAP/Scripts/Enrichment/Genesets/drugbank_hsa.gmt");
			DrugBankParser.parse("mmu", "C:/Users/gardeux/Dropbox/GeneSets/DrugBank/drugbank.xml", "C:/Users/gardeux/Dropbox/ASAP/Scripts/Enrichment/Genesets/drugbank_mmu.gmt");
			KEGGRestApi.generateKEGGDB("hsa", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/kegg_hsa.gmt"); // mmu, rno (rat), dme (droso mela)
			KEGGRestApi.generateKEGGDB("mmu", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/kegg_mmu.gmt"); // mmu, rno (rat), dme (droso mela)
			GeneAtlas.enrichrToGMT("C:/Users/Vincent/Desktop/GeneSets/GeneAtlas/Human_Gene_Atlas.txt", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/gene_atlas.hsa.gmt");
			GeneAtlas.enrichrToGMT("C:/Users/Vincent/Desktop/GeneSets/GeneAtlas/Mouse_Gene_Atlas.txt", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/gene_atlas.mmu.gmt");
			//GODatabase.generateGODB("C:/Users/gardeux/Dropbox/ASAP/Scripts/Enrichment/Genesets/", "hsa");
			//GODatabase.generateGODB("C:/Users/gardeux/Dropbox/ASAP/Scripts/Enrichment/Genesets/", "mmu");*/
			//GODatabase.generateGODB(Parameters.outputFolder, Parameters.taxon);
	}
	
	public static int getTaxonFromOrganismID(int organism_id) 
	{
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "SELECT tax_id,ensembl_db_name FROM organisms WHERE id="+organism_id;
			ResultSet rs = stmt.executeQuery(sql);
			int res = -1;
			String name = null;
			boolean found = rs.next();
			if(found) 
			{
				res = rs.getInt("tax_id");
				name = rs.getString("ensembl_db_name");
			}
			stmt.close();
			if(!found) new ErrorJSON("This organism with id=" + organism_id + " is not in DB");
			if(res == 0) new ErrorJSON("No taxon available for organism id=" + organism_id + ", name=" + name);
			return res;
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
		return -1;
	}
	
	public static int getOrCreateOrganismId(String organism, int typeID)
	{
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "SELECT id FROM organisms WHERE ensembl_db_name='"+organism+"' AND ensembl_subdomain_id='"+typeID+"'";
			ResultSet rs = stmt.executeQuery(sql);
			if(rs.next()) 
			{
				int res = rs.getInt("id");
				rs.close();
				stmt.close();
				return res;
			}
			rs.close();
			System.out.println("No record found. Creating organism " + organism);
			sql = "INSERT INTO organisms (ensembl_db_name, ensembl_subdomain_id) VALUES ('"+organism+"','"+typeID+"')";
			stmt.executeUpdate(sql);
			sql = "SELECT id FROM organisms WHERE ensembl_db_name='"+organism+"' AND ensembl_subdomain_id='"+typeID+"'";
			rs = stmt.executeQuery(sql);
			if(rs.next()) 
			{
				int res = rs.getInt("id");
				rs.close();
				stmt.close();
				return res;
			}
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
		return -1;
	}
	
	public static String getOrganismName(int organism_id)
	{
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "SELECT ensembl_db_name FROM organisms WHERE id="+organism_id;
			ResultSet rs = stmt.executeQuery(sql);
			if(rs.next()) 
			{
				String res = rs.getString("ensembl_db_name");
				stmt.close();
				return res;
			}
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
			catch(SQLException se2){ } // nothing we can do
		}
		return null;
	}
	
	public static void insertOrUpdateGene(String ensembl_id, String name, String biotype, String chr, long gene_length, long sum_exon_length, int organism_id, HashSet<String> alt_names, int latest_ensembl_release)
	{
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = new StringBuilder("SELECT * FROM genes WHERE ensembl_id='").append(ensembl_id).append("' AND organism_id=").append(organism_id).toString();
			ResultSet rs = stmt.executeQuery(sql);
			if(alt_names == null) alt_names = new HashSet<String>();
			if(rs.next()) // If already in DB
			{
				String name_old = rs.getString("name");
				String[] alt_names_old = rs.getString("alt_names").split(",");
				for(String s:alt_names_old) alt_names.add(s);
				if(!name.equals(name_old)) alt_names.add(name_old);
				sql = new StringBuilder("UPDATE genes SET name='").append(name).append("',biotype='").append(biotype).append("',chr='").append(chr).append("',gene_length=").append(gene_length).append(",sum_exon_length=").append(sum_exon_length).append(",alt_names='").append(Gene.buildAltNamesString(alt_names)).append("',latest_ensembl_release=").append(latest_ensembl_release).append(" WHERE ensembl_id='").append(ensembl_id).append("' AND organism_id=").append(organism_id).toString(); // Replace previous entry with the new version
				stmt.executeUpdate(sql);
			}
			else
			{
				//String.format("UPDATE %s SET userid = '%s' WHERE tkey = 1000", TABLE_NAME, adminUserId)
				if(name != null) name = name.replaceAll("'", "\\\\'");
				sql = new StringBuilder("INSERT INTO genes (ensembl_id, name, biotype, chr, gene_length, sum_exon_length, organism_id, alt_names, latest_ensembl_release) VALUES ('").append(ensembl_id).append("','").append(name).append("','").append(biotype).append("','").append(chr).append("',").append(gene_length).append(",").append(sum_exon_length).append(",").append(organism_id).append(",'").append(Gene.buildAltNamesString(alt_names)).append("',").append(latest_ensembl_release).append(")").toString();
				stmt.executeUpdate(sql);
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		finally
		{
			try
			{
				if(stmt != null) stmt.close();
			}
			catch(SQLException se2){ }// nothing we can do
		}
	}	
	
	public static int getNbGenesForOrganism(int organism_id)
	{
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "SELECT count(*) AS total FROM genes WHERE organism_id="+organism_id+";";
			ResultSet rs = stmt.executeQuery(sql);
			if(rs.next()) // Should always work
			{
				int res = rs.getInt("total");
				rs.close();
				stmt.close();
				return res;
			}
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
		return 0;
	}

	public static int getSubdomainID(String subdomain)
	{
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "SELECT id FROM ensembl_subdomains WHERE name='"+subdomain+"'";
			ResultSet rs = stmt.executeQuery(sql);
			if(rs.next())
			{
				int res = rs.getInt("id");
				rs.close();
				stmt.close();
				return res;
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		finally
		{
			try
			{
				if(stmt != null) stmt.close();
			}
			catch(SQLException se2){ }// nothing we can do
		}
		return -1;
	}
	
	public static int getLastRelease(int organism_id)
	{
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "SELECT latest_ensembl_release FROM organisms WHERE id="+organism_id; // Restrict to this organism
			ResultSet rs = stmt.executeQuery(sql);
			if(rs.next())
			{
				int res = rs.getInt("latest_ensembl_release");
				stmt.close();
				return res;
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		finally
		{
			try
			{
				if(stmt != null) stmt.close();
			}
			catch(SQLException se2){ }// nothing we can do
		}
		return 0;
	}
	
	public static void setLastRelease(int organism_id, int release)
	{
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "UPDATE organisms SET latest_ensembl_release = " + release + " WHERE id="+organism_id; // Restrict to genes of this organism
			stmt.executeUpdate(sql);
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
	
	public static void setLastRelease(String type, int release)
	{
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "UPDATE ensembl_subdomains SET latest_ensembl_release = " + release + " WHERE name='" + type + "'"; // Restrict to genes of this organism
			stmt.executeUpdate(sql);
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
	
	public static int getLastRelease(String type)
	{
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "SELECT latest_ensembl_release AS latest FROM ensembl_subdomains WHERE name='"+type+"'";
			ResultSet rs = stmt.executeQuery(sql);
			if(rs.next())
			{
				int res = rs.getInt("latest");
				rs.close();
				stmt.close();
				return res;
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		finally
		{
			try
			{
				if(stmt != null) stmt.close();
			}
			catch(SQLException se2){ }// nothing we can do
		}
		return 0; // first available is 43 (and this is the answer to everything anyway...)
	}
	
	public static HashMap<String, Gene> getGenesInDBByEnsemblID(int organismID)
	{
		HashMap<String, Gene> genes = new HashMap<>();
		
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
		
			String sql = "SELECT * FROM genes WHERE organism_id="+organismID; // Restrict to genes of this organism
			ResultSet rs = stmt.executeQuery(sql);
			 
			// I check the results and compare to my gene list
			while(rs.next())
			{
				Gene g = new Gene();
				g.ensembl_id = rs.getString("ensembl_id");
				g.name = rs.getString("name");
				g.biotype = rs.getString("biotype");
				g.gene_length = rs.getInt("gene_length");
				g.sum_exon_length = rs.getInt("sum_exon_length");
				g.chr = rs.getString("chr");
				g.latest_ensembl_release = rs.getInt("latest_ensembl_release");
				
				g.ensembl_id = g.ensembl_id.toUpperCase();
				g.name = g.name.toUpperCase();
				
				// Dealing with alt_names
				String alt_names = rs.getString("alt_names");
				if(alt_names != null)
				{
					String[] tokens = alt_names.split(",");
					for(String gene:tokens)
					{
						String geneUp = gene.toUpperCase();
						g.alt_names.add(geneUp);
					}
				}
				
				// I add this gene in the HashMap. By ensembl_id ...
				if(genes.get(g.ensembl_id) != null) new ErrorJSON("Two entries in DB for " + g.ensembl_id + " for organism " + organismID);
				genes.put(g.ensembl_id, g);
			}
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
	
	public static void getGenesInDB(int organismID)
	{
		MapGene.init(); // We store the DB in this static object
		
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
		
			String sql = "SELECT * FROM genes WHERE organism_id="+organismID; // Restrict to genes of this organism
			ResultSet rs = stmt.executeQuery(sql);
			 
			// I check the results and compare to my gene list
			while(rs.next())
			{
				Gene g = new Gene();
				g.ensembl_id = rs.getString("ensembl_id");
				g.name = rs.getString("name");
				g.biotype = rs.getString("biotype");
				g.gene_length = rs.getInt("gene_length");
				g.sum_exon_length = rs.getInt("sum_exon_length");
				g.chr = rs.getString("chr");
				g.latest_ensembl_release = rs.getInt("latest_ensembl_release");
				
				// I add this gene in the HashMap. By ensembl_id ...
				ArrayList<Gene> gene_list = MapGene.ensembl_db.get(g.ensembl_id);
				if(gene_list == null) gene_list = new ArrayList<>();
				gene_list.add(g);
				MapGene.ensembl_db.put(g.ensembl_id, gene_list);
				
				// ... and once the gene name
				gene_list = MapGene.gene_db.get(g.name);
				if(gene_list == null) gene_list = new ArrayList<>();
				gene_list.add(g);
				MapGene.gene_db.put(g.name, gene_list);
					
				// ... and for every alt names
				String alt_names = rs.getString("alt_names");
				if(alt_names != null)
				{
					String[] tokens = alt_names.split(",");
					for(String gene:tokens)
					{
						g.alt_names.add(gene);
						gene_list = MapGene.alt_db.get(gene);
						if(gene_list == null) gene_list = new ArrayList<>();
						gene_list.add(g);
						MapGene.alt_db.put(gene, gene_list);
					}
				}
			}
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

	public static void connect()
	{
		try
		{
			Class.forName(JDBC_DRIVER);
			System.out.println("Connecting to database: " + URL);
			conn = DriverManager.getConnection(URL);
			System.out.println("Connected!");
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

	public static void disconnect()
	{
		try
		{
			if(conn!=null) conn.close();
		}
		catch(SQLException se)
		{
			se.printStackTrace();
		}
	}
}

