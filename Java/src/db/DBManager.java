package db;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import com.google.gson.Gson;

import json.ErrorJSON;
import model.GeneSet;
import model.MapGene;
import parsing.model.Gene;

public class DBManager 
{
	public static String JDBC_DRIVER = null;
	public static String URL = null;
	public static Connection conn = null;
	
	public static HashMap<Integer, Integer> listAllOrganisms() 
	{
		HashMap<Integer, Integer> taxToIdMap = new HashMap<>();
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "SELECT id,tax_id FROM organisms WHERE tax_id IS NOT NULL";
			
			ResultSet rs = stmt.executeQuery(sql);
			
			// List the results
			while(rs.next())
			{
				int taxId = rs.getInt("tax_id");
				int id = rs.getInt("id");
				taxToIdMap.put(taxId, id);
			}

			stmt.close();

			return taxToIdMap;
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
		return null;
	}
	
	
	public static long getOrganismFromGeneSets(long geneset_id) 
	{
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			
			String sql = "SELECT organism_id FROM gene_sets WHERE id="+geneset_id;
			ResultSet rs = stmt.executeQuery(sql);
			int res = -1;
			boolean found = rs.next();
			if(found) 
			{
				res = rs.getInt("organism_id");
			}
			stmt.close();
			if(!found) new ErrorJSON("This geneset with id=" + geneset_id + " is not in DB");
			if(res == 0 || res == -1) new ErrorJSON("No organism available for geneset_id id=" + geneset_id);
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
	
	public static void ListGeneSets(int organism_id) 
	{
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "SELECT id FROM gene_sets WHERE organism_id="+organism_id;
			ResultSet rs = stmt.executeQuery(sql);
			while(rs.next()) System.out.println(rs.getInt("id"));
			stmt.close();
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
	
	public static GeneSet getUniqueGeneSet(long id)
	{
		GeneSet res = null;
		Statement stmt = null;
		try
		{	
			stmt = conn.createStatement();
			String sql = "SELECT name,identifier,content,gene_set_id FROM gene_set_items WHERE id="+id;
			ResultSet rs = stmt.executeQuery(sql);
			
			// I go through the results
			while(rs.next())
			{
				res = new GeneSet();
				res.name = rs.getString("name");
				res.identifier = rs.getString("identifier");
				res.gene_set_id = rs.getInt("gene_set_id");
				String[] ids = rs.getString("content").split(",");
				res.content = new HashSet<Long>();
				for(int i = 0;i < ids.length; i++) if(!ids[i].equals("")) res.content.add(Long.parseLong(ids[i])); // They all should be Long
			}
			
			stmt.close();
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
		return res;
	}
	
	public static ArrayList<GeneSet> getGeneSets(long geneset_id) 
	{
		ArrayList<GeneSet> res = new ArrayList<GeneSet>();
		Statement stmt = null;
		try
		{	
			stmt = conn.createStatement();
			String sql = "SELECT identifier,name,content FROM gene_set_items WHERE gene_set_id="+geneset_id;
			ResultSet rs = stmt.executeQuery(sql);
			
			// I go through the results
			while(rs.next())
			{
				GeneSet g = new GeneSet();
				g.name = rs.getString("name");
				g.identifier = rs.getString("identifier");
				
				String[] ids = rs.getString("content").split(",");
				g.content = new HashSet<Long>();
				for(int i = 0;i < ids.length; i++) if(!ids[i].equals("")) g.content.add(Long.parseLong(ids[i])); // They all should be Long

				res.add(g);
			}
			
			stmt.close();
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
		return res;
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
	
	public static HashMap<String, Long> getGenesInDBByID(long organismID)
	{
		HashMap<String, Long> genes = new HashMap<String, Long>();
		
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
		
			String sql = "SELECT id,ensembl_id FROM genes WHERE organism_id="+organismID; // Restrict to genes of this organism
			ResultSet rs = stmt.executeQuery(sql);
			 
			// I check the results and compare to my gene list
			while(rs.next())
			{
				genes.put(rs.getString("ensembl_id"), rs.getLong("id"));
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
				if(g.ensembl_id != null)
				{
					ArrayList<Gene> gene_list = MapGene.ensembl_db.get(g.ensembl_id.toUpperCase());
					if(gene_list == null) gene_list = new ArrayList<>();
					gene_list.add(g);
					MapGene.ensembl_db.put(g.ensembl_id.toUpperCase(), gene_list);
				}
				
				// ... and once the gene name
				if(g.name != null)
				{
					ArrayList<Gene> gene_list = MapGene.gene_db.get(g.name.toUpperCase());
					if(gene_list == null) gene_list = new ArrayList<>();
					gene_list.add(g);
					MapGene.gene_db.put(g.name.toUpperCase(), gene_list);
				}
				
				// ... and for every alt names
				String alt_names = rs.getString("alt_names");
				if(alt_names != null)
				{
					String[] tokens = alt_names.split(",");
					for(String gene:tokens)
					{
						g.alt_names.add(gene);
						ArrayList<Gene> gene_list = MapGene.alt_db.get(gene.toUpperCase());
						if(gene_list == null) gene_list = new ArrayList<>();
						gene_list.add(g);
						MapGene.alt_db.put(gene.toUpperCase(), gene_list);
					}
				}
				
				// ... and for every obsolete alt names
				String obsolete_alt_names = rs.getString("obsolete_alt_names");
				if(obsolete_alt_names != null)
				{
					String[] tokens = obsolete_alt_names.split(",");
					for(String gene:tokens)
					{
						ArrayList<Gene> gene_list = MapGene.obsolete_db.get(gene.toUpperCase());
						if(gene_list == null) gene_list = new ArrayList<>();
						gene_list.add(g);
						MapGene.obsolete_db.put(gene.toUpperCase(), gene_list);
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
	
	public static HashMap<String, Integer> getListCatJSON(long id)
	{
		HashMap<String, Integer> categories = new HashMap<String, Integer>();
		
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
		
			String sql = "SELECT list_cat_json FROM annots WHERE id="+id;
			ResultSet rs = stmt.executeQuery(sql);
			 
			// I check the results and compare to my gene list
			if (!rs.next()) new ErrorJSON("No list_cat_json in annots table with id = " + id);
			else
			{
				String json = rs.getString("list_cat_json");
				Gson gson = new Gson();
				String[] cat = gson.fromJson(json, String[].class);
				for(int i = 0; i < cat.length; i++) categories.put(cat[i], i + 1);
			    if(rs.next()) new ErrorJSON("Too many list_cat_json in annots table with id = " + id);
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
		return categories;
	}
	
	public static int getNbCatJSON(String id)
	{
		int nbCat = -1;
		
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
		
			String sql = "SELECT list_cat_json FROM annots WHERE id="+id;
			ResultSet rs = stmt.executeQuery(sql);
			 
			// I check the results and compare to my gene list
			if (!rs.next()) new ErrorJSON("No list_cat_json in annots table with id = " + id);
			else
			{
				String json = rs.getString("list_cat_json");
				Gson gson = new Gson();
				String[] cat = gson.fromJson(json, String[].class);
				nbCat = cat.length;
			    if(rs.next()) new ErrorJSON("Too many list_cat_json in annots table with id = " + id);
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
		return nbCat;
	}

	public static void connect()
	{
		try
		{
			Class.forName(JDBC_DRIVER);
			conn = DriverManager.getConnection(URL);
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

