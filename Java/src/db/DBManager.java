package db;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;

import model.Parameters;
import parsing.model.Gene;
import tools.Utils;

public class DBManager 
{
	public static String JDBC_DRIVER = null;
	public static String URL = null;
	public static Connection conn = null;
	
	public static final String[] species = new String[] {"ailuropoda_melanoleuca","anas_platyrhynchos","anolis_carolinensis","astyanax_mexicanus","bos_taurus","caenorhabditis_elegans","callithrix_jacchus","canis_familiaris","cavia_porcellus","chlorocebus_sabaeus","choloepus_hoffmanni","ciona_intestinalis","ciona_savignyi","danio_rerio","dasypus_novemcinctus","dipodomys_ordii","drosophila_melanogaster","echinops_telfairi","equus_caballus","erinaceus_europaeus","felis_catus","ficedula_albicollis","gadus_morhua","gallus_gallus","gasterosteus_aculeatus","gorilla_gorilla","homo_sapiens","ictidomys_tridecemlineatus","latimeria_chalumnae","lepisosteus_oculatus","loxodonta_africana","macaca_mulatta","macropus_eugenii","meleagris_gallopavo","microcebus_murinus","monodelphis_domestica","mus_musculus","mustela_putorius_furo","myotis_lucifugus","nomascus_leucogenys","ochotona_princeps","oreochromis_niloticus","ornithorhynchus_anatinus","oryctolagus_cuniculus","oryzias_latipes","otolemur_garnettii","ovis_aries","pan_troglodytes","papio_anubis","pelodiscus_sinensis","petromyzon_marinus","poecilia_formosa","pongo_abelii","procavia_capensis","pteropus_vampyrus","rattus_norvegicus","saccharomyces_cerevisiae","sarcophilus_harrisii","sorex_araneus","sus_scrofa","taeniopygia_guttata","takifugu_rubripes","tarsius_syrichta","tetraodon_nigroviridis","tupaia_belangeri","tursiops_truncatus","vicugna_pacos","xenopus_tropicalis","xiphophorus_maculatus"};
	
	public static void generateEnrichmentDB()
	{
		try 
		{
			DrugBankParser.parse("hsa", "C:/Users/gardeux/Dropbox/GeneSets/DrugBank/drugbank.xml", "C:/Users/gardeux/Dropbox/ASAP/Scripts/Enrichment/Genesets/drugbank_hsa.gmt");
			DrugBankParser.parse("mmu", "C:/Users/gardeux/Dropbox/GeneSets/DrugBank/drugbank.xml", "C:/Users/gardeux/Dropbox/ASAP/Scripts/Enrichment/Genesets/drugbank_mmu.gmt");
			KEGGRestApi.generateKEGGDB("hsa", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/kegg_hsa.gmt"); // mmu, rno (rat), dme (droso mela)
			KEGGRestApi.generateKEGGDB("mmu", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/kegg_mmu.gmt"); // mmu, rno (rat), dme (droso mela)
			GeneAtlas.enrichrToGMT("C:/Users/Vincent/Desktop/GeneSets/GeneAtlas/Human_Gene_Atlas.txt", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/gene_atlas.hsa.gmt");
			GeneAtlas.enrichrToGMT("C:/Users/Vincent/Desktop/GeneSets/GeneAtlas/Mouse_Gene_Atlas.txt", "C://Users/Vincent/Dropbox/ASAP/Scripts/Enrichment/Genesets/gene_atlas.mmu.gmt");
			//GODatabase.generateGODB("C:/Users/gardeux/Dropbox/ASAP/Scripts/Enrichment/Genesets/", "hsa");
			//GODatabase.generateGODB("C:/Users/gardeux/Dropbox/ASAP/Scripts/Enrichment/Genesets/", "mmu");
		} 
		catch (IOException e)
		{
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void createDB()
	{
		long t = System.currentTimeMillis();
		connect();
		System.out.print("Running: DROP TABLE genes...");
		dropTable("genes");
		System.out.println("DONE!");
		
		System.out.print("Running: CREATE TABLE genes...");
		createGeneTable();
		System.out.println("DONE!");
		disconnect();
		System.out.println("Creating DB took " + Utils.toReadableTime(System.currentTimeMillis() - t));
	}
	
	public static void createGeneTable()
	{
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "CREATE TABLE genes " +
					"(id INTEGER NOT NULL AUTO_INCREMENT, " +
					" ensembl_id TEXT, " + 
					" name TEXT, " + 
					" alternate_names TEXT, " +
					" organism_id INTEGER DEFAULT NULL, " +
					" biotype TEXT, " +
					" sum_exon_length INTEGER, " +
					" gene_length INTEGER, " +
					" created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, " +
					" PRIMARY KEY ( id ))";
			stmt.executeUpdate(sql);
			sql = "CREATE TABLE gene_names " +
					"(id INTEGER NOT NULL AUTO_INCREMENT, " +
					" gene_id int references genes, " + 
					" value TEXT, " + 
					" PRIMARY KEY ( id ))";
			stmt.executeUpdate(sql);
		}
		catch(Exception e) { e.printStackTrace(); }
		finally
		{
			try
			{
				if(stmt!=null) stmt.close();
			}
			catch(SQLException se2){ }// nothing we can do
		}
	}
	
	public static void dropTable(String name)
	{
		Statement stmt = null;
		try
		{
			  stmt = conn.createStatement();
			  stmt.executeUpdate("DROP TABLE " + name);
		}
		catch(Exception e) { e.printStackTrace(); }
		finally
		{
			try
			{
				if(stmt!=null) stmt.close();
			}
			catch(SQLException se2){ }// nothing we can do
		}
	}
	
	public static void insertGene(String ensembl_id, String name, String alternateNames, int organism_id, String biotype, int exon_length, int full_length)
	{
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
			String sql = "INSERT INTO genes (ensembl_id, name, alternate_names, organism_id, biotype, exon_length, full_length) VALUES ('"+ensembl_id+"','"+name+"','"+alternateNames+"',"+organism_id+",'"+biotype+"',"+exon_length+","+full_length+")";
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

	public static HashMap<String, ArrayList<Gene>> getGenesInDB()
	{
		HashMap<String, ArrayList<Gene>> genes = new HashMap<>(); // Hopefully it is quick enough to build?? If not it should not be stored.
		connect();
		
		// I perform my request
		Statement stmt = null;
		try
		{
			stmt = conn.createStatement();
		
			String sql = "SELECT * FROM genes WHERE organism_id="+Parameters.organism; // Restrict to genes of this organism
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
				String ensUp = g.ensembl_id.toUpperCase();
				String nameUp = g.name.toUpperCase();
				// I add this gene twice in the HashMap. Once its ensembl_id ...
				ArrayList<Gene> gene_list = genes.get(ensUp);
				if(gene_list == null) gene_list = new ArrayList<>();
				gene_list.add(g);
				genes.put(ensUp, gene_list);
				// ... and once the gene name
				gene_list = genes.get(nameUp);
				if(gene_list == null) gene_list = new ArrayList<>();
				gene_list.add(g);
				genes.put(nameUp, gene_list);
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
		disconnect();
		return genes;
	}

	public static void connect()
	{
		try
		{
			Class.forName(JDBC_DRIVER);
			System.out.print("Connecting to database...");
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

