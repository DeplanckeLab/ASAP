package db;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;

import parsing.model.Gene;

public class EnsemblDB 
{
	public static FTPClient f = null;
	public static HashMap<Integer, HashMap<String, Gene>> genes = null; // By organism
		
	public static int lastRelease(String type)
	{
		switch(type)
		{
			case "vertebrates":
				return lastRelease("ftp.ensembl.org", "pub/");
			case "metazoa":
				return lastRelease("ftp.ensemblgenomes.org", "pub/metazoa/");
			case "bacteria":
				return lastRelease("ftp.ensemblgenomes.org", "pub/bacteria/");
			case "fungi":
				return lastRelease("ftp.ensemblgenomes.org", "pub/fungi/");
			case "plants":
				return lastRelease("ftp.ensemblgenomes.org", "pub/plants/");
			case "protists":
				return lastRelease("ftp.ensemblgenomes.org", "pub/protists/");
		}
		return 0;
	}
	
	public static int lastRelease(String db, String path)
	{
		try
		{
			f = new FTPClient();
			System.out.print("Connecting to "+db+"...");
			f.connect(db);	
			f.setKeepAlive(true);
			f.setControlKeepAliveTimeout(100000);
			f.setControlKeepAliveReplyTimeout(10000);
			f.setConnectTimeout(100000);
			f.setSoTimeout(100000);
		    f.setDefaultTimeout(100000);
		    f.setDataTimeout(100000);
		    f.enterLocalPassiveMode();
		    f.setFileType(FTP.BINARY_FILE_TYPE);
			f.login("anonymous","");
			System.out.println("Connected!");
			FTPFile[] files = f.listFiles(path);
			int maxRelease = 0;
			for(FTPFile fi:files)
			{
				if(fi.getName().startsWith("release-"))
				{
					String[] tokens = fi.getName().split("-");
					int r = Integer.parseInt(tokens[1]);
					if(r > maxRelease) maxRelease = r;
				}
			}
			System.out.println("Last release to date = " + maxRelease);
			f.disconnect();
			return maxRelease;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(-1);
		}
		return -1;
	}
	
	public static void updateDB()
	{
		DBManager.connect();
		String[] subdomains = new String[] {"vertebrates", "metazoa", "bacteria", "fungi", "plants", "protists"}; // Could also be accessed from DB
		for(String s:subdomains)
		{
			int lastInDB = DBManager.getLastRelease(s);
			int lastInEnsembl = lastRelease(s);
			genes = new HashMap<>(); // Reset for every subdomain
			for(int release = lastInDB; release <= lastInEnsembl; release++)
			{
				HashMap<String,String> urls = getGTFURL(release, s);
				if(urls != null) 
				{
					System.out.println(urls.size() + " URLs found for release " + release);
					addReleaseToDB(urls, release, s);
				} else System.out.println("NO URLs found for release " + release);
			}
			// Creating JSON
			for(int organism_id:genes.keySet())
			{
				String organism_name = DBManager.getOrganismName(organism_id);
				System.out.println("Creating JSON for organism " + organism_name + "....");
				Gene.toJSON(genes.get(organism_id), organism_id, organism_name, lastInEnsembl, s);
			}
			// Now updating the DB
			for(int organism_id:genes.keySet())
			{
				System.out.println("Updating organism " + organism_id + " in Database.");
				HashMap<String, Gene> db = genes.get(organism_id);
				for(String ens:db.keySet())
				{
					Gene g = db.get(ens);
					DBManager.insertOrUpdateGene(g.ensembl_id, g.name, g.biotype, g.chr, g.gene_length, g.sum_exon_length, organism_id, g.alt_names, g.latest_ensembl_release);
				}
				System.out.println("Organism " + organism_id + " was updated : " + DBManager.getNbGenesForOrganism(organism_id) + " genes now in database.");
				DBManager.setLastRelease(organism_id, lastInEnsembl);
			}
			DBManager.setLastRelease(s, lastInEnsembl);
		}
		DBManager.disconnect();
	}
	
	public static void addReleaseToDB(HashMap<String,String> urls, int release, String type)
	{
		if(urls.isEmpty()) return;
		int c = 0;
		for(String s:urls.keySet())
		{
			c++;
			int subdomainID = DBManager.getSubdomainID(type);
			System.out.println("Subdomain " + type + " , ID = " + subdomainID);
			int organismId = DBManager.getOrCreateOrganismId(s, subdomainID);
			System.out.println("Organism " + s + " , ID = " + organismId);
			int last = DBManager.getLastRelease(organismId);
			HashMap<String, Gene> genesForOrg = genes.get(organismId);
			if(genesForOrg == null) genesForOrg = DBManager.getGenesInDBByEnsemblID(organismId);// New entry for this organism
			System.out.println("Last release in DB is " + last);
			if(last >= release) System.out.println("Nothing to do. DB is already up-to-date for this species " + s);
			else
			{
				System.out.println("Updating "+ c + "/" + urls.size() + " : " + genesForOrg.size() +" genes for species " + s + " to release " + release);
				String url = urls.get(s);
				if(!url.equals("NA"))
				{
					HashMap<String, Gene> geneInfo = downloadEnsembl(type, url);
					for(String gene_id:geneInfo.keySet())
					{
						Gene g = geneInfo.get(gene_id);
						long l = g.getLength();
						if(l != 0) g.sum_exon_length = g.getLength();
						if(g.chr.equals("dmel_mitochondrion_genome")) g.chr = "MT";
						g.latest_ensembl_release = release; // It was updated to last release
						updateGenes(genesForOrg, g);
						//DBManager.insertOrUpdateGene(gene_id, g.gene_name, g.biotype, g.chr, (g.end - g.start + 1), g.sumExonLength, organismId, g.alternate_names, release);
					}
					genes.put(organismId, genesForOrg);
					System.out.println("Species " + s + " was updated (" + c + "/" + urls.size() + ") to release " + release + " : " + genesForOrg.size() + " genes now in database.");
				}
				else System.out.println("Release " + release + " does not exist for species " + s + ". No change in DB.");
				//DBManager.setLastRelease(organismId, release);
				genes.put(organismId, genesForOrg);
			}
		}
	}
	
	private static void updateGenes(HashMap<String, Gene> genes, Gene g)
	{
		Gene gFound = genes.get(g.ensembl_id); // The old version
		if(gFound != null) // If EnsemblID already in DB, I update it
		{
			// First I merge alt names
			for(String s:gFound.alt_names) g.alt_names.add(s);
			// Then I check if the "name" entry changed
			if(gFound.name!= null && !gFound.name.equals(g.name)) g.alt_names.add(gFound.name);
		}
		genes.put(g.ensembl_id, g);
	}
	
	public static HashMap<String, Gene> downloadEnsembl(String type, String path)
	{
		try
		{
			switch(type)
			{
				case "vertebrates":
					return downloadEnsembl(new URL("ftp://ftp.ensembl.org"+path));
				case "metazoa":
				case "bacteria":
				case "fungi":
				case "plants":
				case "protists":
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(-1);
		}
		return new HashMap<>();
	}
	
	public static HashMap<String, Gene> downloadEnsembl(URL url)
	{
		HashMap<String, Gene> geneInfo = new HashMap<>();
		try 
		{
			boolean foundGeneAnnotation = false;
			
			InputStream is = url.openStream();
	        InputStream gzipStream = new GZIPInputStream(is);
	        BufferedReader br = new BufferedReader(new InputStreamReader(gzipStream)); 
		    
	    	String line = br.readLine();
	    	int l = 1;
	    	while(line != null)
	    	{
	    		if(!line.startsWith("#"))
    			{
    				String[] tokens = line.split("\t");
    				// Parse Line
 					String[] params = tokens[8].split(";");
 					Gene g = new Gene();
					g.start = Long.parseLong(tokens[3]);
					g.end = Long.parseLong(tokens[4]);
					g.chr = tokens[0];
					String type = tokens[2];
					for(String param:params) 
					{
						try
						{
							String value = param.substring(param.indexOf("\"")+1, param.lastIndexOf("\""));
							if(param.contains("gene_name")) g.name = value.toUpperCase();
							if(param.contains("gene_id")) g.ensembl_id = value.toUpperCase();
							if(param.contains("gene_biotype")) g.biotype = value;
						}
						catch(StringIndexOutOfBoundsException oobe)
						{
							System.err.println("Error l." + l + " : parameter '" + param + "' has no associated values in double quotes.");
						}
					}
					if(g.biotype == null) g.biotype = tokens[1]; // earlier versions have biotype in the second column
					if(g.ensembl_id.startsWith("ENS") || g.ensembl_id.startsWith("FB")) // Weird things are not computed
					{
						if(type.equals("gene")) foundGeneAnnotation = true; // Tag
						
						// Is it already existing in the database?
						Gene found = geneInfo.get(g.ensembl_id);
						if(found == null) // Was not found before
						{
							if(type.equals("exon")) g.exon_id.add(g.start, g.end);
							geneInfo.put(g.ensembl_id, g);
						}
						else // Found
						{
							if(type.equals("exon")) found.exon_id.add(g.start, g.end);
							else if(type.equals("gene"))
							{
								found.start = g.start;
								found.end = g.end;
							}
						}
					}
    			}
    			line = br.readLine(); l++;
    		}
    		br.close();
    		if(!foundGeneAnnotation) // If not "gene" (old release)
    		{
				for(String gene_id:geneInfo.keySet())
				{
					Gene g = geneInfo.get(gene_id);
					if(g.exon_id.size() != 0) g.setGeneLengthToExons(); // If this gene is still annotated
				}
    		}
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
		return geneInfo;
	}
	
	private static HashMap<String, String> getGTFURL(int release, String type)
	{
		HashMap<String, String> res = null;
		try
		{
			f = new FTPClient();
			System.out.println("Connecting to ftp.ensembl.org...");
			f.connect("ftp.ensembl.org");	
			f.setKeepAlive(true);
			f.setControlKeepAliveTimeout(100000);
			f.setControlKeepAliveReplyTimeout(10000);
			f.setConnectTimeout(100000);
			f.setSoTimeout(100000);
		    f.setDefaultTimeout(100000);
		    f.setDataTimeout(100000);
		    f.enterLocalPassiveMode();
		    f.setFileType(FTP.BINARY_FILE_TYPE);
			f.login("anonymous","");
			System.out.println("Connected!");
			res = getGTF(release, type);
			f.disconnect();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(-1);
		}
		return res;
	}
	
	private static HashMap<String, String> getGTF(int release, String type) throws Exception // Some disappears from Ensembl => VectorBase (for e.g. some insects)
	{
		switch(type)
		{
			case "vertebrates":
				if(release < 47) return getGTF_43_46(release);
				if(release == 47) return getGTF_47(release);
				return getGTF_48_Max(release);
			case "metazoa": //TODO
			case "bacteria": //TODO
			case "fungi": //TODO
			case "plants": //TODO
			case "protists": //TODO
		}
		return null;
	}
	
	private static HashMap<String, String> getGTF_48_Max(int release) throws Exception // 48-Max // Jusqu'a r80 inclus, un seul fichier gtf
	{
		HashMap<String, String> urls = new HashMap<>();
		FTPFile[] files = f.listDirectories("/pub/release-"+release+"/gtf/");
		for(FTPFile fi:files)
		{
			String url = "/pub/release-"+release+"/gtf/"+fi.getName()+"/";
			FTPFile[] subfiles = f.listFiles(url);
			for(FTPFile fo:subfiles)
			{
				if(fo.getName().endsWith("gtf.gz") && fo.getName().indexOf("abinitio") == -1 && fo.getName().indexOf(".chr") == -1)
				{
					urls.put(fi.getName(), url+fo.getName());
				}
			}
		}
		return urls;
	}
	
	private static HashMap<String, String> getGTF_47(int release) throws Exception // 47
	{
		HashMap<String, String> urls = new HashMap<>();
		String url = "/pub/release-"+release+"/gtf/";
		FTPFile[] files = f.listFiles(url);
		for(FTPFile fi:files)
		{
			String name = fi.getName();
			name = name.substring(0, name.indexOf(".")).toLowerCase();
			urls.put(name, url+fi.getName());
		}
		return urls;
	}
	
	
	private static HashMap<String, String> getGTF_43_46(int release) throws Exception // 43-46
	{
		HashMap<String, String> urls = new HashMap<>();
		FTPFile[] files = f.listDirectories("/pub/release-"+release);
		for(FTPFile fi:files)
		{
			String url = "/pub/release-"+release+"/"+fi.getName()+"/data/gtf/";
			FTPFile[] subfiles = f.listFiles(url);
			for(FTPFile fo:subfiles)
			{
				String name = fi.getName();
				name = name.substring(0, name.indexOf("_", name.indexOf("_") + 1)).toLowerCase();
				urls.put(name, url+fo.getName());
			}
		}
		return urls;
	}
}
