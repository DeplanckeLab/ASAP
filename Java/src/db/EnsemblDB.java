package db;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import org.apache.sis.util.collection.RangeSet;

import model.GeneInfo;
import model.Parameters;

public class EnsemblDB 
{
	public static FTPClient f = null;
	public static HashMap<String, HashMap<Integer, String>> urls = new HashMap<String, HashMap<Integer, String>>();
	public static HashMap<String, GeneInfo> geneInfo = new HashMap<>();
	
	public static void main(String[] args) throws Exception // Run this once
	{
		generateSpeciesURLFromEnsembl();
	}
	
	public static void readSpecies() throws IOException
	{
		BufferedReader br = new BufferedReader(new FileReader("species.txt"));
		String line = br.readLine();
		HashMap<Integer, Integer> lineToRelease = new HashMap<>();
		String[] header = line.split("\t");
		for(int i = 1; i < header.length; i++) lineToRelease.put(i, Integer.parseInt(header[i]));
		line = br.readLine();
		while(line != null) // One line per specie
		{
			String[] tokens = line.split("\t");
			HashMap<Integer, String> urlss = new HashMap<Integer, String>();
			for(int i = 1; i < tokens.length; i++) urlss.put(lineToRelease.get(i), tokens[i]);
			urls.put(tokens[0], urlss);
			line = br.readLine();
		}
		br.close();
	}
	
	public static void createEnsemblDB()
	{
		try
		{
			readSpecies(); // Read our precomputed file
			
			System.out.println("Computing gene info for species " + Parameters.organism_S + " ...");
			HashMap<Integer, String> toLoad = urls.get(Parameters.organism_S);
			for(int i = 43; i <= 87; i++) // order is important
			{
				String url = toLoad.get(i);
				if(!url.equals("NA"))
				{
					downloadEnsembl(toLoad.get(i));
					for(String gene_id:geneInfo.keySet())
					{
						GeneInfo g = geneInfo.get(gene_id);
						long l = g.getLength();
						if(l != 0) g.sumExonLength = g.getLength();
					}
					System.out.println("Release " + i + ": "+geneInfo.size()+" genes now in database.");
				}
				else System.out.println("Release " + i + " does not exist for species " + Parameters.organism_S);
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + Parameters.organism_S + ".txt"));
			bw.write("Ensembl\tName\tAltNames\tBiotype\tGeneLength\tSumExonLength\tChr\n");
			ArrayList<String> tmp = new ArrayList<>();
			for(String gene_id:geneInfo.keySet()) tmp.add(gene_id);
			Collections.sort(tmp);
			for(String gene_id:tmp)
			{
				GeneInfo g = geneInfo.get(gene_id);
				if(g.chr.equals("dmel_mitochondrion_genome")) g.chr = "MT";
				bw.write(gene_id + "\t" + g.gene_name + "\t" + buildAltNamesString(g.alternate_names, g.gene_id.toUpperCase(), g.gene_name.toUpperCase()) + "\t" + g.biotype + "\t" + (g.end - g.start + 1) + "\t" + g.sumExonLength + "\t" + g.chr + "\n");
			}
			bw.close();
		}
		catch(IOException ioe)
		{
			ioe.printStackTrace();
		}
	}
	
	private static String buildAltNamesString(HashSet<String> names, String ensNameUp, String geneNameUp)
	{
		HashSet<String> unique = new HashSet<String>();
		String res = "";
		for(String n:names)
		{
			String nUp = n.toUpperCase();
			if(!unique.contains(nUp))
			{
				if(!nUp.equals(ensNameUp) && !nUp.equals(geneNameUp)) res += n + ","; // Do not take if doublon, or same ens, or same gene name
				unique.add(nUp);
			}
		}
		if(res.endsWith(",")) res = res.substring(0, res.length() - 1); // Remove last comma
		return res;
	}
	
	public static void downloadEnsembl(String path)
	{
		try 
		{
			//int i1 = path.indexOf("release-");
			//int release = Integer.parseInt(path.substring(i1+8, path.indexOf("/", i1)));
			boolean foundGeneAnnotation = false;
			URL url = new URL("ftp://ftp.ensembl.org"+path);
			InputStream is = url.openStream();
	        InputStream gzipStream = new GZIPInputStream(is);
	        BufferedReader br = new BufferedReader(new InputStreamReader(gzipStream)); 
		    
	        for(String key:geneInfo.keySet()) geneInfo.get(key).exon_id = RangeSet.create(Long.class, true, true); // This needs to be reset for every release
	        
	    	String line = br.readLine();
	    	while(line != null)
	    	{
	    		if(!line.startsWith("#"))
    			{
    				String[] tokens = line.split("\t");
    				// Parse Line
 					String[] params = tokens[8].split(";");
					long start = Long.parseLong(tokens[3]);
					long end = Long.parseLong(tokens[4]);
					String gene_name = null;
					String gene_id = null;
					String biotype = null;
					String chr = tokens[0];
					String type = tokens[2];
					for(String param:params) 
					{
						String value = param.substring(param.indexOf("\"")+1, param.lastIndexOf("\""));
						if(param.contains("gene_name")) gene_name = value;
						if(param.contains("gene_id")) gene_id = value;
						if(param.contains("gene_biotype")) biotype = value;
					}
					if(gene_name == null) gene_name = gene_id;
					if(biotype == null) biotype = tokens[1]; // earlier versions have biotype in the second column
					if(gene_id.startsWith("ENS") || gene_id.startsWith("FB")) // Weird things are not computed
					{
						// Is it already existing in the database?
						GeneInfo g = geneInfo.get(gene_id);
						if(g == null) 
						{
							g = new GeneInfo();
							g.gene_id = gene_id; // Will never change
							g.gene_name = gene_name;
						}
						else // If this EnsemblID already in the database
						{
							if(!gene_name.equals(g.gene_name)) // If new name
							{
								if(!g.gene_name.startsWith("ENS") && !g.alternate_names.contains(g.gene_name)) g.alternate_names.add(g.gene_name); // Add the old name
								g.gene_name = gene_name;
								if(g.alternate_names.contains(gene_name)) g.alternate_names.remove(gene_name); // Remove current gene from old names
							}
						}
						g.chr = chr;
						g.biotype = biotype; // Update (do not save previous)
						// Which type is it?
						if(type.equals("gene"))
						{
							foundGeneAnnotation = true;
							g.end = end; 
							g.start = start;
						}
						else if(type.equals("exon")) g.exon_id.add(start, end); // This needs to be flushed for each release
						geneInfo.put(g.gene_id, g);
					}
    			}
    			line = br.readLine();
    		}
    		br.close();
    		if(!foundGeneAnnotation) // If not "gene" (old release)
    		{
				for(String gene_id:geneInfo.keySet())
				{
					GeneInfo g = geneInfo.get(gene_id);
					if(g.exon_id.size() != 0) g.setGeneLengthToExons(); // If this gene is still annotated
				}
    		}
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}	
	}
	
	private static void generateSpeciesURLFromEnsembl() throws Exception
	{
		f = new FTPClient();
		f.connect("ftp.ensembl.org");
		f.login("anonymous","");
		for(int i = 43; i <= 87; i++) getGTF(i);
		ArrayList<String> species = new ArrayList<String>();
		for(String s:urls.keySet()) species.add(s);
		Collections.sort(species);
		BufferedWriter bw = new BufferedWriter(new FileWriter("species.txt"));
		bw.write("species");
		for(int i = 43; i <= 87; i++) bw.write("\t" + i);
		bw.write("\n");
		for(String s:species)
		{
			bw.write(s);
			HashMap<Integer, String> releases = urls.get(s);
			for(int i = 43; i <= 87; i++) 
			{
				String url = releases.get(new Integer(i));
				if(url == null) url = "NA";
				bw.write("\t"+url);
			}
			bw.write("\n");
		}
		bw.close();
	}
	
	private static void getGTF(int release) throws Exception // Some disappears from Ensembl => VectorBase (for e.g. some insects)
	{
		if(release >= 48) getGTF_48_87(release);
		else if(release == 47) getGTF_47(release);
		else if(release >= 43) getGTF_43_46(release);
		else System.exit(-1);
	}
	
	private static void getGTF_48_87(int release) throws Exception // 48-87 // Jusqu'a r80 inclus, un seul fichier gtf
	{

		FTPFile[] files = f.listDirectories("/pub/release-"+release+"/gtf/");
		for(FTPFile fi:files)
		{
			String url = "/pub/release-"+release+"/gtf/"+fi.getName()+"/";
			FTPFile[] subfiles = f.listFiles(url);
			for(FTPFile fo:subfiles)
			{
				if(fo.getName().endsWith("gtf.gz") && fo.getName().indexOf("abinitio") == -1 && fo.getName().indexOf(".chr") == -1)
				{
					HashMap<Integer, String> releases = urls.get(fi.getName());
					if(releases == null) releases = new HashMap<Integer, String>();
					releases.put(release, url+fo.getName());
					urls.put(fi.getName(), releases);
				}
			}
		}
	}
	
	private static void getGTF_47(int release) throws Exception // 47
	{
		String url = "/pub/release-"+release+"/gtf/";
		FTPFile[] files = f.listFiles(url);
		for(FTPFile fi:files)
		{
			String name = fi.getName();
			name= name.substring(0, name.indexOf(".")).toLowerCase();
			HashMap<Integer, String> releases = urls.get(name);
			if(releases == null) releases = new HashMap<Integer, String>();
			releases.put(release, url+fi.getName());
			urls.put(name, releases);
		}
	}
	
	
	private static void getGTF_43_46(int release) throws Exception // 43-46
	{
		FTPFile[] files = f.listDirectories("/pub/release-"+release);
		for(FTPFile fi:files)
		{
			String url = "/pub/release-"+release+"/"+fi.getName()+"/data/gtf/";
			FTPFile[] subfiles = f.listFiles(url);
			for(FTPFile fo:subfiles)
			{
				String name = fi.getName();
				name = name.substring(0, name.indexOf("_", name.indexOf("_") + 1)).toLowerCase();
				HashMap<Integer, String> releases = urls.get(name);
				if(releases == null) releases = new HashMap<Integer, String>();
				releases.put(release, url+fo.getName());
				urls.put(name, releases);
			}
		}
	}
}
