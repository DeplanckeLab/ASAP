package db;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GeneAtlas // ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/atlas-latest-data.tar.gz
{ // http://www.proteinatlas.org/about/download
	public static void enrichrToGMT(String geneAtlasFile, String outputGMTFile) throws IOException
	{
		BufferedReader br = new BufferedReader(new FileReader(geneAtlasFile));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputGMTFile));
		
        String line = null;
        while ((line=br.readLine())!=null) 
        {
        	String[] results = line.split("\t");
        	bw.write(results[0] + "\tnull\tnull");
        	for(int i = 2; i<results.length; i++) bw.write("\t" + results[i].split(",")[0].toUpperCase().trim());
        	bw.write("\n");
        }
		bw.close();
		br.close();
	}
}