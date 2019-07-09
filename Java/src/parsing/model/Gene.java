package parsing.model;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.sis.measure.Range;
import org.apache.sis.util.collection.RangeSet;

import json.ErrorJSON;
import model.Parameters;

public class Gene 
{
	public String ensembl_id;
	public String name;
	public String biotype;
	public long sum_exon_length;
	public long gene_length;
	public long start;
	public long end;
	public String chr;
	public HashSet<String> alt_names;
	public int latest_ensembl_release;
	public RangeSet<Long> exon_id = RangeSet.create(Long.class, true, true);
	
	public Gene()
	{
		alt_names = new HashSet<>();
	}
	
	@Override
	public String toString() 
	{
		return ensembl_id + "\t" + name;
	}
	

	public long getLength()
	{
		long sum = 0;
		for(Range<Long> r:exon_id) sum += ((Long)r.getMaxValue() - (Long)r.getMinValue() + 1);
		return sum;
	}
	
	public void setGeneLengthToExons()
	{
		this.start = Long.MAX_VALUE;
		this.end = Long.MIN_VALUE;
		for(Range<Long> r:exon_id) 
		{
			if((Long)r.getMaxValue() > this.end) this.end = (Long)r.getMaxValue();
			if((Long)r.getMinValue() < this.start) this.start = (Long)r.getMinValue();
		}
	}
	
	public static String buildAltNamesString(HashSet<String> names)
	{
		if(names == null) return "";
		StringBuilder sb = new StringBuilder();
		String prefix = "";
		for(String n:names) 
		{
			sb.append(prefix).append(n);
			prefix = ",";
		}
		return sb.toString();
	}
	
	public static void toJSON(HashMap<String, Gene> genes, int organism_id, String organism_name, int release, String subdomain)
	{
    	try
    	{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output."+organism_id+".json"));
    		bw.write("{\"subdomain\":\"" + subdomain + "\",");
    		bw.write("\"organism_id\":\"" + organism_id + "\",");
    		bw.write("\"organism_name\":\"" + organism_name + "\",");
    		bw.write("\"nber_genes\":" + genes.size() + ",");
    		bw.write("\"release\":" + release + ",");
    		bw.write("\"header\":[\"ensembl_id\",\"name\",\"biotype\",\"chr\",\"gene_length\",\"sum_exon_length\",\"alt_names\",\"latest_ensembl_release\"],");
    		bw.write("\"genes\":");
    		StringBuilder sb = new StringBuilder("[");
    		String prefix = "";
			for(String ens:genes.keySet())
			{
				Gene g = genes.get(ens);
				sb.append(prefix).append("[\"").append(g.ensembl_id).append("\"");
				sb.append(",\"").append(g.name).append("\"");
				sb.append(",\"").append(g.biotype).append("\"");
				sb.append(",\"").append(g.chr).append("\"");
				sb.append(",").append(g.gene_length);
				sb.append(",").append(g.sum_exon_length);
				sb.append(",\"").append(buildAltNamesString(g.alt_names)).append("\"");
				sb.append(",").append(g.latest_ensembl_release);
				sb.append("]");
				prefix = ",";
			}
        	bw.write(sb.append("]").toString());
        	bw.write("}");
        	bw.close();
    	}
    	catch(IOException ioe)
    	{
    		new ErrorJSON(ioe.getMessage());
    	}
	}
}
