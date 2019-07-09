package parsing.model;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import bigarrays.StringArray64;
import model.Metadata;

public class MetaPreparse 
{
	public String name;
	public StringArray64 cellNames;
	public Gene[] geneNames;
	public HashMap<String, Integer> cellMap;
	public HashMap<String, List<Integer>> geneMap;
	public HashSet<String> not_found;
	public HashSet<String> ambiguous;
	public HashMap<Integer, Metadata> metaIndex;
	
	public MetaPreparse(String name) 
	{
		this.name = name;
	}
	
	/*
	@override
	public String toString() 
	{
		StringBuilder sb = new StringBuilder();
		sb.append("\"metadata\":\"").append(this.name).append("\"");
		if(cellNames != null) sb.append(",\"nber_cols\":").append(cellNames.length);
		if(geneNames != null) sb.append(",\"nber_rows\":").append(geneNames.length); 
		if(this.geneNames != null && this.geneNames.length != 0) 
		{ 
			sb.append(",\"genes\":[");
			String prefixT = "";
			for(Gene gene:this.geneNames) 
			{
				sb.append(prefixT).append("\"").append((gene.ensembl_id == null)?gene.name:gene.ensembl_id).append("\"");
				prefixT = ",";
			}
			sb.append("]");
		}
		if(this.cellNames != null && this.cellNames.length != 0) 
		{ 
			sb.append(",\"cells\":[");
			String prefixT = "";
			int cnt = 0;
			for(String cell:this.cellNames) 
			{
				if(cnt == 10) break;
				sb.append(prefixT).append("\"").append(cell).append("\"");
				prefixT = ",";
				cnt++;
			}
			sb.append("]");
		}
		if(this.not_found != null && this.not_found.size() != 0) 
		{ 
			sb.append(",\"not_found\":[");
			String prefixT = "";
			for(String cell:this.not_found) 
			{
				sb.append(prefixT).append("\"").append(cell).append("\"");
				prefixT = ",";
			}
			sb.append("]");
		}
		if(this.ambiguous != null && this.ambiguous.size() != 0) 
		{ 
			sb.append(",\"ambiguous\":[");
			String prefixT = "";
			for(String gene:this.ambiguous) // Only genes can be ambiguous
			{
				sb.append(prefixT).append("\"").append(gene).append("\"");
				prefixT = ",";
			}
			sb.append("]");
		}
		sb.append(",").append(toString(metaIndex.values()));
    	return sb.toString();
	}*/
	
	
}
