package parsing.model;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import bigarrays.StringArray64;
import model.Metadata;

public class MetaPreparse 
{
	public String name;
	public StringArray64 cellNames;
	public Gene[] geneNames;
	public HashMap<String, Integer> cellMap;
	public HashMap<String, Set<Integer>> geneMap;
	public HashSet<String> not_found;
	public HashSet<String> ambiguous;
	public HashMap<Integer, Metadata> metaIndex;
	public HashMap<Integer, Integer> metaIndexMatchParam;
	
	public MetaPreparse(String name) 
	{
		this.name = name;
	}
	
	public String toString() 
	{
		StringBuilder sb = new StringBuilder();
		sb.append("\"metadata\":\"").append(this.name).append("\"");
		if(cellNames != null) sb.append(",\"nber_cols\":").append(cellNames.size());
		if(geneNames != null) sb.append(",\"nber_rows\":").append(geneNames.length); 
		if(this.geneNames != null && this.geneNames.length != 0) 
		{ 
			sb.append(",\"genes\":[");
			String prefixT = "";
			int cnt = 0;
			for(Gene gene:this.geneNames) 
			{
				if(cnt == 10) break;
				sb.append(prefixT).append("\"").append((gene.ensembl_id == null)?gene.name:gene.ensembl_id).append("\"");
				prefixT = ",";
				cnt++;
			}
			sb.append("]");
		}
		if(this.cellNames != null && this.cellNames.size() != 0) 
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
		
		// Check empty metadatas
		sb.append(",\"empty_metadata\":[");
		String prefixT = "";
		for(Metadata me:metaIndex.values())
		{
			if(me.categories != null)
			{
				if(me.categories.size() == 1 && me.categories.iterator().next().equals(""))
				{
					sb.append(prefixT).append("\"").append(me.path).append("\"");
					prefixT = ",";
				}
			}
			else if(me.categoriesMap != null)
			{
				if(me.categoriesMap.size() == 1 && me.categoriesMap.containsKey(""))
				{
					sb.append(prefixT).append("\"").append(me.path).append("\"");
					prefixT = ",";
				}
			}
		}
		sb.append("]");
		
		// Now add all metadata
		sb.append(",").append(Metadata.toString(metaIndex.values()));
    	return sb.toString();
	}
	
	
}
