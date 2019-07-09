package parsing.model;

import java.util.HashSet;
import java.util.Set;

public class GroupPreparse 
{
	public String name;
	public long nbCells;
	public long nbGenes;
	public String[] cellNames;
	public String[] geneNames;
	public float[][] matrix;
	public boolean isCount = true;
	public Set<String> additionalMetadataPath = null;
	
	public GroupPreparse(String name) 
	{
		this.name = name;
		this.additionalMetadataPath = new HashSet<String>();
	}
}
