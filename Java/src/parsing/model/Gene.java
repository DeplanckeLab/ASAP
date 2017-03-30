package parsing.model;

public class Gene 
{
	public String ensembl_id;
	public String name;
	public String biotype;
	public long sum_exon_length;
	public long gene_length;
	public String chr;
	
	@Override
	public String toString() 
	{
		return ensembl_id + "\t" + name;
	}
}
