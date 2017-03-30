package model;

import java.util.HashSet;

import org.apache.sis.measure.Range;
import org.apache.sis.util.collection.RangeSet;

public class GeneInfo
{
	public String gene_id; // Main ID
	public String gene_name;
	public HashSet<String> alternate_names = new HashSet<String>();
	public String biotype = null;
	public long start;
	public long end;
	public String chr;
	public long sumExonLength = 0;
	public RangeSet<Long> exon_id = RangeSet.create(Long.class, true, true);
	
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
	
}

