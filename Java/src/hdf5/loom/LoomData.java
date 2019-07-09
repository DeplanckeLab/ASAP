package hdf5.loom;

import java.util.ArrayList;
import java.util.HashSet;

import bigarrays.DoubleArray64;
import bigarrays.FloatArray64;
import bigarrays.IntArray64;
import bigarrays.LongArray64;
import bigarrays.StringArray64;
import model.ERCC;
import model.Metadata;

public class LoomData 
{
	public ArrayList<Metadata> meta = null;
	public HashSet<Long> removed = null; // If some rows of the main matrix were removed
	
	// Summary infos
	public long nber_cells = 0;
	public long nber_genes = 0;
	public long nber_not_found_genes = 0;
	public long nber_ercc = 0;
	public boolean is_count_table = true;
	public long nber_zeros = 0;
	
	// Main dataset
	public ERCC erccs;
	public StringArray64 gene_names = null; // Not initialized by default
	public StringArray64 ens_names = null; // Not initialized by default
	public StringArray64 cell_names;
	
	// Gene annotations
	public LongArray64 sumExonLength; // Length in bp of exonic portion of the corresponding gene
	public DoubleArray64 sum; // Sum expression for corresponding gene
	public StringArray64 biotypes;
	public StringArray64 chromosomes;
	
	// Cell annotations
	public IntArray64 detected_genes; // Number of genes with >0 reads
	public DoubleArray64 depth; // Sum expression for corresponding cell
	public FloatArray64 ribosomalContent; // Number of read per biotype per cell
	public FloatArray64 proteinCodingContent; // Number of read per biotype per cell
	public FloatArray64 mitochondrialContent; // Number of read per chromosome per cell
	
	public IntArray64 __no_feature = null; // Not initialized
	public IntArray64 __ambiguous = null; // Not initialized
	public IntArray64 __too_low_aQual = null; // Not initialized
	public IntArray64 __not_aligned = null; // Not initialized
	public IntArray64 __alignment_not_unique = null; // Not initialized
	
	public LoomData(long nbGenes, long nbCells) 
	{
		this.removed = new HashSet<Long>();
		this.meta = new ArrayList<Metadata>();	
		this.nber_genes = nbGenes;
		this.nber_cells = nbCells;

		// Gene annotation
		this.sumExonLength = new LongArray64(nbGenes);
		this.sum = new DoubleArray64(nbGenes);
		this.biotypes = new StringArray64(nbGenes);
		this.chromosomes = new StringArray64(nbGenes);
		
		// Cell annotations
		this.cell_names = new StringArray64(nbCells);
		this.depth = new DoubleArray64(nbCells);
		this.detected_genes = new IntArray64(nbCells);
		this.ribosomalContent = new FloatArray64(nbCells);
		this.mitochondrialContent = new FloatArray64(nbCells);
		this.proteinCodingContent = new FloatArray64(nbCells);
	}
}
