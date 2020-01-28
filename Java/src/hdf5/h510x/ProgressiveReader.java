package hdf5.h510x;

import bigarrays.LongArray64;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import json.ErrorJSON;
import json.ParsingJSON;
import model.Parameters;

public class ProgressiveReader 
{
	private IHDF5Reader reader = null;
	private int dataBlockSize = -1;
	private int indicesBlockSize = -1;
	
	LongArray64 indptr = null; // column indexes 'indptr', required for recreating the dense matrix
	
	Block dataBlocks = null;
	Block indicesBlocks = null;
	
	private int currentDataBlockRead = -1;
	private int currentIndicesBlockRead = -1;
	
	private long nbGenes = -1;
	private long nbCells = -1;
	private int blockSizeX = -1;
	public int nbTotalBlocks = -1;
	
	public ProgressiveReader(IHDF5Reader reader, long nbGenes, long nbCells, int blockSizeX, int blockSizeY) 
	{
		this.reader = reader;
		this.blockSizeX = blockSizeX;
		this.dataBlockSize = this.reader.getDataSetInformation("/" + Parameters.selection + "/data").tryGetChunkSizes()[0];
		this.indicesBlockSize = this.reader.getDataSetInformation("/" + Parameters.selection + "/indices").tryGetChunkSizes()[0];
		this.nbGenes = nbGenes;
		this.nbCells = nbCells;
		
		// Read the column indexes 'indptr', required for recreating the dense matrix
		this.indptr = readLong("/" + Parameters.selection + "/indptr"); // This one should have reasonable size (nb cells), so we load it totally
		
		// How many blocks to process?
		this.nbTotalBlocks = (int)Math.ceil((double)this.nbCells / blockSizeX);
		System.out.println("Writing & Parsing " + nbTotalBlocks + " independent blocks");
	}
	
	private LongArray64 readLong(String path)
	{
		if(this.reader == null) new ErrorJSON("Please open the HDF5 file first");
		long length = this.reader.getDataSetInformation(path).getDimensions()[0];
		LongArray64 res = new LongArray64(length); // Know the actual size of the array, and create it sufficiently big
		int nbChunks = (int)(length / LongArray64.chunkSize()) + 1;
		for(int i = 0; i < nbChunks; i++) res.set(i, this.reader.int64().readArrayBlock(path, LongArray64.chunkSize(), i));
		return res;
	}
	
	public float[][] readSubMatrix(long start, long end, ParsingJSON json) // From cols start to end && all rows
	{
		// Create the submatrix to return
		float[][] submatrix = new float[(int)this.nbGenes][blockSizeX];
			
		// Fill the dense submatrix with values != 0
		for(long j = start; j < end; j++) // j varies accross cells ( cols )
		{	
			for(long x = indptr.get(j); x < indptr.get(j+1); x++) // x varies accross genes ( rows )
			{
				float value = getData(x);
				int i = getIndices(x); // Original index
				
				// Put the value in the matrix
				submatrix[i][(int)(j - start)] = value; // i - start should always be in int range (because chunk size is int)
				
				// Process with annotations
				if(json.data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) json.data.is_count_table = false;
				
				// Handle biotype count per cell
				String biotype = json.data.biotypes.get(i);
				if(biotype.equals("protein_coding")) json.data.proteinCodingContent.set(j, json.data.proteinCodingContent.get(j) + value);
				else if(biotype.equals("rRNA")) json.data.ribosomalContent.set(j, json.data.ribosomalContent.get(j) + value);
				if(json.data.chromosomes.get(i).equals("MT")) json.data.mitochondrialContent.set(j, json.data.mitochondrialContent.get(j) + value);
    				
    			// Generate the sums
				json.data.depth.set(j, json.data.depth.get(j) + value);
				json.data.sum.set(i, json.data.sum.get(i) + value);
			}
		}
		
		return submatrix;
	}
	
	private void readNextDataBlock()
	{
		currentDataBlockRead++;
		if(dataBlocks == null) dataBlocks = new Block(0, dataBlockSize - 1);
		else dataBlocks = new Block(dataBlocks.end + 1, dataBlocks.end + dataBlockSize);
		dataBlocks.values = reader.int32().readArrayBlock("/" + Parameters.selection + "/data", dataBlockSize, currentDataBlockRead);
	}
	
	private int getData(long index)
	{
		if(dataBlocks == null || dataBlocks.end < index) readNextDataBlock(); // If this index is not yet available, read an additional block
		return dataBlocks.values[(int)(index - dataBlocks.start)];
	}
	
	private void readNextIndicesBlock()
	{
		currentIndicesBlockRead++;
		if(indicesBlocks == null) indicesBlocks = new Block(0, indicesBlockSize - 1);
		else indicesBlocks = new Block(indicesBlocks.end + 1, indicesBlocks.end + indicesBlockSize);
		indicesBlocks.values = reader.int32().readArrayBlock("/" + Parameters.selection + "/indices", indicesBlockSize, currentIndicesBlockRead);
	}
	
	private int getIndices(long index)
	{
		if(indicesBlocks == null || indicesBlocks.end < index) readNextIndicesBlock(); // If this index is not yet available, read an additional block
		return indicesBlocks.values[(int)(index - indicesBlocks.start)];
	}
}

class Block
{
	public int[] values;
	public long start;
	public long end;
	
	public Block(long start, long end) 
	{
		this.start = start;
		this.end = end;
	}
}