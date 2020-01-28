package normalization;

import bigarrays.DoubleArray64;
import hdf5.loom.LoomFile;
import json.ErrorJSON;
import model.Parameters;
import tools.Utils;

public class Normalization 
{
	/**
	 * Normalization main function (adapted from Seurat method)
	 */	
	public static void runNormalization()
	{
    	// Open Loom file in read/write mode
    	LoomFile loom = new LoomFile("r+", Parameters.loomFile); // TODO I could also create a tmp Loom file to avoid locking the file too long
    	loom.checkLoomFormat();
    	if(!loom.exists("/col_attrs/_Depth")) new ErrorJSON("Error in the Loom file. Path /col_attrs/_Depth should exist!");
    	DoubleArray64 depth = loom.readDoubleArray("/col_attrs/_Depth");
		long[] dim = loom.getDimensions("/matrix");
		//System.out.println("Current /matrix size = " + dim[0] + " x " + dim[1]);
	
		// Process Main Matrix
		int[] blockSize = loom.getChunkSizes("/matrix");
		loom.createEmptyFloat32MatrixDataset(Parameters.oAnnot, dim[0], dim[1], blockSize[0], blockSize[1]);
    	    	
    	// Read the file block per block according to natural storage
		int nbTotalBlocks = (int)Math.ceil((double)dim[0] / blockSize[0]); // dim[0] == nb genes
 		//System.out.println("Writing " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{		
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = loom.readFloatBlock("/matrix", blockSize[0], (int)dim[1], nbBlocks, 0l);

			// Adapting the size to block size if last block
			if(subMatrix.length < blockSize[0])
			{
				float[][] tmpMatrix = new float[blockSize[0]][(int)dim[1]];
				for(int i = 0; i < subMatrix.length; i++) tmpMatrix[i] = subMatrix[i];
				subMatrix = tmpMatrix;
			}
			
			// Parsing Data and generating summary annotations
			for(int x = 0; x < subMatrix.length; x++)
			{
				int i = x + nbBlocks * blockSize[0]; // Original index
				if(i < dim[0]) // In case the block is bigger than the number of genes
				{
					for(int j = 0; j < subMatrix[0].length; j++) // Go through all columns/cells
					{
						subMatrix[x][j] = (float)Utils.log2(1 + (subMatrix[x][j] / depth.get(j)) * Parameters.scale_factor);
					}
				}
			}
			
			// Writing this block to output
			loom.writeFloatBlockDataset(Parameters.oAnnot, subMatrix, nbBlocks, 0);
		}
    	
		// Cause writing fixed-size blocks can extend the matrix size with 0
		loom.resizeDataset(Parameters.oAnnot, dim[0], dim[1]);
		
		// Prepare output String
		StringBuilder sb = new StringBuilder();
    	sb.append("{").append("\"time_idle\":").append(Parameters.idleTime);
    	sb.append(",\"nber_rows\":").append(dim[0]);
    	sb.append(",\"nber_cols\":").append(dim[1]);
		sb.append("}");

		// Close Loom file
		loom.close();
		
    	// Write output.json
    	Utils.writeJSON(sb, Parameters.JSONFileName);
	}
}
