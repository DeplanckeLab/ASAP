package scaling;

import hdf5.loom.LoomFile;
import model.Parameters;
import tools.Utils;

public class Scaling 
{
	/**
	 * Normalization main function (adapted from Seurat method)
	 */	
	public static void runScaling()
	{
    	// Open Loom file in read/write mode
    	LoomFile loom = new LoomFile("r+", Parameters.loomFile); // TODO I could also create a tmp Loom file to avoid locking the file too long
    	loom.checkLoomFormat();
		long[] dim = loom.getDimensions(Parameters.iAnnot);
		//System.out.println("Current /matrix size = " + dim[0] + " x " + dim[1]);
	
		// Process Main Matrix
		int[] blockSize = loom.getChunkSizes(Parameters.iAnnot);
		loom.createEmptyFloat32MatrixDataset(Parameters.oAnnot, dim[0], dim[1], blockSize[0], blockSize[1]);
    	    	
    	// Read the file block per block according to natural storage
		int nbTotalBlocks = (int)Math.ceil((double)dim[0] / blockSize[0]); // dim[0] == nb genes
 		//System.out.println("Writing " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{		
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = loom.readFloatBlock(Parameters.iAnnot, blockSize[0], (int)dim[1], nbBlocks, 0l);

			// Adapting the size to block size if last block
			if(subMatrix.length < blockSize[0])
			{
				float[][] tmpMatrix = new float[blockSize[0]][(int)dim[1]];
				for(int i = 0; i < subMatrix.length; i++) tmpMatrix[i] = subMatrix[i];
				subMatrix = tmpMatrix;
			}
			
			// Parsing Data and generating summary annotations
			for(int x = 0; x < subMatrix.length; x++) // For every row
			{
				int i = x + nbBlocks * blockSize[0]; // Original index
				if(i < dim[0]) // In case the block is bigger than the number of genes
				{
					float rowMean = Utils.mean(subMatrix[x]);
					double rowSdev = 0;
				    if(Parameters.scale)
				    {
						for(int j = 0; j < subMatrix[0].length; j++) // Go through all columns/cells
						{
					    	if(Parameters.center) rowSdev += (subMatrix[x][j] - rowMean) * (subMatrix[x][j] - rowMean);
						    else rowSdev += subMatrix[x][j] * subMatrix[x][j];
						}
						rowSdev = Math.sqrt(rowSdev / (dim[1] - 1));
					} 
				    else rowSdev = 1;
				    if(!Parameters.center) rowMean = 0;

					for(int j = 0; j < subMatrix[0].length; j++) // Go through all columns/cells
					{
						subMatrix[x][j] = (float)((subMatrix[x][j] - rowMean) / rowSdev);
						if(subMatrix[x][j] > Parameters.scale_max) subMatrix[x][j] = Parameters.scale_max;
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
