package differential_expression;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import hdf5.loom.LoomFile;
import json.ErrorJSON;
import model.Metadata;
import model.Parameters;
import tools.Utils;

public class DE 
{
	private static LoomFile loom;
	
	public static void performDE()
	{
		long nber_genes = 0;
		
    	// Open Loom file in read-only
    	loom = new LoomFile("r", Parameters.loomFile);
    	loom.checkLoomFormat();
    	
    	// Run DE
		System.out.println("Performing DE on file : " + Parameters.loomFile);
		switch(Parameters.deModel) {
			case Wilcoxon:
				nber_genes = runWilcox();
				break;
			default:
				new ErrorJSON("This model is not yet implemented");
		}
		
		// Close loom file
		loom.close();
		
		// Write output.json
		StringBuilder sb = new StringBuilder();
    	sb.append("{").append("\"time_idle\":").append(Parameters.idleTime).append(",\"metadata\":[{\"name\":\"").append(Parameters.oAnnot).append("\",\"on\":\"row\",\"type\":\"NUMERIC\",\"nber_cols\":5,\"nber_rows\":");
		sb.append(nber_genes).append("}]}");
		try
		{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
    		bw.write(sb.toString());
        	bw.close();
		}
		catch(IOException ioe) { new ErrorJSON(ioe.getMessage()); }
	}
	
    public static long runWilcox()
    {  	
    	long nber_genes = 0;
    	
    	if(!loom.exists(Parameters.iAnnot)) new ErrorJSON("Error in the Loom file. Path "+Parameters.iAnnot+" should exist!");
    	
    	// Get the groups
    	System.out.print("Parsing Group metadata... ");
    	Metadata groups = loom.readMetadata(Parameters.gAnnot);
    	System.out.println("Group metadata (" + Parameters.gAnnot + ") contains " + groups.categories.size() + " groups");
    	
    	// Handle the case where g2 is complementary group
    	if(Parameters.group_2 == null)
    	{
    		Parameters.group_2 = "Complement_ASAP";
        	for(long i = 0; i < groups.values.size(); i++)
        	{
        		String g = groups.values.get(i);
        		if(!g.equals(Parameters.group_1)) groups.values.set(i, "Complement_ASAP");      		
        	}
    	}
    	
    	// Check that g1 and g2 are in groups
    	int nbG1 = 0, nbG2 = 0;
    	for(long i = 0; i < groups.values.size(); i++)
    	{
    		String g = groups.values.get(i);
    		if(g.equals(Parameters.group_1)) nbG1++;
    		else if(Parameters.group_2 == null) nbG2++;
    		if(Parameters.group_2 != null && g.equals(Parameters.group_2)) nbG2++;
    	}
    	if(nbG1 < 3) { loom.close(); new ErrorJSON("Group 1 should contain at least 3 samples"); }
    	if(nbG2 < 3) { loom.close(); new ErrorJSON("Group 2 should contain at least 3 samples"); }
    	System.out.println("Your reference group (G1) contains " + nbG1 + " cells, while the comparison group (G2) has " + nbG2);
    	
    	// Read the Loom file line by line and perform DE
    	long[] dim = loom.getDimensions();
		int[] blockSize = loom.getChunkSizes();
		
		float[][] results = new float[5][(int)dim[0]]; // 0 logFC", 1 "pval", 2 "FDR", 3 "AveG1", 4 "AveG2"
		
		int nbTotalBlocks = (int)Math.ceil((double)dim[0] / blockSize[0]); // dim[0] === Nb genes ?
		System.out.println("Reading " + nbTotalBlocks + " independent blocks...");
		
		MannWhitneyUTest test = new MannWhitneyUTest();
		nber_genes = 0;
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{		
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = loom.readFloatBlock("/matrix", blockSize[0], (int)dim[1], nbBlocks, 0l);
			
			// Parsing Data and generating summary annotations
			for(int x = 0; x < subMatrix.length; x++)
			{
				int i = x + nbBlocks * blockSize[0]; // Original index of gene
				
				HashMap<String, double[]> toCompute = new HashMap<String, double[]>(); // I should do DoubleArray64 here, but Wilcox test will not handle it
				toCompute.put(Parameters.group_1, new double[nbG1]);
				toCompute.put(Parameters.group_2, new double[nbG2]);
				HashMap<String, Integer> indexes = new HashMap<String, Integer>();
				indexes.put(Parameters.group_1, 0);
				indexes.put(Parameters.group_2, 0);
				HashMap<String, Float> means = new HashMap<String, Float>();
				means.put(Parameters.group_1, 0f);
				means.put(Parameters.group_2, 0f);
				
				for(int j = 0; j < subMatrix[0].length; j++) // Index of cell
				{
					String g = groups.values.get(j);
					double[] array = toCompute.get(g);
					if(array != null) // If this is one of the groups we are interested in
					{
						int index = indexes.get(g);
						array[index] = subMatrix[x][j];
						means.put(g, means.get(g) + (subMatrix[x][j] / array.length));
						indexes.put(g, index + 1);
					}
				}
				//MannWhitneyTest test = new MannWhitneyTest(toCompute.get(Parameters.group_1), toCompute.get(Parameters.group_2));
				//pvals[i] = test.exactSP(); // In principle, in R, if there is ties, it should be the "approx" that is computed instead of the "exact"
				results[0][i] = (float)(Utils.log2(1 + means.get(Parameters.group_1)) - Utils.log2(1 + means.get(Parameters.group_2)));
				results[1][i] = (float)test.mannWhitneyUTest(toCompute.get(Parameters.group_1), toCompute.get(Parameters.group_2)); // This one does not give exactly the same results as R, but is muchhhhhh faster than the JSC on
				results[3][i] = (float)means.get(Parameters.group_1);
				results[4][i] = (float)means.get(Parameters.group_2);
				
				nber_genes++;
			}
		}
		
		// FDR calculation
		results[2] = Utils.p_adjust_F(results[1], "fdr");
		
		// Closing loom Read-only mode
		loom.close();
		
		// Write results
		loom = new LoomFile("r+", Parameters.loomFile);
		loom.writeMatrixMetadata(Parameters.oAnnot, Utils.t(results));
		
		// Return nbGenes (can be filtered)
		return nber_genes;
	}
}


