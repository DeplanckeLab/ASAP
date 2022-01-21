package differential_expression;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import config.Config;
import db.DBManager;
import hdf5.loom.LoomFile;
import json.ErrorJSON;
import json.WarningJSON;
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
    	sb.append("{").append("\"time_idle\":").append(Parameters.idleTime).append(",\"metadata\":[{\"name\":\"").append(Parameters.oAnnot).append("\",\"on\":\"GENE\",\"type\":\"NUMERIC\",\"nber_cols\":5,\"nber_rows\":");
		sb.append(nber_genes).append(",\"headers\":[\"log Fold-Change\",\"p-value\",\"FDR\",\"Avg. Exp. Group 1\",\"Avg. Exp. Group 2\"]").append("}]}");
		try
		{
    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
    		bw.write(sb.toString());
        	bw.close();
		}
		catch(IOException ioe) { new ErrorJSON(ioe.getMessage()); }
	}
	
	public static void findMarkers()
	{
		long nber_genes = 0;
		
		// Recuperate the indexes of all categories for this metadata
		if(Parameters.debugMode) DBManager.URL = Config.ConfigDEV().getURL("asap2_development"); // Annotations/Metadata are not in same DB
		else DBManager.URL = Config.ConfigMAIN().getURL("asap2_development");
		DBManager.connect();
		HashMap<String, Integer> categories = DBManager.getListCatJSON(Parameters.id);
		DBManager.disconnect();
		
    	// Open Loom file in read-only and read metadata
    	loom = new LoomFile("r", Parameters.loomFile);
   	   	if(!loom.exists(Parameters.iAnnot)) new ErrorJSON("Error in the Loom file. Path "+Parameters.iAnnot+" does not exist!");
    	Metadata groups = loom.fillInfoMetadata(Parameters.iAnnot, true);

    	// Error check
    	if(groups.categories.size() < 2) new ErrorJSON("Group metadata (" + Parameters.iAnnot + ") does not contain enough groups");
    	if(!groups.isCategorical()) new ErrorJSON("Group metadata (" + Parameters.iAnnot + ") is not a CATEGORICAL metadata");
    	
    	// Safety check
    	for(String cat:groups.categories) if(categories.get(cat) == null) new ErrorJSON("Categories in database with metadata id = " + Parameters.id + " do NOT match content of Loom file");
    	if(groups.categories.size() != categories.size()) new ErrorJSON("Categories in database with metadata id = " + Parameters.id + " do NOT match content of Loom file");
    	
    	// Read the Loom file line by line and perform DE
    	long[] dim = loom.getDimensions();
		int[] blockSize = loom.getChunkSizes();
		
		// Prepare final results
		HashMap<Integer, float[][]> results = new HashMap<Integer, float[][]>();
		StringBuilder sb = new StringBuilder();
		String prefix = "";
    	for(String g1:groups.categories)
    	{
        	// Check that g1 and g2 are in groups
        	long nbG1 = groups.categoriesMap.get(g1);
        	long nbG2 = groups.values.size() - nbG1;
        	if(nbG1 >= 3 && nbG2 >= 3) // Else I cannot compute
        	{
        		results.put(categories.get(g1), new float[5][(int)dim[0]]); // 1 "logFC", 2 "pval", 3 "FDR", 4 "AveG1", 5 "AveG2"
        	}
        	else
        	{
        		sb.append(prefix).append(g1);
        		prefix = ",";
        	}
    	}
    	if(!prefix.equals("")) WarningJSON.addWarning("Marker genes couldn't be computed for category(ies) "+ sb.toString() +" because there is less than 3 samples in this(ese) categories (or their complement)");
        	
    	// Prepare the indexes of cells/cols, for each group
    	HashMap<Integer, ArrayList<Integer>> indexGroupMap = new HashMap<Integer, ArrayList<Integer>>();
    	for(int i = 0; i < groups.values.size(); i++)
    	{
    	 	int g1 = categories.get(groups.values.get(i));
    	 	
    	 	// Group
    	 	ArrayList<Integer> listIndexes1 = indexGroupMap.get(g1);
    	 	if(listIndexes1 == null) listIndexes1 = new ArrayList<Integer>();
    	 	listIndexes1.add(i);
    	 	indexGroupMap.put(g1, listIndexes1);
    	 	
    	 	// Complement of other groups
	    	for(String g:groups.categories)
	    	{	
	    		int g2 = categories.get(g);
	    		if(g1 != g2)
	    		{
	        	 	ArrayList<Integer> listIndexes2 = indexGroupMap.get(-g2); // We take the negative index for it's complement (starts at 1, so works)
	        	 	if(listIndexes2 == null) listIndexes2 = new ArrayList<Integer>();
	        	 	listIndexes2.add(i);
	        	 	indexGroupMap.put(-g2, listIndexes2);
	    		}
	    	}
    	}
       	
    	// Read by chunks
		int nbTotalBlocks = (int)Math.ceil((double)dim[0] / blockSize[0]); // dim[0] === Nb genes ?

		// Start reading the /matrix (always exists / normalized on the go)
		MannWhitneyUTest test = new MannWhitneyUTest();
		nber_genes = 0;
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{		
			// Retrieve the blocks that will contain all columns (because we compute gene by gene)
			float[][] subMatrix = loom.readFloatBlock("/matrix", blockSize[0], (int)dim[1], nbBlocks, 0l);
			
			// Parsing Data and generating summary annotations
			for(int x = 0; x < subMatrix.length; x++)
			{
				int i = x + nbBlocks * blockSize[0]; // Original index of gene
				
		    	// For each group
		    	for(Integer g1:results.keySet())
		    	{				
		        	// Get indexes
		        	ArrayList<Integer> listIndexes1 = indexGroupMap.get(g1);
		        	ArrayList<Integer> listIndexes2 = indexGroupMap.get(-g1);
		        	double[] array1 = new double[listIndexes1.size()];
		        	for(int j = 0; j < listIndexes1.size(); j++) array1[j] = subMatrix[x][listIndexes1.get(j)];
		        	double[] array2 = new double[listIndexes2.size()];
		        	for(int j = 0; j < listIndexes2.size(); j++) array2[j] = subMatrix[x][listIndexes2.get(j)];

					//MannWhitneyTest test = new MannWhitneyTest(toCompute.get(Parameters.group_1), toCompute.get(Parameters.group_2));
					//pvals[i] = test.exactSP(); // In principle, in R, if there is ties, it should be the "approx" that is computed instead of the "exact"
					
		        	// Update results for this group
		        	float[][] res = results.get(g1);
		        	double mean1 = Utils.mean(array1);
		        	double mean2 = Utils.mean(array2);
		        	res[0][i] = (float)(Utils.log2(1 + mean1) - Utils.log2(1 + mean2));
		        	res[1][i] = (float)test.mannWhitneyUTest(array1, array2); // This one does not give exactly the same results as R, but is muchhhhhh faster than the JSC on
		        	res[3][i] = (float)mean1;
		        	res[4][i] = (float)mean2;
		        	results.put(g1, res);
		    	}
		    	
		    	nber_genes++;
			}
		}
		
		// Closing loom Read-only mode
		loom.close();	
	    
		// FDR calculation
    	for(Integer g1:results.keySet())
    	{
    		float[][] res = results.get(g1);
    		res[2] = Utils.p_adjust_F(res[1], "fdr");
    		results.put(g1, res);
    	}
    			
    	// Write result files
    	ArrayList<String> outputfiles = new ArrayList<String>();
    	for(Integer g1:results.keySet())
    	{
    		float[][] res = results.get(g1);
    		try
    		{
    			outputfiles.add("cat_" + g1 + ".tsv");
	    		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "cat_" + g1 + ".tsv"));
	    		bw.write("logfc\tpval\tfdr\tavg_"+g1+"\tavg_not_"+g1+"\n");
	    		for(int i = 0; i < res[0].length; i++)
	    		{
	    			sb = new StringBuilder();
	    			sb.append(res[0][i]).append("\t").append(res[1][i]).append("\t").append(res[2][i]).append("\t").append(res[3][i]).append("\t").append(res[4][i]).append("\n");
	    			bw.write(sb.toString());
	    		}
	    		bw.close();
    		}
    		catch(IOException ioe)
    		{
    			new ErrorJSON(ioe.getMessage());
    		}
    	}
    	
		// Prepare JSON
		sb = new StringBuilder();
    	sb.append("{").append("\"time_idle\":").append(Parameters.idleTime).append(",\"output_files\":").append(Utils.toString(outputfiles)).append(",\"nber_cols\":5,\"nber_rows\":");
		sb.append(nber_genes).append(",\"headers\":[\"log Fold-Change\",\"p-value\",\"FDR\",\"Avg. Exp. REF group\",\"Avg. Exp. COMPLEMENT\"]").append("}");
		
		// Writing results
    	Utils.writeJSON(sb);
	}
	
    public static long runWilcox()
    {  	
    	long nber_genes = 0;
    	
    	if(!loom.exists(Parameters.iAnnot)) new ErrorJSON("Error in the Loom file. Path "+Parameters.iAnnot+" should exist!");
    	
    	// Get the groups
    	System.out.print("Parsing Group metadata... ");
    	Metadata groups = loom.fillInfoMetadata(Parameters.gAnnot, true);
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


