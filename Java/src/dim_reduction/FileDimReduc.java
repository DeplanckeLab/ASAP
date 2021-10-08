package dim_reduction;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import com.jujutsu.tsne.barneshut.BHTSne;
import com.jujutsu.tsne.barneshut.BarnesHutTSne;
import com.jujutsu.tsne.barneshut.ParallelBHTsne;
import com.jujutsu.tsne.barneshut.TSneConfiguration;
import com.jujutsu.utils.MatrixUtils;
import com.jujutsu.utils.TSneUtils;

import json.DimReducJSON;
import json.ErrorJSON;
import model.Parameters;
import tools.Utils;

public class FileDimReduc 
{
	public static DimReducJSON dimReducJSON = null;
	
    public static void reduceDimension()
    {
    	dimReducJSON = new DimReducJSON(); // Output JSON
        System.out.println("Reducing dimension of file : " + Parameters.fileName);
        try
        {
        	switch(Parameters.dimReducModel)
        	{
        		case TSNE:
        			reduceTSNE();
        	        break;
        		default:
        			System.err.println("Not implemented yet");
        	}
        } 
        catch(IOException ioe)
        {
        	System.err.println(ioe.getMessage());
        	try
        	{
        		new ErrorJSON(ioe.getMessage());
            	BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.json"));
            	bw.write("{\"displayed_error\":\"" + ioe.getMessage() + "\"}");
            	bw.close();
        	}
        	catch(IOException ioe2)
        	{
        		System.err.println(ioe2.getMessage());
        	}
        	System.exit(-1);
        }
        dimReducJSON.writeJSON();
    }
       
    public static void reduceTSNE() throws IOException
    {
    	// Analysis
    	System.out.println("Reading file...");
    	double [][] matrix = MatrixUtils.simpleRead2DMatrix(new File(Parameters.fileName), "\t");
    	System.out.println("Read!");
    	System.gc();
    	System.gc();
    	matrix = Utils.t(matrix);
    	System.out.println("Transposed is computed!");
    	System.gc();
    	System.gc();
    	BarnesHutTSne tsne;
    	boolean parallel = true;
    	if(parallel) {			
    		tsne = new ParallelBHTsne();
    	} else {
    		tsne = new BHTSne();
    	}
    	TSneConfiguration config = TSneUtils.buildConfig(matrix, 2, 55, Parameters.perplexity, 2000, false, 0.5, false, true);

    	System.out.println("Run tSNE!");
    	double [][] Y = tsne.tsne(config);   
    	System.out.println(Y.length);
    	System.out.println(Y[0].length);
    	    // Plot Y or save Y to file and plot with some other tool such as for instance R
    	
    	/*BufferedReader br = new BufferedReader(new FileReader(Parameters.fileName));  // Read twice (because it's 50% faster, writing the file is much faster)
    	String line = br.readLine(); // Skip Header
    	line = br.readLine();
    	int nbGenes = 0;
    	while(line != null)
    	{
 
       		line = br.readLine();
    	}
       	br.close();*/
    }
}
