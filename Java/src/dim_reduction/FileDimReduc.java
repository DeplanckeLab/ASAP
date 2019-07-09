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
    
   /* public static void readHeader(String header) throws IOException
    {
    	String[] cells = header.split("\t");
    	filtJSON.nber_cells = cells.length - 1; // Should be
    	if(filtJSON.nber_cells != parsJSON.nber_cells) new ErrorJSON("Detected different number of cells between parsingJSON("+parsJSON.nber_cells+") and Header("+filtJSON.nber_cells+")");
    	bw.write(header + "\n");
    	cellNames = new String[filtJSON.nber_cells];
    	for(int i = 1; i < cells.length; i++) cellNames[i-1] = cells[i];
    }*/
    
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

    /*public static double[][] getStatsOnFile(boolean getdataset) throws IOException
    {
    	BufferedReader br = new BufferedReader(new FileReader(Parameters.fileName));
    	String line = br.readLine(); // header
    	readHeader(line);
    	double [][] dataset = null;
    	if(Parameters.nbCellsDetected > filtJSON.nber_cells) new ErrorJSON("'Min Detected' should be smaller than the total number of cells/samples i.e. <=" + filtJSON.nber_cells);
    	if(getdataset) dataset = new double[filtJSON.nber_cells][parsJSON.nber_genes];
    	geneNames = new String[parsJSON.nber_genes];
    	colSum = new double[filtJSON.nber_cells];
    	rowSum = new double[parsJSON.nber_genes];
    	rowVar = new double[parsJSON.nber_genes];
    	rowCoeffOfVar = new double[parsJSON.nber_genes];
    	loggeomeans = new double[parsJSON.nber_genes]; // scLVM
    	sizeFactors = new double[filtJSON.nber_cells]; // scLVM
    	int nbGenes = 0;
    	line = br.readLine();
    	while(line != null)
    	{
    		String[] tokens = line.split("\t");
    		double mean = 0;
    	    double M2 = 0;
    		geneNames[nbGenes] = tokens[0];
     		for(int i = 1; i < tokens.length; i++)
    		{
    			double val = Double.parseDouble(tokens[i]);
    			double delta = val - mean;
    			mean = mean + delta/i;
    			M2 = M2 + delta*(val - mean);
    			colSum[i-1] += val;
    			rowSum[nbGenes] += val;
    			loggeomeans[nbGenes] += Math.log(val);
    			if(getdataset) dataset[i-1][nbGenes] = val;
    		}
    		loggeomeans[nbGenes] /= (tokens.length - 1);
    		rowVar[nbGenes] = M2 / (tokens.length - 2);
    		if(mean == 0) rowCoeffOfVar[nbGenes] = 0;
    		else rowCoeffOfVar[nbGenes] = Math.sqrt(rowVar[nbGenes]) / mean;
    		nbGenes++;
            line = br.readLine();
        }
    	br.close();
    	filtJSON.nber_genes = nbGenes;
    	if(filtJSON.nber_genes != parsJSON.nber_genes) new ErrorJSON("Detected different number of genes between parsingJSON("+parsJSON.nber_genes+") and Data Matrix("+filtJSON.nber_genes+")");
    	return dataset;
    }*/
}
