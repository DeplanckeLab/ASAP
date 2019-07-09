package tools;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

import model.Parameters;

public class Stats 
{
	private static Rengine re = null;
	
	public static void startRHandle(boolean showConsole)
	{
		if (!Rengine.versionCheck()) System.err.println("** Version mismatch - Java files don't match library version.");
		if(re == null) 
		{
			if(showConsole) re = new Rengine(new String[] {"--vanilla"}, false, new Console());
			else re = new Rengine(new String[] {"--vanilla"}, false, null);
			//re.eval("if(\"mixtools\" %in% rownames(installed.packages()) == FALSE) {install.packages(\"mixtools\", quiet=TRUE, repos=\"http://cran.us.r-project.org\", lib=\"/gsfs2/home/u30/gardeux/R/x86_64-unknown-linux-gnu-library/3.0\")}");
			//re.eval("library(mixtools)");
			if(!re.waitForR()) System.err.println("Cannot load R");
			else System.out.println("R loaded");
		}
		else System.err.println("R is already running");
	}
		
	public static void scLVM(double[][] normalizedDataset)
	{
		long t = System.currentTimeMillis();
		assignDoubleMatrix(normalizedDataset, "data.parsed.sf");
    	System.out.println("Assign Double Matrix: " + Utils.toReadableTime(System.currentTimeMillis() - t));
		re.eval("require(scLVM)");
		re.eval("png("+Parameters.outputFolder+"\"tech.noise.fit.png\", width=500, height=600, type=\"cairo\")");
		re.eval("data.tech.noise <- fitTechnicalNoise(data.parsed.sf, use_ERCC = F, fit_type = \""+Parameters.fitModel+"\")");
		re.eval("dev.off()");
/*data.plots <- rbind(data.plots, data.frame(name=\"tech.noise.fit.png\", description=\"scLVM fit of technical noise\"))\r\n" + 
				"    png(paste0(output.folder,\"/variable.genes.png\"), width=500, height=600, type=\"cairo\")\r\n" + 
				"    data.variable.genes <- getVariableGenes(data.parsed.sf, data.tech.noise$fit, fit_type = fit.model, threshold = 0.1)\r\n" + 
				"    dev.off()\r\n" + 
				"    data.plots <- rbind(data.plots, data.frame(name=\"variable.genes.png\", description=\"scLVM fit of technical noise + variable genes\"))\r\n" + 
				"    data.out <- as.data.frame(data.parsed[data.variable.genes, ])");*/
		//re.eval("data.count.random$V2 <- data.count.random$V2 / sum(data.count.random$V2)");
	}
	
	public static void stopRHandle()
	{
		re.end();
		System.out.println("R unloaded");
	}
	
    private static REXP assignDoubleMatrix(double[][] matrix, String nameToAssignOn) 
    {
        re.assign(nameToAssignOn, matrix[0]);
        REXP resultMatrix = re.eval(nameToAssignOn + " <- matrix( " + nameToAssignOn + " ,nr=1)");
        for (int i = 1; i < matrix.length; i++) 
        {
            re.assign("temp", matrix[i]);
            resultMatrix = re.eval(nameToAssignOn + " <- rbind(" + nameToAssignOn + ",matrix(temp,nr=1))");
        }
        return resultMatrix;
    }
}

class Console implements RMainLoopCallbacks
{
	public void rBusy(Rengine arg0, int arg1) {
	}

	public String rChooseFile(Rengine arg0, int arg1) {
		return null;
	}

	public void rFlushConsole(Rengine arg0) {
	}

	public void rLoadHistory(Rengine arg0, String arg1) {
	}

	public String rReadConsole(Rengine arg0, String arg1, int arg2) {
		return null;
	}

	public void rSaveHistory(Rengine arg0, String arg1) {
	}

	public void rShowMessage(Rengine arg0, String arg1) {
	}

	public void rWriteConsole(Rengine arg0, String arg1, int arg2) {
		System.out.println(arg1);
	}
}
