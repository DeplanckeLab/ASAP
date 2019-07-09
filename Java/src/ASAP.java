import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import bigarrays.IntArray64;
import bigarrays.StringArray64;
import config.Config;
import db.DBManager;
import db.EnsemblDB;
import db.GODatabase;
import differential_expression.DE;
import dim_reduction.FileDimReduc;
import enrichment.Enrichment;
import filtering.FileFilter;
import hdf5.loom.LoomFile;
import json.ErrorJSON;
import json.FilterCellsJSON;
import model.MetaOn;
import model.Metadata;
import model.Metatype;
import model.Mode;
import model.Parameters;
import parsing.FileParser;
import tools.Utils;

/**
 * @author gardeux
 */
public class ASAP 
{
	public static Mode m = null;

	public static void main(String[] args)
	{
		DBManager.JDBC_DRIVER = Config.getProperty("mDbDriv");
		if(DBManager.JDBC_DRIVER.equals("com.mysql.jdbc.Driver")) DBManager.URL = "jdbc:mysql://";
		else if(DBManager.JDBC_DRIVER.equals("org.postgresql.Driver")) DBManager.URL = "jdbc:postgresql://";
		DBManager.URL += Config.getProperty("mDbHost") + "?user=" + Config.getProperty("mDbUser") + "&password=" + Config.getProperty("mDbPwds");
		
		String[] args2 = readMode(args);
		Parameters.load(args2, m);
		switch(m)
		{
			case CreateEnrichmentDB:
				DBManager.JDBC_DRIVER = "com.mysql.jdbc.Driver";
				DBManager.URL = "jdbc:mysql://mysql.ebi.ac.uk:4085/go_latest?user=go_select&password=amigo";
				DBManager.connect();
				GODatabase.generateGODB(Parameters.outputFolder, Parameters.taxon);
				DBManager.disconnect();
				break;
			case Enrichment:
				long t = System.currentTimeMillis();
				Enrichment.readFiles();
				Enrichment.runEnrichment();
				System.out.println("Enrichment time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case DimensionReduction: 
				t = System.currentTimeMillis();
				FileDimReduc.reduceDimension();
				System.out.println("Dimension reduction time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case DifferentialExpression: 
				t = System.currentTimeMillis();
				DE.performDE();
				System.out.println("Differential expression time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case Preparsing: 
				t = System.currentTimeMillis();
				FileParser.preparse();
				System.out.println("Preparsing time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case Parsing: 
				t = System.currentTimeMillis();
				DBManager.connect();
				DBManager.getGenesInDB(Parameters.organism);
				DBManager.disconnect();
				System.out.println("Accessing DB time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				t = System.currentTimeMillis();
				FileParser.parse();
				System.out.println("Parsing time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case PreparseMetadata: 
				t = System.currentTimeMillis();
				FileParser.preparseMetadata();
				System.out.println("PreparseMetadata time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case ParseMetadata: 
				t = System.currentTimeMillis();
				FileParser.parseMetadata();
				System.out.println("ParsMetadata time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case ExtractRow: 
				LoomFile loom = new LoomFile("r", Parameters.fileName);
				//t = System.currentTimeMillis();
				if(Parameters.index == -1) Parameters.index = loom.getGeneIndex(Parameters.geneName);
				if(Parameters.index == -1) new ErrorJSON("Gene not found");
				float[] row = loom.readRow(Parameters.index, Parameters.iAnnot); // TODO does not work if too big array
				loom.close();
				
				// Creating output string
	    		StringBuilder sb = new StringBuilder();	
	    		sb.append("{\"index\":").append(Parameters.index).append(",").append("\"row\":[");
	    		String prefix = "";
		    	for(int i = 0; i < row.length; i++) 
		    	{
		    		sb.append(prefix).append(Utils.format(row[i]));
		    		prefix = ",";
		    	}
		    	sb.append("]}");
				
				// Writing results
		    	writeJSON(sb, Parameters.JSONFileName);

		    	//System.out.println("ExtractRow time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case ExtractCol: 
				loom = new LoomFile("r", Parameters.fileName);
				//t = System.currentTimeMillis();
				if(Parameters.index == -1) Parameters.index = loom.getCellIndex(Parameters.cellName);
				if(Parameters.index == -1) new ErrorJSON("Cell not found");
				float[] col = loom.readCol(Parameters.index, Parameters.iAnnot); // TODO does not work if too big array
				loom.close();
				
				// Creating output string
	    		sb = new StringBuilder();	
	    		sb.append("{\"index\":").append(Parameters.index).append(",").append("\"col\":[");
	    		prefix = "";
		    	for(int i = 0; i < col.length; i++) 
		    	{
		    		sb.append(prefix).append(Utils.format(col[i]));
		    		prefix = ",";
		    	}
		    	sb.append("]}");
				
				// Writing results
		    	writeJSON(sb, Parameters.JSONFileName);

		    	//System.out.println("ExtractRow time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case ListMetaData:
				//t = System.currentTimeMillis();
				
				// Get Metadata infos
				loom = new LoomFile("r", Parameters.fileName);
				sb = Metadata.toString(loom.listMetadata());
				loom.close();
				
	    		// Writing results
				writeJSON(sb, Parameters.JSONFileName);
				//System.out.println("ListMetaData time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case ExtractMetaData:
				loom = new LoomFile("r", Parameters.fileName);
				//t = System.currentTimeMillis();
				
				// Extract metadata from Loom
				Metadata metadata;
				metadata = loom.readMetadata(Parameters.metaName); // We could in principle simplify this part if we know already the data type / if we want to output the values / etc...
				if(Parameters.metatype != null) metadata.type = Parameters.metatype; 
				
				// Extract cell names if needed
				StringArray64 cellNames = null;
				if(Parameters.details && metadata.on == MetaOn.CELL) cellNames = loom.getCellNames();
				loom.close();
				
				// Prepare the output string
				sb = new StringBuilder();
				metadata.addMeta(sb, !Parameters.evenLessDetails, cellNames, Long.MAX_VALUE);
	    		
	    		// Writing results
				writeJSON(sb, Parameters.JSONFileName);
				//System.out.println("ExtractMetaData time: " + Utils.toReadableTime(System.currentTimeMillis() - t) + " - Idle time = " + Parameters.idleTime + " s");
				break;
			case RemoveMetaData:
				//t = System.currentTimeMillis();

				// Open in writing mode and remove the metadata
				loom = new LoomFile("w+", Parameters.loomFile);
				loom.removeMetadata(Parameters.metaName);
				loom.close();

				// Create output.json
				Metadata meta = new Metadata(Parameters.metaName);
				
				// Create output String
				sb = new StringBuilder();
				meta.addMeta(sb, false, null, Long.MAX_VALUE);
				
	    		// Writing results
				writeJSON(sb, Parameters.outputFolder + "output.json");
				
				//System.out.println("RemoveMetaData time: " + Utils.toReadableTime(System.currentTimeMillis() - t) + " - Idle time = " + Parameters.idleTime + " s");
				break;
			case CreateCellSelection:
				//t = System.currentTimeMillis();
				
				// Read Selected Cells in JSON
				long[] selectedIndexes = FilterCellsJSON.parseJSON(Parameters.JSONFileName).selected_cells; // TODO if too many are selected, this does not work
				
				// Create Metadata and write in Loom file
				loom = new LoomFile("w+", Parameters.loomFile);
				IntArray64 newMetadata = new IntArray64(loom.getDimensions()[1]);
				for(long i:selectedIndexes) newMetadata.set(i, 1);
				loom.writeIntArray(Parameters.metaName, newMetadata);
				loom.close();

				// Create output.json
				HashMap<String, Long> cat = new HashMap<String, Long>();
				cat.put("0", newMetadata.size() - selectedIndexes.length);
				cat.put("1", (long)selectedIndexes.length);
				meta = new Metadata(Parameters.metaName, Metatype.DISCRETE, MetaOn.CELL, newMetadata.size(), 1, cat);
				
				// Create output String
				sb = new StringBuilder();
				meta.addMeta(sb, false, null, Long.MAX_VALUE);

	    		// Writing results
				writeJSON(sb, Parameters.outputFile);
				//System.out.println("CreateCellSelection time: " + Utils.toReadableTime(System.currentTimeMillis() - t) + " - Idle time = " + Parameters.idleTime + " s");
				break;
			case CopyMetaData:
				//t = System.currentTimeMillis();
				LoomFile loomFrom = new LoomFile("r", Parameters.loomFile);
				LoomFile loomTo = new LoomFile("w+", Parameters.loomFile2);

				LoomFile.copyMetadata(Parameters.metaName, loomFrom, loomTo);

				loomFrom.close();
				loomTo.close();
				//System.out.println("CopyMetaData time: " + Utils.toReadableTime(System.currentTimeMillis() - t) + " - Idle time = " + Parameters.idleTime + " s");
				break;
			case UpdateEnsemblDB:
				EnsemblDB.updateDB();
				break;
			case RegenerateNewOrganism:
				//TODO
				new ErrorJSON("Not implemented yet");
				//RegenerateNewOrganism.regenerateJSON();
				break;
			case FilterCells: 
				//t = System.currentTimeMillis();
				FileFilter.filterCells();
				//System.out.println("Filtering Cells time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
			case FilterGenes: 
				//t = System.currentTimeMillis();
				FileFilter.filterGenes();
				//System.out.println("Filtering Genes time: " + Utils.toReadableTime(System.currentTimeMillis() - t));
				break;
		}
	}
	
	public static void writeJSON(StringBuilder content, String outputJSONFile)
	{
		if(outputJSONFile == null) System.out.println(content.toString());
		else
		{
			try
			{
	    		BufferedWriter bw = new BufferedWriter(new FileWriter(outputJSONFile));
	    		bw.write(content.toString());
	        	bw.close();
			}
			catch(IOException ioe) { new ErrorJSON(ioe.getMessage()); }
		}
	}
			
	public static String[] readMode(String[] args)
	{
		String[] args2 = null;
		if(args.length >= 2)
		{
			args2 = new String[args.length - 2];
			int j = 0;
			for(int i = 0; i < args.length; i++) 
			{
				String arg = args[i];
				switch(arg)
				{
					case "-T":
						i++;
						String mode = args[i];
						switch(mode)
						{
							case "CreateEnrichmentDB": m = Mode.CreateEnrichmentDB; break;
							case "Enrichment": m = Mode.Enrichment; break;
							case "Parsing": m = Mode.Parsing; break;
							case "Preparsing": m = Mode.Preparsing; break;
							case "DimensionReduction": m = Mode.DimensionReduction; break;
							case "UpdateEnsemblDB": m = Mode.UpdateEnsemblDB; break;
							case "RegenerateNewOrganism": m = Mode.RegenerateNewOrganism; break;
							case "ExtractRow": m = Mode.ExtractRow; break;
							case "ExtractCol": m = Mode.ExtractCol; break;
							case "ListMetadata": m = Mode.ListMetaData; break;
							case "ExtractMetadata": m = Mode.ExtractMetaData; break;
							case "RemoveMetaData": m = Mode.RemoveMetaData; break;
							case "CopyMetaData": m = Mode.CopyMetaData; break;
							case "CreateCellSelection": m = Mode.CreateCellSelection; break;
							case "PreparseMetadata": m = Mode.PreparseMetadata; break;
							case "ParseMetadata": m = Mode.ParseMetadata; break;
							case "FilterCols": m = Mode.FilterCells; break;
							case "FilterRows": m = Mode.FilterGenes; break;
							case "DifferentialExpression": m = Mode.DifferentialExpression; break;
							default: System.err.println("Mode (-T) " + mode + " does not exist!"); System.out.println("-T %s \t\tMode to run ASAP [Preparsing, Parsing, PreparseMetadata, CreateCellSelection, ParseMetadata, RegenerateOutput, CreateEnrichmentDB, DifferentialExpression, UpdateEnsemblDB, Enrichment, ExtractRow, ExtractCol, ListMetadata, ExtractMetadata, RemoveMetaData, CopyMetaData, FilterCols, FilterRows]."); System.exit(-1);
						}
						break;
					default:
						args2[j] = arg;
						j++;
				}
			}
		}
		if(m == null || args.length < 2)
		{
			System.out.println("Argument -T is mandatory:");
			System.out.println("-T %s \t\tMode to run ASAP [Preparsing, Parsing, PreparseMetadata, CreateCellSelection, ParseMetadata, RegenerateOutput, CreateEnrichmentDB, DifferentialExpression, UpdateEnsemblDB, Enrichment, ExtractRow, ExtractCol, ListMetadata, ExtractMetadata, RemoveMetaData, CopyMetaData, FilterCols, FilterRows].");
			System.exit(-1);
		}
		return args2;
	}
}
