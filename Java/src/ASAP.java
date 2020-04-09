import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import bigarrays.IntArray64;
import bigarrays.LongArray64;
import bigarrays.StringArray64;
import config.Config;
import db.DBManager;
import db.EnsemblDB;
import db.GODatabase;
import db.KEGGRestApi;
import differential_expression.DE;
import dim_reduction.FileDimReduc;
import enrichment.Enrichment;
import filtering.FileFilter;
import hdf5.loom.LoomFile;
import json.CopyMetaJSON;
import json.ErrorJSON;
import json.FilterCellsJSON;
import model.MetaOn;
import model.Metadata;
import model.Metatype;
import model.Mode;
import model.Parameters;
import normalization.Normalization;
import parsing.FileParser;
import scaling.Scaling;
import tools.Utils;

/**
 * @author gardeux
 */
public class ASAP 
{
	public static Mode m = null;
	
	public static void main(String[] args)
	{
		DBManager.JDBC_DRIVER = "org.postgresql.Driver";
		DBManager.URL = "jdbc:postgresql://" + Config.getProperty("mDbHost") + "?user=" + Config.getProperty("mDbUser") + "&password=" + Config.getProperty("mDbPwds");
		
		String[] args2 = readMode(args);
		Parameters.load(args2, m);
		switch(m)
		{
			case CreateKeggDB:
				KEGGRestApi.generateKEGGDB();
				break;
			case CreateGODB:
				GODatabase.generateGODB();
				break;
			case Enrichment:
				Enrichment.runEnrichment();
				break;
			case DimensionReduction: 
				FileDimReduc.reduceDimension();
				break;
			case DifferentialExpression: 
				DE.performDE();
				break;
			case Normalization: 
				Normalization.runNormalization();
				break;
			case Scaling: 
				Scaling.runScaling();
				break;
			case Preparsing: 
				FileParser.preparse();
				break;
			case Parsing: 
				DBManager.connect();
				DBManager.getGenesInDB(Parameters.organism);
				DBManager.disconnect();
				FileParser.parse();
				break;
			case PreparseMetadata: 
				FileParser.preparseMetadata();
				break;
			case ParseMetadata: 
				FileParser.parseMetadata();
				break;
			case ExtractRow: 
				LoomFile loom = new LoomFile("r", Parameters.loomFile);

				// First check the total dimension of the dataset to extract from
				long[] dim = loom.getDimensions("/matrix");
				long nbGenes = dim[0];
				long nbCells = dim[1];
				dim = loom.getDimensions(Parameters.iAnnot); // dim[1] = nb of values to extract for each row
				if(dim.length == 1) { loom.close(); new ErrorJSON("The dataset is not a matrix, it's an array. You should use ExtractMetaData instead"); }
				StringArray64 names = null;
				
				// Retrieve names if asked
				if(Parameters.displayNames)
				{
					if(dim[1] == nbCells) names = loom.getCellNames(); // depends of the metadata/dataset to extract
					else if(dim[1] == nbGenes) names = loom.getGeneHGNC();
					else new ErrorJSON("Not sure what to do here...");
				}
				
				// Retrieve indexes if exist
				if(Parameters.indexes == null) 
				{
					if(Parameters.names != null)
					{
						if(dim[0] == nbCells) Parameters.indexes = loom.getCellIndexes(Parameters.names); // depends of the metadata/dataset to extract
						else if(dim[0] == nbGenes) Parameters.indexes = loom.getGeneIndexes(Parameters.names);
						else new ErrorJSON("Not sure what to do here...");
					}
					else if(Parameters.stable_ids != null)
					{
						if(dim[0] == nbCells) Parameters.indexes = loom.getCellIndexesByStableIds(Parameters.stable_ids); // depends of the metadata/dataset to extract
						else if(dim[0] == nbGenes) Parameters.indexes = loom.getGeneIndexesByStableIds(Parameters.stable_ids);
						else new ErrorJSON("Not sure what to do here...");
					} // the other case should be handled in Parameters
				}
				
				// Creating output string
	    		StringBuilder sb = new StringBuilder();
	    		sb.append("{\"values\":[");
	    		String prefix = "";
	    		if(dim[0] < 100) // TODO replace this by handling CONTIGUOUS / CHUNKED + Handle multiple read in same block
	    		{
	    			// Load the whole thing in memory
	    			float[][] matrix = loom.readFloatMatrix(Parameters.iAnnot);
	    			for(long index:Parameters.indexes)
					{
						sb.append(prefix);
						if(index == -1) sb.append("null");
						else
						{
							sb.append("[");
							String prefix2 = "";
							for(int i = 0; i < matrix[(int)index].length; i++)
							{
								sb.append(prefix2).append(Utils.format(matrix[(int)index][i])); // TODO does not work if too big array
								prefix2 = ",";
							}
							sb.append("]");
						}
						prefix = ",";
					}
	    		}
	    		else
	    		{
					for(long index:Parameters.indexes)
					{
						sb.append(prefix);
						if(index == -1) sb.append("null");
						else
						{
							sb.append("[");
							float[] row = loom.readRow(index, Parameters.iAnnot); // TODO does not work if too big array // TODO handle multiple read in same block
							String prefix2 = "";
							for(int i = 0; i < row.length; i++)
							{
								sb.append(prefix2).append(Utils.format(row[i]));
								prefix2 = ",";
							}
							sb.append("]");
						}
						prefix = ",";
					}
					sb.append("]");
	    		}
				loom.close();

				// If names to be added
				if(Parameters.displayNames)
				{
					sb.append(",\"names\":[");
					prefix = "";
					for(String n:names)
					{
						sb.append(prefix).append("\"").append(n).append("\"");
						prefix = ",";
					}
					sb.append("]");
				}
				sb.append("}");
				
				// Writing results
		    	Utils.writeJSON(sb, Parameters.JSONFileName);
				break;
			case ExtractCol: 
				loom = new LoomFile("r", Parameters.loomFile);

				// First check the total dimension of the dataset to extract from
				dim = loom.getDimensions("/matrix");
				nbGenes = dim[0];
				nbCells = dim[1];
				dim = loom.getDimensions(Parameters.iAnnot); // dim[1] = nb of values to extract for each row
				if(dim.length == 1) { loom.close(); new ErrorJSON("The dataset is not a matrix, it's an array. You should use ExtractMetaData instead"); }
				names = null;
				
				// Retrieve names if asked
				if(Parameters.displayNames)
				{
					if(dim[0] == nbCells) names = loom.getCellNames(); // depends of the metadata/dataset to extract
					else if(dim[0] == nbGenes) names = loom.getGeneHGNC();
					else new ErrorJSON("Not sure what to do here...");
				}
				
				// Retrieve indexes if exist
				if(Parameters.indexes == null) 
				{
					if(dim[1] == nbCells) Parameters.indexes = loom.getCellIndexes(Parameters.names); // depends of the metadata/dataset to extract
					else if(dim[1] == nbGenes) Parameters.indexes = loom.getGeneIndexes(Parameters.names);
					else new ErrorJSON("Not sure what to do here...");
				}
				
				// Creating output string
	    		sb = new StringBuilder();
	    		sb.append("{\"values\":[");
	    		prefix = "";
	    		if(dim[1] < 100) // TODO replace this by handling CONTIGUOUS / CHUNKED + Handle multiple read in same block
	    		{
	    			// Load the whole thing in memory
	    			float[][] matrix = loom.readFloatMatrix(Parameters.iAnnot);
	    			for(long index:Parameters.indexes)
					{
						sb.append(prefix);
						if(index == -1) sb.append("null");
						else
						{
							sb.append("[");
							String prefix2 = "";
							for(int i = 0; i < matrix.length; i++)
							{
								sb.append(prefix2).append(Utils.format(matrix[i][(int)index])); // TODO does not work if too big array
								prefix2 = ",";
							}
							sb.append("]");
						}
						prefix = ",";
					}
	    		}
	    		else
	    		{
					for(long index:Parameters.indexes)
					{
						sb.append(prefix);
						if(index == -1) sb.append("null");
						else
						{
							sb.append("[");
							float[] col = loom.readCol(index, Parameters.iAnnot); // TODO does not work if too big array // TODO handle multiple read in same block
							String prefix2 = "";
							for(int i = 0; i < col.length; i++)
							{
								sb.append(prefix2).append(Utils.format(col[i]));
								prefix2 = ",";
							}
							sb.append("]");
						}
						prefix = ",";
					}
	    		}
	    		sb.append("]");
				loom.close();

				// If names to be added
				if(Parameters.displayNames)
				{
					sb.append(",\"names\":[");
					prefix = "";
					for(String n:names)
					{
						sb.append(prefix).append("\"").append(n).append("\"");
						prefix = ",";
					}
					sb.append("]");
				}
				sb.append("}");
				
				// Writing results
		    	Utils.writeJSON(sb, Parameters.JSONFileName);
				break;
			case ListMetaData:
				// Get Metadata infos
				loom = new LoomFile("r", Parameters.fileName);
				sb = Metadata.toString(loom.listMetadata());
				loom.close();
				
	    		// Writing results
				Utils.writeJSON(sb, Parameters.JSONFileName);
				break;
			case ExtractMetaData:
				loom = new LoomFile("r", Parameters.loomFile);
		
				// Get all metadata to process
				String[] metapath;
				if(Parameters.metaName != null) 
				{
					metapath = new String[1];
					metapath[0] = Parameters.metaName;
				}
				else metapath = CopyMetaJSON.parseJSON(Parameters.JSONFileName);
				if(metapath == null) { loom.close(); new ErrorJSON("The JSON file should contain a meta field, here we don't detect any dataset?"); }
				
				// Prepare the output string
				prefix = "";
				sb = new StringBuilder();
				if(Parameters.metaName == null) sb.append("{\"list_meta\":[");
				
				// Process all metadata
				boolean colDisplayed = false;
				boolean rowDisplayed = false;
				for(String metaToExtract:metapath)
				{
					// Extract metadata from Loom
					Metadata metadata = loom.fillInfoMetadata(metaToExtract, true); // We could in principle simplify this part if we know already the data type / if we want to output the values / etc...
					if(Parameters.metaName == null && Parameters.metatype != null) metadata.type = Parameters.metatype;
					
					sb.append(prefix);
					
					StringArray64 feature_names = null;
					if(Parameters.displayNames && metaToExtract.startsWith("/col_attrs") && !colDisplayed) { feature_names = loom.getCellNames(); colDisplayed = true; } 
					if(Parameters.displayNames && metaToExtract.startsWith("/row_attrs") && !rowDisplayed) { feature_names = loom.getGeneHGNC(); rowDisplayed = true; }
					
					metadata.addMeta(sb, Parameters.displayValues, feature_names, Long.MAX_VALUE);
					
					prefix = ",";
				}
				if(Parameters.metaName == null) sb.append("]");
				
				// Extract cell names if needed // TODO IS NOT WITHIN THE MAIN JSON :/
				/*if(Parameters.metaName != null && Parameters.displayNames) 
				{
					Metadata.addCellNames(sb, prefix, loom.getCellNames(), Long.MAX_VALUE);
					Metadata.addGeneNames(sb, prefix, loom.getGeneHGNC(), Long.MAX_VALUE);
				}*/
					
				// Close Loom file
				if(Parameters.metaName == null) sb.append("}");
				loom.close();
				
	    		// Writing results
				Utils.writeJSON(sb, Parameters.outputFile);
				break;
			case MatchValues:
				loom = new LoomFile("r", Parameters.loomFile);

				// Extract indexes where metadata value matches the input, from Loom
				ArrayList<Long> indexes = loom.getIndexesWhereValueIs(Parameters.iAnnot, Parameters.value);
				loom.close();
								
				// Prepare the output string
				sb = new StringBuilder();
		    	if(Parameters.displayValues)
		    	{
					sb.append("{\"indexes_match\":[");
					prefix = "";
					for(int i=0; i < indexes.size(); i++) 
					{
						sb.append(prefix).append(indexes.get(i));
						prefix = ",";
					}
					sb.append("]}");
		    	}
		    	else
		    	{
		    		sb.append("{\"indexes_match_count\":").append(indexes.size()).append("}");
		    	}

	    		// Writing results
				Utils.writeJSON(sb, Parameters.JSONFileName);
				break;
			case RemoveMetaData:
				// Open in writing mode and remove the metadata
				loom = new LoomFile("r+", Parameters.loomFile);
				loom.removeMetadata(Parameters.metaName);
				loom.close();

				// Create output.json
				Metadata meta = new Metadata(Parameters.metaName);
				
				// Create output String
				sb = new StringBuilder();
				meta.addMeta(sb, false, null, Long.MAX_VALUE);
				
	    		// Writing results
				Utils.writeJSON(sb, Parameters.outputFolder + "output.json");
				break;
			case CreateCellSelection:		
				// Read Selected Cells in JSON
				long[] selectedIndexes = FilterCellsJSON.parseJSON(Parameters.JSONFileName).selected_cells; // TODO if too many are selected, this does not work
				HashSet<Long> selectedStableIDs = new HashSet<>();
				for(long id:selectedIndexes) selectedStableIDs.add(id);
				
				// Create Metadata and write in Loom file
				loom = new LoomFile("r+", Parameters.loomFile);
				LongArray64 stable_ids = loom.readLongArray("/col_attrs/_StableID");
				IntArray64 newMetadata = new IntArray64(loom.getDimensions()[1]);
				for(long i = 0; i < stable_ids.size(); i++) {
					if(selectedStableIDs.remove(stable_ids.get(i))) newMetadata.set(i, 1);
				}
				if(selectedStableIDs.size() > 0) new ErrorJSON("There are " + selectedStableIDs.size() + " ids that are not found.");
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
				Utils.writeJSON(sb, Parameters.outputFile);
				break;
			case CopyMetaData:
				LoomFile loomFrom = new LoomFile("r", Parameters.loomFile);
				LoomFile loomTo = new LoomFile("r+", Parameters.loomFile2);

				if(Parameters.metaName != null) 
				{
					metapath = new String[1];
					metapath[0] = Parameters.metaName;
				}
				else metapath = CopyMetaJSON.parseJSON(Parameters.JSONFileName);
				
				// Create output String
				sb = new StringBuilder("{\"meta\":[");
				prefix = "";
				for(String metaToCopy:metapath) 
				{
					if(LoomFile.copyMetadata(metaToCopy, loomFrom, loomTo))
					{
						sb.append(prefix).append("\"").append(metaToCopy).append("\"");
						prefix = ",";
					}
				}
				sb.append("]}");

				loomFrom.close();
				loomTo.close();

	    		// Writing results
				Utils.writeJSON(sb, Parameters.outputFile);
				
				break;
			case UpdateEnsemblDB:
				EnsemblDB.updateDB();
				break;
			case RegenerateNewOrganism:
				new ErrorJSON("Not implemented yet");
				//RegenerateNewOrganism.regenerateJSON();
				break;
			case FilterCells: 
				FileFilter.filterCells();
				break;
			case FilterGenes: 
				FileFilter.filterGenes();
				break;
			case FilterDEMetadata:
				FileFilter.filterDEMetadata();
				break;
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
							case "CreateKeggDB": m = Mode.CreateKeggDB; break;
							case "CreateGODB": m = Mode.CreateGODB; break;
							case "Enrichment": m = Mode.Enrichment; break;
							case "Parsing": m = Mode.Parsing; break;
							case "Preparsing": m = Mode.Preparsing; break;
							case "DimensionReduction": m = Mode.DimensionReduction; break;
							case "Normalization": m = Mode.Normalization; break;
							case "Scaling": m = Mode.Scaling; break;
							case "UpdateEnsemblDB": m = Mode.UpdateEnsemblDB; break;
							case "RegenerateNewOrganism": m = Mode.RegenerateNewOrganism; break;
							case "ExtractRow": m = Mode.ExtractRow; break;
							case "ExtractCol": m = Mode.ExtractCol; break;
							case "ListMetadata": m = Mode.ListMetaData; break;
							case "ExtractMetadata": m = Mode.ExtractMetaData; break;
							case "MatchValues": m = Mode.MatchValues; break;
							case "RemoveMetaData": m = Mode.RemoveMetaData; break;
							case "CopyMetaData": m = Mode.CopyMetaData; break;
							case "CreateCellSelection": m = Mode.CreateCellSelection; break;
							case "PreparseMetadata": m = Mode.PreparseMetadata; break;
							case "ParseMetadata": m = Mode.ParseMetadata; break;
							case "FilterCols": m = Mode.FilterCells; break;
							case "FilterRows": m = Mode.FilterGenes; break;
							case "FilterDEMetadata": m = Mode.FilterDEMetadata; break;
							case "DifferentialExpression": m = Mode.DifferentialExpression; break;
							default: System.err.println("Mode (-T) " + mode + " does not exist!"); System.out.println("-T %s \t\tMode to run ASAP [Preparsing, Parsing, PreparseMetadata, CreateCellSelection, ParseMetadata, RegenerateOutput, CreateKeggDB, CreateGODB, DifferentialExpression, Normalization, Scaling, UpdateEnsemblDB, Enrichment, ExtractRow, ExtractCol, ListMetadata, ExtractMetadata, MatchValues, RemoveMetaData, CopyMetaData, FilterCols, FilterRows, FilterDEMetadata]."); System.exit(-1);
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
			System.out.println("-T %s \t\tMode to run ASAP [Preparsing, Parsing, PreparseMetadata, CreateCellSelection, ParseMetadata, RegenerateOutput, CreateKeggDB, CreateGODB, DifferentialExpression, Normalization, Scaling, UpdateEnsemblDB, Enrichment, ExtractRow, ExtractCol, ListMetadata, ExtractMetadata, MatchValues, RemoveMetaData, CopyMetaData, FilterCols, FilterRows, FilterDEMetadata].");
			System.exit(-1);
		}
		return args2;
	}
}
