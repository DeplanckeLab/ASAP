import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import bigarrays.DoubleArray64;
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
import model.AnnotationIndex;
import model.MetaOn;
import model.Metadata;
import model.Metatype;
import model.Mode;
import model.Parameters;
import model.Parameters.OutputType;
import model.ResultStats;
import module_score.ModuleScore;
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
		DBManager.JDBC_DRIVER = Config.driver;
		DBManager.URL = Config.ConfigMAIN().getURL("asap2_data_v5");
		
		if(args.length == 0) readMode(args);
		
		// Check if debug mode
		String[] args2 = isDebug(args);
		if(Parameters.debugMode) DBManager.URL = Config.ConfigDEV().getURL("asap2_data_v5");
		
		// Check which tool is called
		args2 = readMode(args2);
		
		// Load associated parameters
		Parameters.load(args2, m);
		
		Utils.initRandomGenerator();
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
			case MarkerEnrichment:
				Enrichment.runMarkerEnrichment();
				break;
			case ModuleScore:
				ModuleScore.runModuleScore();
				break;
			case DimensionReduction: 
				FileDimReduc.reduceDimension();
				break;
			case DifferentialExpression: 
				DE.performDE();
				break;
			case FindMarkers:
				DE.findMarkers();
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
			case IndexByCell:			
				// Recuperate the indexes of all categories for this metadata
				if(Parameters.debugMode) DBManager.URL = Config.ConfigDEV().getURL("asap2_development"); // Annotations/Metadata are not in same DB
				else DBManager.URL = Config.ConfigMAIN().getURL("asap2_development");
				DBManager.connect();
				HashMap<String, Integer> categories = DBManager.getListCatJSON(Parameters.id);
				DBManager.disconnect();
				
				// Reading Loom with metadata to index
				LoomFile loom = new LoomFile("r", Parameters.loomFile);
				LongArray64 stable_ids = loom.getCellStableIds();
				if(stable_ids == null) new ErrorJSON("No StableID in this Loom file: " + Parameters.loomFile);
										
				// Extract metadata from Loom
				Metadata metadata = loom.fillInfoMetadata(Parameters.iAnnot, true);
				if(metadata.values.size() != loom.getNbCells()) new ErrorJSON("Nb of cells do not match this metadata");
				loom.close();
				
				// Opening Loom file containing the cell indexes (Can be the same file, that's why I close it first. I could do a test to not close it if it's the same file.... but meh...)
				if(!new File(Parameters.loomFile2).exists()) new ErrorJSON("This loom file used for storing the indexes, does not exist: " + Parameters.loomFile2);
				LoomFile loomIndex = new LoomFile("r+", Parameters.loomFile2);
				LongArray64 stable_ids_REF = loomIndex.getCellStableIds();
				if(stable_ids_REF == null) new ErrorJSON("No StableID in this Loom file: " + Parameters.loomFile2);
				
				// Extract existing indexes, or create it
				StringArray64 indexesByCell = null;
				if(!loomIndex.exists("/col_attrs/_INDEX_ASAP_CELLS")) indexesByCell = new StringArray64(stable_ids_REF.size());
				else indexesByCell = loomIndex.readStringArray("/col_attrs/_INDEX_ASAP_CELLS");
				
				// Add indexes
				long idx_to_add = 0;
				long sid_to_add = stable_ids.get(idx_to_add);
				for(int i = 0; i < indexesByCell.size(); i++)
				{
					long stable_id = stable_ids_REF.get(i);
					if(stable_id == sid_to_add) // Here I suppose that all stable_ids are sorted. If this is not the case, this will fail.... See error after the loop
					{
						String idx = indexesByCell.get(i);
						StringBuffer sb = new StringBuffer(idx);
						if(!idx.equals("")) sb.append(",");
						Integer val = categories.get(metadata.values.get(idx_to_add));
						if(val == null) new ErrorJSON("This category: " + metadata.values.get(idx_to_add) + " was not found in the DB, for metadata id = " + Parameters.id);
						sb.append(Parameters.id).append(":").append(val);
						indexesByCell.set(i, sb.toString());
						idx_to_add++;
						if(idx_to_add == stable_ids.size()) break;
						sid_to_add = stable_ids.get(idx_to_add);
					}
				}
				if(idx_to_add != stable_ids.size()) new ErrorJSON("Could not find all indexes? Stable_ids are not sorted?");
				
				// Write results in Loom
				loomIndex.writeStringArray("/col_attrs/_INDEX_ASAP_CELLS", indexesByCell);	
				loomIndex.close();
				
				// Creating output string
	    		StringBuilder sb = new StringBuilder();
	    		sb.append("{\"nb_cells_indexed\":").append(idx_to_add).append("}");
	    		
				// Writing results
		    	Utils.writeJSON(sb);

				break;
			case GetIndex:
				// Handle correct DB
				if(Parameters.debugMode) DBManager.URL = Config.ConfigDEV().getURL("asap2_development"); // Annotations/Metadata are not in same DB
				else DBManager.URL = Config.ConfigMAIN().getURL("asap2_development");
				
				// First, recuperate the list of cells (stable_ids)
				long[] cellsToCheck = FilterCellsJSON.parseJSON(Parameters.JSONFileName).selected_cells; // Stable_ids
				HashSet<Long> cellsToCheckMap = new HashSet<>();
				for(Long s:cellsToCheck) cellsToCheckMap.add(s);

				// Opening Loom file containing the cell indexes
				loomIndex = new LoomFile("r", Parameters.loomFile);
				stable_ids_REF = loomIndex.getCellStableIds();
				if(stable_ids_REF == null) new ErrorJSON("No StableID in this Loom file: " + Parameters.loomFile);

				// Extract existing indexes
				indexesByCell = null;
				if(!loomIndex.exists("/col_attrs/_INDEX_ASAP_CELLS")) new ErrorJSON("'/col_attrs/_INDEX_ASAP_CELLS' is not present in loom " + Parameters.loomFile);
				indexesByCell = loomIndex.readStringArray("/col_attrs/_INDEX_ASAP_CELLS");
				loomIndex.close();
				
				// Now compute the stats
				HashMap<String, AnnotationIndex> counts = new HashMap<String, AnnotationIndex>(); // Results
				DBManager.connect(); // Connect to the DB to retrieve the categories
				for(long i = 0; i < indexesByCell.size(); i++)
				{
					long stable_id = stable_ids_REF.get(i);
					if(cellsToCheckMap.remove(stable_id)) // It's a cell to consider
					{
						String row = indexesByCell.get(i);
						// Parsing using fastest method possible (split is slow as hell)
						
						int start = 0;
						HashSet<String> safetyCheck = new HashSet<String>();
						while(start < row.length())
						{
							int columnSep = row.indexOf(':', start);
							int commaSep = row.indexOf(',', columnSep);
							if(commaSep == -1) commaSep = row.length();
							String id = row.substring(start, columnSep);
							if(safetyCheck.contains(id)) new ErrorJSON("Cell with stable_id " + stable_id + " has several time the metadata index " + id);
							safetyCheck.add(id);
							String cat = row.substring(columnSep + 1, commaSep);
							
							AnnotationIndex ai = counts.get(id);
							if(ai == null) ai = new AnnotationIndex(id);
							ai.add(cat);
							counts.put(id, ai);
							start = commaSep + 1;
						}
					}
				}
				DBManager.disconnect();
				
				// Safety check
				if(cellsToCheckMap.size() > 0) new ErrorJSON("Some cells stable_ids are not in the Loom index file: " + Utils.toString(cellsToCheckMap));
				
				// Write output JSON
				sb = new StringBuilder("{");
				String prefix = "";
				for(String metaId:counts.keySet())
				{
					sb.append(prefix).append("\"").append(metaId).append("\":[");
					AnnotationIndex ai = counts.get(metaId);
					String prefix2 = "";
					for(int i = 0; i < ai.categories.length; i++)
					{
						sb.append(prefix2).append(ai.categories[i]);
						prefix2 = ",";
					}
					sb.append("]");
					prefix = ",";
				}		
				sb.append("}");
				
				// Writing results
		    	Utils.writeJSON(sb);
				
				break;
			case GetGeneStats:
				// Cells to extract values from
				LongArray64 cellStableIds = null;		
				if(Parameters.loom_cell_stable_ids != null) {
					loom = new LoomFile("r", Parameters.loom_cell_stable_ids);
					if(!loom.exists("/col_attrs/_StableID")) new ErrorJSON("No StableID in this Loom file: " + Parameters.loom_cell_stable_ids);
					cellStableIds = loom.readLongArray("/col_attrs/_StableID");
					loom.close();
				}
				
				// Alternatively, using the metadata
				LongArray64 cell_indexes = null; // Limit to these cells
				loom = new LoomFile("r", Parameters.loomFile);
				long foundcells = 0;
				if(Parameters.id != -1) // Then Parameters.index is also set
				{
					// Handle correct DB
					if(Parameters.debugMode) DBManager.URL = Config.ConfigDEV().getURL("asap2_development"); // Annotations/Metadata are not in same DB
					else DBManager.URL = Config.ConfigMAIN().getURL("asap2_development");
									
					// Recuperate the indexes of all categories for this metadata
					DBManager.connect();
					String[] res = DBManager.getMetadataCategorie(Parameters.id, (int)Parameters.index);
					Parameters.metaName = res[0];
					Parameters.value = res[1];
					DBManager.disconnect();
					
					// Extract indexes where metadata value matches the input, from Loom
					cell_indexes = LongArray64.convertFrom(loom.getIndexesWhereValueIs(Parameters.metaName, Parameters.value));
					foundcells = cell_indexes.size();
				} 
				else if(cellStableIds != null) 
				{
					cell_indexes = loom.getCellIndexesByStableIds(cellStableIds);
					for(long i = 0; i < cell_indexes.size(); i++) if(cell_indexes.get(i) != -1) foundcells++;
				}
				
				// First check the total dimension of the dataset to extract from
				long[] dim = loom.getDimensions();
				long nbGenes = dim[0];
				long nbCells = dim[1];
				boolean toNormalize = false;
				if(Parameters.iAnnot == null) { toNormalize = true; Parameters.iAnnot = "/matrix"; }
				dim = loom.getDimensions(Parameters.iAnnot); // dim[1] = nb of values to extract for each row
				if(dim[0] != nbGenes && dim[1] != nbCells) { loom.close(); new ErrorJSON("You cannot use GetGeneStats on this dataset: " + Parameters.iAnnot); }
				
		    	// Recuperate depth for normalizing the data
				DoubleArray64 depth = null;
		    	if(toNormalize)
		    	{
					if(!loom.exists("/col_attrs/_Depth")) new ErrorJSON("Error in the Loom file. Path '/col_attrs/_Depth' should exist!");
			    	depth = loom.readDoubleArray("/col_attrs/_Depth");
		    	}
		    	
				// Retrieve indexes if exist
				if(Parameters.indexes == null) 
				{
					if(Parameters.names != null) Parameters.indexes = loom.getGeneIndexes(Parameters.names);
					else if(Parameters.stable_ids != null) Parameters.indexes = loom.getGeneIndexesByStableIds(Parameters.stable_ids);
					// the other case should be handled in Parameters
				}
				
				// Creating output result struct
				ResultStats[] results = new ResultStats[Parameters.indexes.length];
				ResultStats[] results_comp = new ResultStats[Parameters.indexes.length]; // Complementary results
				
	    		// Read the values and compute stats
				for (int i = 0; i < results.length; i++)
				{
					long index = Parameters.indexes[i];
					if(index == -1) 
					{
						results[i] = null;
						results_comp[i] = null;
					}
					else
					{
						// Read the gene expressions
						float[] row = loom.readRow(index, Parameters.iAnnot); // TODO handle multiple read in same block
						
						// In case we need to normalize (by column)
						if(toNormalize)
						{
							for (int j = 0; j < row.length; j++) 
							{
								row[j] = (float)(row[j] / depth.get(j)) * Parameters.scale_factor;
							}
						}
						
						// In case some cells are filtered
						if(cell_indexes != null) 
						{
							// 0. Prepare indexes
							int[] indexes_found = new int[row.length]; // All 0
							
							// 1. Main cells
							float[] filtered_row = new float[(int)foundcells];
							int index_found_cells = 0;
							for(long ci = 0; ci < cell_indexes.size(); ci++)
							{
								long cellI = cell_indexes.get(ci);
								if(cellI != -1) 
								{
									filtered_row[index_found_cells] = row[(int)cellI];
									index_found_cells++;
									indexes_found[(int)cellI] = -1; // This one is already counted (for comp.)
								}
							}
							
							// Compute stats
							ResultStats res = new ResultStats();
							Arrays.sort(filtered_row);
							res.min = filtered_row[0];
							res.median = Utils.median(filtered_row, true);
							res.mean = Utils.mean(filtered_row);
							res.q1 = Utils.q1(filtered_row, true);
							res.q3 = Utils.q3(filtered_row, true);
							res.max = filtered_row[filtered_row.length - 1];
							if(toNormalize) res.log2p(); // log2(1 + x) if normalized on-the-go (AFTER computing the mean)
							results[i] = res;
							
							// 2. Complementary cells
							filtered_row = new float[row.length - (int)foundcells];
							index_found_cells = 0;
							for(int j = 0; j < indexes_found.length; j++)
							{
								if(indexes_found[j] != -1) // If not already counted for the Main cells loop
								{
									filtered_row[index_found_cells] = row[j];
									index_found_cells++;
								}
							}
							
							// Compute stats
							res = new ResultStats();
							Arrays.sort(filtered_row);
							res.min = filtered_row[0];
							res.median = Utils.median(filtered_row, true);
							res.mean = Utils.mean(filtered_row);
							res.q1 = Utils.q1(filtered_row, true);
							res.q3 = Utils.q3(filtered_row, true);
							res.max = filtered_row[filtered_row.length - 1];
							if(toNormalize) res.log2p(); // log2(1 + x) if normalized on-the-go (AFTER computing the mean)						
							results_comp[i] = res;
						}
						else
						{	
							// Compute stats
							ResultStats res = new ResultStats();
							Arrays.sort(row);
							res.min = row[0];
							res.median = Utils.median(row, true);
							res.mean = Utils.mean(row);
							res.q1 = Utils.q1(row, true);
							res.q3 = Utils.q3(row, true);
							res.max = row[row.length - 1];
							if(toNormalize) res.log2p(); // log2(1 + x) if normalized on-the-go (AFTER computing the mean)
							results[i] = res;
							results_comp[i] = null; // No cell filtering
						}
					}
				}
				
				// Close handle
				loom.close();
				
				// Create output string
				sb = new StringBuilder("{");
				prefix = "";
				if(Parameters.names != null)
				{
					sb.append("\"names\":[");
					for(int i = 0; i < Parameters.names.length; i++)
					{
						sb.append(prefix).append(Parameters.names[i]);
						prefix = ",";
					}
				}
				else if(Parameters.stable_ids != null)
				{
					sb.append("\"stable_ids\":[");
					for(int i = 0; i < Parameters.stable_ids.length; i++)
					{
						sb.append(prefix).append(Parameters.stable_ids[i]);
						prefix = ",";
					}
				}
				else
				{
					sb.append("\"indexes\":[");
					for(int i = 0; i < Parameters.indexes.length; i++)
					{
						sb.append(prefix).append(Parameters.indexes[i]);
						prefix = ",";
					}
				}
				
				sb.append("], \"min\":[");
				prefix = "";
				for(int i = 0; i < results.length; i++) 
				{ 
					if(results[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results[i].min); 
					prefix = ",";
				}
				
				sb.append("], \"q1\":[");
				prefix = "";
				for(int i = 0; i < results.length; i++) 
				{ 
					if(results[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results[i].q1); 
					prefix = ",";
				}
				
				sb.append("], \"median\":[");
				prefix = "";
				for(int i = 0; i < results.length; i++) 
				{ 
					if(results[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results[i].median); 
					prefix = ",";
				}

				sb.append("], \"mean\":[");
				prefix = "";
				for(int i = 0; i < results.length; i++) 
				{ 
					if(results[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results[i].mean); 
					prefix = ",";
				}

				sb.append("], \"q3\":[");
				prefix = "";
				for(int i = 0; i < results.length; i++) 
				{ 
					if(results[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results[i].q3); 
					prefix = ",";
				}

				sb.append("], \"max\":[");
				prefix = "";
				for(int i = 0; i < results.length; i++) 
				{ 
					if(results[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results[i].max); 
					prefix = ",";
				}
				
				sb.append("], \"min_comp\":[");
				prefix = "";
				for(int i = 0; i < results_comp.length; i++) 
				{ 
					if(results_comp[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results_comp[i].min); 
					prefix = ",";
				}
				
				sb.append("], \"q1_comp\":[");
				prefix = "";
				for(int i = 0; i < results_comp.length; i++) 
				{ 
					if(results_comp[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results_comp[i].q1); 
					prefix = ",";
				}
				
				sb.append("], \"median_comp\":[");
				prefix = "";
				for(int i = 0; i < results_comp.length; i++) 
				{ 
					if(results_comp[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results_comp[i].median); 
					prefix = ",";
				}

				sb.append("], \"mean_comp\":[");
				prefix = "";
				for(int i = 0; i < results_comp.length; i++) 
				{ 
					if(results_comp[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results_comp[i].mean); 
					prefix = ",";
				}

				sb.append("], \"q3_comp\":[");
				prefix = "";
				for(int i = 0; i < results_comp.length; i++) 
				{ 
					if(results_comp[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results_comp[i].q3); 
					prefix = ",";
				}

				sb.append("], \"max_comp\":[");
				prefix = "";
				for(int i = 0; i < results_comp.length; i++) 
				{ 
					if(results_comp[i] == null) sb.append(prefix).append("null"); 
					else sb.append(prefix).append(results_comp[i].max); 
					prefix = ",";
				}

				sb.append("]}");
				
				// Writing results
		    	Utils.writeJSON(sb, Parameters.JSONFileName);
				break;
			case ExtractRow:
				// Cells to extract values from
				cellStableIds = null;		
				if(Parameters.loom_cell_stable_ids != null) {
					loom = new LoomFile("r", Parameters.loom_cell_stable_ids);
					if(!loom.exists("/col_attrs/_StableID")) new ErrorJSON("No StableID in this Loom file: " + Parameters.loom_cell_stable_ids);
					cellStableIds = loom.readLongArray("/col_attrs/_StableID");
					loom.close();
				}
								
				// Reading Loom with data to extract
				loom = new LoomFile("r", Parameters.loomFile);
				cell_indexes = null; // Limit to these cells
				if(cellStableIds != null) cell_indexes = loom.getCellIndexesByStableIds(cellStableIds);
				
				// First check the total dimension of the dataset to extract from
				dim = loom.getDimensions();
				nbGenes = dim[0];
				nbCells = dim[1];
				dim = loom.getDimensions(Parameters.iAnnot); // dim[1] = nb of values to extract for each row
				if(dim.length == 1) { loom.close(); new ErrorJSON("The dataset is not a matrix, it's an array. You should use ExtractMetadata instead"); }
				if(dim[1] == nbGenes && cellStableIds != null) { loom.close(); new ErrorJSON("You cannot request this dataset and specify -loom_cells"); }
				
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
	    		sb = new StringBuilder();
	    		sb.append("{\"values\":[");
	    		
	    		// In case of sorting
	    		StringBuilder sb2 = new StringBuilder();
	    		if(Parameters.sort) sb2.append(",\"orders\":[");
	    		HashMap<Long, int[]> orders = new HashMap<Long, int[]>(); // Keep the orders in memory
	    		
	    		// Start reading the values
	    		prefix = "";
	    		if(dim[0] < 100) // TODO replace this by handling CONTIGUOUS / CHUNKED + Handle multiple read in same block
	    		{
	    			// Load the whole thing in memory
	    			float[][] matrix = loom.readFloatMatrix(Parameters.iAnnot);
	    			for(long index:Parameters.indexes)
					{
						sb.append(prefix);
						if(Parameters.sort) sb2.append(prefix);
						if(index == -1) 
						{
							sb.append("null");
							if(Parameters.sort) sb2.append("null");
						}
						else
						{
							sb.append("[");
							if(Parameters.sort) sb2.append("[");
							String prefix2 = "";
							if(cell_indexes != null) 
							{
								if(Parameters.sort) new ErrorJSON("Sort is not implemented in this case");
								for(long ci = 0; ci < cell_indexes.size(); ci++)
								{
									long cellI = cell_indexes.get(ci);
									if(cellI == -1) sb.append(prefix2).append("null");
									else sb.append(prefix2).append(Utils.format(matrix[(int)index][(int)cellI])); // TODO does not work if too big array
									prefix2 = ",";
								}
							}
							else 
							{
								float[] row = matrix[(int)index];
								if(Parameters.sort)
								{
									int[] order = Utils.order(row, false);
									orders.put(index, order);
									for(int i:order)
									{
										sb.append(prefix2).append(Utils.format(row[i]));
										sb2.append(prefix2).append(i);
										prefix2 = ",";
									}
								}
								else
								{
									for(int i = 0; i < row.length; i++)
									{
										sb.append(prefix2).append(Utils.format(row[i])); // TODO does not work if too big array
										prefix2 = ",";
									}
								}
							}
							sb.append("]");
							if(Parameters.sort) sb2.append("]");
						}
						prefix = ",";
					}
	    		}
	    		else
	    		{
					for(long index:Parameters.indexes)
					{
						sb.append(prefix);
						if(Parameters.sort) sb2.append(prefix);
						if(index == -1) 
						{
							if(Parameters.sort) sb2.append("null");
							sb.append("null");
						}
						else
						{
							sb.append("[");
							if(Parameters.sort) sb2.append("[");
							float[] row = loom.readRow(index, Parameters.iAnnot); // TODO does not work if too big array // TODO handle multiple read in same block
							String prefix2 = "";
							
							if(cell_indexes != null) 
							{
								if(Parameters.sort) new ErrorJSON("Sort is not implemented in this case");
								for(long ci = 0; ci < cell_indexes.size(); ci++)
								{
									long cellI = cell_indexes.get(ci);
									if(cellI == -1) sb.append(prefix2).append("null");
									else sb.append(prefix2).append(Utils.format(row[(int)cellI])); // TODO does not work if too big array
									prefix2 = ",";
								}
							}
							else 
							{
								if(Parameters.sort)
								{
									int[] order = Utils.order(row, false);
									orders.put(index, order);
									for(int i:order)
									{
										sb.append(prefix2).append(Utils.format(row[i]));
										sb2.append(prefix2).append(i);
										prefix2 = ",";
									}
								}
								else
								{
									for(int i = 0; i < row.length; i++)
									{
										sb.append(prefix2).append(Utils.format(row[i]));
										prefix2 = ",";
									}
								}
							}
							sb.append("]");
							if(Parameters.sort) sb2.append("]");
						}
						prefix = ",";
					}
					sb.append("]");
					if(Parameters.sort) sb2.append("]");
	    		}
				
				// In case we need to output the indexes as well 
				if(Parameters.sort) sb.append(sb2);
				
				// If names to be added
				if(Parameters.displayNames)
				{
					StringArray64 names = null;					
					if(dim[1] == nbCells) names = loom.getCellNames(); // depends of the metadata/dataset to extract
					else if(dim[1] == nbGenes) names = loom.getGeneHGNC();
					else new ErrorJSON("Not sure what to do here...");
					sb.append(",\"names\":[");
					prefix = "";
					if(cell_indexes != null) 
					{
						if(Parameters.sort) new ErrorJSON("Sort is not implemented in this case");
						for(long ci = 0; ci < cell_indexes.size(); ci++)
						{
							long cellI = cell_indexes.get(ci);
							if(cellI == -1) sb.append(prefix).append("null");
							else sb.append(prefix).append("\"").append(names.get(cellI)).append("\""); // TODO does not work if too big array
							prefix = ",";
						}
					}
					else 
					{
						// I don't sort the names....
						for(String n:names)
						{
							sb.append(prefix).append("\"").append(n).append("\"");
							prefix = ",";
						}
					}
					sb.append("]");
				}
				
				// If stable_ids to be added
				if(Parameters.export_stable_ids)
				{
					stable_ids = null;
					if(dim[1] == nbCells) stable_ids = loom.getCellStableIds(); // depends of the metadata/dataset to extract
					else if(dim[1] == nbGenes) stable_ids = loom.getGeneStableIds();
					else new ErrorJSON("Not sure what to do here...");
					sb.append(",\"stable_ids\":[");
					prefix = "";
					
					if(cell_indexes != null) 
					{
						if(Parameters.sort) new ErrorJSON("Sort is not implemented in this case");
						for(long ci = 0; ci < cell_indexes.size(); ci++)
						{
							long cellI = cell_indexes.get(ci);
							if(cellI == -1) sb.append(prefix).append("null");
							else sb.append(prefix).append(stable_ids.get(cellI)); // TODO does not work if too big array
							prefix = ",";
						}
					}
					else 
					{
						// I don't sort the ids....
						for(Long id:stable_ids)
						{
							sb.append(prefix).append(id);
							prefix = ",";
						}
					}
					
					sb.append("]");
				}
				sb.append("}");
				
				// Close handle
				loom.close();
				
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

				// If names to be added
				if(Parameters.displayNames)
				{
					StringArray64 names = null;
					if(dim[0] == nbCells) names = loom.getCellNames(); // depends of the metadata/dataset to extract
					else if(dim[0] == nbGenes) names = loom.getGeneHGNC();
					else new ErrorJSON("Not sure what to do here...");
					sb.append(",\"names\":[");
					prefix = "";
					for(String n:names)
					{
						sb.append(prefix).append("\"").append(n).append("\"");
						prefix = ",";
					}
					sb.append("]");
				}
				
				// If stable_ids to be added
				if(Parameters.export_stable_ids)
				{
					stable_ids = null;
					if(dim[0] == nbCells) stable_ids = loom.getCellStableIds(); // depends of the metadata/dataset to extract
					else if(dim[0] == nbGenes) stable_ids = loom.getGeneStableIds();
					else new ErrorJSON("Not sure what to do here...");
					sb.append(",\"stable_ids\":[");
					prefix = "";
					for(Long id:stable_ids)
					{
						sb.append(prefix).append(id);
						prefix = ",";
					}
					sb.append("]");
				}
				sb.append("}");
				
				// Close Handle
				loom.close();
				
				// Writing results
		    	Utils.writeJSON(sb, Parameters.JSONFileName);
				break;
			case ExtractDataset:
				HashSet<Long> indexesToExtract = null;
				HashSet<Long> stableIdsToExtract = null;
				
				// Check if valid dataset
				if(!Parameters.iAnnot.startsWith("/matrix") && !Parameters.iAnnot.startsWith("/layers/")) new ErrorJSON("Dataset should be /matrix or starts with /layers/");
				
				// Loom to extract stable_id from
				if(Parameters.indexes != null)
				{
					if(Parameters.loom_cell_stable_ids.equals(Parameters.loomFile)) 
					{
						indexesToExtract = new HashSet<Long>();
						for(long l:Parameters.indexes) indexesToExtract.add(l);
					}
				}
				
				//TODO add a case with --stable-ids and same Loom?
				// Read stable_ids to extract in other Loom file
				if(indexesToExtract == null && (Parameters.indexes != null || Parameters.stable_ids != null)) 
				{
					loom = new LoomFile("r", Parameters.loom_cell_stable_ids);
					if(!loom.exists("/col_attrs/_StableID")) new ErrorJSON("No StableID in this Loom file: " + Parameters.loom_cell_stable_ids);
					cellStableIds = loom.readLongArray("/col_attrs/_StableID");
					loom.close();
					
					stableIdsToExtract = new HashSet<Long>();
					if(Parameters.indexes != null)
					{
						for(int i = 0 ; i < Parameters.indexes.length; i++) stableIdsToExtract.add(cellStableIds.get(Parameters.indexes[i]));
					}
					else if(Parameters.stable_ids != null)
					{
						for(int i = 0 ; i < Parameters.stable_ids.length; i++) stableIdsToExtract.add(Parameters.stable_ids[i]);
					}
				}
				
				// Loom to extract data from
				loom = new LoomFile("r", Parameters.loomFile);

				// First check the total dimension of the dataset to extract from
				dim = loom.getDimensions(Parameters.iAnnot);
				if(dim.length == 1) { loom.close(); new ErrorJSON("The dataset is not a matrix, it's an array. It should NOT BE THE CASE!"); }
				nbGenes = dim[0];
				nbCells = dim[1];
				
				// Find indexes to extract in new Loom file
				if(indexesToExtract == null && stableIdsToExtract != null)
				{
					if(!loom.exists("/col_attrs/_StableID")) new ErrorJSON("No StableID in this Loom file: " + Parameters.loomFile);
					cellStableIds = loom.readLongArray("/col_attrs/_StableID");
					
					indexesToExtract = new HashSet<Long>();
					for(long i = 0 ; i < cellStableIds.size(); i++) 
					{
						long stable_id = cellStableIds.get(i);
						if(stableIdsToExtract.contains(stable_id)) indexesToExtract.add(i);
					}
				}
				
				// All values
				if(indexesToExtract == null)
				{
					indexesToExtract = new HashSet<Long>();
					for(long i = 0; i < nbCells; i++) indexesToExtract.add(i);
				}
				
				// Creating output string
	    		sb = new StringBuilder(); // TODO don't create a StringBuilder, it takes too much space in RAM
	    		
				// If names to be added
	    		StringArray64 rownames = null;
	    		StringArray64 colnames = null;
				if(Parameters.displayRowNames) rownames = loom.getGeneHGNC();
				if(Parameters.displayColNames) colnames = loom.getCellNames();
				
				if(Parameters.outputType == OutputType.JSON) 
				{
					sb.append("{");
					if(Parameters.displayRowNames)
					{
						sb.append("\"row_names\":[");
						prefix = "";
						for(String n:rownames)
						{
							sb.append(prefix).append("\"").append(n).append("\"");
							prefix = ",";
						}
						sb.append("],");
					}
					if(Parameters.displayColNames)
					{
						sb.append("\"col_names\":[");
						prefix = "";
						for(String n:colnames)
						{
							sb.append(prefix).append("\"").append(n).append("\"");
							prefix = ",";
						}
						sb.append("],");
					}
					sb.append("\"values\":[");
				}
				else
				{
					if(Parameters.displayColNames)
					{
						prefix = "";
						for(String n:colnames)
						{
							sb.append(prefix).append(n);
							prefix = ",";
						}
						sb.append("\n");
					}
				}

    			// Load the whole thing in memory (not optimal MAN) // TODO don't....
    			float[][] matrix = loom.readFloatMatrix(Parameters.iAnnot);
    			prefix = "";
    			for(long i = 0; i < matrix.length; i++) // genes
				{
					sb.append(prefix);
					String prefix2 = "";
					if(Parameters.outputType == OutputType.JSON) sb.append("[");
					else if(Parameters.displayRowNames) 
					{
						sb.append(rownames.get(i));
						prefix2 = ",";
					}
					for(long j = 0; j < matrix[(int)i].length; j++) // cells
					{
						if(indexesToExtract.contains(j))
						{
							sb.append(prefix2).append(Utils.format(matrix[(int)i][(int)j]));
							prefix2 = ",";
						}
					}
					if(Parameters.outputType == OutputType.JSON) 
					{
						sb.append("]");
						prefix = ",";
					}
					else sb.append("\n");				
				}
	    		if(Parameters.outputType == OutputType.JSON) sb.append("]");
				
				// Close Handle
				loom.close();
				
				// Writing results
		    	if(Parameters.outputType == OutputType.JSON) Utils.writeJSON(sb, Parameters.outputFile);
		    	else Utils.writePlain(sb, Parameters.outputFile);
				break;
			case ListMetadata:
				// Get Metadata infos
				loom = new LoomFile("r", Parameters.fileName);
				sb = new StringBuilder("{").append(Metadata.toString(loom.listMetadata())).append("}");
				loom.close();
				
	    		// Writing results
				Utils.writeJSON(sb, Parameters.JSONFileName);
				break;
			case ExtractMetadata:
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
					metadata = loom.fillInfoMetadata(metaToExtract, true); // We could in principle simplify this part if we know already the data type / if we want to output the values / etc...
					//if(Parameters.metaName != null && Parameters.metatype != null) metadata.type = Parameters.metatype;
					if(Parameters.metatype != null) metadata.type = Parameters.metatype;
					
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
				// Check which metadata(s) to remove
				if(Parameters.metaName != null) 
				{
					metapath = new String[1];
					metapath[0] = Parameters.metaName;
				}
				else metapath = CopyMetaJSON.parseJSON(Parameters.JSONFileName);
				
				// Open in writing mode and remove the metadata
				loom = new LoomFile("r+", Parameters.loomFile);
				
				// Go through metadata to remove
				sb = new StringBuilder("{\"meta\":[");
				prefix = "";
				for(String metaToRemove:metapath) 
				{
					loom.removeMetadata(metaToRemove);
					sb.append(prefix).append("\"").append(metaToRemove).append("\"");
					prefix = ",";
				}
				sb.append("]}");

				// Close Loom file
				loom.close();
				
	    		// Writing results
				Utils.writeJSON(sb, Parameters.outputFile);
				break;
			case RenameMetaData:			
				// Open in writing mode and remove the metadata
				loom = new LoomFile("r+", Parameters.loomFile);
				
				// Rename Metadata
				loom.renameMetadata(Parameters.metaName, Parameters.metaName2);
				
				// Close Loom file
				loom.close();
				
				// Preparing results
				sb = new StringBuilder("{\"meta-from\":");
				sb.append("\"").append(Parameters.metaName).append("\",\"meta-to\":\"").append(Parameters.metaName2).append("\"");
				sb.append("}");

	    		// Writing results
				Utils.writeJSON(sb, Parameters.outputFile);
				break;
			case CreateCellSelection:		
				// Read Selected Cells in JSON
				long[] selectedIndexes = FilterCellsJSON.parseJSON(Parameters.JSONFileName).selected_cells; // TODO if too many are selected, this does not work
				HashSet<Long> selectedStableIDs = new HashSet<>();
				for(long id:selectedIndexes) selectedStableIDs.add(id);
				
				// Create Metadata and write in Loom file
				loom = new LoomFile("r+", Parameters.loomFile);
				stable_ids = loom.readLongArray("/col_attrs/_StableID");
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
				Metadata meta = new Metadata(Parameters.metaName, Metatype.DISCRETE, MetaOn.CELL, newMetadata.size(), 1, cat);
				
				// Create output String
				sb = new StringBuilder();
				meta.addMeta(sb, false, null, Long.MAX_VALUE);

	    		// Writing results
				Utils.writeJSON(sb, Parameters.outputFile);
				break;
			case CopyMetaData:
				LoomFile loomFrom = new LoomFile("r", Parameters.loomFile);
				LoomFile loomTo = new LoomFile("r+", Parameters.loomFile2);

				// Checking dimensions of the files (in case the To file is smaller than the From file)
				long[] dimFrom = loomFrom.getDimensions();
				long[] dimTo = loomTo.getDimensions();
				HashSet<Long> geneIndexesToFilter = null;
				HashSet<Long> cellIndexesToFilter = null;
				if(dimFrom[0] < dimTo[0] || dimFrom[1] < dimTo[1]) new ErrorJSON("Cannot copy from a filtered dataset to a non-filtered one");
				if(dimFrom[0] > dimTo[0]) // Gene Filtering
				{
					LongArray64 stable_ids_from = loomFrom.readLongArray("/row_attrs/_StableID");
					LongArray64 stable_ids_to = loomTo.readLongArray("/row_attrs/_StableID");
					geneIndexesToFilter = new HashSet<Long>();
					long j = 0;
					for(long i = 0; i < stable_ids_to.size(); i ++) // I assume they are sorted... hopefully it is always the case
					{
						long s_to = stable_ids_to.get(i);
						long s_from = stable_ids_from.get(j);
						while(s_from != s_to) 
						{
							geneIndexesToFilter.add(j);
							s_from = stable_ids_from.get(++j);
						}
						j++;
					}
					while(j < stable_ids_from.size()) geneIndexesToFilter.add(j++);
					if(geneIndexesToFilter.size() != (stable_ids_from.size() - stable_ids_to.size())) new ErrorJSON("Could not find all gene stable_ids...");
				}
				if(dimFrom[1] > dimTo[1]) // Cell Filtering
				{
					LongArray64 stable_ids_from = loomFrom.readLongArray("/col_attrs/_StableID");
					LongArray64 stable_ids_to = loomTo.readLongArray("/col_attrs/_StableID");
					cellIndexesToFilter = new HashSet<Long>();
					long j = 0;
					for(long i = 0; i < stable_ids_to.size(); i ++) // I assume they are sorted... hopefully it is always the case
					{
						long s_to = stable_ids_to.get(i);
						long s_from = stable_ids_from.get(j);
						while(s_from != s_to) 
						{
							cellIndexesToFilter.add(j);
							s_from = stable_ids_from.get(++j);
						}
						j++;
					}
					while(j < stable_ids_from.size()) cellIndexesToFilter.add(++j);
					if(cellIndexesToFilter.size() != (stable_ids_from.size() - stable_ids_to.size())) new ErrorJSON("Could not find all cell stable_ids...");
				}
				
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
					if(LoomFile.copyMetadata(metaToCopy, geneIndexesToFilter, cellIndexesToFilter, loomFrom, loomTo))
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
			case FilterCols: 
				FileFilter.filterCells();
				break;
			case FilterRows: 
				FileFilter.filterGenes();
				break;
			case FilterDEMetadata:
				FileFilter.filterDEMetadata();
				break;
		}
	}
			
	public static String[] readMode(String[] args)
	{
		// Find -T
		int index = -1;
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.equals("-T")) {index = i; break;}
		}
		if(index == -1 || args.length < 2) new ErrorJSON("Argument -T is mandatory:\n-T %s \t\tMode to run ASAP. Please select in " + Mode.toArrayString());
			
		String[] args2 = new String[args.length - 2];
		int j = 0;
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			switch(arg)
			{
				case "-T":
					i++;
					try
					{
						m = Mode.valueOf(args[i]);
					}
					catch(IllegalArgumentException iae)
					{
						new ErrorJSON("This mode (-T) '" + args[i] + "' does not exist. Please select in " + Mode.toArrayString());
					}
					break;
				default:
					args2[j] = arg;
					j++;
			}
		}
		return args2;
	}
	
	public static String[] isDebug(String[] args)
	{
		// Find --debug
		String[] args2 = new String[args.length - 1]; // In case --debug exists
		int j = 0;
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.equals("--debug")) 
			{
				if(Parameters.debugMode) new ErrorJSON("There are more than once '--debug' option!");
				Parameters.debugMode = true;
			}
			else
			{
				if(j == args2.length) return args; // No debug
				args2[j] = args[i];
				j++;
			}
		}
		return args2;
	}
}
