package model;

import java.io.File;
import java.util.ArrayList;

import config.Config;
import db.DBManager;
import json.ErrorJSON;
import parsing.model.ColumnName;
import parsing.model.FileType;
import tools.Utils;

public class Parameters 
{
	public static int maxUniverse = 100000; // for FISHER EXACT TEST (max number of genes in background)
	public static enum OutputType{JSON, PLAIN_TEXT};
	
	// Debug
	public static boolean debugMode = false;
	
	// Chunking
	public static final int defaultChunkX = 64;
	public static final int defaultChunkY = 64;
	
	// Shared
	public static boolean isCountMatrix = true;
	public static boolean scientific = false;
	public static String organism_S = null;
	public static String outputFolder = null;
	public static String outputFile = null;
	public static OutputType outputType = OutputType.JSON;
	public static String group = null;
	public static String fileName = null;
	public static String selection = null;
	public static FileType fileType = null;
	public static String loomVersion = "3.0.0";
	public static String fitModel = null;
	public static String erccFile = null;
	public static int organism = 1;
	public static int taxon = -1;
	public static long nCells = -1;
	public static long nGenes = -1;
	public static long[] indexes = null;
	public static long[] stable_ids = null;
	public static boolean export_stable_ids = false;
	public static String loom_cell_stable_ids = null;
	public static String[] names = null;
	public static String metaName = null;
	public static MetaOn which = null;
	public static Metatype metatype = null;
	public static Metatype[] metatypes = null;
	public static boolean displayValues = true;
	public static boolean displayNames = false;
	public static boolean displayRowNames = false;
	public static boolean displayColNames = false;
	public static String loomFile = null;
	public static String loomFile2 = null;
	public static double idleTime = 0;
	public static boolean isIndex = false;
	public static boolean removeAmbiguous = false;
	public static String value = null;
	public static long index = -1;
	public static boolean sort = false;
	public static String defaultMissingValue = "-1";
	public static long randomSeed = 42;
	
	// Module score
	// Seurat
	public static int nBins = 24;
	public static int nBackgroundGenes = 100;
	
	// Dimension Reduction
	public static dim_reduction.model.Model dimReducModel = null;
	
	// t-SNE
	public static int perplexity = -1;
	
	// DE
	public static differential_expression.Model deModel = null;
	public static String iAnnot = null;
	public static String oAnnot = null;
	public static String gAnnot = null;
	public static String group_1 = null;
	public static String group_2 = null;
	public static long id = -1;
	
	// Normalization
	public static long scale_factor = 10000; 
	
	// Scaling
	public static long scale_max = 10; 
	public static boolean scale = true;
	public static boolean center = true;
	
	// Filtering
	public static filtering.model.Model filtModel = null;
	// CPM
	public static int nbCountsPerCell = -1;
	public static int nbCellsDetected = -1;
	
	// CreateDLFile
	public static String JSONFileName = null;
	
	// Enrichment
	public static String adjMethod = "fdr";
	public static long geneset_id = -1;
	public static long[] geneset_ids = null;
	public static enrichment.model.Model enrichModel = enrichment.model.Model.FET;
	public static int maxGenesInPathway = 500;
	public static int minGenesInPathway = 15;
	public static float pThreshold = 1;
	public static float fdrThreshold = 1;
	public static float fcThreshold = 0;
	public static int topThreshold = Integer.MAX_VALUE;
	
	// ModuleScore
	public static module_score.model.Model moduleScoreModel = module_score.model.Model.Seurat;
	public static String metaForComputation = null;
	
	// Parsing
	public static boolean has_header = false;
	public static ColumnName name_column = ColumnName.NONE;
	public static String delimiter = "\t";
	
	public static void load(String[] args, Mode m)
	{
		if(args.length == 0 && m != Mode.UpdateEnsemblDB)
		{
			printHelp(m);
			System.exit(0);
		}
		switch(m)
		{
			case CreateGODB:
				loadCreateGODB(args);
				break;
			case CreateKeggDB:
				loadCreateKeggDB(args);
				break;
			case UpdateEnsemblDB:
				loadUpdateEnsemblDB(args);
				break;
			case Enrichment:
				loadEnrichment(args);
				break;
			case MarkerEnrichment:
				loadMarkerEnrichment(args);
				break;
			case ModuleScore:
				loadModuleScore(args);
				break;
			case Preparsing: 
				loadPreparsing(args);
				break;
			case Parsing: 
				loadParsing(args);
				break;
			case PreparseMetadata: 
				loadPreparseMetadata(args);
				break;
			case ParseMetadata: 
				loadParseMetadata(args);
				break;
			case CreateCellSelection:
				loadCreateCellSelection(args);
				break;
			case FilterRows: 
				loadFilterGenes(args);
				break;
			case DimensionReduction: 
				loadDimensionReduction(args);
				break;
			case DifferentialExpression: 
				loadDifferentialExpression(args);
				break;
			case FindMarkers: 
				loadFindMarkers(args);
				break;
			case Normalization:
				loadNormalization(args);
				break;
			case Scaling:
				loadScaling(args);
				break;
			case RegenerateNewOrganism:
				loadRegenerateNewOrganism(args);
				break;
			case IndexByCell:
				loadIndexByCell(args);
				break;
			case GetIndex:
				loadGetIndex(args);
				break;
			case GetGeneStats:
				loadGetGeneStats(args);
				break;
			case ExtractRow:
				loadExtractRow(args);
				break;
			case ExtractCol:
				loadExtractCol(args);
				break;
			case ExtractDataset:
				loadExtractDataset(args);
				break;
			case ListMetadata:
				loadListMetadata(args);
				break;
			case ExtractMetadata:
				loadExtractMetadata(args);
				break;
			case MatchValues:
				loadMatchValues(args);
				break;
			case RemoveMetaData:
				loadRemoveMetadata(args);
				break;
			case CopyMetaData:
				loadCopyMetadata(args);
				break;
			case FilterCols:
				loadFilterCells(args);
				break;
			case FilterDEMetadata:
				loadFilterDEMetadata(args);
				break;
		}
	}
	
	public static void loadEnrichment(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-m":
						i++;
						switch(args[i])
						{
						case "fet":
							enrichModel = enrichment.model.Model.FET;
							break;
						default:
							new ErrorJSON("The entered model, "+args[i]+ ", does not exist!\nIt should be one of the following: [fet]");
						}
						break;						
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-f":
						i++;
						fileName = args[i];
						fileName = fileName.replaceAll("\\\\", "/");
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-geneset":
						i++;
						try
						{
							geneset_id = Long.parseLong(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-geneset' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-adj":
						i++;
						adjMethod = args[i];
						if(!adjMethod.equals("bonferroni") && !adjMethod.equals("fdr") && !adjMethod.equals("none"))
						{
							new ErrorJSON("The '-adj' option should be followed by one of those values: [bonferroni, fdr, none]. You entered " + args[i]);
						}
						break;
					case "-max":
						i++;
						try
						{
							maxGenesInPathway = Integer.parseInt(args[i]);
							if(maxGenesInPathway <0) throw new NumberFormatException();
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-max' option should be followed by an positive Integer value. You entered " + args[i]);
							
						}
						break;
					case "-min":
						i++;
						try
						{
							minGenesInPathway = Integer.parseInt(args[i]);
							if(minGenesInPathway <0) throw new NumberFormatException();
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-min' option should be followed by an positive Integer value. You entered " + args[i]);
						}
						break;
					case "-h":
						i++;
						DBManager.URL = Config.ConfigMAIN().getURLFromHost(args[i]);
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("No Loom file is specified, please choose one by using the '-loom' option.");
		if(geneset_id == -1) new ErrorJSON("No geneset id is specified, please choose one by using the '-geneset' option.");
		if(fileName == null) new ErrorJSON("No JSON file containing the list of genes to enrich, please choose one by using the '-f' option.");
	}
	
	public static void loadMarkerEnrichment(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{				
					case "-o":
					case "--output-folder":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) f.mkdirs();
						else if(!f.isDirectory()) new ErrorJSON(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "-i":
					case "--input-folder":
						i++;
						fileName = args[i];
						fileName = fileName.replaceAll("\\\\", "/");
						f = new File(fileName);
						if(!f.exists()) new ErrorJSON(fileName + " does not exist.");
						else if(!f.isDirectory()) new ErrorJSON(fileName + " is not a folder.");
						if(!fileName.endsWith("/")) fileName+= "/";
						break;
					case "--genesets":
						i++;
						try
						{
							String[] tokens = args[i].split(",");
							geneset_ids = new long[tokens.length];
							for(int k = 0; k < tokens.length; k++) 
							{
								geneset_ids[k] = Long.parseLong(tokens[k]);
								if(geneset_ids[k] < 0) new ErrorJSON("The '--genesets' option should be followed by positive Long value(s). You entered " + geneset_ids[k]);
							}
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '--genesets' option should be followed by an list of integers (eventually separated by commas if multiple). You entered " + args[i]);
						}
						break;
					case "-h":
						i++;
						DBManager.URL = Config.ConfigMAIN().getURLFromHost(args[i]);
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(geneset_ids == null) new ErrorJSON("No geneset id is specified, please choose one by using the '--genesets' option.");
		if(fileName == null) new ErrorJSON("No JSON file containing the list of genes to enrich, please choose one by using the '-i' option.");
		if(outputFolder == null) new ErrorJSON("No output folder is specified, please choose one by using the '-o' option.");
	}
	
	public static void loadModuleScore(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-m":
						i++;
						switch(args[i])
						{
						case "seurat":
							moduleScoreModel = module_score.model.Model.Seurat;
							break;
						case "pca":
							moduleScoreModel = module_score.model.Model.PCA;
							break;
						default:
							new ErrorJSON("The entered model, "+args[i]+ ", does not exist!\nIt should be one of the following: [seurat, pca]");
						}
						break;						
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-oAnnot":
						i++;
						oAnnot = args[i];
						if(!oAnnot.startsWith("/")) oAnnot = "/" + oAnnot;
						break;
					case "-seed":
						i++;
						try
						{
							randomSeed = Long.parseLong(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-seed' option should be followed by a Long. You entered " + args[i]);
						}
						break;
					case "-nBackgroundGenes":
						i++;
						try
						{
							nBackgroundGenes = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-nBackgroundGenes' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-nBins":
						i++;
						try
						{
							nBins = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-nBins' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-geneset":
						i++;
						try
						{
							geneset_id = Long.parseLong(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-geneset' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-dataset":
						i++;
						metaForComputation = args[i];
						break;
					case "-metadata":
						i++;
						metaName = args[i];
						break;
					case "-sel":
						i++;
						selection = args[i];
						break;
					case "-h":
						i++;
						DBManager.URL = Config.ConfigMAIN().getURLFromHost(args[i]);
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("No Loom file is specified, please choose one by using the '-loom' option.");
		if(geneset_id == -1 && metaName == null) new ErrorJSON("No geneset is specified, please choose one by using the '-geneset' or '-metadata' option.");
		if(geneset_id != -1 && metaName != null) new ErrorJSON("You specified both '-geneset' AND '-metadata' options, pick only ONE");
		if(selection != null && metaName == null) new ErrorJSON("You specified a selection with '-sel' BUT '-metadata' IS NOT SET");
		if(selection == null && metaName != null) new ErrorJSON("You specified a '-metadata' BUT selection with '-sel' IS NOT SET");
		if(metaForComputation == null) new ErrorJSON("You should specify a dataset to use for computing the score, using '-dataset'");
	}

	public static void loadUpdateEnsemblDB(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						new File(outputFolder).mkdirs();
						break;
					case "-h":
						i++;
						DBManager.URL = Config.ConfigMAIN().getURLFromHost(args[i]);
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null)
		{
			printHelp(Mode.UpdateEnsemblDB);
			new ErrorJSON("You need to specify an output folder with -o option.");
		}
	}
	
	public static void loadCreateKeggDB(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						new File(outputFolder).mkdirs();
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null)
		{
			printHelp(Mode.CreateKeggDB);
			new ErrorJSON("Please specify an output folder using -o option");
		}
	}
	
	public static void loadCreateGODB(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						new File(outputFolder).mkdirs();
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null)
		{
			printHelp(Mode.CreateGODB);
			new ErrorJSON("Please specify an output folder using -o option");
		}
	}
	
	public static void loadGetGeneStats(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "--loom":
						i++;
						loomFile = args[i];
						loomFile = loomFile.replaceAll("\\\\", "/");
						if(new File(loomFile).isDirectory()) new ErrorJSON("'-loom' should be followed by a Loom file name, not a folder name.");
						break;
					case "--iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "--id":
						i++;
						try
						{
							id = Long.parseLong(args[i]);
							if(id < 0) new ErrorJSON("The '--id' option should be followed by a Positive Long. You entered " + args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '--id' option should be followed by a Long. You entered " + args[i]);
						}
						break;
					case "--index":
						i++;
						try
						{
							index = Long.parseLong(args[i]);
							if(index < 0) new ErrorJSON("The '--index' option should be followed by a Positive Long. You entered " + args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '--index' option should be followed by a Long. You've entered " + args[i]);
						}
						break;
					case "--indexes":
						i++;
						if(args[i].equals("")) new ErrorJSON("The --indexes option should be followed by positive Long value(s).");
						try
						{
							String[] tokens = args[i].split(",");
							indexes = new long[tokens.length];
							for(int k = 0; k < tokens.length; k++) 
							{
								indexes[k] = Long.parseLong(tokens[k]);
								if(indexes[k] < 0) new ErrorJSON("The '--indexes' option should be followed by positive Long value(s). You entered " + indexes[k]);
							}
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '--indexes' option should be followed by positive Long value(s). You entered " + args[i]);
						}
						break;
					case "--stable-ids":
						i++;
						if(args[i].equals("")) new ErrorJSON("The --stable-ids option should be followed by positive Long value(s).");
						try
						{
							String[] tokens = args[i].split(",");
							stable_ids = new long[tokens.length];
							for(int k = 0; k < tokens.length; k++) 
							{
								stable_ids[k] = Long.parseLong(tokens[k]);
								if(stable_ids[k] < 0) new ErrorJSON("The '--stable-ids' option should be followed by positive Long value(s). You entered " + stable_ids[k]);
							}
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '--stable_ids' option should be followed by positive Long value(s). You entered " + args[i]);
						}
						break;
					case "--names":
						i++;
						if(args[i].equals("")) new ErrorJSON("The '--names' option should be followed by String(s).");
						names = args[i].split(",");
						for(int k = 0; k < names.length; k++) names[k] = names[k].trim().toUpperCase(); 
						break;
					case "--loom-cells":
						i++;
						if(args[i].equals("")) new ErrorJSON("The '--loom-cells' option should be followed by a Loom file name.");
						loom_cell_stable_ids = args[i];
						loom_cell_stable_ids = loom_cell_stable_ids.replaceAll("\\\\", "/");
						if(new File(loom_cell_stable_ids).isDirectory()) new ErrorJSON("'--loom-cells' should be followed by a Loom file name, not a folder name.");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("Loom file parameter is missing");
		if(indexes == null && names == null && stable_ids == null) new ErrorJSON("You need to use at least one of the following options --names, --indexes or --stable-ids");
		if(indexes != null && names != null) new ErrorJSON("Please select whether --indexes or --names options. Not both.");
		if(indexes != null && stable_ids != null) new ErrorJSON("Please select whether --indexes or --stable-ids options. Not both.");
		if(names != null && stable_ids != null) new ErrorJSON("Please select whether --stable-ids or --names options. Not both.");
		if(index != -1 || id != -1)
		{
			if(index == -1) new ErrorJSON("You used --id, please also set --index value.");
			if(id == -1) new ErrorJSON("You used --index, please also set --id value.");
		}
		if(index != -1 && loom_cell_stable_ids != null) new ErrorJSON("Please select whether --id/--index or --loom-cells options. Not both.");
	}
	
	public static void loadExtractRow(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						loomFile = args[i];
						loomFile = loomFile.replaceAll("\\\\", "/");
						if(new File(loomFile).isDirectory()) new ErrorJSON("'-loom' should be followed by a Loom file name, not a folder name.");
						break;
					case "-iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "-prec":
						i++;
						try
						{
							Utils.changeFormatter(Integer.parseInt(args[i]));
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-prec' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "--sort":
						sort = true;
						break;
					case "-indexes":
						i++;
						if(args[i].equals("")) new ErrorJSON("The -indexes option should be followed by String(s).");
						try
						{
							String[] tokens = args[i].split(",");
							indexes = new long[tokens.length];
							for(int k = 0; k < tokens.length; k++) 
							{
								indexes[k] = Long.parseLong(tokens[k]);
								if(indexes[k] < 0) new ErrorJSON("The '-indexes' option should be followed by positive Long value(s). You entered " + indexes[k]);
							}
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-indexes' option should be followed by positive Long value(s). You entered " + args[i]);
						}
						break;
					case "-stable_ids":
						i++;
						if(args[i].equals("")) new ErrorJSON("The -stable_ids option should be followed by String(s).");
						try
						{
							String[] tokens = args[i].split(",");
							stable_ids = new long[tokens.length];
							for(int k = 0; k < tokens.length; k++) 
							{
								stable_ids[k] = Long.parseLong(tokens[k]);
								if(stable_ids[k] < 0) new ErrorJSON("The '-stable_ids' option should be followed by positive Long value(s). You entered " + stable_ids[k]);
							}
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-stable_ids' option should be followed by positive Long value(s). You entered " + args[i]);
						}
						break;
					case "-loom_cells":
						i++;
						if(args[i].equals("")) new ErrorJSON("The -loom_cells option should be followed by String(s).");
						loom_cell_stable_ids = args[i];
						loom_cell_stable_ids = loom_cell_stable_ids.replaceAll("\\\\", "/");
						if(new File(loom_cell_stable_ids).isDirectory()) new ErrorJSON("'-loom_cell_stable_ids' should be followed by a Loom file name, not a folder name.");
						break;
					case "-names":
						i++;
						if(args[i].equals("")) new ErrorJSON("The '-names' option should be followed by String(s).");
						names = args[i].split(",");
						for(int k = 0; k < names.length; k++) names[k] = names[k].trim().toUpperCase(); 
						break;
					case "-display-names":
						displayNames = true;
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null)
		{
			printHelp(Mode.ExtractRow);
			new ErrorJSON("Loom file parameter is missing");
		}
		if(iAnnot == null) iAnnot = "/matrix";
		if(indexes == null && names == null && stable_ids == null) new ErrorJSON("You need to use at least one of the following options -names, -indexes or -stable_ids");
		if(indexes != null && names != null) new ErrorJSON("Please select whether -indexes or -names options. Not both.");
		if(indexes != null && stable_ids != null) new ErrorJSON("Please select whether indexes or -stable_ids options. Not both.");
		if(names != null && stable_ids != null) new ErrorJSON("Please select whether -stable_ids or -names options. Not both.");
	}
	
	public static void loadExtractCol(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						loomFile = args[i];
						loomFile = loomFile.replaceAll("\\\\", "/");
						break;
					case "-iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "-indexes":
						i++;
						if(args[i].equals("")) new ErrorJSON("The -indexes option should be followed by String(s).");
						try
						{
							String[] tokens = args[i].split(",");
							indexes = new long[tokens.length];
							for(int k = 0; k < tokens.length; k++) 
							{
								indexes[k] = Long.parseLong(tokens[k]);
								if(indexes[k] < 0) new ErrorJSON("The '-indexes' option should be followed by positive Long value(s). You entered " + indexes[k]);
							}
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-indexes' option should be followed by positive Long value(s). You entered " + args[i]);
						}
						break;
					case "-prec":
						i++;
						try
						{
							Utils.changeFormatter(Integer.parseInt(args[i]));
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-prec' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-names":
						i++;
						if(args[i].equals("")) new ErrorJSON("The '-names' option should be followed by String(s).");
						names = args[i].split(",");
						for(int k = 0; k < names.length; k++) names[k] = names[k].trim().toUpperCase(); 
						break;
					case "-display-names":
						displayNames = true;
						break;
					case "-display-stable-ids":
						export_stable_ids = true;
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null)
		{
			printHelp(Mode.ExtractCol);
			System.exit(-1);
		}
		if(iAnnot == null) iAnnot = "/matrix";
		if(indexes == null && names == null) new ErrorJSON("You need to use at least one of the two options -names or -indexes");
		if(indexes != null && names != null) new ErrorJSON("You need to use only one of the two options -names or -indexes, not both");
	}
	
	public static void loadExtractDataset(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
					case "--output-file":
						i++;
						outputFile = args[i];
						outputFile = outputFile.replaceAll("\\\\", "/");
						if(new File(outputFile).isDirectory()) new ErrorJSON("'-o' | '--output-file' should be followed by a file name, not a folder name.");
						break;
					case "-t":
					case "--output-type":
						i++;
						try
						{
							outputType = OutputType.valueOf(args[i]);
						}
						catch(IllegalArgumentException iae)
						{
							new ErrorJSON("This output type '" + args[i] + "' does not exist. Please select in [JSON, PLAIN_TEXT]");
						}
						break;
					case "--loom":
						i++;
						loomFile = args[i];
						loomFile = loomFile.replaceAll("\\\\", "/");
						if(new File(loomFile).isDirectory()) new ErrorJSON("'--loom' should be followed by a Loom file name, not a folder name.");
						break;
					case "--iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "--prec":
						i++;
						try
						{
							Utils.changeFormatter(Integer.parseInt(args[i]));
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-prec' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "--row-names":
						displayRowNames = true;
						break;
					case "--col-names":
						displayColNames = true;
						break;
					case "--col-indexes":
						i++;
						if(args[i].equals("")) new ErrorJSON("The --col-indexes option should be followed by String(s).");
						try
						{
							String[] tokens = args[i].split(",");
							indexes = new long[tokens.length];
							for(int k = 0; k < tokens.length; k++) 
							{
								indexes[k] = Long.parseLong(tokens[k]);
								if(indexes[k] < 0) new ErrorJSON("The '--col-indexes' option should be followed by positive Long value(s). You entered " + indexes[k]);
							}
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '--col-indexes' option should be followed by positive Long value(s). You entered " + args[i]);
						}
						break;
					case "--col-stableids":
						i++;
						if(args[i].equals("")) new ErrorJSON("The --col-stableids option should be followed by String(s).");
						try
						{
							String[] tokens = args[i].split(",");
							stable_ids = new long[tokens.length];
							for(int k = 0; k < tokens.length; k++) 
							{
								stable_ids[k] = Long.parseLong(tokens[k]);
								if(stable_ids[k] < 0) new ErrorJSON("The '--col-stableids' option should be followed by positive Long value(s). You entered " + stable_ids[k]);
							}
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '--col-stableids' option should be followed by positive Long value(s). You entered " + args[i]);
						}
						break;
					case "--loom-cells":
						i++;
						if(args[i].equals("")) new ErrorJSON("The --loom-cells option should be followed by String(s).");
						loom_cell_stable_ids = args[i];
						loom_cell_stable_ids = loom_cell_stable_ids.replaceAll("\\\\", "/");
						if(new File(loom_cell_stable_ids).isDirectory()) new ErrorJSON("'--loom-cells' should be followed by a Loom file name, not a folder name.");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null)
		{
			printHelp(Mode.ExtractDataset);
			new ErrorJSON("Loom file parameter is missing", outputFile);
		}
		if(iAnnot == null) iAnnot = "/matrix";
		if(loom_cell_stable_ids == null) loom_cell_stable_ids = loomFile;
	}
	
	public static void loadListMetadata(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-f":
						i++;
						fileName = args[i];
						fileName = fileName.replaceAll("\\\\", "/");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(fileName == null)
		{
			printHelp(Mode.ListMetadata);
			System.exit(-1);
		}
	}
	
	public static void loadExtractMetadata(String[] args)
	{
		int precision = 3; // default
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFile = args[i];
						outputFile = outputFile.replaceAll("\\\\", "/");
						if(new File(outputFile).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						loomFile = args[i];
						loomFile = loomFile.replaceAll("\\\\", "/");
						break;
					case "-meta":
						i++;
						metaName = args[i];
						break;
					case "-metaJSON":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-metaJSON' should be followed by a JSON file containing metadata path in the Loom. Here it is a folder name!");
						if(!new File(JSONFileName).exists()) new ErrorJSON("'-metaJSON' should be followed by a JSON file containing metadata path in the Loom. Here, the file does not exist!");
						break;
					case "-type":
						i++;
						try
						{
							metatype = Metatype.valueOf(args[i]);
						}
						catch(IllegalArgumentException iae)
						{
							new ErrorJSON("This metadata type '" + args[i] + "' does not exist.");
						}
						break;
					case "-prec":
						i++;
						try
						{
							precision = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-prec' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-no-values":
						displayValues = false;
						break;
					case "--scientific":
						scientific = true;
						break;
					case "-names":
						displayNames = true;
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		Utils.changeFormatter(precision);
		if(loomFile == null) new ErrorJSON("Please specify a loom file to extract metadata from using -loom option");
		if(metaName == null && JSONFileName == null) new ErrorJSON("Please specify metadata(s) using -meta or -metaJSON option");
		if(metaName != null && JSONFileName != null) new ErrorJSON("Please specify metadata(s) using -meta OR -metaJSON option. You specified both.");
	}
	
	public static void loadMatchValues(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						loomFile = args[i];
						loomFile = loomFile.replaceAll("\\\\", "/");
						break;
					case "-iAnnot":
						i++;
						iAnnot = args[i];
						break;
					case "-value":
						i++;
						value = args[i];
						break;
					case "-light":
						displayValues = false;
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("You need to specify a Loom file using the -loom option.");
		if(iAnnot == null) new ErrorJSON("You need to specify a metadata path to extract using -iAnnot option.");
		if(value == null) new ErrorJSON("You need to specify a value for the metadata to select indexes equal to this value using -value option.");
	}
	
	public static void loadCreateCellSelection(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFile = args[i];
						outputFile = outputFile.replaceAll("\\\\", "/");
						if(new File(outputFile).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-f":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-f' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-meta":
						i++;
						metaName = args[i];
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(JSONFileName == null)
		{
			printHelp(Mode.CreateCellSelection);
			new ErrorJSON("Please specify the input JSON file (containing list of indexes of cells to filter) using -f option");
		}
		if(loomFile == null)
		{
			printHelp(Mode.CreateCellSelection);
			new ErrorJSON("Please specify a loom file using -loom option");
		}
		if(metaName == null)
		{
			printHelp(Mode.CreateCellSelection);
			new ErrorJSON("Please specify a Metadata Path to create using -meta option");
		}
	}
	
	public static void loadRemoveMetadata(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFile = args[i];
						outputFile = outputFile.replaceAll("\\\\", "/");
						if(new File(outputFile).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-metaJSON":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-metaJSON' should be followed by a JSON file containing metadata path in the Loom. Here it is a folder name!");
						if(!new File(JSONFileName).exists()) new ErrorJSON("'-metaJSON' should be followed by a JSON file containing metadata path in the Loom. Here, the file does not exist!");
						break;
					case "-meta":
						i++;
						metaName = args[i];
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null)
		{
			printHelp(Mode.RemoveMetaData);
			new ErrorJSON("Please specify a loom file using -loom option");
		}
		if(metaName == null && JSONFileName == null)
		{
			printHelp(Mode.RemoveMetaData);
			new ErrorJSON("Please specify a metadata using -meta or -metaJSON option");
		}
	}
	
	public static void loadCopyMetadata(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-loomFrom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-loomTo":
						i++;
						try
						{
							loomFile2 = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile2);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile2);
							if(!c.isFile()) new ErrorJSON(loomFile2 + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-metaJSON":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-metaJSON' should be followed by a JSON file containing metadata path in the Loom. Here it is a folder name!");
						if(!new File(JSONFileName).exists()) new ErrorJSON("'-metaJSON' should be followed by a JSON file containing metadata path in the Loom. Here, the file does not exist!");
						break;
					case "-meta":
						i++;
						metaName = args[i];
						break;
					case "-o":
						i++;
						outputFile = args[i];
						outputFile = outputFile.replaceAll("\\\\", "/");
						if(new File(outputFile).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("Please specify a loom file to copy FROM using -loomFrom option");
		if(loomFile2 == null) new ErrorJSON("Please specify a loom file to copy TO using -loomTo option");
		if(metaName == null && JSONFileName == null) new ErrorJSON("Please specify metadata(s) using -meta or -metaJSON option");
		if(metaName != null && JSONFileName != null) new ErrorJSON("Please specify metadata(s) using -meta OR -metaJSON option. You specified both.");
	}
	
	public static void loadFilterCells(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) f.mkdirs();
						else if(!f.isDirectory()) new ErrorJSON(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-col_indexes_file":
						isIndex = true;
					case "-col_names_file":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-col_names_file' should be followed by a JSON file name containing cell names, not a folder name.");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null)
		{
			printHelp(Mode.FilterCols);
			new ErrorJSON("Please specify an output folder using -o option");
		}
		if(loomFile == null)
		{
			printHelp(Mode.FilterCols);
			new ErrorJSON("Please specify a loom file to filter using -loom option");
		}
		if(JSONFileName == null)
		{
			printHelp(Mode.FilterCols);
			new ErrorJSON("Please specify a JSON file containing Cells to filter out using -col_names_file option");
		}
	}
	
	public static void loadFilterDEMetadata(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "-light":
						displayValues = false;
						break;
					case "-p":
						i++;
						try
						{
							pThreshold = Float.parseFloat(args[i].replaceAll(",", "."));
							if(pThreshold < 0 || pThreshold > 1) new ErrorJSON("The '-p' option should be followed by a float in [0,1]. You entered " + args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-p' option should be followed by a Float. You entered " + args[i]);
						}
						break;
					case "-fdr":
						i++;
						try
						{
							fdrThreshold = Float.parseFloat(args[i].replaceAll(",", "."));
							if(fdrThreshold < 0 || fdrThreshold > 1) new ErrorJSON("The '-fdr' option should be followed by a float in [0,1]. You entered " + args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-fdr' option should be followed by a Float. You entered " + args[i]);
						}
						break;
					case "-fc":
						i++;
						try
						{
							fcThreshold = Float.parseFloat(args[i].replaceAll(",", "."));
							if(fcThreshold < 0 ) new ErrorJSON("The '-fc' option should be followed by a POSITIVE Float. You entered " + args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-fc' option should be followed by a Float. You entered " + args[i]);
						}
						break;
					case "-top":
						i++;
						try
						{
							topThreshold = Integer.parseInt(args[i]);
							if(topThreshold <= 0 ) new ErrorJSON("The '-top' option should be followed by a POSITIVE Float. You entered " + args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-top' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("No Loom file is specified, please use the '-loom' option.\n");
		if(iAnnot == null) new ErrorJSON("No metadata (for DE results) is specified, please use the '-iAnnot' option.\n");
	}
		
	public static void loadRegenerateNewOrganism(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						new File(outputFolder).mkdirs();
						break;
					case "-organism":
						i++;
						try
						{
							organism = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-organism' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-j":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(JSONFileName == null || outputFolder == null)
		{
			printHelp(Mode.RegenerateNewOrganism);
			System.exit(-1);
		}
	}
	
	public static void loadIndexByCell(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
					case "--output":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) f.mkdirs();
						else if(!f.isDirectory()) new ErrorJSON(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "--loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--loomIndex":
						i++;
						try
						{
							loomFile2 = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile2);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile2);
							if(!c.isFile()) new ErrorJSON(loomFile2 + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '--loomIndex' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "--id":
						i++;
						try
						{
							id = Long.parseLong(args[i]);
							if(id < 0) new ErrorJSON("The '--id' option should be followed by a Positive Long. You entered " + args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '--id' option should be followed by a Long. You entered " + args[i]);
						}
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("No Loom file is specified, please use the '--loom' option.\n");
		if(loomFile2 == null) new ErrorJSON("No Index Loom file is specified, please use the '--loomIndex' option.\n");
		if(iAnnot == null) new ErrorJSON("No input dataset is specified, please use the '--iAnnot' option.\n");
		if(id == -1) new ErrorJSON("No table id is specified for the metadata, please use the '--id' option.\n");
		if(outputFolder != null) new File(outputFolder).mkdirs();
	}
	
	public static void loadGetIndex(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
					case "--output":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) f.mkdirs();
						else if(!f.isDirectory()) new ErrorJSON(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "--loomIndex":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '--loomIndex' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--jsonListCells":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'--jsonListCells' should be followed by a JSON file name containing cell names, not a folder name.");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("No Index Loom file is specified, please use the '--loomIndex' option.\n");
		if(JSONFileName == null) new ErrorJSON("No JSON file name (with cell names) is specified, please use the '--jsonListCells' option.\n");
		if(outputFolder != null) new File(outputFolder).mkdirs();
	}
	
	public static void loadParsing(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) f.mkdirs();
						else if(!f.isDirectory()) new ErrorJSON(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "-type":
						i++;
						try
						{
							fileType = FileType.valueOf(args[i]);
						}
						catch(IllegalArgumentException iae)
						{
							new ErrorJSON("This file type '" + args[i] + "' does not exist.");
						}
						break;
					case "-ncells":
						i++;
						try
						{
							nCells = Long.parseLong(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-ncells' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-ngenes":
						i++;
						try
						{
							nGenes = Long.parseLong(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-ngenes' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-f":
						i++;
						try
						{
							fileName = args[i].replaceAll("\\\\", "/");
							File c = new File(fileName);
							if(!c.exists()) new ErrorJSON("No file at path " + fileName);
							if(!c.isFile()) new ErrorJSON(fileName + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-f' option should be followed by dataset file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-header":
						i++;
						has_header = Boolean.parseBoolean(args[i]);
						break;
					case "-d":
						i++;
						delimiter = args[i];
						break;
					case "-organism":
						i++;
						try
						{
							organism = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-organism' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-col":
						i++;
						switch(args[i])
						{
							case "first": name_column = ColumnName.FIRST; break;
							case "last": name_column = ColumnName.LAST; break;
							case "none": name_column = ColumnName.NONE; break;
							default:
								new ErrorJSON("The '-col' option should be followed by [last, first, none]. You entered " + args[i]);
						}
						break;
					case "-sel":
						i++;
						selection = args[i];
						selection = selection.replaceAll("\\\\", "/");
						break;
					case "-h":
						i++;
						DBManager.URL = Config.ConfigMAIN().getURLFromHost(args[i]);
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null) new ErrorJSON("-o is mandatory option for parsing.");
		if(fileName == null) new ErrorJSON("-f is mandatory option for parsing.");
		if(nCells == -1) new ErrorJSON("-ncells is mandatory option for parsing.");
		if(nGenes == -1) new ErrorJSON("-ngenes is mandatory option for parsing.");
		if(fileType == null) new ErrorJSON("-type is mandatory option for parsing.");
		if(fileType == FileType.H5_10x || fileType == FileType.LOOM)
		{
			if(has_header) new ErrorJSON("The header should not be set to true if filetype = H5 or LOOM");
			if(!delimiter.equals("\t")) new ErrorJSON("The Delimiter should not be set when filetype = H5 or LOOM");
			if(name_column != ColumnName.NONE) new ErrorJSON("There should not be a '-col' option when filetype = H5 or LOOM");
		}
	}
	
	public static void loadParseMetadata(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-type":
						i++;
						try
						{
							fileType = FileType.valueOf(args[i]);
						}
						catch(IllegalArgumentException iae)
						{
							new ErrorJSON("This file type '" + args[i] + "' does not exist.");
						}
						break;
					case "-f":
						i++;
						try
						{
							fileName = args[i].replaceAll("\\\\", "/");
							File c = new File(fileName);
							if(!c.exists()) new ErrorJSON("No file at path " + fileName);
							if(!c.isFile()) new ErrorJSON(fileName + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-f' option should be followed by dataset file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-header":
						i++;
						has_header = Boolean.parseBoolean(args[i]);
						break;
					case "-d":
						i++;
						delimiter = args[i];
						break;
					case "-col":
						i++;
						switch(args[i])
						{
							case "first": name_column = ColumnName.FIRST; break;
							case "last": name_column = ColumnName.LAST; break;
							case "none": name_column = ColumnName.NONE; break;
							default:
								new ErrorJSON("The '-col' option should be followed by [last, first, none]. You entered " + args[i]);
						}
						break;
					case "-sel":
						i++;
						selection = args[i];
						selection = selection.replaceAll("\\\\", "/");
						break;
					case "-which":
						i++;
						switch(args[i].toLowerCase())
						{
							case "gene":
								Parameters.which = MetaOn.GENE;
								break;
							case "cell":
								Parameters.which = MetaOn.CELL;
								break;
							default:
								new ErrorJSON(args[i] + " is incompatible with the 'which' option. Please choose between [cell, gene].");
						}
						break;
					case "-metadataType":
						i++;
						String[] mTypes = args[i].split(",");
						try
						{
							ArrayList<Metatype> mT = new ArrayList<Metatype>();
							for(String m:mTypes)
							{
								m = m.trim();
								if(!m.equals("")) mT.add(Metatype.valueOf(m));
							}
							metatypes = new Metatype[mT.size()];
							for(int k = 0; k < metatypes.length; k++) metatypes[k] = mT.get(k);
						}
						catch(IllegalArgumentException iae)
						{
							new ErrorJSON("Unknown metadata type in '" + args[i] + "'");
						}
						break;
					case "-removeAmbiguous":
						removeAmbiguous = true;
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		//if(removeAmbiguous) System.out.println("You chose to remove ambiguous genes.");
		//if(which != null) System.out.println("You specified that the metadata was applied on " + which);
		//else 
		//{
		//	System.out.println("You DID NOT specify if the metadata was applied on [cell,gene] using the -which option. Set default to CELL.");
		//	which = MetaOn.CELL;
		//}
		if(fileName == null) new ErrorJSON("Please specify a metadata file to preparse with '-f'");
		if(loomFile == null) new ErrorJSON("Please specify a Loom file to which adding metadata with '-loom'");
		if(fileType == null) new ErrorJSON("-type is mandatory option for metadata parsing.");
		if(fileType == FileType.H5_10x || fileType == FileType.LOOM) new ErrorJSON("Cannot extract metadata from these files atm.");
	}
	
	public static void loadPreparsing(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						new File(outputFolder).mkdirs();
						break;
					case "-f":
						i++;
						fileName = args[i];
						fileName = fileName.replaceAll("\\\\", "/");
						break;
					case "-header":
						i++;
						has_header = Boolean.parseBoolean(args[i]);
						break;
					case "-d":
						i++;
						delimiter = args[i];
						break;
					case "-col":
						i++;
						switch(args[i])
						{
							case "first": name_column = ColumnName.FIRST; break;
							case "last": name_column = ColumnName.LAST; break;
							case "none": name_column = ColumnName.NONE; break;
							default:
								new ErrorJSON("The '-col' option should be followed by [last, first, none]. You entered " + args[i]);
						}
						break;
					case "-organism":
						i++;
						try
						{
							organism = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-organism' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-sel":
						i++;
						selection = args[i];
						selection = selection.replaceAll("\\\\", "/");
						break;
					case "-h":
						i++;
						DBManager.URL = Config.ConfigMAIN().getURLFromHost(args[i]);
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null) new ErrorJSON("The output folder '-o' is mandatory for Preparsing.");
		if(fileName == null) new ErrorJSON("Please specify a file to preparse with 'f'");
	}
	
	public static void loadPreparseMetadata(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-f":
						i++;
						fileName = args[i];
						fileName = fileName.replaceAll("\\\\", "/");
						break;
					case "-header":
						i++;
						has_header = Boolean.parseBoolean(args[i]);
						break;
					case "-d":
						i++;
						delimiter = args[i];
						break;
					case "-col":
						i++;
						switch(args[i])
						{
							case "first": name_column = ColumnName.FIRST; break;
							case "last": name_column = ColumnName.LAST; break;
							case "none": name_column = ColumnName.NONE; break;
							default:
								new ErrorJSON("The '-col' option should be followed by [last, first, none]. You entered " + args[i]);
						}
						break;
					case "-sel":
						i++;
						selection = args[i];
						selection = selection.replaceAll("\\\\", "/");
						break;
					case "-which":
						i++;
						switch(args[i].toLowerCase())
						{
							case "gene":
								Parameters.which = MetaOn.GENE;
								break;
							case "cell":
								Parameters.which = MetaOn.CELL;
								break;
							default:
								new ErrorJSON(args[i] + " is incompatible with the 'which' option. Please choose between [cell, gene].");
						}
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		//if(which != null) System.out.println("You specified that the metadata was applied on " + which);
		//else 
		//{
		//	System.out.println("You DID NOT specify if the metadata was applied on [cell,gene] using the -which option. Set default to CELL.");
		//	which = MetaOn.CELL;
		//}
		if(fileName == null) new ErrorJSON("Please specify a metadata file to preparse with '-f'");
		if(loomFile == null) new ErrorJSON("Please specify a Loom file to which adding metadata with '-loom'");
	}
	
	public static void loadFilterGenes(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						new File(outputFolder).mkdirs();
						break;
					case "-loom":
						i++;
						loomFile = args[i];
						loomFile = loomFile.replaceAll("\\\\", "/");
						break;
					case "-row_names_file":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-row_names_file' should be followed by a JSON file name containing gene IDs to keep, not a folder name.");
						break;
					case "-m":
						i++;
						switch(args[i])
						{
						case "basic":
							filtModel = filtering.model.Model.BASIC;
							break;
						case "cpm":
							filtModel = filtering.model.Model.CPM;
							i++;
							try
							{
								nbCountsPerCell = Integer.parseInt(args[i]);
							}
							catch(NumberFormatException nfe)
							{
								new ErrorJSON("The '-m cpm' model should be followed by two Integers: 'nbCountsPerCell nbCellsDetected'. The first one you entered is not an Integer: " + args[i]);
							}
							catch(ArrayIndexOutOfBoundsException aioobe)
							{
								new ErrorJSON("The '-m cpm' model should be followed by two Integers: 'nbCountsPerCell nbCellsDetected'. The first one is missing!");
							}
							i++;
							try
							{
								nbCellsDetected = Integer.parseInt(args[i]);
							}
							catch(NumberFormatException nfe)
							{
								new ErrorJSON("The '-m cpm' model should be followed by two Integers: 'nbCountsPerCell nbCellsDetected'. The second one you entered is not an Integer: " + args[i]);
							}
							catch(ArrayIndexOutOfBoundsException aioobe)
							{
								new ErrorJSON("The '-m cpm' model should be followed by two Integers: 'nbCountsPerCell nbCellsDetected'. The second one is missing!");
							}
							if(nbCountsPerCell < 0 || nbCellsDetected < 0) new ErrorJSON("The '-m cpm' model should be followed by two positive Integers. You entered '" + nbCountsPerCell + " " + nbCellsDetected + "'." );
							break;
						case "keep":
							filtModel = filtering.model.Model.KEEP;
							break;
						default:
							new ErrorJSON("The entered model, "+args[i]+ ", does not exist!\nIt should be one of the following: [basic, cpm, keep]");
						}
						break;				
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(filtModel== null || outputFolder == null || loomFile == null)
		{
			printHelp(Mode.FilterRows);
			String error = "Filtering cannot be run because parameters are missing:\n";
			if(filtModel == null) error += "No model is specified, please choose a model by using the '-m' option.\n";
			if(loomFile == null) error += "No file is specified, please choose a data file by using the '-f' option.\n";
			if(outputFolder == null) error += "No output folder is specified, please choose an output file by using the '-o' option.\n";
			new ErrorJSON(error); // This stops the program
		}
		if(filtModel == filtering.model.Model.KEEP && JSONFileName == null)
		{
			printHelp(Mode.FilterRows);
			new ErrorJSON("Please specify a JSON file containing Cells to filter out after filtering using -row_names_file option");
		}
		if(filtModel == filtering.model.Model.CPM && (nbCountsPerCell == -1 || nbCellsDetected == -1))
		{
			new ErrorJSON("The '-m cpm' model should be followed by two Integers: 'nbCountsPerCell nbCellsDetected'.");
		}
	}
	
	public static void loadDimensionReduction(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						new File(outputFolder).mkdirs();
						break;
					case "-f":
						i++;
						fileName = args[i];
						fileName = fileName.replaceAll("\\\\", "/");
						break;
					case "-m":
						i++;
						switch(args[i])
						{
						case "pca":
							dimReducModel = dim_reduction.model.Model.PCA;
							break;
						case "mds":
							dimReducModel = dim_reduction.model.Model.MDS;
							break;
						case "zifa":
							dimReducModel = dim_reduction.model.Model.ZIFA;
							break;
						case "tsne":
							dimReducModel = dim_reduction.model.Model.TSNE;
							i++;
							try
							{
								perplexity = Integer.parseInt(args[i]);
							}
							catch(NumberFormatException nfe)
							{
								new ErrorJSON("The '-m tsne' model should be followed by an Integer: 'perplexity'. The value you entered is not an Integer: " + args[i]);
							}
							catch(ArrayIndexOutOfBoundsException aioobe)
							{
								new ErrorJSON("The '-m tsne' model should be followed by an Integer: 'perplexity'. This parameter is missing!");
							}
							if(perplexity <= 0) new ErrorJSON("The '-m tsne' model should be followed by an Integer: 'perplexity' greater than 0. You entered '" + perplexity + "'." );
							break;
						default:
							new ErrorJSON("The entered model, "+args[i]+ ", does not exist!\nIt should be one of the following: [pca, mds, tsne, zifa]");
						}
						break;				
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(dimReducModel == null || outputFolder == null || fileName == null)
		{
			printHelp(Mode.DimensionReduction);
			String error = "Filtering cannot be run because parameters are missing:\n";
			if(dimReducModel == null) error += "No model is specified, please choose a model by using the '-m' option.\n";
			if(fileName == null) error += "No file is specified, please choose a data file by using the '-f' option.\n";
			if(outputFolder == null) error += "No output folder is specified, please choose an output file by using the '-o' option.\n";
			new ErrorJSON(error); // This stops the program
		}
		new File(outputFolder).mkdirs();
		if(dimReducModel == dim_reduction.model.Model.TSNE && perplexity == -1)
		{
			new ErrorJSON("The '-m tsne' model should be followed by an Integer: 'perplexity'. This parameter is missing!");
		}
	}
	
	public static void loadDifferentialExpression(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) f.mkdirs();
						else if(!f.isDirectory()) new ErrorJSON(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "-oAnnot":
						i++;
						oAnnot = args[i];
						if(!oAnnot.startsWith("/")) oAnnot = "/" + oAnnot;
						break;
					case "-gAnnot":
						i++;
						gAnnot = args[i];
						if(!gAnnot.startsWith("/")) gAnnot = "/" + gAnnot;
						break;
					case "-g1":
						i++;
						group_1 = args[i];
						break;
					case "-g2":
						i++;
						if(args[i].equals("null")) group_2 = null;
						else group_2 = args[i];
						break;
					case "-m":
						i++;
						switch(args[i])
						{
						case "wilcox_asap":
							deModel = differential_expression.Model.Wilcoxon;
							break;
						default:
							new ErrorJSON("The entered model, "+args[i]+ ", does not exist!\nIt should be one of the following: [wilcox-asap]");
						}
						break;
					case "--is_count_table":
						i++;
						if(args[i].equals("false")) isCountMatrix = false;
						else if(args[i].equals("true")) isCountMatrix = true;
						else new ErrorJSON("The parameter '--is_count_table' should be followed by one of the following: [true, false]");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(deModel == null) new ErrorJSON("No model is specified, please choose a model by using the '-m' option.\n");
		if(loomFile == null) new ErrorJSON("No Loom file is specified, please use the '-loom' option.\n");
		if(outputFolder == null) new ErrorJSON("No Output folder is specified, please use the '-o' option.\n");
		if(iAnnot == null) new ErrorJSON("No input dataset is specified, please use the '-iAnnot' option.\n");
		if(oAnnot == null) new ErrorJSON("No output metadata for results is specified, please use the '-oAnnot' option.\n");
		if(gAnnot == null) new ErrorJSON("No group metadata is specified, please use the '-gAnnot' option.\n");
		if(group_1 == null) new ErrorJSON("No reference group is specified, please use the '-g1' option.\n");
		if(group_2 == null) System.out.println("No comparison group is specified, will perform reference group vs the rest of the groups");
		if(group_1 == group_2) new ErrorJSON("Cannot compute DE from the same group.");
		new File(outputFolder).mkdirs();
	}
	
	public static void loadFindMarkers(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
					case "--output":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) f.mkdirs();
						else if(!f.isDirectory()) new ErrorJSON(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "--loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "--id":
						i++;
						try
						{
							id = Long.parseLong(args[i]);
							if(id < 0) new ErrorJSON("The '--id' option should be followed by a Positive Long. You entered " + args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '--id' option should be followed by a Long. You entered " + args[i]);
						}
						break;
					case "--is_count_table":
						i++;
						if(args[i].equals("false")) isCountMatrix = false;
						else if(args[i].equals("true")) isCountMatrix = true;
						else new ErrorJSON("The parameter '--is_count_table' should be followed by one of the following: [true, false]");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("No Loom file is specified, please use the '--loom' option");
		if(iAnnot == null) new ErrorJSON("No input dataset is specified, please use the '--iAnnot' option");
		if(id == -1) new ErrorJSON("No table id is specified for the metadata, please use the '--id' option");
		if(outputFolder == null) new ErrorJSON("No output folder was specified with '-o' option. Please specify an output folder.");
		new File(outputFolder).mkdirs();
	}
	
	public static void loadNormalization(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-oAnnot":
						i++;
						oAnnot = args[i];
						if(!oAnnot.startsWith("/")) oAnnot = "/" + oAnnot;
						break;
					case "-scaleFactor":
						i++;
						try
						{
							scale_factor = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-scaleFactor' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					default:
						new ErrorJSON("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("No Loom file is specified, please use the '-loom' option.\n");
		if(oAnnot == null) new ErrorJSON("No output metadata for results is specified, please use the '-oAnnot' option.\n");
	}
	
	public static void loadScaling(String[] args)
	{
		for(int i = 0; i < args.length; i++) 
		{
			String arg = args[i];
			if(arg.startsWith("-"))
			{
				switch(arg)
				{
					case "-o":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						if(new File(JSONFileName).isDirectory()) new ErrorJSON("'-o' should be followed by a JSON file name, not a folder name.");
						break;
					case "-loom":
						i++;
						try
						{
							loomFile = args[i].replaceAll("\\\\", "/");
							File c = new File(loomFile);
							if(!c.exists()) new ErrorJSON("No file at path " + loomFile);
							if(!c.isFile()) new ErrorJSON(loomFile + " is not a file");
						}
						catch(Exception e)
						{
							new ErrorJSON("The '-loom' option should be followed by Loom file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "-oAnnot":
						i++;
						oAnnot = args[i];
						if(!oAnnot.startsWith("/")) oAnnot = "/" + oAnnot;
						break;
					case "-scaleMax":
						i++;
						try
						{
							scale_max = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-scaleMax' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-scale":
						i++;
						scale = Boolean.parseBoolean(args[i]);
						break;
					case "-center":
						i++;
						center = Boolean.parseBoolean(args[i]);
						break;
					default:
						new ErrorJSON("Unused argument: " + arg);
				}
			}
		}
		if(loomFile == null) new ErrorJSON("No Loom file is specified, please use the '-loom' option.\n");
		if(oAnnot == null) new ErrorJSON("No output metadata for results is specified, please use the '-oAnnot' option.\n");
		if(iAnnot == null) new ErrorJSON("No input metadata for results is specified, please use the '-iAnnot' option.\n");
	}
	
	public static void printHelp(Mode m)
	{
		switch(m)
		{
			case Enrichment: 
				System.out.println("Enrichment Mode\n\nOptions:");
				System.out.println("-m %s \t\tChoose a model amongst [fet, default = fet].");
				System.out.println("-o %s \t\t[Optional] Output JSON file.");
				System.out.println("-loom %s \t[Required] Loom file containing the genes to enrich.");
				System.out.println("-f %s \t\t[Required] The JSON file containing the selected cell indexes.");
				System.out.println("-geneset %i \t[Required]Id of the gene_sets to use for enrichment in the DB.");
				System.out.println("-adj %s \tStatitical adjustment method for multiple comparision [bonferroni, fdr, or none, default = fdr].");
				System.out.println("-min %i \tMinimum number of genes in a pathway for being taken into consideration [default = 15].");
				System.out.println("-max %i \tMaximum number of genes in a pathway for being taken into consideration [default = 500].");
				if(Parameters.debugMode) System.out.println("-h %s \t\tTo change specific host (default is " + Config.ConfigDEV().getProperty("mDbHost") + ")");
				else System.out.println("-h %s \t\tTo change specific host (default is " + Config.ConfigMAIN().getProperty("mDbHost") + ")");
				break;
			case MarkerEnrichment: 
				System.out.println("MarkerEnrichment Mode\n\nOptions:");
				System.out.println("-o | --output-folder %s \t[Required] Output folder");
				System.out.println("-i | --input-folder %s \t\t[Required] The folder path containing the FindMarker's output");
				System.out.println("--genesets %i,%i,... \t\t[Required]Id of the gene_sets to use for enrichment in the DB (separated by comma, if multiple).");
				if(Parameters.debugMode) System.out.println("-h %s \t\t\t\tTo change specific host (default is " + Config.ConfigDEV().getProperty("mDbHost") + ")");
				else System.out.println("-h %s \t\t\t\tTo change specific host (default is " + Config.ConfigMAIN().getProperty("mDbHost") + ")");
				break;
			case ModuleScore: 
				System.out.println("ModuleScore Mode\n\nOptions:");
				System.out.println("-m %s \t\tChoose a model amongst [seurat, pca] (default = seurat)");
				System.out.println("-o %s \t\t[Optional] Output JSON file.");
				System.out.println("-loom %s \t[Required] Loom file to compute the module score on.");
				System.out.println("-geneset %i \t[Required, or -metadata]Id of the gene set (in the DB) to use for selecting the genes used for the module score");
				System.out.println("-metadata %s \t[Required, or -geneset]Name of the metadata (in the Loom) to use for selecting the genes used for the module score");
				System.out.println("-dataset %s \t[Required]Name of the dataset (in the Loom) to use for computing the module score on its values");
				System.out.println("-oAnnot %s \t[Optional] Output metadata for storing the scores (instead of JSON file)");
				System.out.println("-seed %i \tValue of the seed for the pseudo-random generator (default = 42)");
				System.out.println("-nBackgroundGenes %i \tNumber of background gene to take randomly (with replacement) from one features' same bin (default = 100)");
				System.out.println("-nBins %i \tNumber of bins to generate (default = 24)");
				System.out.println("-sel %s \tIn case of a metadata, name of the reference selection");
				if(Parameters.debugMode) System.out.println("-h %s \t\tTo change specific host (default is " + Config.ConfigDEV().getProperty("mDbHost") + ")");
				else System.out.println("-h %s \t\tTo change specific host (default is " + Config.ConfigMAIN().getProperty("mDbHost") + ")");
				break;
			case Preparsing: 
				System.out.println("Preparsing Mode\n\nOptions:");
				System.out.println("-col %s \tName Column [none, first, last].");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-f %s \t\tFile to parse.");
				System.out.println("-organism %i \tId of the organism.");
				System.out.println("-header %b \tThe file has a header [true, false].");
				System.out.println("-sel %s \tIn case of an archive, or a h5 with multiple groups, name of entry to load as a dataset.");
				System.out.println("-d %s \t\tDelimiter.");
				if(Parameters.debugMode) System.out.println("-h %s \t\tTo change specific host (default is " + Config.ConfigDEV().getProperty("mDbHost") + ")");
				else System.out.println("-h %s \t\tTo change specific host (default is " + Config.ConfigMAIN().getProperty("mDbHost") + ")");
				break;
			case Parsing: 
				System.out.println("Parsing Mode\n\nOptions:");
				System.out.println("-ncells %s \t[Required] Number of cells (from preparsing)");
				System.out.println("-ngenes %s \t[Required] Number of genes (from preparsing)");
				System.out.println("-f %s \t\t[Required] File to parse.");
				System.out.println("-type %s \tFile type [RAW_TEXT, LOOM, H5_10x, ARCHIVE, ARCHIVE_COMPRESSED, COMPRESSED]");
				System.out.println("-col %s \tName Column [none, first, last].");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-organism %i \tId of the organism.");
				System.out.println("-header %b \tThe file has a header [true, false].");
				System.out.println("-sel %s \tIn case of an archive, or a h5 with multiple groups, name of entry to load as a dataset.");
				System.out.println("-d %s \t\tDelimiter.");
				break;
			case PreparseMetadata: 
				System.out.println("Preparse Metadata Mode\n\nOptions:");
				System.out.println("-col %s \tName Column [none, first, last].");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-loom %s \t\tLoom file to annotate.");
				System.out.println("-f %s \t\tFile to parse.");
				System.out.println("-header %b \tThe file has a header [true, false].");
				System.out.println("-sel %s \tIn case of an archive, or a h5 with multiple groups, name of entry to load as a dataset.");
				System.out.println("-which %s \tIf you know the metadata type of annotation [gene, cell] (default is CELL).");
				System.out.println("-d %s \t\tDelimiter.");
				break;
			case ParseMetadata: 
				System.out.println("Parse Metadata Mode\n\nOptions:");
				System.out.println("-col %s \tName Column [none, first, last].");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-loom %s \t\tLoom file to annotate.");
				System.out.println("-f %s \t\tFile to parse.");
				System.out.println("-type %s \t\tFile type [RAW_TEXT, LOOM, H5_10x, ARCHIVE, ARCHIVE_COMPRESSED, COMPRESSED] (of note, LOOM and H5_10x is not handled for now)");
				System.out.println("-header %b \tThe file has a header [true, false].");
				System.out.println("-sel %s \tIn case of an archive, or a h5 with multiple groups, name of entry to load as a dataset.");
				System.out.println("-d %s \t\tDelimiter.");
				System.out.println("-which %s \tMetadata type of annotation [gene, cell] (default is CELL).");
				System.out.println("-metadataType %s \tMetadata types [DISCRETE, NUMERIC, STRING], separated by comma if multiple metadata to load (default is to impute them for each metadata).");
				System.out.println("-removeAmbiguous\tWith this option, all ambiguous genes are discarded (default is to take all of them)");
				break;
			case CreateCellSelection:
				System.out.println("CreateCellSelection Mode\n\nOptions:");
				System.out.println("-o %s \t\t[Optional] Output JSON file.");
				System.out.println("-loom %s \tLoom file to annotate.");
				System.out.println("-f %s \t\tThe JSON file containing the selected cell indexes.");
				System.out.println("-meta %s \tPath of the metadata to create.");
				break;
			case GetGeneStats:
				System.out.println("GetGeneStats Mode\n\nOptions:");
				System.out.println("-o %s \t\t[Optional] Output JSON file.");
				System.out.println("--loom %s \t[Required] Loom file to read.");
				System.out.println("--iAnnot %s \t[Optional] Input dataset e.g. '/matrix' (default)");
				System.out.println("--stable-ids %i [One is required] Stable_id(s) of the row to read (separated by comma if multiple)");
				System.out.println("--indexes %i \t[One is required] Index(es) of the row to read (separated by comma if multiple)");
				System.out.println("--names %s \t[One is required] Name(s) of the gene to extract (separated by comma if multiple)");
				System.out.println("--loom-cells %s [Optional] Loom file of the columns/cells to export");
				System.out.println("--metadata %s \t[Optional] Path of input metadata for selecting cells");
				System.out.println("--index %s\t[Optional] Extract only cells of that value.");
				System.out.println("--id %i \t[Optional] Id of the metadata in the table annot (database asap)");
				break;
			case ExtractRow: 
				System.out.println("ExtractRow Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-loom %s \t[Required] Loom file to read.");
				System.out.println("-iAnnot %s \t\tInput dataset e.g. '/matrix' (default)");
				System.out.println("-stable_ids %s \t[One is required] Stable_id(s) of the row to read (separated by comma if multiple)");
				System.out.println("-loom_cells %s \t[Optional] Loom file of the columns/cells to export");
				System.out.println("-indexes %s \t[One is required] Index(es) of the row to read (separated by comma if multiple)");
				System.out.println("-names %s \t[One is required] Name(s) of the gene to extract (separated by comma if multiple)");
				System.out.println("-display-names %s \t[Optional] For adding the name(s) of the cell/genes extracted in the output JSON");
				System.out.println("-display-stable-ids %s \t[Optional] For adding the stable_id(s) of the cell/genes extracted in the output JSON");
				System.out.println("-prec \t[Optional] For trimming the values to this number of significant number after the decimal (default = 3).");
				System.out.println("--sort \t[Optional] For sorting the returned values. Returns an additional vector with the original indexes order.");
				break;
			case ExtractCol: 
				System.out.println("ExtractCol Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-loom %s \t[Required] Loom file to read.");
				System.out.println("-iAnnot %s \tInput dataset e.g. '/matrix' (default)");
				System.out.println("-indexes %s \t[One is required] Index(es) of the col to read (separated by comma if multiple)");
				System.out.println("-names %s \t[One is required] Name(s) of the cell to extract (separated by comma if multiple)");
				System.out.println("-display-names %s \t[Optional] For adding the name(s) of the cell/genes extracted in the output JSON");
				System.out.println("-display-stable-ids %s \t[Optional] For adding the stable_id(s) of the cell/genes extracted in the output JSON");
				System.out.println("-prec %i \t[Optional] For trimming the values to this number of significant number after the decimal (default = 3).");
				break;
			case ExtractDataset: 
				System.out.println("ExtractDataset Mode\n\nOptions:");
				System.out.println("--loom %s \tLoom file to read");
				System.out.println("-o | --output-file %s \tOutput file (default = standard output)");
				System.out.println("-t | --output-type %s \t[JSON, PLAIN_TEXT] Output type (default = JSON)");
				System.out.println("--iAnnot %s \tPath of input dataset (default = '/matrix')");
				System.out.println("--row-names \tFor adding the row names in the output file");
				System.out.println("--col-names \tFor adding the col names in the output file");
				System.out.println("--prec %i \tFor trimming the values to this number of significant number after the decimal (default = 3).");
				System.out.println("--col-indexes %s\tIndex(es) of the columns to export (separated by comma if multiple)");
				System.out.println("--col-stableids %s\tStable_id(s) of the columns to export (separated by comma if multiple)");
				System.out.println("--loom-cells %s\tLoom file of the columns/cells to export (default: same loom)");
				break;
			case DimensionReduction: 
				System.out.println("Dimension Reduction Mode\n\nOptions:");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-f %s \t\tFile to parse.");
				System.out.println("-m %s \t\tModel to use for dimension reduction. It should be one of the following: [pca, mds, tsne, zifa]");
				break;
			case DifferentialExpression: 
				System.out.println("Differential Expression Mode\n\nOptions:");
				System.out.println("-o %s \t\t\tOutput folder.");
				System.out.println("-loom %s \t\tLoom file to annotate.");
				System.out.println("-m %s \t\t\tModel to use for differential expression. It should be one of the following: [wilcox-asap]");
				System.out.println("-iAnnot %s \t\tInput dataset e.g. '/matrix'");
				System.out.println("-oAnnot %s \t\tOutput metadata for results");
				System.out.println("-gAnnot %s \t\tGroup metadata.");
				System.out.println("-g1 %s \t\t\tReference group.");
				System.out.println("-g2 %s \t\t\tComparison group.");
				System.out.println("--is_count_table %s \tSay if the count matrix is raw counts (integer) and should be logged: [true, false]");
				break;
			case Normalization: 
				System.out.println("Normalization Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-loom %s \t\tLoom file to annotate.");
				System.out.println("-oAnnot %s \t\tOutput metadata for normalized dataset");
				break;
			case Scaling: 
				System.out.println("Scaling Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-loom %s \t\tLoom file to annotate.");
				System.out.println("-iAnnot %s \t\tInput dataset e.g. '/layers/toto'");
				System.out.println("-oAnnot %s \t\tOutput metadata for normalized dataset");
				System.out.println("-scaleMax %s \t\tMax scaling factor.");
				System.out.println("-scale %s \t\t[false, true] whether to scale the dataset");
				System.out.println("-center %s \\t\t[false, true] whether to center the dataset");
				break;
			case UpdateEnsemblDB:
				System.out.println("UpdateEnsemblDB\n\nOptions:");
				System.out.println("-o %s \t\tOutput folder.");
				if(Parameters.debugMode) System.out.println("-h %s \t\tTo change specific host (default is " + Config.ConfigDEV().getProperty("mDbHost") + ")");
				else System.out.println("-h %s \t\tTo change specific host (default is " + Config.ConfigMAIN().getProperty("mDbHost") + ")");
				break;
			case CreateKeggDB:
				System.out.println("CreateKeggDB\n\nOptions:");
				System.out.println("-o %s \tOutput folder.");
				break;
			case CreateGODB:
				System.out.println("CreateGODB\n\nOptions:");
				System.out.println("-o %s \tOutput folder.");
				break;
			case RegenerateNewOrganism:
				System.out.println("RegenerateNewOrganism Mode\n\nOptions:");
				System.out.println("-organism %i \tId of the organism.");
				System.out.println("-o %s \t\tOutput folder where are the 'not_found_genes.txt', 'output.json' and 'gene_names.json' files to modify.");
				System.out.println("-j %s \t\tThe JSON file containing the gene names.");
				break;
			case ListMetadata: 
				System.out.println("ListMetadata Mode\n\nOptions:");
				System.out.println("-o %s \t[Required] Output JSON file.");
				System.out.println("-f %s \t[Required] Loom file to read.");
				break;
			case ExtractMetadata: 
				System.out.println("ExtractMetadata Mode\n\nOptions:");
				System.out.println("-o %s\t[Optional] Output JSON file.");
				System.out.println("-loom %s\t[Required] Loom file to read.");
				System.out.println("-meta %s\t[Required or -metaJSON] Name of the metadata to extract.");
				System.out.println("-metaJSON %s\t[Required or -meta] JSON file containing full path of metadata(s) to extract.");
				System.out.println("-no-values\t[Optional] For not displaying the values.");
				System.out.println("-names\t[Optional] For displaying the gene/cell names.");
				System.out.println("-prec %i\t[Optional] For trimming the values to this number of significant number after the decimal (default = 3).");
				System.out.println("--scientific \t[Optional] For outputting the numeric results in scientific notation.");
				System.out.println("-type %s\t\t[Optional] Type of the metadata (for not guessing it) [DISCRETE, NUMERIC, STRING]");
				break;
			case MatchValues: 
				System.out.println("MatchValues Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-loom %s \t[Required] Loom file to read.");
				System.out.println("-iAnnot %s \t[Required] Name of the metadata to extract ids from.");
				System.out.println("-light \t[Optional] For not outputting the id of the genes/cells (just the number).");
				System.out.println("-value %s\t[Required] Extract only rows/cols indexes of that value.");
				break;
			case RemoveMetaData:
				System.out.println("RemoveMetaData Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-loom %s \t[Required] Loom file to modify.");
				System.out.println("-meta %s \t[Required or -metaJSON] Full path of the metadata to remove.");
				System.out.println("-metaJSON %s \t[Required or -meta] JSON file containing full path of metadata(s) to remove.");
				break;
			case CopyMetaData:
				System.out.println("CopyMetaData Mode\n\nOptions:");
				System.out.println("-loomFrom %s \t[Required] Loom file to read from.");
				System.out.println("-loomTo %s \t[Required] Loom file to add the metadata to.");
				System.out.println("-meta %s \t[Required or -metaJSON] Full path of the metadata to copy.");
				System.out.println("-metaJSON %s \t[Required or -meta] JSON file containing full path of metadata(s) to copy.");
				System.out.println("-o %s \t[Optional] Output JSON file containing metadata that were actually copied.");
				break;
			case FilterCols: 
				System.out.println("FilterCols Mode\n\nOptions:");
				System.out.println("-o %s \t[Required] Output folder for generated Loom and JSON files.");
				System.out.println("-loom %s \t[Required] Loom file to read.");
				System.out.println("-col_names_file %s \t[Required one of both] JSON file containing Cells to filter out.");
				System.out.println("-col_indexes_file %s \t[Required one of both] JSON file containing Cells to filter out.");
				break;
			case FilterRows: 
				System.out.println("FilterRows Mode\n\nOptions:");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-loom %s \tFile to parse.");
				System.out.println("-m %s \t\tModel to use for filtering. It should be one of the following: [basic, cpm, keep]");
				System.out.println("-row_names_file %s \t[Required] JSON file containing Genes to Keep.");
				break;
			case FilterDEMetadata:
				System.out.println("FilterDEMetadata Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-loom %s \t[Required] Input Loom file");
				System.out.println("-iAnnot %s \t[Required] Input metadata e.g. '/row_attrs/toto'");
				System.out.println("-p %f \t\tThreshold for nominal p-value [defaut=1]");
				System.out.println("-fdr %f \tThreshold for FDR [defaut=1]");
				System.out.println("-fc %f \t\tThreshold for Fold Change [default=0]");
				System.out.println("-top %i \tSelecting only the top X genes (ranked by FC)");
				System.out.println("-light \tFor not outputting the genes ids (just the number).");
				break;
			case FindMarkers:
				System.out.println("FindMarkers Mode\n\nOptions:");
				System.out.println("-o | --output %s \t[Required] Output folder");
				System.out.println("--loom %s \t\t[Required] Input Loom file");
				System.out.println("--iAnnot %s \t\t[Required] Input metadata e.g. '/row_attrs/toto'");
				System.out.println("--id %i \t\t[Required] Id of the metadata in the table annot (database asap)");
				System.out.println("--is_count_table %s \tSay if the count matrix is raw counts (integer) and should be logged: [true, false]");
				break;
			case IndexByCell:
				System.out.println("IndexByCell Mode\n\nOptions:");
				System.out.println("-o | --output %s \t[Optional] Output JSON file.");
				System.out.println("--loom %s \t[Required] Input Loom file");
				System.out.println("--loomIndex %s \t[Required] Loom file containing index");
				System.out.println("--iAnnot %s \t[Required] Input metadata e.g. '/row_attrs/toto'");
				System.out.println("--id %i \t[Required] Id of the metadata in the table annot (database asap)");
				break;
			case GetIndex: 
				System.out.println("GetIndex Mode\n\nOptions:");
				System.out.println("-o | --output %s \t[Optional] Output JSON file.");
				System.out.println("--loomIndex %s \t[Required] Loom file containing index");
				System.out.println("--jsonListCells %s \t[Required] JSON file containing cells to consider");
				break;
		}
		System.out.println();
	}
}
