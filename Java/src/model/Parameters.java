package model;

import java.io.File;

import config.Config;
import db.DBManager;
import json.ErrorJSON;
import parsing.model.ColumnName;
import parsing.model.FileType;
import tools.Utils;

public class Parameters 
{
	// Shared
	public static String organism_S = null;
	public static String outputFolder = null;
	public static String outputFile = null;
	public static String group = null;
	public static String fileName = null;
	public static String selection = null;
	public static FileType fileType = null;
	public static String fitModel = null;
	public static String erccFile = null;
	public static int organism = 1;
	public static int taxon = -1;
	public static String configFile = "asap.conf";
	public static long nCells = -1;
	public static long nGenes = -1;
	public static long index = -1;
	public static String geneName = null;
	public static String cellName = null;
	public static String metaName = null;
	public static MetaOn which = null;
	public static Metatype metatype = null;
	public static boolean details = true;
	public static boolean evenLessDetails = false;
	public static String loomFile = null;
	public static String loomFile2 = null;
	public static double idleTime = 0;
	public static boolean isIndex = false;
	
	// Dimension Reduction
	public static dim_reduction.model.Model dimReducModel = null;
	// TSNE
	public static int perplexity = -1;
	
	// DE
	public static differential_expression.Model deModel = null;
	public static String iAnnot = null;
	public static String oAnnot = null;
	public static String gAnnot = null;
	public static String group_1 = null;
	public static String group_2 = null;
	
	// Filtering
	public static filtering.model.Model filtModel = null;
	// CPM
	public static int nbCountsPerCell = -1;
	public static int nbCellsDetected = -1;
	
	// CreateDLFile
	public static String JSONFileName = null;
	
	// Enrichment
	public static int nbRepeat = -1;
	public static int randomSeed = 42;
	public static boolean isSilent = false;
	
	public static String adjMethod = "fdr";
	public static double probaCutoff = -1;
	public static enrichment.model.Model enrichModel = null;
	
	public static int maxGenesInPathway = 500;
	public static int minGenesInPathway = 15;

	public static String pathwayFile = null;
	public static String backgroundFile = null;
	public static String listGenesFile = null;
	
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
			case CreateEnrichmentDB:
				loadCreateEnrichmentDB(args);
				break;
			case UpdateEnsemblDB:
				loadUpdateEnsemblDB(args);
				break;
			case Enrichment:
				loadEnrichment(args);
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
			case FilterGenes: 
				loadFilterGenes(args);
				break;
			case DimensionReduction: 
				loadDimensionReduction(args);
				break;
			case DifferentialExpression: 
				loadDifferentialExpression(args);
				break;
			case RegenerateNewOrganism:
				loadRegenerateNewOrganism(args);
				break;
			case ExtractRow:
				loadExtractRow(args);
				break;
			case ExtractCol:
				loadExtractCol(args);
				break;
			case ListMetaData:
				loadListMetadata(args);
				break;
			case ExtractMetaData:
				loadExtractMetadata(args);
				break;
			case RemoveMetaData:
				loadRemoveMetadata(args);
				break;
			case CopyMetaData:
				loadCopyMetadata(args);
				break;
			case FilterCells:
				loadFilterCells(args);
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
						case "gsea":
							enrichModel = enrichment.model.Model.GSEA;
							break;
						case "hypergeo":
							enrichModel = enrichment.model.Model.HyperGeometric;
							break;
						case "fet":
							enrichModel = enrichment.model.Model.FET;
							break;
						default:
							new ErrorJSON("The entered model, "+args[i]+ ", does not exist!\nIt should be one of the following: [gsea, hypergeo, fet]");
						}
						break;
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						new File(outputFolder).mkdirs();
						break;
					case "-path":
						i++;
						pathwayFile = args[i];
						break;
					case "-background":
						i++;
						backgroundFile = args[i];
						break;
					case "-test":
						i++;
						listGenesFile = args[i];
						break;
					case "-n":
						i++;
						try
						{
							nbRepeat = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-n' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-p":
						i++;
						try
						{
							probaCutoff = Double.parseDouble(args[i]);
							if(probaCutoff <0 || probaCutoff >1) throw new NumberFormatException();
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-p' option should be followed by an Double value in [0, 1]. You entered " + args[i]);
						}
						break;
					case "-s":
						i++;
						try
						{
							randomSeed = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-s' option should be followed by an Integer value. You entered " + args[i]);
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
					case "-silent":
						isSilent = true;
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
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(enrichModel == null) new ErrorJSON("No model is specified, please choose a model by using the '-m' option.");
		if(pathwayFile == null) new ErrorJSON("No pathway file is specified, please choose a data file by using the '-path' option.");
		if(outputFolder == null) new ErrorJSON("No output folder is specified, please choose an output file by using the '-o' option.");
		if(backgroundFile == null) new ErrorJSON("No background file is specified, please choose a background file by using the '-background' option.");
		if(listGenesFile == null) new ErrorJSON("No gene list file is specified, please choose a gene list file by using the '-test' option.");
		
		new File(outputFolder).mkdirs();
		if(enrichModel == enrichment.model.Model.GSEA && nbRepeat == -1)
		{
			System.out.println("The model '" + enrichModel + "' is using a permutation resampling, but you did not specify the number of permutation to perform using the '-n' option. Using default number: 10000.");
			nbRepeat = 10000;
		}
		if((enrichModel == enrichment.model.Model.FET || enrichModel == enrichment.model.Model.HyperGeometric) && probaCutoff == -1)
		{
			System.out.println("The model '" + enrichModel + "' requires a probability cutoff for considering a gene as deregulated, but you did not specify this number by using the '-p' option. Using default value: 0.05.");
			probaCutoff = 0.05;
		}
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
						if(DBManager.JDBC_DRIVER.equals("com.mysql.jdbc.Driver")) DBManager.URL = "jdbc:mysql://";
						else if(DBManager.JDBC_DRIVER.equals("org.postgresql.Driver")) DBManager.URL = "jdbc:postgresql://";
						DBManager.URL += args[i] + "?user=" + Config.getProperty("mDbUser") + "&password=" + Config.getProperty("mDbPwds");
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
	
	public static void loadCreateEnrichmentDB(String[] args)
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
							DBManager.connect();
							taxon = DBManager.getTaxonFromOrganismID(organism);
							DBManager.disconnect();
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-organism' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null)
		{
			printHelp(Mode.CreateEnrichmentDB);
			new ErrorJSON("Please specify an output folder using -o option");
		}
		if(taxon == -1)
		{
			printHelp(Mode.CreateEnrichmentDB);
			new ErrorJSON("Please specify an organism_id using -organism option");
		}
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
					case "-f":
						i++;
						fileName = args[i];
						fileName = fileName.replaceAll("\\\\", "/");
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
					case "-i":
						i++;
						try
						{
							index = Long.parseLong(args[i]);
							if(index <0) throw new NumberFormatException();
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-i' option should be followed by an positive Integer value. You entered " + args[i]);
						}
						break;
					case "-gene":
						i++;
						geneName = args[i].trim().toUpperCase();
						if(geneName.equals("")) new ErrorJSON("The -gene option should be followed by a gene name. You entered " + args[i]);
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(fileName == null)
		{
			printHelp(Mode.ExtractRow);
			System.exit(-1);
		}
		if(iAnnot == null) iAnnot = "/matrix";
		if(index == -1 && geneName == null) new ErrorJSON("You need to use at least one of the two options -gene or -i");
		if(index != -1 && geneName != null) new ErrorJSON("You need to use only one of the two options -gene or -i, not both");
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
					case "-f":
						i++;
						fileName = args[i];
						fileName = fileName.replaceAll("\\\\", "/");
						break;
					case "-iAnnot":
						i++;
						iAnnot = args[i];
						if(!iAnnot.startsWith("/")) iAnnot = "/" + iAnnot;
						break;
					case "-i":
						i++;
						try
						{
							index = Long.parseLong(args[i]);
							if(index <0) throw new NumberFormatException();
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-i' option should be followed by an positive Integer value. You entered " + args[i]);
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
					case "-cell":
						i++;
						cellName = args[i].trim();
						if(cellName.equals("")) new ErrorJSON("The -cell option should be followed by a cell name. You entered " + args[i]);
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(fileName == null)
		{
			printHelp(Mode.ExtractCol);
			System.exit(-1);
		}
		if(iAnnot == null) iAnnot = "/matrix";
		if(index == -1 && cellName == null) new ErrorJSON("You need to use at least one of the two options -cell or -i");
		if(index != -1 && cellName != null) new ErrorJSON("You need to use only one of the two options -cell or -i, not both");
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
			printHelp(Mode.ListMetaData);
			System.exit(-1);
		}
	}
	
	public static void loadExtractMetadata(String[] args)
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
					case "-meta":
						i++;
						metaName = args[i];
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
							Utils.changeFormatter(Integer.parseInt(args[i]));
						}
						catch(NumberFormatException nfe)
						{
							new ErrorJSON("The '-prec' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "-light":
						details = false;
						break;
					case "-ultralight":
						evenLessDetails = true;
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(fileName == null || metaName == null)
		{
			printHelp(Mode.ExtractMetaData);
			System.exit(-1);
		}
		//if(!details) System.out.println("'-light' option activated. In case it is a cell metadata, the cell names will NOT be outputted.");
		//if(evenLessDetails) System.out.println("'-ultralight' option activated. Neither the values nor the cells names will be outputted.");
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
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
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
					case "-meta":
						i++;
						metaName = args[i];
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null)
		{
			printHelp(Mode.RemoveMetaData);
			new ErrorJSON("Please specify an output folder using -o option");
		}
		if(loomFile == null)
		{
			printHelp(Mode.RemoveMetaData);
			new ErrorJSON("Please specify a loom file using -loom option");
		}
		if(metaName == null)
		{
			printHelp(Mode.RemoveMetaData);
			new ErrorJSON("Please specify a dataset using -meta option");
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
			printHelp(Mode.CopyMetaData);
			new ErrorJSON("Please specify a loom file to copy FROM using -loomFrom option");
		}
		if(loomFile == null)
		{
			printHelp(Mode.CopyMetaData);
			new ErrorJSON("Please specify a loom file to copy TO using -loomTo option");
		}
		if(metaName == null)
		{
			printHelp(Mode.CopyMetaData);
			new ErrorJSON("Please specify a dataset using -meta option");
		}
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
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
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
			printHelp(Mode.FilterCells);
			new ErrorJSON("Please specify an output folder using -o option");
		}
		if(loomFile == null)
		{
			printHelp(Mode.FilterCells);
			new ErrorJSON("Please specify a loom file to filter using -loom option");
		}
		if(JSONFileName == null)
		{
			printHelp(Mode.FilterCells);
			new ErrorJSON("Please specify a JSON file containing Cells to filter out using -col_names_file option");
		}
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
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
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
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null || fileName == null)
		{
			printHelp(Mode.Parsing);
			System.exit(-1);
		}
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
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
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
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null || loomFile == null || fileName == null)
		{
			printHelp(Mode.Parsing);
			System.exit(-1);
		}
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
						if(DBManager.JDBC_DRIVER.equals("com.mysql.jdbc.Driver")) DBManager.URL = "jdbc:mysql://";
						else if(DBManager.JDBC_DRIVER.equals("org.postgresql.Driver")) DBManager.URL = "jdbc:postgresql://";
						DBManager.URL += args[i] + "?user=" + Config.getProperty("mDbUser") + "&password=" + Config.getProperty("mDbPwds");
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
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						new File(outputFolder).mkdirs();
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
		if(which != null) System.out.println("You specified that the metadata was applied on " + which);
		else 
		{
			System.out.println("You DID NOT specify if the metadata was applied on [cell,gene] using the -which option. Set default to CELL.");
			which = MetaOn.CELL;
		}
		if(outputFolder == null) new ErrorJSON("The output folder '-o' is mandatory for PreparseMetadata.");
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
			printHelp(Mode.FilterGenes);
			String error = "Filtering cannot be run because parameters are missing:\n";
			if(filtModel == null) error += "No model is specified, please choose a model by using the '-m' option.\n";
			if(loomFile == null) error += "No file is specified, please choose a data file by using the '-f' option.\n";
			if(outputFolder == null) error += "No output folder is specified, please choose an output file by using the '-o' option.\n";
			new ErrorJSON(error); // This stops the program
		}
		if(filtModel == filtering.model.Model.KEEP && JSONFileName == null)
		{
			printHelp(Mode.FilterGenes);
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
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
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
						case "wilcox-asap":
							deModel = differential_expression.Model.Wilcoxon;
							break;
						default:
							new ErrorJSON("The entered model, "+args[i]+ ", does not exist!\nIt should be one of the following: [wilcox-asap]");
						}
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
	
	public static void printHelp(Mode m)
	{
		switch(m)
		{
			case Enrichment: 
				System.out.println("Enrichment Mode\n\nOptions:");
				System.out.println("-m %s \t\tChoose a model among [gsea, hypergeo, fet].");
				System.out.println("-o %s \t\tOutput folder");
				System.out.println("-path %s \tPathway/Gene Mapping file.");
				System.out.println("-background %s \tBackground file [matrix].");
				System.out.println("-test %s \tList genes to enrich file [JSON].");
				System.out.println("-n %i \t\tNumber of permutation resampling to perform for models [gsea].");
				System.out.println("-p %f \t\tProbability threshold for considering a gene as deregulated for models [fet, hypergeo].");
				System.out.println("-s %i \t\tRandom seed for the generator of pseudo-random numbers [default = 42].");
				System.out.println("-adj %s \tStatitical adjustment method for multiple comparision [bonferroni, fdr, or none, default = fdr].");
				System.out.println("-min %i \tMinimum number of genes in a pathway for being taken into consideration [default = 15].");
				System.out.println("-max %i \tMaximum number of genes in a pathway for being taken into consideration [default = 500].");
				System.out.println("-silent \tDo not print message in the standard output.");
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
				System.out.println("-h %s \t\tTo change specific host (default is " + Config.getProperty("mDbHost") + ")");
				break;
			case Parsing: 
				System.out.println("Parsing Mode\n\nOptions:");
				System.out.println("-col %s \tName Column [none, first, last].");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-f %s \t\tFile to parse.");
				System.out.println("-type %s \t\tFile type [RAW_TEXT, LOOM, H5_10x, ARCHIVE, ARCHIVE_COMPRESSED, COMPRESSED]");
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
				System.out.println("-type %s \t\tFile type [RAW_TEXT, LOOM, H5_10x, ARCHIVE, ARCHIVE_COMPRESSED, COMPRESSED]");
				System.out.println("-header %b \tThe file has a header [true, false].");
				System.out.println("-sel %s \tIn case of an archive, or a h5 with multiple groups, name of entry to load as a dataset.");
				System.out.println("-d %s \t\tDelimiter.");
				break;
			case CreateCellSelection:
				System.out.println("CreateCellSelection Mode\n\nOptions:");
				System.out.println("-o %s \t\t[Optional] Output JSON file.");
				System.out.println("-loom %s \tLoom file to annotate.");
				System.out.println("-f %s \t\tThe JSON file containing the selected cell indexes.");
				System.out.println("-meta %s \tPath of the metadata to create.");
				break;
			case ExtractRow: 
				System.out.println("ExtractRow Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-f %s \t[Required] Loom file to read.");
				System.out.println("-iAnnot %s \t\tInput dataset e.g. '/matrix' (default)");
				System.out.println("-i %s \t[One is required] Index of the row to read");
				System.out.println("-gene %s \t[One is required] Name of the gene to extract");
				System.out.println("-prec \t[Optional] For trimming the values to this number of significant number after the decimal (default = 3).");
				break;
			case ExtractCol: 
				System.out.println("ExtractCol Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-f %s \t[Required] Loom file to read.");
				System.out.println("-iAnnot %s \tInput dataset e.g. '/matrix' (default)");
				System.out.println("-i %s \t[One is required] Index of the col to read");
				System.out.println("-cell %s \t[One is required] Name of the cell to extract");
				System.out.println("-prec \t[Optional] For trimming the values to this number of significant number after the decimal (default = 3).");
				break;
			case DimensionReduction: 
				System.out.println("Dimension Reduction Mode\n\nOptions:");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-f %s \t\tFile to parse.");
				System.out.println("-m %s \t\tModel to use for dimension reduction. It should be one of the following: [pca, mds, tsne, zifa]");
				break;
			case DifferentialExpression: 
				System.out.println("Differential Expression Mode\n\nOptions:");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-loom %s \t\tLoom file to annotate.");
				System.out.println("-m %s \t\tModel to use for differential expression. It should be one of the following: [wilcox-asap]");
				System.out.println("-iAnnot %s \t\tInput dataset e.g. '/matrix'");
				System.out.println("-oAnnot %s \t\tOutput metadata for results");
				System.out.println("-gAnnot %s \t\tGroup metadata.");
				System.out.println("-g1 %s \t\tReference group.");
				System.out.println("-g2 %s \t\tComparison group.");
				break;
			case UpdateEnsemblDB:
				System.out.println("UpdateEnsemblDB\n\nOptions:");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-h %s \t\tTo change specific host (default is " + Config.getProperty("mDbHost") + ")");
				break;
			case CreateEnrichmentDB:
				System.out.println("CreateEnrichmentDB\n\nOptions:");
				System.out.println("-o %s \tOutput folder.");
				System.out.println("-organism %i \tId of the organism.");
				break;
			case RegenerateNewOrganism:
				System.out.println("RegenerateNewOrganism Mode\n\nOptions:");
				System.out.println("-organism %i \tId of the organism.");
				System.out.println("-o %s \t\tOutput folder where are the 'not_found_genes.txt', 'output.json' and 'gene_names.json' files to modify.");
				System.out.println("-j %s \t\tThe JSON file containing the gene names.");
				break;
			case ListMetaData: 
				System.out.println("ListMetadata Mode\n\nOptions:");
				System.out.println("-o %s \t[Required] Output JSON file.");
				System.out.println("-f %s \t[Required] Loom file to read.");
				break;
			case ExtractMetaData: 
				System.out.println("ExtractMetaData Mode\n\nOptions:");
				System.out.println("-o %s \t[Optional] Output JSON file.");
				System.out.println("-f %s \t[Required] Loom file to read.");
				System.out.println("-meta %s \t[Required] Name of the metadata to extract.");
				System.out.println("-light \t[Optional] For not outputting the gene/cell names.");
				System.out.println("-ultralight \t[Optional] For not outputting the gene/cell names nor the values.");
				System.out.println("-prec \t[Optional] For trimming the values to this number of significant number after the decimal (default = 3).");
				System.out.println("-type %s \t\t[Optional] Type of the metadata (for not guessing it) [DISCRETE, NUMERIC, STRING]");
				break;
			case RemoveMetaData:
				System.out.println("RemoveMetaData Mode\n\nOptions:");
				System.out.println("-o %s \t[Required] Output folder for JSON file.");
				System.out.println("-loom %s \t[Required] Loom file to modify.");
				System.out.println("-meta %s \t[Required] Full path of the metadata to remove.");
				break;
			case CopyMetaData:
				System.out.println("CopyMetaData Mode\n\nOptions:");
				System.out.println("-loomFrom %s \t[Required] Loom file to read from.");
				System.out.println("-loomTo %s \t[Required] Loom file to add the metadata to.");
				System.out.println("-meta %s \t[Required] Full path of the metadata to copy.");
				break;
			case FilterCells: 
				System.out.println("FilterCells Mode\n\nOptions:");
				System.out.println("-o %s \t[Required] Output folder for generated Loom and JSON files.");
				System.out.println("-loom %s \t[Required] Loom file to read.");
				System.out.println("-col_names_file %s \t[Required one of both] JSON file containing Cells to filter out.");
				System.out.println("-col_indexes_file %s \t[Required one of both] JSON file containing Cells to filter out.");
				break;
			case FilterGenes: 
				System.out.println("FilterGenes Mode\n\nOptions:");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-loom %s \tFile to parse.");
				System.out.println("-m %s \t\tModel to use for filtering. It should be one of the following: [basic, cpm, keep]");
				System.out.println("-row_names_file %s \t[Required] JSON file containing Genes to Keep.");
				break;
		}
		System.out.println();
	}
}
