package model;

import java.io.File;

import db.DBManager;
import enrichment.model.Model;
import parsing.model.ColumnName;

public class Parameters 
{
	// Shared
	public static String organism_S = null;
	public static String outputFolder = null;
	public static String fileName = null;
	public static int organism = 1;
	
	// CreateDLFile
	public static String outputFile= null;
	public static String JSONFileName = null;
	
	// Enrichment
	public static int nbRepeat = -1;
	public static int randomSeed = 42;
	public static boolean isSilent = false;
	
	public static String adjMethod = "fdr";
	public static double probaCutoff = -1;
	public static Model model = null;

	public static int maxGenesInPathway = 500;
	public static int minGenesInPathway = 15;

	public static String pathwayFile = null;
	public static String backgroundFile = null;
	public static String listGenesFile = null;
	
	// Parsing
	public static boolean has_header = true;
	public static int skip_line = 0;
	public static ColumnName name_column = ColumnName.FIRST;
	public static String delimiter = "\t";
	
	public static void load(String[] args, Mode m)
	{
		if(args.length == 0)
		{
			printHelp(m);
			System.exit(0);
		}
		switch(m)
		{
			case CreateEnrichmentDB:
			case CreateGenesDB:
				// Nothing to do
				break;
			case CreateEnsemblDB:
				loadCreateEnsemblDB(args);
				break;
			case Enrichment:
				loadEnrichment(args);
				break;
			case Parsing: 
				loadParsing(args);
				break;
			case RegenerateNewOrganism:
				loadRegenerateNewOrganism(args);
				break;
			case CreateDLFile:
				loadCreateDLFile(args);
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
							model = Model.GSEA;
							break;
						case "hypergeo":
							model = Model.HyperGeometric;
							break;
						case "fet":
							model = Model.FET;
							break;
						default:
							System.err.println("The entered model, "+args[i]+ ", does not exist!\nIt should be one of the following: [gsea, hypergeo, fet]");
							System.exit(-1);
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
							System.err.println("The '-n' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
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
							System.err.println("The '-p' option should be followed by an Double value in [0, 1]. You entered " + args[i]);
							System.exit(-1);
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
							System.err.println("The '-s' option should be followed by an Integer value. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-adj":
						i++;
						adjMethod = args[i];
						if(!adjMethod.equals("bonferroni") && !adjMethod.equals("fdr") && !adjMethod.equals("none"))
						{
							System.err.println("The '-adj' option should be followed by one of those values: [bonferroni, fdr, none]. You entered " + args[i]);
							System.exit(-1);
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
							System.err.println("The '-max' option should be followed by an positive Integer value. You entered " + args[i]);
							System.exit(-1);
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
							System.err.println("The '-min' option should be followed by an positive Integer value. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		
		if(model == null || pathwayFile == null || outputFolder == null || backgroundFile == null || listGenesFile == null)
		{
			if(model == null) System.err.println("No model is specified, please choose a model by using the '-m' option.");
			if(pathwayFile == null) System.err.println("No pathway file is specified, please choose a data file by using the '-path' option.");
			if(outputFolder == null) System.err.println("No output folder is specified, please choose an output file by using the '-o' option.");
			if(backgroundFile == null) System.err.println("No background file is specified, please choose a background file by using the '-background' option.");
			if(listGenesFile == null) System.err.println("No gene list file is specified, please choose a gene list file by using the '-test' option.");
			System.out.println();
			if(!Parameters.isSilent) System.err.println("Note: Enrichment cannot be run because files are missing.");
			System.exit(-1);
		}
		new File(outputFolder).mkdirs();
		if(model == Model.GSEA && nbRepeat == -1)
		{
			System.out.println("The model '" + model + "' is using a permutation resampling, but you did not specify the number of permutation to perform using the '-n' option. Using default number: 10000.");
			nbRepeat = 10000;
		}
		if((model == Model.FET || model == Model.HyperGeometric) && probaCutoff == -1)
		{
			System.out.println("The model '" + model + "' requires a probability cutoff for considering a gene as deregulated, but you did not specify this number by using the '-p' option. Using default value: 0.05.");
			probaCutoff = 0.05;
		}
	}
	
	public static void loadCreateEnsemblDB(String[] args)
	{
		boolean found = false;
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
						organism_S = args[i];
						for(String spe:DBManager.species) if(spe.equals(organism_S)) {found = true; break;}
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFolder == null || organism_S == null)
		{
			printHelp(Mode.CreateEnsemblDB);
			System.exit(-1);
		}
		if(!found)
		{
			System.err.println("This organism (" + organism_S + ") was not found in our database.");
			System.exit(-1);
		}
	}
	
	public static void loadCreateDLFile(String[] args)
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
						break;
					case "-f":
						i++;
						fileName = args[i];
						fileName = fileName.replaceAll("\\\\", "/");
						break;
					case "-j":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
						break;
					default:
						System.err.println("Unused argument: " + arg);
				}
			}
		}
		if(outputFile == null || fileName == null || JSONFileName == null)
		{
			printHelp(Mode.CreateDLFile);
			System.exit(-1);
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
							System.err.println("The '-organism' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-j":
						i++;
						JSONFileName = args[i];
						JSONFileName = JSONFileName.replaceAll("\\\\", "/");
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
					case "-skip":
						i++;
						try
						{
							skip_line = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-s' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
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
							System.err.println("The '-organism' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
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
								System.err.println("The '-col' option should be followed by [last, first, none]. You entered " + args[i]);
								System.exit(-1);
						}
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
	}
	
	public static void printHelp(Mode m)
	{
		switch(m)
		{
			case CreateEnrichmentDB:
				break;
			case CreateGenesDB:
				break;
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
			case Parsing: 
				System.out.println("Parsing Mode\n\nOptions:");
				System.out.println("-col %s \tName Column [none, first, last].");
				System.out.println("-o %s \t\tOutput folder.");
				System.out.println("-f %s \t\tFile to parse.");
				System.out.println("-organism %i \tId of the organism.");
				System.out.println("-header %b \tThe file has a header [true, false].");
				System.out.println("-d %s \t\tDelimiter.");
				System.out.println("-skip %i \tNumber of lines to skip at the beginning of the file.");
				break;
			case CreateEnsemblDB:
				System.out.println("CreateEnsemblDB\n\nOptions:");
				System.out.println("-o %s \tOutput folder.");
				System.out.println("-organism %s \tEnsembl name of the organism [mus_musculus, homo_sapiens, ...].");
				break;
			case RegenerateNewOrganism:
				System.out.println("RegenerateNewOrganism Mode\n\nOptions:");
				System.out.println("-organism %i \tId of the organism.");
				System.out.println("-o %s \t\tOutput folder where are the 'not_found_genes.txt', 'output.json' and 'gene_names.json' files to modify.");
				System.out.println("-j %s \t\tThe JSON file containing the gene names.");
				break;
			case CreateDLFile:
				System.out.println("CreateDLFile Mode\n\nOptions:");
				System.out.println("-f %s \t\tThe matrix file with id of genes.");
				System.out.println("-o %s \t\tThe output file to create.");
				System.out.println("-j %s \t\tThe JSON file containing the gene names.");
				break;
		}
		System.out.println();
	}
}
