package db;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class DrugBankParser // Get XML from here https://www.drugbank.ca/releases/latest
{
	//public static ArrayList<String> organisms = new ArrayList<String>();
	public static HashSet<String> accepted_organisms;
	public static HashMap<String, String> drug_description = new HashMap<>();
	
	public static void setOrganism(String organism)
	{
		accepted_organisms = new HashSet<String>();
		switch(organism)
		{
		case "hsa":
			accepted_organisms.add("Human");
			accepted_organisms.add("human");
			break;
		case "mmu":
			accepted_organisms.add("Mouse");
			break;
		case "ecoli": // TODO correct name
			accepted_organisms.add("Escherichia coli");
			accepted_organisms.add("Escherichia coli (strain K12)");
			break;
		case "rat": // TODO correct name
			accepted_organisms.add("Rat");
			break;
		default:
			System.err.println("This organism does not exist");
			System.exit(-1);
		}
	}
	
	public static void parse(String organism, String associationFile, String outputGMTFile) throws IOException 
	{
		setOrganism(organism);
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputGMTFile));
		HashMap<String, ArrayList<String>> assoc = loadDataAssoc(associationFile); // Parse XML
		//System.out.print(organisms.size() + " organisms were found in XML file\n[\"organisms\":");
		//for(String o:organisms) System.out.print(o + ",");
		//System.out.println("]");
		// Create GMT from parsed XML
		for(String geneset:assoc.keySet())
		{
			ArrayList<String> genes = assoc.get(geneset);
			if(genes.size() > 0) // Do not write empty genesets
			{
				bw.write(geneset+"\t"+drug_description.get(geneset)+"\thttp://www.drugbank.ca/drugs/"+geneset);
				for (String gene:genes) bw.write("\t"+gene); // List gene
				bw.write("\n");
			}
		}
		bw.close();
	}
	
	private static HashMap<String, ArrayList<String>> loadDataAssoc(String filename)
	{
		if(filename == null) return null;
		HashMap<String, ArrayList<String>> data_assoc = new HashMap<>();
		try
		{
			File xmlFile = new File(filename);
			DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
			Document doc = dBuilder.parse(xmlFile);
			doc.getDocumentElement().normalize();
			NodeList drugs = doc.getFirstChild().getChildNodes();
			for (int d = 0; d < drugs.getLength(); d++) 
			{
				Node dNode = drugs.item(d);
				if (dNode.getNodeName().equals("drug") && dNode.getNodeType() == Node.ELEMENT_NODE) 
				{
					ArrayList<String> genes = new ArrayList<String>();
					Element drug = (Element) dNode;
					String drugId = drug.getElementsByTagName("drugbank-id").item(0).getTextContent();
					drug_description.put(drugId, drug.getElementsByTagName("name").item(0).getTextContent());
					NodeList subnodes = drug.getChildNodes();
					for (int s = 0; s < subnodes.getLength(); s++) 
					{
						Node sNode = subnodes.item(s);
						if (sNode.getNodeName().equals("targets") && sNode.getNodeType() == Node.ELEMENT_NODE) 
						{
							Element targets = (Element)sNode;
							NodeList targetSet = targets.getChildNodes();
							for (int t = 0; t < targetSet.getLength(); t++) 
							{
								Node tNode = targetSet.item(t);
								if (tNode.getNodeName().equals("target") && tNode.getNodeType() == Node.ELEMENT_NODE) 
								{
									Element target = (Element)tNode;
									String organism = target.getElementsByTagName("organism").item(0).getTextContent();
									//if(!organisms.contains(organism) && organism != "") organisms.add(organism);
									if(accepted_organisms.contains(organism.trim())) // Restrict to accepted genes
									{
										NodeList nodes = target.getChildNodes();
										for (int tt = 0; tt < nodes.getLength(); tt++) 
										{
											Node n = nodes.item(tt);
											if (n.getNodeName().equals("polypeptide") && n.getNodeType() == Node.ELEMENT_NODE) 
											{
												Element polypeptide = (Element)n;
												String geneId = polypeptide.getElementsByTagName("gene-name").item(0).getTextContent().trim();
												if(!genes.contains(geneId) && !geneId.equals("")) genes.add(geneId);
											}	
										}
									}
								}
							}
						}
					}
					data_assoc.put(drugId, genes);
				}
			}
		}
		catch(FileNotFoundException nfe)
		{
			nfe.printStackTrace();
			System.exit(-1);
		} 
		catch (ParserConfigurationException e) 
		{
			e.printStackTrace();
			System.exit(-1);
		} 
		catch (SAXException e) 
		{
			e.printStackTrace();
			System.exit(-1);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
			System.exit(-1);
		}
		return data_assoc;
	}
}