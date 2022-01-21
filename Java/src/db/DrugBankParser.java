package db;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

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
	public static HashMap<Integer, StringBuilder> output = new HashMap<>();
	
	public static void main(String[] args)
	{
		//DBManager.JDBC_DRIVER = Config.getProperty("mDbDriv");
		//if(DBManager.JDBC_DRIVER.equals("com.mysql.jdbc.Driver")) DBManager.URL = "jdbc:mysql://";
		//else if(DBManager.JDBC_DRIVER.equals("org.postgresql.Driver")) DBManager.URL = "jdbc:postgresql://";
		//DBManager.URL += Config.getProperty("mDbHost") + "?user=" + Config.getProperty("mDbUser") + "&password=" + Config.getProperty("mDbPwds");
		
		
		DBManager.connect();
		HashMap<Integer, Integer> taxonsInDB = DBManager.listAllOrganisms();
		System.out.println("There are " + taxonsInDB.size() + " organisms in ASAP DB.");
		DBManager.disconnect();
		
		for(Integer tax:taxonsInDB.keySet()) output.put(tax, new StringBuilder());
		
		parse("C:/Users/gardeux/Downloads/full database.xml");
		
		// Create GMT from found species
		try
		{
			for(Integer tax:output.keySet())
			{
				String content = output.get(tax).toString();
				if(!content.equals(""))
				{
					BufferedWriter bw = new BufferedWriter(new FileWriter("C:/Users/gardeux/Desktop/Drugbank/drugbank."+tax+".gmt"));
					bw.write(content);
					bw.close();
				}
			}
		}
		catch(IOException ioe)
		{
			System.err.println(ioe.getMessage());
		}
	}
	
	private static void parse(String filename)
	{
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
					HashMap<Integer, String> taxonXgenes = new HashMap<Integer, String>();
					Element drug = (Element) dNode;
					String drugId = drug.getElementsByTagName("drugbank-id").item(0).getTextContent();
					String drug_description = drug.getElementsByTagName("name").item(0).getTextContent();
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
									NodeList nodes = target.getChildNodes();
									for (int tt = 0; tt < nodes.getLength(); tt++) 
									{
										Node n = nodes.item(tt);
										if (n.getNodeName().equals("polypeptide") && n.getNodeType() == Node.ELEMENT_NODE) 
										{
											Element polypeptide = (Element)n;
											
											String taxId = polypeptide.getElementsByTagName("organism").item(0).getAttributes().getNamedItem("ncbi-taxonomy-id").getTextContent().trim();
											if(!taxId.equals(""))
											{
												int taxon = Integer.parseInt(taxId);
												String geneId = polypeptide.getElementsByTagName("gene-name").item(0).getTextContent().trim();
												
												String genes = taxonXgenes.get(taxon);
												if(genes == null) taxonXgenes.put(taxon, geneId);
												else taxonXgenes.put(taxon, genes + "\t" + geneId);
											}
										}	
									}
								}
							}
						}
					}
					
					for(Integer tax:taxonXgenes.keySet())
					{
						if(output.get(tax) != null) output.get(tax).append(drugId).append("\t").append(drug_description).append("\thttp://www.drugbank.ca/drugs/").append(drugId).append("\t").append(taxonXgenes.get(tax)).append("\n");
					}
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
	}
}