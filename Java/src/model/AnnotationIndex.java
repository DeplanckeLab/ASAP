package model;

import db.DBManager;
import json.ErrorJSON;

public class AnnotationIndex 
{
	public String idMeta;
	public int[] categories;
	
	public AnnotationIndex(String idMeta)
	{
		this.idMeta = idMeta;
		int nbCategories = DBManager.getNbCatJSON(this.idMeta);
		if(nbCategories == -1) new ErrorJSON("Error accessing the categories with metadata id = " + this.idMeta);
		this.categories = new int[nbCategories]; // with indexes 1, 2, 3, 4...
	}
	
	public void add(String cat)
	{
		int catIndex = Integer.parseInt(cat) - 1;
		this.categories[catIndex]++;
	}
}
