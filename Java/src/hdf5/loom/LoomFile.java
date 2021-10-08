package hdf5.loom;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import bigarrays.DoubleArray64;
import bigarrays.FloatArray64;
import bigarrays.IntArray64;
import bigarrays.LongArray64;
import bigarrays.StringArray64;
import ch.systemsx.cisd.hdf5.HDF5DataClass;
import ch.systemsx.cisd.hdf5.HDF5DataSetInformation;
import ch.systemsx.cisd.hdf5.HDF5DataTypeInformation;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5FloatStorageFeatures.HDF5FloatStorageFeatureBuilder;
import ch.systemsx.cisd.hdf5.HDF5GenericStorageFeatures.HDF5GenericStorageFeatureBuilder;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures.HDF5IntStorageFeatureBuilder;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import json.ErrorJSON;
import json.WarningJSON;
import model.MetaOn;
import model.Metadata;
import model.Metatype;
import model.Parameters;
import parsing.model.FileType;
import parsing.model.Gene;
import tools.Utils;

public class LoomFile 
{
	private static HashSet<LoomFile> all_open_handles = new HashSet<LoomFile>();
	
	private IHDF5Reader handle = null;
	private String loomPath = null;
	
	public boolean readOnly = false;

	public LoomFile(String type, String filename) 
	{
		loomPath = filename;
		if(type.equals("w")) // Always create new
		{
			createEmptyLoomFile(filename);
			readOnly = false;
		}
		else if(type.equals("r+")) // Modify if existing 
		{
			if(!new File(filename).exists()) createEmptyLoomFile(filename);

			// Handle the LOCK
			int nbLocked = 0;
			while(true)
			{
				boolean isLocked = false;
				
				try
				{				
					this.handle = HDF5Factory.open(filename);
					checkLoomFormat();
				}
				catch(Exception ex)
				{
					isLocked = true;
					nbLocked++;
					if(nbLocked >= 3600) // One hour
					{
						new ErrorJSON("The file cannot be unlocked");
					}
				}
				/*
				 Another exception is thrown randomly: UnsupportedOperationException
				 So I catch everything and pray :(
				catch(HDF5LibraryException ex)
				{
					if(ex.getMajorErrorNumber() == HDF5Constants.H5E_FILE && ex.getMinorErrorNumber() == HDF5Constants.H5E_BADFILE) isLocked = true;
					//if(ex.getMinorErrorNumber() == HDF5Constants.H5E_CANTLOCK) isLocked = true;				
					else new ErrorJSON(ex.getMessage());
				}*/
				if(!isLocked) break;
				else 
				{
					System.out.println("The file is locked.");
					try
					{
						Thread.sleep(1000);
					}
					catch(InterruptedException ie) { new ErrorJSON(ie.getMessage()); }
					Parameters.idleTime += 1;
				}
			}
			
			readOnly = false;
		}
		else if(type.equals("r")) 
		{
			if(!new File(filename).exists()) new ErrorJSON("This Loom file does not exists: " + filename);
			
			// Handle the LOCK
			int nbLocked = 0;
			while(true)
			{
				boolean isLocked = false;
				try
				{
					this.handle = HDF5Factory.openForReading(filename);
					checkLoomFormat();
				}
				catch(Exception ex)
				{
					isLocked = true;
					nbLocked++;
					if(nbLocked >= 3600) // One hour
					{
						new ErrorJSON("The file cannot be unlocked");
					}
				}
				/*
				 Another exception is thrown randomly: UnsupportedOperationException
				 So I catch everything and pray :(
				catch(HDF5LibraryException ex)
				{
					if(ex.getMajorErrorNumber() == HDF5Constants.H5E_FILE && ex.getMinorErrorNumber() == HDF5Constants.H5E_BADFILE) isLocked = true;
					//if(ex.getMinorErrorNumber() == HDF5Constants.H5E_CANTLOCK) isLocked = true;				
					else new ErrorJSON(ex.getMessage());
				}*/
				if(!isLocked) break;
				else 
				{
					System.out.println("The file is locked.");
					try
					{
						Thread.sleep(1000);
					}
					catch(InterruptedException ie) { new ErrorJSON(ie.getMessage()); }
					Parameters.idleTime += 1;
				}
			}
			
			readOnly = true;
		}
		else new ErrorJSON("Unknown type: " + type);
		all_open_handles.add(this);
	}
	
	public String getLoomPath()
	{
		return loomPath;
	}
	
	public void checkLoomFormat() // Following http://linnarssonlab.org/loompy/format/index.html#specification
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		StringBuffer sb = new StringBuffer();
		if(!this.handle.object().isDataSet("/matrix")) sb.append("Loom format is not correct, there is no /matrix dataset\n");
		if(!this.handle.isGroup("col_attrs")) sb.append("Loom format is not correct, there is no /col_attrs group");
		if(!this.handle.isGroup("col_graphs")) sb.append("Loom format is not correct, there is no /col_graphs group");
		if(!this.handle.isGroup("row_attrs")) sb.append("Loom format is not correct, there is no /row_attrs group");
		if(!this.handle.isGroup("row_graphs")) sb.append("Loom format is not correct, there is no /row_graphs group");
		String res = sb.toString();
		if(res.length() != 0) { this.close(); new ErrorJSON(res); }
	}
	
	public FileType isLoom() // Following http://linnarssonlab.org/loompy/format/index.html#specification
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(!this.handle.object().isDataSet("/matrix")) return FileType.UNKNOWN;
		if(!this.handle.isGroup("col_attrs")) return FileType.UNKNOWN;
		if(!this.handle.isGroup("col_graphs")) return FileType.UNKNOWN;
		if(!this.handle.isGroup("row_attrs")) return FileType.UNKNOWN;
		if(!this.handle.isGroup("row_graphs")) return FileType.UNKNOWN;
		Parameters.loomVersion = "2.0.0";
		if(!this.handle.isGroup("attrs")) return FileType.LOOM;
		if(!this.handle.exists("attrs/LOOM_SPEC_VERSION")) return FileType.LOOM;
		Parameters.loomVersion = this.handle.readString("attrs/LOOM_SPEC_VERSION");
		return FileType.LOOM;
	}
	
	public float[] readRow(long index, String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		long[] dim = this.handle.getDataSetInformation(path).getDimensions(); // Know the actual size of the array
		float[][] res = this.handle.float32().readMatrixBlock(path, 1, (int)dim[1], index, 0l); // TODO does not work if too big array
		return res[0];
	}
	
	public float[] readCol(long index, String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		long[] dim = this.handle.getDataSetInformation(path).getDimensions(); // Know the actual size of the array
		float[][] res = this.handle.float32().readMatrixBlock(path, (int)dim[0], 1, 0l, index); // TODO does not work if too big array
		return Utils.t(res)[0];
	}
	
	public float[][] readFloatBlock(String path, int blockSizeX, int blockSizeY, long nbBlocksX, long nbBlocksY)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		return this.handle.float32().readMatrixBlock(path, blockSizeX, blockSizeY, nbBlocksX, nbBlocksY); // TODO does not work if too big array
	}
	
	public int[][] readIntBlock(String path, int blockSizeX, int blockSizeY, long nbBlocksX, long nbBlocksY)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		return this.handle.int32().readMatrixBlock(path, blockSizeX, blockSizeY, nbBlocksX, nbBlocksY); // TODO does not work if too big array
	}
	
	public List<Metadata> listMetadata()
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		List<Metadata> res = new ArrayList<>();
		List<String> m = this.handle.getGroupMembers("/row_attrs");
		for(String mem:m)
		{
			if(!this.handle.isGroup("/row_attrs/" + mem))
			{
				Metadata meta = new Metadata();
				meta.path = "/row_attrs/" + mem;
				meta.on = MetaOn.GENE;
				meta.size = this.getSizeInBytes("/row_attrs/" + mem);
				long[] dims = this.getDimensions("/row_attrs/" + mem); // TODO HANDLE THE 1D AND 2D ARRAYS
				if(dims.length == 0) // Not an array => single value
				{
					meta.nbrow = 1;
					meta.nbcol = 1;
				}
				else if(dims.length == 1)
				{
					meta.nbrow = dims[0];
					meta.nbcol = 1;
				}
				else
				{
					meta.nbrow = dims[0];
					meta.nbcol = dims[1];
				}
				res.add(meta);
			}
		}
		m = this.handle.getGroupMembers("/col_attrs");
		for(String mem:m)
		{
			if(!this.handle.isGroup("/col_attrs/" + mem))
			{
				Metadata meta = new Metadata();
				meta.path = "/col_attrs/" + mem;
				meta.on = MetaOn.CELL;
				meta.size = this.getSizeInBytes("/col_attrs/" + mem);
				long[] dims = this.getDimensions("/col_attrs/" + mem);
				if(dims.length == 0) // Not an array => single value
				{
					meta.nbrow = 1;
					meta.nbcol = 1;
				}
				else if(dims.length == 1)
				{
					meta.nbrow = 1;
					meta.nbcol = dims[0];
				}
				else
				{
					meta.nbrow = dims[1];
					meta.nbcol = dims[0];
				}
				res.add(meta);
			}
		}
		m = this.handle.getGroupMembers("/attrs");
		for(String mem:m)
		{
			if(!this.handle.isGroup("/attrs/" + mem))
			{
				Metadata meta = new Metadata();
				meta.path = "/attrs/" + mem;
				meta.on = MetaOn.GLOBAL;
				meta.size = this.getSizeInBytes(meta.path);
				long[] dims = this.getDimensions(meta.path);
				
				if(dims.length == 0) // Not an array => single value
				{
					meta.nbrow = 1;
					meta.nbcol = 1;
				}
				else if(dims.length == 1)
				{
					meta.nbrow = 1;
					meta.nbcol = dims[0];
				}
				else // matrix
				{
					meta.nbrow = dims[1];
					meta.nbcol = dims[0];
				}
				res.add(meta);
			}
		}
		m = this.handle.getGroupMembers("/layers");
		for(String mem:m)
		{
			if(!this.handle.isGroup("/layers/" + mem))
			{
				Metadata meta = new Metadata();
				meta.path = "/layers/" + mem;
				meta.on = MetaOn.EXPRESSION_MATRIX;
				meta.size = this.getSizeInBytes(meta.path);
				long[] dims = this.getDimensions(meta.path);
				if(dims.length < 2) // Not an array => single value
				{
					new ErrorJSON("/layers group should contain only 2D matrices");
				}
				else
				{
					meta.nbrow = dims[0];
					meta.nbcol = dims[1];
				}
				res.add(meta);
			}
		}
		return res;
	}
	
	public int[] getChunkSizes()
	{
		return getChunkSizes("/matrix");
	}
	
	public int[] getChunkSizes(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		return this.handle.getDataSetInformation(path).tryGetChunkSizes();
	}
	
	public void removeMetadata(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) { this.close(); new ErrorJSON("Cannot delete from Read-Only file"); }
		if(this.handle.object().exists(path)) ((IHDF5Writer)this.handle).object().delete(path);
		else { this.close(); new ErrorJSON("This metadata (" + path + ") is not found in the Loom file"); }
	}

	public ArrayList<Long> getIndexesWhereValueIs(String path, String value)
	{
		// First check if loom file is opened
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		
		// Then check if the path corresponds to a metadata path
		if(!path.contains("col_attrs/") && !path.contains("row_attrs/")) new ErrorJSON("Your path is not a proper CELL or GENE metadata");
		
		// Finally check if the path exist in the Loom
		if(!this.handle.exists(path)) new ErrorJSON("Error in the Loom file. Path " + path + " does not exist!");
			
		// Check dimensions
		long[] dim = this.getDimensions(path);
		if(dim.length > 1 && dim[1] > 0) new ErrorJSON("Not implemented for matrices..."); // Matrix
		
		// Variables
		ArrayList<Long> indexesMatching = new ArrayList<Long>();
		
		// Check the metadata type and create the output
		try
		{
			HDF5DataClass h5_class = this.getDataClass(path);
			switch(h5_class)
			{
				case FLOAT:
					FloatArray64 fValues = readFloatArray(path);
					float fToCompare= Float.parseFloat(value.replaceAll(",", "."));
					for (long i = 0; i < fValues.size(); i++) if(fValues.get(i) == fToCompare) indexesMatching.add(i);
					break;		
				case INTEGER:
					IntArray64 iValues = readIntArray(path);
					int iToCompare= Integer.parseInt(value.replaceAll(",", "."));
					for (long i = 0; i < iValues.size(); i++) if(iValues.get(i) == iToCompare) indexesMatching.add(i);
					break;
				case STRING:
					StringArray64 sValues = readStringArray(path);
					for (long i = 0; i < sValues.size(); i++) if(sValues.get(i).equals(value)) indexesMatching.add(i);
					break;
				default:
					new ErrorJSON("Cannot handle this type of metadata yet.");
					break;
			}
		}
		catch(NumberFormatException nfe)
		{
			new ErrorJSON("The value you entered with -sel (" + value + ") is incompatible with the metadata you selected (" + path + ")");
		}
		
		// Return the indexes matching the value
		return indexesMatching;
	}
	
	public Metadata fillInfoMetadata(String path, boolean checkCategories)
	{
		Metadata out = new Metadata();
		out.path = path;
		out.size = this.getSizeInBytes(path);
		
		// Check if loom file is opened
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		// Check if the path exist in the Loom
		if(!this.handle.exists(path)) new ErrorJSON("Error in the Loom file. Path " + path + " does not exist!");
		
		// Then check if the path corresponds to a metadata path
		if(path.startsWith("/col_attrs/")) out.on = MetaOn.CELL;
		else if(path.startsWith("/row_attrs/")) out.on = MetaOn.GENE;
		else if(path.equals("/matrix") || path.startsWith("/layers")) out.on = MetaOn.EXPRESSION_MATRIX;
		else if(path.startsWith("/attrs/")) out.on = MetaOn.GLOBAL;
		else new ErrorJSON("Your path is not a proper CELL, GENE or GLOBAL metadata, nor a EXPRESSION_MATRIX");
		
		// For count occurrences of each value
		out.categories = new HashSet<>();
		out.categoriesMap = new HashMap<String, Long>();
		
		// Prepare reading	
		long[] dim = this.getDimensions(path);

		// Check the metadata type
		HDF5DataClass h5_class = this.getDataClass(path);
		switch(h5_class)
		{
			case FLOAT:
				out.type = Metatype.NUMERIC; // a priori
				if(dim.length == 0) // Not an array => single value
				{
					out.value = "" + readFloat(path);
					out.categories = null;
					out.nbcol = 1;
					out.nbrow = 1;
					out.type = Metatype.NUMERIC;
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					if(out.on != MetaOn.EXPRESSION_MATRIX)
					{
						if(checkCategories)
						{
							float[][] values = readFloatMatrix(path);
							out.matrixValues = new String[values.length][values[0].length];
							for (int i = 0; i < out.matrixValues.length; i++)
							{
								for (int j = 0; j < out.matrixValues[i].length; j++) 
								{
									float v = values[i][j];
									String vs = "" + v;
									out.categories.add(vs); // Checking if not discrete
									Long count = out.categoriesMap.get(vs);
									if(count == null) count = 0L;
									out.categoriesMap.put(vs, count + 1);
									out.matrixValues[i][j] =  vs;
								}
							}
							if(out.isCategorical(values.length * values[0].length)) out.type = Metatype.DISCRETE;
							else out.categories = null;
						}
					}
					else
					{
						out.nbcol = dim[1];
						out.nbrow = dim[0];
					}
				}
				else // Vector
				{
					if(checkCategories)
					{
						FloatArray64 values = readFloatArray(path);
						out.values = new StringArray64(values.size());
						for (long i = 0; i < out.values.size(); i++)
						{
							float v = values.get(i);
							String vs = "" + v;
							out.categories.add(vs); // Checking if not discrete
							Long count = out.categoriesMap.get(vs);
							if(count == null) count = 0L;
							out.categoriesMap.put(vs, count + 1);
							out.values.set(i, vs);
						}
						if(out.isCategorical()) out.type = Metatype.DISCRETE;
						else out.categories = null;
					}
				}
				break;		
			case INTEGER:
				out.type = Metatype.NUMERIC; // a priori
				if(dim.length == 0) // Not an array => single value
				{
					out.value = "" + readInt(path);
					out.categories = null;
					out.nbcol = 1;
					out.nbrow = 1;
					out.type = Metatype.NUMERIC;
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					if(out.on != MetaOn.EXPRESSION_MATRIX)
					{
						if(checkCategories)
						{
							int[][] values = readIntMatrix(path);
							out.matrixValues = new String[values.length][values[0].length];
							for (int i = 0; i < out.matrixValues.length; i++)
							{
								for (int j = 0; j < out.matrixValues[i].length; j++) 
								{
									int v = values[i][j];
									String vs = "" + v;
									out.categories.add(vs); // Checking if not discrete
									Long count = out.categoriesMap.get(vs);
									if(count == null) count = 0L;
									out.categoriesMap.put(vs, count + 1);
									out.matrixValues[i][j] =  vs;
								}
							}
							if(out.isCategorical(values.length * values[0].length)) out.type = Metatype.DISCRETE;
							else out.categories = null;
						}
					}
					else
					{
						out.nbcol = dim[1];
						out.nbrow = dim[0];
					}
				}
				else // Vector
				{
					if(checkCategories)
					{
						IntArray64 values = readIntArray(path);
						out.values = new StringArray64(values.size());
						for (long i = 0; i < out.values.size(); i++)
						{
							int v = values.get(i);
							String vs = "" + v;
							out.categories.add(vs); // Checking if not discrete
							Long count = out.categoriesMap.get(vs);
							if(count == null) count = 0L;
							out.categoriesMap.put(vs, count + 1);
							out.values.set(i, vs);
						}
						if(out.isCategorical()) out.type = Metatype.DISCRETE;
						else out.categories = null;
					}
				}
				break;	
			case STRING:
				out.type = Metatype.STRING; // a priori
				if(dim.length == 0) // Not an array => single value
				{
					out.value = readString(path);
					out.value = Utils.handleSpecialCharacters(out.value);
					out.stringWasCorrected = true;
					out.categories = null;
					out.nbcol = 1;
					out.nbrow = 1;
					out.type = Metatype.STRING;
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					new ErrorJSON("A 2D array of STRING is not allowed");
				}
				else // Vector
				{
					if(checkCategories)
					{
						out.values = readStringArray(path);
						for (long i = 0; i < out.values.size(); i++)
						{
							String v = out.values.get(i);
							v = Utils.handleSpecialCharacters(v);
							out.categories.add(v); // Still checking if not discrete
							Long count = out.categoriesMap.get(v);
							if(count == null) count = 0L;
							out.categoriesMap.put(v, count + 1);
						}
						out.stringWasCorrected = true;
						if(out.isCategorical()) out.type = Metatype.DISCRETE;
						else out.categories = null;
					}
				}
				break;
			default:
				out.type = Metatype.NOT_HANDLED;
				break;
		}
		return out;
	}
	
	public long[] getDimensions() // 0 = Nb Cells, 1 = Nb Genes
	{
		return this.getDimensions("/matrix");
	}
	
	public long[] getDimensions(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(!this.handle.object().exists(path)) { this.close(); new ErrorJSON("Dataset '"+path+"' does not exist"); }
		return this.handle.getDataSetInformation(path).getDimensions();
	}
	
	public long getSizeInBytes(String path) // TODO if the dataste is GZIP, I cannot retrieve the true value
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(!this.handle.object().exists(path)) { this.close(); new ErrorJSON("Dataset '"+path+"' does not exist"); }
		HDF5DataSetInformation info = this.handle.getDataSetInformation(path);
		long nb = info.getNumberOfElements();;
		HDF5DataTypeInformation type = info.getTypeInformation();
		if(nb == 0) nb = 1;
		long esize = type.getElementSize();
		if(type.isVariableLengthString()) esize = type.getElementSizeForPadding() * 2; // 8 * 2 = 16
		return nb * esize;
	}
	
	public float[][] readRows(long indexS, int nb)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		long[] dim = this.handle.getDataSetInformation("/matrix").getDimensions(); // Know the actual size of the array
		float[][] res = this.handle.float32().readMatrixBlock("/matrix", nb, (int)dim[1], indexS, 0l); // TODO does not work if too big array
		return res;
	}
	
	public float[][] readRows(ArrayList<Long> indexes, String dataset)
	{
		if(indexes == null || indexes.isEmpty()) return null;
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(!this.exists(dataset)) new ErrorJSON("This dataset does not exist in the Loom file: " + dataset);
		// First, let's sort the indexes and remove duplicates
		TreeSet<Long> unique = new TreeSet<Long>();
		unique.addAll(indexes);
		Iterator<Long> it = unique.iterator();
		// Some info
		long[] dim = this.getDimensions(dataset); // Know the actual size of the array
		int[] blockSize = this.getChunkSizes(dataset);
		int nbTotalBlocks = (int)Math.ceil((float)dim[0] / blockSize[0]);
		// Output
		float[][] res = new float[indexes.size()][(int)dim[1]]; // TODO does not work if too many cells
		// Now, I'll go through the blocks of the dataset, and check if it contains our indexes of interest
		long currentIndex = it.next();
		int currentRow = 0;
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{
			// Check if our index is within this block range
			float[][] block = null;
			while(currentIndex >= nbBlocks * blockSize[0] && currentIndex < (nbBlocks + 1) * blockSize[0])
			{
				// Retrieve the required block with all columns (because we read row by row)
				if(block == null) block = this.handle.float32().readMatrixBlock(dataset, blockSize[0], (int)dim[1], nbBlocks, 0l);
				
				// Assign the result matrix (we need to find the index value in the current submatrix)
				res[currentRow++] = block[(int)(currentIndex - nbBlocks * blockSize[0])];
				
				// Next index to recuperate
				if(it.hasNext()) currentIndex = it.next();
				else currentIndex = -1;
			}

			// Note: For the last row, block.length will be < blockSize[0]
		}
		
		return res;
	}
	
	public long[] getGeneIndexes(String[] names)
	{
		int found = 0;
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		HashMap<String, Long> tmp = new HashMap<>();
		for(String s:names) tmp.put(s.toUpperCase(), -1l);
		if(this.handle.object().isDataSet("/row_attrs/Gene"))
		{
			StringArray64 genes = this.readStringArray("/row_attrs/Gene");
			for(long i = 0; i < genes.size(); i++) 
			{
				String g = genes.get(i).toUpperCase();
				if(tmp.containsKey(g)) 
				{
					tmp.put(g, i);
					found++;
					if(found == names.length) break;
				}
			}
		}
		if(this.handle.object().isDataSet("/row_attrs/Accession") && found != names.length)
		{
			StringArray64 genes = this.readStringArray("/row_attrs/Accession");
			for(long i = 0; i < genes.size(); i++) 
			{
				String g = genes.get(i).toUpperCase();
				if(tmp.containsKey(g)) 
				{
					tmp.put(g, i);
					found++;
					if(found == names.length) break;
				}
			}
		}
		// We do not store the Alt Names in the Loom file
		
		// Prepare output res
		long[] res = new long[names.length];
		for (int i = 0; i < names.length; i++) res[i] = tmp.get(names[i].toUpperCase());
		
		return res;
	}
	
	public long[] getGeneIndexesByStableIds(long[] stable_ids)
	{
		int found = 0;
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		HashMap<Long, Long> tmp = new HashMap<>();
		for(long stable_id:stable_ids) tmp.put(stable_id, -1l);
		if(this.handle.object().isDataSet("/row_attrs/_StableID"))
		{
			LongArray64 genes = this.readLongArray("/row_attrs/_StableID");
			for(long i = 0; i < genes.size(); i++) 
			{
				long g = genes.get(i);
				if(tmp.containsKey(g)) 
				{
					tmp.put(g, i);
					found++;
					if(found == stable_ids.length) break;
				}
			}
		}
		// Prepare output res
		long[] res = new long[stable_ids.length];
		for (int i = 0; i < stable_ids.length; i++) res[i] = tmp.get(stable_ids[i]);
		
		return res;
	}
	
	public long[] getCellIndexes(String[] names)
	{
		int found = 0;
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		HashMap<String, Long> tmp = new HashMap<>();
		for(String s:names) tmp.put(s.toUpperCase(), -1l);
		if(this.handle.object().isDataSet("/col_attrs/CellID"))
		{
			StringArray64 cells = this.readStringArray("/col_attrs/CellID");
			for(long i = 0; i < cells.size(); i++) 
			{
				String c = cells.get(i).toUpperCase();
				if(tmp.containsKey(c)) 
				{
					tmp.put(c, i);
					found++;
					if(found == names.length) break;
				}
			}
		}
		// Prepare output res
		long[] res = new long[names.length];
		for (int i = 0; i < names.length; i++) res[i] = tmp.get(names[i].toUpperCase());
		
		return res;
	}
	
	public long[] getCellIndexesByStableIds(long[] stable_ids)
	{
		int found = 0;
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		HashMap<Long, Long> tmp = new HashMap<>();
		for(long stable_id:stable_ids) tmp.put(stable_id, -1l);
		if(this.handle.object().isDataSet("/col_attrs/_StableID"))
		{
			LongArray64 cells = this.readLongArray("/col_attrs/_StableID");
			for(long i = 0; i < cells.size(); i++) 
			{
				long c = cells.get(i);
				if(tmp.containsKey(c)) 
				{
					tmp.put(c, i);
					found++;
					if(found == stable_ids.length) break;
				}
			}
		} else new ErrorJSON("No StableID in this Loom file: " + this.loomPath);
		// Prepare output res
		long[] res = new long[stable_ids.length];
		for (int i = 0; i < stable_ids.length; i++) res[i] = tmp.get(stable_ids[i]);
		
		return res;
	}
	
	public LongArray64 getCellIndexesByStableIds(LongArray64 stable_ids)
	{
		int found = 0;
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		HashMap<Long, Long> tmp = new HashMap<>(); // I'm afraid this could be enormous?
		for(long sid = 0; sid < stable_ids.size(); sid ++) tmp.put(stable_ids.get(sid), -1l);
			
		if(this.handle.object().isDataSet("/col_attrs/_StableID"))
		{
			LongArray64 cells = this.readLongArray("/col_attrs/_StableID");
			for(long i = 0; i < cells.size(); i++) 
			{
				long c = cells.get(i);
				if(tmp.containsKey(c)) 
				{
					tmp.put(c, i);
					found++;
					if(found == stable_ids.size()) break;
				}
			}
		} else new ErrorJSON("No StableID in this Loom file: " + this.loomPath);
		// Prepare output res
		LongArray64 res = new LongArray64(stable_ids.size());
		for (int i = 0; i < stable_ids.size(); i++) res.set(i, tmp.get(stable_ids.get(i)));
		
		return res;
	}

	public float[] readCol(long index)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		long[] dim = this.handle.getDataSetInformation("/matrix").getDimensions(); // Know the actual size of the array
		float[][] res = this.handle.float32().readMatrixBlock("/matrix", (int)dim[0], 1, 0l, index); // TODO does not work if too big array
		float[] row = new float[(int)dim[0]];
		for(int i = 0; i < dim[0]; i++) row[i] = res[i][0];
		return row;
	}
	
	public LongArray64 readLongArray(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		long length = this.handle.getDataSetInformation(path).getDimensions()[0];
		LongArray64 res = new LongArray64(length); // Know the actual size of the array, and create it sufficiently big
		int nbChunks = (int)(length / LongArray64.chunkSize()) + 1;
		for(int i = 0; i < nbChunks; i++) res.set(i, this.handle.int64().readArrayBlock(path, LongArray64.chunkSize(), i));
		return res;
	}
	
	public DoubleArray64 readDoubleArray(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		long length = this.handle.getDataSetInformation(path).getDimensions()[0];
		DoubleArray64 res = new DoubleArray64(length); // Know the actual size of the array, and create it sufficiently big
		int nbChunks = (int)(length / DoubleArray64.chunkSize()) + 1;
		for(int i = 0; i < nbChunks; i++) res.set(i, this.handle.float64().readArrayBlock(path, DoubleArray64.chunkSize(), i));
		return res;
	}
	
	public FloatArray64 readFloatArray(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		long length = this.handle.getDataSetInformation(path).getDimensions()[0];
		FloatArray64 res = new FloatArray64(length); // Know the actual size of the array, and create it sufficiently big
		int nbChunks = (int)(length / FloatArray64.chunkSize()) + 1;
		for(int i = 0; i < nbChunks; i++) res.set(i, this.handle.float32().readArrayBlock(path, FloatArray64.chunkSize(), i));
		return res;
	}
	
	public float[][] readFloatMatrix(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		float[][] res = this.handle.float32().readMatrix(path); // TODO does not work if too big array
		return res;
	}
	
	public IntArray64 readIntArray(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		long length = this.handle.getDataSetInformation(path).getDimensions()[0];
		IntArray64 res = new IntArray64(length); // Know the actual size of the array, and create it sufficiently big
		int nbChunks = (int)(length / IntArray64.chunkSize()) + 1;
		for(int i = 0; i < nbChunks; i++) res.set(i, this.handle.int32().readArrayBlock(path, IntArray64.chunkSize(), i));
		return res;
	}
	
	public int[][] readIntMatrix(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		int[][] res = this.handle.int32().readMatrix(path); // TODO does not work if too big array
		return res;
	}
	
	public StringArray64 readStringArray(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(!this.handle.object().exists(path)) new ErrorJSON("Dataset '"+path+"' does not exist");
		long length = this.handle.getDataSetInformation(path).getDimensions()[0];
		StringArray64 res = new StringArray64(length); // Know the actual size of the array, and create it sufficiently big
		int nbChunks = (int)(length / StringArray64.chunkSize()) + 1;
		for(int i = 0; i < nbChunks; i++) res.set(i, this.handle.string().readArrayBlock(path, StringArray64.chunkSize(), i));
		return res;
	}
	
	public String readString(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(!this.handle.object().exists(path)) new ErrorJSON("Dataset '"+path+"' does not exist");
		return this.handle.string().read(path);
	}
	
	public float readFloat(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(!this.handle.object().exists(path)) new ErrorJSON("Dataset '"+path+"' does not exist");
		return this.handle.float32().read(path);
	}
	
	public int readInt(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(!this.handle.object().exists(path)) new ErrorJSON("Dataset '"+path+"' does not exist");
		return this.handle.int32().read(path);
	}
	
	public StringArray64 getCellNames()
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.exists("/col_attrs/CellID")) return this.readStringArray("/col_attrs/CellID");
		return null;
	}
	
	public LongArray64 getCellStableIds()
	{
		return this.readLongArray("/col_attrs/_StableID");
	}
	
	public LongArray64 getGeneStableIds()
	{
		return this.readLongArray("/row_attrs/_StableID");
	}
	
	public boolean exists(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		return this.handle.exists(path);
	}
	
	public HDF5DataClass getDataClass(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		return this.handle.getDataSetInformation(path).getTypeInformation().getDataClass();
	}
	
	public StringArray64 getGeneHGNC()
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.exists("/row_attrs/Gene")) return this.readStringArray("/row_attrs/Gene");
		return null;
	}
	
	public Gene[] getGeneNames()
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		Gene[] genes = null;
		if(this.exists("/row_attrs/Gene"))
		{
			String[] tmp = this.handle.string().readArray("/row_attrs/Gene");
			if(genes == null) genes = new Gene[tmp.length];
			for(int i = 0; i < tmp.length; i++) 
			{
				if(genes[i] == null) genes[i] = new Gene();
				genes[i].name= tmp[i];
			}
		}
		if(this.exists("/row_attrs/Accession"))
		{
			String[] tmp = this.handle.string().readArray("/row_attrs/Accession");
			if(genes == null) genes = new Gene[tmp.length];
			for(int i = 0; i < tmp.length; i++) 
			{
				if(genes[i] == null) genes[i] = new Gene();
				genes[i].ensembl_id = tmp[i];
			}
		}
		// We do not store the AltNames in the Loom file
		return genes;
	}
	
	public void close()
	{
		this.close(true);
	}
	
	public void close(boolean remove)
	{
		if(this.handle == null) System.err.println("Loom file is closed already");
		if(this.handle.object().exists("/__DATA_TYPES__") && !this.readOnly) ((IHDF5Writer)this.handle).object().delete("/__DATA_TYPES__");
		this.handle.close();
		if(remove) all_open_handles.remove(this);
		this.handle = null;
	}
	
	public static void close_all()
	{
		for(LoomFile f:all_open_handles) f.close(false);
		all_open_handles.clear();
	}
	
	/***** Writing *****/
	
	private void createEmptyLoomFile(String filename) // matrix.length == nb genes
    {
		if(handle != null) new ErrorJSON("Cannot open two Loom in the same Loom object, please close the first handle first");
		
    	// Create Loom file
    	IHDF5WriterConfigurator config = HDF5Factory.configure(filename);
    	this.handle = config.dontUseExtendableDataTypes().overwrite().useSimpleDataSpaceForAttributes().useUTF8CharacterEncoding().writer();
	    
        // Required annotations
    	((IHDF5Writer)this.handle).object().createGroup("/attrs"); // Loom v3
    	((IHDF5Writer)this.handle).string().write("/attrs/LOOM_SPEC_VERSION", "3.0.0"); // Loom v3
    	((IHDF5Writer)this.handle).object().createGroup("/col_attrs");
    	((IHDF5Writer)this.handle).object().createGroup("/row_attrs");
    	((IHDF5Writer)this.handle).object().createGroup("/layers"); // this is not required but we use it a lot
    	((IHDF5Writer)this.handle).object().createGroup("/row_graphs"); // Required group, even if empty
    	((IHDF5Writer)this.handle).object().createGroup("/col_graphs"); // Required group, even if empty
    }
	
	public void writeStringArray(String path, StringArray64 array)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly Loom file");
		int nbChunks = (int)(array.size() / StringArray64.chunkSize()) + 1;
		((IHDF5Writer)this.handle).string().createArrayVL(path, array.size(), (int)Math.min(StringArray64.chunkSize(), array.size()));
		for(int i = 0; i < nbChunks; i++) ((IHDF5Writer)this.handle).string().writeArrayBlock(path, array.getByChunk(i), i);
	}
	
	public void writeFloatArray(String path, FloatArray64 array)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly Loom file");
		int nbChunks = (int)(array.size() / FloatArray64.chunkSize()) + 1;
		((IHDF5Writer)this.handle).float32().createArray(path, array.size(), (int)Math.min(FloatArray64.chunkSize(), array.size()));
		for(int i = 0; i < nbChunks; i++) ((IHDF5Writer)this.handle).float32().writeArrayBlock(path, array.getByChunk(i), i);
	}
	
	public void writeDoubleArray(String path, DoubleArray64 array)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly Loom file");
		int nbChunks = (int)(array.size() / DoubleArray64.chunkSize()) + 1;
		((IHDF5Writer)this.handle).float64().createArray(path, array.size(), (int)Math.min(FloatArray64.chunkSize(), array.size()));
		for(int i = 0; i < nbChunks; i++) ((IHDF5Writer)this.handle).float64().writeArrayBlock(path, array.getByChunk(i), i);
	}
	
	public void writeIntArray(String path, IntArray64 array)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly Loom file");
		int nbChunks = (int)(array.size() / IntArray64.chunkSize()) + 1;
		((IHDF5Writer)this.handle).int32().createArray(path, array.size(), (int)Math.min(IntArray64.chunkSize(), array.size()));
		for(int i = 0; i < nbChunks; i++) ((IHDF5Writer)this.handle).int32().writeArrayBlock(path, array.getByChunk(i), i);
	}
	
	public void writeLongArray(String path, LongArray64 array)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly Loom file");
		int nbChunks = (int)(array.size() / LongArray64.chunkSize()) + 1;
		((IHDF5Writer)this.handle).int64().createArray(path, array.size(), (int)Math.min(LongArray64.chunkSize(), array.size()));
		for(int i = 0; i < nbChunks; i++) ((IHDF5Writer)this.handle).int64().writeArrayBlock(path, array.getByChunk(i), i);
	}
	
	public void writeFloat(String path, float value)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly Loom file");
		((IHDF5Writer)this.handle).float32().write(path, value);
	}
	
	public void writeInt(String path, int value)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly Loom file");
		((IHDF5Writer)this.handle).int32().write(path, value);
	}
	
	public void writeString(String path, String value)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly Loom file");
		((IHDF5Writer)this.handle).string().write(path, value);
	}
	
	public void writeMatrixMetadata(String path, float[][] data)
    {
		writeFloatMatrix(path, data);
    }
	
	public void writeFloatMatrix(String path, float[][] data)
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
		((IHDF5Writer)this.handle).float32().writeMatrix(path, data);
    }
	
	public void writeFloatArray(String path, float[] data)
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
		((IHDF5Writer)this.handle).float32().writeArray(path, data);
    }
	
	public void writeIntMatrix(String path, int[][] data)
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
		((IHDF5Writer)this.handle).int32().writeMatrix(path, data);
    }
	
	public void writeIntArray(String path, int[] data)
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
		((IHDF5Writer)this.handle).int32().writeArray(path, data);
    }
	
	public void writeStringArray(String path, String[] data)
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
		((IHDF5Writer)this.handle).string().writeArray(path, data);
    }
	
	/*public void writeBlockInFloatMatrix(float[][] data, int startIndex, HashSet<Long> filteredGenes) // matrix.length == nb genes
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        // Write all matrix stored in RAM, by chunks of [1, nbCells] i.e. gene by gene
	    if(writtenRow == -1) writtenRow = 0;
        for(int i = 0; i < data.length; i++)
	    {
	    	if(!filteredGenes.contains((long)startIndex + i))
	    	{
		    	float[][] m = new float[1][data[0].length];
		    	for(int k = 0; k < data[0].length; k++) m[0][k] = data[i][k];
		    	((IHDF5Writer)this.handle).float32().writeMatrixBlock("/matrix", m, writtenRow, 0);
		    	writtenRow++;
	    	}
	    }
    }
	
	public void writeBlockInMatrix(float[][] data, long nber_cells, int startIndex, HashSet<Long> filteredCells) // matrix.length == nb genes
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        // Write all matrix stored in RAM, by chunks of [1, nbCells] i.e. gene by gene
        for(int i = 0; i < data.length; i++)
	    {
		    float[][] m = new float[1][(int)nber_cells]; // TODO if big array
		    int pos = 0;
		    for(int k = 0; k < data[i].length; k++) 
		    {
		    	if(!filteredCells.contains((long)k))
		    	{
		    		m[0][pos] = data[i][k];
		    		pos++;
		    	}
		    }
		    ((IHDF5Writer)this.handle).float32().writeMatrixBlock("/matrix", m, startIndex + i, 0);
	    }
    }
	
	public void writeBlockInCustomChunkedMatrix(float[][] data, int startIndex) // matrix.length == nb genes
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        // Write all matrix stored in RAM, by chunks of [nbGenes, chunkLength] i.e. block of cells by block of cells
		((IHDF5Writer)this.handle).float32().writeMatrixBlock("/matrix", data, 0, startIndex);
    }*/
	
	public void writeFloatBlockDataset(String path, float[][] data, int blockNumberX, int blockNumberY) // matrix.length == nb genes
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
		((IHDF5Writer)this.handle).float32().writeMatrixBlock(path, data, blockNumberX, blockNumberY);
    }
	
	public void writeIntBlockDataset(String path, int[][] data, int blockNumberX, int blockNumberY) // matrix.length == nb genes
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
		((IHDF5Writer)this.handle).int32().writeMatrixBlock(path, data, blockNumberX, blockNumberY);
    }
	
	public void createEmptyFloat32MatrixDataset(String path, long nbGenes, long nbCells, int chunkSizeX, int chunkSizeY)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        HDF5FloatStorageFeatureBuilder features = new HDF5FloatStorageFeatureBuilder();
        features.noScaling().chunkedStorageLayout().datasetReplacementEnforceKeepExisting().deflateLevel((byte)2);
        ((IHDF5Writer)this.handle).float32().createMatrix(path, nbGenes, nbCells, chunkSizeX, chunkSizeY, features.features());
	}
	
	public void createEmptyInt32MatrixDataset(String path, long nbGenes, long nbCells, int chunkSizeX, int chunkSizeY)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        HDF5IntStorageFeatureBuilder features = new HDF5IntStorageFeatureBuilder();
        features.noScaling().chunkedStorageLayout().datasetReplacementEnforceKeepExisting().deflateLevel((byte)2);
        ((IHDF5Writer)this.handle).int32().createMatrix(path, nbGenes, nbCells, chunkSizeX, chunkSizeY, features.features());
	}
	
	public void createEmptyFloat32VectorDataset(String path, long length, int chunkSize)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        HDF5FloatStorageFeatureBuilder features = new HDF5FloatStorageFeatureBuilder();
        features.noScaling().chunkedStorageLayout().datasetReplacementEnforceKeepExisting().deflateLevel((byte)2);
        ((IHDF5Writer)this.handle).float32().createArray(path, length, chunkSize, features.features());
	}
	
	public void createEmptyInt32VectorDataset(String path, long length, int chunkSize)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        HDF5IntStorageFeatureBuilder features = new HDF5IntStorageFeatureBuilder();
        features.noScaling().chunkedStorageLayout().datasetReplacementEnforceKeepExisting().deflateLevel((byte)2);
        ((IHDF5Writer)this.handle).int32().createArray(path, length, chunkSize, features.features());
	}
	
	public void createEmptyStringVectorDataset(String path, int length, int chunkSize)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        HDF5GenericStorageFeatureBuilder features = new HDF5GenericStorageFeatureBuilder();
        features.noScaling().chunkedStorageLayout().datasetReplacementEnforceKeepExisting().deflateLevel((byte)2);
        ((IHDF5Writer)this.handle).string().createArray(path, length, chunkSize, features.features());
	}
			
	public void resizeDataset(String path, long x, long y)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly");
		((IHDF5Writer)this.handle).object().setDataSetDimensions(path, new long[] { x, y });
	}
    
    /**** Static methods ******/
	public static void copyFloatMatrixByCell(String path, LoomFile from, LoomFile to)
	{
    	LoomFile.copyFloatMatrixByCell(path, from, to, null, null);
	}
	
    public static void copyFloatMatrixByCell(String path, LoomFile from, LoomFile to, LoomData data)
	{
    	LoomFile.copyFloatMatrixByCell(path, from, to, null, data);
	}
    
    public static void copyFloatMatrixByCell(String path, LoomFile from, LoomFile to, HashSet<Long> toFilterGenes, LoomData data)
	{
    	if(toFilterGenes == null) toFilterGenes = new HashSet<Long>();

		long[] dim = from.getDimensions(); // dim[0] == genes
    	if(data != null && data.nber_genes != (dim[0] - toFilterGenes.size())) new ErrorJSON("data size does not match!");
    	if(data != null && data.nber_cells != dim[1]) new ErrorJSON("data size does not match!");	
		
    	int[] blockSize = from.getChunkSizes();
		blockSize[0] = 64; // We take all rows anyways, so we can change to whatever size we want
		if(blockSize[1] < 64 && 64 % blockSize[1] == 0) blockSize[1] = 64; // If chunk < 64 and multiple of 64, then better to pass to 64. If not, then we keep the original stuff
		int nbTotalBlocks = (int)Math.ceil((float)dim[1] / blockSize[1]);
    	
    	// Initialize main matrix in Loom with 0s
    	to.createEmptyFloat32MatrixDataset(path, dim[0] - toFilterGenes.size(), nbTotalBlocks * blockSize[1], blockSize[0], blockSize[1]); // I create a bit more cols, that will be filtered later on
    	
    	// Read the original file blockSize x totalNbGenes
		//System.out.println("Writing & Parsing " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = from.readFloatBlock(path, (int)dim[0], blockSize[1], 0l, nbBlocks); // TODO handle bigger datasets? Why blockSize is int? Am I doing the correct stuff?
			
			// Compute stats if needed (for parsing mainly)
			if(data != null)
			{
				// Parsing Data and generating summary annotations
				int i_filt_index = -1;
				for(int i = 0; i < subMatrix.length; i++) // Original index of gene
				{
					if(!toFilterGenes.contains((long)i)) // If not gene filtered
					{
						i_filt_index++;
						for(int y = 0; y < subMatrix[0].length; y++) // cells
						{
							long j = y + nbBlocks * blockSize[1]; // Original index of cells
							
			    			if(j < dim[1]) // In case the block is bigger than the number of cells
			    			{
								float value = subMatrix[i][y];
								
								if(data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) data.is_count_table = false;
									
				    			// Handle biotype count per cell
				    			String biotype = data.biotypes.get(i_filt_index);
				    			if(biotype.equals("protein_coding")) data.proteinCodingContent.set(j, data.proteinCodingContent.get(j) + value);
					    		else if(biotype.equals("rRNA")) data.ribosomalContent.set(j, data.ribosomalContent.get(j) + value);
				    			if(data.chromosomes.get(i_filt_index).equals("MT")) data.mitochondrialContent.set(j, data.mitochondrialContent.get(j) + value);
				        			
				        		// Generate the sums
				    			data.depth.set(j, data.depth.get(j) + value);
				    			data.sum.set(i_filt_index, data.sum.get(i_filt_index) + value);
									
								// Number of zeroes / detected genes 
								if(value == 0) data.nber_zeros++;
								else data.detected_genes.set(j, data.detected_genes.get(j) + 1);
							}
		    			}
					}
				}
			}
			
			// Restraining to non-filtered cells
			float[][] tmpMatrix = subMatrix;
			if(!toFilterGenes.isEmpty())
			{
				tmpMatrix = new float[(int)dim[0] - toFilterGenes.size()][blockSize[1]]; // TODO handle bigger datasets?
				int rowIndex = 0;
				for(int i = 0; i < subMatrix.length; i++)
				{
					if(!toFilterGenes.contains((long)i)) 
					{
						for(int j = 0; j < subMatrix[i].length; j++) tmpMatrix[rowIndex][j] = subMatrix[i][j];
						rowIndex++;
					}
				}
			}
			else if(subMatrix[0].length < blockSize[1]) // Adapting the size to block size if last block
			{
				tmpMatrix = new float[(int)dim[0]][blockSize[1]];
				for(int i = 0; i < subMatrix.length; i++) for(int j = 0; j < subMatrix[i].length; j++) tmpMatrix[i][j] = subMatrix[i][j];
			}
			
			// Writing this block to output
			to.writeFloatBlockDataset(path, tmpMatrix, 0, nbBlocks);
		}
		to.resizeDataset(path, (int)dim[0] - toFilterGenes.size(), (int)dim[1]); // Cause writing fixed-size blocks can extend the matrix size with 0
	}
    
    public static void copyFloatMatrixByGene(String path, LoomFile from, LoomFile to, LoomData data)
	{
    	LoomFile.copyFloatMatrixByGene(path, from, to, null, data);
	}
    
    public static void copyFloatMatrixByGene(String path, LoomFile from, LoomFile to, HashSet<Long> toFilterCells, LoomData data)
	{
    	if(toFilterCells == null) toFilterCells = new HashSet<Long>();

		long[] dim = from.getDimensions(); // dim[0] == genes
    	if(data != null && data.nber_genes != dim[0]) new ErrorJSON("data size does not match!");
    	if(data != null && data.nber_cells != (dim[1] - toFilterCells.size())) new ErrorJSON("data size does not match!");	
		
    	int[] blockSize = from.getChunkSizes();
		blockSize[1] = 64; // We take all rows anyways, so we can change to whatever size we want
		if(blockSize[0] < 64 && 64 % blockSize[0] == 0) blockSize[0] = 64; // If chunk < 64 and multiple of 64, then better to pass to 64. If not, then we keep the original stuff
		int nbTotalBlocks = (int)Math.ceil((float)dim[0] / blockSize[0]);
    	
    	// Initialize main matrix in Loom with 0s
    	to.createEmptyFloat32MatrixDataset(path, nbTotalBlocks * blockSize[0], dim[1] - toFilterCells.size(), blockSize[0], blockSize[1]); // I create a bit more cols, that will be filtered later on
    	
    	// Read the original file blockSize x totalNbGenes
		//System.out.println("Writing & Parsing " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			float[][] subMatrix = from.readFloatBlock(path, blockSize[0], (int)dim[1], nbBlocks, 0l); // TODO handle bigger datasets? Why blockSize is int? Am I doing the correct stuff?
			
			// First, generate some stats if required (for parsing mainly)
			if(data != null)
			{
				// Parsing Data and generating summary annotations
				for(int x = 0; x < subMatrix.length; x++)
				{
					long i = x + nbBlocks * blockSize[0]; // Original index of genes
	    			if(i < dim[0]) // In case the block is bigger than the number of genes
	    			{
	    				int j_filt_index = -1;
						for(int j = 0; j < subMatrix[0].length; j++) // cells
						{
							if(!toFilterCells.contains((long)j)) // If this cell is not filtered
							{
								j_filt_index++;
								
								float value = subMatrix[x][j];
								
								if(data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) data.is_count_table = false;
									
				    			// Handle biotype count per cell
				    			String biotype = data.biotypes.get(i);
				    			if(biotype.equals("protein_coding")) data.proteinCodingContent.set(j_filt_index, data.proteinCodingContent.get(j_filt_index) + value);
				    			else if(biotype.equals("rRNA")) data.ribosomalContent.set(j_filt_index, data.ribosomalContent.get(j_filt_index) + value);
				    			if(data.chromosomes.get(i).equals("MT")) data.mitochondrialContent.set(j_filt_index, data.mitochondrialContent.get(j_filt_index) + value);
				        			
				        		// Generate the sums
				    			data.depth.set(j_filt_index, data.depth.get(j_filt_index) + value);
				    			data.sum.set(i, data.sum.get(i) + value);
									
								// Number of zeroes / detected genes 
								if(value == 0) data.nber_zeros++;
								else data.detected_genes.set(j_filt_index, data.detected_genes.get(j_filt_index) + 1);
							}
						}
	    			}
				}
			}
			
			// Preparing to write data: restraining to non-filtered cells
			float[][] tmpMatrix = subMatrix;
			if(!toFilterCells.isEmpty())
			{
				tmpMatrix = new float[blockSize[0]][(int)dim[1] - toFilterCells.size()]; // TODO handle bigger datasets?
				int colIndex = 0;
				for(int j = 0; j < subMatrix[0].length; j++)
				{
					if(!toFilterCells.contains((long)j)) 
					{
						for(int i = 0; i < subMatrix.length; i++) tmpMatrix[i][colIndex] = subMatrix[i][j];
						colIndex++;
					}
				}
			}
			else if(subMatrix.length < blockSize[0]) // Adapting the size to block size if last block
			{
				tmpMatrix = new float[blockSize[0]][(int)dim[1]];
				for(int i = 0; i < subMatrix.length; i++) tmpMatrix[i] = subMatrix[i];
			}
		
			// Writing this block to output
			to.writeFloatBlockDataset(path, tmpMatrix, nbBlocks, 0);
		}
		to.resizeDataset(path, (int)dim[0], (int)dim[1] - toFilterCells.size()); // Cause writing fixed-size blocks can extend the matrix size with 0
	}
    
    public static void copyIntMatrixByCell(String path, LoomFile from, LoomFile to)
	{
    	LoomFile.copyIntMatrixByCell(path, from, to, null, null);
	}
    
    public static void copyIntMatrixByCell(String path, LoomFile from, LoomFile to, LoomData data)
	{
    	LoomFile.copyIntMatrixByCell(path, from, to, null, data);
	}
    
    public static void copyIntMatrixByCell(String path, LoomFile from, LoomFile to, HashSet<Long> toFilterGenes, LoomData data)
	{
    	if(toFilterGenes == null) toFilterGenes = new HashSet<Long>();

		long[] dim = from.getDimensions(); // dim[0] == genes
    	if(data != null && data.nber_genes != (dim[0] - toFilterGenes.size())) new ErrorJSON("data size does not match!");
    	if(data != null && data.nber_cells != dim[1]) new ErrorJSON("data size does not match!");	
		
    	int[] blockSize = from.getChunkSizes();
		blockSize[0] = 64; // We take all rows anyways, so we can change to whatever size we want
		if(blockSize[1] < 64 && 64 % blockSize[1] == 0) blockSize[1] = 64; // If chunk < 64 and multiple of 64, then better to pass to 64. If not, then we keep the original stuff
		int nbTotalBlocks = (int)Math.ceil((float)dim[1] / blockSize[1]);
    	
    	// Initialize main matrix in Loom with 0s
    	to.createEmptyInt32MatrixDataset(path, dim[0] - toFilterGenes.size(), nbTotalBlocks * blockSize[1], blockSize[0], blockSize[1]); // I create a bit more cols, that will be filtered later on
    	
    	// Read the original file blockSize x totalNbGenes
		//System.out.println("Writing & Parsing " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			int[][] subMatrix = from.readIntBlock(path, (int)dim[0], blockSize[1], 0l, nbBlocks); // TODO handle bigger datasets? Why blockSize is int? Am I doing the correct stuff?
			
			// Compute stats if needed (for parsing mainly)
			if(data != null)
			{
				int i_filt_index = -1;
				// Parsing Data and generating summary annotations
				for(int i = 0; i < subMatrix.length; i++) // Original index of gene
				{
					if(!toFilterGenes.contains((long)i)) // If not gene filtered
					{
						i_filt_index++;
						
						for(int y = 0; y < subMatrix[0].length; y++) // cells
						{
							long j = y + nbBlocks * blockSize[1]; // Original index of cells
							
			    			if(j < dim[1]) // In case the block is bigger than the number of cells
			    			{
								int value = subMatrix[i][y];
								
								if(data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) data.is_count_table = false;
									
				    			// Handle biotype count per cell
				    			String biotype = data.biotypes.get(i_filt_index);
				    			if(biotype.equals("protein_coding")) data.proteinCodingContent.set(j, data.proteinCodingContent.get(j) + value);
				    			else if(biotype.equals("rRNA")) data.ribosomalContent.set(j, data.ribosomalContent.get(j) + value);
				    			if(data.chromosomes.get(i_filt_index).equals("MT")) data.mitochondrialContent.set(j, data.mitochondrialContent.get(j) + value);
				        			
				        		// Generate the sums
				    			data.depth.set(j, data.depth.get(j) + value);
				    			data.sum.set(i_filt_index, data.sum.get(i_filt_index) + value);
									
								// Number of zeroes / detected genes 
								if(value == 0) data.nber_zeros++;
								else data.detected_genes.set(j, data.detected_genes.get(j) + 1);
							}
		    			}
					}
				}
			}
			
			// Restraining to non-filtered cells
			int[][] tmpMatrix = subMatrix;
			if(!toFilterGenes.isEmpty())
			{
				tmpMatrix = new int[(int)dim[0] - toFilterGenes.size()][blockSize[1]]; // TODO handle bigger datasets?
				int rowIndex = 0;
				for(int i = 0; i < subMatrix.length; i++)
				{
					if(!toFilterGenes.contains((long)i)) 
					{
						for(int j = 0; j < subMatrix[i].length; j++) tmpMatrix[rowIndex][j] = subMatrix[i][j];
						rowIndex++;
					}
				}
			}
			else if(subMatrix[0].length < blockSize[1]) // Adapting the size to block size if last block
			{
				tmpMatrix = new int[(int)dim[0]][blockSize[1]];
				for(int i = 0; i < subMatrix.length; i++) for(int j = 0; j < subMatrix[i].length; j++) tmpMatrix[i][j] = subMatrix[i][j];
			}
			
			// Writing this block to output
			to.writeIntBlockDataset(path, tmpMatrix, 0, nbBlocks);
		}
		to.resizeDataset(path, (int)dim[0] - toFilterGenes.size(), (int)dim[1]); // Cause writing fixed-size blocks can extend the matrix size with 0
	}
    
    public static void copyIntMatrixByGene(String path, LoomFile from, LoomFile to)
	{
    	LoomFile.copyIntMatrixByGene(path, from, to, null, null);
	}
    
    public static void copyIntMatrixByGene(String path, LoomFile from, LoomFile to, LoomData data)
	{
    	LoomFile.copyIntMatrixByGene(path, from, to, null, data);
	}
    
    public static void copyIntMatrixByGene(String path, LoomFile from, LoomFile to, HashSet<Long> toFilterCells, LoomData data)
	{
    	if(toFilterCells == null) toFilterCells = new HashSet<Long>();

		long[] dim = from.getDimensions(); // dim[0] == genes
    	if(data != null && data.nber_genes != dim[0]) new ErrorJSON("data size does not match!");
    	if(data != null && data.nber_cells != (dim[1] - toFilterCells.size())) new ErrorJSON("data size does not match!");
    	
    	int[] blockSize = from.getChunkSizes();
		blockSize[1] = 64; // We take all rows anyways, so we can change to whatever size we want
		if(blockSize[0] < 64 && 64 % blockSize[0] == 0) blockSize[0] = 64; // If chunk < 64 and multiple of 64, then better to pass to 64. If not, then we keep the original stuff
		int nbTotalBlocks = (int)Math.ceil((float)dim[0] / blockSize[0]);
    	
    	// Initialize main matrix in Loom with 0s
    	to.createEmptyInt32MatrixDataset(path, nbTotalBlocks * blockSize[0], dim[1] - toFilterCells.size(), blockSize[0], blockSize[1]); // I create a bit more cols, that will be filtered later on
    	
    	// Read the original file blockSize x totalNbGenes
		//System.out.println("Writing & Parsing " + nbTotalBlocks + " independent blocks...");
		for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
		{
			// Retrieve the blocks that will contain all columns (because we write gene by gene)
			int[][] subMatrix = from.readIntBlock(path, blockSize[0], (int)dim[1], nbBlocks, 0l); // TODO handle bigger datasets? Why blockSize is int? Am I doing the correct stuff?
			
			// First, generate some stats if required (for parsing mainly)
			if(data != null)
			{
				// Parsing Data and generating summary annotations
				for(int x = 0; x < subMatrix.length; x++)
				{
					long i = x + nbBlocks * blockSize[0]; // Original index of genes
	    			if(i < dim[0]) // In case the block is bigger than the number of genes
	    			{
	    				int j_filt_index = -1;
						for(int j = 0; j < subMatrix[0].length; j++) // cells
						{
							if(!toFilterCells.contains((long)j)) // If this cell is not filtered
							{
								j_filt_index++;
								
								int value = subMatrix[x][j];
								
								if(data.is_count_table && Math.abs(value - Math.round(value)) > 1E-5) data.is_count_table = false;
									
				    			// Handle biotype count per cell
				    			String biotype = data.biotypes.get(i);
				    			if(biotype.equals("protein_coding")) data.proteinCodingContent.set(j_filt_index, data.proteinCodingContent.get(j_filt_index) + value);
				    			else if(biotype.equals("rRNA")) data.ribosomalContent.set(j_filt_index, data.ribosomalContent.get(j_filt_index) + value);
				    			if(data.chromosomes.get(i).equals("MT")) data.mitochondrialContent.set(j_filt_index, data.mitochondrialContent.get(j_filt_index) + value);
				        			
				        		// Generate the sums
				    			data.depth.set(j_filt_index, data.depth.get(j_filt_index) + value);
				    			data.sum.set(i, data.sum.get(i) + value);
									
								// Number of zeroes / detected genes 
								if(value == 0) data.nber_zeros++;
								else data.detected_genes.set(j_filt_index, data.detected_genes.get(j_filt_index) + 1);
							}
						}
	    			}
				}
			}
			
			// Preparing to write data: restraining to non-filtered cells
			int[][] tmpMatrix = subMatrix;
			if(!toFilterCells.isEmpty())
			{
				tmpMatrix = new int[blockSize[0]][(int)dim[1] - toFilterCells.size()]; // TODO handle bigger datasets?
				int colIndex = 0;
				for(int j = 0; j < subMatrix[0].length; j++)
				{
					if(!toFilterCells.contains((long)j)) 
					{
						for(int i = 0; i < subMatrix.length; i++) tmpMatrix[i][colIndex] = subMatrix[i][j];
						colIndex++;
					}
				}
			}
			else if(subMatrix.length < blockSize[0]) // Adapting the size to block size if last block
			{
				tmpMatrix = new int[blockSize[0]][(int)dim[1]];
				for(int i = 0; i < subMatrix.length; i++) tmpMatrix[i] = subMatrix[i];
			}
		
			// Writing this block to output
			to.writeIntBlockDataset(path, tmpMatrix, nbBlocks, 0);
		}
		to.resizeDataset(path, (int)dim[0], (int)dim[1] - toFilterCells.size()); // Cause writing fixed-size blocks can extend the matrix size with 0
	}
    
    public static boolean copyMetadata(String path, LoomFile from, LoomFile to)
	{
    	Metadata meta = from.fillInfoMetadata(path, false);
    	return copyMetadata(meta, null, null, from, to);
	}
    
    public static boolean copyMetadata(Metadata m, LoomFile from, LoomFile to)
	{
    	return copyMetadata(m, null, null, from, to);
	}
    
    public static boolean copyMetadata(String path, HashSet<Long> toFilterGenes, HashSet<Long> toFilterCells, LoomFile from, LoomFile to)
   	{
       	Metadata meta = from.fillInfoMetadata(path, false);
       	return copyMetadata(meta, toFilterGenes, toFilterCells, from, to);
   	}
    
    public static boolean copyMetadata(Metadata m, HashSet<Long> toFilterGenes, HashSet<Long> toFilterCells, LoomFile from, LoomFile to)
	{
    	// Handling potential filters to perform
    	if(toFilterCells == null) toFilterCells = new HashSet<Long>();
    	if(toFilterGenes == null) toFilterGenes = new HashSet<Long>();
    	
		// Check if destination can be written
    	if(to.readOnly) new ErrorJSON("Cannot write metadata in 'to' Loom File. Please open it for writing.");
    	
    	// Get dataset infos
		long[] dim = from.getDimensions(m.path);
		int[] blockSize = from.getChunkSizes(m.path);
		boolean flag_copied = false;
		if(to.exists(m.path)) to.removeMetadata(m.path); // Remove if existing
		
		// Depending on the data type, we handle it differently
		HDF5DataClass h5_class = from.getDataClass(m.path);
    	switch (h5_class) 
		{
			case FLOAT:	
				if(dim.length == 0) // Not an array => single value
				{
					to.writeFloat(m.path, from.readFloat(m.path));
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					if(m.on == MetaOn.CELL)
					{
						float[][] values = from.readFloatMatrix(m.path);
						if(toFilterCells.isEmpty()) to.writeFloatMatrix(m.path, values);
						else
						{
							m.nbcol = dim[1];
							m.nbrow = values.length - toFilterCells.size();
							float[][] towrite = new float[(int)m.nbrow][(int)m.nbcol];
							int k = 0;
							for (int i = 0; i < values.length; i++) 
							{
								if(!toFilterCells.contains((long)i))
								{
									towrite[k] = values[i];
									k++;
								}
							}
							to.writeFloatMatrix(m.path, towrite);
						}
						m.size = to.getSizeInBytes(m.path);
					}
					else if(m.on == MetaOn.GENE)
					{
						float[][] values = from.readFloatMatrix(m.path);
						if(toFilterGenes.isEmpty()) to.writeFloatMatrix(m.path, values);
						else
						{
							m.nbcol = dim[1];
							m.nbrow = values.length - toFilterGenes.size();
							float[][] towrite = new float[(int)m.nbrow][(int)m.nbcol];
							int k = 0;
							for (int i = 0; i < values.length; i++) 
							{
								if(!toFilterGenes.contains((long)i))
								{
									towrite[k] = values[i];
									k++;
								}
							}
							to.writeFloatMatrix(m.path, towrite);
						}
						m.size = to.getSizeInBytes(m.path);
					}
					else if(m.on == MetaOn.EXPRESSION_MATRIX && blockSize == null) // blocksize == null => contiguous storage
					{
						float[][] values = from.readFloatMatrix(m.path);
						m.nbrow = values.length - toFilterGenes.size();
						m.nbcol = values[0].length - toFilterCells.size();
						float[][] towrite = new float[(int)m.nbrow][(int)m.nbcol];
						int k_i = 0;
						for (int i = 0; i < values.length; i++) 
						{
							if(!toFilterGenes.contains((long)i))
							{
								if(toFilterCells.isEmpty()) towrite[k_i] = values[i];
								else
								{
									int k_j = 0;
									for(int j = 0; j < values.length; j++) 
									{
										if(!toFilterCells.contains((long)j)) towrite[k_i][k_j++] = values[i][j];
									}
								}
								k_i++;
							}
						}
						to.writeFloatMatrix(m.path, towrite);
						m.size = to.getSizeInBytes(m.path);
					} 
					else if(m.on == MetaOn.EXPRESSION_MATRIX && blockSize != null) 
					{
						if(!toFilterCells.isEmpty() && !toFilterGenes.isEmpty()) new ErrorJSON("Filtering both GENEs AND CELLs in an EXPRESSION_MATRIX is not implemented yet");
						if(!toFilterCells.isEmpty()) copyFloatMatrixByGene(m.path, from, to, toFilterCells, null);
						else if(!toFilterGenes.isEmpty()) copyFloatMatrixByCell(m.path, from, to, toFilterGenes, null);
						else 
						{
					    	if(blockSize[0] / blockSize[1] > 2) copyFloatMatrixByCell(m.path, from, to, null); // Handling uneven chunks
					    	else LoomFile.copyFloatMatrixByGene(m.path, from, to, null);
						}
						m.size = to.getSizeInBytes(m.path);
					}
					else new ErrorJSON("Why is this Matrix dataset stored here??" + m.path);
				}
				else // Vector
				{
					FloatArray64 values = from.readFloatArray(m.path);
					m.nbrow = values.size() - toFilterCells.size();
					m.nbcol = values.size() - toFilterGenes.size(); // We don't really change the nb of cols. Just for ASAP.
					
					if(m.on == MetaOn.CELL)
					{
						if(toFilterCells.isEmpty()) to.writeFloatArray(m.path, values);
						else
						{
							FloatArray64 towrite = new FloatArray64(values.size() - toFilterCells.size());
							long k = 0;
							for (long i = 0; i < values.size(); i++) 
							{
								if(!toFilterCells.contains(i))
								{
									towrite.set(k, values.get(i));
									k++;
								}
							}
							to.writeFloatArray(m.path, towrite);
							m.size = to.getSizeInBytes(m.path);
						}
					}
										
					if(m.on == MetaOn.GENE)
					{
						if(toFilterGenes.isEmpty()) to.writeFloatArray(m.path, values);
						else
						{
							FloatArray64 towrite = new FloatArray64(values.size() - toFilterGenes.size());
							long k = 0;
							for (long i = 0; i < values.size(); i++) 
							{
								if(!toFilterGenes.contains(i))
								{
									towrite.set(k, values.get(i));
									k++;
								}
							}
							to.writeFloatArray(m.path, towrite);
							m.size = to.getSizeInBytes(m.path);
						}
					}
				}
				flag_copied = true;
				break;
			case INTEGER:
				if(dim.length == 0) // Not an array => single value
				{
					to.writeInt(m.path, from.readInt(m.path));
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					if(m.on == MetaOn.CELL)
					{
						int[][] values = from.readIntMatrix(m.path);
						if(toFilterCells.isEmpty()) to.writeIntMatrix(m.path, values);
						else
						{
							m.nbcol = dim[1];
							m.nbrow = values.length - toFilterCells.size();
							int[][] towrite = new int[(int)m.nbrow][(int)m.nbcol];
							int k = 0;
							for (int i = 0; i < values.length; i++) 
							{
								if(!toFilterCells.contains((long)i))
								{
									towrite[k] = values[i];
									k++;
								}
							}
							to.writeIntMatrix(m.path, towrite);
						}
						m.size = to.getSizeInBytes(m.path);
					}
					else if(m.on == MetaOn.GENE)
					{
						int[][] values = from.readIntMatrix(m.path);
						if(toFilterGenes.isEmpty()) to.writeIntMatrix(m.path, values);
						else
						{
							m.nbcol = dim[1];
							m.nbrow = values.length - toFilterGenes.size();
							int[][] towrite = new int[(int)m.nbrow][(int)m.nbcol];
							int k = 0;
							for (int i = 0; i < values.length; i++) 
							{
								if(!toFilterGenes.contains((long)i))
								{
									towrite[k] = values[i];
									k++;
								}
							}
							to.writeIntMatrix(m.path, towrite);
						}
						m.size = to.getSizeInBytes(m.path);
					}
					else if(m.on == MetaOn.EXPRESSION_MATRIX && blockSize == null) // blocksize == null => contiguous storage
					{
						int[][] values = from.readIntMatrix(m.path);
						m.nbrow = values.length - toFilterGenes.size();
						m.nbcol = values[0].length - toFilterCells.size();
						int[][] towrite = new int[(int)m.nbrow][(int)m.nbcol];
						int k_i = 0;
						for (int i = 0; i < values.length; i++) 
						{
							if(!toFilterGenes.contains((long)i))
							{
								if(toFilterCells.isEmpty()) towrite[k_i] = values[i];
								else
								{
									int k_j = 0;
									for(int j = 0; j < values.length; j++) 
									{
										if(!toFilterCells.contains((long)j)) towrite[k_i][k_j++] = values[i][j];
									}
								}
								k_i++;
							}
						}
						to.writeIntMatrix(m.path, towrite);
						m.size = to.getSizeInBytes(m.path);
					} 
					else if(m.on == MetaOn.EXPRESSION_MATRIX && blockSize != null) 
					{
						if(!toFilterCells.isEmpty() && !toFilterGenes.isEmpty()) new ErrorJSON("Filtering both GENEs AND CELLs in an EXPRESSION_MATRIX is not implemented yet");
						if(!toFilterCells.isEmpty()) copyIntMatrixByGene(m.path, from, to, toFilterCells, null);
						else if(!toFilterGenes.isEmpty()) copyIntMatrixByCell(m.path, from, to, toFilterGenes, null);
						else 
						{
					    	if(blockSize[0] / blockSize[1] > 2) copyIntMatrixByCell(m.path, from, to, null); // Handling uneven chunks
					    	else LoomFile.copyIntMatrixByGene(m.path, from, to, null);
						}
						m.size = to.getSizeInBytes(m.path);
					}
					else new ErrorJSON("Why is this Matrix dataset stored here??" + m.path);
				}
				else // Vector
				{
					IntArray64 values = from.readIntArray(m.path);
					m.nbrow = values.size() - toFilterCells.size();
					m.nbcol = values.size() - toFilterGenes.size(); // We don't really change the nb of cols. Just for ASAP.
					
					if(m.on == MetaOn.CELL)
					{
						if(toFilterCells.isEmpty()) to.writeIntArray(m.path, values);
						else
						{
							IntArray64 towrite = new IntArray64(values.size() - toFilterCells.size());
							long k = 0;
							for (long i = 0; i < values.size(); i++) 
							{
								if(!toFilterCells.contains(i))
								{
									towrite.set(k, values.get(i));
									k++;
								}
							}
							to.writeIntArray(m.path, towrite);
							m.size = to.getSizeInBytes(m.path);
						}
					}
										
					if(m.on == MetaOn.GENE)
					{
						if(toFilterGenes.isEmpty()) to.writeIntArray(m.path, values);
						else
						{
							IntArray64 towrite = new IntArray64(values.size() - toFilterGenes.size());
							long k = 0;
							for (long i = 0; i < values.size(); i++) 
							{
								if(!toFilterGenes.contains(i))
								{
									towrite.set(k, values.get(i));
									k++;
								}
							}
							to.writeIntArray(m.path, towrite);
							m.size = to.getSizeInBytes(m.path);
						}
					}
				}
				flag_copied = true;
				break;
			case STRING:
				if(dim.length == 0) // Not an array => single value
				{
					to.writeString(m.path, from.readString(m.path));
				}
				else if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					new ErrorJSON("Cannot create a 2D array of STRING");
				}
				else // Vector
				{
					StringArray64 values = from.readStringArray(m.path);
					m.nbrow = values.size() - toFilterCells.size();
					m.nbcol = values.size() - toFilterGenes.size(); // We don't really change the nb of cols. Just for ASAP.
					
					if(m.on == MetaOn.CELL)
					{
						if(toFilterCells.isEmpty()) to.writeStringArray(m.path, values);
						else
						{
							StringArray64 towrite = new StringArray64(values.size() - toFilterCells.size());
							long k = 0;
							for (long i = 0; i < values.size(); i++) 
							{
								if(!toFilterCells.contains(i))
								{
									towrite.set(k, values.get(i));
									k++;
								}
							}
							to.writeStringArray(m.path, towrite);
							m.size = to.getSizeInBytes(m.path);
						}
					}
										
					if(m.on == MetaOn.GENE)
					{
						if(toFilterGenes.isEmpty()) to.writeStringArray(m.path, values);
						else
						{
							StringArray64 towrite = new StringArray64(values.size() - toFilterGenes.size());
							long k = 0;
							for (long i = 0; i < values.size(); i++) 
							{
								if(!toFilterGenes.contains(i))
								{
									towrite.set(k, values.get(i));
									k++;
								}
							}
							to.writeStringArray(m.path, towrite);
							m.size = to.getSizeInBytes(m.path);
						}
					}
				}
				flag_copied = true;
				break;
			default:
				m.type = Metatype.NOT_HANDLED;
				break;
		}
		return flag_copied;
	}
    
    public void writeMetadata(String path, float[] array)
	{  	
		// Check if destination can be written
    	if(this.readOnly) new ErrorJSON("Cannot write metadata in Loom File. Please open it for writing.");
    	
    	// Get dataset infos
		if(this.exists(path)) 
		{
			WarningJSON.addWarning(path + " was already in the LOOM, it was overwriten");
			this.removeMetadata(path); // Remove if existing
		}
		
		// Write it
		this.writeFloatArray(path, array);
	}
    
    public boolean writeMetadata(Metadata m)
	{
    	boolean flag_copied = false;
    	
		// Check if destination can be written
    	if(this.readOnly) new ErrorJSON("Cannot write metadata in Loom File. Please open it for writing.");
    	
    	// Get dataset infos
		if(this.exists(m.path)) 
		{
			WarningJSON.addWarning(m.path + " was already in the LOOM, it was overwriten");
			this.removeMetadata(m.path); // Remove if existing
		}
		
		// Depending on the data type, we handle it differently
    	switch (m.type) 
		{
    		case NUMERIC:
    			FloatArray64 array = new FloatArray64(m.values.size());
    			boolean isInteger = true;
    			long missing = 0;
    			for(long i = 0; i < array.size(); i++)
    			{
    				String tmp = m.values.get(i);
    				if(tmp.equals("")) 
    				{
    					missing++;
    					tmp = Parameters.defaultMissingValue; // In case there is a missing value
    				}
    				float value = Float.parseFloat(tmp.replaceAll(",", "."));
    				if(isInteger && (int)value != value) isInteger = false; 
    				array.set(i, value);
    			}
    			if(missing != 0) WarningJSON.addWarning(missing + " missing values replaced by " + Parameters.defaultMissingValue + " for metadata " + m.path);
    			if(isInteger) this.writeIntArray(m.path, array.toIntArray());
    			else this.writeFloatArray(m.path, array);
    			flag_copied = true;
    			break;
    		case DISCRETE:
    		case STRING:
    			this.writeStringArray(m.path, m.values);
    			flag_copied = true;
    			break;
    		default:
    			WarningJSON.addWarning(m.path + " was NOT written in the LOOM, this type is NOT HANDLED");
		}
    	
		return flag_copied;
	}
    
    public static void fillLoomFile(LoomFile loom, LoomData data) // matrix.length == nb genes
    {
    	if(loom.readOnly) new ErrorJSON("Cannot write in this Loom file.");
    	
    	// Create Metadata array
    	data.meta = new ArrayList<>();
    	data.meta.add(new Metadata("/attrs/LOOM_SPEC_VERSION", Metatype.STRING, MetaOn.GLOBAL, 1, 1)); // Loom v3
    	
        // Cell annotations
    	loom.writeStringArray("/col_attrs/CellID", data.cell_names);
    	data.meta.add(new Metadata("/col_attrs/CellID", Metatype.STRING, MetaOn.CELL, data.cell_names.size(), 1));
    	
    	// Stable Ids
    	LongArray64 stable_cell_ids = new LongArray64(data.cell_names.size());
    	for(long i = 0; i < stable_cell_ids.size(); i++) stable_cell_ids.set(i, i);
    	loom.writeLongArray("/col_attrs/_StableID", stable_cell_ids);
    	data.meta.add(new Metadata("/col_attrs/_StableID", Metatype.NUMERIC, MetaOn.CELL, stable_cell_ids.size(), 1));
    	
        // Other generated stuff
        loom.writeDoubleArray("/col_attrs/_Depth", data.depth);
        data.meta.add(new Metadata("/col_attrs/_Depth", Metatype.NUMERIC, MetaOn.CELL, data.depth.size(), 1));
        loom.writeIntArray("/col_attrs/_Detected_Genes", data.detected_genes);
        data.meta.add(new Metadata("/col_attrs/_Detected_Genes", Metatype.NUMERIC, MetaOn.CELL, data.detected_genes.size(), 1));
        
        // Biotypes/chromosomes
        for(long j = 0; j < data.mitochondrialContent.size(); j++) 
        {
        	if(data.depth.get(j) == 0) data.mitochondrialContent.set(j, 0f);
        	else data.mitochondrialContent.set(j, (float)(data.mitochondrialContent.get(j) / data.depth.get(j) * 100));
        }
        loom.writeFloatArray("/col_attrs/_Mitochondrial_Content", data.mitochondrialContent);
 		data.meta.add(new Metadata("/col_attrs/_Mitochondrial_Content", Metatype.NUMERIC, MetaOn.CELL, data.mitochondrialContent.size(), 1));
        
        for(long j = 0; j < data.proteinCodingContent.size(); j++)        
        {
        	if(data.depth.get(j) == 0)  data.proteinCodingContent.set(j, 0f);
        	else data.proteinCodingContent.set(j, (float)(data.proteinCodingContent.get(j) / data.depth.get(j) * 100));
        } 
        loom.writeFloatArray("/col_attrs/_Protein_Coding_Content", data.proteinCodingContent);
 		data.meta.add(new Metadata("/col_attrs/_Protein_Coding_Content", Metatype.NUMERIC, MetaOn.CELL, data.proteinCodingContent.size(), 1));
 		
        for(long j = 0; j < data.ribosomalContent.size(); j++) 
        {
        	if(data.depth.get(j) == 0)  data.ribosomalContent.set(j, 0f);
        	else data.ribosomalContent.set(j, (float)(data.ribosomalContent.get(j) / data.depth.get(j) * 100));
        }
        loom.writeFloatArray("/col_attrs/_Ribosomal_Content", data.ribosomalContent);
 		data.meta.add(new Metadata("/col_attrs/_Ribosomal_Content", Metatype.NUMERIC, MetaOn.CELL, data.ribosomalContent.size(), 1));
     	
        // Genes and Gene annotations     
 		loom.writeStringArray("/row_attrs/Original_Gene", data.original_gene_names);
    	data.meta.add(new Metadata("/row_attrs/Original_Gene", Metatype.STRING, MetaOn.GENE, 1, data.original_gene_names.size()));
 		loom.writeStringArray("/row_attrs/Gene", data.gene_names); // HGNC names
    	data.meta.add(new Metadata("/row_attrs/Gene", Metatype.STRING, MetaOn.GENE, 1, data.gene_names.size()));
    	loom.writeStringArray("/row_attrs/Accession", data.ens_names);  // EnsemblIDs
    	data.meta.add(new Metadata("/row_attrs/Accession", Metatype.STRING, MetaOn.GENE, 1, data.ens_names.size()));
    	// Stable Ids
    	LongArray64 stable_ids = new LongArray64(data.ens_names.size());
    	for(long i = 0; i < stable_ids.size(); i++) stable_ids.set(i, i);
    	loom.writeLongArray("/row_attrs/_StableID", stable_ids);
    	data.meta.add(new Metadata("/row_attrs/_StableID", Metatype.NUMERIC, MetaOn.GENE, 1, stable_ids.size()));
    	// We do not store alt names in the Loom file 
    	loom.writeDoubleArray("/row_attrs/_Sum", data.sum);
    	data.meta.add(new Metadata("/row_attrs/_Sum", Metatype.NUMERIC, MetaOn.GENE, 1, data.sum.size()));
    	loom.writeLongArray("/row_attrs/_SumExonLength", data.sumExonLength);
    	data.meta.add(new Metadata("/row_attrs/_SumExonLength", Metatype.NUMERIC, MetaOn.GENE, 1, data.sumExonLength.size()));
    	loom.writeStringArray("/row_attrs/_Biotypes", data.biotypes);
    	Metadata m = new Metadata("/row_attrs/_Biotypes", Metatype.DISCRETE, MetaOn.GENE, 1, data.biotypes.size(), data.biotypes);
    	m.fillMap();
    	m.values = null;
    	data.meta.add(m);
    	loom.writeStringArray("/row_attrs/_Chromosomes", data.chromosomes);
    	m = new Metadata("/row_attrs/_Chromosomes", Metatype.DISCRETE, MetaOn.GENE, 1, data.chromosomes.size(), data.chromosomes);
    	m.fillMap();
    	m.values = null;
    	data.meta.add(m);
    }
}
