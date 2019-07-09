package hdf5.loom;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

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
import hdf.hdf5lib.HDF5Constants;
import hdf.hdf5lib.exceptions.HDF5LibraryException;
import json.ErrorJSON;
import model.MetaOn;
import model.Metadata;
import model.Metatype;
import model.Parameters;
import parsing.model.Gene;
import tools.Utils;

public class LoomFile 
{
	private IHDF5Reader handle = null;
	
	public boolean readOnly = false;
	
	public int writtenRow; // For writing matrix

	public LoomFile(String type, String filename) 
	{
		writtenRow = -1;
		
		if(type.equals("w")) // Always create new
		{
			createEmptyLoomFile(filename);
			readOnly = false;
		}
		else if(type.equals("w+")) // Modify if existing 
		{
			if(!new File(filename).exists()) createEmptyLoomFile(filename);

			// Handle the LOCK
			while(true)
			{
				boolean isLocked = false;
				try
				{				
					this.handle = HDF5Factory.open(filename);
					if(!isLoomFormatOK()) new ErrorJSON("Loom format is not OK: " + filename);
				}
				catch(HDF5LibraryException ex)
				{
					if(ex.getMajorErrorNumber() == HDF5Constants.H5E_FILE && ex.getMinorErrorNumber() == HDF5Constants.H5E_BADFILE) isLocked = true;
					//if(ex.getMinorErrorNumber() == HDF5Constants.H5E_CANTLOCK) isLocked = true;				
					else new ErrorJSON(ex.getMessage());
				}
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
			// Handle the LOCK
			while(true)
			{
				boolean isLocked = false;
				try
				{
					this.handle = HDF5Factory.openForReading(filename);
					if(!isLoomFormatOK()) new ErrorJSON("Loom format is not OK: " + filename);
				}
				catch(HDF5LibraryException ex)
				{
					if(ex.getMajorErrorNumber() == HDF5Constants.H5E_FILE && ex.getMinorErrorNumber() == HDF5Constants.H5E_BADFILE) isLocked = true;
					//if(ex.getMinorErrorNumber() == HDF5Constants.H5E_CANTLOCK) isLocked = true;				
					else new ErrorJSON(ex.getMessage());
				}
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
		
	}
	
	public boolean isLoomFormatOK() // Following http://linnarssonlab.org/loompy/format/index.html#specification
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(!this.handle.object().isDataSet("/matrix")) return false;
		if(!this.handle.isGroup("col_attrs")) return false;
		if(!this.handle.isGroup("col_graphs")) return false;
		if(!this.handle.isGroup("row_attrs")) return false;
		if(!this.handle.isGroup("row_graphs")) return false;
		return true;
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
	
	public float[][] readBlock(int blockSizeX, int blockSizeY, long nbBlocks)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		return this.handle.float32().readMatrixBlock("/matrix", blockSizeX, blockSizeY, nbBlocks, 0l); // TODO does not work if too big array
	}
	
	public float[][] readFloatBlock(String path, int blockSizeX, int blockSizeY, long nbBlocks)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		return this.handle.float32().readMatrixBlock(path, blockSizeX, blockSizeY, nbBlocks, 0l); // TODO does not work if too big array
	}
	
	public int[][] readIntBlock(String path, int blockSizeX, int blockSizeY, long nbBlocks)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		return this.handle.int32().readMatrixBlock(path, blockSizeX, blockSizeY, nbBlocks, 0l); // TODO does not work if too big array
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
				if(dims.length == 1)
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
				long[] dims = this.getDimensions("/col_attrs/" + mem); // TODO HANDLE THE 1D AND 2D ARRAYS
				if(dims.length == 1)
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
		if(this.readOnly) new ErrorJSON("Cannot delete from Read-Only file");
		if(this.handle.object().exists(path)) ((IHDF5Writer)this.handle).object().delete(path);
		else new ErrorJSON("This metadata (" + path + ") is not found in the Loom file");
	}
	
	public HDF5DataClass getDatasetType(String path)
	{
		return this.handle.getDataSetInformation(path).getTypeInformation().getDataClass();
	}
	
	public Metadata readMetadata(String path)
	{
		Metadata out = new Metadata();
		out.path = path;
		out.size = this.getSizeInBytes(path);
		
		// First check if loom file is opened
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		
		// Then check if the path corresponds to a metadata path
		if(path.contains("col_attrs/")) out.on = MetaOn.CELL;
		else if(path.contains("row_attrs/")) out.on = MetaOn.GENE;
		else new ErrorJSON("Your path is not a proper CELL or GENE metadata");
		
		// Finally check if the path exist in the Loom
		if(!this.handle.exists(path)) new ErrorJSON("Error in the Loom file. Path " + path + " does not exist!");
		
		// For count occurrences of each value
		out.categories = new HashSet<>();
		out.categoriesMap = new HashMap<String, Long>();
		
		// Prepare reading	
		long[] dim = this.getDimensions(path);
		//int[] blockSize = this.getChunkSizes(path);
		//int nbTotalBlocks = (int)Math.ceil((double)dim[0] / blockSize[0]); // dim[0] === Nb genes ?
		//System.out.println("Reading " + nbTotalBlocks + " independent blocks...");

		// Check the metadata type
		HDF5DataClass h5_class = this.getDataClass(path);
		switch(h5_class)
		{
			case FLOAT:
				out.type = Metatype.NUMERIC; // a priori
				if(dim.length > 1 && dim[1] > 0) // Matrix
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
				else // Vector
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
				break;		
			case INTEGER:
				out.type = Metatype.NUMERIC; // a priori
				if(dim.length > 1 && dim[1] > 0) // Matrix
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
				else // Vector
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
				break;	
			case STRING:
				out.type = Metatype.STRING; // a priori
				if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					new ErrorJSON("A 2D array of STRING is not allowed");
				}
				else // Vector
				{
					out.values = readStringArray(path);
					for (long i = 0; i < out.values.size(); i++)
					{
						String v = out.values.get(i);
						out.categories.add(v); // Still checking if not discrete
						Long count = out.categoriesMap.get(v);
						if(count == null) count = 0L;
						out.categoriesMap.put(v, count + 1);
					}
					if(out.isCategorical()) out.type = Metatype.DISCRETE;
					else out.categories = null;
				}
				break;
			default:
				out.type = Metatype.NOT_HANDLED;
				break;
		}
		return out;
	}
	
	public long getNbERCCs()
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.handle.object().exists("/col_attrs/_ERCCs"))
		{
			long[] dim = this.handle.getDataSetInformation("/col_attrs/_ERCCs").getDimensions();
			if(dim.length > 1) return dim[1];
			return 1;
		}
		return 0;
	}
	
	public long[] getDimensions() // 0 = Nb Cells, 1 = Nb Genes
	{
		return this.getDimensions("/matrix");
	}
	
	public long[] getDimensions(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		return this.handle.getDataSetInformation(path).getDimensions();
	}
	
	public long getSizeInBytes(String path)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		//return this.handle.getDataSetInformation(path).getSize();
		HDF5DataSetInformation info = this.handle.getDataSetInformation(path);
		long nb = info.getNumberOfElements();
		HDF5DataTypeInformation type = info.getTypeInformation();
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
	
	public int getGeneIndex(String geneName)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.handle.object().isDataSet("/row_attrs/Gene"))
		{
			String[] genes = this.handle.string().readArray("/row_attrs/Gene");
			for(int i = 0; i < genes.length; i++) if(genes[i].equals(geneName)) return i; // CASE should match ?
		}
		if(this.handle.object().isDataSet("/row_attrs/Accession"))
		{
			String[] genes = this.handle.string().readArray("/row_attrs/Accession");
			for(int i = 0; i < genes.length; i++) if(genes[i].equals(geneName)) return i; // CASE should match ?
		}
		// We do not store the Alt Names in the Loom file
		return -1;
	}
	
	public int getCellIndex(String cellName)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.handle.object().isDataSet("/col_attrs/CellID"))
		{
			String[] cells = this.handle.string().readArray("/col_attrs/CellID");
			for(int i = 0; i < cells.length; i++) if(cells[i].equals(cellName)) return i;
		}
		return -1; // Not found
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
		long length = this.handle.getDataSetInformation(path).getDimensions()[0];
		StringArray64 res = new StringArray64(length); // Know the actual size of the array, and create it sufficiently big
		int nbChunks = (int)(length / StringArray64.chunkSize()) + 1;
		for(int i = 0; i < nbChunks; i++) res.set(i, this.handle.string().readArrayBlock(path, StringArray64.chunkSize(), i));
		return res;
	}
	
	public StringArray64 getCellNames()
	{
		return this.readStringArray("/col_attrs/CellID");
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
	
	public Gene[] getGeneNames()
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		Gene[] genes = null;
		if(this.handle.object().isDataSet("/row_attrs/Gene"))
		{
			String[] tmp = this.handle.string().readArray("/row_attrs/Gene");
			if(genes == null) genes = new Gene[tmp.length];
			for(int i = 0; i < tmp.length; i++) 
			{
				if(genes[i] == null) genes[i] = new Gene();
				genes[i].name= tmp[i];
			}
		}
		if(this.handle.object().isDataSet("/row_attrs/Accession"))
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
		if(this.handle == null) System.err.println("Loom file is closed already");
		if(this.handle.object().exists("/__DATA_TYPES__") && !this.readOnly) ((IHDF5Writer)this.handle).object().delete("/__DATA_TYPES__");
		this.handle.close();
		this.handle = null;
	}
	
	/***** Writing *****/
	
	private void createEmptyLoomFile(String filename) // matrix.length == nb genes
    {
		if(handle != null) new ErrorJSON("Cannot open two Loom in the same Loom object, please close the first handle first");
		
    	// Create Loom file
    	IHDF5WriterConfigurator config = HDF5Factory.configure(filename);
    	this.handle = config.dontUseExtendableDataTypes().overwrite().useSimpleDataSpaceForAttributes().writer();
	    
        // Required annotations
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
	
	public void writeBlockInMatrix(float[][] data, int startIndex, HashSet<Long> filteredGenes) // matrix.length == nb genes
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
    }
	
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
	
	public void writeFloatRow(String path, float[] data, long index)
    {
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        float[][] m = new float[1][data.length];
        m[0] = data;
        ((IHDF5Writer)this.handle).float32().writeMatrixBlock(path, m, index, 0);
    }
	
	public void createEmptyMatrix(long nbGenes, long nbCells)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        HDF5FloatStorageFeatureBuilder features = new HDF5FloatStorageFeatureBuilder();
        features.noScaling().chunkedStorageLayout().datasetReplacementEnforceKeepExisting().deflateLevel((byte)2);
        ((IHDF5Writer)this.handle).float32().createMatrix("/matrix", nbGenes, nbCells, 1, (int)nbCells, features.features());
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
	
	public void createEmptyMatrix(long nbGenes, long nbCells, int chunkSize)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly matrix");
        HDF5FloatStorageFeatureBuilder features = new HDF5FloatStorageFeatureBuilder();
        features.noScaling().chunkedStorageLayout().datasetReplacementEnforceKeepExisting().deflateLevel((byte)2);
        ((IHDF5Writer)this.handle).float32().createMatrix("/matrix", nbGenes, nbCells, 64, chunkSize, features.features());
	}
		
	public void resizeDataset(String path, long x, long y)
	{
		if(this.handle == null) new ErrorJSON("Please open the Loom file first");
		if(this.readOnly) new ErrorJSON("Cannot write in readOnly");
		((IHDF5Writer)this.handle).object().setDataSetDimensions(path, new long[] { x, y });
	}
    
    /**** Static methods ******/
    
    public static void copyMetadata(String path, LoomFile from, LoomFile to)
	{
    	if(to.readOnly) new ErrorJSON("Cannot write metadata in 'to' Loom File. Please open it for writing.");
		
		long[] dim = from.getDimensions(path);
		int[] blockSize = from.getChunkSizes(path);
		int nbTotalBlocks = (int)Math.ceil((double)dim[0] / blockSize[0]); // dim[0] === Nb genes ?
		System.out.println("Reading " + nbTotalBlocks + " independent blocks...");
		boolean resizeNeeded = false;
		
		if(to.exists(path)) to.removeMetadata(path); // Remove if existing
		HDF5DataClass h5_class = from.getDataClass(path);
		switch (h5_class) 
		{
			case FLOAT:	
				if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					to.createEmptyFloat32MatrixDataset(path, dim[0], dim[1], blockSize[0], blockSize[1]); // Create matrix with same chunking
					for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
					{		
						// Retrieve the blocks that will contain all columns
						float[][] subMatrix = from.readFloatBlock(path, blockSize[0], (int)dim[1], nbBlocks); // Always read all blocks in one full line
						
						if(subMatrix.length < blockSize[0]) // Because it does not write the end of the array, if the size is less than blocksize
						{
							float[][] subMatrixTmp = new float[blockSize[0]][(int)dim[1]];
							for (int i = 0; i < subMatrix.length; i++) for (int j = 0; j < subMatrix[i].length; j++) subMatrixTmp[i][j] = subMatrix[i][j];
							subMatrix = subMatrixTmp;
							resizeNeeded = true;
						}
						
						// Writing data
						to.writeFloatBlockDataset(path, subMatrix, nbBlocks, 0);
					}
					
					// Because we added extra rows
					if(resizeNeeded) to.resizeDataset(path, dim[0], dim[1]);
				}
				else // Vector
				{
					to.writeFloatArray(path, from.readFloatArray(path));
				}
				break;
			case INTEGER:
				if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					to.createEmptyInt32MatrixDataset(path, dim[0], dim[1], blockSize[0], blockSize[1]); // Create matrix with same chunking
					for(int nbBlocks = 0; nbBlocks < nbTotalBlocks; nbBlocks++)
					{		
						// Retrieve the blocks that will contain all columns
						int[][] subMatrix = from.readIntBlock(path, blockSize[0], (int)dim[1], nbBlocks); // Always read all blocks in one full line
						
						if(subMatrix.length < blockSize[0]) // Because it does not write the end of the array, if the size is less than blocksize
						{
							int[][] subMatrixTmp = new int[blockSize[0]][(int)dim[1]];
							for (int i = 0; i < subMatrix.length; i++) for (int j = 0; j < subMatrix[i].length; j++) subMatrixTmp[i][j] = subMatrix[i][j];
							subMatrix = subMatrixTmp;
							resizeNeeded = true;
						}
						
						// Writing data
						to.writeIntBlockDataset(path, subMatrix, nbBlocks, 0);
					}
					
					// Because we added extra rows
					if(resizeNeeded) to.resizeDataset(path, dim[0], dim[1]);
				}
				else // Vector
				{
					to.writeIntArray(path, from.readIntArray(path));
				}
				break;
			case STRING:
				if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					new ErrorJSON("Cannot create a 2D array of STRING");
				}
				else // Vector
				{
					to.writeStringArray(path, from.readStringArray(path));
				}
				break;
			default:
				new ErrorJSON("This metadata type " + h5_class + " is not handled.");
				break;
		}
	}
    
    public static void copyMetadata(Metadata m, HashSet<Long> toFilter, boolean filterCells, LoomFile from, LoomFile to)
	{
    	if(to.readOnly) new ErrorJSON("Cannot write metadata in 'to' Loom File. Please open it for writing.");
		HDF5DataClass h5_class = from.getDataClass(m.path);
		long[] dim = from.getDimensions(m.path);
		boolean writeAll = true;
		if(m.on == MetaOn.CELL && filterCells) writeAll = false;
		if(m.on == MetaOn.GENE && !filterCells) writeAll = false;
		switch (h5_class) 
		{
			case FLOAT:	
				if(dim.length > 1 && dim[1] > 0) // Matrix // TODO HANDLE BIG SIZE 2D
				{
					float[][] values = from.handle.float32().readMatrix(m.path);
					if(writeAll)((IHDF5Writer)to.handle).float32().writeMatrix(m.path, values);
					else
					{
						float[][] towrite = new float[values.length - toFilter.size()][(int)dim[1]];
						int k = 0;
						for (int i = 0; i < values.length; i++) 
						{
							if(!toFilter.contains((long)i))
							{
								towrite[k] = values[i];
								k++;
							}
						}
						((IHDF5Writer)to.handle).float32().writeMatrix(m.path, towrite);
					}
				}
				else // Vector
				{
					FloatArray64 values = from.readFloatArray(m.path);
					if(writeAll) to.writeFloatArray(m.path, values);
					else
					{
						FloatArray64 towrite = new FloatArray64(values.size() - toFilter.size());
						long k = 0;
						for (long i = 0; i < values.size(); i++) 
						{
							if(!toFilter.contains(i))
							{
								towrite.set(k, values.get(i));
								k++;
							}
						}
						to.writeFloatArray(m.path, towrite);
					}
				}
				break;
			case INTEGER:
				if(dim.length > 1 && dim[1] > 0) // Matrix // TODO HANDLE BIG SIZE 2D
				{
					int[][] values = from.handle.int32().readMatrix(m.path);
					if(writeAll)((IHDF5Writer)to.handle).int32().writeMatrix(m.path, values);
					else
					{
						int[][] towrite = new int[values.length - toFilter.size()][(int)dim[1]];
						int k = 0;
						for (int i = 0; i < values.length; i++) 
						{
							if(!toFilter.contains((long)i))
							{
								towrite[k] = values[i];
								k++;
							}
						}
						((IHDF5Writer)to.handle).int32().writeMatrix(m.path, towrite);
					}
				}
				else // Vector
				{
					IntArray64 values = from.readIntArray(m.path);
					if(writeAll) to.writeIntArray(m.path, values);
					else
					{
						IntArray64 towrite = new IntArray64(values.size() - toFilter.size());
						long k = 0;
						for (long i = 0; i < values.size(); i++) 
						{
							if(!toFilter.contains(i))
							{
								towrite.set(k, values.get(i));
								k++;
							}
						}
						to.writeIntArray(m.path, towrite);
					}
				}
				break;
			case STRING:
				if(dim.length > 1 && dim[1] > 0) // Matrix
				{
					new ErrorJSON("Cannot create a 2D array of STRING");
				}
				else // Vector
				{
					StringArray64 values = from.readStringArray(m.path);
					if(writeAll) to.writeStringArray(m.path, values);
					else
					{				
						StringArray64 towrite = new StringArray64(values.size() - toFilter.size());
						long k = 0;
						for (long i = 0; i < values.size(); i++) 
						{
							if(!toFilter.contains(i))
							{
								towrite.set(k, values.get(i));
								k++;
							}
						}
						to.writeStringArray(m.path, towrite);
					}
				}
				break;
			default:
				m.type = Metatype.NOT_HANDLED;
				break;
		}
	}
    
    public static void fillLoomFile(LoomFile loom, LoomData data) // matrix.length == nb genes
    {
    	if(loom.readOnly) new ErrorJSON("Cannot write in this Loom file.");
    	
    	// Create Metadata array
    	data.meta = new ArrayList<>();
    	
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
        if(data.__no_feature != null) {
        	loom.writeIntArray("/col_attrs/_No_feature", data.__no_feature);
        	data.meta.add(new Metadata("/col_attrs/_No_feature", Metatype.NUMERIC, MetaOn.CELL, data.__no_feature.size(), 1));
        }
        if(data.__ambiguous != null) {
        	loom.writeIntArray("/col_attrs/_Ambiguous", data.__ambiguous);
        	data.meta.add(new Metadata("/col_attrs/_Ambiguous", Metatype.NUMERIC, MetaOn.CELL, data.__ambiguous.size(), 1));
        }
        if(data.__too_low_aQual != null) {
        	loom.writeIntArray("/col_attrs/_Too_low_aQual", data.__too_low_aQual);
        	data.meta.add(new Metadata("/col_attrs/_Too_low_aQual", Metatype.NUMERIC, MetaOn.CELL, data.__too_low_aQual.size(), 1));
        }
        if(data.__not_aligned != null) {
        	loom.writeIntArray("/col_attrs/_Not_aligned", data.__not_aligned);
        	data.meta.add(new Metadata("/col_attrs/_Not_aligned", Metatype.NUMERIC, MetaOn.CELL, data.__not_aligned.size(), 1));
        }
        if(data.__alignment_not_unique != null) {
        	loom.writeIntArray("/col_attrs/_Alignment_not_unique", data.__alignment_not_unique);
        	data.meta.add(new Metadata("/col_attrs/_Alignment_not_unique", Metatype.NUMERIC, MetaOn.CELL, data.__alignment_not_unique.size(), 1));
        }
        
        // Biotypes/chromosomes
        for(long j = 0; j < data.mitochondrialContent.size(); j++) data.mitochondrialContent.set(j, (float)(data.mitochondrialContent.get(j) / data.depth.get(j) * 100));
        loom.writeFloatArray("/col_attrs/_Mitochondrial_Content", data.mitochondrialContent);
 		data.meta.add(new Metadata("/col_attrs/_Mitochondrial_Content", Metatype.NUMERIC, MetaOn.CELL, data.mitochondrialContent.size(), 1));
        
        for(long j = 0; j < data.proteinCodingContent.size(); j++) data.proteinCodingContent.set(j, (float)(data.proteinCodingContent.get(j) / data.depth.get(j) * 100));
        loom.writeFloatArray("/col_attrs/_Protein_Coding_Content", data.proteinCodingContent);
 		data.meta.add(new Metadata("/col_attrs/_Protein_Coding_Content", Metatype.NUMERIC, MetaOn.CELL, data.proteinCodingContent.size(), 1));
 		
        for(long j = 0; j < data.ribosomalContent.size(); j++) data.ribosomalContent.set(j, (float)(data.ribosomalContent.get(j) / data.depth.get(j) * 100));
        loom.writeFloatArray("/col_attrs/_Ribosomal_Content", data.ribosomalContent);
 		data.meta.add(new Metadata("/col_attrs/_Ribosomal_Content", Metatype.NUMERIC, MetaOn.CELL, data.ribosomalContent.size(), 1));
     	
        // ERCCs
        if(data.erccs != null) 
        {
        	float[][] ercc_matrix = Utils.t(data.erccs.getArray());
        	((IHDF5Writer)loom.handle).float32().writeMatrix("/col_attrs/_ERCCs", ercc_matrix);
        	((IHDF5Writer)loom.handle).string().setArrayAttr("/col_attrs/_ERCCs", "colnames", data.erccs.getNames()); // Save the column names as an attribute       	
        	data.meta.add(new Metadata("/col_attrs/_ERCCs", Metatype.NUMERIC, MetaOn.CELL, ercc_matrix.length, ercc_matrix[0].length));
        }
        
        // Genes and Gene annotations     
        if(!data.removed.isEmpty()) // If gene annotations to remove from final matrix
        {
        	long writtenRows = 0;
        	long nbRows = loom.getDimensions()[0]; // Know the actual size of the array
        	StringArray64 genes = new StringArray64(nbRows);
        	StringArray64 ens = new StringArray64(nbRows);
        	DoubleArray64 sums = new DoubleArray64(nbRows);
        	LongArray64 _sumExonLength = new LongArray64(nbRows);
        	StringArray64 _biotypes = new StringArray64(nbRows);
        	StringArray64 _chromosomes = new StringArray64(nbRows);
        	for(long i = 0; i < Parameters.nGenes; i++)
        	{
        		if(!data.removed.contains(i))
        		{
        			genes.set(writtenRows, data.gene_names.get(i));
        			ens.set(writtenRows, data.ens_names.get(i));
        			sums.set(writtenRows, data.sum.get(i));
                	_sumExonLength.set(writtenRows, data.sumExonLength.get(i));
                	_biotypes.set(writtenRows, data.biotypes.get(i));
                	_chromosomes.set(writtenRows, data.chromosomes.get(i));
        			writtenRows++;
        		}
        	}
        	//HDF5GenericStorageFeatureBuilder features = new HDF5GenericStorageFeatureBuilder();
            //features.noScaling().datasetReplacementEnforceKeepExisting().deflateLevel((byte)2);
        	loom.writeStringArray("/row_attrs/Gene", genes); // HGNC names
        	data.meta.add(new Metadata("/row_attrs/Gene", Metatype.STRING, MetaOn.GENE, 1, genes.size())); 	
        	loom.writeStringArray("/row_attrs/Accession", ens); // EnsemblIDs
            data.meta.add(new Metadata("/row_attrs/Accession", Metatype.STRING, MetaOn.GENE, 1, ens.size()));
        	// Stable Ids
        	LongArray64 stable_ids = new LongArray64(ens.size());
        	for(long i = 0; i < stable_ids.size(); i++) stable_ids.set(i, i);
        	loom.writeLongArray("/row_attrs/_StableID", stable_ids);
        	data.meta.add(new Metadata("/row_attrs/_StableID", Metatype.NUMERIC, MetaOn.GENE, 1, stable_ids.size()));
            // We do not store alt names in the Loom file 
            loom.writeDoubleArray("/row_attrs/_Sum", sums);
        	data.meta.add(new Metadata("/row_attrs/_Sum", Metatype.NUMERIC, MetaOn.GENE, 1, sums.size()));
        	loom.writeLongArray("/row_attrs/_SumExonLength", _sumExonLength);
        	data.meta.add(new Metadata("/row_attrs/_SumExonLength", Metatype.NUMERIC, MetaOn.GENE, 1, _sumExonLength.size()));	
        	loom.writeStringArray("/row_attrs/_Biotypes", _biotypes);
        	Metadata m = new Metadata("/row_attrs/_Biotypes", Metatype.DISCRETE, MetaOn.GENE, 1, _biotypes.size(), _biotypes);
        	m.fillMap();
        	m.values = null;
        	data.meta.add(m);    	
        	loom.writeStringArray("/row_attrs/_Chromosomes", _chromosomes);
        	m = new Metadata("/row_attrs/_Chromosomes", Metatype.DISCRETE, MetaOn.GENE, 1, _chromosomes.size(), _chromosomes);
        	m.fillMap();
        	m.values = null;
        	data.meta.add(m);
        }
        else // Simply write the full thing
        {
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
}
