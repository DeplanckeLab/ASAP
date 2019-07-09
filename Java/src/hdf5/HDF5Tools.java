package hdf5;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import hdf5.h510x.H510xHandler;
import hdf5.loom.LoomFile;
import json.ErrorJSON;
import model.Parameters;
import parsing.model.FileType;

public class HDF5Tools 
{
	public static void checkIfHDF5orLoom()
	{
		Parameters.fileType = null;
        try
        {
        	if(HDF5Factory.isHDF5File(Parameters.fileName))
        	{
        		if(H510xHandler.is10xFormatOK()) Parameters.fileType = FileType.H5_10x;
            	else 
            	{ 
            		LoomFile loom = new LoomFile("r", Parameters.fileName);
            		if(loom.isLoomFormatOK()) Parameters.fileType = FileType.LOOM;
            		loom.close();
            	}
        	}
        }
	    catch(Exception e)
	    {
	    	new ErrorJSON(e.getMessage());
	    }
	}
}


