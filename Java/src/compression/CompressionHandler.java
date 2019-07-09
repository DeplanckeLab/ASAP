package compression;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.apache.commons.compress.archivers.ArchiveEntry;
import org.apache.commons.compress.archivers.ArchiveException;
import org.apache.commons.compress.archivers.ArchiveInputStream;
import org.apache.commons.compress.archivers.ArchiveStreamFactory;
import org.apache.commons.compress.compressors.CompressorException;
import org.apache.commons.compress.compressors.CompressorInputStream;
import org.apache.commons.compress.compressors.CompressorStreamFactory;

import json.ErrorJSON;
import model.Parameters;
import parsing.model.FileType;

public class CompressionHandler 
{
	public static InputStream getReader()
	{
		try
		{
			InputStream in = new BufferedInputStream(Files.newInputStream(Paths.get(Parameters.fileName)));
			if(Parameters.fileType != null)
			{
				switch(Parameters.fileType)
				{
					case ARCHIVE: return getReaderArchived(in);
					case COMPRESSED: return getReaderCompressed(in);
					case ARCHIVE_COMPRESSED: return getReaderArchiveCompressed(in);
					default: return null; // This should not happen
				}
			}
			else
			{
				InputStream reader = getReaderArchiveCompressed(in);
				if(reader != null) return reader;
				in = new BufferedInputStream(Files.newInputStream(Paths.get(Parameters.fileName)));
				reader = getReaderCompressed(in);
				if(reader != null) return reader;
				reader = getReaderArchived(in);
				return reader;
			}

		}
		catch (IOException ioe) 
		{
			new ErrorJSON("Error while opening the file: " + ioe.getMessage());
		}
		return null;
	}
	
	public static ArrayList<String> getList()
	{
		ArrayList<String> res = null;
		try
		{
			InputStream in = new BufferedInputStream(Files.newInputStream(Paths.get(Parameters.fileName)));
			res = getListArchiveCompressed(in);
			in = new BufferedInputStream(Files.newInputStream(Paths.get(Parameters.fileName)));
			if(res == null) res = getListCompressed(in);
			if(res == null) res = getListArchived(in);
		}
		catch (IOException ioe) 
		{
			new ErrorJSON("Error while opening the file: " + ioe.getMessage());
		}
		return res;
	}
	
	private static InputStream getReaderArchived(InputStream in)
	{
		if(Parameters.selection == null) { System.err.println("Cannot get a reader from an archive without a selected item"); return null; }
		try
		{
			ArchiveInputStream i = new ArchiveStreamFactory().createArchiveInputStream(in);
			ArchiveEntry entry = null;
		    while ((entry = i.getNextEntry()) != null) 
		    {
		        if(!entry.isDirectory() && entry.getName().equals(Parameters.selection)) 
		        {
				    System.out.println("Detected as archive file");
				    Parameters.fileType = FileType.ARCHIVE;
		        	return i;
		        }
		    }
		    new ErrorJSON("Did not found your entry in the archive: " + Parameters.selection);
		}
		catch (IOException ioe) { new ErrorJSON("Error while listing the content of the file: " + ioe.getMessage()); }
		catch (ArchiveException ae) { return null; }
		return null;
	}
	
	private static ArrayList<String> getListArchived(InputStream in)
	{
		ArrayList<String> entries = new ArrayList<>();
		try
		{
			ArchiveInputStream i = new ArchiveStreamFactory().createArchiveInputStream(in);
			ArchiveEntry entry = null;
		    while ((entry = i.getNextEntry()) != null) 
		    {
		        if(!i.canReadEntryData(entry)) System.err.println("[IGNORED]Archive contains a non-readable entry " + entry.getName());
		        if(!entry.isDirectory()) entries.add(entry.getName());
		    }
		    System.out.println("Detected as archive file (list of items to be given)");
		    Parameters.fileType = FileType.ARCHIVE;
		}
		catch (IOException ioe) { new ErrorJSON("Error while listing the content of the file: " + ioe.getMessage()); }
		catch (ArchiveException ae) { return null; }
		return entries;
	}
	
	private static InputStream getReaderCompressed(InputStream in)
	{
		try
		{
			CompressorInputStream i = new CompressorStreamFactory().createCompressorInputStream(in);
		    System.out.println("Detected as compressed file");
		    Parameters.fileType = FileType.COMPRESSED;
			return i;
		}
		catch (CompressorException ae) { return null; }
	}
	
	private static ArrayList<String> getListCompressed(InputStream in) // Actually only one file to return
	{
		ArrayList<String> entries = new ArrayList<>();
		try
		{
			new CompressorStreamFactory().createCompressorInputStream(in);
			entries.add(Parameters.fileName.substring(Parameters.fileName.lastIndexOf("/") + 1, Parameters.fileName.lastIndexOf(".")));
			System.out.println("Detected as compressed file (only one item compressed)");
			Parameters.fileType = FileType.COMPRESSED;
		}
		catch (CompressorException ce) { return null; }
		return entries;
	}
	
	private static InputStream getReaderArchiveCompressed(InputStream in)
	{
		if(Parameters.selection == null) { System.err.println("Cannot get a reader from an archive without a selected item"); return null; }
		try
		{
			CompressorInputStream cis = new CompressorStreamFactory().createCompressorInputStream(in);
			ArchiveInputStream i = new ArchiveStreamFactory().createArchiveInputStream(new BufferedInputStream(cis));
			ArchiveEntry entry = null;
		    while ((entry = i.getNextEntry()) != null) 
		    {
		        if(!entry.isDirectory() && entry.getName().equals(Parameters.selection)) 
		        {
				    System.out.println("Detected as archive-compressed file");
				    Parameters.fileType = FileType.ARCHIVE_COMPRESSED;
		        	return i;
		        }
		    }
		    new ErrorJSON("Did not found your entry in the archive: " + Parameters.selection);
		}
		catch (IOException ioe) { new ErrorJSON("Error while listing the content of the file: " + ioe.getMessage()); }
		catch (ArchiveException ae) { return null; }
		catch (CompressorException ce) { return null; }
		return null;
	}
	
	private static ArrayList<String> getListArchiveCompressed(InputStream in)
	{
		ArrayList<String> entries = new ArrayList<>();
		try
		{
			CompressorInputStream cis = new CompressorStreamFactory().createCompressorInputStream(in);
			ArchiveInputStream i = new ArchiveStreamFactory().createArchiveInputStream(new BufferedInputStream(cis));
			ArchiveEntry entry = null;
		    while ((entry = i.getNextEntry()) != null) 
		    {
		        if (!i.canReadEntryData(entry)) System.err.println("[IGNORED]Archive contains a non-readable entry " + entry.getName());
		        if(!entry.isDirectory()) entries.add(entry.getName());
		    }
		    System.out.println("Detected as archive compressed file (list of items to be given)");
		    Parameters.fileType = FileType.ARCHIVE_COMPRESSED;
		}
		catch (IOException ioe) { new ErrorJSON("Error while listing the content of the file: " + ioe.getMessage()); }
		catch (ArchiveException ae) { return null; }
		catch (CompressorException ce) { return null; }
		return entries;
	}
}
