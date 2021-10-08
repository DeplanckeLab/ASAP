package bigarrays;

import java.util.Arrays;
import java.util.Iterator;

/**
 *  @author vincent.gardeux@epfl.ch
 *
 * Class for representing a String static array requiring address space larger than 32 bits. 
 */

public class StringArray64 implements Iterable<String>
{
	private static final int CHUNK_SIZE = 1024*1024*512; //1GiB

    private long size;
    private String[][] data;

    public StringArray64(String[] array) // TODO remove this
    {
        this.size = array.length;
        if(size == 0) data = null;
        else 
        {
            data = new String[1][array.length];
            data[0] = array;
        }
    }
    
    public StringArray64(long size)
    {
        this.size = size;
        if(size == 0) data = null;
        else 
        {
            int chunks = (int)(size/CHUNK_SIZE);
            int remainder = (int)(size - ((long)chunks)*CHUNK_SIZE);
            data = new String[chunks+(remainder==0?0:1)][];
            for(int idx=chunks; --idx>=0; ) { data[idx] = new String[(int)CHUNK_SIZE]; Arrays.fill(data[idx], ""); }
            if(remainder != 0) { data[chunks] = new String[remainder]; Arrays.fill(data[chunks], ""); }
        }
    }
    
    public static int chunkSize()
    {
    	return CHUNK_SIZE;
    }
    
    public String get(long index) 
    {
        if(index < 0 || index >= size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index / CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk) * CHUNK_SIZE));
        return data[chunk][offset];
    }
    
    public String[] getByChunk(int chunk) 
    {
    	if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        return data[chunk];
    }
    
    public void set(long index, String b) 
    {
        if(index<0 || index>=size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index/CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk)*CHUNK_SIZE));
        data[chunk][offset] = b;
    }
    
    public void set(int chunk, String[] b) 
    {
        if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        data[chunk] = b;
    }

    public StringArray64 copy()
    {
    	StringArray64 copy = new StringArray64(this.size);
    	for(long index = 0; index < this.size; index++) copy.set(index, this.get(index));
    	return copy;
    }
    
    public long size() 
    {
        return this.size;
    }

	@Override
	public Iterator<String> iterator() 
	{
		return new StringArray64Iterator();
	}
	
	private class StringArray64Iterator implements Iterator<String>
	{
        private long position = 0;
 
        public boolean hasNext() 
        {
            if (position < size) return true;
            return false;
        }
 
        public String next() 
        {
            if(this.hasNext()) return get(position++);
            return null;
        }
    }
}