package bigarrays;

import java.util.ArrayList;
import java.util.Iterator;

/**
 *  @author vincent.gardeux@epfl.ch
 *
 * Class for representing a long static array requiring address space larger than 32 bits. 
 */

public class LongArray64 implements Iterable<Long>
{
	private static final int CHUNK_SIZE = 1024*1024*1024; //1GiB

    private long size;
    private long[][] data;

    public LongArray64(long size)
    {
        this.size = size;
        if(size == 0) data = null;
        else 
        {
            int chunks = (int)(size/CHUNK_SIZE);
            int remainder = (int)(size - ((long)chunks)*CHUNK_SIZE);
            data = new long[chunks+(remainder==0?0:1)][];
            for(int idx=chunks; --idx>=0; ) data[idx] = new long[(int)CHUNK_SIZE];
            if(remainder != 0) data[chunks] = new long[remainder];
        }
    }
    
    public static int chunkSize()
    {
    	return CHUNK_SIZE;
    }
    
    public long get(long index) 
    {
        if(index < 0 || index >= size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index / CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk) * CHUNK_SIZE));
        return data[chunk][offset];
    }
    
    public long[] getByChunk(int chunk) 
    {
    	if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        return data[chunk];
    }
    
    public void set(long index, long b) 
    {
        if(index<0 || index>=size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index/CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk)*CHUNK_SIZE));
        data[chunk][offset] = b;
    }
    
    public void set(int chunk, long[] b) 
    {
        if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        data[chunk] = b;
    }

    public long size() 
    {
        return this.size;
    }
    
    public static LongArray64 convertFrom(ArrayList<Long> list)
    {
    	LongArray64 array = new LongArray64(list.size());
    	for(int i = 0; i < list.size(); i++) array.set(i, list.get(i));
    	return array;
    }
    
	@Override
	public Iterator<Long> iterator() 
	{
		return new LongArray64Iterator();
	}
	
	private class LongArray64Iterator implements Iterator<Long>
	{
        private long position = 0;
 
        public boolean hasNext() 
        {
            if (position < size) return true;
            return false;
        }
 
        public Long next() 
        {
            if(this.hasNext()) return get(position++);
            return null;
        }
    }
}
