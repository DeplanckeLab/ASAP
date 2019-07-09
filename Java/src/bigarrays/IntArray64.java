package bigarrays;

/**
 *  @author vincent.gardeux@epfl.ch
 *
 * Class for representing a integer static array requiring address space larger than 32 bits. 
 */

public class IntArray64 
{
	private static final int CHUNK_SIZE = 1024*1024*512;

    private long size;
    private int[][] data;

    public IntArray64(long size)
    {
        this.size = size;
        if(size == 0) data = null;
        else 
        {
            int chunks = (int)(size/CHUNK_SIZE);
            int remainder = (int)(size - ((long)chunks)*CHUNK_SIZE);
            data = new int[chunks+(remainder==0?0:1)][];
            for(int idx=chunks; --idx>=0; ) data[idx] = new int[(int)CHUNK_SIZE];
            if(remainder != 0) data[chunks] = new int[remainder];
        }
    }
    
    public static int chunkSize()
    {
    	return CHUNK_SIZE;
    }
    
    public int get(long index) 
    {
        if(index < 0 || index >= size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index / CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk) * CHUNK_SIZE));
        return data[chunk][offset];
    }
    
    public int[] getByChunk(int chunk) 
    {
    	if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        return data[chunk];
    }
    
    public void set(long index, int b) 
    {
        if(index<0 || index>=size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index/CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk)*CHUNK_SIZE));
        data[chunk][offset] = b;
    }
    
    public void set(int chunk, int[] b) 
    {
        if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        data[chunk] = b;
    }

    public long size() 
    {
        return this.size;
    }
}
