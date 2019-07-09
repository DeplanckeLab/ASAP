package bigarrays;

/**
 *  @author vincent.gardeux@epfl.ch
 *
 * Class for representing a double static array requiring address space larger than 32 bits. 
 */

public class DoubleArray64
{
	private static final int CHUNK_SIZE = 1024*1024*1024; //1GiB

    private long size;
    private double[][] data;

    public DoubleArray64(long size)
    {
        this.size = size;
        if(size == 0) data = null;
        else 
        {
            int chunks = (int)(size/CHUNK_SIZE);
            int remainder = (int)(size - ((long)chunks)*CHUNK_SIZE);
            data = new double[chunks+(remainder==0?0:1)][];
            for(int idx=chunks; --idx>=0; ) data[idx] = new double[(int)CHUNK_SIZE];
            if(remainder != 0) data[chunks] = new double[remainder];
        }
    }
    
    public static int chunkSize()
    {
    	return CHUNK_SIZE;
    }
    
    public double get(long index) 
    {
        if(index < 0 || index >= size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index / CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk) * CHUNK_SIZE));
        return data[chunk][offset];
    }
    
    public double[] getByChunk(int chunk) 
    {
    	if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        return data[chunk];
    }
    
    public void set(long index, double b) 
    {
        if(index<0 || index>=size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index/CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk)*CHUNK_SIZE));
        data[chunk][offset] = b;
    }
    
    public void set(int chunk, double[] b) 
    {
        if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        data[chunk] = b;
    }

    public long size() 
    {
        return this.size;
    }
}