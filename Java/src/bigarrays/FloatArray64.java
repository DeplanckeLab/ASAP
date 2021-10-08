package bigarrays;

/**
 *  @author vincent.gardeux@epfl.ch
 *
 * Class for representing a float static array requiring address space larger than 32 bits. 
 */

public class FloatArray64
{
	private static final int CHUNK_SIZE = 1024*1024*512;

    private long size;
    private float[][] data;

    public FloatArray64(long size)
    {
        this.size = size;
        if(size == 0) data = null;
        else 
        {
            int chunks = (int)(size/CHUNK_SIZE);
            int remainder = (int)(size - ((long)chunks)*CHUNK_SIZE);
            data = new float[chunks+(remainder==0?0:1)][];
            for(int idx=chunks; --idx>=0; ) data[idx] = new float[(int)CHUNK_SIZE];
            if(remainder != 0) data[chunks] = new float[remainder];
        }
    }
       
    public static int chunkSize()
    {
    	return CHUNK_SIZE;
    }
    
    public float[] toArray()
    {
    	if(this.data.length == 1) return this.data[0];
    	return null;
    }
    
    public float get(long index) 
    {
        if(index < 0 || index >= size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index / CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk) * CHUNK_SIZE));
        return data[chunk][offset];
    }
    
    public float[] getByChunk(int chunk) 
    {
    	if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        return data[chunk];
    }
    
    public void set(long index, float b) 
    {
        if(index<0 || index>=size) throw new IndexOutOfBoundsException("Error attempting to access data element "+index+".  Array is "+size+" elements long.");
        int chunk = (int)(index/CHUNK_SIZE);
        int offset = (int)(index - (((long)chunk)*CHUNK_SIZE));
        data[chunk][offset] = b;
    }
    
    public void set(int chunk, float[] b) 
    {
        if(chunk >= data.length) throw new IndexOutOfBoundsException("Error attempting to access chunk "+chunk+".  Array is "+size+" elements long [" + data.length + " chunks]");
        data[chunk] = b;
    }

    public long size() 
    {
        return this.size;
    }
    
    public IntArray64 toIntArray()
    {
    	IntArray64 array = new IntArray64(this.size);
    	for(long i = 0; i < size; i++) array.set(i, (int)this.get(i));
    	return array;
    }
}