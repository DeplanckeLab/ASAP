package tools;


import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

public class Utils 
{
	private static Random rand = new Random();
	
	public static void setSeed(int seed)
	{
		rand = new Random(seed);
	}
	
	public static void listdirs(String directoryName, ArrayList<File> files) 
	{
	    File directory = new File(directoryName);
	    File[] fList = directory.listFiles();
	    for (File file : fList) if (file.isDirectory()) files.add(file);
	}
	
	public static void listfiles(String directoryName, ArrayList<File> files) 
	{
	    File directory = new File(directoryName);
	    File[] fList = directory.listFiles();
	    for (File file : fList) if (file.isFile()) files.add(file);
	}
    
	public static Random getRandomGenerator()
	{
		return rand;
	}
	
	public static double mean(double[] data)
	{
		double sum = 0.0;
		for(double a : data) sum += a;
		return sum/data.length;
	}
	
	public static double mean(ArrayList<Integer> data)
	{
		double sum = 0.0;
		for(Integer a : data) sum += a;
		return sum/(double)data.size();
	}
	
	public static int sum(int[] array)
	{
		int sum = 0;
		for(int i:array) sum+=i;
		return sum;
	}
	
	public static double mean(int[] data)
	{
		double sum = 0.0;
		for(double a : data) sum += a;
		return sum/data.length;
	}

	public static double median(int[] data)
	{
		Arrays.sort(data);
		if (data.length % 2 == 0) return ((double)data[data.length/2] + (double)data[data.length/2 - 1])/2;
		return (double) data[data.length/2];
	}
	
	public static double var(double[] data, double mean)
	{
		double var = 0;
		for(double a :data) var += (a-mean)*(a-mean);
		return var/(data.length-1);
	}
	
	public static double var(int[] data, double mean)
	{
		double var = 0;
		for(double a :data) var += (a-mean)*(a-mean);
		return var/(data.length-1);
	}
	
	public static double sd(double[] data, double mean)
	{
		double var = var(data, mean);
		return Math.sqrt(var);
	}
	
	public static double sd(int[] data, double mean)
	{
		double var = var(data, mean);
		return Math.sqrt(var);
	}
	
	public static double var(double[] data)
	{
	    double mean = 0;
	    double M2 = 0;
	    if(data.length < 2) return 0;
	    for (int i = 0; i < data.length; i++) 
	    {
	        double delta = data[i] - mean;
	        mean = mean + delta/(i+1);
	        M2 = M2 + delta*(data[i] - mean);
		}
	    return M2/(data.length - 1);
	}
	
	public static double sd(double[] data)
	{
		double var = var(data);
		return Math.sqrt(var);
	}
	
	public static double cv(double[] data) // Coefficient of Variation
	{
		double mu = mean(data);
		double theta = sd(data);
		return theta / mu;
	}
	
	public static double[] toArray(ArrayList<Double> data)
	{
		double[] res = new double[data.size()];
		for (int i = 0; i < res.length; i++) res[i] = data.get(i);
		return res;
	}
	
    public static double cov(double[] x, double[] y)
    {
        double result = 0d;
        int length = x.length;
        double xMean = mean(x);
        double yMean = mean(y);
        for (int i = 0; i < length; i++) 
        {
        	double xDev = x[i] - xMean;
        	double yDev = y[i] - yMean;
        	result += (xDev * yDev - result) / (i + 1);
        }
        return result * ((double) length / (double)(length - 1)); // Bias correction
    }
    
    public static double cov(double[] x, double[] y, double xMean, double yMean)
    {
        double result = 0d;
        int length = x.length;
        for (int i = 0; i < length; i++) 
        {
        	double xDev = x[i] - xMean;
        	double yDev = y[i] - yMean;
        	result += (xDev * yDev - result) / (i + 1);
        }
        return result * ((double) length / (double)(length - 1)); // Bias correction
    }
	
	public static String[] sortD(Map<String, Double> map)
	{
		List<Double> values = new ArrayList<>(map.values());
		Collections.sort(values, Collections.reverseOrder());
		String[] sortedIndexes = new String[values.size()];
		for(String key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	public static boolean contains(String[] array, String value)
	{
		for(String val:array) if(val.equals(value)) return true;
		return false;
	}

	public static String[] sortI(Map<String, Integer> map)
	{
		List<Integer> values = new ArrayList<>(map.values());
		Collections.sort(values, Collections.reverseOrder());
		String[] sortedIndexes = new String[values.size()];
		for(String key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	public static String[] sortL(Map<String, Long> map)
	{
		List<Long> values = new ArrayList<>(map.values());
		Collections.sort(values, Collections.reverseOrder());
		String[] sortedIndexes = new String[values.size()];
		for(String key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	public static String[] sortKeys(Set<String> keySet)
	{
		String[] keys = keySet.toArray(new String[keySet.size()]);
		Arrays.sort(keys);
		return keys;
	}
	
	public static String toReadableTime(long ms)
	{
		if(ms < 1000) return ""+ms+" ms";
		long s = ms / 1000;
		ms = ms % 1000;
		if(s < 60) return s+" s "+ms+" ms";
		long mn = s / 60;
		s = s % 60;
		if(mn < 60) return mn +" mn "+s+" s "+ms+" ms";
		long h = mn / 60;
		mn = mn % 60;
		if(h < 24) return h +" h "+ mn +" mn "+s+" s "+ms+" ms";
		long d = h / 24;
		h = h % 24;
		return d+ " d " + h +" h "+ mn +" mn "+s+" s "+ms+" ms";
	}
}
