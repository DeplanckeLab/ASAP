package tools;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.IntStream;

import json.ErrorJSON;
import json.WarningJSON;
import model.Parameters;

public class Utils 
{
	private static Random rand;
	private static DecimalFormat df;
	private static final double LOG2 = Math.log(2);
	
	static 
	{
		changeFormatter(3);
	}
	
	public static void initRandomGenerator()
	{
		rand = new Random(Parameters.randomSeed);
	}
	
	public static int[] sample(List<Integer> array, int nbTimes) // with replacement
	{
		int[] res = new int[nbTimes];
		for (int i = 0; i < nbTimes; i++) res[i] = array.get(getRandomGenerator().nextInt(array.size()));
		return res;
	}
	
	public static void changeFormatter(int nbSignificantDigits)
	{
		StringBuilder sb = new StringBuilder("0");
		if(nbSignificantDigits > 0) sb.append(".");
		for(int i = 0; i < nbSignificantDigits; i++) sb.append("#");
		if(Parameters.scientific) sb.append("E0");
		DecimalFormatSymbols symbols = DecimalFormatSymbols.getInstance();
		symbols.setDecimalSeparator('.');
		df = new DecimalFormat(sb.toString(), symbols);
	}
	
	public static void writeJSON(StringBuilder content, String outputJSONFile)
	{
		String out = content.toString();
		if(WarningJSON.isAnyWarning())
		{
			int index = out.lastIndexOf("}");
			out = out.substring(0, index);
			out = out + "," + WarningJSON.getJSON() + "}";
		}
		if(outputJSONFile == null) System.out.println(out);
		else
		{
			try
			{
	    		BufferedWriter bw = new BufferedWriter(new FileWriter(outputJSONFile));
	    		bw.write(out);
	        	bw.close();
			}
			catch(IOException ioe) { new ErrorJSON(ioe.getMessage()); }
		}
	}
	
	public static void writePlain(StringBuilder content, String outputFile)
	{
		String out = content.toString();
		if(WarningJSON.isAnyWarning())
		{
			System.err.println(WarningJSON.getJSON());
		}
		if(outputFile == null) System.out.println(out);
		else
		{
			try
			{
	    		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
	    		bw.write(out);
	        	bw.close();
			}
			catch(IOException ioe) { new ErrorJSON(ioe.getMessage()); }
		}
	}
	
	public static String handleSpecialCharacters(String content)
	{
		// Remove return characters
		if(content.contains("\r")) content = content.replaceAll("\r", "");
		if(content.contains("\n")) content = content.replaceAll("", "");
		
		// Handle special characters
		// First replace all \ by \\
		content = content.replaceAll("\\\\", "\\\\\\\\");
		// Then replace " by \" 
		content = content.replaceAll("\"", "\\\\\"");
		return content;
	}
	
	public static double[] p_adjust(final double[] pvalues, String adjMethod)
    {
        if (!adjMethod.equals("BH") && !adjMethod.equals("fdr") && !adjMethod.equals("bonferroni") && !adjMethod.equals("none"))
        {
            System.out.println("This adjustment method is not implemented.");
            System.exit(0);
        }
        if (adjMethod.equals("fdr")) adjMethod = "BH";
        if (adjMethod.equals("none")) return pvalues;

        Integer[] idx = new Integer[pvalues.length];
        for (int i = 0; i < idx.length; i++) idx[i] = i;
        double[] adj_p_value = new double[pvalues.length];

        Arrays.sort(idx, new Comparator<Integer>()
        { // I sort the indexes with respect to the double values
            @Override
            public int compare(Integer o1, Integer o2)
            {
                return Double.compare(pvalues[o2], pvalues[o1]);
            }
        });
        int[] index = new int[idx.length];
        for (int i = 0; i < index.length; i++)
        {
            index[idx[i]] = i;
        }

        for (int i = 0; i < pvalues.length; i++)
        {
            if (adjMethod.equals("BH"))
            {
                adj_p_value[i] = pvalues.length / (pvalues.length - (double) index[i]) * pvalues[i];
            }
            if (adjMethod.equals("bonferroni"))
            {
                adj_p_value[i] = pvalues.length * pvalues[i];
            }
        }

        double min = Double.MAX_VALUE;
        for (int i = 0; i < index.length; i++) // cummin
        {
            double adjP = adj_p_value[idx[i]];
            if (adjP < min)
            {
                min = adjP;
            } else
            {
                adjP = min;
            }
            adj_p_value[idx[i]] = Math.min(1, adjP);
        }

        return adj_p_value;
    }
	
	public static float[] p_adjust_F(final float[] pvalues, String adjMethod)
    {
        if (!adjMethod.equals("BH") && !adjMethod.equals("fdr") && !adjMethod.equals("bonferroni") && !adjMethod.equals("none"))
        {
            System.out.println("This adjustment method is not implemented.");
            System.exit(0);
        }
        if (adjMethod.equals("fdr")) adjMethod = "BH";
        if (adjMethod.equals("none")) return pvalues;

        Integer[] idx = new Integer[pvalues.length];
        for (int i = 0; i < idx.length; i++) idx[i] = i;
        float[] adj_p_value = new float[pvalues.length];

        Arrays.sort(idx, new Comparator<Integer>()
        { // I sort the indexes with respect to the double values
            @Override
            public int compare(Integer o1, Integer o2)
            {
                return Double.compare(pvalues[o2], pvalues[o1]);
            }
        });
        int[] index = new int[idx.length];
        for (int i = 0; i < index.length; i++)
        {
            index[idx[i]] = i;
        }

        for (int i = 0; i < pvalues.length; i++)
        {
            if (adjMethod.equals("BH"))
            {
                adj_p_value[i] = pvalues.length / (pvalues.length - (float) index[i]) * pvalues[i];
            }
            if (adjMethod.equals("bonferroni"))
            {
                adj_p_value[i] = pvalues.length * pvalues[i];
            }
        }

        double min = Double.MAX_VALUE;
        for (int i = 0; i < index.length; i++) // cummin
        {
            double adjP = adj_p_value[idx[i]];
            if (adjP < min)
            {
                min = adjP;
            } else
            {
                adjP = min;
            }
            adj_p_value[idx[i]] = (float)Math.min(1, adjP);
        }

        return adj_p_value;
    }
	
	public static double log2(double x)
	{
	    return Math.log(x) / LOG2;
	}
	
	public static void setSeed(int seed)
	{
		rand = new Random(seed);
	}
	
	public static String format(double n)
	{
		return df.format(n);
	}
	
	public static String format(String n)
	{
		if(n.equals("NaN")) return "null";
		n = n.replaceAll(",", ".");
		return Utils.format(Double.parseDouble(n));
	}
	
	public static double [][] t(double [][] matrix)
	{
		double[][] output = new double[matrix[0].length][matrix.length];
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[i].length; j++)
			{
				output[j][i] = matrix[i][j];
			}
		}
		return output;
	}
	
	public static float [][] t(float [][] matrix)
	{
		float[][] output = new float[matrix[0].length][matrix.length];
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[i].length; j++)
			{
				output[j][i] = matrix[i][j];
			}
		}
		return output;
	}
	
	/**
	 * Fastest implementation to count non-empty lines of a text file
	 * @param is InputStream of the text file to read
	 * @return number of non-empty lines
	 * @throws IOException
	 */
	public static long countNonEmptyLines(InputStream is) throws IOException
	{
	    byte[] c = new byte[1024];
	    int count = 0;
	    byte lastc = '\n';
	    int readChars = is.read(c);
	    while(readChars != -1) // In some cases, it can return something else than 1024. Returns -1 if empty.
	    {
	    	for(int i=0; i<readChars; i++)
	    	{
	    		if((c[i] == '\n' || c[i] == '\r') && lastc != '\n' && lastc != '\r') count++;
	    		lastc = c[i];
	    	}
	    	readChars = is.read(c);
	    } 
	    is.close();
	    if(lastc == '\n' || lastc == '\r') return count;
	    return count + 1;
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
	
	public static float mean(float[] data)
	{
		float sum = 0;
		for(float a : data) sum += a;
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
	
	public static double median(ArrayList<Double> data)
	{
		Collections.sort(data);
		if (data.size() % 2 == 0) return ((double)data.get(data.size()/2) + (double)(data.get(data.size()/2 - 1)))/2;
		return (double) data.get(data.size()/2);
	}
	
	public static Double quartile(double[] data, float quartile)
	{
	    Arrays.sort(data);
	    if(quartile < 0 || quartile > 1) return null;
	    int q = (int) Math.round(data.length * quartile);
	    if(q >= data.length) q = data.length - 1;
	    return data[q];
	}
	
	public static int[] quantiles(float[] array, int nbin)
	{
		int[] ranks = Utils.rank(array);
		int[] quantiles = new int[ranks.length];
		int binLength = ranks.length / nbin;
		for (int i = 0; i < quantiles.length; i++) quantiles[i] = (ranks[i] / binLength) + 1;
		return quantiles;
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
	
	public static StringBuilder toString(float[] array)			
	{
		StringBuilder res = new StringBuilder("[");
		if(array != null)
		{
			String prefix = "";
			for(float f:array)
			{
				res.append(prefix).append(f);
				prefix = ",";
			}
		}
		res.append("]");
		return res;
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
	
	// Be careful, order != rank
	public static int[] order(int[] array) // return array of indexes
	{
		int[] sortedIndexes = IntStream.range(0, array.length).boxed().sorted((i, j) -> array[i] - array[j]).mapToInt(ele -> ele).toArray();
		return sortedIndexes;
	}
	
	public static int[] order(double[] array, boolean reversed) // return array of indexes
	{
		HashMap<Integer, Double> map = new HashMap<>(array.length);
		for (int i = 0; i < array.length; i++) map.put(i, array[i]); // Copy input
		return sortD(map, reversed);
	}
	
	public static int[] order(float[] array, boolean reversed) // return array of indexes
	{
		HashMap<Integer, Float> map = new HashMap<>(array.length);
		for (int i = 0; i < array.length; i++) map.put(i, array[i]); // Copy input
		return sortF(map, reversed);
	}
	
	public static boolean contains(String[] array, String value)
	{
		for(String val:array) if(val.equals(value)) return true;
		return false;
	}
	
	public static String[] sortI(Map<String, Integer> map)
	{
		return sortI(map, false);
	}
	
	public static String[] sortL(Map<String, Long> map)
	{
		return sortL(map, false);
	}

	public static String[] sortI(Map<String, Integer> map, boolean reversed)
	{
		List<Integer> values = new ArrayList<>(map.values());
		if(reversed) Collections.sort(values, Collections.reverseOrder());
		else Collections.sort(values);
		String[] sortedIndexes = new String[values.size()];
		for(String key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	public static String[] sortL(Map<String, Long> map, boolean reversed)
	{
		List<Long> values = new ArrayList<>(map.values());
		if(reversed) Collections.sort(values, Collections.reverseOrder());
		else Collections.sort(values);
		String[] sortedIndexes = new String[values.size()];
		for(String key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	public static int[] sortD(HashMap<Integer, Double> map, boolean reversed)
	{
		List<Double> values = new ArrayList<>(map.values());
		if(reversed) Collections.sort(values, Collections.reverseOrder());
		else Collections.sort(values);
		int[] sortedIndexes = new int[values.size()];
		for(Integer key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	public static int[] sortF(HashMap<Integer, Float> map, boolean reversed)
	{
		List<Float> values = new ArrayList<>(map.values());
		if(reversed) Collections.sort(values, Collections.reverseOrder());
		else Collections.sort(values);
		int[] sortedIndexes = new int[values.size()];
		for(Integer key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	
	public static int[] rank(float[] array)
	{
	    int[] R = new int[array.length];
	    Integer [] I = new Integer[array.length];
	    for(int i = 0; i < array.length; i++) I[i] = i;
	    Arrays.sort(I, (i0, i1) -> (int) Math.signum(array[i0]-array[i1]));
	    int j = 0;
	    for(int i = 0; i < array.length; i++)
	    {
	        if(array[I[i]] != array[I[j]]) j = i;
	        R[I[i]] = j;
	    }
	    return R;
	}
	
	/*public static int[] rank(float[] array, boolean reversed)
	{
		HashMap<Float, List<Integer>> map = new HashMap<>(array.length);
		for (int i = 0; i < array.length; i++) 
		{
			List<Integer> indexes = map.get(array[i]);
			if(indexes == null) indexes = new ArrayList<>();
			indexes.add(i);
			map.put(array[i], indexes);
		}
		return rankF(map, reversed);
	}
	
	public static int[] rankF(HashMap<Float, List<Integer>> map, boolean reversed)
	{
		Set<Float> values = new HashSet<>(map.keySet());
		if(reversed) Collections.sort(values, Collections.reverseOrder());
		else Collections.sort(values);
		int[] rankedIndexes = new int[values.size()];
		for(Float value:values) 
		{
			int index = values.indexOf(map.get(key));
			rankedIndexes[index] = key;
			values.set(index, null);
		}
		return rankedIndexes;
	}*/
	
	public static String[] sortKeys(Set<String> keySet)
	{
		String[] keys = keySet.toArray(new String[keySet.size()]);
		Arrays.sort(keys);
		return keys;
	}
	
	public static float[] colMeans(float[][] input)
	{
		if(input == null || input.length == 0) return null;
		int nCols = input[0].length;
		float[] res = new float[nCols];
		for (int j = 0; j < nCols; j++) 
		{
			for (int i = 0; i < input.length; i++) res[j] += input[i][j];
			res[j] = res[j] / input.length;
		}
		return res;
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
