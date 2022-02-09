package tools;

import model.Parameters;

/**
 * <p>
 * The code used to calculate a Fisher p-value comes originally from a
 * <a href="http://infofarm.affrc.go.jp/~kadasowa/fishertest.htm">JavaScript program</a>
 * by T. Kadosawa (kadosawa@niaes.affrc.go.jp) and was adapted from a Java program developed by David Hopwood
 */
public class FisherExactTest_2X2
{
	public enum Alternative{GREATER,LESS,TWO_SIDED};
	private static double[] logFactorial;
	private static FisherExactTest_2X2 singleton = null;
	
	public static FisherExactTest_2X2 getSingleton()
	{
		if(singleton == null) singleton = new FisherExactTest_2X2(Parameters.maxUniverse);
		return singleton;
	}
	
	private FisherExactTest_2X2(int universe)
	{
		logFactorial = new double[universe]; // Precompute this for faster computations afterwards (max should be >= size universe)
		logFactorial[0] = 0.0;
        for(int i = 1; i < logFactorial.length; i++) logFactorial[i] = logFactorial[i-1] + Math.log(i);
	}
	    
    public double fisher(int a, int b, int c, int d, Alternative alt) // oneTailed = less, not greater
    {
        if (a * d > b * c || (alt == Alternative.GREATER && a * d < b * c)) // Odds Ratio > 1 or < 1
        {
            a = a + b; b = a - b; a = a - b; 
            c = c + d; d = c - d; c = c - d;
        }
        if (a > d) { a = a + d; d = a - d; a = a - d; }
        if (b > c) { b = b + c; c = b - c; b = b - c; }

        int a_org = a;
        double p_sum = 0.0d;

        double p = fisherSub(a, b, c, d);
        double p_1 = p;

        while (a >= 0)
        {
            p_sum += p;
            if (a == 0) break;
            --a; ++b; ++c; --d;
            p = fisherSub(a, b, c, d);
        }
        if(alt != Alternative.TWO_SIDED) return p_sum;

        a = b; b = 0; c = c - a; d = d + a;
        p = fisherSub(a, b, c, d);

        while (p < p_1) 
        {
            if (a == a_org) break;
            p_sum += p;
            --a; ++b; ++c; --d;
            p = fisherSub(a, b, c, d);
        }
        return p_sum;
    }

    public double fisher2(int a, int b, int c, int d, boolean oneTailed) // oneTailed = less, not greater
    {
    	System.out.println((double)(a * d) / (b * c));
        if (a * d < b * c) // Odds Ratio < 1
        {
            a = a + b; b = a - b; a = a - b; 
            c = c + d; d = c - d; c = c - d;
        }
        if (a > d) { a = a + d; d = a - d; a = a - d; }
        if (b > c) { b = b + c; c = b - c; b = b - c; }

        int a_org = a;
        double p_sum = 0.0d;

        double p = fisherSub(a, b, c, d);
        double p_1 = p;

        while (a >= 0)
        {
            p_sum += p;
            if (a == 0) break;
            --a; ++b; ++c; --d;
            p = fisherSub(a, b, c, d);
        }
        if(oneTailed) return p_sum;

        a = b; b = 0; c = c - a; d = d + a;
        p = fisherSub(a, b, c, d);

        while (p < p_1) 
        {
            if (a == a_org) break;
            p_sum += p;
            --a; ++b; ++c; --d;
            p = fisherSub(a, b, c, d);
        }
        return p_sum;
    }
    
    private double fisherSub(int a, int b, int c, int d)
    {
        return Math.exp(logFactorial[a + b] + logFactorial[c + d] + logFactorial[a + c] + logFactorial[b + d] - logFactorial[a + b + c + d] - logFactorial[a] - logFactorial[b] - logFactorial[c] - logFactorial[d]);
    }
}

