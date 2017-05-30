package lkr74.mathgenerics;

public class VisitBitArray {
	
	// Class for fast & compact flagging/status changes
	// Leonard Krylov 2017
	
	private long[] array = null;
	private int[] resetStack = null;
	private int bitSets = 0, resets = 0, maxResets = 6;
	private static int criterionResets = -1;

	public VisitBitArray(int items) {
		array = new long[bitSets = (items >> 6) + 1];
		// TODO: microbenchmarking for maxResets should be done on repeated basis since allocation becomes sluggish during low memory periods,
		// 6 longs is fine for low-congestion environment
		if (criterionResets > 0) maxResets = criterionResets;
		if (maxResets > 0) { resetStack = new int[maxResets + 1]; resetStack[0] = -1; }
	}

	private VisitBitArray(int items, int maxResets) {
		array = new long[bitSets = (items >> 6) + 1];
		this.maxResets = maxResets;
		resetStack = new int[maxResets + 1];
		resetStack[0] = -1;
	}

	public void clearVisits() {
		if (resets >= maxResets)		{ array = new long[bitSets];  resets = 0; }
		else if (resets > 0)			{
			for (int i = 1; i <= resets; i++) array[resetStack[i]] = 0; resets = 0; resetStack[0] = -1; }
	}
	
	public boolean visited(int i) { return (array[i >> 6] & (0x1L << (i & 63))) != 0; }
	
	public void visit(int i) {
		int iD64 = i >> 6;
		array[iD64] |= (0x1L << (i & 63));
		// simple heuristic for eliminating duplicate bitfield resettings: check if it's same field as the previous
		if (resets < maxResets && resetStack[resets] != iD64)
			resetStack[++resets] = iD64;
	}
	
	public static int criterionmaxResets(boolean verbose) {
		int optimumMaxReset = 0;
		for (int maxR = 1; maxR < 4000; maxR += 1) {
			//if (maxR % 100 == 0) System.gc();
			//System.gc();
			long tTime = 0, tTimeR = 0;
			for (int j = 0; j < 5000; j++) {
				VisitBitArray vba = new VisitBitArray(800000, 0);
				for (int i = 0; i < maxR; i++) vba.visit((int)(Math.random()*799999));
				long tStart = System.nanoTime();
				vba.clearVisits();
				if (j >= 1000) tTime += System.nanoTime() - tStart;
			}
			if (verbose) System.out.println("vba: " + tTime/4000.0);
			//System.gc();
			for (int j = 0; j < 5000; j++) {
				VisitBitArray vbaR = new VisitBitArray(800000, maxR);
				for (int i = 0; i < maxR; i++) vbaR.visit((int)(Math.random()*799999));
				long tStartR = System.nanoTime();
				vbaR.clearVisits();
				if (j >= 1000) tTimeR += System.nanoTime() - tStartR;
			}
			if (verbose) System.out.println("vbaR: " + tTimeR/4000.0);
			if (tTimeR > tTime) { optimumMaxReset = maxR; break; }
		}
		if (verbose) System.out.println("optimum: " + optimumMaxReset);
		return criterionResets = optimumMaxReset;
	}

}
