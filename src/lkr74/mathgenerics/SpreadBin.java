package lkr74.mathgenerics;

import java.security.InvalidParameterException;

public class SpreadBin {
	
	// Class for constructing and outputting single-value quantifiable distributions
	// Leonard Krylov 2017

	private int bins = 0, barScale = 100, illegalCases = 0;
	public double[] div = null;
	public long[] bin = null;
	
	public SpreadBin(int n, double min, double max, int barScale) {		// instantiates distribution over n bins
		if (min >= max) throw new InvalidParameterException("SpreadBin.SpreadBin(): invalid min/max parameters");
		if (n <= 0) throw new InvalidParameterException("SpreadBin.SpreadBin(): invalid bin count");
		if (barScale <= 0) throw new InvalidParameterException("SpreadBin.SpreadBin(): invalid bar scale");
		bins = n;
		this.barScale = barScale;
		div = new double[bins + 1];
		bin = new long[bins];
		double dlt = (max - min) / (double)bins;
		div[0] = min;
		for (int b = 1; b <= bins; b++) div[b] = min + dlt * b;
	}
	
	public void add(double v) {
		if (v < div[0] || v > div[bins]) {
			illegalCases++; return; }						// illegal case: value was outside the specified min - max range
		int b = 1;
		while (b <= bins && v > div[b]) b++;
		bin[--b]++;
	}
	
	public int maxBin() {									// gets largest bin of all data clusters
		long vMax = 0;
		int maxBin = -1;
		for (int b = 0; b < bins; b++) if (bin[b] > vMax) { vMax = bin[b]; maxBin = b; }
		return maxBin;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (illegalCases > 0) System.out.println("Illegal cases!: " + illegalCases);
		sb.append("Distribution across " + bins + " samples:\n");
		long vMax = 0;
		for (int b = 0; b < bins; b++) if (bin[b] > vMax) vMax = bin[b];
		if (vMax == 0) return sb.toString();
		for (int b = 0; b < bins; b++) {					// generate horizontal bars
			int barL = (int)((barScale * bin[b]) / vMax);
			sb.append(String.format("%-12d%-12.4f", bin[b], div[b+1]));
			for (int bar = 0; bar < barL; bar++) sb.append("|");
			sb.append("\n");
		}
		return sb.toString();
	}
}
