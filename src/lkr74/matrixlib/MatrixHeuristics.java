package lkr74.matrixlib;

import lkr74.mathgenerics.RandFill;

public class MatrixHeuristics {
	
	// TODO: the Matrix heuristics class is meant to do some automated microbenchmark tests to determine
	// certain boundaries of matrixlib types utilisation and for picking algorithms strategically

	public static final int USE_RC = 0, USE_CSR = 1, USE_NSP = 2;
	public static int[][] multiplyAlgoPrio;
	
	private double[] getRandomData(int r, int c, int fillrate) {
		double[] data = new double[r * c];
		RandFill rfill = new RandFill(r * c);
		
		for (int fill = 0; fill < (r * c) * ((float)fillrate / 100.0); fillrate++)
			data[rfill.getRandom()] = Math.random();
		return data;
	}
	
	// TODO: test multiplying speeds of the three matrix types on different dimensions and fill rates
	public void measureMultiplySpeed(int runs) {
		
		for (int scale = 4; scale <= 60; scale++) {							// test different matrix scales

			for (int fillrate = 10; fillrate <= 100; fillrate++) {			// test different matrix fillrates

				double[] data = getRandomData(scale, scale, fillrate);
				Matrix CR = new Matrix("CR", scale, scale, data, null);
				CSRMatrix CSR = new CSRMatrix("CSR", scale, scale, data, null);
				NSPMatrix NSP = new NSPMatrix("NSP", scale, scale, data, null);
				data = getRandomData(scale, scale, fillrate);
				Matrix CR2 = new Matrix("CR2", scale, scale, data, null);
				CSRMatrix CSR2 = new CSRMatrix("CSR2", scale, scale, data, null);
				NSPMatrix NSP2 = new NSPMatrix("NSP2", scale, scale, data, null);
				
				long tStart, tDelta;
				for (int tests = 0; tests < runs; tests++) {
					tStart = System.nanoTime();
					Matrix CR3 = CR.multiply(CR2);
					tDelta = System.nanoTime() - tStart;
				}
			}

		}
	}
	
}
