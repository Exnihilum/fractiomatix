package lkr74.matrixlib;

import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

import lkr74.mathgenerics.XYLineChart_AWT;

public class MatrixApp {
	
	
	// create Chart data from timings of 2 compared tested functions, showing variation in one parameter (ex: size/fillrate)
	static XYDataset createStatisticSet(double[][] timingLists, String[] testNames, int[] testCases)
	{
		final XYSeriesCollection dataset = new XYSeriesCollection();          
		for (int j = 0; j < timingLists.length; j++) {
			final XYSeries testSeries = new XYSeries(testNames[testCases[j]+4], false);

			for (int i = 0; i < timingLists[j].length; i++)
				testSeries.add(i+2, timingLists[j][i]);
	
			dataset.addSeries(testSeries);          
		}
		return dataset;
	}
	
	final static int TESTRUNS = 200, ITERATIONS = 400;

	// several matrix tests baked into one method. Supply negative index for an RC test and positive for a CRS test
	private static void matrixMultiTest(int[] testCases) {

		double[][] testRuns = new double[testCases.length][TESTRUNS];
		boolean testTrueEquality = true;	// set true to try testing equality when the matrices ARE equal
		String[] testNames = {"RC add", "RC equality", "RC multiply", "RC copy", "n/a", "CRS copy", "CRS multiply", "CRS equality", "CRS add"};

		// the preruns are supposed to warm up the JVM
		int preruns = 500;
		for (int r = 0; r < TESTRUNS; r++)
		{
			long start;
			// create matrices to copy between, we'll test copying from RC to CRS and back
			// matrix size will keep increasing from 2x2 to 50x50
			Matrix M0 = new Matrix("M0", r+2, r+2,  Matrix.Type.Null);
			Matrix M1 = new Matrix("M1", r+2, r+2,  Matrix.Type.Null);
			Matrix M2 = new CSRMatrix("M2", r+2, r+2,  Matrix.Type.Null);
			Matrix M3 = new CSRMatrix("M3", r+2, r+2,  Matrix.Type.Null);
			
			// fill approx. 1/8 of matrix fields with random numbers
			for (int j = 0; j < ((r+2)*(r+2))/8; j++) {	
				M0.valueTo(((int)(Math.random()*100)) % (r+1), ((int)(Math.random()*100)) % (r+1), Math.random());
				M2.valueTo(((int)(Math.random()*100)) % (r+1), ((int)(Math.random()*100)) % (r+1), Math.random());
				if (!testTrueEquality) {
					M1.valueTo(((int)(Math.random()*100)) % (r+1), ((int)(Math.random()*100)) % (r+1), Math.random());
					M3.valueTo(((int)(Math.random()*100)) % (r+1), ((int)(Math.random()*100)) % (r+1), Math.random());
				}
			}
			if (testTrueEquality) { M1 = M0.clone(); M3 = M2.clone(); }

			// iterate through test cases
			int tcnt = 0;
			boolean result;
			for (int testtype : testCases) {
			    start = System.nanoTime();
				for (int i = 0; i < ITERATIONS; i++)
					switch (testtype) {
					case -1:	M2 = M0.clone(); break;							// copy from RC to CRS
					case -2:	Matrix.multiply(M0, M1); break;					// RC multiply
					case -3:	result = M0.equals(M1); break;					// RC equality
					case -4:	Matrix.add(M0, M1); break;						// RC add
					case 1:		M0 = M3.clone(); break;							// copy from CRS to RC
					case 2:		CSRMatrix.multiply(M2, M3); break;				// CRS multiply
					case 3:		M2.equals(M3); break;							// CRS equality
					case 4:		CSRMatrix.add(M2, M3); break;					// CRS add
					}
			    testRuns[tcnt++][r] = (System.nanoTime() - start) / ITERATIONS;
			}

		    if (--preruns > 0) r--;
		}

		// chart for comparison timing runs
		XYDataset bestFitChartSet;
		bestFitChartSet = createStatisticSet(testRuns, testNames, testCases);
		XYLineChart_AWT bestFitChart = new XYLineChart_AWT("Matrix stress test", "matrix size", "nanosecs/matrix op.", bestFitChartSet, 1024, 768);
		bestFitChart.pack( );          
		RefineryUtilities.centerFrameOnScreen(bestFitChart);   
		bestFitChart.setVisible( true ); 
	}
	
	
	// test client
	public static void main(String[] args) {
				
		double[] d = {1, 0, 3, 0, 5, 0, 0, 1, 3};
		double[] d8 = {1,2,3,1,2}, d9 = {12,24,36};
		double[] d2 = {	5,1,2,0,4,
						1,4,2,1,3,
						2,2,5,4,0,
						0,1,4,1,3,
						4,3,0,3,4,
						0,0,0,0,0};
		double[] d3 = {	0,3,0,7,0,6,0,6,3,
						0,4,8,0,0,2,0,7,7,
						9,0,3,9,0,4,0,9,0,
						0,0,0,0,0,0,7,3,1,
						0,1,3,0,9,9,0,0,0,
						0,4,0,2,1,0,0,9,9,
						3,0,3,4,0,3,0,1,3,
						0,9,0,0,0,0,3,0,0,
						0,6,0,9,0,6,7,7,6};
		double[] d4 = {	0,2,0,0,0,6,0,
						5,0,8,0,0,0,0,
						1,0,3,0,0,0,0,
						0,0,3,0,0,0,7,
						0,0,3,0,0,0,0,
						0,0,3,0,0,0,0,
						0,0,3,0,0,0,0,
						0,0,3,0,0,0,3,
						0,6,0,0,0,0,0};
		double[] d5 = {	1,2,3,
						3,2,1,
						2,1,3};
		double[] d10 = {1,2,3,4,5,6,7,8,7,6,
						2,2,3,4,5,6,7,8,7,6,
						3,3,3,4,5,6,7,8,7,6,
						4,4,4,4,5,6,7,8,7,6,
						5,5,5,5,5,6,7,8,7,6,
						6,6,6,6,6,6,7,8,7,6,
						7,7,7,7,7,7,7,8,7,6,
						8,8,8,8,8,8,8,8,7,6,
						7,7,7,7,7,7,7,7,7,6,
						6,6,6,6,6,6,6,6,6,6};
		double[] d20 = {1,0,0,0,0,0,0,0,0,0,
						0,0,0,0,0,0,0,0,0,0,
						0,0,1,0,0,0,0,0,0,0,
						0,0,0,2,0,0,0,0,0,0,
						0,0,0,0,3,0,0,0,0,0,
						0,0,0,2,0,0,0,0,0,0,
						0,0,0,0,0,0,1,0,0,0,
						0,0,0,0,0,0,0,2,1,0,
						0,0,0,0,0,0,0,0,7,0,
						0,0,0,0,0,0,0,0,0,8};
		double[] d11 = {1,2,3,4,5,6,7,8,
						2,2,3,4,5,6,7,8,
						3,3,3,4,5,6,7,8,
						4,4,4,4,5,6,7,8,
						5,5,5,5,5,6,7,8,
						6,6,6,6,6,6,7,8,
						7,7,7,7,7,7,7,8,
						8,8,8,8,8,8,8,8};
			
//		int[] tests = {-1, 2, -2, -4};
//		matrixMultiTest(tests);
//		if(1==1) return;

//		Matrix D8 = new Matrix("D8", 8, 8, d11);
//		D8.convergent();
//		Matrix d7 = new Matrix("d8", 8, 1, d);
//		D8.conditionDiagonal(d7, false, false);		// test method with add method and don't create a bitimage
//		D8.convergent();
//		System.out.println(d7.toString());

//		
//		Matrix D2 = new Matrix("D2", 5, 5, d2);
//		Matrix x5 = new Matrix("x5", 5, 1, d8);
//		double eigenvalue = Matrix.eigenPowerMethod(D2, x5, 0.01, 100);
//		System.out.println("Eigenvalue: " + eigenvalue + "\n");
		
		Matrix D20 = new Matrix("Dd", 10, 10, d20);
		System.out.println(D20.diagonality());
		
		int iters = 1;
		long tstart, tend;

		Matrix.DEBUG_LEVEL = 2;

		Matrix D5i = new Matrix("D5i", 3, 3);
		Matrix c2 = new Matrix("c2", 3, 1, d9), x6 = null;
		Matrix[] XA = null;
		Matrix D5 = new Matrix("A", 3, 3, d5);
		Matrix[] UVl = D5.factorise();
		x6 = D5.solve(c2);
		if (UVl != null)	x6 = c2.solveCrout(UVl[0], UVl[1]);
		else				System.out.println("Crout solver: nonfactorisable matrix.");

		x6 = D5.solveGaussJordan(c2, D5i);
		
		tstart = System.nanoTime();
		System.out.println(D5.toString());
		for (int i = 0; i < iters; i++)
			XA = D5.solveGaussJordan2(c2, false);
		tend = System.nanoTime();
		System.out.printf("solveGaussJordan2() averaged %.1f ns\n", (double)(tend - tstart)/iters);
		if(1==1) return;

		Matrix Q1 = new Matrix("Q1", 128, 128, Matrix.Type.Random), Q11;
		Matrix Q2 = new Matrix("Q2", 128, 128, Matrix.Type.Random), Q22;
			
		for (int tests = 0; tests < 5; tests++) {
			tstart = System.nanoTime();
			for (int i = 0; i < iters; i++)
				Q11 = Matrix.multiply(Q1, Q1);
			tend = System.nanoTime();
		    System.out.printf("multiply() averaged %.1f ns\n", (double)(tend - tstart)/iters);
		    System.out.println("Mults: " + Matrix.mulFlops_DEBUG);
			tstart = System.nanoTime();
			for (int i = 0; i < iters; i++)
				Q22 = Matrix.multiplyStrasWin(Q2, Q2, 8);
			tend = System.nanoTime();
		    System.out.printf("multiplyStrasWin() averaged %.1f ns\n", (double)(tend - tstart)/iters);
		    System.out.println("Mults: " + Matrix.mulFlopsSW_DEBUG);
		    System.out.println("Adds: " + Matrix.mulAdopsSW_DEBUG);
		    System.out.println("Recurses: " + Matrix.mulSW_DEBUG);
		}
		Matrix.DEBUG_LEVEL = 2;
		if(1==1) return;
		
		// test if multiplication of rescaled matrices produce same results as unexpanded ones
		Matrix S1 = new Matrix("S1", 6, 5, d2);
		S1 = Matrix.multiply(S1, S1.transpose(true));
		S1 = new Matrix("S1", 6, 5, d2);
		S1 = S1.rescale(0, 0, 6, 6, true);
		S1 = Matrix.multiply(S1, S1.transpose(true));
		if(1==1) return;
		
		Matrix A1 = new Matrix("A1", 9, 9, d3);
		A1 = Matrix.multiply(2, A1);
		Matrix A2 = new Matrix("A2", 9, 7, d4);
		//A2.transpose();
		Matrix A3 = Matrix.multiply(A1, A2);

		Matrix G = new CSRMatrix("G", 9, 9, d3);

		Matrix.DEBUG_LEVEL = 0;
		for (int tests = 0; tests < 3; tests++) {
			G.makeNonThreaded();
			tstart = System.nanoTime();
			for (int i = 0; i < iters; i++)
				Matrix.determinantLaplace(G, 2);
			tend = System.nanoTime();
		    System.out.printf("determinantLaplaceR2() singlethread averaged %.1f ns\n", (double)(tend - tstart)/iters);
			System.out.println("determinant: " + G.det);
			System.out.println("recursions: " + Matrix.detL_DEBUG + "\n");
			G.makeThreaded();
			tstart = System.nanoTime();
			for (int i = 0; i < iters; i++)
				Matrix.determinantLaplace(G, 3);
			tend = System.nanoTime();
		    System.out.printf("determinantLaplaceR3() multithread averaged of %.1f ns\n", (double)(tend - tstart)/iters);
			System.out.println("determinant: " + G.det);
			System.out.println("recursions: " + Matrix.detL_DEBUG + "\n");
		}
		Matrix.DEBUG_LEVEL = 1;

		CSRMatrix x1 = new CSRMatrix("x1", 1, 9, d);
		CSRMatrix x2 = new CSRMatrix("x2", 7, 1, d);
		CSRMatrix V = new CSRMatrix("V", 9, 7, d4);
		x2.swap(6, 2);
		CSRMatrix S = CSRMatrix.multiply(x1, V);
		S = CSRMatrix.multiply(S, x2);
		if(1==1) return;

		// test out CRSMatrix conversion and access
		Matrix O = new Matrix("O", 9, 7, d3);
		O.transpose(false);
		Matrix R = new Matrix("R", 9, 7, d3);
		//G = Matrix.multiply(O, R);
		G = G.eliminateRowColumn(6, 6, true);

		// test centering method for Matrix
		Matrix P = Matrix.center(G);

		BinBitImage.compact(d2, O.bitImage.data[0], 0);
		System.out.println(BinBitImage.binBitToString());
		CSRMatrix I = new CSRMatrix("I", O.M, O.N, O.getDataRef());
		System.out.println("O equals I: " + O.equals(I));
		O.valueTo(2, 3, 0.777);
		System.out.println("O equals I: " + O.equals(I));
		
		System.out.println("Create CSR matrix, read particular values:");
		CSRMatrix csrO = new CSRMatrix("O2", O.M, O.N, O.getDataRef());
		System.out.println("(3,1): " + csrO.valueOf(3, 1));
		System.out.println("(1,0): " + csrO.valueOf(1, 0));
		
		Matrix D = new Matrix("D", 3, 3, d);

		D.doGaussElimination();
		CSRMatrix csrA = new CSRMatrix("A2", 5, 5,  Matrix.Type.Random);
		
		System.out.println("Create CSR matrix, copy to RC matrix and transpose:");
		Matrix B = new Matrix("B", csrA.M, csrA.N, Matrix.Type.Null);
		B = csrA.clone();
		B.transpose(false);

		CSRMatrix csrS = new CSRMatrix("S", 5, 5, Matrix.Type.Null);
		CSRMatrix csrT = new CSRMatrix("T", 5, 5, Matrix.Type.Identity);
		csrT.valueTo(1, 3, 9.9);
		csrS.valueTo(2, 2, 9.9);
		csrS.valueTo(0, 4, -8);
		B = csrT.clone();
		CSRMatrix csrU = CSRMatrix.add(csrS, csrT);
		B = csrU.clone();
		CSRMatrix csrD = CSRMatrix.multiply(csrU, csrT);
		B = csrD.clone();
		
		Matrix C = new Matrix("C", 5, 5, Matrix.Type.Identity);
		System.out.println("Identity matrix\n" + C.toString());
		System.out.println();

		//System.out.println("Multiplied matrices\n" + Matrix.multiply(csrA, B).toString());
		System.out.println();

		// shouldn't be equal since AB != BA in general
		//System.out.println(Matrix.multiply(csrA, B).equals(Matrix.multiply(B, csrA)));
		System.out.println();

		Matrix b = new Matrix("b", 5, 1, Matrix.Type.Random);

		Matrix A = new Matrix("A3", csrA.M, csrA.N, Matrix.Type.Null);
		A = csrA.clone();
		Matrix x = A.solve(b);

		Matrix.multiply(A, x);
		
	}

}
