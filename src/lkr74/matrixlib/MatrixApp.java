package lkr74.matrixlib;

import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

import com.jogamp.gluegen.DebugEmitter;

import lkr74.mathgenerics.MiscMath;
import lkr74.mathgenerics.MiscMath.RandFill;
import lkr74.mathgenerics.XYLineChart_AWT;
import lkr74.matrixlib.Matrix.Type;

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

	// several matrix tests baked into one method. Supply negative index for an RC test and positive for a CSR test
	private static void matrixMultiTest(int[] testCases) {

		double[][] testRuns = new double[testCases.length][TESTRUNS];
		boolean testTrueEquality = true;	// set true to try testing equality when the matrices ARE equal
		String[] testNames = {"RC add", "RC equality", "RC multiply", "RC copy", "n/a", "CSR copy", "CSR multiply", "CSR equality", "CSR add"};
		
		int debuglevel = Matrix.DEBUG_LEVEL;
		Matrix.DEBUG_LEVEL = 0;			// no output to console during tests
		
		// the preruns are supposed to warm up the JVM
		int preruns = 500;
		for (int r = 0; r < TESTRUNS; r++)
		{
			int mWidth = r + 2, mSize = mWidth * mWidth;
			long start;
			// create matrices to copy between, we'll test copying from RC to CSR and back
			// matrix size will keep increasing from 2x2 to 50x50
			Matrix M0 = new Matrix("M0", mWidth, mWidth,  Matrix.Type.Null);
			Matrix M1 = new Matrix("M1", mWidth, mWidth,  Matrix.Type.Null);
			CSRMatrix M2 = new CSRMatrix("M2", mWidth, mWidth,  Matrix.Type.Null);
			CSRMatrix M3 = new CSRMatrix("M3", mWidth, mWidth,  Matrix.Type.Null);
			Matrix M4;
			CSRMatrix M5;
	
			RandFill rfill = new MiscMath().new RandFill(mSize);

			// fill approx. 1/8 of matrix fields with random numbers
			for (int j = 0; j < mSize/8; j++) {
				int rndPos = rfill.getRandom();
				M0.valueTo(rndPos / M0.N, rndPos % M0.M, Math.random());
				M2.valueTo(rndPos / M0.N, rndPos % M0.M, Math.random());
				if (!testTrueEquality) {
					M1.valueTo(rndPos / M0.N, rndPos % M0.M, Math.random());
					M3.valueTo(rndPos / M0.N, rndPos % M0.M, Math.random());
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
					case -1:	M4 = M0.clone(); break;							// copy from RC to CSR
					case -2:	M4 = M0.multiply(M1); break;					// RC multiply
					case -3:	result = M0.equals(M1); break;					// RC equality
					case -4:	M4 = M0.add(M1, true); break;					// RC add
					case 1:		M5 = M3.clone(); break;							// copy from CSR to RC
					case 2:		M5 = M2.multiply(M3); break;					// CSR multiply
					case 3:		M2.equals(M3); break;							// CSR equality
					case 4:		M5 = M2.add(M3, true); break;					// CSR add
					}
			    testRuns[tcnt++][r] = (System.nanoTime() - start) / ITERATIONS;
			}

		    if (--preruns > 0) r--;
		}

		Matrix.DEBUG_LEVEL = debuglevel;
		
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
		double[] d8 = {1,2,3,1,2}, d9 = {12,24,36,5,6,7}, d9b = {12,5,24,6,36,7}, v1 = {8,15,19};
		double[] d1i = {1,2,3,4,5,6,7,8};
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
						0,2,1,
						2,1,3};
		double[] gs = {	12,-5,4,
						6,16,8,
						-4,2,34};
		double[] cn = {	3,1,1,
						1,4,2,
						2,1,5};
		double[] testQR = {	4,1,3,2,
							5,4,-1,5,
							7,3,0,7,
							1,2,5,8,
							-3,-5,2,8,
							5,7,4,9};
		double[] testHH = {	 4, 1,-2, 2,
							 1, 2, 0, 1,
							-2, 0, 3,-2,
							 2, 1,-2,-1 };
		double[] testCh = {	 4,12,-16,12,37,-43,-16,-43,98};
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
		double[] d20 = { 1,-3, 0, 0, 0, 0, 0, 0, 0, 0,
						 0,-4,-3, 0, 0, 0, 0, 0, 0, 0,
						-1, 7, 1, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 2,-2, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 3, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 9, 2, 1, 0,
						 0, 0, 0, 0,-3, 0, 0, 2, 7, 0,
						 0, 0, 0, 0, 0, 0, 0, 0, 0, 8};
		double[] d11 = {1,2,3,4,5,6,7,8,
						2,2,3,4,5,6,7,8,
						3,3,3,4,5,6,7,8,
						4,4,4,4,5,6,7,8,
						5,5,5,5,5,6,7,8,
						6,6,6,6,6,6,7,8,
						7,7,7,7,7,7,7,8,
						8,8,8,8,8,8,8,8};
			
//		int[] mtests = {4, 2, -2, -4};
//		matrixMultiTest(mtests);
//		if(1==1) return;

		// Test Cholesky factorisation
		//Matrix Ch = new Matrix("C", 3, 3, testCh, null);
		Matrix Ch = new Matrix("C", 10, 10, d20, null);
		Ch = Ch.transpose(true).multiply(Ch);
		Ch.factoriseCholesky();
		if(1==1) return;
		
		// Test Householder reduction form
		Matrix HH = new Matrix("HH", 4, 4, testHH, null);
		HH = HH.reduceHouseholder(true);
		
		// Test MarixMarket file loading routine
		MatrixMarketIO mmIO = new MatrixMarketIO("data/mcca.mtx");
		Matrix MM = mmIO.getMatrix();
		System.out.println(MM.toString());
		// Test conversion of MatrixMarket data to CSR sparse format
		CSRMatrix MMcsr = CSRMatrix.convert(MM);
		System.out.println(MMcsr.toString());
		
		Matrix MM_Cholesky = MM.factoriseCholesky();
		
		if (MM_Cholesky != null) {
			MatrixBMPImage MM_image = new MatrixBMPImage(MM_Cholesky);
			MM_image.write();
		}
		
		Matrix D8 = new Matrix("D8", 8, 8, d11, null);
		D8.convergent();

		// Test addition
		Matrix D2 = new Matrix("D2", 6, 5, d2, null);
		System.out.println(D2.toString());
		D2 = D2.add(D2, true);
		// Test multiplication
		D2 = D2.multiply(D2.transpose(true));

		// Test Gram-Schmidt orthonormalisation for non-tranposed & transposed
		Matrix O = new Matrix("O", 5, 5, d2, null);
		System.out.println(O.toString());
		Matrix Q = O.orthogonalise(true);
		O.transpose(true).orthogonalise(true);
		
		// test LU = A decomposition with permutation and A = LU recomposition with depermutation
		Matrix[] lLU = Q.decomposeLU(true, true);
		Matrix LU = lLU[0].multiply(lLU[1]);
		LU.unpermute(lLU[0].mutator);
		System.out.println(LU);
		
		// Test eigenvalue & eigenvector finder findEigenQR(), which takes an orthonormalised matrix
		Matrix E = O.findEigenQR(Q, 150, 30, 0.005);

		// test QR decomposition
		O.decomposeQR(Q);
		
		// Test Power Method for finding largest eigenvalue
		D2 = new Matrix("D2", 5, 5, d2, null);
		Matrix x5 = new Matrix("x5", 5, 1, d8, null);
		System.out.println("Largest eigenvalue: " + D2.eigenPowerValue(x5, 0.01, 100) + "\n");
		
		// Test method that finds diagonal bandwidth of a sparse matrix
		System.out.println("half bandwidth: " + new Matrix("Dd", 10, 10, d20, null).getHalfBandwidth());
		
		int iters = 100;
		long tstart, tend;

		// Test partial pivoting Gauss solver
		Matrix D5 = new Matrix("A", 3, 3, d5, null);
		Matrix x6 = D5.solveGaussPP(new Matrix("c2", 3, 1, d9, null));
		// test remultiplying it
		D5.multiply(x6);
		
		Matrix N = new Matrix("N", 3, 3, cn, null);
		N.convergent();
		x6 = N.solveGaussSeidel(new Matrix("v1", 3, 1, v1, null), 100, 0.1);
		
		// Test pivoting LU decomposer and back substitution methods
		Matrix[] bLU = new Matrix("c2", 3, 1, d9, null).backSubstituteLU(D5, null, true);
		if (bLU == null) System.out.println("backSubstituteLU2() received singular matrix.");

		// Test partial pivoting Gauss Jordan with inverse matrix creation and sparse diagonal matrix optimisation
		Matrix D20 = new Matrix("D", 10, 10, d20, null);
		x6 = D20.solveGaussJordanPPDO(
				new Matrix("b", 10, 1, d10, null),
				new Matrix("D", 10, 10, Matrix.Type.Null));
		if (x6 == null) System.out.println("solveGaussJordan() returned singular matrix.");
		System.out.println("Determinant: " + D20.det);
		
		// Test full pivoting Gauss-Jordan with in-situ matrix & input vector -> solution vector transformation
		Matrix c2 = new Matrix("c2", 3, 2, d9b, null), X;
		Matrix.DEBUG_LEVEL--;
		tstart = System.nanoTime();
		for (int i = 0; i < iters; i++)
			X = D5.solveGaussJordanFP(c2, true);
		tend = System.nanoTime();
		Matrix.DEBUG_LEVEL++;
		System.out.printf("solveGaussJordanFP() averaged %.1f ns\n", (double)(tend - tstart)/iters);

		Matrix Q1 = new Matrix("Q1", 128, 128, Matrix.Type.Random), Q11;
		Matrix Q2 = new Matrix("Q2", 128, 128, Matrix.Type.Random), Q22;
		
		// Test ordinary multiply versus Strassen-Winograd multiply algorithm
//		Matrix.DEBUG_LEVEL--;
//		for (int tests = 0; tests < 5; tests++) {
//
//			tstart = System.nanoTime();
//			for (int i = 0; i < iters; i++) Q11 = Q1.multiply(Q1);
//			tend = System.nanoTime();
//
//			System.out.printf("multiply() averaged %.1f ns\n", (double)(tend - tstart)/iters);
//		    System.out.println("Mults: " + Matrix.mulFlops_DEBUG);
//
//		    tstart = System.nanoTime();
//			for (int i = 0; i < iters; i++) Q22 = Matrix.multiplyStrasWin(Q2, Q2, 8);
//			tend = System.nanoTime();
//
//			System.out.printf("multiplyStrasWin() averaged %.1f ns\n", (double)(tend - tstart)/iters);
//		    System.out.println("Mults: " + Matrix.mulFlopsSW_DEBUG);
//		    System.out.println("Adds: " + Matrix.mulAdopsSW_DEBUG);
//		    System.out.println("Recurses: " + Matrix.mulSW_DEBUG);
//		}
//		Matrix.DEBUG_LEVEL++;
		
		// Test if multiplication of rescaled matrices produce same results as unexpanded ones
		Matrix S1 = new Matrix("S1", 6, 5, d2, null);
		S1 = S1.multiply(S1.transpose(true));
		S1 = new Matrix("S1", 6, 5, d2, null);
		S1 = S1.rescale(0, 0, 6, 6, true);
		S1 = S1.multiply(S1.transpose(true));
		
		// Test value x matrix multiplication
		Matrix A1 = new Matrix("A1", 9, 9, d3, null);
		A1 = A1.multiply(-2000, false);

		Matrix G = new CSRMatrix("G", 9, 9, d3, null);

		// Test Laplacian determinant finding method
//		Matrix.DEBUG_LEVEL--;
//		for (int tests = 0; tests < 3; tests++) {
//			
//			tstart = System.nanoTime();
//			for (int i = 0; i < iters; i++) Matrix.determinantLaplace(G, 2);
//			tend = System.nanoTime();
//			
//		    System.out.printf("determinantLaplaceR2() averaged %.1f ns\n", (double)(tend - tstart)/iters);
//			System.out.println("determinant: " + G.det);
//			System.out.println("recursions: " + Matrix.detL_DEBUG + "\n");
//		}
//		Matrix.DEBUG_LEVEL++;

		// Test CSR matrices, test vAv style multiplication
		CSRMatrix x1 = new CSRMatrix("x1", 1, 9, d, null);
		CSRMatrix x2 = new CSRMatrix("x2", 7, 1, d, null);
		CSRMatrix V = new CSRMatrix("V", 9, 7, d4, null);
		x2.swap(6, 2);
		CSRMatrix S = x1.multiply(V);
		S = S.multiply(x2);

		// Test CSRMatrix conversion and polymorphic access
		CSRMatrix O2 = new CSRMatrix("O", 9, 7, d3, null);
		O2.transpose(false);
		G = O2.multiply(new CSRMatrix("R", 9, 7, d3, null));
		G = G.eliminateRowColumn(6, 6, true);

		// test centering method for Matrix
		Matrix P = G.center(true);

		BinBitImage.compact(d2, O2.bitImage.data[0], 0);
		System.out.println(BinBitImage.binBitToString());
		CSRMatrix I = new CSRMatrix("I", O2.M, O2.N, O2.getDataRef()[0], O.getDataRef()[1]);
		System.out.println(O2.name + " equals " + I.name + ": " + O2.equals(I));
		O2.valueTo(2, 3, 0.777);
		System.out.println(O2.name + " equals " + I.name + ": " + O2.equals(I));
		
		System.out.println("Create CSR matrix, read particular values:");
		CSRMatrix csrO = new CSRMatrix("O2", O2.M, O2.N, O2.getDataRef()[0], O2.getDataRef()[1]);
		System.out.println("(3,1): " + csrO.valueOf(3, 1));
		System.out.println("(1,0): " + csrO.valueOf(1, 0));
		
		Matrix D = new Matrix("D", 3, 3, d, null);

		D.doGaussEliminationPP();
		CSRMatrix csrA = new CSRMatrix("A2", 5, 5,  Matrix.Type.Random);

		CSRMatrix csrS = new CSRMatrix("S", 5, 5, Matrix.Type.Null);
		CSRMatrix csrT = new CSRMatrix("T", 5, 5, Matrix.Type.Identity);
		csrT.valueTo(1, 3, 9.9);
		csrS.valueTo(2, 2, 9.9);
		csrS.valueTo(0, 4, -8);
		CSRMatrix csrU = csrS.add(csrT, true);
		CSRMatrix csrD = csrU.multiply(csrT);
		CSRMatrix B = csrD.clone();
		
		Matrix C = new Matrix("C", 5, 5, Matrix.Type.Identity);
		System.out.println("Identity matrix\n" + C.toString());
		System.out.println();

		System.out.println("Multiplied matrices polymorphically\n" + csrA.multiply(B));
		System.out.println();

		// shouldn't be equal since AB != BA in general
		csrS = new CSRMatrix("S", 5, 5, d2, null);
		csrT = new CSRMatrix("T", 5, 5, d2, null);
		csrT = csrT.add(5, false);
		System.out.println(csrS.multiply(csrT).equals(csrT.multiply(csrS)));
		System.out.println();

		Matrix b = new Matrix("b", 5, 1, Matrix.Type.Random);

		Matrix A = new Matrix("A3", csrA.M, csrA.N, Matrix.Type.Null);
		A = csrA.clone();
		Matrix x = A.solveGaussPP(b);
		A = A.multiply(x);
		
	}

}
