package lkr74.matrixlib;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;
import lkr74.fem1.FEM1;
import lkr74.fem1.FEM1Element;
import lkr74.fem1.FEM1Octant;
import lkr74.fem1.FEM1Octree;
import lkr74.mathgenerics.RandFill;
import lkr74.mathgenerics.XYLineChart_AWT;
import lkr74.sizeof.SizeOf;

public class MatrixApp {
	
	final static boolean COPY = true, NO_COPY = false;
	
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
			Matrix M0 = new Matrix("M0", mWidth, mWidth,  Matrix.Type.Null, 1);
			Matrix M1 = new Matrix("M1", mWidth, mWidth,  Matrix.Type.Null, 1);
			CSRMatrix M2 = new CSRMatrix("M2", mWidth, mWidth,  Matrix.Type.Null, 1);
			CSRMatrix M3 = new CSRMatrix("M3", mWidth, mWidth,  Matrix.Type.Null, 1);
			Matrix M4;
			CSRMatrix M5;
	
			RandFill rfill = new RandFill(mSize);

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
					case -4:	M4 = M0.add(M1, COPY); break;					// RC add
					case 1:		M5 = M3.clone(); break;							// copy from CSR to RC
					case 2:		M5 = M2.multiply(M3); break;					// CSR multiply
					case 3:		M2.equals(M3); break;							// CSR equality
					case 4:		M5 = M2.add(M3, COPY); break;					// CSR add
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
		bestFitChart.pack();          
		RefineryUtilities.centerFrameOnScreen(bestFitChart);   
		bestFitChart.setVisible( true ); 
	}
	
	
	
	
	public static void testLUdecomposure(String fileName, boolean toImage, boolean toFile, boolean toGraphViz, int benchRuns) {
		
		double[] vB = {
			0,0,0,1,2,3,4,0,0,0,6,5,4,3,2,1,0,0,0,0,8,0,0,0,0,0,0,0,0,5,4,3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,4,3,0,0,0,0,0,9,1,1,1,1,0,9,4,4,4,5,5,5,1,1,1,0,0,9,8,7,6,5,4,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,2,3,0,0,0,0,0,1,1,0,0,0,0,0,0,9,9,9,0,0,0,0,0,0,0,0,0,0,6,5,4,0,0,0,0,0,0,2,0,3,0,0,0,0,0,0,0,0,0,0,3,4,5,6,7,8,0,9,4,4,4,5,5,5,1,1,1,0,0,9,8,7,6,5,4,0,0,0,0,0,0,0};
		double[] d20 = {  1,-3, 0, 0, 0, 0, 0, 0, 0,
				 		  0,-4, 0, 7, 0, 0, 0, 1, 0,
				 		 -1, 1, 1, 1, 0, 0, 0, 1, 0,
				 		  0, 0, 0, 2, 0, 5, 0, 0, 0,
				 		  0, 0, 3, 1, 3, 1, 0, 1, 0,
				 		  0, 0, 0, 0, 0, 1, 0, 7, 0,
				 		  0, 0, 0, 0,-1, 1, 2, 1, 0,
				 		  0, 0, 0, 0, 0, 0, 0,-1, 4,
				 		  0, 0, 0, 0, 0, 0, 3, 7, 3};
		
		// if filename supplied, load in that matrix of MatrixMarket datatype
		Matrix MM = null;
		if (fileName != "") {
			MatrixMarketIO mmIO = new MatrixMarketIO(fileName, 0, true);
			MM = mmIO.getMatrix();
			//MM.toFile(0);
		} else {
			MM = new Matrix("MM", 9, 9, d20, null);
			if (benchRuns < 2) System.out.println(MM.toString());
		}
		MM.analyse(true);

		Matrix[] LU = null;
		long tStart = 0, tEnd = 0, preRuns = benchRuns/100;
		
		if (benchRuns > 1) Matrix.DEBUG_LEVEL--;			// print matrix debug data if we're not benchmarking
		
		for (int r = 0; r < benchRuns; r++) {
			if (r == preRuns) tStart = System.nanoTime();
			LU = MM.decomposeLU(COPY, false, true);			// decompose to a copy, don't split into L & U, use row permutation	
		}
		tEnd = System.nanoTime();
		
		if (benchRuns > 1)  Matrix.DEBUG_LEVEL++;
		
		System.out.printf("Matrix.decomposeLU() averaged %.1f ns\n", (double)(tEnd - tStart)/(benchRuns - preRuns));

		if (LU[0] != null) {
			LU[0].name = "LU";
			if (benchRuns < 2) System.out.println(LU[0].toString());
			LU[0].analyse(false);
			System.out.println(Arrays.toString(LU[0].mutator));
			if (toFile) LU[0].toFile(0);
			if (toImage) {
				MatrixBMPImage MM_image = new MatrixBMPImage(LU[0]);
				MM_image.write();
			}
		}
		Matrix vectorB = new Matrix("b", fileName != "" ? MM.M : 9, 1, vB, null);
		Matrix bLU[] = vectorB.backSubstituteLU(null, LU[0], COPY);
		bLU[0].name = "x";
		if (toFile) bLU[0].toFile(5);

		// try loading, toString, LU decomposing and writing to BMP of a NSPMatix type
		NSPMatrix MM2 = null;
		if (fileName != "") {
			MatrixMarketIO mmIO = new MatrixMarketIO(fileName, 2, true);
			MM2 = mmIO.getNSPMatrix();				
		} else {
			MM2 = new NSPMatrix("MM", 9, 9, d20, null);
			MM2.name = "MM_NSP";
			if (benchRuns < 2) System.out.println(MM2.toString());
			if (toGraphViz) MM2.toGraphViz(true, true, true);
		}
		System.out.println("MM2 non-zeroes: " + MM2.nNZ + ", percent: " + (100f * (float)MM2.nNZ / (MM2.M*MM2.N)) + "%");
		
		NSPMatrix[] LU2 = null;
		
		if (benchRuns > 1) Matrix.DEBUG_LEVEL--;
		
		for (int r = 0; r < benchRuns; r++) {
			if (r == preRuns) tStart = System.nanoTime();
			LU2 = MM2.decomposeLU(COPY, false);				// decompose to a copy, don't split into L & U
		}
		tEnd = System.nanoTime();
		
		if (benchRuns > 1) Matrix.DEBUG_LEVEL++;
		
		System.out.printf("NSPMatrix.decomposeLU() averaged %.1f ns\n", (double)(tEnd - tStart)/(benchRuns - preRuns));

		if (LU2[0] != null) {
			LU2[0].name = "LU_NSP";
			System.out.println("LU non-zeroes: " + LU2[0].nNZ + ", percent: " + (100f * (float)LU2[0].nNZ / (LU2[0].M*LU2[0].N)) + "%");
			System.out.println(Arrays.toString(LU2[0].mutator));
			if (toFile) LU2[0].toFile(0);
			if (toGraphViz) LU2[0].toGraphViz(true, true, true);
			if (benchRuns < 2) System.out.println(LU2[0].toString());
			if (toImage) {
				NSPMatrixBMPImage MM_image = new NSPMatrixBMPImage(LU2[0]);
				MM_image.write();
			}
		}
		NSPMatrix vectorB2 = new NSPMatrix("b", fileName != "" ? MM.M : 9, 1, vB, null);
		NSPMatrix bLU2[] = vectorB2.backSubstituteLU(null, LU2[0], true);
		bLU2[0].name = "x_NSP";
		if (toFile) bLU2[0].toFile(5);
	}
	
	
	
	
	// test NspNode finder with three search algorithms (iterative, linear & binary) at different heuristic switch levels	
	public static int[] testFindHVSpNode(int arraySize, int iters, int filterIter, boolean verbose) {
		
		long tstart, tend;
		double[][] testRuns = new double[3][arraySize];
		
		int offs = (int)(Math.random()*arraySize) / 4, sizeStep = (int)Math.sqrt(arraySize);
		int itrLim = NSPMatrix.nodeSearch_iterateLim, linLim = NSPMatrix.nodeSearch_linearLim;
		
		NSPNode[] nodes = new NSPNode[arraySize];
		//for (int i = 0, r = offs; i < arraySize; i++, r += 1 + (int)(Math.random()*31))
		for (int i = 0, rstep = offs; i < arraySize; i++, rstep += 1 + (int)(Math.random() * sizeStep))
			nodes[i] = new NSPNode(0, rstep, i, i);
		
		// prerun
		for (int k = 0; k < iters * 5; k++) {
			int csought = (int)(Math.random()*(nodes[arraySize-1].c - nodes[0].c));
			int found = NSPMatrix.findHVspNode(nodes, 0, arraySize-1, -1, csought);
			found += csought; csought += found;
		}
		
		if (verbose) System.out.println("NSPMatrix.findHVspNode() search test for increasing seeking lengths in a randomly ascending array");
		//RandFill rFill = new RandFill(sizeStep);
		int[] randSeek = new int[arraySize];
		
		for (int sLen = 1; sLen < arraySize; sLen++) {
			
			// rerandomise node array on every test run
			offs = (int)(Math.random()*arraySize);								// add random amount of zeroes at start of dataset
			for (int i = 0, rstep = offs; i < sLen; i++) {
			//for (int i = 0, rstep = offs; i < sLen; i++, rstep++)
				randSeek[i] = nodes[i].c = rstep;								// memorise what column values we randomly generated
				rstep += 1 + (int)(Math.random() * sizeStep);
			}
			for (int i = 0; i < sLen; i++) {									// switch around sought column values randomly
				int r = (int)(Math.random() * sLen); int tmp = randSeek[i]; randSeek[i] = randSeek[r]; randSeek[r] = tmp; }
			
			//NSPMatrix.nodeSearch_iterateLim = arraySize + 1; NSPMatrix.nodeSearch_linearLim = arraySize + 1;
			NSPMatrix.nodeSearch_iterateLim = Integer.MAX_VALUE; NSPMatrix.nodeSearch_linearLim = Integer.MAX_VALUE;
			System.gc();
			tstart = System.nanoTime();
			for (int k = 0, s = 0; k < iters; k++) {
				//int csought = (int)(Math.random()*(nodes[sLen - 1].c - nodes[0].c));
				//int csought = (int)(Math.random() * cMax);
				int found = NSPMatrix.findHVspNode(nodes, 0, sLen - 1, -1, randSeek[s++]);
				if (s >= sLen) s = 0;
				// show the finding, or the return of what was the nearest element
				//int cfound = (found < 0 ? nodes[-found-1].c : nodes[found].c);
				//System.out.println("sought: " + csought + (found < 0 ? ", nearest: " : ", found: ") + cfound);
				found += found;													// fake variable usage to avoid compiler optimisations
			}
			tend = System.nanoTime();
			testRuns[0][sLen-1] = (double)(tend - tstart)/iters;
			if (verbose && (sLen % sizeStep == 0 || sLen < sizeStep)) System.out.printf("%d elems, iter: %5.1f ", sLen, (double)(tend - tstart)/iters);
			
			NSPMatrix.nodeSearch_iterateLim = 0;
			System.gc();
			tstart = System.nanoTime();
			for (int k = 0, s = 0; k < iters; k++) {
				//int csought = (int)(Math.random()*(nodes[sLen - 1].c - nodes[0].c));
				//int csought = (int)(Math.random() * cMax);
				int found = NSPMatrix.findHVspNode(nodes, 0, sLen - 1, -1, randSeek[s++]);
				if (s >= sLen) s = 0;
				found += found;													// fake variable usage to avoid compiler optimisations
			}
			tend = System.nanoTime();
			testRuns[1][sLen-1] = (double)(tend - tstart)/iters;
			if (verbose && (sLen % sizeStep == 0 || sLen < sizeStep)) System.out.printf("linr: %5.1f ", (double)(tend - tstart)/iters);

			NSPMatrix.nodeSearch_linearLim = 0;
			System.gc();
			tstart = System.nanoTime();
			for (int k = 0, s = 0; k < iters; k++) {
				//int csought = (int)(Math.random()*(nodes[sLen - 1].c - nodes[0].c));
				//int csought = (int)(Math.random() * cMax);
				int found = NSPMatrix.findHVspNode(nodes, 0, sLen - 1, -1, randSeek[s++]);
				if (s >= sLen) s = 0;
				found += found;													// fake variable usage to avoid compiler optimisations
			}
			tend = System.nanoTime();
			testRuns[2][sLen-1] = (double)(tend - tstart)/iters;
			if (verbose && (sLen % sizeStep == 0 || sLen < sizeStep)) System.out.printf("bina: %5.1f\n", (double)(tend - tstart)/iters);
		}
		NSPMatrix.nodeSearch_iterateLim = itrLim; NSPMatrix.nodeSearch_linearLim = linLim;
		

		for (int f = 0; f < filterIter; f++) {									// do simple multipass filtering of the data curves
			double[][] testRuns2 = testRuns.clone();
			for (int i = 0; i < arraySize; i++) {
				testRuns2[0][i] = (testRuns[0][i-1<0?0:i-1] + testRuns[0][i] + testRuns[0][i+1>=arraySize?arraySize-1:i+1]) * 0.3333;
				testRuns2[1][i] = (testRuns[1][i-1<0?0:i-1] + testRuns[1][i] + testRuns[1][i+1>=arraySize?arraySize-1:i+1]) * 0.3333;
				testRuns2[2][i] = (testRuns[2][i-1<0?0:i-1] + testRuns[2][i] + testRuns[2][i+1>=arraySize?arraySize-1:i+1]) * 0.3333;
			}
			testRuns = testRuns2;
		}
		
		// find out which order and at what element counts to implement the methods (demands proper amount of data filtering)
		int[] prioChoice = {0,Integer.MAX_VALUE,1,Integer.MAX_VALUE,2};
		if (filterIter < 10) { int[] prioChoice2 = {0,50,1,Integer.MAX_VALUE,2}; prioChoice = prioChoice2;		// the default assignment
		} else  {
			if (testRuns[0][0] > testRuns[1][0]) { int tmp = prioChoice[0]; prioChoice[0] = prioChoice[2]; prioChoice[2] = tmp; }	// sort priorities
			if (testRuns[1][0] > testRuns[2][0]) { int tmp = prioChoice[2]; prioChoice[2] = prioChoice[4]; prioChoice[4] = tmp; }
			if (testRuns[0][0] > testRuns[1][0]) { int tmp = prioChoice[0]; prioChoice[0] = prioChoice[2]; prioChoice[2] = tmp; }
			int i = 0;															// find breakeven element counts for priorities
			for (; i < arraySize; i++) { if (testRuns[prioChoice[0]][i] > testRuns[prioChoice[2]][i]) { prioChoice[1] = i - 1; break; } }
			for (; i < arraySize; i++) { if (testRuns[prioChoice[2]][i] > testRuns[prioChoice[4]][i]) { prioChoice[3] = i - 1; break; } }
		}

		
		// chart for comparison timing runs
		if (verbose) {
			XYDataset bestFitChartSet;
			String[] testNames = {"","","","","iterative search","linear search","binary search"};
			int[] testCases = {0, 1, 2};
			bestFitChartSet = createStatisticSet(testRuns, testNames, testCases);
			XYLineChart_AWT bestFitChart = new XYLineChart_AWT(
					"NSPMatrix.findHVspNode() test", "elements searched", "nanosecs/search op.", bestFitChartSet, 1024, 768);
			bestFitChart.pack();          
			RefineryUtilities.centerFrameOnScreen(bestFitChart);   
			bestFitChart.setVisible( true );
		}
		return prioChoice;
	}
	
	
	
	// test client
	public static void main(String[] args) {
				
		double[] d = {1, 0, 3, 0, 5, 0, 0, 1, 3};
		double[] d8 = {1,2,3,1,2}, d9 = {12,24,36,5,6,7}, d9b = {12,5,24,6,36,7}, v1 = {8,15,19};
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
						0,0,3,4,0,3,0,1,3,
						3,9,0,0,0,0,3,0,0,
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
		double[] d21 = {	 1,-3, 0, 0, 0, 0, 0, 0, 0,
							 0,-4, 0, 7, 0, 0, 0, 1, 0,
							-1, 1, 1, 1, 0, 0, 0, 1, 0,
							 0, 0, 0, 2, 0, 5, 0, 0, 0,
							 0, 0, 3, 1, 3, 1, 0, 1, 0,
							 0, 0, 0, 0, 0, 1, 0, 7, 0,
							 0, 0, 0, 0,-1, 1, 2, 1, 0,
							 0, 0, 0, 0, 0, 0, 0,-1, 4,
							 0, 0, 0, 0, 0, 0, 3, 7, 3};
		double[] d22 = {	 10, 2, 3, 4, 5, 0, 0, 0, 0, 0,
							  1,-.1, 3, 4, 5, 0, 0, 0, 0, 0,
							  1, 2,10, 4, 5, 0, 0, 0, 0, 0,
							  1, 2, 3,10, 5, 0, 0, 0, 0, 0,
							  1, 2, 9, 4,10, 0, 0, 0, 0, 0,
							  0, 0, 0, 0, 0, 1, 2, 2, 0, 0,
							  0, 0, 0, 0, 0, 2,10, 2, 0, 0,
							  0, 0, 0, 0, 0, 2, 2,10, 0, 0,
							  0, 0, 0, 0, 0, 0, 0, 0,10, 7,
							  0, 0, 0, 0, 0, 0, 0, 0, 7,10};
		double[] d11 = {1,2,3,4,5,6,7,8,
						2,2,3,4,5,6,7,8,
						3,3,3,4,5,6,7,8,
						4,4,4,4,5,6,7,8,
						5,5,5,5,5,6,7,8,
						6,6,6,6,6,6,7,8,
						7,7,7,7,7,7,7,8,
						8,8,8,8,8,8,8,8};
			
		// runs matrix multiplication comparison
//		int[] mtests = {4, 2, -2, -4};
//		matrixMultiTest(mtests);
//		if(1==1) return;

		int iters = 100000;
		long tStart, tEnd;
		
		
		// test element search heuristics of triple-method element finder testFindHVSpNode()
		// this test will tell what element limits to set for each type of search approach: iterative/linear/binary search
//		testFindHVSpNode(300, 150000, 10, true);

//		double[] w1D = {.55555555,.88888888,.55555555 }, w = new double[9*3];
//		for (int j = 0; j < 3; j++)
//			for (int i = 0; i < 3; i++) w[3*j+i] = w1D[i]*w1D[j];
		
		double[] iP = FEM1.integrationPoints(8, 8);
		iP = FEM1.integrationPoints(8, 27);					// test integration point coordinates calculators
		iP = FEM1.integrationWeights(8, 27);
		double[] matProps = {.1, 1};
		double[] tensor = FEM1.matStiffness_1D(matProps);
		
		// test what size of maximal node clustering per octant leaf gives fastest search of two closest nodes			
//		int optimalCluster = FEM1Benchmark.testOctreeBuildClosestNodes("data/landscape.obj", 2, 500, 20, true);	
		// test octree build single & multitasked, with variable thread spawn capping values
//		FEM1Benchmark.testBuildingOctree("data/landscape.obj", true, false, true, true);
//		if(1==1) return;
		
		// generates the permutations of the stencils used for matching in the IST algorithm
//		int[] perms = FEM1.stencilPermutations();
//		System.out.println(FEM1.stencilPermToJavaCode(perms));
		//VisitBitArray.criterionmaxResets(true);
		
		// test two example tetrahedra, one perfect and one semi-bad
//		double[] perfectTetrahedron = {-0.2309,0,-0.3266, 0.2309,-0.3266,0, -0.2309,0,0.3266, 0.2309,0.3266,0 };
//		double[] badTetrahedron = {-.2525,-0.0679,-1.9709, -0.2894,0.0215,-1.97, -0.264,0.003,-1.9, -.36,.054,-1.9 };
//		double ptVal = fem.tetraVolumePositivity(perfectTetrahedron, 0, 1, 2, 3);
//		ptVal = fem.tetraVolumePositivity(badTetrahedron, 0, 1, 2, 3);
//		ptVal = fem.tetraQuality(perfectTetrahedron, 0, 1, 2, 3);
//		ptVal = fem.tetraQuality(badTetrahedron, 0, 1, 2, 3);
//		ptVal = fem.tetraSmallestDihedral(fem.tetraDihedralAngles(perfectTetrahedron, 0, 1, 2, 3));
//		ptVal = fem.tetraSmallestDihedral(fem.tetraDihedralAngles(badTetrahedron, 0, 1, 2, 3));	
		// print out the smallest dihedral angles and qualities of perfectly shaped BCC tetrahedra of the IST stuffer
//		System.out.println(FEM1.toStringBCCtetraInfo());
		
		FEM1 fem = null;
		BufferedReader br = null;

		// test loading a single object mesh from OBJ file with smoothing groups and normals
		fem.setDebugLevel(2);
		try {	br = new BufferedReader(new FileReader("data/Beethoven.obj"));
				fem = new FEM1(br, FEM1.MESH_PSC); br.close();
		} catch (FileNotFoundException e) { e.printStackTrace();
		} catch (IOException e) { e.printStackTrace(); }
		
		// test creating an octree from  FEM1 system elements: initialising the first octant with all data, then calling subdivider buildOctree()
		//FEM1Octree octree = new FEM1Octree(fem, FEM1Octree.DO_NODES|FEM1Octree.DO_EDGES|FEM1Octree.DO_FACETS);
		tStart = System.nanoTime();
		FEM1Octree octree = new FEM1Octree(fem, FEM1Octree.DO_FACETS, 0);			// DO_FACETS activates automatic calculation of maxLevel
		octree.root.build(octree, 0, true);	
		fem.setDebugLevel(2);
		FEM1Octree latticeTree = fem.volumeMeshIST(octree, 0.2, 5, true);			// subdivide to leaf octant size of 0.2m (= width of smallest element)
		latticeTree.toOBJ(3, 7, false, true, 0);									// test octree output to OBJ
		System.out.printf("FEM1.volumeMeshIST() took %d ns\n", System.nanoTime() - tStart);
		fem.toOBJ(true, true);													// test outputting tetrahedral data to OBJ file
		
		// microbenchmark test of IST volume generator
//		fem.setDebugLevel(1);
//		fem.clearISTsolution();
//		tStart = System.nanoTime();
//		long timeOT = 0;
//		int iNum = 50;
//		for (int i = 0; i < iNum; i++) {
//			long timerOT = System.nanoTime();
//			octree = new FEM1Octree(fem, FEM1Octree.DO_FACETS, 0);					// DO_FACETS activates automatic calculation of maxLevel
//			octree.root.build(octree, 0, true);
//			timeOT += System.nanoTime() - timerOT;
//			//if (i==iNum-1) fem.setDebugLevel(2);
//			fem.volumeMeshIST(octree, 0.1, 0, true);								// subdivide to leaf octant size of 0.1m
//		}
//		System.out.printf("FEM1OCtant.build() took %d ns\n", timeOT/iNum);
//		System.out.printf("GeocTree + IST generation took %d ns from average of %d iterations\n", (System.nanoTime() - tStart)/iNum, iNum);
				
		FEM1.getExecutor(0).shutdown();
		try { while (!FEM1.getExecutor(0).awaitTermination(10, TimeUnit.SECONDS));
		} catch (InterruptedException e) { e.printStackTrace(); }				
		if(1==1) return;
		
		FEM1Octree lTree = fem.latticeTree(octree, 57, false);
		// test recursive octant collector versus nonrecursive one (nonrecursive seems faster)
//		FEM1Octant[] leafArray = null, leafArray2 = {null, null}; tStart = System.nanoTime();
//		for (int i = 0; i < 5000; i++) {
//			leafArray = lTree.root.leafOctantArray(lTree, lTree.maxLevel);	// nonrecursive
//			leafArray = lTree.root.leafOctantArrayR(lTree, lTree.maxLevel);	// recursive
//			leafArray2[0] = leafArray[0]; }
//		long tTime = System.nanoTime() - tStart; System.out.println((double)tTime/5000);
//		leafArray2[0] = leafArray2[1];
		
		// test to track down an octant and get it's 1-neighbourhood or 2-neighbourhoos
		FEM1Octant oct = lTree.root.locateCoordinate(2.2341, -0.1424, 3.0697, lTree.topLevel);
		FEM1Octant[] octA = new FEM1Octant[6];
		oct.getFaceNeighbours(lTree, octA, null);
		octA = oct.getFaceEdgeNeighbours(lTree, null, 0, -1);
		
		// test collection functions for different levels of octants
//		FEM1Octant[] leaves = lTree.root.leafOctantArray(lTree);
//		FEM1Octant[][] layers = lTree.root.layerOctantArray(lTree);
		lTree.toOBJ(true, false, 0);
		

		// test doing a 2:1 Weak Condition subdivision of octree, where every neighbour can only have a neighbour twice the size but not larger
		// a criterion important for differential equation fields or construction of graded mesh volumes
		//lTree.root.split2to1(lTree, true);
		//lTree.root.toOBJ(lTree, false, 0);
		//if(1==1) return;
		
		tStart = tEnd = 0;
		for (int i = 0; i < 600; i++) {
			lTree = fem.latticeTree(octree, 64, true);
			tStart = System.nanoTime();
			lTree.root.split2to1(lTree, false);
			tEnd += System.nanoTime() - tStart;
		}
		System.out.printf("FEM1.split2to1() %.1f ns\n", (double)tEnd/600.0);
		FEM1.getExecutor(0).shutdown();
		try { while (!FEM1.getExecutor(0).awaitTermination(10, TimeUnit.SECONDS));
		} catch (InterruptedException e) { e.printStackTrace(); }				

		if(1==1) return;

		tEnd = 0;
		long tEnd2 = 0;
		for (int i = 0; i < 100; i++) {
			long tStart2 = System.nanoTime();
			lTree = fem.latticeTree(octree, 64, false);
			tEnd2 += System.nanoTime() - tStart2;
			tStart = System.nanoTime();
			lTree.root.split2to1(lTree, false);
			tEnd += System.nanoTime() - tStart;
		}
		System.out.printf("FEM1.latticeTree() averaged %.1f ns\n", (double)(tEnd2) / 100.0);
		System.out.printf("FEM1Octree.split2to1_2() averaged %.1f ns\n", (double)(tEnd) / 100.0);
		//lTree.root.toOBJ(lTree, true, FEM1Octree.DO_NODES);
		if(1==1) return;
		
		// test getting bounding box from an edge and a facet
		double[] bbox = fem.edgeBBox(0);
		// test getting octree leaves overlapped by a bounding box
		for (int f = 0; f < fem.polygons; f++) {
			bbox = fem.facetBBox(f);
			FEM1Octant[] overlap = octree.root.octantArrayByBBox(octree, bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5], false, false);
			overlap[0] = overlap[1];
		}
		
		//octree.root.addEdge(octree, 0, -43.8261, 3.9602, 0, -66.2773, 26.6212, 0, octree.maxEdges, 0);
		//octree.root.addEdge(octree, -84, 77.8848, 3.9602, 0, -66.2773, 26.6212, 0, octree.maxEdges, 0);
		//octree.root.toOBJ(octree, true, FEM1Octree.DO_NODES);
//		octree.root.toOBJ(octree, true, FEM1Octree.DO_EDGES);
//		octree.root.toOBJ(octree, true, FEM1Octree.DO_FACETS);
		// export closestNodes() visualisation to OBJ file
		fem.closestNodesToOBJ(octree);
		
//		FEM1Octree[] octantsClosest = new FEM1Octree[2];
//		double[] distanceClosest = new double[1];
//		for (int n = 0; n < fem.nodes; n++) {
//			int nClosest = octree.closestNode(fem, n, distanceClosest, octantsClosest);
//			System.out.println("Closest pair " + n + "|" + nClosest + ", d: " + distanceClosest[0]);
//		}
		
		// test loading a single object from OBJ file with handmade tetrahedral data into element objects		
		try {
			br = new BufferedReader(new FileReader("data/tetrahedral1.obj"));
			fem = new FEM1(br, FEM1.MESH_HANDMADE_OBJECTIFIED);
			br.close();
		} catch (IOException e) { e.printStackTrace(); }
		
		// test enclosure of a coordinate within tetrahedron 0
		int enclosure = fem.getElement2(0).tetraVertexEnclosure(0.1791, -0.3527, 2.7746, true);
		double[] isect = new double[12];
		// test intersection of a segment with tetrahedron 0
		enclosure = fem.getElement2(0).tetraSegmentIntersection(0.571, -0.2194, 1.918, -0.181, -0.5808, 0.9791, true, isect);
		System.out.println(FEM1Element.toBitsInteger(enclosure, "                  iiiiFFFFeFFFFe\n                  43214321e4321e\n"));
		
		// determine what nodes are internal (not neighbouring boundary) of FEM-system
		int[] internalNodes = fem.internalNodes();
		// test node deletion and resulting neighbourhood readjustment
//		int[] deletionRef = {0, 1, 2};
//		fem.deleteElements(deletionRef);
		
		// test method that locates optimal position of a fourth tetrahedral node relative to 3 supplied previous nodes
		FEM1Element optElem = fem.optimalTetrahedron(0, 1, 2, fem.addNode(0, 0, 0, (byte)0));
		
		// test sorting the FEM-system's tetrahedra according to quality
		tStart = System.nanoTime();
		int[] tWorst = fem.tetraWorstSort(10, .8, .2, false);
		tEnd = System.nanoTime();
		System.out.printf("FEM1.tetraWorstSort() averaged %.1f ns\n", (double)(tEnd - tStart));
		
		// run test of every test criterion appliable to a tetrahedron, for every tetrahedron of FEM-system
		for (int e = 0; e < fem.elements2; e++) {
			FEM1Element elem = fem.getElement2(e);
			if (elem == null) continue;
			System.out.println("Element " + e + ":");
			System.out.println("Circumscribed: " + elem.tetraCircumRadius());
			double[] cCenter = elem.tetraCircumcenter();
			System.out.println("Circumcenter: " + cCenter[0] + ", " + cCenter[1] + ", " + cCenter[2]);
			double[] cBary = elem.tetraBarycentric();
			System.out.println("Barycentric: " + cBary[0] + ", " + cBary[1] + ", " + cBary[2] + ", " + cBary[3]);
			elem.propagateInterfaces(fem, FEM1Element.PROPAGATE_EDGES, false);
			System.out.println("Inscribed: " + elem.tetraInscribedRadius());
			elem.propagateInterfaces(fem, FEM1Element.PROPAGATE_AREAS, false);
			System.out.println("Quality: " + elem.tetraQuality());
			double[] massC = elem.tetraMassCenter();
			System.out.println("Mass centre: " + massC[0] + ", " + massC[1] + ", " + massC[2]);
			elem.tetraVolumeGradient(true, false, true, false);
			System.out.println("Volume gradient node0: " + elem.data[10] + ", " + elem.data[11] + ", " + elem.data[12]);
			System.out.println("Volume gradient node1: " + elem.data[13] + ", " + elem.data[14] + ", " + elem.data[15]);
			System.out.println("Volume gradient node2: " + elem.data[16] + ", " + elem.data[17] + ", " + elem.data[18]);
			System.out.println("Volume gradient node3: " + elem.data[19] + ", " + elem.data[20] + ", " + elem.data[21]);
		}
		if(1==1) return;
		
		// test first type of FEM-system (Long Chen, "Programming of Finite Elements Method in Matlab")
		// TODO: unclear how Matlab treats partial matrix indexation, or a defect in Chen's code gives unpredicted results
		NSPMatrix Mfem = fem.assembleNSPMatrix();			// test FEM1 tetrahedral stiffness matrix assembly into an NSPMatrix
		new SizeOf(Mfem, true);								// test SizeOf() method, determining total memory usage of an arbitrary object
		Matrix Mfem2 = fem.assembleMatrix();				// test FEM1 tetrahedral stiffness matrix assembly into a Matrix
		new SizeOf(Mfem2, true);
			
		
		boolean testAllSolvers = false;
		Matrix[] LU2 = null;
		if (testAllSolvers) {
			Matrix A10x10 = new Matrix("Asubp", 10, 10, d22, null);
			Matrix v10x1 = new Matrix("b", 10, 1, d2, null);
			//A10x10.decomposeLU(NO_COPY, false, false);
			v10x1.backSubstituteLU(A10x10, null, COPY);

			NSPMatrix N10x10 = new NSPMatrix("Asubp", 10, 10, d22, null);
			NSPMatrix vN10x1 = new NSPMatrix("b", 10, 1, d2, null);
			vN10x1.backSubstituteLU(N10x10, null, COPY);

			A10x10.solveGaussJordanPP(v10x1, null);
			A10x10.solveGaussJordanFP(v10x1, COPY);
			A10x10.solveGaussPP(v10x1);
			A10x10.solveGaussSeidel(v10x1, 10, 0.02);
			
			// test direct construction of a FrontalMatrix, it's LU-decomposition and conversion to a Matrix, then backsubstitution with a problem vector
			FrontalMatrix FA10x10 = new FrontalMatrix(10, 10, d22);
			FA10x10.decomposeLU(true);
			System.out.println(FA10x10.toString());
			Matrix LU10x10f = new Matrix("Afm", FA10x10);
			v10x1 = new Matrix("b", 10, 1, d2, null);
			v10x1.backSubstituteLU(null, LU10x10f, false);
			LU2 = A10x10.decomposeLU();
		}
		
		MatrixMarketIO mmIO = new MatrixMarketIO("data/d_dyn.mtx", 2, true);
		//MatrixMarketIO mmIO = new MatrixMarketIO("data/mcca.mtx", 2, true);
		NSPMatrix MM2 = mmIO.getNSPMatrix();
		//MM2.toFile(0);
		new SizeOf(MM2, true);
		if(1==1) return;
		
		boolean doMaxMatchingPermutation = true;
		if (doMaxMatchingPermutation) {
			// test maximum matching finding with Hopcroft-Karp algorithm
			BipartiteDM biDM = new BipartiteDM(MM2);
			//biDM.toString();
			biDM.maximumMatchingHK();
			biDM.findBTF();
			// test permuting the paired rows-columns to the diagonal
			MM2 = MM2.permuteMaximumMatching(biDM, NO_COPY, true);
			//MM2.toFile(0);
		}

		// test construction of multifrontal unsymmetric DAG with transitive reduction of symbolic factorisation graph
		FrontalDAG dataDAG = null;
		//MM2 = new NSPMatrix("MM2", 10, 10, d22, null);
		tStart = System.nanoTime();
		for (int i = 0; i < 1; i++) {
			FrontalDAG taskDAG = FrontalDAG.taskDAG(MM2, false);
			//taskDAG.toGraphViz(true, true, true);
			dataDAG = taskDAG.dataDAG(MM2);
			NSPMatrix LU = taskDAG.decomposeLU(MM2);
			System.out.println(LU.toString());
		}
		tEnd = System.nanoTime();
		System.out.printf("taskDAG.decomposeLU() averaged %.1f ns\n", (double)(tEnd - tStart)/10000);
		dataDAG.clearVisits();
		if(1==1) return;

		// test LU decomposure for normal & NSP sparse matrix types, no output to image, nor file, nor graphviz
		testLUdecomposure("data/mcca.mtx", false, false, false, 10000);
		//testLUdecomposure("", false, false, false, 1);
		
		// test sparse dynamic NSPMatrix, creation, multiplying, printout,
		// zeroes purging, value setting and zeroing, row/column swapping
		NSPMatrix G8 = new NSPMatrix("G", 9, 7, d4, null);
		NSPMatrix G8b = G8.transpose(COPY);
		G8 = G8.multiply(G8b);
		System.out.println(G8.toString());
		
		// test zero value purging of NSPMatrix by setting some values to zero and calling purgeZeroes()
		G8.Hsp[1].array[0].v = 0; G8.Hsp[2].array[0].v = 0; G8.Hsp[3].array[0].v = 0;
		System.out.println(G8.purgeZeroes());
		System.out.println(G8.toString());
		//if(1==1) return;

		// test the valueTo() & valueOf() getters and setters of NSPMattix
		G8 = new NSPMatrix("G", 9, 9, d3, null);
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) { double v = G8.valueOf(i, j); System.out.print((v != 0 ? v : " - ") + "  "); }
			System.out.println("\n");
		}
		System.out.println("\n");
		G8.valueTo(8, 0, -10000); G8.valueTo(2, 0, 0); G8.valueTo(6, 0, 0); G8.valueTo(6, 0, 7);
		G8.swapHVspArrays(8, 0, 1);
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) { double v = G8.valueOf(i, j); System.out.print((v != 0 ? v : " - ") + "  "); }
			System.out.println("\n");
		}
		
		// test sparse row inner product method multiplyHVsp() of NSPMatrix
		double n = NSPMatrix.multiplyHVsp(G8.Hsp[1], G8.Hsp[0], 0, 0);
		System.out.println("multiplyHVsp: " + n);
		NSPArray nd = NSPMatrix.addHVsp(G8.Hsp[0], G8.Vsp[6], 2.0, 0, 1, 1);
		System.out.println(nd.toString());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) { double v = G8.valueOf(i, j); System.out.print((v != 0 ? v : " - ") + "  "); }
			System.out.println("\n");
		}
		//if(1==1) return;

		
		// Test Cholesky factorisation
		//Matrix Ch = new Matrix("C", 3, 3, testCh, null);
		Matrix Ch = new Matrix("C", 10, 10, d20, null);
		// convert unsymmetric matrix to symmetric, A^T.A style, since Cholesky only works for symmetric
		Ch = Ch.transpose(COPY).multiply(Ch);
		Ch.factoriseCholesky();
		
		// Test Householder reduction form
		Matrix HH = new Matrix("HH", 4, 4, testHH, null);
		HH = HH.reduceHouseholder(COPY);
		
		// Test MatrixMarket file loading routine
		MatrixMarketIO mmIO2 = new MatrixMarketIO("data/mcca.mtx", 0, true);
		Matrix MM = mmIO2.getMatrix();
		System.out.println(MM.toString());
		// Test conversion of Matrix format to CSR sparse format
		CSRMatrix MMcsr = CSRMatrix.convert(MM);
		System.out.println(MMcsr.toString());	
//		Matrix MM_Cholesky = MM.factoriseCholesky();	
//		if (MM_Cholesky != null) {
//			MatrixBMPImage MM_image = new MatrixBMPImage(MM_Cholesky);
//			MM_image.write();
//		}
		
		// test if a matrix will pass convergence criterion for a Gauss-Seidel type iteration
		Matrix D8 = new Matrix("D8", 8, 8, d11, null);
		D8.isConvergent();

		// Test addition
		Matrix D2 = new Matrix("D2", 6, 5, d2, null);
		System.out.println(D2.toString());
		D2 = D2.add(D2, COPY);
		// Test multiplication
		D2 = D2.multiply(D2.transpose(COPY));

		// Test Gram-Schmidt orthonormalisation for non-tranposed & transposed
		Matrix O = new Matrix("O", 5, 5, d2, null);
		System.out.println(O.toString());
		Matrix Q = O.orthogonalise(COPY);
		O.transpose(COPY).orthogonalise(COPY);
		
		// test LU = A decomposition with permutation and A = LU recomposition with unpermutation
		Matrix[] lLU = Q.decomposeLU(COPY, true, true);
		Matrix LU = lLU[0].multiply(lLU[1]);
		LU.unpermute();
		System.out.println(LU);
		
		// Test eigenvalue & eigenvector finder findEigenQR(), which takes an orthonormalised matrix
		Matrix E = O.findEigenQR(Q, 150, 30, 0.005);

		// test QR decomposition
		O.decomposeQR(Q);
		
		// Test Power Method for finding largest eigenvalue
		D2 = new Matrix("D2", 5, 5, d2, null);
		Matrix x5 = new Matrix("x5", 5, 1, d8, null);
		System.out.println("Largest eigenvalue: " + D2.eigenPowerValue(x5, 0.01, 100) + "\n");
		
		// Test method that finds out how much sparse matrix data is "sticking out" from the diagonal
		System.out.println("half bandwidth: " + new Matrix("Dd", 10, 10, d20, null).getHalfBandwidth());
		
		iters = 100;

		// Test partial pivoting Gauss solver
		Matrix D5 = new Matrix("A", 3, 3, d5, null);
		Matrix x6 = D5.solveGaussPP(new Matrix("c2", 3, 1, d9, null));
		// test remultiplying it
		D5.multiply(x6);
		
		Matrix N = new Matrix("N", 3, 3, cn, null);
		N.isConvergent();
		x6 = N.solveGaussSeidel(new Matrix("v1", 3, 1, v1, null), 100, 0.1);
		
		// Test pivoting LU decomposer and back substitution methods
		Matrix[] bLU = new Matrix("c2", 3, 1, d9, null).backSubstituteLU(D5, null, COPY);
		if (bLU == null) System.out.println("backSubstituteLU2() received singular matrix.");

		// Test partial pivoting Gauss Jordan with inverse matrix creation and sparse diagonal matrix optimisation
		Matrix D20 = new Matrix("D", 10, 10, d20, null);
		x6 = D20.solveGaussJordanPPDO(
				new Matrix("b", 10, 1, d10, null),
				new Matrix("D", 10, 10, Matrix.Type.Null, 1));
		if (x6 == null) System.out.println("solveGaussJordan() returned singular matrix.");
		System.out.println("Determinant: " + D20.det);
		
		// Test full pivoting Gauss-Jordan with in-situ matrix & input vector -> solution vector transformation
		Matrix c2 = new Matrix("c2", 3, 2, d9b, null), X;
		Matrix.DEBUG_LEVEL--;
		tStart = System.nanoTime();
		for (int i = 0; i < iters; i++)
			X = D5.solveGaussJordanFP(c2, COPY);
		tEnd = System.nanoTime();
		Matrix.DEBUG_LEVEL++;
		System.out.printf("solveGaussJordanFP() averaged %.1f ns\n", (double)(tEnd - tStart)/iters);

		// Test ordinary multiply versus Strassen-Winograd multiply algorithm
		Matrix Q1 = new Matrix("Q1", 128, 128, Matrix.Type.Random, 1), Q11;
		Matrix Q2 = new Matrix("Q2", 128, 128, Matrix.Type.Random, 1), Q22;
		Matrix.DEBUG_LEVEL--;
		for (int tests = 0; tests < 5; tests++) {

			tStart = System.nanoTime();
			for (int i = 0; i < iters; i++) Q11 = Q1.multiply(Q1);
			tEnd = System.nanoTime();

			System.out.printf("multiply() averaged %.1f ns\n", (double)(tEnd - tStart)/iters);
		    System.out.println("Mults: " + Matrix.mulFlops_DEBUG);
		    System.out.println("Adds: " + Matrix.mulAdops_DEBUG);

		    tStart = System.nanoTime();
			for (int i = 0; i < iters; i++) Q22 = Matrix.multiplyStrasWin(Q2, Q2, 8);
			tEnd = System.nanoTime();

			System.out.printf("multiplyStrasWin() averaged %.1f ns\n", (double)(tEnd - tStart)/iters);
		    System.out.println("Mults: " + Matrix.mulFlopsSW_DEBUG);
		    System.out.println("Adds: " + Matrix.mulAdopsSW_DEBUG);
		    System.out.println("Recurses: " + Matrix.mulSW_DEBUG);
		}
		Matrix.DEBUG_LEVEL++;
		
		// Test if multiplication of rescaled matrices produce same results as unexpanded ones
		Matrix S1 = new Matrix("S1", 6, 5, d2, null);
		S1 = S1.multiply(S1.transpose(COPY));
		S1 = new Matrix("S1", 6, 5, d2, null);
		S1 = S1.rescale(0, 0, 6, 6, COPY);
		S1 = S1.multiply(S1.transpose(COPY));
		
		// Test multiplying matrix with a value
		Matrix A1 = new Matrix("A1", 9, 9, d3, null);
		A1 = A1.multiply(-2000, NO_COPY);

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

		// Test CSR matrices, test vAv style inner product
		CSRMatrix x1 = new CSRMatrix("x1", 1, 9, d, null);
		CSRMatrix x2 = new CSRMatrix("x2", 7, 1, d, null);
		CSRMatrix V = new CSRMatrix("V", 9, 7, d4, null);
		x2.swap(6, 2);
		CSRMatrix S = x1.multiply(V).multiply(x2);

		// Test CSRMatrix conversion and polymorphic access
		CSRMatrix O2 = new CSRMatrix("O", 9, 7, d3, null);
		O2.transpose(COPY);
		G = O2.multiply(new CSRMatrix("R", 9, 7, d3, null));
		G = G.eliminateRowColumn(6, 6, true);

		// test centering method for Matrix
		Matrix P = G.center(COPY);

		// test BinBitImage methods of bitflagging nonzero values, comparing matrices through their bit images
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
		CSRMatrix csrA = new CSRMatrix("A2", 5, 5,  Matrix.Type.Random, 1);

		CSRMatrix csrS = new CSRMatrix("S", 5, 5, Matrix.Type.Null, 1);
		CSRMatrix csrT = new CSRMatrix("T", 5, 5, Matrix.Type.Identity, 1);
		csrT.valueTo(1, 3, 9.9);
		csrS.valueTo(2, 2, 9.9);
		csrS.valueTo(0, 4, -8);
		CSRMatrix csrU = csrS.add(csrT, COPY);
		CSRMatrix csrD = csrU.multiply(csrT);
		CSRMatrix B = csrD.clone();
		
		Matrix C = new Matrix("C", 5, 5, Matrix.Type.Identity, 1);
		System.out.println("Identity matrix\n" + C.toString());
		System.out.println();

		System.out.println("Multiplied matrices polymorphically\n" + csrA.multiply(B));
		System.out.println();

		// shouldn't be equal since AB != BA in general
		csrS = new CSRMatrix("S", 5, 5, d2, null);
		csrT = new CSRMatrix("T", 5, 5, d2, null);
		csrT = csrT.add(5, NO_COPY);
		System.out.println(csrS.multiply(csrT).equals(csrT.multiply(csrS)));
		System.out.println();

		Matrix b = new Matrix("b", 5, 1, Matrix.Type.Random, 1);

		Matrix A = new Matrix("A3", csrA.M, csrA.N, Matrix.Type.Null, 1);
		A = csrA.clone();
		Matrix x = A.solveGaussPP(b);
		A = A.multiply(x);
		
	}

}
