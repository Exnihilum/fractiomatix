package lkr74.fem1;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.concurrent.TimeUnit;
import java.awt.BasicStroke; 
import org.jfree.chart.ChartPanel; 
import org.jfree.chart.JFreeChart; 
import org.jfree.data.xy.XYDataset; 
import org.jfree.data.xy.XYSeries; 
import org.jfree.ui.ApplicationFrame; 
import org.jfree.chart.plot.XYPlot; 
import org.jfree.chart.ChartFactory; 
import org.jfree.chart.plot.PlotOrientation; 
import org.jfree.data.xy.XYSeriesCollection; 
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;

import org.jfree.ui.RefineryUtilities;
import lkr74.fem1.FEM1;
import lkr74.fem1.FEM1Octant;
import lkr74.mathgenerics.XYLineChart_AWT;

public class FEM1Benchmark extends ApplicationFrame {
	private static final long serialVersionUID = -1795220774825920995L;

	public FEM1Benchmark( String title, String xtag, String ytag, XYDataset dataset, int width, int height) {
		super(title);
		JFreeChart xylineChart = ChartFactory.createXYLineChart(title,xtag,ytag,dataset,PlotOrientation.VERTICAL,true,true,false);
		ChartPanel chartPanel = new ChartPanel( xylineChart );
		chartPanel.setPreferredSize( new java.awt.Dimension( width , height ) );
		final XYPlot plot = xylineChart.getXYPlot();
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
		renderer.setSeriesPaint(0 , Color.BLACK );	renderer.setSeriesPaint(1 , Color.BLUE );
		renderer.setSeriesPaint(2 , Color.RED );	renderer.setSeriesPaint(3 , Color.GREEN );
		renderer.setSeriesStroke(0 , new BasicStroke(2));	renderer.setSeriesStroke(1 , new BasicStroke(2));
		renderer.setSeriesStroke(2 , new BasicStroke(2));	renderer.setSeriesStroke(3 , new BasicStroke(2));
		plot.setRenderer( renderer ); 
		setContentPane( chartPanel ); 
	}

	// create Chart data from timings of 2 compared tested functions, showing variation in one parameter (ex: size/fillrate)
	static XYDataset createStatisticSet(double[][] timingLists, String[] testNames, int[] testCases) {
		final XYSeriesCollection dataset = new XYSeriesCollection();          
		for (int j = 0; j < timingLists.length; j++) {
			final XYSeries testSeries = new XYSeries(testNames[testCases[j]+4], false);
			for (int i = 0; i < timingLists[j].length; i++) testSeries.add(i+2, timingLists[j][i]);
			dataset.addSeries(testSeries); }
		return dataset;
	}

	
	// test what size of maximal node clustering per octant leaf gives fastest search of two closest nodes, method returns optimal cluster size for the data
	// inputs: iterate creating and node searching octrees from clusterSize1 to clusterSize2
	// filterIter tells how many filtering rounds to execute on data curve
	final static int DO_NODES = FEM1Octant.DO_NODES, DO_EDGES = FEM1Octant.DO_EDGES, DO_FACETS = FEM1Octant.DO_FACETS;
	public static int testOctreeBuildClosestNodes(String dataFile, int clusterSize1, int clusterSize2, int filterIter, boolean verbose) {
		
		if (clusterSize1 < 2) clusterSize1 = 2;
		if (clusterSize2 < clusterSize1) clusterSize2 = clusterSize1 + 1;
		
		double[][] testRuns = new double[3][clusterSize2 - clusterSize1];
		FEM1 fem = null;
		// OBJ data loading block
		try {	BufferedReader br = new BufferedReader(new FileReader(dataFile));
				fem = new FEM1("benchmark", br, FEM1.MESH_PSC); br.close();
		} catch (FileNotFoundException e) { e.printStackTrace(); } catch (IOException e) { e.printStackTrace(); }
			
		// test what size of maximal node clustering per octant leaf gives fastest search of two closest nodes
		long tStart, tEnd, tEnd2;
		for (int i = clusterSize1, warmup = 10; i < clusterSize2;) {
			System.gc();
			tStart = System.nanoTime();
			tEnd = 0; tEnd2 = 0;
			FEM1Octree octree = new FEM1Octree(fem, DO_NODES|DO_EDGES|DO_FACETS, 0);
			octree.maxNodes = i;
			octree.root.build(octree, 0, false);
			tEnd = System.nanoTime() - tStart;
			int[] closest = new int[fem.nodes * 2];					// worst case allocation
			double[] distance = new double[fem.nodes];				// worst case allocation
			int nPairs = fem.closestNodePairs(octree, closest, distance, null);
			tEnd2 = System.nanoTime() - tStart;
			if (warmup-- > 0) continue; else {
				testRuns[0][i - clusterSize1] = tEnd *.001;
				testRuns[1][i - clusterSize1] = tEnd2 *.001;
				testRuns[2][i - clusterSize1] = (tEnd + tEnd2) *.001;
				if (verbose)	System.out.printf("Max node cluster size " + i++ + " took %.1f ns, " + nPairs + " closest pairs took %.1f ns, sum: %.1f\n",
													(double)tEnd, (double)tEnd2, (double)(tEnd+tEnd2));
			}
		}
		int clSizeRange = clusterSize2 - clusterSize1 - 1;
		for (int f = 0; f < filterIter; f++) {					// apply simple averaging filtering
			double[][] testRuns2 = testRuns.clone();
			for (int i = 1; i < clSizeRange; i++) {
				//testRuns2[0][i] = (testRuns[0][i-1] + testRuns[0][i] + testRuns[0][i+1]) * 0.3333;
				//testRuns2[1][i] = (testRuns[1][i-1] + testRuns[1][i] + testRuns[1][i+1]) * 0.3333;
				testRuns2[2][i] = (testRuns[2][i-1] + testRuns[2][i] + testRuns[2][i+1]) * 0.3333;
			}
			testRuns = testRuns2;
		}
		double lowest = testRuns[2][0];
		int optimalSize = clusterSize1;
		for (int i = 1; i <= clSizeRange; i++) if (testRuns[2][i] < lowest) { lowest = testRuns[2][i]; optimalSize = i + clusterSize1; }
		if (verbose) {
			System.out.println("Optimal cluster size: " + optimalSize);
			XYDataset bestFitChartSet;
			String[] testNames = {"","","","","octree build","closest node search","build + search"};
			int[] testCases = {0, 1, 2};
			bestFitChartSet = createStatisticSet(testRuns, testNames, testCases);
			XYLineChart_AWT bestFitChart = new XYLineChart_AWT(
					"Octree varying cluster size build and closest node search test", "exec. time", "msecs/op.", bestFitChartSet, 1024, 768);
			bestFitChart.pack();          
			RefineryUtilities.centerFrameOnScreen(bestFitChart);   
			bestFitChart.setVisible( true ); 
		}	
		return optimalSize;
	}
	
	
	public static void testBuildingOctree(String dataFile, boolean testSpeedCompOnly, boolean testMaxItems, boolean testThreadSplitF, boolean verbose) {
		// test creating an octree from elements of fem system
		FEM1 fem = null;
		// OBJ data loading block
		try {	BufferedReader br = new BufferedReader(new FileReader(dataFile));
				fem = new FEM1("benchmark", br, FEM1.MESH_PSC); br.close();
		} catch (FileNotFoundException e) { e.printStackTrace(); } catch (IOException e) { e.printStackTrace(); }

		double time_thread=0, time_no_thread=0;
		int altNo = 120 + 8;
		for (int alt = 0; alt < altNo; alt++) {					// do alternations of threaded/nonthreaded building, summing process times
			System.gc();
			FEM1Octree octree = new FEM1Octree(fem, DO_NODES|DO_EDGES|DO_FACETS, 0);
			long tStart = System.nanoTime();
			octree.root.build(octree, 0, false);
			//System.out.println(octree.splitNtiming + " " + octree.splitEtiming + " " + octree.splitFtiming + " ");
			if (alt < 8) time_no_thread += (System.nanoTime() - tStart);
	
			System.gc();
			octree = new FEM1Octree(fem, DO_NODES|DO_EDGES|DO_FACETS, 0);
			tStart = System.nanoTime();
			octree.root.build(octree, 0, true);
			if (alt < 8) time_thread += (System.nanoTime() - tStart);
		}
		if (verbose) {
			double tnt, tt;
			altNo -= 8;
			System.out.printf("FEM1.buildOctree() single task took average %.1f ns\n", tnt = time_no_thread / (double)altNo);
			System.out.printf("FEM1.buildOctree() multitask took average %.1f ns\n", tt = time_thread / (double)altNo);
			System.out.printf("Factor: %.1f\n", tnt / tt);
		}
		if (testSpeedCompOnly) {
			FEM1.getExecutor(0).shutdown();
			try { FEM1.getExecutor(0).awaitTermination(10, TimeUnit.SECONDS);
			} catch (InterruptedException e) { e.printStackTrace(); }
			return;
		}

		double[][] testRuns = new double[1][200];
		double avgRuns = 4;
		int samples = 200;
		
		// either check what amount of branch processing is spedient for motivating spawning of a new thread
		// or check how maxItems limiter of branching disbalance/growth influences multitasking
		for (int iC = 1, t = 0; iC < samples; iC++, t++) {				// do test of multithreaded with increasing item count cutoff for triggering new threads
			System.gc();
			long tStart = System.nanoTime();
			for (int i = 0; i < avgRuns; i++) {						// average from 20 multithread runs
				FEM1Octree octree = new FEM1Octree(fem, DO_NODES|DO_EDGES|DO_FACETS, 0);
				if (testThreadSplitF) octree.threadSplitF = iC;
				else {
					if (!testMaxItems) octree.threadSpawnItemsCap = iC*5;
					else { octree.maxNodes = octree.maxEdges = octree.maxFacets = iC; }
				}
				octree.root.build(octree, 0, true);
			}
			testRuns[0][t] = (double)(System.nanoTime() - tStart) / avgRuns;
			if (testThreadSplitF)
					System.out.printf("FEM1.buildOctree() multitask with thread split factor" + iC + " took %.1f ns\n", testRuns[0][t]);
			else	System.out.printf("FEM1.buildOctree() multitask with " + (testMaxItems ? "items max "+iC : "item cap "+iC*5) + " took %.1f ns\n", testRuns[0][t]);
		}
//		for (int f = 0; f < 10; f++) {								// apply simple averaging filtering
//			double[] testRuns0 = new double[samples];
//			testRuns0[0] = (testRuns[0][1] + testRuns[0][0] + testRuns[0][2]) * 0.3333;
//			for (int i = 1; i < samples - 1; i++) {
//				testRuns0[i] = (testRuns[0][i-1] + testRuns[0][i] + testRuns[0][i+1]) * 0.3333;
//			}
//			testRuns0[samples - 1] = (testRuns[0][samples-2] + testRuns[0][samples - 1] + testRuns[0][samples - 3]) * 0.3333;
//			testRuns[0] = testRuns0;
//		}
		int optimalSize = 0;										// find optimal items cap for recursive thread spawning
		if (!testThreadSplitF) {
			double lowest = testRuns[0][0];
			for (int i = 1; i < samples; i++) if (testRuns[0][i] < lowest) { lowest = testRuns[0][i]; optimalSize = i; }
		}
		if (verbose) {
			if (optimalSize != 0) System.out.println("Optimal minimum branch items size for thread spawn: " + optimalSize);
			XYDataset bestFitChartSet;
			String[] testNames = {"","","","","octree build"};
			int[] testCases = {0, 1};
			bestFitChartSet = createStatisticSet(testRuns, testNames, testCases);
			XYLineChart_AWT bestFitChart = new XYLineChart_AWT(
					"", "item count thread spawn trigger cap", "nsecs/op.", bestFitChartSet, 1024, 768);
			bestFitChart.pack();          
			RefineryUtilities.centerFrameOnScreen(bestFitChart);   
			bestFitChart.setVisible( true ); 
		}	

		FEM1.getExecutor(0).shutdown();
		try { FEM1.getExecutor(0).awaitTermination(10, TimeUnit.SECONDS);
		} catch (InterruptedException e) { e.printStackTrace(); }
	}
	
	
	// test different nodeWork allocation factors, comparing speed of execution (nwfNo = nodeWork factor number)
	public static int testAllocFactor(int nwfNo, int size, boolean verbose) {
		double[][] testRuns = new double[3][nwfNo];
		boolean warmup = true;
		for (int nwf = 1; nwf < nwfNo; nwf++) {
			FEM1 fem2 = new FEM1("test");
			fem2.nodeworkFactor = nwf;
			long tStart = System.nanoTime();
			for (int n = 0; n < size; n++) fem2.addNode(0, 0, 0, (byte)0);
			for (int n = 0; n < size; n++) fem2.deleteNode((int)(Math.random() * fem2.nodes));
			long tEnd = System.nanoTime();
			if (warmup) { nwf--; warmup = false; continue; }
			testRuns[0][nwf - 1] = (double)(tEnd - tStart) / 1000000.0;
			testRuns[1][nwf - 1] = 0.1 * (double)(tEnd - tStart) / (double)(fem2.node.length + fem2.nodeWork.length);
			tStart = System.nanoTime();
			for (int n = 0; n < size; n++)
				if (Math.random() > 0.2)	fem2.addNode(0, 0, 0, (byte)0);
				else						fem2.deleteNode((int)(Math.random() * fem2.nodes));
			tEnd = System.nanoTime();
			testRuns[2][nwf - 1] = (double)(tEnd - tStart) / 1000000.0;
			if (verbose)
				System.out.printf("factor " + fem2.nodeworkFactor + "add*N+delete*N: %.1f ms, random add/delete: %.1f ms\n",
						testRuns[0][nwf - 1], testRuns[2][nwf - 1]);
		}
		
		if (verbose) {
			XYDataset bestFitChartSet;
			String[] testNames = {"","","","","add*N+delete*N","add*N+delete*N / mem.usage","random add/delete"};
			int[] testCases = {0, 1, 2};
			bestFitChartSet = createStatisticSet(testRuns, testNames, testCases);
			XYLineChart_AWT bestFitChart = new XYLineChart_AWT(
					"nodeWork allocation factor", "exec. time", "nsecs/test", bestFitChartSet, 1024, 768);
			bestFitChart.pack();          
			RefineryUtilities.centerFrameOnScreen(bestFitChart);   
			bestFitChart.setVisible( true ); 
		}
		return 0;
	}

}
