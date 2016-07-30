package lkr74.mathgenerics;

import java.awt.Color; 
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

public class XYLineChart_AWT extends ApplicationFrame 
{
	private static final long serialVersionUID = -3401891241593865248L;

	public XYLineChart_AWT( String title, String xtag, String ytag, XYDataset dataset, int width, int height)
	{
		super(title);
		JFreeChart xylineChart = ChartFactory.createXYLineChart(
				title, xtag, ytag, dataset, PlotOrientation.VERTICAL, true, true, false);

		ChartPanel chartPanel = new ChartPanel( xylineChart );
		chartPanel.setPreferredSize( new java.awt.Dimension( width , height ) );
		final XYPlot plot = xylineChart.getXYPlot();
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
		renderer.setSeriesPaint(0 , Color.BLACK );
		renderer.setSeriesPaint(1 , Color.BLUE );
		renderer.setSeriesPaint(2 , Color.RED );
		renderer.setSeriesPaint(3 , Color.GREEN );
		renderer.setSeriesStroke(0 , new BasicStroke(2));
		renderer.setSeriesStroke(1 , new BasicStroke(2));
		renderer.setSeriesStroke(2 , new BasicStroke(2));
		renderer.setSeriesStroke(3 , new BasicStroke(2));
		plot.setRenderer( renderer ); 
		setContentPane( chartPanel ); 
	}

	static XYDataset createElevatorDataset(double[] dlist, double[] tlist, double tstart, double dt)
	{
		final XYSeries elevator1 = new XYSeries("Elevator ascending (t->d)", false);
		final XYSeries elevator1t = new XYSeries("Elevator ascending (d->t)", false);
		int i = 0;
		for(double d: dlist) {
			elevator1.add(tstart, d);
			elevator1t.add(tlist[i++], d);
			//elevator1t.add(tstart, tlist[i++]);
			tstart += dt;
		}
		final XYSeriesCollection dataset = new XYSeriesCollection();          
		dataset.addSeries(elevator1);          
		dataset.addSeries(elevator1t);          
		return dataset;
	}


	static XYDataset create2ndOrderSystemDataset(double[] rlist, double step)
	{
		final XYSeries response = new XYSeries("System Laplacian response", false);
		double stepoffset = 0;
		for(double r: rlist) {
			response.add(stepoffset, r);
			stepoffset += step;
		}
		final XYSeriesCollection dataset = new XYSeriesCollection();          
		dataset.addSeries(response);          
		return dataset;
	}


	static XYDataset createBestFitDataset(double[] xList, double[] yList, double[] bestFit, double[] bestFitP)
	{
		final XYSeries dataSeries = new XYSeries("Datapoints", false);
		final XYSeries fitSeries = new XYSeries("Best linear fit to datapoints", false);
		final XYSeries fitPSeries = new XYSeries("Interpolation from a set of partial best-fit linears", false);
		int i = 0;
		
		// construct the datapoint series
		for(double x: xList)  dataSeries.add(x, yList[i++]);
		
		// construct the best-fit linear f() series
		fitSeries.add(xList[0], bestFit[0] + bestFit[1] * xList[0]);
		fitSeries.add(xList[xList.length-1], bestFit[0] + bestFit[1] * xList[xList.length-1]);
		
		// construct the interpolation of linear f() partials series, sample it at 100 points
		if (bestFitP.length < 4) {
			fitPSeries.add(xList[0], bestFitP[1] + bestFitP[2] * xList[0]);
			fitPSeries.add(xList[xList.length-1], bestFitP[1] + bestFitP[2] * xList[xList.length-1]);			
		} else {
			double dstep = (bestFitP[bestFitP.length - 3] - bestFitP[0] + 8) / 100.0;
			for (double d = bestFitP[0]-4; d < bestFitP[bestFitP.length-3] + 4; d += dstep)
				fitPSeries.add(d, MiscMath.getInterpolationOfLinPartials(bestFitP, d));
		}
		
		final XYSeriesCollection dataset = new XYSeriesCollection();          
		dataset.addSeries(dataSeries);          
		dataset.addSeries(fitSeries);          
		dataset.addSeries(fitPSeries);          
		return dataset;
	}
}