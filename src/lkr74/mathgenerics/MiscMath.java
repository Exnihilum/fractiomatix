package lkr74.mathgenerics;

import org.jfree.data.xy.XYDataset;
import org.jfree.ui.RefineryUtilities;

import lkr74.mathgenerics.XYLineChart_AWT;




public class MiscMath {
	
	// the hat basis function ramps up to 1 and ramps down to 0 again, reaches zero on reaching span of hat basis
	public class HatBasis {
		// default span is 1
		private double span = 1, xtip = 0.5;
		
		// Validate & construct a hat basis method
		HatBasis(double span) {
			if (span == 0) throw new IllegalArgumentException("Hat basis f() span is 0");
			if (span < 0) throw new IllegalArgumentException("Hat basis f() span must be positive");
			this.span = span;
			this.xtip = span / 2;
		}
		
		// get the result from instantiated hat basis
		public double get(double x) {
			if (x < 0 || x > span)	return 0;
			if (x < xtip)			return x / xtip;
			return 2 - x / xtip;
			
		}
	}
	
	
	
	// RandFill class will get an (integer) range and will guarantee to return a random value within that range
	// but never the same value as previously, thus populating the range evenly with random values,
	// a useful class for iteratively test-filling arrays & matrices with values in a linearly asymptotic fashion
	// "sectors" holds integer pairs, the first is start index of a sector, the second an end index
	// there is always at least one occupied slot between each sector, every new occupied slot can create a new sector partition
	public class RandFill {
		private int sectorCnt = 1, slotCount;
		private int[] sectors;
		
		public RandFill(int range) {
			if (range < 1) throw new RuntimeException("RandFill(): Invalid range.");
			this.slotCount = range;
			sectors = new int[range + 1];		// holds definitions of sectors between occupied random value slots
			sectors[1] = range - 1;				// define the first default sector as covering entire range
		}
		
		public int remainingSlots() { return slotCount; }
		
		public int getRandom() {
			
			// zero sectors left means we've exhausted all random slots, return 0
			if (sectorCnt == 0) return 0;
			
			int rndSec = (int)(Math.random() * sectorCnt) * 2;			// select a random sector
			int secLen = sectors[rndSec + 1] - sectors[rndSec] + 1;		// get it's length
			int rnd;
			
			// did we find a single-element sector?
			if (secLen == 1) {
				sectorCnt--;
				slotCount--;
				rnd = sectors[rndSec];
				sectors[rndSec] = sectors[sectorCnt * 2];				// delete this 1-element sector by copying last-in-array sector over it
				sectors[rndSec + 1] = sectors[sectorCnt * 2 + 1];
				return rnd;												// return this random slot
			}
			// select a random slot within this sector
			rnd = (int)(Math.random() * secLen) + sectors[rndSec];
			if (rnd == sectors[rndSec])	sectors[rndSec]++;				// slot is at start of sector, increment sector start point
			else if (rnd == sectors[rndSec + 1]) sectors[rndSec + 1]--;	// slot is at end of sector, decrement sector end point
			else {														// slot is in middle of sector, partition into two sectors
				sectors[sectorCnt * 2] = rnd + 1;						// create new sector at end of sector array
				sectors[sectorCnt * 2 + 1] = sectors[rndSec + 1];		// let it start above the slot
				sectors[rndSec + 1] = rnd - 1;							// decrement current sector below the slot
				sectorCnt++;
			}
			slotCount--;
			return rnd;
		}
	}

	
	public static double[] solveCubicPolynomial(double A, double B, double C) {
		// max 4 values to return: result1, result1, complexpart(+/-)
		// if solution is noncomplex, complexpart will return zero
		if (A == 0) throw new RuntimeException("solveCubicPolynomial(): divide by zero.");

		double[] solution = new double[3];
		double rootpart = B * B - 4 * A * C, rootpart2 = Math.sqrt(Math.abs(rootpart)) / (2 * A);
		solution[0] = solution[1] = -B / (2 * A);
		// check if we'll have a complex solution
		if (rootpart < 0) {
			solution[2] = rootpart2;	// the complex part
			return solution;
		}
		solution[0] += rootpart2;
		solution[1] -= rootpart2;
		return solution;
	}
	
	
	
	// given a list of x,y data pairs, fit a linear function to them
	// returns in a new buffer the b0, b1 from linear function y = b0 + b1 * x
	public static double[] findLinearF(double[] xdata, double[] ydata) {

		if (xdata.length != ydata.length)
			throw new RuntimeException("findLinearF(): (x,y) value pair counts don't match.");
		int n = xdata.length;
		if (n <= 0)
			throw new RuntimeException("findLinearF(): invalid element count.");

		double sumxy = 0, sumx = 0, sumxsq = 0, sumy = 0;
		for (double x : xdata) { sumx += x; sumxsq += x * x; }
		for (double y : ydata) sumy += y;
		for (int i = 0; i < xdata.length; i++) sumxy += xdata[i] * ydata[i];

		double div = (sumxsq - sumx * sumx / n);
		if (div == 0) throw new RuntimeException("findLinearF(): divide by zero.");

		double[] b = new double[2];
		b[1] = (sumxy - (sumx * sumy) / n) / div;
		b[0] = (sumy - b[1] * sumx) / n;
		return b;
	}
	
	
	// works like findLinearF() except it calculates partial linear functions on sequential parts of the (x,y)-data
	// "plength" tells the length of every partial sequence for which to find a linear approximation
	// method returns (x0, b00, b10), (x1, b01, b11), ... etc sequentially in a flat array
	// method expects an increasing, SORTED array where x0 <= x1 <= ... <= xn
	public static double[] findLinearPartialF(double[] xdata, double[] ydata, int plength) {
		if (xdata.length != ydata.length || xdata.length == 0)
			throw new RuntimeException("findLinearPartialF(): (x,y) value pair counts invalid.");
		if (plength <= 0)
			throw new RuntimeException("findLinearPartialF(): invalid partial sequence length value supplied.");

		double sumxy, sumx, sumxsq, sumy, b1;
		int bstep = 0, n = xdata.length, npartials = n / plength + (n % plength != 0 ? 1 : 0);
		// this array holds (x, b0, b1) triples: centroid of sequence, bias, factor
		if (plength >= n) npartials = 1;
		double[] xb0b1array = new double[npartials * 3];		

		// for every partial sequence of the (x,y)-array
		for (int s = 0; s < npartials; s++) {
			
			int offset = s * plength, pl = plength;
			sumxy = 0; sumx = 0; sumxsq = 0; sumy = 0;
			
			// offset into and do summations over a partial sequence of the (x,y)-array
			for (int i = 0; i < plength; i++) {
				// the last partial sequence will obviously often be too short, stop iterating in that case
				if (n <= offset + i) { pl = i; break; }
				double x = xdata[offset + i], y = ydata[offset + i]; 
				sumx += x; sumxsq += x * x;
				sumy += y;
				sumxy += x * y;
			}

			double div = (sumxsq - sumx * sumx / (double)pl);
			if (div == 0) throw new RuntimeException("findLinearPartialF(): divide by zero.");

			// insert the centroid x between the start x and end x of this sequence
			xb0b1array[bstep++] = (xdata[offset] + xdata[offset + pl - 1]) / 2.0;
			
			// insert calculated (b0, b1) pair into the array of partial bias + factor values 
			xb0b1array[++bstep] = b1 = (sumxy - (sumx * sumy) / (double)pl) / div;
			xb0b1array[--bstep] = (sumy - b1 * sumx) / (double)pl;
			bstep += 2;
		}
		return xb0b1array;
	}
	
	
	
	// getInterpolationOfLinPartials() uses the set of partial linear functions gotten from findLinearPartialF()
	// to get a value interpolated from these functions, which lie in a series,
	// defined by the (x, b0, b1) center-of-partial + bias + factor values
	public static double getInterpolationOfLinPartials(double[] xb0b1array, double x) {
		
		int npartials = xb0b1array.length / 3;
		
		// if there's just one linear f() stored, return it's linear value
		if (npartials == 1)
			return xb0b1array[1] + xb0b1array[2] * x;

		// walk through the (x, b0, b1) triplet array of partial linear f() factors, the x-es expected to be sorted
		int t = 0;
		for (int i = 0; i < npartials - 1; i++, t += 3) {
			double xa = xb0b1array[t];
			
			// if x lies past this interpolation domain of the linear f() partials
			if (x >= xa) {
				double xb = xb0b1array[t + 3];
				// if x fits into this partial linear f()'s domain
				if (x < xb) {
					double itpfa = (xb - x) / (xb - xa), itpfb = (x - xa) / (xb - xa);
					return	(itpfa * xb0b1array[t + 1] + itpfb * xb0b1array[t + 4]) +	// = b0 interpolated
							(itpfa * xb0b1array[t + 2] + itpfb * xb0b1array[t + 5]) *	// = b1 interpolated
							x;
				}
			// we're before even the first linear f() domain, so retract the linear value of the first f()
			} else return xb0b1array[t + 1] + xb0b1array[t + 2] * x;
		}
		// we're past every linear f() domain, so continue the linear value of the last f()
		return xb0b1array[t + 1] + xb0b1array[t + 2] * x;
	}
	
	
	
	static double getLaplacian2ndORderSystem(double s, double sigma, double omegad) {
		
		if (s == 0 || omegad == 0) throw new RuntimeException("solveCubicPolynomial(): divide by zero.");
		
		double spsigma = s + sigma, ssigomeghyp = (spsigma * spsigma + omegad * omegad);
		if (ssigomeghyp == 0) throw new RuntimeException("solveCubicPolynomial(): divide by zero.");
		return (1 / s) - spsigma / ssigomeghyp - (sigma / omegad) * omegad / ssigomeghyp;
	}
	
	
	public static double triangle3DArea(double[] a, double[] b, double[] c, boolean root) {
		double x1 = b[0] - a[0], x2 = b[1] - a[1], x3 = b[2] - a[2];
		double y1 = c[0] - a[0], y2 = c[1] - a[1], y3 = c[2] - a[2];
		double sq1 = x2 * y3 - x3 * y2, sq2 = x3 * y1 - x1 * y3, sq3 = x1 * y2 - x2 * y1;
		if (root)	return 0.5 * Math.sqrt(sq1 * sq1 + sq2 * sq2 + sq3 * sq3);
		else		return sq1 * sq1 + sq2 * sq2 + sq3 * sq3;
	}
	
	public static double triangle3DArea2(double[] r1, double[] r2, double[] r3, boolean root) {
		double dr1 = r1[0] - r2[0], dr2 = r1[1] - r2[1], dr3 = r1[2] - r2[2];
		double a = Math.sqrt(dr1 * dr1 + dr2 * dr2 + dr3 * dr3);
		dr1 = r2[0] - r3[0]; dr2 = r2[1] - r3[1]; dr3 = r2[2] - r3[2];
		double b = Math.sqrt(dr1 * dr1 + dr2 * dr2 + dr3 * dr3);
		dr1 = r3[0] - r1[0]; dr2 = r3[1] - r1[1]; dr3 = r3[2] - r1[2];
		double c = Math.sqrt(dr1 * dr1 + dr2 * dr2 + dr3 * dr3);
		double s = (a + b + c) / 2;
		if (root)	return Math.sqrt(s * (s - a) * (s - b) * (s - c));
		else		return s * (s - a) * (s - b) * (s - c);
	}
	
	

	public static class ElevatorMovement {
		
		// ElevatorMovement class forms an elevation pattern for a constantly accelerating & deccelerating elevator.
		// Given a time it will return correct elevator position, with base position at zero.
		// Elevation pattern is enveloped by start acceleration, end acceleration and max attainable speed inbetween.
		private boolean notmoving = false;
		// set to defaults: d = 1m (elevation distance), v = 1m/s (max attained speed) and (a1, a2) = 1m/s2 (constant acceleration)
		private double d = 1, v = 1, a1 = 1, a2 = 1;
		// t = total expected ascend/descend time
		// t1 = time when speed v is attained with accel. a1
		// t2 = time when elevator must begin deccelerate with accel. a2
		private double t, t1, t2;
		// elevator height attained at times t1 & t2
		private double d1, d2;
		
		// Construct the elevator elevation case.
		// the Elevator class features recalculation of v if it isn't attainanle for a given d & a1 & a2
		public ElevatorMovement(double d, double v, double a1, double a2) {

			// Parameter sanity check, any of (d, v, a1, a2) = 0 will immobilise ElevatorMovement
			if (d == 0 || v == 0 || a1 == 0 || a2 == 0) { notmoving = true; return; }
			if (d < 0 || a1 < 0 || a2 < 0)
				throw new IllegalArgumentException("ElevatorMovement: distance or accelerations cannot be negative");
	
			this.d = d;
			this.a1 = a1;
			this.a2 = a2;
			double dd1 = v * v / (2 * a1);	// delta d1 is distance attained by reaching speed (0 -> v) at acceleration a1
			double dd2 = v * v / (2 * a2);	// delta d2 is distance attained by reaching speed (0 -> v) at acceleration a2

			// Check if speed v is attainable within elevator's ascend/descend envelope
			// by calculating delta d1 & d2 for a/decceleration phases a1 & a2
			if ((dd1 + dd2) > d) {
				// v unattainable, rescale dd1 & dd2 into elevator envelope
				double ddfac = 1 / (dd1 + dd2);
				dd1 *= d * ddfac;
				dd2 *= d * ddfac;
				// recalculate max attainable speed by elevator from dd1 & a1 perspective
				v = Math.sqrt(2 * a1 * dd1);
				System.out.printf("Elevator speed unattainable, recalculated to %.2f\n", v);
			}
			this.v = v;
			
			// get timepoints t1 and t2, will be important for getting elevation at arbitrary timepoint
			t1 = Math.sqrt(2 * dd1 / a1);
			t2 = Math.sqrt(2 * dd2 / a2);
			double t12 = (d - dd1 - dd2) / v;
			t = t1 + t12 + t2;
			t2 = t1 + t12;

			// get elevation points for t1 & t2
			d1 = dd1;
			d2 = d - dd2;
			System.out.printf("t = %.2f, t1 = %.2f, t2 = %.2f\n", t, t1, t2);
		}
		
		// get elevator distance at point t when ascending
		public double getAscendingDistance(double t) {
			
			if (notmoving || t <= 0) return 0;
			if (t >= this.t) return d;
			if (t < t1) return a1 * t * t / 2;
			if (t > t2) {
				double dt = this.t - t;
				return d - a2 * dt * dt / 2;
			}
			// the straight part of the elevation
			return d1 + v * (t - t1);
		}

		public double getDescendingDistance(double t) {
			return getAscendingDistance(this.t - t);
		}
		
		// calculates the time from the displacement of the elevator
		public double getTime(double d) {
			
			if (notmoving || d <= 0) return 0;
			if (d >= this.d) return t;
			if (d < d1) return Math.sqrt(2 * d / a1);
			if (d > d2) return t - Math.sqrt(2 * (this.d - d) / a2);
			return t1 + (d - d1) / v;
		}
	}
	
	

	public static void main(String[] args) {
		
		// set true to print out the simulated data chart
		boolean drawchart = true;
		double[] dlist, tlist;
		// tpad = flat timestretch padded before & after dataset span
		// dt = delta t (time step)
		double tpad = 0.2, dt = 0.02;
		
		double[][] tv = {{-1,0,0}, {-1,1,0}, {0,0,0}};
		double[][] tv2 = {{0,0,1}, {0,1,1}, {0,0,0}};
		
		RandFill rfill = new MiscMath().new RandFill(100);
		for (int i = 0; i < 110; i++)
			System.out.println(rfill.getRandom());
		
		System.out.println("3D Triangle area: " + triangle3DArea(tv[0], tv[1], tv[2], true));
		System.out.println("3D Triangle area: " + triangle3DArea(tv2[0], tv2[1], tv2[2], true));
		System.out.println("3D Triangle area: " + triangle3DArea2(tv[0], tv[1], tv[2], true));
		System.out.println("3D Triangle area: " + triangle3DArea2(tv2[0], tv2[1], tv2[2], true));

		System.out.println();
		double[] s1 = solveCubicPolynomial(1, 3f/2f, 10f/2f);
		if (s1[2] == 0)
				System.out.println("cubic solution: x1: " + s1[0] + ", x2: " + s1[1] + "\n");
		else	System.out.println("cubic complex solution: x1: " + s1[0] + ", x2: " + s1[1] + ", c: " + s1[2] + " * sqrt(-1)\n");
		
		ElevatorMovement elevator1m = new ElevatorMovement(100, 80, 1, 1);

		System.out.printf("Elevator profile: t:%.2f vmax:%.2f a1:%.2f a2:%.2f\n\n", elevator1m.t, elevator1m.v, elevator1m.a1, elevator1m.a2);
		
		// text-only putput:
		if (!drawchart) {
			for (double time=-tpad; time < elevator1m.t + tpad; time += dt)
				System.out.printf("Elevator: t:%.2f ascend:%.2f ascend-by-t:%.2f\n",
						time, elevator1m.getAscendingDistance(time), elevator1m.getTime(elevator1m.getAscendingDistance(time)));
		} else {
			
			// calculate necessary size of datalist for elevator displacement sampling, including time padding
			dlist = new double[(int)((elevator1m.t + tpad * 2) / dt)];
			// converserly, this is the time samples list
			tlist = new double[dlist.length];
			double t = -tpad;
			// samples of elevator displacement calculated & stored here
			for (int i = 0; i < dlist.length; i++) {
				// get displacement distance from a given time
				dlist[i] = elevator1m.getAscendingDistance(t);
				// get a time from a given displacement distance
				tlist[i] = elevator1m.getTime(dlist[i]);
				t += dt;
			}
			
			// sample the Laplacian 2nd Order System
			double step = 1f/20f, s = step;
			double[] rlist = new double[200];
			for (int i = 0; i < rlist.length; i++) {
				// sigma = rate of decay, omegad = oscillation freq.
				rlist[i] = getLaplacian2ndORderSystem(s, 0.9, 500);
				s += step;
			}

			
			// sample (x,y) values for best-fit linear function findLinearF()
			double[] xList = new double[40], yList = new double[40], bestFit, bestPartialFit;
			for (int x = 0; x < xList.length; x++) {
				xList[x] = x;
				// randomise some noisy function
				yList[x] = Math.sin(x/6.0) + Math.random() * 0.1;
				//yList[x] = 4 - Math.pow(x, 3.5) / 20f + Math.pow(x+5, 3) / 3f + Math.random()*400;
			}
			bestFit = findLinearF(xList, yList);
			
			
			// calculate the linear partials set for interpolation
			bestPartialFit = findLinearPartialF(xList, yList, 10);
			System.out.println("Best Partial Linear Fits:");
			for (int i = 0; i < bestPartialFit.length/3; i++)
				System.out.println("x: " + bestPartialFit[i*3] + " b0: " + bestPartialFit[i*3+1] + " b1: " + bestPartialFit[i*3+2]);

			
			// chart for testing elevator function
			XYDataset elevatorChartSet = XYLineChart_AWT.createElevatorDataset(dlist, tlist, -tpad, dt);
			XYLineChart_AWT chart = new XYLineChart_AWT("Elevator displacement profile", "time (s)", "displacement (m)", elevatorChartSet, 800, 600);
			chart.pack( );          
			RefineryUtilities.centerFrameOnScreen( chart );   
			chart.setVisible( true ); 

			// chart for testing best-fit function
			XYDataset bestFitChartSet = XYLineChart_AWT.createBestFitDataset(xList, yList, bestFit, bestPartialFit);
			XYLineChart_AWT bestFitChart = new XYLineChart_AWT("Best fit profile", "random data", "fitted linear", bestFitChartSet, 1024, 768);
			bestFitChart.pack( );          
			RefineryUtilities.centerFrameOnScreen(bestFitChart);   
			//bestFitChart.setVisible( true ); 

			// chart for testing Laplacian
//			XYDataset secondOrderSysChartSet = XYLineChart_AWT.create2ndOrderSystemDataset(rlist, step);
//			XYLineChart_AWT chart2 = new XYLineChart_AWT("2nd Order System Laplacian Response", secondOrderSysChartSet, 800, 600);
//			chart2.pack( );          
//			RefineryUtilities.centerFrameOnScreen( chart2 );          
			//chart2.setVisible( true );
			
		}
	}
}
