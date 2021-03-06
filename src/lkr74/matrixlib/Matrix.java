package lkr74.matrixlib;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.InvalidParameterException;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;


public class Matrix implements Cloneable {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			FLAT ARRAY MATRIX OPERATIONS																//
	//			Leonard Krylov 2016																			//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	public static final double ROUNDOFF_ERROR = 1e-8;
	static double pivotCritFactor = 0.5;						// pivotCritFactor * candidate pivot will elect candidate if > than current pivot
	static double pivotFailFactor = Matrix.ROUNDOFF_ERROR*100;	// the smallest value allowed for a pivot to categorise as failed
	static double pivotValid = Matrix.ROUNDOFF_ERROR*101;		// the smallest value allowed for a pivot to categorise as valid

	static final int NULL_MATRIX = 1;
	static final int COMPLEX_MATRIX = 1<<1;
	static final int SINGULAR_MATRIX = 1<<8;
	static final boolean REAL = false;
	static final boolean IMAGINARY = true;
	static final int L1_CACHE_BYTES = 32, DOUBLE_BYTES = 8;
	static final int CACHED_DOUBLES = L1_CACHE_BYTES / DOUBLE_BYTES;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			INSTANCE-LEVEL VALUES
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public String name;
	public int M, N; 								// number of rows & columns
	protected double[] data, idata;					// row-column matrix uses a flat array, data = real part, idata = imaginary part
	protected int status;							// status bits
	public double det = 0;							// last calculated determinant value
	public double norm1 = -1;						// last calculated norm-1 value

	// optimisation hierarchies & variables
	// binBitLength = no. of values gathered in binBit, procNum = number of concurrent processors available
	protected boolean rowAspectSparsest = true;		// indicates whether row or column aspect is most sparse for determinant analysis
	public int halfBandwidth = -1;					// tells how far a sparse matrix deviates from it's diagonal band, p = |i - j|
	protected int detSign = 1;						// tells sign of the determinant, according to heuristics of analyseRowsColumns() 
	public int nNZ = 0;								// tells how many nonzeroes matrix contains
	protected int[][] analysis = null;				// used by analyseRowsColumns() for determinant finding heuristics
	protected int[][] mutator = null;				// used to index mutated rows & columns for certain algorithms
	protected BinBitImage bitImage;
	
	// multithreading variables
	static ExecutorService executor;				// thread pool for matrix ops
	static List<Future<Double>> futureDoubleList;	// list of return doubles from Callable threads
	//protected boolean threaded = false;			// indicates if matrix should/can be run in multithreaded mode	
	//static int taskNum, procNum;
	static Thread taskList[];
	volatile double vdet = 0;						// threadsafe determinant accumulator for multithreaded methods
	
	// global variables
	protected static int DEBUG_LEVEL = 2;
	protected volatile static int nameCount = 1;
	// debugging global variables
	protected static volatile int detL_DEBUG;
	protected static volatile int mulSW_DEBUG;
	protected static volatile int mulFlops_DEBUG, mulAdops_DEBUG;
	protected static volatile int mulFlopsSW_DEBUG, mulAdopsSW_DEBUG;
	
	// skeleton instantiator of a matrix
	public Matrix(String name, int r, int c) {
		// initialise status flags, all bits = 0 means an empty matrix of type Matrix
		if (r < 1 || c < 1) throw new InvalidParameterException("Matrix(): Illegal matrix dimensions.");
		this.name = name + nameCount++;
		this.M = r;
		this.N = c;
		mutator = new int[2][];
		setNull();
	}

	//	instantiates Matrix with data of choice
	public Matrix(String name, int r, int c, Type type, double v) {
		// initialise status flags, all bits = 0 means an empty matrix of type Matrix
		if (r < 1 || c < 1) throw new InvalidParameterException("Matrix(): Illegal matrix dimensions.");
		this.name = name + nameCount++;
		this.M = r;
		this.N = c;
		this.generateData(type, v);
		if (type == Type.Null_Complex) { setNull(); setComplex(); }
		if (type == Type.Null) setNull();
		//for (int i = 0; i < r; i++) mutator[i] = i;
		bitImage = new BinBitImage(this);
		mutator = new int[2][];
	}
	
	//	instantiates a matrix with a provided dataset, cloning the dataset into this matrix
	public Matrix(String name, int M, int N, double[] data, double[] idata) {
		if (M < 1 || N < 1)							throw new RuntimeException("Matrix(): Illegal matrix dimensions.");
		if (data.length < M * N) 					throw new RuntimeException("Matrix(): Illegal data array.");
		if (idata != null && idata.length < M * N) 	throw new RuntimeException("Matrix(): Illegal complex data array.");
		this.name = name + nameCount++;
		this.M = M;
		this.N = N;
		this.putData(data, idata);
		bitImage = new BinBitImage(this);
		status = 0;
		mutator = new int[2][];
		if (idata != null) setComplex();
	}
	
	// instantiates a Matrix from a FrontalMatrix
	public Matrix(String name, FrontalMatrix fm) {
		M = fm.r; N = fm.c;
		this.name = name;
		
		// copy over the data from frontal matrix, taking care that the frontal can have extra data padding (cTot > c & rTot > r)
		data = new double[fm.r * fm.c];
		for (int i = 0; i < fm.r; i++)
			for (int j = 0; j < fm.c; j++)	data[i * N + j] = fm.data[i * fm.cTot + j];
		
		// convert the frontal's inherent column/row indexes into the mutator arrays of Matrix class
		// note: there might be no permutations at all, in which case we have a linearly increasing indexation
		mutator = new int[2][];
		mutator[0] = new int[fm.r]; mutator[1] = new int[fm.c];
		for (int i = fm.cTot, iEnd = i + fm.r, i2 = 0; i < iEnd; i++, i2++) mutator[0][i2] = fm.idxCR[i];
		for (int j = 0; j < fm.c; j++) mutator[1][j] = fm.idxCR[j];
		
		bitImage = new BinBitImage(this);
	}


/*	
	// initialise Executor (left to the user, since one can't get Runtime's processor count during class loading)
	public static void initExecutor(int processors) {
        // define the Callable executor pool of reusable threads
		if (processors == 0) 	executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		else					executor = Executors.newFixedThreadPool(processors);
		procNum = (processors == 0 ? Runtime.getRuntime().availableProcessors() : processors);
  		// define a Future list of Callable return values constituting doubles (it will expand dynamically)
		futureDoubleList = new ArrayList<Future<Double>>();	
	}

	// initialise Executor (left to the user, since one can't get Runtime's processor count during class loading)
	public static void initTaskList(int processors) {
        // define the Callable executor pool of reusable threads
		if (processors == 0) 	taskList = new Thread[Runtime.getRuntime().availableProcessors()];
		else					taskList = new Thread[processors];
		procNum = (processors == 0 ? Runtime.getRuntime().availableProcessors() : processors);
		taskNum = 0;
	}

	
	public static void finishTaskList() {
		for (Thread task: taskList) {
			if (task == null) continue;
			try { task.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		taskNum = 0;
	}
*/


	public double valueOf(int r, int c) { return data[r * N + c]; }
	
	public double[] valueOfC(int r, int c) {
		if (data == null || idata == null) throw new RuntimeException("Matrix.valueOfC(): Unallocated datafields.");
		double[] iv = new double[2];
		int rNc = r * N + c;
		iv[0] = data[rNc]; iv[1] = idata[rNc];
		return iv;
	}

	public void valueTo(int r, int c, double v) {
		if (data == null) throw new RuntimeException("Matrix.valueTo(): Unallocated datafield.");
		data[r * N + c] = v;
		if (!nearZero(v)) {
			bitImage.setBit(r, c);					// set the corresponding bit in bitImage if value is nonzero
			clearNull();
			int width = r - c;						// readjust half bandwidth if value coordinates are "sticking out" from diagonal
			if (width < 0) width = -width;
			if (halfBandwidth < width) halfBandwidth = width;
		}
	}
	public void valueToC(int r, int c, double[] iv) {
		if (data == null || idata == null) throw new RuntimeException("Matrix.valueToC(): Unallocated datafields.");
		int rNc = r * N + c;
		data[rNc] = iv[0]; idata[rNc] = iv[1];
		if (!nearZero(iv[0]) || !nearZero(iv[1])) {
			bitImage.setBit(r, c);					// set the corresponding bit in bitImage if value is nonzero
			clearNull();
			int width = r - c;						// readjust half bandwidth if value coordinates are "sticking out" from diagonal
			if (width < 0) width = -width;
			if (halfBandwidth < width) halfBandwidth = width;
		}
	}

	
	// generate identity matrix with scalar v in the diagonal
	public Matrix identity(int s, double v) { Matrix M = new Matrix("I", s, s); M.generateData(Type.Identity, v); return M; }
	
	
	// calculates 1-norm of a matrix
	private double norm1() {
		int MN = M * N;
		double sum;
		for (int j = 0; j < N; j++) {
			sum = 0;
			for (int i = j; i < MN; i += N) { double v = data[i]; sum += v < 0 ? -v : v; }
			if (sum > norm1) norm1 = sum;
		}
		return norm1;
	}

	
	
	// calculates Frobenius norm of a matrix
	private double normFrobenius() {
		int MN = M * N;
		double sum = 0;
		for (int i = 0; i < MN; i++) { double v = data[i]; sum += v * v; }
		return Math.sqrt(sum);
	}


	
	// polymorphic matrix centering method will numerically center the matrix, columnwise
	public Matrix center(boolean copy) {
		Matrix Ac = this;
		if (copy) Ac = clone();
		double[][] dataSet = Ac.getDataRef();
		double[] dataAc = dataSet[0];
		int MN = M * N;
		for (int j = 0; j < N; j++) {
			double avg = 0;
			for (int iNj = j; iNj < MN; iNj += N) avg += dataAc[iNj];
			avg /= (double)M;
			for (int iNj = j; iNj < MN; iNj += N) dataAc[iNj] -= avg;
		}
		dataAc = dataSet[1];
		if (isComplex()) {
			double[] idataAc = Ac.idata;
			for (int j = 0; j < N; j++) {
				double avg = 0;
				for (int iNj = j; iNj < MN; iNj += N) avg += idataAc[iNj];
				avg /= (double)M;
				for (int iNj = j; iNj < MN; iNj += N) idataAc[iNj] -= avg;
			}
		}
		Ac.halfBandwidth = -1;
		Ac.bitImage.make();
		return Ac;
	}

	
	// generates different types of matrices, a supplied value can be used to alter the matrix modulo
	public void generateData(Type type, double v) {
		if (M < 1 || N < 1)
			throw new InvalidParameterException("Matrix.generateData(): Illegal matrix dimensions.");
		
		double[] data = null, idata = null;
		boolean findHalfbandWidth = false;
		
		switch (type) {
			case Null_Complex:
				idata = new double[M * N];
				setComplex();
			case Null:
				data = new double[M * N];
				halfBandwidth = -1;
				setNull();
				break;
			case Unit:
				data = new double[M * N];
				for (int i = 0; i < M * N; i++) data[i] = 1;
				if (M != N) halfBandwidth = -1; else halfBandwidth = M - 1;
				break;
			case Identity:
				if (M != N) 	throw new InvalidParameterException("Matrix.generateData(): Identity matrix must be square.");
				data = new double[M * N];
				for (int i = 0; i < N; i++) data[i * N + i] = v;
				halfBandwidth = 0;
				break;
			case Random_Complex:
				idata = new double[M * N];
				for (int i = 0; i < M * N; i++) idata[i] = Math.random();
				setComplex();
			case Random:
				data = new double[M * N];
				for (int i = 0; i < M * N; i++) data[i] = Math.random();
				if (M != N) halfBandwidth = -1; else halfBandwidth = M - 1;
				break;
			case Centering:
				if (M != N) 	throw new InvalidParameterException("Matrix.generateData(): Centering matrix must be square.");
				data = new double[M * N];
				for (int i = 0; i < M; i++)
					for (int j = M * N, jEnd = j + N; j < jEnd; j++) {
						double d = 1.0 / (double)M;
						if (i == j)		data[j] = 1.0 - d;
						else			data[j] = - d;
					}
				halfBandwidth = -1;
				break;
			case UpperHessenberg:
				data = new double[M * N];
				int m = (int)v;
				for (int i = 0; i < M; i++)
					for (int j = 0; j < N; j++)
						if (i - 1 <= j && j <= m + i - 1)
							data[i * N + j] = 1;
				if (M != N) halfBandwidth = -1; else findHalfbandWidth = true;
				break;
			default:
				throw new InvalidParameterException("Matrix.generateData(): Illegal matrix type.");
		}
		putDataRef(data, idata);
		if (findHalfbandWidth) getHalfBandwidth();
	}

	
	public Matrix rescale(int Mi, int Ni, int Mo, int No, boolean doBitImage) {
		
		if (Mi > Mo || Ni > No) throw new InvalidParameterException("Matrix.rescale(): Invalid data subset size.");
		int newM = Mo - Mi, newN = No - Ni;
		// if rescaled field ends up outside the matrix datafield, set matrix to Null matrix
		
		Matrix A = new Matrix("R", newM, newN, Matrix.Type.Null, 1);
		A.name = (DEBUG_LEVEL > 1 ? name + "(r:" + Mi + "," + Ni + "," + Mo + "," + No + ")" : name + "(r)");

		double[] newdata = A.data;

		for (int nr = Mi, r = 0; nr < Mo; nr++, r++) {
			// store a row only if indexing stays inside row bounds of original matrix
			if (nr >= 0 && nr < M ) {
				int nrofs = r * newN, rofs = nr * N;
				for (int nc = Ni, c = 0; nc < No; nc++, c++) {
					// store a value only if indexing stays inside column bounds of original matrix
					if (nc >= 0 && nc < N)
						newdata[nrofs + c] = data[rofs + nc];
				}
			}
		}
		if (!(Mi >= M || Ni >= N || Mo < 0 || No < 0)) A.clearNull();
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("Rescaling matrix " + A.name);
			System.out.println(A.toString());
		}
		
		if (doBitImage) A.bitImage = new BinBitImage(A);
		return A;
	}

	
	
	// the Gram-Schmidt orthonormalisation of column vectors constituting a matrix
	// method also finds out if vectors are linearly independent
	// if copy = false, use supplied matrix-form vectors for calculations
	public Matrix orthogonalise(boolean copy) {
		
		Matrix An = this;
		// if copy = false, An will be used to store intermediate q-hat vectors
		if (copy) An = this.clone();
		else An.halfBandwidth = -1;
		double r11 = An.absLength(0);										// r11 := ||x1||
		if (nearZero(r11)) return null;										// null vector found, abort

		int M = An.M, N = An.N;
		Matrix Q = new Matrix("Q", M, N, Matrix.Type.Null, 1);
		double[] dataAn = An.data, dataQ = Q.data;
		r11 = 1.0 / r11;
		for (int r = 0, rN = 0; r < M; r++, rN += N)
			dataQ[rN] = dataAn[rN] * r11;									// q1 := x1 / r11
		
		// main loop steps column by column
		// note: should be optimised for cache efficiency to run along rows, not columns
		for (int j = 1; j < N; j++) {
			for (int i = 0; i < j; i++) {
				// use vector j in An as q^, multiply by columns
				double rij = multiplyVectors(An, Q, j, i, false);			// rij := (q^,qi)
				for (int k = 0, kNi = i, kNj = j; k < M; k++, kNi += N, kNj += N)
					dataAn[kNj] -= rij * dataQ[kNi];						// q^ := q^ - rij*qi
			}
			double rjj = An.absLength(j);									// rjj := ||q^||
			if (nearZero(rjj)) return null;									// null vector found, abort
			rjj = 1.0 / rjj;
			for (int r = 0, rN = j; r < M; r++, rN += N)
				dataQ[rN] = dataAn[rN] * rjj;								// qj := q^ / rjj
		}

		if (DEBUG_LEVEL > 1) {
			System.out.println("Gram-Schmidt orthonormalisation" + (copy ? ": " : "(input matrix wasn't saved): "));
			System.out.println(Q.toString());
		}
		return Q;
	}
	
	
	
	// Cholesky factoriser conditions symmetric and positive definite matrices for fast O(4*p*n) linear system solving
	public Matrix factoriseCholesky() {
		
		if (M < 1 || N < 1 || M != N)
			throw new InvalidParameterException("Matrix.factoriseCholesky(): Invalid matrix dimensions.");

		Matrix R = new Matrix("R", M, N, Matrix.Type.Null, 1);	// factorisation goes into this matrix
		double[] dataR = R.data;
		
		int p = 0;
		//if (halfBandwidth < 0) p = getHalfBandwidth();
		
		for (int i = 0; i < M; i++) {
			
			double rki = 0;
			int im1 = i - 1, kStart = i - p;
			if (kStart < 0) kStart = 0;
			for (int k = kStart * N + i, k1 = kStart; k1 <= im1; k += N, k1++)
				rki += dataR[k] * dataR[k];

			int iNi = i * N + i;
			dataR[iNi] = Math.sqrt(data[iNi] - rki);
			
			int jEnd = i + p;
			if (jEnd >= N) jEnd = N - 1;
			if (nearZero(dataR[iNi])) return null;					// divide by zero, zeroes on diagonal not allowed
			double riiD = 1.0 / dataR[iNi];
			
			for (int j = i + 1; j <= jEnd; j++) {
				double rkirkj = 0;
				for (int kN = kStart * N, k = kStart; k <= im1; k++, kN += N)
					rkirkj += dataR[kN + i] * dataR[kN + j];
				int iNj = i * N + j;
				dataR[iNj] = riiD * (data[iNj] - rkirkj);
			}				
		}
		if (DEBUG_LEVEL > 1) {
			System.out.println("Cholesky factorisation:");
			System.out.println(R.toString());
		}
		return R;
	}

	
	
	// given a Gram-Schmidt orthonormalised matrix Q and original matrix A
	// returns the matrix R, such that QR = A
	public Matrix decomposeQR(Matrix Q) {
	
		if (M < 1 || N < 1 || !(N==Q.M && M==Q.N))
			throw new InvalidParameterException("Matrix.decomposeQR(): Invalid matrix dimensions.");

		Matrix R = new Matrix("R", Q.M, N, Matrix.Type.Null, 1);	// matrix R will hold dot products of matrix elements in A & Q
		double[] dataR = R.data;
		
		for (int i = 0; i < Q.M; i++) {
			int iN = i * N;
			// matrix D will be upper triangular
			for (int j = i, iNj = iN + j; j < N; j++, iNj++)
				dataR[iNj] =  multiplyVectors(this, Q, j, i, false);
		}
		if (DEBUG_LEVEL > 1) {
			System.out.println("QR decomposition:");
			System.out.println(Q.toString());
			System.out.println(R.toString());
		}
		return R;
	}

	
	
	// reduces matrix to Householder form
	public Matrix reduceHouseholder(boolean copy) {

		if (M < 1 || N < 1) throw new InvalidParameterException("Matrix.reduceHouseholder(): Invalid matrix dimensions.");
		if (M != N) throw new InvalidParameterException("Matrix.reduceHouseholder(): Matrix not square.");

		DEBUG_LEVEL--;
		Matrix A = this;
		if (copy) A = this.clone();
		Matrix I = new Matrix("I", M, N, Matrix.Type.Identity, 1);
		Matrix v = new Matrix("v", M, 1, Matrix.Type.Null, 1), V;
		Matrix vT = new Matrix("vT", 1, N, Matrix.Type.Null, 1);
		double[] dataA = A.data;
		double a;
		
		// proceed columnwise
		for (int k = 0, MN = M * N; k < N - 1; k++) {
			double colSum = 0;
			int j1 = N * (k + 1) + k;
			for (int j = j1; j < MN; j += N) {
				a = dataA[j];
				colSum += a * a;					
			}
			a = dataA[j1];
			double alpha = (a < 0 ? 1.0 : -1.0) * Math.sqrt(colSum);	// alpha = -sign(a(k+1,k) * sqrt(sum(a(j,k)^2)
			double r = 0.5 / Math.sqrt(0.5 * alpha * (alpha - a));
			
			for (int kv = 0; kv <= k; kv++)								// set up v vector that will be multiplied up to a matrix
				v.data[kv] = vT.data[kv] = 0;
			v.data[k+1] = vT.data[k+1] = (a - alpha) * r;
			for (int j = k + 2, kvN = j * N + k; j < N; j++, kvN += N)
				v.data[j] = vT.data[j] = dataA[kvN] * r;
			
			V = v.multiply(vT).multiply(2, false);						// 2*v.v^T -> V
			Matrix P = I.subtract(V, true);								// P = I - 2*v.v^T
			A = P.multiply(A).multiply(P);								// A(k+1) = P(k).A(k).P(k)
			dataA = A.data;
		}
		
		DEBUG_LEVEL++;
		if (DEBUG_LEVEL > 1) {
			System.out.println("Householder reduction:");
			System.out.println(A.toString());
		}
		A.halfBandwidth = -1;
		return A;
	}
	
	

	// finds out number of zeroes in each row & column, sums them up into an analysis array
	public void analyse(boolean countNonZeroes) {
		
		if (M < 2 || N < 2) {
			if (M < 1 || N < 1) throw new InvalidParameterException("Matrix.analyse(): Invalid matrix dimensions.");
			return;		// 1x1 matrix, nothing to analyse
		}

		double[][] dataSet = getDataRef();
		double[] data = dataSet[0], idata = dataSet[1];		// get real & imaginary parts of data
		detSign = 1;
		rowAspectSparsest = true;
		int signR = 1, signC = 1;
		// zRCount & zCCount store number of zeroes in the rows and columns
		if (analysis == null) analysis = new int[4][];
		analysis[0] = new int[M];			// holds zeroes counts for each row
		analysis[1] = new int[N];			// holds zeroes counts for each column
		analysis[2] = new int[M];			// holds indexes into zCount[0], sorted according to most zeroes
		analysis[3] = new int[N];			// holds indexes into zCount[1], sorted according to most zeroes
		
		int[] zRCount = analysis[0], zCCount = analysis[1], zRIdx = analysis[2], zCIdx = analysis[3];
		// initialise the row & column zero-count indexation arrays
		for (int i = 0, iend = zRIdx.length; i < iend; i++) zRIdx[i] = i;
		for (int i = 0, iend = zCIdx.length; i < iend; i++) zCIdx[i] = i;
		
		nNZ = 0;							// reset nonzero counter, we're resumming nonzeroes
		for (int i = 0; i < M; i++) {
			int iN = N * i;
			if (idata != null) {
				for (int j = 0; j < N; j++)
					if (nearZero(data[iN + j]) && nearZero(idata[iN + j])) { zRCount[i]++; zCCount[j]++; }
					else nNZ++;
			} else {
				for (int j = 0; j < N; j++)
					if (nearZero(data[iN + j])) { zRCount[i]++; zCCount[j]++; }
					else nNZ++;
			}
		}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("Matrix.analyse() of " + this.name + ":");
			System.out.println("non-zeroes: " + nNZ + ", percent: " + (100f * (float)nNZ / (M*N)) + "%");
		}
		if (countNonZeroes) return;				// if only nonzero count requested, return here
		
		// do bubble sort of row indexes array (highest first), tracking the sign value
		for (boolean sorted = false; !sorted;) {
			sorted = true;
			for (int i = 0, iend = zRCount.length - 1; i < iend; i++)
				if (zRCount[zRIdx[i]] < zRCount[zRIdx[i + 1]]) {
					int temp = zRIdx[i]; 		// swap around indexes
					zRIdx[i] = zRIdx[i + 1];
					zRIdx[i + 1] = temp;
					signR = -signR;				// every row swap changes sign of determinant
					sorted = false;
				}
		}
		// do bubble sort of column indexes array (highest first), tracking the sign value
		for (boolean sorted = false; !sorted;) {
			sorted = true;
			for (int i = 0, iend = zCCount.length - 1; i < iend; i++)
				if (zCCount[zCIdx[i]] < zCCount[zCIdx[i + 1]]) {
					int temp = zCIdx[i]; 		// swap around indexes
					zCIdx[i] = zCIdx[i + 1];
					zCIdx[i + 1] = temp;
					signC = -signC;				// every column swap changes sign of determinant
					sorted = false;
				}
		}
		// find out which aspect (row vs column) gave most zeroes
		if (zRCount[zRIdx[0]] > zCCount[zCIdx[0]]) detSign = signR;
		// row aspect gave most zeroes, swap around array indexes
		else	{
			detSign = signC;
			rowAspectSparsest = false;
			analysis[0] = analysis[1]; analysis[1] = zRCount;
			analysis[2] = analysis[3]; analysis[3] = zRIdx;
		}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("matrix " + this.name + "'s row & column analysis:");
			System.out.println(Arrays.toString(analysis[0]));
			System.out.println(Arrays.toString(analysis[1]));
			System.out.println(Arrays.toString(analysis[2]));
			System.out.println(Arrays.toString(analysis[3]) + "\n");
		}
	}
	
	
	
	// swap two rows of matrix
	public void swapRows(int r1, int r2) {
		double temp;
		for (int i = 0, or1 = r1 * N, or2 = r2 * N; i < N; i++, or1++, or2++) {
			temp = data[or1]; data[or1] = data[or2]; data[or2] = temp;
			if (isComplex()) { temp = idata[or1]; idata[or1] = idata[or2]; idata[or2] = temp; }
		}
		det = -det;						// determinant is inverted by row swapping
		if (DEBUG_LEVEL > 2) System.out.println("Matrix.swapRows:\n" + this.toString());
		if (bitImage != null)
			bitImage.swapRows(r1, r2);		// swap bitImage rows
		halfBandwidth = -1;
	}


	
	
	// transpose the matrix, returning data within same matrix or in a copied matrix
	public Matrix transpose(boolean copy) {
		
		if (M < 1 || N < 1) throw new InvalidParameterException("Matrix.transpose(): Invalid matrix dimensions.");

		Matrix T = this;
		
		// given a non-square matrix, reallocate new data field, as values cannot be just exchanged across diagonal
		if (M != N) {
			if (copy) T = new Matrix(name + "^T", N, M);
			else name = name + "^T"; 
					
			double[] newdata = new double[N * M], inewdata = null;
			if (isComplex()) inewdata = new double[N * M];
			
			for (int i = 0; i < M; i++) {
				int iN = i * N;
				for (int j = 0, jMi = i, iNj = iN; j < N; j++, jMi += M, iNj++)
					newdata[jMi] = data[iNj];
				if (isComplex())
					for (int j = 0, jMi = i, iNj = iN; j < N; j++, jMi += M, iNj++) inewdata[jMi] = idata[iNj];
			}
			
			T.data = newdata;
			T.status = status;
			T.idata = inewdata;
			if (!copy) { int temp = T.M; T.M = T.N; T.N = temp; } // swap row & column count if we're not making a new matrix
			T.bitImage = new BinBitImage(T);
		} else {
			if (copy) { T = this.clone(); T.name = T.name + "^T"; }
			else name = name + "^T";
			double[] dataT = T.data, idataT = T.idata;
			
			// the approach for a square matrix
			double temp;
			for (int i = 0; i < M; i++) {
				int iN = i * N;
				for (int j = i + 1, jMpi, iNj = iN + j; j < N; j++, iNj++) {
					jMpi = j * M + i;
					temp = dataT[jMpi]; dataT[jMpi] = dataT[iNj]; dataT[iNj] = temp;
					if (isComplex()) { temp = idataT[jMpi]; idataT[jMpi] = idataT[iNj]; idataT[iNj] = temp; }
					T.bitImage.transposeBit(i, j);
				}
			}		
		}
		if (DEBUG_LEVEL > 1) System.out.println(T.toString());
		return T;
	}
	
	
	// normalises the matrix/vector/transposevector A, for a matrix the columns are normalised individually
	public Matrix normalise(boolean copy) {

		if (M < 1 || N < 1) throw new InvalidParameterException("Matrix.normalise(): Invalid matrix dimensions.");
		Matrix C;
		double[] dataC = null;
		if (copy) C = this.clone(); else C = this;
		double[][] dataSet = C.getDataRef();
		dataC = dataSet[0];
		int MN = M * N;

		if (M > 1)
			// if we're not dealing with a transposed vector
			for (int c = 0; c < N; c++) {
				double norm = 0;
				for (int r = c; r < MN; r+=N)
					norm += dataC[r] * dataC[r];							// Euclidean norm ||A||2
				norm = 1.0 / Math.sqrt(norm);
				for (int r = c; r < MN; r += N) dataC[r] *= norm;			// normalise
			}
		else {
			double norm = 0;
			for (int c = 0; c < N; c++)
				norm += dataC[c] * dataC[c];								// Euclidean norm ||A||2

			norm = 1.0 / Math.sqrt(norm);
			for (int c = 0; c < N; c++) dataC[c] *= norm;					// normalise
		}
		
		if (DEBUG_LEVEL > 1) System.out.println(this.toString());
		return C;
	}

	
	public Matrix add(Matrix B, boolean copy) { return addSub(B, false, copy); }
	public Matrix subtract(Matrix B, boolean copy) { return addSub(B, true, copy); }
	
	// polymorphic addSub will do A+B or A-B returning resultant C of matrix type A, or add/subtract directly to A if clone = false
	// the data from both matrices are gotten through the polymorphic method getDataRef()
	public Matrix addSub(Matrix B, boolean subtract, boolean copy) {
		
		if (B.M != M || B.N != N) throw new InvalidParameterException("Matrix.addSub(): Nonmatching matrix dimensions.");
		
		// bring in data from input matrices to the superclass format for adding
		double[][] dataSet = getDataRef();
		double[] dataA = dataSet[0], idataA = dataSet[1];
		dataSet = B.getDataRef();
		double[] dataB = dataSet[0], idataB = dataSet[1], dataC, idataC;
		Matrix C;
		if (copy) {
			C = new Matrix("C", M, N, Matrix.Type.Null, 1);
			dataC = C.data; idataC = C.idata;
		} else {
			C = this; 
			dataC = dataA; idataC = idataA;
		}
		
		C.name = (Matrix.DEBUG_LEVEL > 1 ? "(" + name + (subtract?"-":"+") + B.name + ")" : (subtract?"D":"S") + nameCount++);
		
		if (subtract)	for (int i = 0, MN = M * N; i < MN; i++) dataC[i] = dataA[i] - dataB[i];
		else			for (int i = 0, MN = M * N; i < MN; i++) dataC[i] = dataA[i] + dataB[i];
		// case tree takes care of situation when any of dataA/dataB/dataC can be null array
		if (idataB != null) {
			if (idataA != null) {
				// there is imaginary data in A and B, make sure C has buffer allocated
				if (idataC == null) { C.idata = new double[M * N]; C.setComplex(); }
				if (subtract)	for (int i = 0, MN = M * N; i < MN; i++) idataC[i] = idataA[i] - idataB[i];
				else			for (int i = 0, MN = M * N; i < MN; i++) idataC[i] = idataA[i] + idataB[i];
			} else idataC = idataB.clone();
		} else
			if (idataA != null && idataC != idataA) idataC = idataA.clone();
		
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		clearNull();
		// add/sub. simply increases half bandwidth to the bandwidth of the matrix of it's largest extent
		C.halfBandwidth = halfBandwidth > B.halfBandwidth ? halfBandwidth : B.halfBandwidth;
		C.bitImage.make();
		return C;		
	}

	
	
	// get absolute length for a vector or transposevector or a column vector in a matrix
	// "i" tells which column in a matrix to calculate length for
	public double absLength(int i) {
		double sq, vnorm = 0;
		boolean store = false;											// store length in determinant if we're deaing with a vector-form matrix
		int mN, vLen;
		if (N == 1)			{ i = 0; mN = 1; vLen = M; store = true; }	// column vector?
		else if (M == 1)	{ i = 0; mN = 1; vLen = N; store = true; }	// transposed vector?
		else {
			if (i >= N) 	throw new InvalidParameterException("Matrix.absLength(): Column index out of bounds.");
			mN = N; vLen = M;											// column vector from a matrix?
		}
		
		if (isComplex())	for (int c = 0, o = i; c < vLen; c++, o += mN) { sq = data[o] + idata[o]; vnorm += sq * sq; }
		else				for (int c = 0, o = i; c < vLen; c++, o += mN) { sq = data[o]; vnorm += sq * sq; }
		vnorm = Math.sqrt(vnorm);
		if (store) det = vnorm;
		return vnorm;
	}

	
	
	// takes two matrix-form vectors and returns an inner product value
	// ia & ib tell which two column/row vectors to multiply from the two matrices supplied
	// byRows = true tells to multiply row vectors, otherwise multiplication by columns happens
	public static double multiplyVectors(Matrix A, Matrix B, int ia, int ib, boolean byRows) {
		
		int aM = A.M, aN = A.N, bM = B.M, bN = B.N;
		if (aM < 1 || aN < 1 || bM < 1 || bN < 1)
			throw new InvalidParameterException("Matrix.multiplyVectors(): Invalid vector dimensions.");	
		
		double p = 0;
		double[] dataA = A.data, dataB = B.data;
		
		// case 1: two vector-form matrices, either one can be in row/column form
		if ((aM == 1 && (bM == 1 && aN == bN || bN == 1 && aN == bM)) ||
			(aN == 1 && (bM == 1 && aM == bN || bN == 1 && aM == bM))) {
			for (int i = 0, N = aN; i < N; i++) p += A.data[i] * B.data[i];
			return p;
		}
		
		// case 2: two matrices with arbitrary numbers of column vectors of same length
		if (aM == bM) {
			if (byRows) {
				if (ia >= aM || ib >= bM) throw new InvalidParameterException("Matrix.multiplyVectors(): Vector indexes out of bounds.");

				for (int i = bN, ca = ia * aN, cb = ib * aN; i > 0; i--, ca++, cb++)
					p += dataA[ca] * dataB[cb];
			} else {
				if (ia >= aN || ib >= bN) throw new InvalidParameterException("Matrix.multiplyVectors(): Vector indexes out of bounds.");

				for (int i = bM, ca = ia, cb = ib; i > 0; i--, ca += aN, cb += bN)
					p += dataA[ca] * dataB[cb];
			}
			return p;
		}
		throw new InvalidParameterException("Matrix.multiplyVectors(): Nonmatching matrix-form vector dimensions.");
	}
	
	
	
	// polymorphic matrix product will do AB = C, with zero-ignoring optimisation for sparse matrices
	// this is a regular nested loop implementation, with the cache-friendliness of iterating over neighbouring elements in the innermost loop
	public Matrix multiply(Matrix B) {

		int aM = M, aN = N, bM = B.M, bN = B.N;
		if (aM < 1 || aN < 1 || bM < 1 || bN < 1) 	throw new InvalidParameterException("Matrix.multiply(): Invalid matrix dimensions.");
		if (aN != bM) 								throw new InvalidParameterException("Matrix.multiply(): Nonmatching matrix dimensions.");
//		mulFlops_DEBUG = aM * 2 + aM * aN;			// counts number of multiplications
		
		double[] dataA = this.getDataRef()[0], dataB = B.getDataRef()[0];
		
		String newname = (Matrix.DEBUG_LEVEL > 1 ? "(" + name + "." + B.name + ")" : "P");
		Matrix C = new Matrix(newname, aM, bN);
		double[] dataC = C.data = new double[aM * bN];

		for (int i = 0; i < aM; i++) {
			int iN = i * aN, iCN = C.N * i;
			for (int j = 0; j < aN; j++) {
				int jBN = j * bN;
				double v = dataA[iN + j];
				if (!nearZero(v)) { 					// optimisation for sparse matrices
					for (int k1 = iCN, k2 = jBN, k1End = iCN + bN; k1 < k1End; k1++, k2++)
						dataC[k1] += v * dataB[k2];
//					mulFlops_DEBUG += bN;
//					mulAdops_DEBUG += bN;
				}
			}
		}		
	
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		C.clearNull();
		C.halfBandwidth = -1;
		C.bitImage = new BinBitImage(C);				// bitImage needs full reconstitution
		return C;
	}

		
	
	// polymorphic matrix product will do AB = C, with zero-ignoring optimisation for sparse matrices
	// this is attempt at a cache-optimised nested blocked-loops implementation with unrolled innermost adding to matrix C
	// thus the processed added elements of C are completely blocked in the cache, the outer iteration will scan over CONCURRENT_LOOPS number
	// of rows of B (which hopefully fit into L2) and CACHED_DOUBLES number of elements (affixed to 4) of A are fixed in static variables
	// the unrolling (affixed to 4 elements) eliminates redundant tests and the looping overhead for adding to the CACHED_DOUBLES of elements of C
	// note1: contains two code blocks (uncomment whichever), one for unrolled inner loop and one for regular inner loop
	// note2: this is in fact much slower (!) than regular multiply(), which tells to trust the compiler and that handoptimising can be a waste of time
	final static int CONCURRENT_LOOPS = 4;
	public Matrix multiply_cacheopt(Matrix B) {

		int aM = M, aN = N, bM = B.M, bN = B.N;
		if (aM < 1 || aN < 1 || bM < 1 || bN < 1) 	throw new InvalidParameterException("Matrix.multiply(): Invalid matrix dimensions.");
		if (aN != bM) 								throw new InvalidParameterException("Matrix.multiply(): Nonmatching matrix dimensions.");
		mulFlops_DEBUG = aM * 2 + aM * aN;			// counts number of multiplications
		
		double[] dataA = this.getDataRef()[0], dB = B.getDataRef()[0];
		
		String newname = (Matrix.DEBUG_LEVEL > 1 ? "(" + name + "." + B.name + ")" : "P");
		Matrix C = new Matrix(newname, aM, bN);
		double[] dC = C.data = new double[aM * bN];

		int andBitsBlk = 0xFFFFFFFF - (CACHED_DOUBLES - 1), kcEnd = bN & andBitsBlk, bNRem = bN & (CONCURRENT_LOOPS - 1);
		int jBlock = CONCURRENT_LOOPS * bN;
		
		for (int i = 0; i < bN; i++) {												// iterate in cacheblocks over rows of A
			int iN = i * aN, iCN = bN * i;
			for (int jc = 0; jc < aN; jc += CONCURRENT_LOOPS) {						// iterate in concurrent summations over columns of B
				
				int iNjc2 = iN + jc, k22 = jc * bN;
				int k2Blk, aMRem;
				if (jc + CONCURRENT_LOOPS > bM) { aMRem = bM - jc; k2Blk = aMRem * bN; }
				else { k2Blk = jBlock; aMRem = CONCURRENT_LOOPS; }
				
				int nnzA = 0;
				double v1 = 0, v2 = 0, v3 = 0, v4 = 0;
				switch (aMRem) {													// prepare quick test for zero status of values from A
					case 1:	if(!nearZero(v1=dataA[iNjc2])) nnzA++; break;
					case 2: if(!nearZero(v1=dataA[iNjc2])) nnzA++; if(!nearZero(v2=dataA[iNjc2+1])) nnzA+=2;
					case 3: if(!nearZero(v1=dataA[iNjc2++])) nnzA++; if(!nearZero(v2=dataA[iNjc2++])) nnzA+=2;
							if(!nearZero(v3=dataA[iNjc2])) nnzA+=4; iNjc2-=2; break;
					case 4: if(!nearZero(v1=dataA[iNjc2++])) nnzA++; if(!nearZero(v2=dataA[iNjc2++])) nnzA+=2;
							if(!nearZero(v3=dataA[iNjc2++])) nnzA+=4; if(!nearZero(v4=dataA[iNjc2])) nnzA+=8; iNjc2-=3;
				}
				if (nnzA == 0) continue;
				
				for (int kc = 0; kc < bN; kc += CACHED_DOUBLES) {					// iterate in cacheblocks over columns of C
					
					// k1Blk will switch from iterating over full cacheline to only remaining elements at end of row
					int k1Blk = kc < kcEnd ? CACHED_DOUBLES : bNRem, k2Add = bN - k1Blk;
					int k12 = iCN + kc, k2 = k22 + kc, k2End = k2 + k2Blk, k1;	
					// DEBUG: the inner loop unrolled option
					switch (nnzA) {
					case 0: 	k2 += bN * 4; break;
					case 1:	k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v1*dB[k2++];case 3:dC[k1++]+=v1*dB[k2++];case 2:dC[k1++]+=v1*dB[k2++];case 1:dC[k1++]+=v1*dB[k2++];}break;
					case 2:		k2 += bN;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v2*dB[k2++];case 3:dC[k1++]+=v2*dB[k2++];case 2:dC[k1++]+=v2*dB[k2++];case 1:dC[k1++]+=v2*dB[k2++];}break;
					case 3:	k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v1*dB[k2++];case 3:dC[k1++]+=v1*dB[k2++];case 2:dC[k1++]+=v1*dB[k2++];case 1:dC[k1++]+=v1*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v2*dB[k2++];case 3:dC[k1++]+=v2*dB[k2++];case 2:dC[k1++]+=v2*dB[k2++];case 1:dC[k1++]+=v2*dB[k2++];}break;
					case 4:		k2 += bN * 2;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v3*dB[k2++];case 3:dC[k1++]+=v3*dB[k2++];case 2:dC[k1++]+=v3*dB[k2++];case 1:dC[k1++]+=v3*dB[k2++];}break;
					case 5:	k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v1*dB[k2++];case 3:dC[k1++]+=v1*dB[k2++];case 2:dC[k1++]+=v1*dB[k2++];case 1:dC[k1++]+=v1*dB[k2++];}k2+=k2Add+bN;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v3*dB[k2++];case 3:dC[k1++]+=v3*dB[k2++];case 2:dC[k1++]+=v3*dB[k2++];case 1:dC[k1++]+=v3*dB[k2++];}break;
					case 6:		k2 += bN;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v2*dB[k2++];case 3:dC[k1++]+=v2*dB[k2++];case 2:dC[k1++]+=v2*dB[k2++];case 1:dC[k1++]+=v2*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v3*dB[k2++];case 3:dC[k1++]+=v3*dB[k2++];case 2:dC[k1++]+=v3*dB[k2++];case 1:dC[k1++]+=v3*dB[k2++];}break;
					case 7:	k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v1*dB[k2++];case 3:dC[k1++]+=v1*dB[k2++];case 2:dC[k1++]+=v1*dB[k2++];case 1:dC[k1++]+=v1*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v2*dB[k2++];case 3:dC[k1++]+=v2*dB[k2++];case 2:dC[k1++]+=v2*dB[k2++];case 1:dC[k1++]+=v2*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v3*dB[k2++];case 3:dC[k1++]+=v3*dB[k2++];case 2:dC[k1++]+=v3*dB[k2++];case 1:dC[k1++]+=v3*dB[k2++];}break;
					case 8:		k2 += bN * 3;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v4*dB[k2++];case 3:dC[k1++]+=v4*dB[k2++];case 2:dC[k1++]+=v4*dB[k2++];case 1:dC[k1++]+=v4*dB[k2++];}break;
					case 9:	k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v1*dB[k2++];case 3:dC[k1++]+=v1*dB[k2++];case 2:dC[k1++]+=v1*dB[k2++];case 1:dC[k1++]+=v1*dB[k2++];}k2+=k2Add+bN*2;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v4*dB[k2++];case 3:dC[k1++]+=v4*dB[k2++];case 2:dC[k1++]+=v4*dB[k2++];case 1:dC[k1++]+=v4*dB[k2++];}break;
					case 10:	k2 += bN;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v2*dB[k2++];case 3:dC[k1++]+=v2*dB[k2++];case 2:dC[k1++]+=v2*dB[k2++];case 1:dC[k1++]+=v2*dB[k2++];}k2+=k2Add+bN;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v4*dB[k2++];case 3:dC[k1++]+=v4*dB[k2++];case 2:dC[k1++]+=v4*dB[k2++];case 1:dC[k1++]+=v4*dB[k2++];}break;
					case 11:k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v1*dB[k2++];case 3:dC[k1++]+=v1*dB[k2++];case 2:dC[k1++]+=v1*dB[k2++];case 1:dC[k1++]+=v1*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v2*dB[k2++];case 3:dC[k1++]+=v2*dB[k2++];case 2:dC[k1++]+=v2*dB[k2++];case 1:dC[k1++]+=v2*dB[k2++];}k2+=k2Add+bN;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v4*dB[k2++];case 3:dC[k1++]+=v4*dB[k2++];case 2:dC[k1++]+=v4*dB[k2++];case 1:dC[k1++]+=v4*dB[k2++];}break;
					case 12:	k2 += bN * 2;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v3*dB[k2++];case 3:dC[k1++]+=v3*dB[k2++];case 2:dC[k1++]+=v3*dB[k2++];case 1:dC[k1++]+=v3*dB[k2++];}k2+=k2Add+bN;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v4*dB[k2++];case 3:dC[k1++]+=v4*dB[k2++];case 2:dC[k1++]+=v4*dB[k2++];case 1:dC[k1++]+=v4*dB[k2++];}break;
					case 13:k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v1*dB[k2++];case 3:dC[k1++]+=v1*dB[k2++];case 2:dC[k1++]+=v1*dB[k2++];case 1:dC[k1++]+=v1*dB[k2++];}k2+=k2Add+bN;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v3*dB[k2++];case 3:dC[k1++]+=v3*dB[k2++];case 2:dC[k1++]+=v3*dB[k2++];case 1:dC[k1++]+=v3*dB[k2++];}k2+=k2Add+bN;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v4*dB[k2++];case 3:dC[k1++]+=v4*dB[k2++];case 2:dC[k1++]+=v4*dB[k2++];case 1:dC[k1++]+=v4*dB[k2++];}break;		
					case 14:	k2 += bN;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v2*dB[k2++];case 3:dC[k1++]+=v2*dB[k2++];case 2:dC[k1++]+=v2*dB[k2++];case 1:dC[k1++]+=v2*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v3*dB[k2++];case 3:dC[k1++]+=v3*dB[k2++];case 2:dC[k1++]+=v3*dB[k2++];case 1:dC[k1++]+=v3*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v4*dB[k2++];case 3:dC[k1++]+=v4*dB[k2++];case 2:dC[k1++]+=v4*dB[k2++];case 1:dC[k1++]+=v4*dB[k2++];}break;		
					case 15:k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v1*dB[k2++];case 3:dC[k1++]+=v1*dB[k2++];case 2:dC[k1++]+=v1*dB[k2++];case 1:dC[k1++]+=v1*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v2*dB[k2++];case 3:dC[k1++]+=v2*dB[k2++];case 2:dC[k1++]+=v2*dB[k2++];case 1:dC[k1++]+=v2*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v3*dB[k2++];case 3:dC[k1++]+=v3*dB[k2++];case 2:dC[k1++]+=v3*dB[k2++];case 1:dC[k1++]+=v3*dB[k2++];}k2+=k2Add;if(k2>k2End)break;
							k1=k12;switch (k1Blk) {case 4:dC[k1++]+=v4*dB[k2++];case 3:dC[k1++]+=v4*dB[k2++];case 2:dC[k1++]+=v4*dB[k2++];case 1:dC[k1++]+=v4*dB[k2++];}break;		
					}
					// DEBUG: the inner unrolled loop option

					// DEBUG: the inner regular loop option
					/*for (int nnzA2 = nnzA, iNjc = iNjc2; k2 < k2End;) {
						double v = dataA[iNjc++];
						if ((nnzA2 & 1) != 0) { 									// optimisation for sparse matrices: check for zero mult.
							k1 = iCN + kc;
							switch (k1Blk) {	case 4: dC[k1++] += v * dB[k2++]; case 3: dC[k1++] += v * dB[k2++];
												case 2: dC[k1++] += v * dB[k2++]; case 1: dC[k1++] += v * dB[k2++]; }
							k2 += k2Add;
						} else k2 += bN;
						nnzA2 >>= 1;
					}*/
					// DEBUG: the inner regular loop option
				}
			}
		}		

		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		C.clearNull();
		C.halfBandwidth = -1;
		C.bitImage = new BinBitImage(C);				// bitImage needs full reconstitution
		return C;
	}
	
	
	
	// multiplies/scales matrix with a value (bitimage will be unchanged)
	public Matrix multiply(double v, boolean copy) {
		
		if (M < 1 || N < 1) throw new InvalidParameterException("Matrix.multiply(): Invalid matrix dimensions.");
	
		Matrix C = this;
		double[] dataA = this.getDataRef()[0], dataC = dataA;
		if (copy) {
			String newname = (Matrix.DEBUG_LEVEL > 1 ? "s*" + name : "S" + nameCount++);
			C = new Matrix(newname, M, N, Matrix.Type.Null, 1);
			dataC = C.data;
		}		
		
		for (int i = 0, MN = M * N; i < MN; i++) dataC[i] = v * dataA[i];
		
		C.data = dataC;
		if (M == N) C.det = Math.pow(v, M) * det;
		C.clearNull();
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		return C;
	}
	
	
	

	static double[][] buffer2x2StrasWin;				// preallocated buffer for the 2x2 submatrix cases
	
	// multiplier uses the Strassen-Winograd algorithm for matrix multiplication, adapted for flat data arrays
	// method accepts only power of 2 matrices, buffering of matrices with rescale(oldM, oldN, newM, newN, true) is recommended
	// truncP = the size of submatrix at which recursion stops and the normal multiplicator takes over
	// the algorithm uses the least amount of buffers possible, by consistently reusing them instead of allocating new ones
	public static Matrix multiplyStrasWin(Matrix A, Matrix B, int truncP) {
		
		if (A.M < 1 || A.N < 1 || B.M < 1 || B.N < 1) throw new InvalidParameterException("Matrix.multiply(): Invalid matrix dimensions.");
		if (A.N != B.M) throw new InvalidParameterException("Matrix.multiply(): Nonmatching matrix dimensions.");
		mulSW_DEBUG = mulFlopsSW_DEBUG = mulAdopsSW_DEBUG = 0;
		
		// squarify and expand matrices to nearest 2^n size
		int majorM = A.M > A.N ? A.M : A.N, newM;
		for (newM = 2; newM < majorM; newM <<= 1);		// do 2^(n+1) until we reach or excel major size of the matrices
		if (A.M != A.N || newM != majorM) {
			A = A.rescale(0, 0, newM, newM, false);		// no BitImage needed for transient data
			B = B.rescale(0, 0, newM, newM, false);
		}
		if (truncP > newM) truncP = newM;				// truncation point can't be larger than the matrix

		String newname = (DEBUG_LEVEL > 1 ? "(" + A.name + "." + B.name + ")" : "W" + Matrix.nameCount++);
		Matrix C = new Matrix(newname, newM, newM, Type.Null, 1);
		
		buffer2x2StrasWin = new double[8][truncP > 16 ? 16*16 : truncP*truncP];
		
		// subInfo carries along following data:
		// (truncP) the matrix size truncation point when standard multiplicator is used instead
		// (dim) the current submatrix dimension 
		// 3x offsets into the flat datafields, 3x the matrix width of the datafields
		int[] subInfo = {truncP, newM, 0, 0, 0, newM, newM, newM};
		multiplyStrasWin2(A.data, B.data, C.data, subInfo);
		
		C.clearNull();
		C.halfBandwidth = -1;
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		return C;

	}
	
	public static void multiplyStrasWin2(double[] dA, double[] dB, double[] dC, int[] subInfo) {

		int truncP = subInfo[0], dim = subInfo[1];
		int offsA = subInfo[2], offsB = subInfo[3], offsC = subInfo[4];
		int dimA = subInfo[5], dimB = subInfo[6], dimC = subInfo[7];
//		Matrix.mulSW_DEBUG++;
//		Matrix.mulFlopsSW_DEBUG += 4;
				
		// if we reached base case of 2x2, return straight 2x2 multiplication
		if (dim == 2) {
			double a11 = dA[offsA], a12 = dA[offsA + 1], a21 = dA[offsA + dimA], a22 = dA[offsA + dimA + 1];
			double b11 = dB[offsB], b12 = dB[offsB + 1], b21 = dB[offsB + dimB], b22 = dB[offsB + dimB + 1];
			dC[offsC] = a11 * b11 + a12 * b21;
			dC[offsC++ + dimC] = a21 * b11 + a22 * b21;
			dC[offsC] = a11 * b12 + a12 * b22;
			dC[offsC + dimC] = a21 * b12 + a22 * b22;
//			Matrix.mulFlopsSW_DEBUG += 8;
//			Matrix.mulAdopsSW_DEBUG += 11;
			return;
		}
		
		// if we reached the recursion truncation point, multiply current submatrix with "regular" algorithm
		if (truncP >= dim) {
//			Matrix.mulFlopsSW_DEBUG += dim * 2 + dim * dim + dim * dim * dim;
//			Matrix.mulAdopsSW_DEBUG += dim * 2 + dim * dim + dim * dim * dim * 2;
			for (int i = 0; i < dim; i++) {
				int ioffsA = offsA + i * dimA, ioffsC = offsC + i * dimC;
				for (int j = 0; j < dim; j++) {
					int joffsB = offsB + j * dimB;
					double v = dA[ioffsA + j];
					if (!nearZero(v))
						for (int k = 0; k < dim; k++)
							dC[ioffsC + k] += v * dB[joffsB + k];
				}
			}		
			return;
		}
		
		int sdim2 = dim>>1, ssize = sdim2 * sdim2;		// next sublevel dimensions
		// offsets to middle row (half-way into submatrix) of current recursive level
		int hoffsA = dimA * sdim2, hoffsB = dimB * sdim2, hoffsC = dimC * sdim2;
		// the skips tell how many steps to skip to get to next row in the data during incremental operations
		int skipA = dimA - sdim2, skipB = dimB - sdim2, skipC = dimC - sdim2;

		// A21 + A22 -> S1
		double[] dS1 = (sdim2 == 2 ? buffer2x2StrasWin[0] : new double[ssize]);
		for (int i = 0, offsA21 = offsA + hoffsA, offsA22 = offsA21 + sdim2, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	 dS1[offs] = dA[offsA21++] + dA[offsA22++];
			offsA21 += skipA; offsA22 += skipA;
		}

		// B12 - B11 -> T1
		double[] dT1 = (sdim2 == 2 ? buffer2x2StrasWin[1] : new double[ssize]);
		for (int i = 0, offsB11 = offsB, offsB12 = offsB + sdim2, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	 dT1[offs] = dB[offsB12++] - dB[offsB11++];
			offsB12 += skipB; offsB11 += skipB;
		}

		// S1 . T1 -> P5
		double[] dP5 = (sdim2 == 2 ? buffer2x2StrasWin[2] : new double[ssize]);
		int[] subInfo3 = {truncP, sdim2, 0, 0, 0, sdim2, sdim2, sdim2};
		multiplyStrasWin2(dS1, dT1, dP5, subInfo3);

		// B22 - T1 -> T2 (=T1)
		double[] dT2 = dT1;
		for (int i = 0, offsB22 = offsB + hoffsB + sdim2, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	 dT2[offs] = dB[offsB22++] - dT1[offs];
			offsB22 += skipB;
		}

		// S1 - A11 -> S2 (=S1)
		double[] dS2 = dS1;
		for (int i = 0, offsA11 = offsA, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	 dS2[offs] = dS1[offs] - dA[offsA11++];
			offsA11 += skipA;
		}
			
		// S2 . T2 -> P6
		double[] dP6 = (sdim2 == 2 ? buffer2x2StrasWin[3] : new double[ssize]);
		multiplyStrasWin2(dS2, dT2, dP6, subInfo3);

		// A12 - S2 -> S4
		double[] dS4 = (sdim2 == 2 ? buffer2x2StrasWin[4] : new double[ssize]);
		for (int i = 0, offsA12 = offsA + sdim2, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	 dS4[offs] = dA[offsA12++] - dS2[offs];
			offsA12 += skipA;
		}

		// S4 . B22 -> P3
		double[] dP3 = (sdim2 == 2 ? buffer2x2StrasWin[5] : new double[ssize]);
		int[] subInfo5 = {truncP, sdim2, 0, offsB + hoffsB + sdim2, 0, sdim2, dimB, sdim2};
		multiplyStrasWin2(dS4, dB, dP3, subInfo5);

		// A11 - A21 -> S3 (=S1=S2)
		double[] dS3 = dS2;
		for (int i = 0, offsA11 = offsA, offsA21 = offsA + hoffsA, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	 dS3[offs] = dA[offsA11++] - dA[offsA21++];
			offsA11 += skipA; offsA21 += skipA;
		}
		
		// B22 - B12 -> T3 (=S4)
		double[] dT3 = dS4;
		for (int i = 0, offsB12 = offsB + sdim2, offsB22 = offsB12 + hoffsB, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	 dT3[offs] = dB[offsB22++] - dB[offsB12++];
			offsB22 += skipB; offsB12 += skipB;
		}

		// T2 - B21 -> T4 (=T1=T2)
		double[] dT4 = dT2;
		for (int i = 0, offsB21 = offsB + hoffsB, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	 dT4[offs] = dT2[offs] - dB[offsB21++];
			offsB21 += skipB;
		}

		// A11 . B11 -> P1
		int[] subInfo2 = {truncP, sdim2, offsA, offsB, 0, dimA, dimB, sdim2};
		double[] dP1 = (sdim2 == 2 ? buffer2x2StrasWin[6] : new double[ssize]);
		multiplyStrasWin2(dA, dB, dP1, subInfo2);

		// S3 . T3 -> P7
		double[] dP7 = (sdim2 == 2 ? buffer2x2StrasWin[7] : new double[ssize]);
		multiplyStrasWin2(dS3, dT3, dP7, subInfo3);

		// P1 + P6 -> U2 (=T3=S4)
		double[] dU2 = dT3;
		for (int i = 0, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	dU2[offs] = dP1[offs] + dP6[offs];
		}

		// A12 . B21 -> P2 (=P6)
		subInfo2[2] += sdim2;
		subInfo2[3] += hoffsB;
		double[] dP2 = dP6;
		multiplyStrasWin2(dA, dB, dP2, subInfo2);

		// P1 + P2 -> C11
		for (int i = 0, offsC11 = offsC, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	dC[offsC11++] = dP1[offs] + dP2[offs];
			offsC11 += skipC;
		}

		// A22 . T4 -> P4 (=P2)
		double[] dP4 = dP2;
		subInfo2[2] = offsA + hoffsA + sdim2;
		subInfo2[3] = 0;
		subInfo2[5] = dimA;
		subInfo2[6] = sdim2;
		multiplyStrasWin2(dA, dT4, dP4, subInfo2);
							
		// U2 + P7 -> U3 (=T1=T2=T4)
		double[] dU3 = dT4;
		for (int i = 0, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	dU3[offs] = dU2[offs] + dP7[offs];
		}

		// U2 + P5 -> U4 (=S1=S2=S3)
		double[] dU4 = dS3;
		for (int i = 0, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	dU4[offs] = dU2[offs] + dP5[offs];
		}

		// U4 + P3 -> U5 (=T3=S4=U2), C12
		double[] dU5 = dU2;
		for (int i = 0, offsC12 = offsC + sdim2, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	 dU5[offs] = dC[offsC12++] = dU4[offs] + dP3[offs];
			offsC12 += skipC; 
		}

		// U3 - P4 -> C21
		for (int i = 0, offsC21 = offsC + hoffsC, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	dC[offsC21++] = dU3[offs] - dP4[offs];
			offsC21 += skipC;
		}

		// U3 + P5 -> C22
		for (int i = 0, offsC22 = offsC + hoffsC + sdim2, offs = 0; i < sdim2; i++)  {
			for (int j = 0; j < sdim2; j++, offs++)	dC[offsC22++] = dU3[offs] + dP5[offs];
			offsC22 += skipC; 
		}
		
//		Matrix.mulAdopsSW_DEBUG += (23 + sdim2 * 31);
	}
	

	
	
	// polymorphic comparer
	public boolean equals(Matrix B) {
		
		if (B.M != M || B.N != N) return false;
		
		// compare bitImages directly for bitwise inequalities
		if (!bitImage.equals(B.bitImage)) return false;
		
		double[][] dataSet = getDataRef();
		double[] dataA = dataSet[0], idataA = dataSet[1];
		dataSet = B.getDataRef();
		double[] dataB = dataSet[0], idataB = dataSet[1];

		for (int i = 0; i < bitImage.data.length; i++)
			if (!BinBitImage.compare(dataA, dataB, bitImage.data[i])) return false;

		// ordinary cell-by-cell comparison
		for (int i = 0; i < B.M; i++) {
			int iN = i * N;
			for (int j = iN, jEnd = j + N; j < jEnd; j++)
				if (dataA[j] != dataB[j]) return false;
		}

		// binary subdivision comparison with aid from the BinBitImage bitflags
		if (idataA != null) {
			if (idataB != null) {
				for (int i = 0; i < bitImage.data.length; i++)
					if (!BinBitImage.compare(idataA, idataB, bitImage.data[i])) return false;
			// one matrix having imaginary data while the other not constitutes inequality
			} else return false;
		}
		// both matrices lacking imaginary data (while passing real-value test) constitutes equality
		if (idataB == null) return true;
		// one matrix having imaginary data while the other not constitutes inequality
		return false;		
	}

	
	
	// method tests how far the non-zero elements extend from the diagonal, returning a value of the range
	// this allows elimination-type solvers to stop elimination beyond the width of the diagonal data band
	// method runs from left/right edges towards the diagonal in a symmetrical fashion, looking for first non-zero element
	public int getHalfBandwidth() {

		if (M < 1 || N < 1) throw new InvalidParameterException("Matrix.getHalfBandwidth(): Invalid matrix dimansions.");
		if (M != N) throw new InvalidParameterException("Matrix.getHalfBandwidth(): Matrix not square.");

		halfBandwidth = 0;
		int symOffs = M * M - 1;
		for (int r = 1; r < M; r++) {
			int rN = r * N, rNsym = symOffs - rN;
			for (int i = 0, rNi = rN, rNsymi = rNsym; i < r - halfBandwidth; i++, rNi++, rNsymi--) {
				// is element non-zero?
				if (!nearZero(data[rNi]) || !nearZero(data[rNsymi])) {
					// found the first/farthest-out element of this row, break out of loop
					if (r - i > halfBandwidth) halfBandwidth = r - i;
					break;
				}
				if (idata != null) {
					// is imaginary element non-zero?
					if (!nearZero(idata[rNi]) || !nearZero(idata[rNsymi])) {
						// found the first/farthest-out element of this row, break out of loop
						if (r - i > halfBandwidth) halfBandwidth = r - i;
						break;
					}				
				}
			}
		}
		if (DEBUG_LEVEL > 1) System.out.println("half bandwidth: " + halfBandwidth);
		return halfBandwidth;
	}
	
	

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			LINEAR SYSTEM SOLVING METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	// applies iterative convergence criterion on matrix and gives true if it's guaranteed to converge
	// criterion is that the sum of quotient of row's elements and diagonal element should be below 1 for every row
	// iterative convergency is strongly influenced by diagonal elements being as large as possible
	public boolean isConvergent() {

		boolean converges = true;
		if (DEBUG_LEVEL > 1) System.out.println("matrix " + this.name + " convergence criterion coefficients:");
		for (int r = 0; r < M; r++) {
			int rN = r * N;
			double critSum = 0;
			for (int c = 0; c < N; c++) {
				if (c != r) {
					double v = data[rN + c] / data[rN + r];
					critSum += (v > 0 ? v : -v);
					if (critSum >= 1) {
						converges = false;
						if (DEBUG_LEVEL == 0) return false;
					}
				}
			}
			if (DEBUG_LEVEL > 1) System.out.print(String.format("%5.2f", critSum) + (r < M - 1 ? ", " : " "));
		}
		if (DEBUG_LEVEL > 1) System.out.println("\nConvergence" + (converges ? "" : " not") + " guaranteed\n");
		return converges;
	}
	
	
	
	// does Gauss elimination, returns false if matrix was singular
	// taken from http://introcs.cs.princeton.edu/java/95linear/Matrix.java.html
	public boolean doGaussEliminationPP() {

		if (M != N) throw new InvalidParameterException("Matrix.doGaussEliminationPP(): Matrix not square.");

		// get a copy of matrix data field
		double[] data = getDataRef()[0];

		// Gaussian elimination with partial pivoting
		for (int i = 0; i < N; i++) {
			int iN = i * N;

			// find pivot row and swap
			int max = i;
			double vabs = Math.abs(data[max * N + i]);
			for (int j = i + 1; j < N; j++) {
				if (Math.abs(data[j * N + i]) > vabs)	max = j;
			}
			// swap rows i & max
			swapRows(i, max);

			// singular
			if (nearZero(data[iN + i])) { status |= SINGULAR_MATRIX; return false; }

			// pivot within A
			double iival = data[iN + i];
			for (int j = i + 1; j < N; j++) {
				double p = data[j * N + i] / iival;
				int jN = j * N;
				for (int k = i + 1; k < N; k++) {
					data[jN + k] -= data[iN + k] * p;
				}
				data[jN + i] = 0.0;
			}
		}
		// the matrix will have a changed data field & bitImage
		putDataRef(data, null);
		if (DEBUG_LEVEL > 1) {
			System.out.println("Gaussian elimination:");
			System.out.println(this.toString());
		}
		
		bitImage.make();
		return true;
	}

	
	// return x = A^-1 b, assuming A is square and has full rank
	// taken from http://introcs.cs.princeton.edu/java/95linear/Matrix.java.html
	public Matrix solveGaussPP(Matrix b1) {
		if (M != N || b1.M != N || b1.N != 1)
			throw new InvalidParameterException("Matrix.solveGaussPP(): Invalid matrix/vector dimensions.");

		// create copies of the data, bring it into native matrix format
		Matrix A = this.clone();
		Matrix b = new Matrix("b", N, 1, b1.getDataRef()[0], null);
		double[] dataA = A.data, datab = b.data;

		// Gaussian elimination with partial pivoting
		for (int i = 0; i < N; i++) {
			int iN = i * N;

			// find pivot row and swap
			int max = i;
			double vAbs = dataA[iN + i];
			for (int i1 = i + 1, jNi = iN + N + i, jNiEnd = N * M; jNi < jNiEnd; jNi += N, i1++) {
				double v = dataA[jNi];
				if (vAbs < v || -vAbs > v) { max = i1; vAbs = v < 0 ? -v : v; }
			}
			if (i < max) {
				A.swapRows(i, max);
				b.swapRows(i, max);
			}

			// if matrix is singular, return null
			double dAi = dataA[iN + i];
			if (nearZero(dAi)) { status |= SINGULAR_MATRIX; return null; }

			// pivot within b
			dAi = 1.0 / dAi;
			double db = datab[i];
			for (int j = i + 1, jNi = N * j + i; j < N; j++, jNi += N)
				datab[j] -= db * dataA[jNi] * dAi;


			// pivot within A
			for (int j = i + 1, jNi = N * j + i; j < N; j++, jNi += N) {
				double p = dataA[jNi] * dAi;
				int jN = j * N;
				for (int k = i + 1, jNk = jN + k, iNk = iN + k; k < N; k++)
					dataA[jNk++] -= dataA[iNk++] * p;
				dataA[jN + i] = 0.0;
			}
		}

		String newname;
		if (Matrix.DEBUG_LEVEL > 1) newname = "x(" + this.name + "^-1*" + b1.name + ")";
		else						newname = "x" + nameCount++;
		Matrix x = new Matrix(newname, N, 1, Matrix.Type.Null, 1);
		
		// back substitution
		for (int j = N - 1; j >= 0; j--) {
			double t = 0.0;
			int jM = j * M;
			for (int k = j + 1, jMk = jM + k; k < N; k++) t += dataA[jMk++] * x.data[k];
			x.data[j] = (datab[j] - t) / dataA[jM + j];
		}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("Gaussian with partial pivoting solver:");
			System.out.println(x.toString());
		}
		x.bitImage = new BinBitImage(x);
		return x;
	}

	
	// return x = A^-1 b applying Gauss-Jordan method with full pivoting
	// adapted from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
	// added determinant computation into this algorithm
	// b[] supplies a list of constant vectors to produce results for
	// if copy = true, result matrix copy returned in Matrix array, supplied matrix is preserved
	// otherwise the supplied constant vectors transform into result vectors and supplied matrix is destroyed
	public Matrix solveGaussJordanFP(Matrix b, boolean copy) {
		
		if (M != N || b.M != N)	throw new InvalidParameterException("Matrix.solveGaussJordanFP(): Invalid matrix/vector dimensions.");

		Matrix A = this, x = b;
		if (copy) {
			x = b.clone();
			A = this.clone();
		} else
			A.halfBandwidth = x.halfBandwidth = -1;
		
		double[] datax = x.data, dataA = A.data;
		// integer arrays indxc, indxr and ipiv do the bookkeeping on the pivoting
		int[] indxc = new int[N], indxr = new int[N], iPiv = new int[N];
		int iCol = 0, iRow = 0, XN = x.N;
		double det = 1;
		
		// major loop of column reduction
		for(int i = 0; i < N; i++) {
			double vMax = 0;
			for (int j = 0; j < N; j++) {
				int jN = N * j;
				if (iPiv[j] != 1)
					for(int k = 0; k < N; k++) {
						double v = dataA[jN++];
						if (v < 0) v = -v;
						if (iPiv[k] == 0 && v >= vMax) {
							vMax = v;
							iRow = j;
							iCol = k;
						}
					}
			}
			iPiv[iCol]++;
			// if a better pivot element is found, then we swap in it's row to the diagonal,
			// columns are not swapped but relabelled: indxc[i], the column of ith pivot element
			// is the one reduced, indxr[i] is the original row location of that pivot element,
			// indxr[i] =/= indxc[i] implies column swap, with this bookkeeping the solution vectors
			// come out correct while the inverse matrix is permutated
			if (iRow != iCol) {
				swapRows(iRow, iCol);
				x.swapRows(iRow, iCol);
			}
			indxc[i] = iRow;
			indxr[i] = iCol;
			
			double pivInv = dataA[N * iCol + iCol];
			// signal a singular matrix by returning null
			if (nearZero(pivInv)) { status |= SINGULAR_MATRIX; return null; }
			// the determinat is built up by multiplying in all the pivot divisors
			det *= pivInv;
			pivInv = 1.0 / pivInv;
			int iColN = iCol * N, iColXN = iCol * x.N;
			// divide pivot row by pivot element
			for (int l = iColN, lEnd = l + N; l < lEnd; l++)		dataA[l] *= pivInv;
			for (int l = iColXN, lEnd = l + XN; l < lEnd; l++)		datax[l] *= pivInv;

			if (DEBUG_LEVEL > 2) {
				System.out.println("iCol: " + iCol + " pivot: " + 1.0/pivInv);
				System.out.println(A.toString());
			}

			for (int r = 0; r < N; r++)
				// subtract pivot row from all rows except itself
				if (r != iCol) {
					int rN = r * N, rXN = r * XN;
					double v = dataA[rN + iCol];
					if (DEBUG_LEVEL > 2) System.out.println("v:" + v);
					dataA[rN + iCol] = 0;				// the column of the pivot gets zeroed out (except on pivot row)
					for (int c1 = rN, c2 = iColN, cEnd = c1 + N; c1 < cEnd;)		dataA[c1++] -= dataA[c2++] * v;
					for (int c1 = rXN, c2 = iColXN, cEnd = c1 + XN; c1 < cEnd;)		datax[c1++] -= datax[c2++] * v;
				}
			if (DEBUG_LEVEL > 2) System.out.println(A.toString());
		}
		
		// reverse-permutate the result vectors from the column swaps, by swapping columns in reverse order of their buildup
		// (note: the original code seems to contain a mistake here)
		for (int l = N - 1; l >= 0; l--) {
			if (indxr[l] != indxc[l]) {
				for (int k = 0; k < XN; k++) {
					int kNrl = k * N + indxr[l], kNrc = k * N + indxc[l];
					double v = datax[kNrl];
					datax[kNrl] = datax[kNrc];
					datax[kNrc] = v;
				}
				det = -det;
			}
		}

		if (copy)	this.det = det;				// the determinant belongs to the untransformed matrix, which is preserved if duplicate=true
		else			x.det = det;			// if input matrix is transformed, return determinant in solution vectors matrix
		x.name = "X" + nameCount++;
		if (DEBUG_LEVEL > 1) {
			System.out.println("Gauss-Jordan full pivoting solver, result vectors:");
			System.out.println(x.toString());
		}
		return x;
	}
	
	
	// return x = A^-1 b with Gauss-Jordan partial pivoting method, inverse of A returned in Ai, returns solution vector x
	// this method checks if matrix is sparse to a diagonal band and changes row processing width to the
	// maximum overlap of the diagonal bands
	public Matrix solveGaussJordanPPDO(Matrix c, Matrix Ai) {

		if (M != N || c.M != N || c.N != 1)
							throw new InvalidParameterException("Matrix.solveGaussJordanPPDO(): Invalid matrix/vector dimensions.");
		if (Ai == null) 	throw new InvalidParameterException("Matrix.solveGaussJordanPPDO(): Invalid inverse reference supplied.");

		Matrix A = this.clone();
		det = 1;
		
		Matrix x = c.clone();
		x.name = Matrix.DEBUG_LEVEL > 1 ? "x(" + name + "^-1*" + c.name + ")" : "x" + nameCount++;
		
		Ai.M = Ai.N = M;
		Ai.generateData(Matrix.Type.Identity, 1);
		double[] dataA = A.data, dataAi = Ai.data, datax = x.data;

		// find out how far from the diagonal the sparse data stretches at maximum
		if (halfBandwidth < 0) getHalfBandwidth();
		// turn it into bandwidth number, telling how many rows above & below to process for a given row
		int bandwidth = halfBandwidth * 2;	
		
		// loop handles every case of subtracting current unitised row from all other rows
		for (int r = 0; r < M; r++) {			

			// find highest value for current diagonal element and swap in it's row
			double vMax = ROUNDOFF_ERROR;
			int rPivot = r;
			// for a matrix sparse/zeroed outside a diagonal band we'll only process the available row data above & below pivot
			// that is overlapping the pivot datawidth, the degree of overlap was calculated from return of diagonality() method
			for (int r2 = r, bandwidth1 = bandwidth; r2 < M && bandwidth1 > 0; r2++, bandwidth1--) {
				double v = dataA[N * r2 + r];
				if (v < 0) v = -v;
				if (v > vMax) { vMax = v; rPivot = r2; }
				if (DEBUG_LEVEL > 2) System.out.println(r2 + " " + v + " " + vMax);
			}
			if (vMax <= ROUNDOFF_ERROR)
				{ status |= SINGULAR_MATRIX; return null; }			// got a singular matrix, abort
			
			if (r != rPivot) {										// another row had highest value for current diagonal element, swap it in
				A.swapRows(r, rPivot);
				x.swapRows(r, rPivot);
				Ai.swapRows(r, rPivot);
				det = -det;
				if (DEBUG_LEVEL > 2) {
					System.out.println("Swapped in row " + rPivot + " into pivot row " + r);
					System.out.println(A.toString());
				}
			}

			int iStart = (r - bandwidth < 0 ? 0 : r - bandwidth);		// find out bounds of sparse datawidth of pivot row
			int iEnd = (r + bandwidth >= M ? M - 1 : r + bandwidth);
			
			for (int i = iStart; i <= iEnd; i++) {						// the diagonal element column unitising loop
				int iN = i * N;
				double div = dataA[iN + r];								// get divisor for current lined-up row of A, c, Ai structures
				if (!nearZero(div)) {
					det *= div;
					div = 1.0 / div;
					// go along the lined-up row of A, c, Ai structures, dividing by current element
					if (i <= r) {
						for (int j = iN, jEnd = j + N; j < jEnd; j++) {
							dataA[j] *= div;								// divide row in A (which moves towards Identity matrix)
							dataAi[j] *= div;								// divide row in Ai (which moves towards an inverse of A)
						}
					} else {
						// we can ignore the already zeroed lower-left triangle of matrix
						for (int j = iN + r, jEnd = j + N - r; j < jEnd; j++)	dataA[j] *= div;
						for (int j = iN, jEnd = j + N; j < jEnd; j++)			dataAi[j] *= div;
					}
					datax[i] *= div;									// divide input vector c
				} else {
					// divisor was zero, operate on current row with the pivot row added to it
					int rN = r * N;
					det *= dataA[rN + r] + dataA[iN + r];
					div = 1.0 / dataA[rN + r];
					if (i <= r)
						for (int j = 0, iNj = iN, rNj = rN; j < N; j++, iNj++, rNj++) {
							dataA[iNj] = (dataA[iNj] + dataA[rNj]) * div;		// divide row in A (which moves towards Identity matrix)
							dataAi[iNj] = (dataAi[iNj] + dataAi[rNj]) * div;	// divide row in Ai (which moves towards an inverse of A)
						}
					else {
						// we can ignore the already zeroed lower-left triangle of matrix
						for (int j = r, iNj = iN+r, rNj = rN+r; j < N; j++, iNj++, rNj++)
							dataA[iNj] = (dataA[iNj] + dataA[rNj]) * div;
						for (int j = 0, iNj = iN, rNj = rN; j < N; j++, iNj++, rNj++)
							dataAi[iNj] = (dataAi[iNj] + dataAi[rNj]) * div;
					}
					datax[i] = (datax[i] + datax[r]) * div;				// divide input vector c
				}
			}
			if (DEBUG_LEVEL > 2) System.out.println(A.toString());

			int rN = r * N;
			// for every row i in the lined-up A, c, Ai structures, subtract row (of diagonal element) r from row i, both above and below r
			for (int i = iStart; i <= iEnd; i++) {
				int iN = i * N;				
				// don't subtract row r from itself
				if (i != r) {
					if (i < r)
						for (int j = 0, iNj = iN, rNj = rN; j < N; j++, iNj++, rNj++) {
							dataA[iNj] -= dataA[rNj];					// subtract from A
							dataAi[iNj] -= dataAi[rNj];					// subtract from Ai
						}
					else {
						for (int j = r, iNj = iN+r, rNj = rN+r; j < N; j++, iNj++, rNj++)
							dataA[iNj] -= dataA[rNj];					// subtract from A
						for (int j = 0, iNj = iN, rNj = rN; j < N; j++, iNj++, rNj++)
							dataAi[iNj] -= dataAi[rNj];					// subtract from Ai
					}
					datax[i] -= datax[r];							// subtract from c
				}
			}
			if (DEBUG_LEVEL > 2) System.out.println(A.toString());
		}
		
		// finally, divide all rows by a'(ii) for a final normalising, this finally turns A into an Identity matrix
		for (int r = 0; r < M; r++) {
			int rN = r * N;
			double div = dataA[rN + r];
			det *= div;
			div = 1.0 / div;
			
			dataA[rN + r] *= div;									// A is almost an identity matrix, just one element to divide
			for (int j = rN, jEnd = j + N; j < jEnd; j++)
				dataAi[j] *= div;									// divide row r of Ai by a(ii), turning it into the inverse of A

			datax[r] *= div;										// divide c(i) by a(ii), it will turn into solution vector x
		}
		A.det = 1.0;												// |I| = 1
		Ai.det = 1.0 / det;											// |A^-1| = 1 / |A|

		if (DEBUG_LEVEL > 1) {
			System.out.println("Gauss-Jordan partial pivoting solver with sparse inverse matrix calculation, result vector:");
			System.out.println(x.toString());
			System.out.println("input matrix transformed to identity:");
			System.out.println(A.toString());
			System.out.println("input matrix inverse:");
			System.out.println(Ai.toString());
			System.out.println("half bandwidth: " + halfBandwidth);
		}
		return x;
	}

	
	
	// iterative Gauss-Seidel solver will not guranteedly succeed to converge to a vector,
	// but could be used with certainty if convergence() criterion is fulfilled, and fall back
	// to standard Gaussian elimination solvers
	public Matrix solveGaussSeidel(Matrix c, int itersMax, double maxError) {
		
		if (M != N || c.M != N || c.N != 1)
			throw new InvalidParameterException("Matrix.solveGaussJordan(): Invalid matrix/vector dimensions.");

		DEBUG_LEVEL--;
		
		// create initial solution vector
		Matrix x = new Matrix("x", M, 1, Matrix.Type.Random, 1), xOld = null;
		double[] datac = c.data, datax = x.data;
		// current iteration error & previous iteration error
		double cuErr = 0, prErr = 0;
		boolean converges = true;
		
		// do the refining iterations
		int i = 0;
		for (; i < itersMax; i++) {
			
			xOld = x.clone();
			// do Gauss-Seidel iteration step over every row
			for (int r = 0; r < M; r++) {

				double sum = 0;
				int rN = r * N;
				for (int c1 = 0,rNc = rN; c1 < M; c1++, rNc++)
					if (c1 != r)	sum += data[rNc] * datax[c1];
	
				datax[r] = (datac[r] - sum) / data[rN + r];
			}
			// check error vector length between previous & current solution vector
			System.out.println("Error: " + cuErr + ", error diff: " + (cuErr - prErr));
			prErr = cuErr;
			cuErr = xOld.subtract(x, false).absLength(0);
			//if (i > 3 && cuErr - prErr > 0) { converges = false; break; }	// latest error larger than previous, won't converge, stop iteration
			if (cuErr <= maxError && cuErr >= -maxError) break;		// error below maxError, stop iteration, quit with result
		}
		
		DEBUG_LEVEL++;
		if (DEBUG_LEVEL > 1) {
			if (!converges || i >= itersMax)
							System.out.println("Gauss-Seidel iterator solver failed converging.");
			else			System.out.println("Gauss-Seidel iterator solver done in " + i + " iterations.");
			System.out.println(x.toString());
		}	
		if (converges) return x; else return null;
	}
	
	
	
	// return x = A^-1 b with Gauss-Jordan partial pivoting method, inverse of A returned in Ai (if supplied), returns solution vector x
	// this method is like solveGaussJordanPP() minus optimisations for diagonal-sparse matrices
	public Matrix solveGaussJordanPP(Matrix b, Matrix Ai) {

		if (M != N || b.M != N || b.N != 1)
							throw new InvalidParameterException("Matrix.solveGaussJordan(): Invalid matrix/vector dimensions.");

		Matrix A = this.clone();
		det = 1;
		
		Matrix x = b.clone();
		if (Matrix.DEBUG_LEVEL > 1) x.name = "x(" + name + "^-1*" + b.name + ")";
		else						x.name = "x" + nameCount++;
		
		double[] dataA = A.data, dataAi = null, datax = x.data;
		if (Ai != null) {
			Ai.M = Ai.N = M;
			Ai.generateData(Matrix.Type.Identity, 1);
			dataAi = Ai.data;
		}
		
		// loop handles every case of subtracting current unitised row from all other rows
		for (int r = 0; r < M; r++) {			

			// find highest value for current diagonal element and swap in it's row
			double vMax = ROUNDOFF_ERROR;
			int rPivot = r;
			//for (int r2 = r, offDiag = r + diag; r2 < offDiag; r2++) {
			for (int r2 = r; r2 < M; r2++) {
				double v = dataA[N * r2 + r];
				if (v < 0) v = -v;
				if (v > vMax) { vMax = v; rPivot = r2; }
				if (DEBUG_LEVEL > 2) System.out.println(r2 + " " + v + " " + vMax);
			}
			if (vMax <= ROUNDOFF_ERROR)
				{ status |= SINGULAR_MATRIX; return null; }			// got a singular matrix, abort
			if (r != rPivot) {										// another row had highest value for current diagonal element, swap it in
				A.swapRows(r, rPivot);
				x.swapRows(r, rPivot);
				if (Ai != null) Ai.swapRows(r, rPivot);
				det = -det;
				if (DEBUG_LEVEL > 2) {
					System.out.println("Swapped in row " + rPivot + " into pivot row " + r);
					System.out.println(A.toString());
				}
			}

			for (int i = 0; i < M; i++) {							// the diagonal element column unitising loop
				int iN = i * N;
				double div = dataA[iN + r];							// get divisor for current lined-up row of A, c, Ai structures
				if (!nearZero(div)) {
					det *= div;
					div = 1.0 / div;
					// go along the lined-up row of A, c, Ai structures, dividing by current element
					if (Ai != null)
						for (int j = iN, jEnd = j + N; j < jEnd; j++) {
							dataA[j] *= div;							// divide row in A (which moves towards Identity matrix)
							dataAi[j] *= div;							// divide row in Ai (which moves towards an inverse of A)
						}
					else for (int j = iN, jEnd = j + N; j < jEnd; j++)
						dataA[j] *= div;
					
					datax[i] *= div;								// divide input vector c
					
				} else {											// divisor was zero, operate on current row with pivot row added to it
					
					int rN = r * N;
					det *= dataA[rN + r] + dataA[iN + r];
					div = 1.0 / dataA[rN + r];
					if (Ai != null)
						for (int rNj = rN, iNj = iN, iNjEnd = iN + N; iNj < iNjEnd; iNj++, rNj++) {
							dataA[iNj] = (dataA[iNj] + dataA[rNj]) * div;		// divide row in A (which moves towards Identity matrix)
							dataAi[iNj] = (dataAi[iNj] + dataAi[rNj]) * div;	// divide row in Ai (which moves towards an inverse of A)
						}
					else for (int rNj = rN, iNj = iN, iNjEnd = iN + N; iNj < iNjEnd; iNj++, rNj++)
						dataA[iNj] = (dataA[iNj] + dataA[rNj]) * div;
					
					datax[i] = (datax[i] + datax[r]) * div;			// divide input vector c
				}
			}
			if (DEBUG_LEVEL > 2) System.out.println(A.toString());

			int rN = r * N;
			// for every row i in the lined-up A, c, Ai structures, subtract row (of diagonal element) r from row i, both above and below r
			for (int i = 0; i < M; i++) {
				int iN = i * N;				
				// don't subtract row r from itself
				if (i != r) {
					if (Ai != null)
						for (int rNj = rN, iNj = iN, iNjEnd = iN + N; iNj < iNjEnd; iNj++, rNj++) {
							dataA[iNj] -= dataA[rNj];					// subtract from A
							dataAi[iNj] -= dataAi[rNj];					// subtract from Ai
						}
					else for (int rNj = rN, iNj = iN, iNjEnd = iN + N; iNj < iNjEnd; iNj++, rNj++)
						dataA[iNj] -= dataA[rNj];
					
					datax[i] -= datax[r];							// subtract from c
				}
			}
			if (DEBUG_LEVEL > 2) System.out.println(A.toString());
		}
		
		// finally, divide all rows by a'(ii) for a final normalising, this finally turns A into an Identity matrix
		for (int r = 0; r < M; r++) {
			int rN = r * N;
			det *= dataA[rN + r];
			double div = 1.0 / dataA[rN + r];
			if (Ai != null)
				for (int rNj = rN, rNjEnd = rNj + N; rNj < rNjEnd; rNj++)	{
					dataA[rNj] *= div;
					dataAi[rNj] *= div;								// divide row r of Ai by a(ii), turning it into the inverse of A
				}
			else for (int rNj = rN, rNjEnd = rNj + N; rNj < rNjEnd; rNj++)
				dataA[rNj] *= div;
			
			datax[r] *= div;									// divide c(i) by a(ii), it will turn into solution vector x
		}
		A.det = 1;												// |I| = 1
		if (Ai != null) Ai.det = 1.0 / det;						// |A^-1| = 1 / |A|

		if (DEBUG_LEVEL > 1) {
			System.out.println("Gauss-Jordan partial pivoting solver with inverse matrix calculation, result vector:");
			System.out.println(x.toString());
			System.out.println("input matrix transformed to identity:");
			System.out.println(A.toString());
			if (Ai != null) {
				System.out.println("input matrix inverse:");
				System.out.println(Ai.toString());
			}
		}
		return x;
	}


	
	
	// LU combined single-matrix LU decompose with pivoting and recording of permutations, which also calculates determinant
	// adapted from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
	// the combined decompose is stored directly in current matrix, unless copy = true, then a copy is returned
	// if generateLandU = true, a normal two-matrix L & U result returned, otherwise a combined LU matrix returned
	// also the row mutation result is stored in result matrixes mutator field
	public Matrix[] decomposeLU(boolean copy, boolean generateLandU, boolean permute) {
		
		if (M != N)	throw new InvalidParameterException("Matrix.decomposeLU2(): Matrix not square.");
		if (M < 1)	throw new InvalidParameterException("Matrix.decomposeLU2(): Invalid matrix.");
		
		Matrix[] lLU = new Matrix[2];
		int[] mutatorLU = null;
		int swaps = 0, MN = M * N;
		Matrix LU = this;
		double[] dataLU = data;
		
		if (copy) {
			LU = clone();
			LU.name = "LU" + nameCount++;
			dataLU = LU.data;
			if (mutator[0] == null) mutatorLU = LU.mutator[0] = new int[M];
		} else {
			name = "LU" + nameCount++;			// this matrix will be modified, change it's name
			mutatorLU = mutator[0] = new int[M];
		}
		
		double[] vv = new double[M];
		double d = 1;							// every permutation toggles sign of d
		
		for (int i = 0; i < M; i++) {			// looping over rows getting implicit scaling factor of every row, putting it in vv[]
			int iN = i * N;
			double vMax = 0;
			for (int j = iN, jEnd = j + M; j < jEnd; j++) {
				double v = dataLU[j];
				if (v < -vMax || v > vMax) vMax = (v < 0 ? -v : v);
			}
			if (nearZero(vMax)) { status |= SINGULAR_MATRIX; return null; }	// if singular matrix found, return null
			vv[i] = 1.0 / vMax;												// save scaling factor
			mutatorLU[i] = i;
		}
		
		for (int j = 0; j < M; j++) {
			int jN = j * N;
			for (int i = 0; i < j; i++) {
				int iN = i * N;
				double sum = dataLU[iN + j];
				for (int k = iN, kNj = j, kEnd = k + i; k < kEnd; k++, kNj += N)
					sum -= dataLU[k] * dataLU[kNj];
				dataLU[iN + j] = sum;
			}
			
			int iMax = 0;
			double vMax = 0;											// vMax will store the largest pivot element found
			for (int i = j; i < M; i++) {
				int iN = i * N;
				double sum = dataLU[iN + j];
				for (int k = iN, kNj = j, kEnd = iN + j; k < kEnd; k++, kNj += N)
					sum -= dataLU[k] * dataLU[kNj];
				dataLU[iN + j] = sum;
				double v = vv[i] * (sum < 0 ? -sum : sum);
				// have we found a better pivot than the best so far?
				if (v > vMax) { vMax = v; iMax = i; }
			}
			if (j != iMax && permute) {									// test if there's a better pivot and we need to swap rows (if permute flag allows)
				swaps++;
				LU.swapRows(iMax, j);									// yes
				//System.out.println("Swapped " + iMax + " and " + j);
				d = -d;
				vv[iMax] = vv[j];										// also change scale factor
				int temp = mutatorLU[j];
				mutatorLU[j] = mutatorLU[iMax];
				mutatorLU[iMax] = temp;
			}
			
			// zero at pivot element implies a singular matrix, but some applications can get by with tweaking away the zero
			if (nearZero(dataLU[jN + j])) {
				if (copy) status |= SINGULAR_MATRIX;					// mark original (unchanged) matrix as singular
				dataLU[jN + j] = ROUNDOFF_ERROR * (norm1 < 0 ? this.norm1() : norm1);
			}
			if (j < M - 1) {
				double v = 1.0 / dataLU[jN + j];
				for (int ij = M * (j + 1) + j; ij < MN; ij += N)
					dataLU[ij] *= v;							// divide column by pivot element
			}
		}

		lLU[0] = LU;

		// use d to calculate determinant, put it in LU
		for (int i = 0, M1 = M + 1; i < MN; i += M1) d *= dataLU[i];
		if (copy) LU.det = d; else this.det = d;
		
		// was generation of individual L & U matrices requested?
		if (generateLandU) {
			// convert LU into L, moving it's upper triangle data into U
			lLU[0].name = "L";									// here we treat dataLU as if containing data of L
			lLU[1] = new Matrix("U", M, N, Matrix.Type.Null, 1);
			double[] dataU = lLU[1].data;
			for (int i = 0; i < M; i++)
				for (int j = i; j < N; j++) {					// iterate over the top triangle to copy data into U
					dataU[i * M + j] = dataLU[i * M + j];
					if (i == j) dataLU[i * M + j] = 1;			// L's diagonal has only 1:s
					else dataLU[i * M + j] = 0;					// L's top triangle is zero
				}
			lLU[0].det = lLU[1].det = d;
			lLU[1].mutator[0] = mutatorLU;
		}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("Matrix.decomposeLU() result:");
			System.out.println(LU.toString());
			if (generateLandU) System.out.println(lLU[1].toString());
			System.out.println("LU row mutator:");
			System.out.println(Arrays.toString(mutatorLU));
		}
		return lLU;
	}
	
	

	// unpermute a matrix with mutator values altered
	public void unpermute() {
		if (mutator[0] != null) {
			int[] mutatorR = mutator[0];
			for (int i = M - 1; i >= 0; i--) if (mutatorR[i] != i)	swapRows(i, mutatorR[i]);
			halfBandwidth = -1;
		}
		if (mutator[1] != null) {
			int[] mutatorC = mutator[1];
			for (int i = M - 1; i >= 0; i--) if (mutatorC[i] != i)	swapRows(i, mutatorC[i]);
			halfBandwidth = -1;
		}
	}
	
	
	// method takes the combined factorised LU matrix using it to solve a linear system
	// where only the constant vector is changing, the matrix coefficients assumed to be fixed
	// adapted from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
	// the triangular matrices L & U come here in a COMBINED form in matrix LU, split by l(ii)+u(ii) in the diagonal; l(ii) = 1
	// matrix A, LU are applied on this vector, array containing solution vector and LU matrix is the return value
	// if LU is null, then the method calls decomposeLU2() with A, if A is null then LU is used
	// method optimisation takes into account that b might contain many leading zero elements
	public Matrix[] backSubstituteLU(Matrix A, Matrix LU, boolean copy) {
		
		if (A == null && LU == null)
			throw new InvalidParameterException("Matrix.backSubstituteLU2(): Null matrices.");
		if (N != 1 || M != (A == null ? LU.M : A.M))
			throw new InvalidParameterException("Matrix.backSubstituteLU2(): Invalid constant vector.");			
		Matrix[] bLU = new Matrix[2];
		
		// LU = null tells algorithm to create a new LU decomposition, which will be returned in the matrix array
		if (LU == null) {
			if (A.M != A.N)	throw new InvalidParameterException("Matrix.decomposeLU(): Matrix A not square.");			
			if (A.M < 1)	throw new InvalidParameterException("Matrix.decomposeLU(): Invalid matrix.");
			// copy A into a LU matrix and decompose LU, if copy=false it will destroy matrix A returning it as LU
			bLU[1] = A.decomposeLU(copy, false, true)[0];
			// decompose failed on a singular matrix, return null
			if (bLU[1] == null)	return null;	
		} else
			bLU[1] = LU;
		
		// if copy=true, do not return solution vector in b vector, but in a copy
		if (copy) bLU[0] = this.clone(); else bLU[0] = this;
		if (Matrix.DEBUG_LEVEL > 1) bLU[0].name = "x((" + bLU[1].name + ")^-1*" + name + ")";
		else						bLU[0].name = "x" + nameCount++;

		double[] dataLU = bLU[1].data, datab = bLU[0].data;
		
		int ii = -1, N = bLU[1].N;							// when ii > 0 it will become an index of first nonvanishing element of b
		int[] mutatorLU = bLU[1].mutator[0];
		double sum = 0;
		
		for (int i = 0; i < N; i++) {						
			int iN = i * N, ip = mutatorLU[i];				// the permutations of decomposition stage need to be unpermuted during the passes
			sum = datab[ip];								// depermute: take a value from it's permuted position
			datab[ip] = datab[i];							// depermute: return former value at this position to it's place
			if (ii >= 0)
				for (int jb = ii, j = iN + ii, jEnd = j + i - 1; j <= jEnd; j++, jb++)
					sum -= dataLU[j] * datab[jb];
			else
				if (!nearZero(sum)) ii = i;					// nonzero element found, we'll have to do the sums in loop above from now on
			datab[i] = sum;
		}
		for (int i = N - 1; i >= 0; i--) {					// backsubstitution pass
			int iN = i * N;
			sum = datab[i];
			for (int j = i + 1; j < N; j++)
				sum -= dataLU[iN + j] * datab[j];
			datab[i] = sum / dataLU[iN + i];				// store an element of solution vector X
		}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("Backsubstitution LU solver:");
			System.out.println(bLU[0].toString());
			//if (A != null) System.out.println(bLU[1].toString());
		}
		return bLU;
	}
	
	
	// factorises/decomposes matrix into two diagonal matrices U and V, useful for the LU backsubstitution linear solving method:
	//		uuu			v..
	// U:	.uu		L:	vv.
	//		..u			vvv
	public Matrix[] decomposeLU() {
		
		if (M != N)	throw new RuntimeException("Matrix.decomposeLU(): Matrix not square.");
		if (M < 1)	throw new RuntimeException("Matrix.decomposeLU(): Invalid matrix dimansions.");

		// all diagonal elements must be nonzero
		for (int d = 0; d < M; d++)
			if (nearZero(data[d * M + d])) { status |= SINGULAR_MATRIX; return null; }
			
		Matrix[] lLU = new Matrix[2];
		lLU[0] = new Matrix("L", M, N, Matrix.Type.Null, 1);
		lLU[1] = new Matrix("U", M, N, Matrix.Type.Null, 1);
		double[] dataL = lLU[0].data, dataU = lLU[1].data;

		for (int c1 = 0; c1 < N; c1++) {
			
			// calculate u(r,c)
			if (c1 == 0)	dataU[0] = data[0];
			else
				for (int r = 0; r <= c1; r++) {
					int rN = r * N;
					double luSum = 0;
					for (int k = 0; k < r; k++)
						luSum += dataL[rN + k] * dataU[k * N + c1];
					dataU[rN + c1] = data[rN + c1] - luSum;
				}
			
			// calculate l(r,c)
			double uDiag = 1.0 / dataU[c1 * N + c1];
			for (int r = c1; r < N; r++) {
				int rN = r * N;
				double luSum = 0;
				if (c1 == 0) dataL[rN] = data[rN];
				else {
					for (int k = 0; k < c1; k++)
						luSum += dataL[rN + k] * dataU[k * N + c1];
					dataL[rN + c1] = (data[rN + c1] - luSum) * uDiag;
				}
			}
		}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("decomposeLU() result:");
			System.out.println(lLU[0].toString());
			System.out.println(lLU[1].toString());
		}
		return lLU;
	}

	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			DETERMINANT HEURISTICS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	// helper method that does quick zero determinant heuristics on the matrix/submatrix
	public boolean hasZeroDeterminant(boolean thoroughTest) {
		// Cacc accumulates sums along the columns during the horisontal passes
		double[] data = getDataRef()[0], Cacc = new double[N], Racc = new double[M];
		// Racc accumulates sum along the row
		double v;
		for (int i = 0; i < M; i++) {
			int iN = N * i;
			for (int j = 0; j < N; j++) {
				v = data[iN + j];
				Racc[i] += v;
				Cacc[j] += v;
			}
			// found zero row, determinant = 0
			if (nearZero(Racc[i])) return true;
		}
		
		// check presense of zero columns
		for (int j = 0; j < N; j++) {
			// found zero column, determinant = 0
			if (nearZero(Cacc[j])) return true;
			double caccj = Cacc[j];
			
			if (thoroughTest) {
				// test finding two columns that are equal
				for (int i = j+1; i < N; i++) {
					double diff = Cacc[i] - caccj;
					if (nearZero(diff)) {
						// sums of two columns were equal, test individual values
						for (int ct = 0; ct < M;) {
							int ctN = ct * N;
							double diff2 = data[ctN + j] - data[ctN + i];
							if (!nearZero(diff2))
								break;
							if (++ct >= M) return true;		// columns were equal, determinant = 0
						}
					}
				}
			}
		}
		if (!thoroughTest) return false;

		// test finding two rows that are equal
		for (int j = 0; j < M; j++) {
			double raccj = Racc[j];
			for (int i = j + 1; i < M; i++) {
				double diff = Racc[i] - raccj;
				if (nearZero(diff)) {
					// sums of two rows were equal, test individual values
					for (int ct = 0; ct < N;) {
						double diff2 = data[i * N + ct] - data[j * N + ct];
						if (!nearZero(diff2)) break;
						if (++ct >= N) return true;		// columns were equal, determinant = 0
					}
				}
			}
		}
		return false;
	}
	
	
	// quick zero determinant heuristics helper method for version 2 of determinant Laplace Expansion
	// uses the relative increment list activec to jump past recursively eliminated columns
	public boolean hasZeroDeterminant2(int[][] activec, int columns) {

		int topb = M - columns;
		double[] data = getDataRef()[0];
		// Cacc accumulates sums along the columns during the horisontal passes, Racc accumulates sums along the rows
		double[] Cacc = new double[N], Racc = new double[M];
		double v;

		for (int i = topb; i < M; i++) {
			// get row index from detAnalysis if has been created
			int iN = N * (analysis == null ? i : analysis[2][i]);
			for (int j = 0, ac = activec[0][0]; j < columns; j++, ac += activec[0][ac]) {
				v = data[iN + ac - 1];
				Racc[i] += v;
				Cacc[j] += v;
			}
			// found zero row, determinant = 0
			if (nearZero(Racc[i])) return true;
		}
		
		// check presense of zero columns
		for (int j = 0; j < columns; j++) if (nearZero(Cacc[j])) return true; // found zero column, determinant = 0
		return false;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			IMPLEMENTATION OF LAPLACIAN DETERMINANT FINDER (of no particular usefulness!)
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	// multimode implementation of the Laplace Expansion Recursive algorithm for finding the determinant
	public static double determinantLaplace(Matrix M, int version) {
		
		if (M.M != M.N) throw new InvalidParameterException("Matrix.determinantLaplace(): Nonsquare matrix.");
		if (M.M < 1 || M.N < 1) throw new InvalidParameterException("Matrix.determinantLaplace(): Invalid matrix.");
		if (M.M == 1 && M.N == 1) return M.getDataRef()[0][0];
		Matrix.detL_DEBUG = 0;
		M.det = 0;
		
		if (M.hasZeroDeterminant(true)) return 0;
		
		// prepare activec: 
		// 0th array: the pseudo-linked list for fast column skipping
		// 1st array: storing the current n-tuple of column indexes for searches in the tree structure
		// 2nd-Nth arrays: different branch levels of tree structure for flagging visited branches
		int[][] activec = new int[M.M][];
		activec[0] = new int[M.N + 1];
		activec[1] = new int[M.N];
		for (int i = 0; i < activec[0].length; i++) activec[0][i] = 1;	// initial column(n)->column(n+1) step is 1
		
		switch (version) {
		case 2:
			M.analyse(false);
			// find out if the zero-count of the row or column aspect is best for determinant calculation, transpose if column aspect is best
			// note: the determinant is not affected by transposing
			if (!M.rowAspectSparsest) M.transpose(false);
			M.det = determinantLaplaceR2(M, activec, M.N);
			if (!M.rowAspectSparsest) M.transpose(false);	// transpose back matrix
			break;
		
		case 1:	M.det = determinantLaplaceR1(M, false); break;
		
		default:
		}
		return M.det;
	}

	
	// determinant Laplace Expansion Recursive algorithm with row-column elimination
	// (a slow algorithm, but simple to analyse as it creates a new reduced matrix on each recursion)
	public static double determinantLaplaceR1(Matrix M, boolean doZeroTest) {
		
		// counts number of times we've entered the method recursively
		Matrix.detL_DEBUG++;
		double[] data = M.getDataRef()[0];
		// base case: 2x2 matrix solves directly, end of recursion
		if (M.N == 2)
			return data[0] * data[3] - data[1] * data[2];
		
		double det = 0;		// determinant is accumulated here
		// basecase: 3x3 matrix solves directly, end of recursion
		if (M.N == 3) {
			if (!nearZero(data[0])) det += data[0] * (data[4] * data[8] - data[5] * data[7]);
			if (!nearZero(data[1])) det -= data[1] * (data[3] * data[8] - data[5] * data[6]);
			if (!nearZero(data[2])) det += data[2] * (data[3] * data[7] - data[4] * data[6]);
			return det;
		}

		// check if submatrix can return zero determinant immediately, but don't do this test at first entry
		if (doZeroTest && M.hasZeroDeterminant(false)) return 0;

		// for every column
		for (int j = 0; j < M.N; j++) {
			// if multiplicand is not zero (this test makes sparse matrices very fast to calculate)
			if (!nearZero(data[j])) {
				// recursively call determinantLaplace() on row/column-eliminated submatrices
				det +=  (((j + 2) & 1) == 0 ? 1 : -1) * data[j] * determinantLaplaceR1(M.eliminateRowColumn(0, j, false), true);
			}
		}
		return det;
	}

	
	
	// version two of determinant Laplace Expansion does not use the costly matrix reconstruction method
	// eliminateRowColumn(), but instead readjusts a relative increment list of active columns of each recursive submatrix,
	// activec contains a list of relative increment indexes, telling how much to increment to get to next active column
	public static double determinantLaplaceR2(Matrix M, int[][] activec, int columns) {
		
		Matrix.detL_DEBUG++;
		double[] data = M.getDataRef()[0];
		int tofs1, tofs2, topb = M.M - columns;
		
		// get 2x2 base case row offsets tofs1 & tofs2 from detAnalysis if it has been created
		if (M.analysis != null) {		
			int[] zRIdx = M.analysis[2];
			tofs1 = zRIdx[topb] * M.N;
			tofs2 = zRIdx[topb + 1] * M.N;
			
			if (columns == 2) {
				int ac = activec[0][0], tofsac = tofs1 + ac - 1, tofsac2 = tofs2 + ac - 1;
				return data[tofsac] * data[tofsac2 + activec[0][ac]] - data[tofsac + activec[0][ac]] * data[tofsac2];
			}
		} else {
			tofs1 = topb * M.N;	
			// base case: 2x2 matrix solves directly, end of recursion
			if (columns == 2) {
				int ac = activec[0][0], tofsac1 = tofs1 + ac - 1, tofsac2 = tofs1 + ac + activec[0][ac] - 1;
				return data[tofsac1] * data[M.N + tofsac2] - data[tofsac2] * data[M.N + tofsac1];
			}
		}
			
		// don't do the test at first entry (already done), otherwise test if submatrix conforms to zero determinant heuristics
		if (topb > 0 && M.hasZeroDeterminant2(activec, columns)) return 0;
		
		double det = 0;
				
		// walk the relative increment list that skips inactive columns, ac_prev = previous index, ac = current index
		// the j is independently counting up the column number
		for (int j = 0, ac_prev = 0, ac = activec[0][0]; j < columns; j++, ac_prev = ac, ac += activec[0][ac]) {

			// readjust increment list to skip current column in next recursion
			int oldac = activec[0][ac_prev];
			activec[0][ac_prev] += activec[0][ac];
			
			double v = data[tofs1 + ac - 1];
			// don't care recursing if we're multiplying with a zero value
			if (!nearZero(v))
				det += (((j + 2) & 1) == 0 ? 1 : -1) * v * determinantLaplaceR2(M, activec, columns - 1);

			// restore increment list for the higher recursive stage
			activec[0][ac_prev] = oldac;
		}
		return det;
	}
	
	

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			EIGENMETHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	// the Power Method will find the dominating eigenvalue and eigenvector of matrix A, given initial guess vector y
	// eigenvalue is returned by method, eigenvector is returned in y
	public double eigenPowerValue(Matrix y, double thetaError, int iters) {
		
		if (y.M < 2 || y.N != 1) 	throw new InvalidParameterException("Matrix.eigenPowerValue(): Invalid guess vector.");
		if (M != N)					throw new InvalidParameterException("Matrix.eigenPowerValue(): Matrix not square.");
		if (N != y.M)				throw new InvalidParameterException("Matrix.eigenPowerValue(): Matrix-vector dimensionality mismatch.");

		if (DEBUG_LEVEL > 1) System.out.println("Matrix to search: \n" + toString());
		
		int debug_lvl = DEBUG_LEVEL, k;
		DEBUG_LEVEL = 0;
		
		Matrix v = y.clone();
		double theta = 0, oldtheta = 0;
		iters = iters < 1 ? 1 : iters;
		
		// do maximum "iters" iterations
		for (k = 0; k < iters; k++) {
	
			v = this.multiply(y);
			oldtheta = theta;
			theta = v.absLength(0);
			if (DEBUG_LEVEL > 1) System.out.println(theta);
			y = v.multiply(1.0 / theta, false);

			if (oldtheta - theta >= -thetaError && oldtheta - theta <= thetaError)
				break;		// attained eigenvalue precision?
		}
		DEBUG_LEVEL = debug_lvl;
		if (DEBUG_LEVEL > 1) {
			System.out.println("eigenPowerMethod() converged in " + k + " iterations.");
			System.out.println(y.toString());
		}
		return theta;
	}
	
	
	
	// findEigenQR() uses the common A(k+1) = Q*AQ iterative diagonalisation routine, and Q is an orthonormal matrix
	public Matrix findEigenQR(Matrix Q, int iters, int vIters, double maxErr) {

		DEBUG_LEVEL--;
		Matrix Alam = this.clone(), QT = Q.transpose(true);
		DEBUG_LEVEL++;
		Alam = Alam.reduceHouseholder(false);		// A matrix reduced to Householder form reduces the iterations
		DEBUG_LEVEL--;
		
		double currErr = 0, lastErr = 0;
		double[] dataAlam = Alam.data, lastD = new double[M];
		maxErr *= maxErr;							// compare squared values directly to avoid unnecessary square root
		
		// QR iteration will reduce matrix A to eigenvalues on the diagonal
		int i = 0;
		for (; i < iters; i++) {

			Alam = (QT.multiply(Alam)).multiply(Q);			// A(k+1) = Q*(k)A(k)Q(k)
			Q = Alam.orthogonalise(true);
			QT = Q.transpose(true);
			
			dataAlam = Alam.data;
			lastErr = currErr;
			currErr = 0;
			for (int l = 0, j = 0; j < M; j++, l+= N + 1) {
				double d = lastD[j] - dataAlam[l];
				currErr += d * d;							// get residual vector length
				lastD[j] = dataAlam[l];						// save diagonal for lambda values error comparison
			}
			//currAbs = Math.sqrt(currAbs);
			System.out.println("Error: " + String.format("%.4f, delta: %.4f", currErr, lastErr-currErr));
			if (currErr < maxErr) break;
		}
		int itersDone = i;
		
		// the eigenvectors to solve start out filled with ones
		Matrix[] ev = new Matrix[M];
		for (i = 0; i < M; i++)
			ev[i] = new Matrix("z", M, 1, Matrix.Type.Unit, 1);

		// solve the eigenvectors
		Matrix Ilam = identity(M, 1);
		for(int j = 0; j < M; j++) {
			// this is the I*lambda matrix
			for (int iM = M+1, d = 0, jMj = j*M+j, dEnd = d*M; d < dEnd; d += iM) Ilam.data[d] = dataAlam[jMj];
			// LU decompose (A - lambda*I)
			Matrix[] lLU = this.subtract(Ilam, true).decomposeLU(true, true, true);
			Matrix[] xl = {ev[j], null};
			// iterate 4 times: LU*x(l+1) = x(l)
			for (int l = 0; l < 4; l++) {
				xl = xl[0].backSubstituteLU(null, lLU[1], false);
				if (xl == null) { ev[j] = null; break; }			// if backsubstitution fails, set that eigenvector to null
				xl[0].normalise(false);
			}
		}
		
		DEBUG_LEVEL++;

		if (DEBUG_LEVEL > 1) {
			System.out.println("QR eigenfinder took " + itersDone + " iterations, results:");
			for(i = 0; i < M; i++) {
				System.out.println("eigenvalue" + (i+1) + ": " + Alam.data[i * M + i]);
				System.out.println("eigenvector" + (i+1) + (ev[i] != null ? ":" : " failed converging."));
				if (ev[i] != null) {
					for(int r = 0; r < M; r++) {
						double v = ev[i].data[r];
						System.out.println((v < 0 ? "|" : "| ") + String.format("%.5f", v) + " |");
					}
				}
			}
			if (DEBUG_LEVEL > 2) {
				System.out.println(Alam.toString());
				System.out.println(Q.toString());
			}
		}
		return Alam;
	}

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			DATA HANDLING METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	public double[][] getDataRef() {
		double[][] dataSet = new double[2][];
		dataSet[0] = data;
		dataSet[1] = idata;
		return dataSet;
	}
	
	public double[][] getData() {
		double[][] dataSet = new double[2][];
		dataSet[0] = (data != null ? data.clone() : null);
		dataSet[1] = (idata != null ? idata.clone() : null);
		return dataSet;
	}
	
	public void putDataRef(double[] data, double[] idata) { this.data = data; this.idata = idata; clearNull(); }
	
	public void putData(double[] data, double[] idata) {
		this.data = (data != null ? data.clone() : null);
		this.idata = (idata != null ? idata.clone() : null);
		if (data != null) clearNull();
	}

	
	@Override
	public Matrix clone() {
		Object O = null;
		try { O = super.clone(); } catch (CloneNotSupportedException e) { e.printStackTrace(); }
		Matrix A = (Matrix) O;
		if (data != null) { A.data = data.clone(); clearNull(); }
		if (idata != null) A.idata = idata.clone();
		if (bitImage != null) A.bitImage = bitImage.clone(A);
		A.mutator = new int[2][];
		if (mutator[0] != null) A.mutator[0] = mutator[0].clone();
		if (mutator[1] != null) A.mutator[1] = mutator[1].clone();
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
		return A;
	}
	

	// this clone operation will deconvert from inheritor classes (CSR & NSP) into row-column through overridden getData() method
	public Matrix toMatrix() {
		double[][] dataSet = getDataRef();
		Matrix A = new Matrix(name, M, N, dataSet[0], dataSet[1]);
		A.status = status;
		if (bitImage != null) A.bitImage = bitImage.clone(A);
		A.mutator = new int[2][];
		if (mutator[0] != null) A.mutator[0] = mutator[0].clone();
		if (mutator[1] != null) A.mutator[1] = mutator[1].clone();
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(A.toString());
		return A;
	}
	

	public void zero() {
		if (data != null) for (int i = 0, MN = M * N; i < MN; i++) data[i] = 0;
		if (idata != null) for (int i = 0, MN = M * N; i < MN; i++) idata[i] = 0;
		bitImage.zero();
		if (mutator[0] != null) mutator[0] = new int[M];
		if (mutator[1] != null) mutator[1] = new int[M];
		setNull();
		halfBandwidth = -1;
	}

	// method wipes out a specific row and column
	public Matrix eliminateRowColumn(int r, int c, boolean makeBitImage) {
		
		if (r < 0 || r > M - 1 || c < 0 || c > N - 1)
			throw new InvalidParameterException("Matrix.eliminateRowColumn(): Row or column out of bounds.");
		if (M < 2 || N < 2)
			throw new InvalidParameterException("Matrix.eliminateRowColumn(): Invalid matrix size.");

		String newname;
		if (DEBUG_LEVEL > 1) 	newname = new String(name + "(M-" + r + ",N-" + c + ")");
		else 					newname = name;
		Matrix A = null;
		if (isComplex())	A = new Matrix(newname, M - 1, N - 1, Matrix.Type.Null_Complex, 1);
		else				A = new Matrix(newname, M - 1, N - 1, Matrix.Type.Null, 1);
		
		for (int i = 0, ii = 0; i < M; i++) {
			if (i != r) {
				int iiAN = ii * A.N, iN = i * N;
				// real-valued matrix case
				if (!isComplex()) {
					for (int j = 0, jj = 0; j < N; j++)
						if (j != c) A.data[iiAN + jj++] = data[iN + j];
				// complex-valued matrix case
				} else {
					for (int j = 0, jj = 0; j < N; j++)
						if (j != c) {
							A.data[iiAN + jj] = data[iN + j];	
							A.idata[iiAN + jj++] = idata[iN + j]; }
				}
				ii++;
			}
		}

		if (Matrix.DEBUG_LEVEL > 1) System.out.println(this.toString());
		if (makeBitImage) A.bitImage = new BinBitImage(A);
		return A;
	}

	
	// method directly inverts 3x3 matrix with determinant returned at end of array
	public static double[] invert3x3(double[] m3x3) {
		double[] m3x3i = new double[10];
		double v3 = m3x3[3], v4 = m3x3[4], v5 = m3x3[5], v6 = m3x3[6], v7 = m3x3[7], v8 = m3x3[8];
		m3x3i[0] = v4 * v8 - v5 * v7;
		m3x3i[3] = -(v3 * v8 - v5 * v6);
		double det = m3x3i[6] = v3 * v7 - v4 * v6;
		double v0 = m3x3[0], v1 = m3x3[1], v2 = m3x3[2];
		det = v0 * m3x3i[0] + v1 * m3x3i[3] + v2 * det;	
		if (Matrix.nearZero(det)) return null;										// signal singular matrix with null return
		double detD  = 1.0 / det;
		m3x3i[0] *= detD; m3x3i[3] *= detD; m3x3i[6] *= detD;
		m3x3i[1] = -(v1 * v8 - v2 * v7) * detD;
		m3x3i[2] = (v1 * v5 - v2 * v4) * detD;
		m3x3i[4] = (v0 * v8 - v2 * v6) * detD;
		m3x3i[5] = -(v0 * v5 - v2 * v3) * detD;
		m3x3i[7] = -(v0 * v7 - v1 * v6) * detD;
		m3x3i[8] = (v0 * v4 - v1 * v3) * detD;
		m3x3i[9] = det;																// also stores down determinant at end of data
		return m3x3i;
	}
	
	// this method expects stepwise calling, from two methods processing same matrix, where one only needs determinant, the other a complete inverse
	// preserving the initial data of determinant calculation saves on FLOPs, so the method assumes that either:
	//			a) no calculations have been done (m3x3i = null), calculate determinant, return it, expecting second call for the rest
	//			b) part of calculation pertaining to finding determinant was already done, the inverse needs finishing
	public static double[] invert3x3dualCall(double[] m3x3, double[] m3x3i) {
		double v0, v1, v2;
		double v3 = m3x3[3], v4 = m3x3[4], v5 = m3x3[5], v6 = m3x3[6], v7 = m3x3[7], v8 = m3x3[8];
		if (m3x3i == null) {
			m3x3i = new double[11];
			m3x3i[0] = v4 * v8 - v5 * v7;
			m3x3i[3] = -(v3 * v8 - v5 * v6);
			double det = m3x3i[6] = v3 * v7 - v4 * v6;
			v0 = m3x3[0]; v1 = m3x3[1]; v2 = m3x3[2];
			det = v0 * m3x3i[0] + v1 * m3x3i[3] + v2 * det;	
			if (Matrix.nearZero(det)) return null;						// signal singular matrix with null return			
			m3x3i[9] = det;												// also stores down determinant at end of data
			return m3x3i;
		} else {
			v0 = m3x3[0]; v1 = m3x3[1]; v2 = m3x3[2];
		}
		double detD  = 1.0 / m3x3i[9];
		m3x3i[0] *= detD; m3x3i[3] *= detD; m3x3i[6] *= detD;
		m3x3i[1] = -(v1 * v8 - v2 * v7) * detD;
		m3x3i[2] = (v1 * v5 - v2 * v4) * detD;
		m3x3i[4] = (v0 * v8 - v2 * v6) * detD;
		m3x3i[5] = -(v0 * v5 - v2 * v3) * detD;
		m3x3i[7] = -(v0 * v7 - v1 * v6) * detD;
		m3x3i[8] = (v0 * v4 - v1 * v3) * detD;
		m3x3i[10] = 1;													// signal that matrix is fully inverted
		return m3x3i;
	}


	// method returns determinant of 3x3 matrix
	public static double determinant3x3(double[] m3x3) {
		double v3 = m3x3[3], v4 = m3x3[4], v5 = m3x3[5], v6 = m3x3[6], v7 = m3x3[7], v8 = m3x3[8];
		double m0 = v4 * v8 - v5 * v7;
		double m3 = -(v3 * v8 - v5 * v6);
		double det = v3 * v7 - v4 * v6;
		return m3x3[0] * m0 + m3x3[1] * m3 + m3x3[2] * det;
	}

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			INLINE/SHORTFORM METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	public static boolean nearZero(double v) { return (v < -ROUNDOFF_ERROR || v > ROUNDOFF_ERROR ? false : true); }
	public static boolean nearZeroC(double[] iv) {
		return (	iv[0] < -ROUNDOFF_ERROR || iv[0] > ROUNDOFF_ERROR ||
					iv[1] < -ROUNDOFF_ERROR || iv[1] > ROUNDOFF_ERROR ? false : true);
	}
	public static boolean nearZeroC(double v, double iv) {
		return (	v < -ROUNDOFF_ERROR || v > ROUNDOFF_ERROR ||
					iv < -ROUNDOFF_ERROR || iv > ROUNDOFF_ERROR ? false : true);
	}
	boolean isNull() { return (status & NULL_MATRIX) != 0; }
	void clearNull() { status ^= NULL_MATRIX & status; }
	void setNull() { status |= NULL_MATRIX; }
	boolean isComplex() { return (status & COMPLEX_MATRIX) != 0; }
	void setComplex() { status |= COMPLEX_MATRIX; }
	
	// get complex modulus of imaginary value
	double complexM(double r, double i) { return Math.sqrt(r * r + i * i); };
	
	public void readjustHalfBandwidth(int r, int c) {
		int width = r - c;								// readjust half bandwidth if value coordinates are "sticking out" from diagonal
		if (width < 0) width = -width;
		if (halfBandwidth < width) halfBandwidth = width;
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// limits output of a double to max 5 characters + sign, for following cases:
	// -XXeEi		v >= 10000
	// -XXXXi		v >= 100
	// -XX.Di		v >= 10
	// -X.DDi		v < 10
	static String to5chars(double v, boolean complex) {
		int decimals = 0;
		if (v >= 10000 || v <= -10000) {
			for (; v >= 100 || v <= -100; v /= 10.0) decimals++;
			if (decimals >= 10) { v /= 10.0; decimals++; }
			if (complex)
					return (v < 0 ? " " : "  ") + (int)v + "e" + decimals + "  ";
			else	return (v < 0 ? " " : "  ") + (int)v + "e" + decimals + (complex ? "i " : " ");
		}
		
		boolean isInt = nearZero(v - (int)v);
			
		if (v >= 100 || v <= -100) {
			return complex ? String.format("%6di", (int)v) : String.format("%6d ", (int)v);
			//return complex ? String.format("%6.0fi", v) : String.format("%6.0f ", v);
		}
			
		if (v >= 10 || v <= -10) {
			if (isInt) return complex ? String.format("%5di ", (int)v) : String.format("%5d  ", (int)v);
			return complex ? String.format("%6.1fi", v) : String.format("%6.1f ", v);
		}
			
		if (isInt) return complex ? String.format("%4di  ", (int)v) : String.format("%4d   ", (int)v);
		return complex ? String.format("%6.2fi", v) : String.format("%6.2f ", v);
	}
	
	
	
	static int MAX_PRINTEXTENT = 50;
	
	@Override
	public String toString() {
		
		double[][] dataSet = this.getDataRef();
		double[] data = dataSet[0], idata = dataSet[1];
		
		StringBuilder sb = new StringBuilder();
		int maxM = M > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : M;
		int maxN = N > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : N;
		
		sb.append(((M == 1 || N == 1) ? "vector: " :"matrix: ") + name + (M == 1 ? "^T\n" : "\n"));
		sb.append("size: " + M + ", " + N);
		
		if (data != null) {
			
			// write out column indexes, if they were permuted
			if (mutator[1] != null) {
				sb.append(mutator[0] != null ? "\n        " : "\n ");
				for (int j = 0; j < maxN; j++) sb.append(to5chars(mutator[1][j], false));
			}
			
			for (int i = 0; i < maxM; i++) {
				int iN = i * N;
				if (sb.length() > 0) sb.append("\n");		// don't add newline at start of matrix printout
				
				// write out row index, if it was permuted
				sb.append(mutator[0] != null ? to5chars(mutator[0][i], false) + "|" : "|");
				
				for (int j = 0; j < maxN; j++)
					if (nearZero(data[iN + j])) 	sb.append("   -   ");
					else							sb.append(to5chars(data[iN + j], false));

				// if matrix was bigger than allowed printout bounds, indicate the continuation
				if (maxN < N) if (i % 4 == 0) sb.append(" ..."); else sb.append("    ");
				
				// any imaginary data comes as a second line under the real data line
				if (idata != null) {
					sb.append("     |\n|");
					for (int j = 0; j < maxN; j++)
						if (nearZero(data[iN + j]))	sb.append("       ");
						else						sb.append(to5chars(idata[iN + j], true));
					
					sb.append(" |\n|");
					for (int j = 0; j < maxN; j++) sb.append("       ");
					if (maxN < N) sb.append("    ");
				}
				
				sb.append(" |");
			}
			// if matrix was bigger than allowed printout bounds, indicate the continuation
			if (maxM < M) {
				for (int i = 0; i < 3; i++) {
					sb.append("\n|");
					for (int j = 0; j < maxN; j++)
						if (j%2 == 0) sb.append("   .   "); else sb.append("       ");
					sb.append("     |");
				}
			}
		} else sb.append("[null]\n");
		if (M == N) {
			sb.append("\ndeterminant: " + det + "\n");
			if (halfBandwidth > 0)
				sb.append("half bandwidth: " + halfBandwidth + "\n");
		}
		else		sb.append("\n");
			
		return sb.toString();	
	}
	
	
	// method writes matrix to file in full rectangular form
	// if precision > 0 the no. of decimals is = precision, otherwise values are written in condensed 5 char mode
	public void toFile(int precision) {
		
		File file = new File(name + ".txt");
		if (!file.exists()) {
			try {	file.createNewFile();
			} catch (IOException e) { e.printStackTrace(); }
		}
		BufferedWriter bw = null;
		try {		bw = new BufferedWriter(new FileWriter(file));
		} catch (IOException e) { e.printStackTrace(); }

		double[][] dataSet = this.getDataRef();
		double[] data = dataSet[0], idata = dataSet[1];
		
		StringBuilder sb = new StringBuilder();
		
		sb.append(((M == 1 || N == 1) ? "vector: " :"matrix: ") + name + (M == 1 ? "^T\n" : "\n"));
		sb.append("size: " + M + ", " + N);
		
		if (data != null) {
			
			// write out column indexes, if they were permuted
			if (mutator[1] != null) {
				sb.append(mutator[0] != null ? "\n        " : "\n ");
				for (int j = 0; j < N; j++) sb.append(to5chars(mutator[1][j], false));
			}

			for (int i = 0; i < M; i++) {
				int iN = i * N;
				if (sb.length() > 0) sb.append("\n");		// don't add newline at start of matrix printout
				
				// write out row index, if it was permuted
				sb.append(mutator[0] != null ? to5chars(mutator[0][i], false) + "|" : "|");

				for (int j = 0; j < N; j++)
					if (nearZero(data[iN + j])) sb.append("   -   ");
					else	sb.append(precision > 0 ?
								String.format("%." + precision + "f", data[iN + j]) :
								to5chars(data[iN + j], false));

				// any imaginary data comes as a second line under the real data line
				if (idata != null) {
					sb.append("     |\n|");
					for (int j = 0; j < N; j++)
						if (nearZero(data[iN + j])) sb.append("       ");
						else	sb.append(precision > 0 ?
									String.format("%." + precision + "f", data[iN + j]) :
									to5chars(idata[iN + j], true));
					
					sb.append(" |\n|");
					for (int j = 0; j < N; j++) sb.append("       ");
				}
				
				sb.append(" |");
			}
		} else sb.append("[null]\n");
		if (M == N) {
			sb.append("\ndeterminant: " + det + "\n");
			if (halfBandwidth > 0)
				sb.append("half bandwidth: " + halfBandwidth + "\n");
		}
		else		sb.append("\n");

		try {
			bw.write(sb.toString());
			bw.flush();
			bw.close();
		} catch (IOException e) { e.printStackTrace(); }

	}

	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			MATRIX TYPES
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	public enum Type {
		Null(0), Identity(1), Random(2), Centering(3), Unit(4), UpperHessenberg(5), Null_Complex(16), Random_Complex(17);
		private int type;
		private Type(int type) { this.type = type; }
	}
}

