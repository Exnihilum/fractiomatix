package lkr74.matrixlib;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;


public class Matrix implements Cloneable {

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			FIXED VALUES
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	static final double ROUNDOFF_ERROR = 1e-8;
	static final int IDENTITY_MATRIX = 1;
	static final int RC_MATRIX = 1<<1;
	static final int CSR_MATRIX = 1<<2;
	static final int SINGULAR_MATRIX = 1<<8;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			INSTANCE-LEVEL VALUES
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	protected String name;
	protected int M, N; 						// number of rows & columns
	protected double[] data;					// row-column matrix uses a flat array
	protected int status;
	public double det = 0;						// last calculated determinant value

	// optimisation hierarchies & variables
	// binBitLength = no. of values gathered in binBit, procNum = number of concurrent processors available
	protected boolean rowAspectSparsest = true;	// indicates whether row or column aspect is most sparse for determinant analysis
	protected int detSign = 1;					// tells sign of the determinant, according to heuristics of analyseRowsColumns() 
	protected int[][] detAnalysis = null;		// used by analyseRowsColumns() for determinant finding heuristics
	protected BinBitImage bitImage;
	
	// multithreading variables
	static ExecutorService executor;				// thread pool for matrix ops
	static List<Future<Double>> futureDoubleList;	// list of return doubles from Callable threads
	protected boolean threaded = false;				// indicates if matrix should/can be run in multithreaded mode	
	static int taskNum, procNum;
	static Thread taskList[];
	volatile double vdet = 0;						// threadsafe determinant accumulator for multithreaded methods
	
	// global variables
	protected static int DEBUG_LEVEL = 2;
	protected volatile static int nameCount = 1;
	// debugging global variables
	protected static volatile int detL_DEBUG;
	protected static volatile int mulSW_DEBUG;
	protected static volatile int mulFlops_DEBUG;
	protected static volatile int mulFlopsSW_DEBUG;
	protected static volatile int mulAdopsSW_DEBUG;
	
	// skeleton instantiator of a matrix
	public Matrix(String name, int r, int c) {
		// initialise status flags, all bits = 0 means an empty matrix of type Matrix
		if (r < 1 || c < 1) throw new RuntimeException("Matrix(): Illegal matrix dimensions.");
		this.name = name;
		this.M = r;
		this.N = c;
	}

	//	instantiates Matrix with data of choice
	public Matrix(String name, int r, int c, Type type) {
		// initialise status flags, all bits = 0 means an empty matrix of type Matrix
		if (r < 1 || c < 1) throw new RuntimeException("Matrix(): Illegal matrix dimensions.");
		this.name = name;
		this.M = r;
		this.N = c;
		data = Matrix.generateData(M, N, type);
		bitImage = new BinBitImage(this);
	}
	
	//	instantiates a matrix with a provided dataset, cloning the dataset into this matrix
	public Matrix(String name, int M, int N, double[] data) {
		if (M < 1 || N < 1 || data.length < M * N)
			throw new RuntimeException("Matrix(): Illegal matrix dimensions.");
		this.name = name;
		this.M = M;
		this.N = N;
		this.putData(data);
		bitImage = new BinBitImage(this);
	}

	
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

	public void makeThreaded() { this.threaded = true; }		// activate multithreading flag for this matrix
	public void makeNonThreaded() { this.threaded = false; }	// deactivate multithreading flag for this matrix

	public double[] getDataRef() { return data; }
	public double[] getData() { return data.clone(); }
	public void putDataRef(double[] data) { this.data = data; }
	public void putData(double[] data) { this.data = data.clone(); }
	
	
	// polymorphic clone operation will deconvert from other formats into row-column through overridden getData() method
	@Override
	public Matrix clone() {
		Object O = null;
		try { O = super.clone(); } catch (CloneNotSupportedException e) { e.printStackTrace(); }
		Matrix A = (Matrix) O;
		A.name = name;
		A.M = M;
		A.N = N;
		A.status = status;
		if (data != null) A.putDataRef(getData());
		A.bitImage = bitImage.clone(A);
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
		return A;
	}

	public void zero() {
		for (int i = 0, MN = M * N; i < MN; i++) data[i] = 0;
		bitImage.zero();
	}

	public Matrix eliminateRowColumn(int r, int c, boolean makeBitImage) {
		
		if (r < 0 || r > M - 1 || c < 0 || c > N - 1)
			throw new RuntimeException("Matrix.eliminateRowColumn(): Row or column out of bounds.");
		if (M < 2 || N < 2)
			throw new RuntimeException("Matrix.eliminateRowColumn(): Invalid matrix size.");

		String newname;
		if (DEBUG_LEVEL > 1) 	newname = new String(name + "(M-" + r + ",N-" + c);
		else 					newname = name + nameCount++;
		Matrix A = new Matrix(newname, M - 1, N - 1, Matrix.Type.Null);
		
		for (int i = 0, ii = 0; i < M; i++) {
			if (i != r) {
				for (int j = 0, jj = 0; j < N; j++)
					if (j != c)
						A.data[ii * A.N + jj++] = data[i * N + j];
				ii++;
			}
		}
		if (Matrix.DEBUG_LEVEL > 1) System.out.println(this.toString());
		A.bitImage = new BinBitImage(A);
		return A;
	}

	public double valueOf(int r, int c) { return data[r * N + c]; }

	public void valueTo(int r, int c, double v) {
		data[r * N + c] = v;
		bitImage.setBit(r, c);				// set the corresponding bit in bitImage
	}

	public Matrix identity(int s) { return new Matrix("I", s, s, Type.Identity); }
	
	public static Matrix center(Matrix A) {
		Matrix At = A.clone();
		At.transpose(false);
		Matrix C = new Matrix("C", At.M, At.N,  Matrix.Type.Centering);
		Matrix Ac = Matrix.multiply(C, At);
		Ac.transpose(false);
		Ac.bitImage = new BinBitImage(Ac);
		return Ac;
	}

	
	public static double[] generateData(int r, int c, Type type) {
		if (r < 1 || c < 1)
			throw new RuntimeException("Matrix.generateData(): Illegal matrix dimensions.");
		double[] data = new double[r*c];
		
		switch (type) {
			case Null:
				break;
			case Identity:
				if (r != c) 	throw new RuntimeException("Matrix.generateData(): Identity matrix must be square.");
				for (int i = 0; i < c; i++)
					data[i * c + i] = 1;
				break;
			case Random:
				for (int i = 0; i < r; i++)
					for (int j = 0; j < c; j++)
						data[i * c + j] = Math.random();
				break;
			case Centering:
				if (r != c) 	throw new RuntimeException("Matrix.generateData(): Centering matrix must be square.");
				for (int i = 0; i < r; i++)
					for (int j = 0; j < c; j++) {
						if (i == j)		data[i * c + j] = 1.0 - 1.0 / (double)r;
						else			data[i * c + j] = -1.0 / (double)r;
					}
				break;
			default:
				throw new RuntimeException("Matrix.generateData(): Illegal matrix type.");
		}
		return data;
	}

	
	public Matrix rescale(int Mi, int Ni, int Mo, int No, boolean doBitImage) {
		
		if (Mi > Mo || Ni > No) throw new RuntimeException("Matrix.getData(): Invalid data subset size.");
		int newM = Mo - Mi, newN = No - Ni;
		
		String newname;
		if (DEBUG_LEVEL > 1)	newname = name + "(r:" + Mi + "," + Ni + "," + Mo + "," + No + ")";
		else					newname = name + "(r)";
		int dbg = Matrix.DEBUG_LEVEL;
		Matrix.DEBUG_LEVEL = 0;
		Matrix A = new Matrix(newname, newM, newN, Matrix.Type.Null);
		Matrix.DEBUG_LEVEL = dbg;
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
		if (DEBUG_LEVEL > 1) System.out.println(A.toString());
		if (doBitImage) A.bitImage = new BinBitImage(A);
		return A;
	}

	
	
	// finds out number of zeroes in each row & column, sums them up into an analysis array, useful for finding determinant
	public void analyseRowsColumns() {
		
		if (M < 2 || N < 2) {
			if (M < 1 || N < 1) throw new RuntimeException("Matrix.analyseRowsColumns(): Invalid matrix dimensions.");
			return;		// 1x1 matrix, nothing to analyse
		}

		double[] data = getDataRef();
		detSign = 1;
		rowAspectSparsest = true;
		int signR = 1, signC = 1;
		// zRCount & zCCount store number of zeroes in the rows and columns
		if (detAnalysis == null) detAnalysis = new int[4][];
		detAnalysis[0] = new int[M];			// holds zeroes counts for each row
		detAnalysis[1] = new int[N];			// holds zeroes counts for each column
		detAnalysis[2] = new int[M];			// holds indexes into zCount[0], sorted according to most zeroes
		detAnalysis[3] = new int[N];			// holds indexes into zCount[1], sorted according to most zeroes
		
		int[] zRCount = detAnalysis[0], zCCount = detAnalysis[1], zRIdx = detAnalysis[2], zCIdx = detAnalysis[3];
		// initialise the row & column zero-count indexation arrays
		for (int i = 0, iend = zRIdx.length; i < iend; i++) zRIdx[i] = i;
		for (int i = 0, iend = zCIdx.length; i < iend; i++) zCIdx[i] = i;
		
		double v;
		for (int i = 0; i < M; i++) {
			int iN = N * i;
			for (int j = 0; j < N; j++) {
				v = data[iN + j];
				if (v <= ROUNDOFF_ERROR && v >= -ROUNDOFF_ERROR) { zRCount[i]++; zCCount[j]++; }
			}
		}
		
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
			detAnalysis[0] = detAnalysis[1]; detAnalysis[1] = zRCount;
			detAnalysis[2] = detAnalysis[3]; detAnalysis[3] = zRIdx;
		}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("matrix " + this.name + "'s row & column analysis:");
			System.out.println(Arrays.toString(detAnalysis[0]));
			System.out.println(Arrays.toString(detAnalysis[1]));
			System.out.println(Arrays.toString(detAnalysis[2]));
			System.out.println(Arrays.toString(detAnalysis[3]) + "\n");
		}
	}
	
	
	
	// swap two rows of matrix
	public void swap(int r1, int r2) {
		double temp;
		for (int i = 0, or1 = r1 * N, or2 = r2 * N; i < N; i++) {
			temp = data[or1 + i];
			data[or1 + i] = data[or2 + i];
			data[or2 + i] = temp;
		}
		det = -det;						// determinant is inverted by row swapping
		if (DEBUG_LEVEL > 2) System.out.println("swap:\n" + this.toString());
		if (bitImage != null)
			bitImage.swapRows(r1, r2);		// swap bitImage rows
	}


	
	
	// transpose the matrix, returning data within same matrix or in a copied matrix
	public Matrix transpose(boolean copy) {
		
		if (M < 1 || N < 1) throw new RuntimeException("Matrix.transpose(): Invalid matrix dimensions.");

		Matrix T = this;
		
		// given a non-square matrix, reallocate new data field
		if (M != N) {
			if (copy) T = new Matrix(name + "^T", N, M);
			double[] newdata = new double[N * M];
			
			for (int i = 0; i < M; i++) {
				int iN = i * N;
				for (int j = 0; j < N; j++) newdata[j * M + i] = data[iN + j];
			}
			
			T.data = newdata;
			if (DEBUG_LEVEL > 1) System.out.println(T.toString());
			T.bitImage = new BinBitImage(T);
		} else {
			if (copy) { T = this.clone(); T.name = T.name + "^T"; }
			double[] dataT = T.data;
			// the approach for a square matrix
			double temp;
			for (int i = 0; i < M; i++) {
				int iN = i * N;
				for (int j = i + 1, jMpi; j < N; j++) {
					jMpi = j * M + i;
					temp = dataT[jMpi];
					dataT[jMpi] = dataT[iN + j];
					dataT[iN + j] = temp;
					T.bitImage.transposeBit(i, j);
				}
			}
			if (DEBUG_LEVEL > 1) System.out.println(T.toString());
		}	
		return T;
	}
	
	
	// normalises the matrix/vector/transposevector A, for a matrix the columns are normalised individually
	public static void normalise(Matrix A) {

		if (A.M < 1 || A.N < 1) throw new RuntimeException("Matrix.normalise(): Invalid matrix dimensions.");

		double[] data = A.data;
		if (A.M > 1)
			// if we're not dealing with a transposed vector
			for (int c = 0, N = A.N, M = A.M; c < N; c++) {
				double norm = 0, sq;
				for (int r = 0; r < M; r++) {
					sq = data[r * N + c];
					norm += sq * sq;	// Euclidean norm ||A||2
				}
				norm = Math.sqrt(norm);
				for (int r = 0; r < M; r++) data[r * N + c] /= norm;	// normalise
			}
		else {
			double norm = 0, sq;
			for (int c = 0, N = A.N; c < N; c++) {
				sq = data[c];
				norm += sq * sq;	// Euclidean norm ||A||2
			}
			norm = Math.sqrt(norm);
			for (int c = 0, N = A.N; c < N; c++) data[c] /= norm;		// normalise
		}
		if (DEBUG_LEVEL > 1) System.out.println(A.toString());
	}

	
	
	// polymorphic add will do A+B returning resultant C of matrix type A
	// the data from both matrices are gotten through the polymorphic method getDataRef()
	public static Matrix add(Matrix A, Matrix B) {
		
		if (B.M != A.M || B.N != A.N)
			throw new RuntimeException("Matrix.add(): Nonmatching matrix dimensions.");
		double[] dataA = A.getDataRef(), dataB = B.getDataRef(), dataC;
		dataC = new double[A.M * A.N];
		
		Matrix C = A.clone();
		if (Matrix.DEBUG_LEVEL > 1) C.name = "(" + A.name + "+" + B.name + ")";
		else						C.name = "A" + nameCount++;
		
		for (int i = 0; i < C.M; i++) {
			int iN = i * A.N;
			for (int j = 0; j < C.N; j++)
				dataC[iN + j] = dataA[iN + j] + dataB[iN + j];
		}
		C.putDataRef(dataC);
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		C.bitImage.make();
		return C;		
	}

	
	
	// polymorphic subtract will do A-B returning resultant C of matrix type A
	public static Matrix subtract(Matrix A, Matrix B) {
		
		if (B.M != A.M || B.N != A.N)
			throw new RuntimeException("Matrix.subtract(): Nonmatching matrix dimensions.");
		double[] dataA = A.getDataRef(), dataB = B.getDataRef(), dataC;
		dataC = new double[A.M * A.N];
		
		Matrix C = A.clone();
		if (Matrix.DEBUG_LEVEL > 1) C.name = "(" + A.name + "-" + B.name + ")";
		else						C.name = "S" + nameCount++;

		for (int i = 0; i < C.M; i++) {
			int iN = i * A.N;
			for (int j = 0; j < C.N; j++)
				dataC[iN + j] = dataA[iN + j] - dataB[iN + j];
		}
		C.putDataRef(dataC);
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		C.bitImage.make();
		return C;		
	}

	
	
	// get absolute length for a vector or transposevector
	public static double absLength(Matrix v) {
		double sq, vnorm = 0;
		if (v.N == 1) for (int r = 0, M = v.M; r < M; r++) { sq = v.data[r]; vnorm += sq * sq; }
		else if (v.M == 1) for (int c = 0, N = v.N; c < N; c++) { sq = v.data[c]; vnorm += sq * sq; }
		else throw new RuntimeException("Matrix.vnorm(): Invalid vector.");
		return Math.sqrt(vnorm);
	}

	
	
	// takes two vectors and returns a value
	public static double vmultiply(Matrix a, Matrix b) {
		
		if ((a.M > 1 && a.N != 1) || (a.N > 1 && a.M != 1))
			throw new RuntimeException("Matrix.vmultiply(): Invalid vector dimensions.");
		if ((a.M == b.M && a.N != b.N) || (a.M == b.N && a.N != b.M))
			throw new RuntimeException("Matrix.vmultiply(): Vectors do not match.");
		
		double c = 0;
		for (int i = 0, N = a.N; i < N; i++) c += a.data[i] * b.data[i];
		return c;
	}
	
	// polymorphic matrix product will do AB = C
	public static Matrix multiply(Matrix A, Matrix B) {

		if (A.M < 1 || A.N < 1 || B.M < 1 || B.N < 1) 	throw new RuntimeException("Matrix.multiply(): Invalid matrix dimensions.");
		if (A.N != B.M) 								throw new RuntimeException("Matrix.multiply(): Nonmatching matrix dimensions.");
		mulFlops_DEBUG = A.M * 2 + A.M * A.N;	// counts number of multiplications
		
		double[] dataA = A.getDataRef(), dataB = B.getDataRef();
		double[] dataC = new double[A.M * B.N];
		
		String newname;
		if (Matrix.DEBUG_LEVEL > 1) newname = "(" + A.name + "." + B.name + ")";
		else						newname = "P" + nameCount++;
		Matrix C = new Matrix(newname, A.M, B.N);

		for (int i = 0; i < A.M; i++) {
			int iN = i * A.N, iCN = C.N * i;
			for (int j = 0; j < A.N; j++) {
				int jBN = j * B.N;
				double v = dataA[iN + j];
				if (v < -ROUNDOFF_ERROR || v > ROUNDOFF_ERROR) { // optimisation for sparse matrices
					for (int k = 0; k < B.N; k++)
						dataC[iCN + k] += v * dataB[jBN + k];
					mulFlops_DEBUG += B.N;
				}
			}
		}		
		C.putDataRef(dataC);
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		C.bitImage = new BinBitImage(C);		// bitImage needs full reconstitution
		return C;
	}
	
	// multiplies/scales matrix with a value, bitimage will be unchanged (multiplying with zero is meaningless)
	public static Matrix multiply(double v, Matrix A) {
		if (A.M < 1 || A.N < 1)
			throw new RuntimeException("Matrix.multiply(): Invalid matrix dimensions.");
	
		double[] dataA = A.getDataRef();
		double[] dataC = new double[A.M * A.N];
		
		String newname;
		if (Matrix.DEBUG_LEVEL > 1) newname = "s*" + A.name;
		else						newname = "S" + nameCount++;
		Matrix C = new Matrix(newname, A.M, A.N);

		for (int i = 0, M = A.M; i < M; i++) {
			int iN = i * A.N;
			for (int j = 0, N = A.N; j < N; j++) dataC[iN + j] = v * dataA[iN + j];
		}
		
		C.putDataRef(dataC);
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		return C;
	}
	
	

	static double[][] buffer2x2StrasWin;				// preallocated buffer for the 2x2 submatrix cases
	
	// multiplier uses the Strassen-Winograd algorithm for matrix multiplication, adapted for flat data arrays
	// truncP = the size of submatrix at which recursion stops and the normal multiplicator takes over
	public static Matrix multiplyStrasWin(Matrix A, Matrix B, int truncP) {
		
		if (A.M < 1 || A.N < 1 || B.M < 1 || B.N < 1) throw new RuntimeException("Matrix.multiply(): Invalid matrix dimensions.");
		if (A.N != B.M) throw new RuntimeException("Matrix.multiply(): Nonmatching matrix dimensions.");
		mulSW_DEBUG = mulFlopsSW_DEBUG = mulAdopsSW_DEBUG = 0;
				
		// squarify and expand matrices to nearest 2^n size
		int majorM = A.M > A.N ? A.M : A.N, newM;
		for (newM = 2; newM < majorM; newM <<= 1);		// do 2^(n+1) until we reach or excel major size of the matrices
		if (A.M != A.N || newM != majorM) {
			A = A.rescale(0, 0, newM, newM, false);		// no BitImage needed for transient data
			B = B.rescale(0, 0, newM, newM, false);
		}
		if (truncP > newM) truncP = newM;				// truncation point can't be larger than the matrix

		String newname;
		if (Matrix.DEBUG_LEVEL > 1) newname = "(" + A.name + "." + B.name + ")";
		else						newname = "W" + Matrix.nameCount++;
		Matrix C = new Matrix(newname, newM, newM, Matrix.Type.Null);
		
		buffer2x2StrasWin = new double[8][truncP > 16 ? 16*16 : truncP*truncP];
		
		// subInfo carries along following data:
		// (truncP) the matrix size truncation point when standard multiplicator is used instead
		// (dim) the current submatrix dimension 
		// 3x offsets into the flat datafields, 3x the matrix width of the datafields
		int[] subInfo = {truncP, newM, 0, 0, 0, newM, newM, newM};
		multiplyStrasWin2(A.data, B.data, C.data, subInfo);
			
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		return C;

	}
	
	public static void multiplyStrasWin2(double[] dA, double[] dB, double[] dC, int[] subInfo) {

		int truncP = subInfo[0], dim = subInfo[1];
		int offsA = subInfo[2], offsB = subInfo[3], offsC = subInfo[4];
		int dimA = subInfo[5], dimB = subInfo[6], dimC = subInfo[7];
		Matrix.mulSW_DEBUG++;
		Matrix.mulFlopsSW_DEBUG += 4;
				
		// if we reached base case of 2x2, return straight 2x2 multiplication
		if (dim == 2) {
			double a11 = dA[offsA], a12 = dA[offsA + 1], a21 = dA[offsA + dimA], a22 = dA[offsA + dimA + 1];
			double b11 = dB[offsB], b12 = dB[offsB + 1], b21 = dB[offsB + dimB], b22 = dB[offsB + dimB + 1];
			dC[offsC] = a11 * b11 + a12 * b21;
			dC[offsC++ + dimC] = a21 * b11 + a22 * b21;
			dC[offsC] = a11 * b12 + a12 * b22;
			dC[offsC + dimC] = a21 * b12 + a22 * b22;
			Matrix.mulFlopsSW_DEBUG += 8;
			Matrix.mulAdopsSW_DEBUG += 11;
			return;
		}
		
		// if we reached the recursion truncation point, multiply current submatrix with standard algorithm
		if (truncP >= dim) {
			Matrix.mulFlopsSW_DEBUG += dim * 2 + dim * dim + dim * dim * dim;
			Matrix.mulAdopsSW_DEBUG += dim * 2 + dim * dim + dim * dim * dim * 2;
			for (int i = 0; i < dim; i++) {
				int ioffsA = offsA + i * dimA, ioffsC = offsC + i * dimC;
				for (int j = 0; j < dim; j++) {
					int joffsB = offsB + j * dimB;
					double v = dA[ioffsA + j];
					if (v < -Matrix.ROUNDOFF_ERROR || v > Matrix.ROUNDOFF_ERROR)
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
		
		Matrix.mulAdopsSW_DEBUG += (23 + sdim2 * 31);
	}
	

	
	
	// polymorphic comparer
	public boolean equals(Matrix B) {
		
		if (B.M != this.M || B.N != this.N) return false;
		
		// compare bitImages directly to mask out inequalities
		if (!this.bitImage.equals(B.bitImage)) return false;
		
		double[] data1 = getDataRef(), data2 = B.getDataRef();
		for (int i = 0; i < bitImage.data.length; i++) {
			if (!BinBitImage.compare(data1, data2, bitImage.data[i])) return false;
		}
		// ordinary cell-by-cell comparison
//		for (int i = 0; i < A.M; i++)
//			for (int j = 0; j < A.N; j++)
//				if (dataA[i * A.N + j] != dataB[i * A.N + j]) return false;
		return true;		
	}

	
	
	// method tests how far the non-zero elements extend from the diagonal, returning a value of the range
	// this allows elimination-type solvers to stop elimination beyond the width of the diagonal data band
	// method runs from left/right edges towards the diagonal in a symmetrical fashion, looking for first non-zero element
	public int diagonality() {
		int width = 0, symOffs = M * M - 1;
		for (int r = 1; r < M; r++) {
			int rN = r * N, rNsym = symOffs - rN;
			for (int i = 0; i < r - width; i++)
				// is element non-zero?
				if (	data[rN + i] > ROUNDOFF_ERROR || data[rN + i] < -ROUNDOFF_ERROR ||
						data[rNsym - i] > ROUNDOFF_ERROR || data[rNsym - i] <-ROUNDOFF_ERROR) {
					// found the first/farthest-out element of this row, break out of loop
					if (r - i > width) width = r - i;
					break;
				}
		}
		return width;
	}
	
	
	
	// factorises/decomposes matrix into two diagonal matrices U and V, useful for the Crout linear solving method:
	//		uuu			v..
	// U:	.uu		v:	vv.
	//		..u			vvv
	public Matrix[] decomposeLU() {
		
		if (M != N)	throw new RuntimeException("Matrix.factorise(): Matrix not square.");
		if (M < 1)	throw new RuntimeException("Matrix.factorise(): Invalid matrix.");

		// all diagonal elements must be nonzero
		for (int d = 0; d < M; d++) {
			double v = data[d * M + d];
			if (v > -ROUNDOFF_ERROR && v < ROUNDOFF_ERROR) return null;
		}
			
		Matrix[] lUV = new Matrix[2];
		lUV[0] = new Matrix("U" + nameCount, M, N, Matrix.Type.Null);
		lUV[1] = new Matrix("V" + nameCount, M, N, Matrix.Type.Null);
		double[] dataU = lUV[0].data, dataV = lUV[1].data;

		for (int c = 0; c < N; c++) {
			
			// calculate u(r,c)
			if (c == 0)	dataU[0] = data[0];
			else
				for (int r = 0; r <= c; r++) {
					int rN = r * N;
					double vuSum = 0;
					for (int k = 0; k < r; k++)
						vuSum += dataV[rN + k] * dataU[k * N + c];
					dataU[rN + c] = data[rN + c] - vuSum;
				}
			
			// calculate v(r,c)
			double uDiag = 1.0 / dataU[c * N + c];
			for (int r = c; r < N; r++) {
				int rN = r * N;
				double vuSum = 0;
				if (c == 0) dataV[rN] = data[rN];
				else {
					for (int k = 0; k < c; k++)
						vuSum += dataV[rN + k] * dataU[k * N + c];
					dataV[rN + c] = (data[rN + c] - vuSum) * uDiag;
				}
			}
		}
		return lUV;
	}
	
	
//	public Matrix[] decomposeLUpivoting() {
//		double[] v = new double[N];
//		
//	}
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			LINEAR SYSTEM SOLVING METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	// applies iterative convergence criterion on matrix and gives true if it's guaranteed to converge
	// criterion is that the sum of quotient of row's elements and diagonal element should be below 1 for every row
	// iterative convergency is stronly influenced by diagonal elements being as large as possible
	public boolean convergent() {

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
		if (DEBUG_LEVEL > 1) System.out.println(converges ? "\n" : "\nconvergence unlikely\n");
		return converges;
	}
	
	
	// conditioning for precision and iterative methods: swap in the largest elements of rows into diagonal positions
	// the constant vector c needs to be swapped as well
	// method returns array "mutator" with the mutated indexes for facilitating reconstitution of the solved vector
	public int[] conditionDiagonal(Matrix c, boolean swapMethod, boolean doBitImage) {
		
		if (M != N)					throw new RuntimeException("Matrix.conditionDiagonal(): Matrix not square.");
		if (M < 1)					throw new RuntimeException("Matrix.conditionDiagonal(): Invalid matrix.");
		if (M != c.M || c.N != 1)	throw new RuntimeException("Matrix.conditionDiagonal(): Invalid constant vector.");
				
		double[] data = this.data, newdata = new double[M*N];
		int[] mutator = null;
		if (swapMethod) {
			mutator = new int[M];
			for (int i = 0; i < M; i++) mutator[i] = i;
		}
		
		for (int r1 = 0, m; r1 < M; r1++) {
			int r1N = r1 * N;
			double a_r1r1 = data[r1N + r1];
			if (a_r1r1 < 0) a_r1r1 = -a_r1r1;
			// find largest (absolute) element in current row (skip checking previous rows for swapMethod)
			if (swapMethod) {
				// apply method that swaps row r1 & row r2 if diagonals of both rows get larger on swapping 
				for (int r2 = 0; r2 < N; r2++) {
					if (r1 != r2) {
					int r2N = r2 * N;

					double a_r2r1 = data[r2N + r1], a_r1r2 = data[r1N + r2], a_r2r2 = data[r2N + r2];
					if ((a_r1r1 < a_r2r1 || a_r1r1 > -a_r2r1) && (a_r1r2 >= a_r2r2 || a_r1r2 <= -a_r2r2)) {
						swap(r1, r2);
						c.swap(r1, r2);					// swap matching elements in vector c
						m = mutator[r1];				// swap indexes in mutator array
						mutator[r1] = mutator[r2];
						mutator[r2] = m;
					}
					}
				}
			} else {
				// apply adding method which adds row r2 to r1 if it provides a larger value for diagonal element of row r1
				int rMaxDiag = 0;
				double vMaxDiag = 0;
				// find the largest diagonal from all other rows than r1
				for (int r2 = 0; r2 < N; r2++) {
					if (r2 != r1) {
						double a_r2r1 = data[r2 * N + r2];
						if (a_r2r1 < 0) a_r2r1 = -a_r2r1;
						if (vMaxDiag < a_r2r1) {
							rMaxDiag = r2;
							vMaxDiag = a_r2r1;
						}
					}
				}
				// add in that row into row r1 only if it had a larger diagonal
				if (vMaxDiag > a_r1r1) {
					int rMax = rMaxDiag * N;
					for (int c1 = 0; c1 < N; c1++) newdata[r1N + c1] = data[r1N + c1] + data[rMax + c1];
				}
			}
		}
		if (!swapMethod) this.data = newdata;
		if (DEBUG_LEVEL > 1) System.out.println("conditionDiagonal():\n" + toString());
		if (doBitImage) bitImage.make();
		return mutator;
	}
	
	
	
	// does Gauss elimination, returns false if matrix was singular
	// taken from http://introcs.cs.princeton.edu/java/95linear/Matrix.java.html
	public boolean doGaussElimination() {

		if (M != N) throw new RuntimeException("Matrix.doGaussElimination(): Matrix not square.");

		// get a copy of matrix data field
		double[] data = getData();

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
			swap(i, max);

			// singular
			if (data[iN + i] <= ROUNDOFF_ERROR && data[iN + i] >= -ROUNDOFF_ERROR) {
				status |= SINGULAR_MATRIX;
				return false;
			}

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
		putDataRef(data);
		if (DEBUG_LEVEL > 1) {
			System.out.println("Gaussian elimination:");
			System.out.println(this.toString());
		}
		bitImage.make();
		return true;
	}

	
	// return x = A^-1 b, assuming A is square and has full rank
	// taken from http://introcs.cs.princeton.edu/java/95linear/Matrix.java.html
	public Matrix solveGaussPartPivoting(Matrix rhs) {
		if (M != N || rhs.M != N || rhs.N != 1)
			throw new RuntimeException("Matrix.solve(): Invalid matrix/vector dimensions.");

		// create copies of the data
		Matrix A = this.clone();
		Matrix b = new Matrix("b", N, 1, rhs.data);
		double[] dataA = A.data, datab = b.data;

		// Gaussian elimination with partial pivoting
		for (int i = 0; i < N; i++) {
			int iN = i * N;

			// find pivot row and swap
			int max = i;
			double vAbs = dataA[iN + i];
			for (int j = i + 1; j < N; j++) {
				double v = dataA[j * N + i];
				if (v > vAbs || v < -vAbs) max = j;
			}
			A.swap(i, max);
			b.swap(i, max);

			// if matrix is singular, return null
			double dAi = dataA[iN + i];
			if (dAi < ROUNDOFF_ERROR && dAi > -ROUNDOFF_ERROR) return null;

			// pivot within b
			dAi = 1.0 / dAi;
			for (int j = i + 1; j < N; j++)
				datab[j] -= datab[i] * dataA[j * N + i] * dAi;

			// pivot within A
			for (int j = i + 1; j < N; j++) {
				double p = dataA[j * M + i] * dAi;
				int jN = j * N;
				for (int k = i + 1; k < N; k++) {
					dataA[jN + k] -= dataA[iN + k] * p;
				}
				dataA[jN + i] = 0.0;
			}
		}

		String newname;
		if (Matrix.DEBUG_LEVEL > 1) newname = "x(" + this.name + "^-1*" + rhs.name + ")";
		else						newname = "x" + nameCount++;
		Matrix x = new Matrix(newname, N, 1, Matrix.Type.Null);
		
		// back substitution
		for (int j = N - 1; j >= 0; j--) {
			double t = 0.0;
			int jM = j * M;
			for (int k = j + 1; k < N; k++)
				t += dataA[jM + k] * x.data[k];
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
	// if duplicate=true, result matrix and inverse matrix returned in Matrix array
	// otherwise this matrix will transform into inverse, and supplied constant vectors into result vectors
	public Matrix[] solveGaussJordanFullPivoting(Matrix B, boolean duplicate) {
		
		if (M != N || B.M != N)	throw new RuntimeException("Matrix.solve(): Invalid matrix/vector dimensions.");

		Matrix A = this, X = B;
		Matrix[] XA = new Matrix[2];
		if (duplicate) {
			X = B.clone();
			A = this.clone();
		}
		XA[0] = X;
		XA[1] = A;
		
		double[] dataX = X.data, dataA = A.data;
		// integer arrays indxc, indxr and ipiv do the bookkeeping on the pivoting
		int[] indxc = new int[N], indxr = new int[N], iPiv = new int[N];
		int iCol = 0, iRow = 0, XN = X.N;
		double det = 1;
		
		// major loop of column reduction
		for(int i = 0; i < N; i++) {
			double vMax = 0;
			for (int j = 0; j < N; j++) {
				int jN = N * j;
				if (iPiv[j] != 1)
					for(int k = 0; k < N; k++) {
						double v = dataA[jN + k];
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
				swap(iRow, iCol);
				X.swap(iRow, iCol);
			}
			indxc[i] = iRow;
			indxr[i] = iCol;
			
			double pivInv = dataA[N * iCol + iCol];
			// signal a singular matrix by returning null
			if (pivInv > -ROUNDOFF_ERROR && pivInv < ROUNDOFF_ERROR) return null;
			// the determinat is built up by multiplying in all the pivot divisors
			det *= pivInv;
			pivInv = 1.0 / pivInv;
			int iColN = iCol * N, iColXN = iCol * X.N;
			// divide pivot row by pivot element
			for (int l = 0; l < N; l++)		dataA[iColN + l] *= pivInv;
			for (int l = 0; l < XN; l++)	dataX[iColXN + l] *= pivInv;

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
					for (int c = 0; c < N; c++)		dataA[rN + c] -= dataA[iColN + c] * v;
					for (int c = 0; c < XN; c++)	dataX[rXN + c] -= dataX[iColXN + c] * v;
				}
			if (DEBUG_LEVEL > 2) System.out.println(A.toString());
		}
		
		// reverse-permutate the result vectors from the column swaps, by swapping columns in reverse order of their buildup
		// (note: the original code seems to contain a mistake here)
		for (int l = N - 1; l >= 0; l--) {
			if (indxr[l] != indxc[l]) {
				for (int k = 0; k < XN; k++) {
					int kNrl = k * N + indxr[l], kNrc = k * N + indxc[l];
					double v = dataX[kNrl];
					dataX[kNrl] = dataX[kNrc];
					dataX[kNrc] = v;
				}
				det = -det;
			}
		}

		if (duplicate)	this.det = det;			// the determinant belongs to the untransformed matrix, which is preserved if duplicate=true
		else			XA[0].det = det;		// if input matrix is transformed, return determinant in solution vectors matrix
		XA[0].name = "X" + nameCount;
		XA[1].name = this.name + "^-1";
		if (DEBUG_LEVEL > 1) {
			System.out.println("Gauss-Jordan full pivoting solver, result vectors:");
			System.out.println(XA[0].toString());
		}
		return XA;
	}
	
	
	// return x = A^-1 b applying Gauss-Jordan method, inverse of A returned in Ai, solution vector x returned by method
	public Matrix solveGaussJordan(Matrix c, Matrix Ai) {

		if (M != N || c.M != N || c.N != 1)
							throw new RuntimeException("Matrix.solve(): Invalid matrix/vector dimensions.");
		if (Ai == null) 	throw new RuntimeException("Matrix.solve(): Invalid inverse reference supplied.");

		Matrix A = this.clone();
		det = 1;
		
		Matrix x = c.clone();
		if (Matrix.DEBUG_LEVEL > 1) x.name = "x(" + name + "^-1*" + c.name + ")";
		else						x.name = "x" + nameCount++;
		
		Ai.data = generateData(M, N, Matrix.Type.Identity);
		Ai.M = Ai.N = M;
		double[] dataA = A.data, dataAi = Ai.data, datax = x.data;
		
		// loop handles every case of subtracting current unitised row from all other rows
		for (int r = 0; r < M; r++) {			

			// find highest value for current diagonal element and swap in it's row
			double vMax = ROUNDOFF_ERROR;
			int rPivot = r;
			for (int r2 = r; r2 < M; r2++) {
				double v = dataA[N * r2 + r];
				if (DEBUG_LEVEL > 2) System.out.println(r+ " " + v + " " + vMax);
				if (v < 0) v = -v;
				if (v > vMax) { vMax = v; rPivot = r2; }
			}
			if (vMax <= ROUNDOFF_ERROR) return null;				// got a singular matrix, abort
			if (r != rPivot) {										// another row had highest value for current diagonal element, swap it in
				A.swap(r, rPivot);
				x.swap(r, rPivot);
				Ai.swap(r, rPivot);
				det = -det;
				if (DEBUG_LEVEL > 2) {
					System.out.println("Swapped in row " + rPivot + " into pivot row " + r);
					System.out.println(A.toString());
				}
			}


			for (int i = 0; i < M; i++) {							// the diagonal element column unitising loop
				int iN = i * N;
				double div = dataA[iN + r];							// get divisor for current lined-up row of A, c, Ai structures
				if (div < -ROUNDOFF_ERROR || div > ROUNDOFF_ERROR) {
					det *= div;
					div = 1.0 / div;
					// go along the lined-up row of A, c, Ai structures, dividing by current element
					for (int j = 0; j < N; j++) {
						dataA[iN + j] *= div;							// divide row in A (which moves towards Identity matrix)
						dataAi[iN + j] *= div;							// divide row in Ai (which moves towards an inverse of A)
					}	
					datax[i] *= div;									// divide input vector c
				} else {
					int rN = r * N;
					det *= dataA[rN + r] + dataA[iN + r];
					div = 1.0 / dataA[rN + r];
					for (int j = 0; j < N; j++) {
						dataA[iN + j] = (dataA[iN + j] + dataA[rN + j]) * div;	// divide row in A (which moves towards Identity matrix)
						dataAi[iN + j] = (dataAi[iN + j] + dataAi[rN + j]) * div;	// divide row in Ai (which moves towards an inverse of A)
					}	
					datax[i] += datax[r] * div;		// divide input vector c
				}
			}
			if (DEBUG_LEVEL > 2) System.out.println(A.toString());

			int rN = r * N;
			// for every row i in the lined-up A, c, Ai structures, subtract row (of diagonal element) r from row i, both above and below r
			for (int i = 0; i < M; i++) {
				int iN = i * N;				
				// don't subtract row r from itself
				if (i != r) {
					for (int j = 0; j < N; j++) {
						dataA[iN + j] -= dataA[rN + j];				// subtract from A
						dataAi[iN + j] -= dataAi[rN + j];			// subtract from Ai
					}
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
			for (int j = 0; j < N; j++)	{
				dataA[rN + j] *= div;
				dataAi[rN + j] *= div;	// divide row r of Ai by a(ii), turning it into the inverse of A
			}
			datax[r] *= div;									// divide c(i) by a(ii), it will turn into solution vector x
		}
		A.det = 1;					// |I| = 1
		Ai.det = 1.0 / det;			// |A^-1| = 1 / |A|

		if (DEBUG_LEVEL > 1) {
			System.out.println("Gauss-Jordan partial pivoting solver with inverse matrix calculation, result vector:");
			System.out.println(x.toString());
			System.out.println("input matrix transformed to identity:");
			System.out.println(A.toString());
			System.out.println("input matrix inverse:");
			System.out.println(Ai.toString());
		}
		return x;
	}
	
	
	
	public Matrix solveCrout(Matrix U, Matrix V) {
		
		Matrix cp = new Matrix("c'", U.M, 1, Matrix.Type.Null);					// c' needs to be generated from input constant vector c
		String newname;
		if (Matrix.DEBUG_LEVEL > 1) newname = "x((VU)^-1*" + this.name + ")";
		else						newname = "x" + nameCount++;
		Matrix x = new Matrix(newname, U.M, 1, Matrix.Type.Null);									// the solution vector x
		double[] dataU = U.data, dataV = V.data, datacp = cp.data, datax = x.data;

		// loop handles every element of c, c', V, U and x
		// c'(r) = (c(r) - sum(v(r,k) * c'(k))) / v(r,r)
		datacp[0] = data[0] / dataV[0];
		for (int r = 1; r < M; r++) {
			int rN = V.N * r;
			double sum = 0;
			for (int k = 0; k < r; k++) sum += dataV[rN + k] * datacp[k];
			datacp[r] = (data[r] - sum) / dataV[rN + r];
		}

		// x(r) = (c'(r) - sum(u(r,k) * x(k))) / u(r,r)
		int r = M - 1;
		datax[r] = datacp[r] / dataU[r * U.N + r--];
		for (; r >= 0; r--) {
			int rN = U.N * r;
			double sum = 0;
			for (int k = r + 1; k < M; k++) sum += dataU[rN + k] * datax[k];
			datax[r] = (datacp[r] - sum) / dataU[rN + r];
		}

		if (DEBUG_LEVEL > 1) {
			System.out.println("Crout solver:");
			System.out.println(x.toString());
		}
		return x;	
	}
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			DETERMINANT/INVERSE MATRIX METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	// helper method that does quick zero determinant heuristics on the matrix/submatrix
	public boolean hasZeroDeterminant(boolean thoroughTest) {
		// Cacc accumulates sums along the columns during the horisontal passes
		double[] data = getDataRef(), Cacc = new double[N], Racc = new double[M];
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
			if (Racc[i] <= ROUNDOFF_ERROR && Racc[i] >= -ROUNDOFF_ERROR) return true;
		}
		
		// check presense of zero columns
		for (int j = 0; j < N; j++) {
			// found zero column, determinant = 0
			if (Cacc[j] <= ROUNDOFF_ERROR && Cacc[j] >= -ROUNDOFF_ERROR) return true;
			double caccj = Cacc[j];
			
			if (thoroughTest) {
				// test finding two columns that are equal
				for (int i = j+1; i < N; i++) {
					if (Math.abs(Cacc[i] - caccj) <= ROUNDOFF_ERROR) {
						// sums of two columns were equal, test individual values
						for (int ct = 0; ct < M;) {
							int ctN = ct * N;
							if (Math.abs(data[ctN + j] - data[ctN + i]) > ROUNDOFF_ERROR)
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
			for (int i = j + 1; i < M; i++)
				if (Math.abs(Racc[i] - raccj) <= ROUNDOFF_ERROR) {
					// sums of two rows were equal, test individual values
					for (int ct = 0; ct < N;) {
						if (Math.abs(data[i * N + ct] - data[j * N + ct]) > ROUNDOFF_ERROR)
							break;
						if (++ct >= N) return true;		// columns were equal, determinant = 0
					}
				}
		}
		return false;
	}
	
	
	// multimode implementation of the Laplace Expansion Recursive algorithm for finding the determinant
	public static double determinantLaplace(Matrix M, int version) {
		
		if (M.M != M.N) throw new RuntimeException("Matrix.determinantLaplace(): Nonsquare matrix.");
		if (M.M < 1 || M.N < 1) throw new RuntimeException("Matrix.determinantLaplace(): Invalid matrix.");
		if (M.M == 1 && M.N == 1) return M.getDataRef()[0];
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
			M.analyseRowsColumns();
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
		double[] data = M.getDataRef();
		// base case: 2x2 matrix solves directly, end of recursion
		if (M.N == 2)
			return data[0] * data[3] - data[1] * data[2];
		
		double det = 0;		// determinant is accumulated here
		// basecase: 3x3 matrix solves directly, end of recursion
		if (M.N == 3) {
			if (data[0] > ROUNDOFF_ERROR || data[0] < -ROUNDOFF_ERROR)
				det += data[0] * (data[4] * data[8] - data[5] * data[7]);
			if (data[1] > ROUNDOFF_ERROR || data[1] < -ROUNDOFF_ERROR)
				det -= data[1] * (data[3] * data[8] - data[5] * data[6]);
			if (data[2] > ROUNDOFF_ERROR || data[2] < -ROUNDOFF_ERROR)
				det += data[2] * (data[3] * data[7] - data[4] * data[6]);
			return det;
		}

		// check if submatrix can return zero determinant immediately, but don't do this test at first entry
		if (doZeroTest && M.hasZeroDeterminant(false)) return 0;

		// for every column
		for (int j = 0; j < M.N; j++) {
			// if multiplicand is not zero (this test makes sparse matrices very fast to calculate)
			if (data[j] > ROUNDOFF_ERROR || data[j] < -ROUNDOFF_ERROR) {
				// recursively call determinantLaplace() on row/column-eliminated submatrices
				det +=  (((j + 2) & 1) == 0 ? 1 : -1) * data[j] * determinantLaplaceR1(M.eliminateRowColumn(0, j, false), true);
			}
		}

		return det;
	}

	
	// quick zero determinant heuristics helper method for version 2 of determinant Laplace Expansion
	// uses the relative increment list activec to jump past recursively eliminated columns
	public boolean hasZeroDeterminant2(int[][] activec, int columns) {

		int topb = M - columns;
		double[] data = getDataRef();
		// Cacc accumulates sums along the columns during the horisontal passes, Racc accumulates sums along the rows
		double[] Cacc = new double[N], Racc = new double[M];
		double v;

		for (int i = topb; i < M; i++) {
			// get row index from detAnalysis if has been created
			int iN = N * (detAnalysis == null ? i : detAnalysis[2][i]);
			for (int j = 0, ac = activec[0][0]; j < columns; j++, ac += activec[0][ac]) {
				v = data[iN + ac - 1];
				Racc[i] += v;
				Cacc[j] += v;
			}
			// found zero row, determinant = 0
			if (Racc[i] <= ROUNDOFF_ERROR && Racc[i] >= -ROUNDOFF_ERROR) return true;
		}
		
		// check presense of zero columns
		for (int j = 0; j < columns; j++) {
			// found zero column, determinant = 0
			if (Cacc[j] <= ROUNDOFF_ERROR && Cacc[j] >= -ROUNDOFF_ERROR) return true;	
		}
		return false;
	}
	
	
	// version two of determinant Laplace Expansion does not use the costly matrix reconstruction method
	// eliminateRowColumn(), but instead readjusts a relative increment list of active columns of each recursive submatrix,
	// activec contains a list of relative increment indexes, telling how much to increment to get to next active column
	public static double determinantLaplaceR2(Matrix M, int[][] activec, int columns) {
		
		Matrix.detL_DEBUG++;
		double[] data = M.getDataRef();
		int tofs1, tofs2, topb = M.M - columns;
		
		// get 2x2 base case row offsets tofs1 & tofs2 from detAnalysis if it has been created
		if (M.detAnalysis != null) {		
			int[] zRIdx = M.detAnalysis[2];
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
			if (v < -ROUNDOFF_ERROR || v > ROUNDOFF_ERROR)
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
	public static double eigenPowerMethod(Matrix A, Matrix y, double thetaError, int iters) {
		
		if (y.M < 2 || y.N != 1) 	throw new RuntimeException("Matrix.eigenPowerMethod(): Invalid guess vector.");
		if (A.M != A.N)				throw new RuntimeException("Matrix.eigenPowerMethod(): Matrix not square.");
		if (A.N != y.M)				throw new RuntimeException("Matrix.eigenPowerMethod(): Matrix-vector dimensionality mismatch.");

		if (DEBUG_LEVEL > 1) System.out.println("Matrix to search: \n" + A.toString());
		
		int debug_lvl = DEBUG_LEVEL, k;
		DEBUG_LEVEL = 0;
		
		Matrix v = y.clone();
		double theta = 0, oldtheta = 0;
		iters = iters < 1 ? 1 : iters;
		
		// do maximum "iters" iterations
		for (k = 0; k < iters; k++) {
	
			v = Matrix.multiply(A, y);
			oldtheta = theta;
			theta = absLength(v);
			if (DEBUG_LEVEL > 1) System.out.println(theta);
			y = Matrix.multiply(1.0 / theta, v);

			if (oldtheta - theta >= -thetaError && oldtheta - theta <= thetaError) break;		// attained eigenvalue precision?
		}
		DEBUG_LEVEL = debug_lvl;
		if (DEBUG_LEVEL > 1) {
			System.out.println("eigenPowerMethod() converged in " + k + " iterations.");
			System.out.println(y.toString());
		}
		return theta;
	}
	
	

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		System.out.println("matrix: " + name);
		if (data != null)
			for (int i = 0; i < M; i++) {
				if (sb.length() > 0) sb.append("\n");
				sb.append("| ");
				for (int j = 0; j < N; j++)
					if (data[i * N + j] > ROUNDOFF_ERROR || data[i * N + j] < -ROUNDOFF_ERROR)
							sb.append(String.format("%5.2f ", data[i * N + j]));
					else	sb.append("      ");
				sb.append("|");
			}
		else sb.append("[null]\n");
		if (M == N) sb.append("\nDeterminant: " + det + "\n");
		else		sb.append("\n");
		return sb.toString();	
	}
	
	
	public enum Type {
		Null(0), Identity(1), Random(2), Centering(3);
		private int type;
		private Type(int type) { this.type = type; }
	}


}

