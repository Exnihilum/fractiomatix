package lkr74.matrixlib;

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
	static final int NULL_MATRIX = 1;
	static final int COMPLEX_MATRIX = 1<<1;
	static final int SINGULAR_MATRIX = 1<<8;
	static final boolean REAL = false;
	static final boolean IMAGINARY = true;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			INSTANCE-LEVEL VALUES
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	protected String name;
	protected int M, N; 						// number of rows & columns
	protected double[] data, idata;				// row-column matrix uses a flat array, data = real part, idata = imaginary part
	protected int status;
	public double det = 0;						// last calculated determinant value

	// optimisation hierarchies & variables
	// binBitLength = no. of values gathered in binBit, procNum = number of concurrent processors available
	protected boolean rowAspectSparsest = true;	// indicates whether row or column aspect is most sparse for determinant analysis
	protected int detSign = 1;					// tells sign of the determinant, according to heuristics of analyseRowsColumns() 
	protected int[][] detAnalysis = null;		// used by analyseRowsColumns() for determinant finding heuristics
	protected int[] mutator = null;				// used to index mutated rows for certain algorithms
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
		this.name = name + nameCount++;
		this.M = r;
		this.N = c;
		setNull();
	}

	//	instantiates Matrix with data of choice
	public Matrix(String name, int r, int c, Type type) {
		// initialise status flags, all bits = 0 means an empty matrix of type Matrix
		if (r < 1 || c < 1) throw new RuntimeException("Matrix(): Illegal matrix dimensions.");
		this.name = name + nameCount++;
		this.M = r;
		this.N = c;
		double[][] dataSet = Matrix.generateData(M, N, type, 1);
		data = dataSet[0];
		idata = dataSet[1];
		if (type == Type.Null || type == Type.Null_Complex) setNull();
		bitImage = new BinBitImage(this);
	}
	
	//	instantiates a matrix with a provided dataset, cloning the dataset into this matrix
	public Matrix(String name, int M, int N, double[] data, double[] idata) {
		if (M < 1 || N < 1 || data.length < M * N) throw new RuntimeException("Matrix(): Illegal matrix dimensions.");
		if (idata != null && idata.length < M * N) throw new RuntimeException("Matrix(): Illegal imaginary array.");
		this.name = name + nameCount++;
		this.M = M;
		this.N = N;
		this.putData(data, idata);
		bitImage = new BinBitImage(this);
		status = 0;
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
		if (mutator != null) A.mutator = mutator.clone();
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
		return A;
	}

	// this clone operation will deconvert from other formats into row-column through overridden getData() method
	public Matrix cloneMatrix() {
		double[][] dataSet = getDataRef();
		Matrix A = new Matrix(name, M, N, dataSet[0], dataSet[1]);
		A.status = status;
		if (bitImage != null) A.bitImage = bitImage.clone(A);
		if (mutator != null) A.mutator = mutator.clone();
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(A.toString());
		return A;
	}
	

	public void zero() {
		if (data != null) for (int i = 0, MN = M * N; i < MN; i++) data[i] = 0;
		if (idata != null) for (int i = 0, MN = M * N; i < MN; i++) idata[i] = 0;
		bitImage.zero();
		if (mutator != null) for (int i = 0; i < M; i++) mutator[i] = 0;
		setNull();
	}

	public Matrix eliminateRowColumn(int r, int c, boolean makeBitImage) {
		
		if (r < 0 || r > M - 1 || c < 0 || c > N - 1)
			throw new RuntimeException("Matrix.eliminateRowColumn(): Row or column out of bounds.");
		if (M < 2 || N < 2)
			throw new RuntimeException("Matrix.eliminateRowColumn(): Invalid matrix size.");

		String newname;
		if (DEBUG_LEVEL > 1) 	newname = new String(name + "(M-" + r + ",N-" + c + ")");
		else 					newname = name;
		Matrix A = null;
		if (idata != null)	A = new Matrix(newname, M - 1, N - 1, Matrix.Type.Null_Complex);
		else				A = new Matrix(newname, M - 1, N - 1, Matrix.Type.Null);
		
		for (int i = 0, ii = 0; i < M; i++) {
			if (i != r) {
				int iiAN = ii * A.N, iN = i * N;
				// real-valued matrix case
				if (idata == null) {
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
		A.bitImage = new BinBitImage(A);
		return A;
	}

	public double valueOf(int r, int c) { return data[r * N + c]; }
	
	public double[] valueOfC(int r, int c) {
		if (data == null || idata == null) throw new RuntimeException("Matrix.valueOf(): Unallocated datafields.");
		double[] iv = new double[2];
		int rNc = r * N + c;
		iv[0] = data[rNc]; iv[1] = idata[rNc];
		return iv;
	}

	public void valueTo(int r, int c, double v) {
		if (data == null) throw new RuntimeException("Matrix.ivalueTo(): Unallocated datafield.");
		data[r * N + c] = v;
		if (!nearZero(v)) {
			bitImage.setBit(r, c);				// set the corresponding bit in bitImage if value is nonzero
			clearNull();
		}
	}
	public void valueToC(int r, int c, double[] iv) {
		if (data == null || idata == null) throw new RuntimeException("Matrix.ivalueTo(): Unallocated datafields.");
		int rNc = r * N + c;
		data[rNc] = iv[0]; idata[rNc] = iv[1];
		if (!nearZero(iv[0]) || !nearZero(iv[1])) {
			bitImage.setBit(r, c);					// set the corresponding bit in bitImage if value is nonzero
			clearNull();
		}
	}

	
	// generate identity matrix with scalar v in the diagonal
	public Matrix identity(int s, double v) {
		Matrix M = new Matrix("I", s, s);
		double[][] dataSet = generateData(s, s, Type.Identity, v);
		M.data = dataSet[0];
		return M;
		}

	
	// polymorphic matrix centering method will numerically center the matrix, columnwise
	public Matrix center(boolean copy) {
		Matrix Ac = this;
		if (copy) Ac = clone();
		double[][] dataSet = Ac.getDataRef();
		double[] dataAc = dataSet[0];
		for (int j = 0; j < N; j++) {
			double avg = 0;
			for (int i = 0, iNj = j; i < M; i++, iNj += N) avg += dataAc[iNj];
			avg /= (double)M;
			for (int i = 0, iNj = j; i < M; i++, iNj += N) dataAc[iNj] -= avg;
		}
		dataAc = dataSet[1];
		if (dataAc != null) {
			double[] idataAc = Ac.idata;
			for (int j = 0; j < N; j++) {
				double avg = 0;
				for (int i = 0, iNj = j; i < M; i++, iNj += N) avg += idataAc[iNj];
				avg /= (double)M;
				for (int i = 0, iNj = j; i < M; i++, iNj += N) idataAc[iNj] -= avg;
			}
		}
		Ac.bitImage.make();
		return Ac;
	}

	
	// generates different types of matrices, a supplied value can be used to alter the matrix modulo
	public static double[][] generateData(int r, int c, Type type, double v) {
		if (r < 1 || c < 1)
			throw new RuntimeException("Matrix.generateData(): Illegal matrix dimensions.");
		double dataSet[][] = new double[2][];
		
		switch (type) {
			case Null_Complex:
				dataSet[1] = new double[r * c];
			case Null:
				dataSet[0] = new double[r * c];
				break;
			case Identity:
				if (r != c) 	throw new RuntimeException("Matrix.generateData(): Identity matrix must be square.");
				dataSet[0] = new double[r * c];
				for (int i = 0; i < c; i++)
					dataSet[0][i * c + i] = v;
				break;
			case Random_Complex:
				dataSet[1] = new double[r * c];
				for (int i = 0; i < r * c; i++) dataSet[1][i] = Math.random();
			case Random:
				dataSet[0] = new double[r * c];
				for (int i = 0; i < r * c; i++) dataSet[0][i] = Math.random();
				break;
			case Centering:
				if (r != c) 	throw new RuntimeException("Matrix.generateData(): Centering matrix must be square.");
				for (int i = 0; i < r; i++)
					for (int j = i * c, jEnd = j + c; j < jEnd; j++) {
						double d = 1.0 / (double)r;
						if (i == j)		dataSet[0][j] = 1.0 - d;
						else			dataSet[0][j] = - d;
					}
				break;
			default:
				throw new RuntimeException("Matrix.generateData(): Illegal matrix type.");
		}
		return dataSet;
	}

	
	public Matrix rescale(int Mi, int Ni, int Mo, int No, boolean doBitImage) {
		
		if (Mi > Mo || Ni > No) throw new RuntimeException("Matrix.rescale(): Invalid data subset size.");
		int newM = Mo - Mi, newN = No - Ni;
		// if rescaled field ends up outside the matrix datafield, set matrix to Null matrix
		if (Mi >= M || Ni >= N || Mo < 0 || No < 0) setNull();
		
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

	
	
	// the Gram-Schmidt orthonormalisation of column vectors constituting a matrix
	// method also finds out if vectors are linearly independent
	// if copy = false, use supplied matrix-form vectors for calculations
	public Matrix orthogonalise(boolean copy) {
		
		Matrix An = this;
		// if copy = false, An will be used to store intermediate q-hat vectors
		if (copy) An = this.clone();
		double r11 = An.absLength(0);										// r11 := ||x1||
		if (nearZero(r11)) return null;										// null vector found, abort

		int M = An.M, N = An.N;
		Matrix Q = new Matrix("Q", M, N, Matrix.Type.Null);
		r11 = 1.0 / r11;
		for (int r = 0, rN = 0; r < M; r++, rN += N)
			Q.data[rN] = An.data[rN] * r11;									// q1 := x1 / r11
		
		// main loop steps column by column
		for (int j = 1; j < N; j++) {
			for (int i = 0; i < j; i++) {
				// use vector j in An as q^, multiply by columns
				double rij = multiplyVectors(An, Q, j, i, false);			// rij := (q^,qi)
				for (int k = 0, kNi = i, kNj = j; k < M; k++, kNi += N, kNj += N)
					An.data[kNj] -= rij * Q.data[kNi];						// q^ := q^ - rij*qi
			}
			double rjj = An.absLength(j);									// rjj := ||q^||
			if (nearZero(rjj)) return null;									// null vector found, abort
			rjj = 1.0 / rjj;
			for (int r = 0, rN = j; r < M; r++, rN += N)
				Q.data[rN] = An.data[rN] * rjj;								// qj := q^ / rjj
		}

		if (DEBUG_LEVEL > 1) {
			System.out.println("Gram-Schmidt orthonormalisation" + (copy ? ": " : "(input matrix wasn't saved): "));
			System.out.println(Q.toString());
		}
		return Q;
	}
	
	
	
	public Matrix decomposeQR(Matrix Q) {
	
		if (M < 1 || N < 1 || !(N==Q.M && M==Q.N))
			throw new RuntimeException("Matrix.decomposeQR(): Invalid matrix dimensions.");

		Matrix R = new Matrix("R", Q.M, N, Matrix.Type.Null);	// matrix D will hold dot products of matrix elements in A & Q
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

	
	
	public Matrix toHouseholder(boolean copy) {
		
		Matrix A = this;
		if (copy) A = this.clone();
		Matrix H = new Matrix("H", M, N, Matrix.Type.Null);
		Matrix I = new Matrix("I", M, N, Matrix.Type.Identity);
		Matrix v = new Matrix("v", M, 1, Matrix.Type.Null), V;
		Matrix vT = new Matrix("vT", 1, N, Matrix.Type.Null);
		
		// proceed columnwise
		for (int j = 0; j < N; j++) {
			double colSum = 0;
			for (int i = N; i < M * N; i+= N)
				colSum += data[i] * data[i];
			double alpha = (data[N] < 0 ? 1 : -1) * Math.sqrt(colSum);
			double r = 0.5 / Math.sqrt(0.5 * alpha * (alpha - data[N]));
			
			v.data[0] = vT.data[0] = 0;
			v.data[1] = vT.data[1] = (data[N] - alpha) * r;
			for (int i = N * 2; i < M * N; i+= N) v.data[i] = vT.data[i] = data[i] * r;
			
			V = v.multiply(vT).multiply(2, false);
			Matrix P = I.subtract(v, true);				// P = I - 2v.v^T
			A = P.multiply(A).multiply(P);				// A(n+1) = PAP
		}
		return A;
	}
	
	

	// finds out number of zeroes in each row & column, sums them up into an analysis array, useful for finding determinant
	public void analyseRowsColumns() {
		
		if (M < 2 || N < 2) {
			if (M < 1 || N < 1) throw new RuntimeException("Matrix.analyseRowsColumns(): Invalid matrix dimensions.");
			return;		// 1x1 matrix, nothing to analyse
		}

		double[][] dataSet = getDataRef();
		double[] data = dataSet[0], idata = dataSet[1];		// get real & imaginary parts of data
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
		
		for (int i = 0; i < M; i++) {
			int iN = N * i;
			if (idata != null) {
				for (int j = 0; j < N; j++)
					if (nearZero(data[iN + j]) || nearZero(idata[iN + j])) { zRCount[i]++; zCCount[j]++; }
			} else {
				for (int j = 0; j < N; j++)
					if (nearZero(data[iN + j])) { zRCount[i]++; zCCount[j]++; }
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
	public void swapRows(int r1, int r2) {
		double temp;
		for (int i = 0, or1 = r1 * N, or2 = r2 * N; i < N; i++, or1++, or2++) {
			temp = data[or1]; data[or1] = data[or2]; data[or2] = temp;
			if (idata != null) { temp = idata[or1]; idata[or1] = idata[or2]; idata[or2] = temp; }
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
			else name = name + "^T"; 
					
			double[] newdata = new double[N * M], inewdata = null;
			if (idata != null) inewdata = new double[N * M];
			
			for (int i = 0; i < M; i++) {
				int iN = i * N;
				for (int j = 0, jMi = i, iNj = iN; j < N; j++, jMi += M, iNj++)
					newdata[jMi] = data[iNj];
				if (idata != null)
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
					if (idata != null) { temp = idataT[jMpi]; idataT[jMpi] = idataT[iNj]; idataT[iNj] = temp; }
					T.bitImage.transposeBit(i, j);
				}
			}		
		}
		if (DEBUG_LEVEL > 1) System.out.println(T.toString());
		return T;
	}
	
	
	// normalises the matrix/vector/transposevector A, for a matrix the columns are normalised individually
	public Matrix normalise(boolean copy) {

		if (M < 1 || N < 1) throw new RuntimeException("Matrix.normalise(): Invalid matrix dimensions.");
		Matrix C;
		double[] dataC = null;
		if (copy) C = this.clone(); else C = this;
		double[][] dataSet = C.getDataRef();
		dataC = dataSet[0];

		if (M > 1)
			// if we're not dealing with a transposed vector
			for (int c = 0; c < N; c++) {
				double norm = 0, sq;
				for (int r = 0; r < M; r++) {
					sq = dataC[r * N + c];
					norm += sq * sq;									// Euclidean norm ||A||2
				}
				norm = 1.0 / Math.sqrt(norm);
				for (int r = 0; r < M; r++) dataC[r * N + c] *= norm;	// normalise
			}
		else {
			double norm = 0, sq;
			for (int c = 0; c < N; c++) {
				sq = dataC[c];
				norm += sq * sq;										// Euclidean norm ||A||2
			}
			norm = Math.sqrt(norm);
			for (int c = 0; c < N; c++) dataC[c] /= norm;		// normalise
		}
		
		if (DEBUG_LEVEL > 1) System.out.println(this.toString());
		return C;
	}

	
	public Matrix add(Matrix B, boolean copy) { return addSub(B, false, copy); }
	public Matrix subtract(Matrix B, boolean copy) { return addSub(B, true, copy); }
	
	// polymorphic addSub will do A+B or A-B returning resultant C of matrix type A, or add/subtract directly to A if clone = false
	// the data from both matrices are gotten through the polymorphic method getDataRef()
	public Matrix addSub(Matrix B, boolean subtract, boolean copy) {
		
		if (B.M != M || B.N != N) throw new RuntimeException("Matrix.addSub(): Nonmatching matrix dimensions.");
		
		// bring in data from input matrices to the superclass format for adding
		double[][] dataSet = getDataRef();
		double[] dataA = dataSet[0], idataA = dataSet[1];
		dataSet = B.getDataRef();
		double[] dataB = dataSet[0], idataB = dataSet[1], dataC, idataC;
		Matrix C;
		if (copy) {
			C = new Matrix("C", M, N, Matrix.Type.Null);
			dataC = C.data; idataC = C.idata;
		} else {
			C = this; 
			dataC = dataA; idataC = idataA;
		}
		
		C.name = (Matrix.DEBUG_LEVEL > 1 ? "(" + name + (subtract?"-":"+") + B.name + ")" : "A" + nameCount++);
		
		if (subtract)	for (int i = 0, MN = M * N; i < MN; i++) dataC[i] = dataA[i] - dataB[i];
		else			for (int i = 0, MN = M * N; i < MN; i++) dataC[i] = dataA[i] + dataB[i];
		// case tree takes care of situation when any of dataA/dataB/dataC can be null array
		if (idataB != null) {
			if (idataA != null) {
				// there is imaginary data in A and B, make sure C has buffer allocated
				if (idataC == null) C.idata = new double[M * N];
				if (subtract)	for (int i = 0, MN = M * N; i < MN; i++) idataC[i] = idataA[i] - idataB[i];
				else			for (int i = 0, MN = M * N; i < MN; i++) idataC[i] = idataA[i] + idataB[i];
			} else idataC = idataB.clone();
		} else
			if (idataA != null && idataC != idataA) idataC = idataA.clone();
		
		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		clearNull();
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
			if (i >= N) 	throw new RuntimeException("Matrix.absLength(): Column index out of bounds.");
			mN = N; vLen = M;						// column vector from a matrix?
		}
		
		if (idata != null)	for (int c = 0, o = i; c < vLen; c++, o += mN) { sq = data[o] + idata[o]; vnorm += sq * sq; }
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
			throw new RuntimeException("Matrix.multiplyVectors(): Invalid vector dimensions.");	
		
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
				if (ia >= aM || ib >= bM) throw new RuntimeException("Matrix.multiplyVectors(): Vector indexes out of bounds.");

				for (int i = bN, ca = ia * aN, cb = ib * aN; i > 0; i--, ca++, cb++)
					p += dataA[ca] * dataB[cb];
			} else {
				if (ia >= aN || ib >= bN) throw new RuntimeException("Matrix.multiplyVectors(): Vector indexes out of bounds.");

				for (int i = bM, ca = ia, cb = ib; i > 0; i--, ca += aN, cb += bN)
					p += dataA[ca] * dataB[cb];
			}
			return p;
		}
		throw new RuntimeException("Matrix.multiplyVectors(): Nonmatching matrix-form vector dimensions.");
	}
	
	
	
	
	// polymorphic matrix product will do AB = C
	public Matrix multiply(Matrix B) {

		int aM = M, aN = N, bM = B.M, bN = B.N;
		if (aM < 1 || aN < 1 || bM < 1 || bN < 1) 	throw new RuntimeException("Matrix.multiply(): Invalid matrix dimensions.");
		if (aN != bM) 								throw new RuntimeException("Matrix.multiply(): Nonmatching matrix dimensions.");
		mulFlops_DEBUG = aM * 2 + aM * aN;			// counts number of multiplications
		
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
//					for (int k = 0; k < B.N; k++)
//						dataC[iCN + k] += v * dataB[jBN + k];
					mulFlops_DEBUG += bN;
				}
			}
		}		

		if (DEBUG_LEVEL > 1) System.out.println(C.toString());
		clearNull();
		C.bitImage = new BinBitImage(C);				// bitImage needs full reconstitution
		return C;
	}
	
	
	
	// multiplies/scales matrix with a value, bitimage will be unchanged (multiplying with zero is meaningless)
	public Matrix multiply(double v, boolean copy) {
		
		if (M < 1 || N < 1) throw new RuntimeException("Matrix.multiply(): Invalid matrix dimensions.");
	
		Matrix C = this;
		double[] dataA = this.getDataRef()[0], dataC = dataA;
		if (copy) {
			String newname = (Matrix.DEBUG_LEVEL > 1 ? "s*" + name : "S" + nameCount++);
			C = new Matrix(newname, M, N, Matrix.Type.Null);
			dataC = C.data;
		}		
		
		for (int i = 0, MN = M * N; i < MN; i++) dataC[i] = v * dataA[i];
		
		C.data = dataC;
		clearNull();
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

		String newname = (DEBUG_LEVEL > 1 ? "(" + A.name + "." + B.name + ")" : "W" + Matrix.nameCount++);
		Matrix C = new Matrix(newname, newM, newM, Type.Null);
		
		buffer2x2StrasWin = new double[8][truncP > 16 ? 16*16 : truncP*truncP];
		
		// subInfo carries along following data:
		// (truncP) the matrix size truncation point when standard multiplicator is used instead
		// (dim) the current submatrix dimension 
		// 3x offsets into the flat datafields, 3x the matrix width of the datafields
		int[] subInfo = {truncP, newM, 0, 0, 0, newM, newM, newM};
		multiplyStrasWin2(A.data, B.data, C.data, subInfo);
		
		C.clearNull();
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
		
		Matrix.mulAdopsSW_DEBUG += (23 + sdim2 * 31);
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
//		for (int i = 0; i < A.M; i++)
//			for (int j = 0; j < A.N; j++)
//				if (dataA[i * A.N + j] != dataB[i * A.N + j]) return false;

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
	public int diagonality() {
		int width = 0, symOffs = M * M - 1;
		for (int r = 1; r < M; r++) {
			int rN = r * N, rNsym = symOffs - rN;
			for (int i = 0, rNi = rN, rNsymi = rNsym; i < r - width; i++, rNi++, rNsymi--) {
				// is element non-zero?
				if (!nearZero(data[rNi]) || !nearZero(data[rNsymi])) {
					// found the first/farthest-out element of this row, break out of loop
					if (r - i > width) width = r - i;
					break;
				}
				if (idata != null) {
					// is imaginary element non-zero?
					if (!nearZero(idata[rNi]) || !nearZero(idata[rNsymi])) {
						// found the first/farthest-out element of this row, break out of loop
						if (r - i > width) width = r - i;
						break;
					}				
				}
			}
		}
		return width;
	}
	
	

	
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
		if (DEBUG_LEVEL > 1) System.out.println("\nConvergence" + (converges ? "" : " not") + " guaranteed\n");
		return converges;
	}
	
	
	
	// does Gauss elimination, returns false if matrix was singular
	// taken from http://introcs.cs.princeton.edu/java/95linear/Matrix.java.html
	public boolean doGaussEliminationPP() {

		if (M != N) throw new RuntimeException("Matrix.doGaussEliminationPP(): Matrix not square.");

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
	public Matrix solveGaussPP(Matrix rhs) {
		if (M != N || rhs.M != N || rhs.N != 1)
			throw new RuntimeException("Matrix.solveGaussPP(): Invalid matrix/vector dimensions.");

		// create copies of the data, bring it into native matrix format
		Matrix A = this.cloneMatrix();
		Matrix b = new Matrix("b", N, 1, rhs.getDataRef()[0], null);
		double[] dataA = A.data, datab = b.data;

		// Gaussian elimination with partial pivoting
		for (int i = 0; i < N; i++) {
			int iN = i * N;

			// find pivot row and swap
			int max = i;
			double vAbs = dataA[iN + i];
			for (int j = i + 1, jNi = i; j < N; j++, jNi += N) {
				double v = dataA[jNi];
				if (v > vAbs || v < -vAbs) max = j;
			}
			A.swapRows(i, max);
			b.swapRows(i, max);

			// if matrix is singular, return null
			double dAi = dataA[iN + i];
			if (nearZero(dAi)) { status |= SINGULAR_MATRIX; return null; }

			// pivot within b
			dAi = 1.0 / dAi;
			for (int j = i + 1, jNi = N * j + i; j < N; j++, jNi += N)
				datab[j] -= datab[i] * dataA[jNi] * dAi;


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
		if (Matrix.DEBUG_LEVEL > 1) newname = "x(" + this.name + "^-1*" + rhs.name + ")";
		else						newname = "x" + nameCount++;
		Matrix x = new Matrix(newname, N, 1, Matrix.Type.Null);
		
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
	public Matrix solveGaussJordanFP(Matrix B, boolean copy) {
		
		if (M != N || B.M != N)	throw new RuntimeException("Matrix.solveGaussJordanFP(): Invalid matrix/vector dimensions.");

		Matrix A = this, X = B;
		if (copy) {
			X = B.clone();
			A = this.clone();
		}
		
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
				X.swapRows(iRow, iCol);
			}
			indxc[i] = iRow;
			indxr[i] = iCol;
			
			double pivInv = dataA[N * iCol + iCol];
			// signal a singular matrix by returning null
			if (nearZero(pivInv)) { status |= SINGULAR_MATRIX; return null; }
			// the determinat is built up by multiplying in all the pivot divisors
			det *= pivInv;
			pivInv = 1.0 / pivInv;
			int iColN = iCol * N, iColXN = iCol * X.N;
			// divide pivot row by pivot element
			for (int l = iColN, lEnd = l + N; l < lEnd; l++)		dataA[l] *= pivInv;
			for (int l = iColXN, lEnd = l + XN; l < lEnd; l++)		dataX[l] *= pivInv;

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
					for (int c1 = rXN, c2 = iColXN, cEnd = c1 + XN; c1 < cEnd;)		dataX[c1++] -= dataX[c2++] * v;
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

		if (copy)	this.det = det;				// the determinant belongs to the untransformed matrix, which is preserved if duplicate=true
		else			X.det = det;			// if input matrix is transformed, return determinant in solution vectors matrix
		X.name = "X" + nameCount++;
		if (DEBUG_LEVEL > 1) {
			System.out.println("Gauss-Jordan full pivoting solver, result vectors:");
			System.out.println(X.toString());
		}
		return X;
	}
	
	
	// return x = A^-1 b with Gauss-Jordan partial pivoting method, inverse of A returned in Ai, returns solution vector x
	// this method checks if matrix is sparse to a diagonal band and changes row processing width to the
	// maximum overlap of the diagonal bands
	public Matrix solveGaussJordanPPDO(Matrix c, Matrix Ai) {

		if (M != N || c.M != N || c.N != 1)
							throw new RuntimeException("Matrix.solveGaussJordanPPDO(): Invalid matrix/vector dimensions.");
		if (Ai == null) 	throw new RuntimeException("Matrix.solveGaussJordanPPDO(): Invalid inverse reference supplied.");

		Matrix A = this.clone();
		det = 1;
		
		Matrix x = c.clone();
		if (Matrix.DEBUG_LEVEL > 1) x.name = "x(" + name + "^-1*" + c.name + ")";
		else						x.name = "x" + nameCount++;
		
		Ai.data = generateData(M, N, Matrix.Type.Identity, 1)[0];
		Ai.M = Ai.N = M;
		double[] dataA = A.data, dataAi = Ai.data, datax = x.data;

		int diag = diagonality();	// find out how far from the diagonal the sparse data stretches at maximum
		int overlap = diag * 2;		// tur it into overlap number, telling how many rows above & below to process for a given row
		
		// loop handles every case of subtracting current unitised row from all other rows
		for (int r = 0; r < M; r++) {			

			// find highest value for current diagonal element and swap in it's row
			double vMax = ROUNDOFF_ERROR;
			int rPivot = r;
			// for a matrix sparse/zeroed outside a diagonal band we'll only process the available row data above & below pivot
			// that is overlapping the pivot datawidth, the degree of overlap was calculated from return of diagonality() method
			for (int r2 = r, overlap1 = overlap; r2 < M && overlap1 > 0; r2++, overlap1--) {
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

			int iStart = (r - overlap < 0 ? 0 : r - overlap);			// find out bounds of sparse datawidth of pivot row
			int iEnd = (r + overlap >= M ? M - 1 : r + overlap);
			
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
		A.det = 1;													// |I| = 1
		Ai.det = 1.0 / det;											// |A^-1| = 1 / |A|

		if (DEBUG_LEVEL > 1) {
			System.out.println("Gauss-Jordan partial pivoting solver with sparse inverse matrix calculation, result vector:");
			System.out.println(x.toString());
			System.out.println("input matrix transformed to identity:");
			System.out.println(A.toString());
			System.out.println("input matrix inverse:");
			System.out.println(Ai.toString());
			System.out.println("Diagonality: " + diag);
		}
		return x;
	}

	
	
	public Matrix solveGaussSeidel(Matrix c, int itersMax, double maxError) {
		
		if (M != N || c.M != N || c.N != 1)
			throw new RuntimeException("Matrix.solveGaussJordan(): Invalid matrix/vector dimensions.");

		DEBUG_LEVEL--;
		
		// create initial solution vector
		Matrix x = new Matrix("x", M, 1, Matrix.Type.Random), xOld = null;
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
							System.out.println("Gauss-Seidel iterator solver won't converge.");
			else			System.out.println("Gauss-Seidel iterator solver done in " + i + " iterations.");
			System.out.println(x.toString());
		}	
		if (converges) return x; else return null;
	}
	
	
	
	// return x = A^-1 b with Gauss-Jordan partial pivoting method, inverse of A returned in Ai, returns solution vector x
	// this method is like solveGaussJordanPP() minus optimisations for diagonal-sparse matrices
	public Matrix solveGaussJordanPP(Matrix c, Matrix Ai) {

		if (M != N || c.M != N || c.N != 1)
							throw new RuntimeException("Matrix.solveGaussJordan(): Invalid matrix/vector dimensions.");
		if (Ai == null) 	throw new RuntimeException("Matrix.solveGaussJordan(): Invalid inverse reference supplied.");

		Matrix A = this.clone();
		det = 1;
		
		Matrix x = c.clone();
		if (Matrix.DEBUG_LEVEL > 1) x.name = "x(" + name + "^-1*" + c.name + ")";
		else						x.name = "x" + nameCount++;
		
		Ai.data = generateData(M, N, Matrix.Type.Identity, 1)[0];
		Ai.M = Ai.N = M;
		double[] dataA = A.data, dataAi = Ai.data, datax = x.data;
		
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
				Ai.swapRows(r, rPivot);
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
					for (int j = 0; j < N; j++) {
						dataA[iN + j] *= div;							// divide row in A (which moves towards Identity matrix)
						dataAi[iN + j] *= div;							// divide row in Ai (which moves towards an inverse of A)
					}	
					datax[i] *= div;									// divide input vector c
				} else {												// divisor was zero, operate on current row with pivot row added to it
					int rN = r * N;
					det *= dataA[rN + r] + dataA[iN + r];
					div = 1.0 / dataA[rN + r];
					for (int j = 0; j < N; j++) {
						dataA[iN + j] = (dataA[iN + j] + dataA[rN + j]) * div;	// divide row in A (which moves towards Identity matrix)
						dataAi[iN + j] = (dataAi[iN + j] + dataAi[rN + j]) * div;	// divide row in Ai (which moves towards an inverse of A)
					}	
					datax[i] = (datax[i] + datax[r]) * div;		// divide input vector c
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
				dataAi[rN + j] *= div;							// divide row r of Ai by a(ii), turning it into the inverse of A
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

	
	
	// factorises/decomposes matrix into two diagonal matrices U and V, useful for the LU backsubstitution linear solving method:
	//		uuu			v..
	// U:	.uu		L:	vv.
	//		..u			vvv
	public Matrix[] decomposeLU() {
		
		if (M != N)	throw new RuntimeException("Matrix.decomposeLU(): Matrix not square.");
		if (M < 1)	throw new RuntimeException("Matrix.decomposeLU(): Invalid matrix.");

		// all diagonal elements must be nonzero
		for (int d = 0; d < M; d++)
			if (nearZero(data[d * M + d])) { status |= SINGULAR_MATRIX; return null; }
			
		Matrix[] lLU = new Matrix[2];
		lLU[0] = new Matrix("L", M, N, Matrix.Type.Null);
		lLU[1] = new Matrix("U", M, N, Matrix.Type.Null);
		double[] dataL = lLU[0].data, dataU = lLU[1].data;

		for (int c1 = 0; c1 < N; c1++) {
			
			// calculate u(r,c)
			if (c1 == 0)	dataU[0] = data[0];
			else
				for (int r = 0; r <= c1; r++) {
					int rN = r * N;
					double vuSum = 0;
					for (int k = 0; k < r; k++)
						vuSum += dataL[rN + k] * dataU[k * N + c1];
					dataU[rN + c1] = data[rN + c1] - vuSum;
				}
			
			// calculate v(r,c)
			double uDiag = 1.0 / dataU[c1 * N + c1];
			for (int r = c1; r < N; r++) {
				int rN = r * N;
				double vuSum = 0;
				if (c1 == 0) dataL[rN] = data[rN];
				else {
					for (int k = 0; k < c1; k++)
						vuSum += dataL[rN + k] * dataU[k * N + c1];
					dataL[rN + c1] = (data[rN + c1] - vuSum) * uDiag;
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

	
	// method takes the factorised L & U triangular matrices applying them on this constant vector to solve
	// a linear system where only the constant vector is changing, the matrix coefficients assumed to be fixed
	// method returns the solution vector x
	public Matrix backSubstituteLU(Matrix L, Matrix U) {
		
		Matrix cp = new Matrix("c'", U.M, 1, Matrix.Type.Null);					// c' needs to be generated from input constant vector c
		String newname;
		if (Matrix.DEBUG_LEVEL > 1) newname = "x((LU)^-1*" + this.name + ")";
		else						newname = "x";
		Matrix x = new Matrix(newname, U.M, 1, Matrix.Type.Null);									// the solution vector x
		double[] dataU = U.data, dataV = L.data, datacp = cp.data, datax = x.data;

		// loop handles every element of c, c', V, U and x
		// c'(r) = (c(r) - sum(v(r,k) * c'(k))) / v(r,r)
		datacp[0] = data[0] / dataV[0];
		for (int r = 1; r < M; r++) {
			int rN = L.N * r;
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
			System.out.println("backSubstituteLU solver:");
			System.out.println(x.toString());
		}
		return x;	
	}
	
	
	// LU decompose with pivoting and recording of permutations, which also calculates determinant
	// adapted from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
	// the decompose is stored directly in current matrix, unless copy = true, then a copy is returned
	Matrix decomposeLU2(boolean copy) {
		
		if (M != N)	throw new RuntimeException("Matrix.decomposeLU2(): Matrix not square.");
		if (M < 1)	throw new RuntimeException("Matrix.decomposeLU2(): Invalid matrix.");
		int[] mutatorLU = null;
		Matrix LU = null;
		double[] dataLU = data;
		if (copy) {
			LU = clone();
			LU.name = "LU" + nameCount++;
			dataLU = LU.data;
			if (mutator == null) mutatorLU = LU.mutator = new int[M];
		} else {
			name = "LU" + nameCount++;			// this matrix will be modified, change it's name
			mutatorLU = mutator = new int[M];
		}
		
		double[] vv = new double[M];
		int d = 1;								// every permutation toggles sign of d
		
		for (int i = 0; i < M; i++) {			// looping over rows getting scaling factor of every row, putting it in vv
			int iN = i * N;
			double vMax = 0;
			for (int j = 0; j < M; j++) {
				double v = dataLU[iN + j];
				if (v < -vMax || v > vMax) vMax = (v < 0 ? -v : v);
			}
			if (nearZero(vMax)) { status |= SINGULAR_MATRIX; return null; }	// if singular matrix found, return null
			vv[i] = 1.0 / vMax;												// save scaling factor
		}
		
		for (int j = 0; j < M; j++) {
			int jN = j * N;
			for (int i = 0; i < j; i++) {
				int iN = i * N;
				double sum = dataLU[iN + j];
				for (int k = 0; k < i; k++) sum -= dataLU[iN + k] * dataLU[k * N + j];
				dataLU[iN + j] = sum;
			}
			
			int iMax = 0;
			double vMax = 0;											// vMax will store the largest pivot element found
			for (int i = j; i < M; i++) {
				int iN = i * N;
				double sum = dataLU[iN + j];
				for (int k = 0; k < j; k++) sum -= dataLU[iN + k] * dataLU[k * N + j];
				dataLU[iN + j] = sum;
				double v = vv[i] * (sum < 0 ? -sum : sum);
				// have we found a better pivot than the best so far?
				if (v >= vMax) { vMax = v; iMax = i; }
			}
			if (j != iMax) {											// test if there's a better pivot and we need to swap rows
				if (copy) LU.swapRows(iMax, j); else swapRows(iMax, j);			// yes
				d = -d;
				vv[iMax] = vv[j];										// also change scale factor
			}
			mutatorLU[j] = iMax;
			// zero at pivot element implies a singular matrix, but some applications can get by with tweaking away the zero
			if (dataLU[jN + j] == 0) {
				status |= SINGULAR_MATRIX;
				dataLU[jN + j] = ROUNDOFF_ERROR;
			}
			if (j < M - 1) {
				double v = 1.0 / dataLU[jN + j];
				for (int i = j + 1; i < N; i++) dataLU[i * N + j] *= v;	// divide column by pivot element
			}
		}
		
		// use d to calculate determinant, put it in LU
		for (int i = 0; i < M; i++) d *= dataLU[i * M + i];
		if (copy) LU.det = d; else this.det = d;
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("decomposeLU2() result:");
			System.out.println(copy ? LU.toString() : this.toString());
			System.out.println(Arrays.toString(mutatorLU));
		}
		if (copy) return LU; else return this;
	}
	
	
	// method takes the factorised L & U triangular matrices using them to solve a linear system
	// where only the constant vector is changing, the matrix coefficients assumed to be fixed
	// adapted from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
	// the triangular matrices L & U come here in a COMBINED form in matrix LU, split by 1:s in the diagonal
	// matrix A, LU & vector b are the input, matrix array with solution vector and LU matrix is the return value
	// if LU is null, then the method calls decomposeLU2() with A, if A is null then LU is used
	// method optimisation takes into account that b might contain many leading zero elements
	public static Matrix[] backSubstituteLU2(Matrix A, Matrix LU, Matrix b, boolean copy) {
		
		if (A == null && LU == null || b == null) throw new RuntimeException("Matrix.backSubstituteLU2(): Null matrices/vector.");
		if (b.N != 1 || b.M != (A == null ? LU.M : A.M))	throw new RuntimeException("Matrix.backSubstituteLU2(): Invalid constant vector.");			
		Matrix[] bLU = new Matrix[2];
		
		// a LU = null tells algorithm to create a new LU decompose, which will be returned in the matrix array
		if (LU == null) {
			if (A.M != A.N)	throw new RuntimeException("Matrix.decomposeLU2(): Matrix A not square.");			
			if (A.M < 1)	throw new RuntimeException("Matrix.decomposeLU2(): Invalid matrix.");
			// copy A into a LU matrix and decompose LU, if copy=false it will destroy matrix A returning it as LU
			bLU[1] = A.decomposeLU2(copy);
			// decompose failed on a singular matrix, return null
			if (bLU[1] == null)	return null;	
		} else
			bLU[1] = LU;
		
		// if copy=true, do not return solution vector in b vector, but in a copy
		if (copy) bLU[0] = b.clone(); else bLU[0] = b;
		if (Matrix.DEBUG_LEVEL > 1) bLU[0].name = "x((" + bLU[1].name + ")^-1*" + b.name + ")";
		else						bLU[0].name = "x" + nameCount++;

		double[] dataLU = bLU[1].data, datab = bLU[0].data, databi = bLU[0].idata;
		
		int ii = -1, N = bLU[1].N;							// when ii > 0 it will become an index of first nonvanishing element of b
		int[] mutatorLU = bLU[1].mutator;
		double sum = 0;
		
		for (int i = 0; i < N; i++) {						
			int iN = i * N, ip = mutatorLU[i];				// the permutations of decomposition stage need to be unpermuted during the passes
			sum = datab[ip];
			datab[ip] = datab[i];
			if (ii >= 0)
				for (int j = ii; j <= i - 1; j++) sum -= dataLU[iN + j] * datab[j];
			else
				if (!nearZero(sum)) ii = i;					// nonzero element found, we'll have to do the sums in loop above from now on
			datab[i] = sum;
		}
		for (int i = N - 1; i >= 0; i--) {					// backsubstitution pass
			int iN = i * N;
			sum = datab[i];
			for (int j = i + 1; j < N; j++) sum -= dataLU[iN + j] * datab[j];
			datab[i] = sum / dataLU[iN + i];				// store an element of solution vector X
		}
		
		// do the imaginary part of vector b if it exists
		if (databi != null) {
			for (int i = 0; i < N; i++) {						
				int iN = i * N, ip = mutatorLU[i];			// the permutations of decomposition stage need to be unpermuted during the passes
				sum = databi[ip];
				databi[ip] = databi[i];
				if (ii >= 0)
					for (int j = ii; j <= i - 1; j++) sum -= dataLU[iN + j] * databi[j];
				else
					if (!nearZero(sum)) ii = i;				// nonzero element found, we'll have to do the sums in loop above from now on
				databi[i] = sum;
			}
			for (int i = N - 1; i >= 0; i--) {				// backsubstitution pass
				int iN = i * N;
				sum = databi[i];
				for (int j = i + 1; j < N; j++) sum -= dataLU[iN + j] * databi[j];
				databi[i] = sum / dataLU[iN + i];			// store an element of solution vector X
			}
		}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("backSubstituteLU2 solver:");
			System.out.println(bLU[0].toString());
			System.out.println(bLU[1].toString());
		}
		return bLU;
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

	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			IMPLEMENTATION OF LAPLACIAN DETERMINAN FINDER (only of academic interest)
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	// multimode implementation of the Laplace Expansion Recursive algorithm for finding the determinant
	public static double determinantLaplace(Matrix M, int version) {
		
		if (M.M != M.N) throw new RuntimeException("Matrix.determinantLaplace(): Nonsquare matrix.");
		if (M.M < 1 || M.N < 1) throw new RuntimeException("Matrix.determinantLaplace(): Invalid matrix.");
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
			int iN = N * (detAnalysis == null ? i : detAnalysis[2][i]);
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
	
	
	// version two of determinant Laplace Expansion does not use the costly matrix reconstruction method
	// eliminateRowColumn(), but instead readjusts a relative increment list of active columns of each recursive submatrix,
	// activec contains a list of relative increment indexes, telling how much to increment to get to next active column
	public static double determinantLaplaceR2(Matrix M, int[][] activec, int columns) {
		
		Matrix.detL_DEBUG++;
		double[] data = M.getDataRef()[0];
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
		
		if (y.M < 2 || y.N != 1) 	throw new RuntimeException("Matrix.eigenPowerValue(): Invalid guess vector.");
		if (M != N)					throw new RuntimeException("Matrix.eigenPowerValue(): Matrix not square.");
		if (N != y.M)				throw new RuntimeException("Matrix.eigenPowerValue(): Matrix-vector dimensionality mismatch.");

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

			if (oldtheta - theta >= -thetaError && oldtheta - theta <= thetaError) break;		// attained eigenvalue precision?
		}
		DEBUG_LEVEL = debug_lvl;
		if (DEBUG_LEVEL > 1) {
			System.out.println("eigenPowerMethod() converged in " + k + " iterations.");
			System.out.println(y.toString());
		}
		return theta;
	}
	
	
	
	// findEigenQR() uses the common A(k+1) = Q*AQ iteration routine
	public Matrix findEigenQR(Matrix Q, int iters, int vIters, double maxErr) {

		DEBUG_LEVEL--;
		Matrix Alam = this.clone(), QT = Q.transpose(true);
		// the null matrix for eigenvector solving
		Matrix z = new Matrix("z", M, 1, Matrix.Type.Null);
		Matrix[] v = new Matrix[M];
		double currErr = 0, lastErr = 0;
		double[] dataAlam = Alam.data, lastD = new double[M];
		maxErr *= maxErr;					// compare squared values to avoid unnecessary square root
		
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
		DEBUG_LEVEL++;
		
		// time to find the eigenvectors through v(i) = (lambda(i)*I - A)^-1 * zero-vector
		for(int j = 0; j < M; j++) {
			Matrix Ilam = identity(M, dataAlam[j * M + j]);
			Matrix AmIlam = this.subtract(Ilam, true);
			//v[j] = (this.subtract(Ilam, true)).solveGaussSeidel(z, vIters, maxErr);
			v[j] = AmIlam.solveGaussJordanFP(z, true);
		}

		if (DEBUG_LEVEL > 1) {
			System.out.println("QR eigenfinder took " + i + " iterations, results:");
			for(i = 0; i < M; i++) {
				System.out.println("eigenvalue" + (i+1) + ": " + Alam.data[i * M + i]);
				//System.out.println("eigenvector" + (i+1) + (v[i] != null ? ":" : " failed converging."));
				//if (v[i] != null) System.out.println(v[i].toString());
			}
			System.out.println(Alam.toString());
			System.out.println(Q.toString());
		}
		return Alam;
	}

	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			INLINE/SHORTFORM METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	public static boolean nearZero(double v) { return (v < -ROUNDOFF_ERROR || v > ROUNDOFF_ERROR ? false : true); }
	boolean isNull() { return (status & NULL_MATRIX) != 0; }
	void clearNull() { status ^= NULL_MATRIX & status; }
	void setNull() { status |= NULL_MATRIX; }
	
	// get complex modulus of imaginary value
	double complexM(double r, double i) { return Math.sqrt(r * r + i * i); };
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// limits output of a double to max 5 characters, for following cases:
	// -XXeE		v >= 10000
	// -XXXX		v >= 100
	// -XX.D		v >= 10
	// -X.DD		v < 10
	private static String to5chars(double v, boolean imaginary) {
		int decimals = 0;
		if (v >= 10000 || v <= -10000) {
			for (; v >= 100 || v <= -100; v /= 10.0) decimals++;
			if (imaginary)
					return (v < 0 ? " " : "  ") + (int)v + "e" + decimals + "  ";
			else	return (v < 0 ? " " : "  ") + (int)v + "e" + decimals + (imaginary ? "i " : " ");
		}
		
		boolean isInt = nearZero(v - (int)v);
			
		if (v >= 100 || v <= -100) {
			return imaginary ? String.format("%6di", (int)v) : String.format("%6d ", (int)v);
			//return imaginary ? String.format("%6.0fi", v) : String.format("%6.0f ", v);
		}
			
		if (v >= 10 || v <= -10) {
			if (isInt) return imaginary ? String.format("%5di ", (int)v) : String.format("%5d  ", (int)v);
			return imaginary ? String.format("%6.1fi", v) : String.format("%6.1f ", v);
		}
			
		if (isInt) return imaginary ? String.format("%4di  ", (int)v) : String.format("%4d   ", (int)v);
		return imaginary ? String.format("%6.2fi", v) : String.format("%6.2f ", v);
	}
	
	private static int MAX_PRINTEXTENT = 20;
	
	@Override
	public String toString() {
		
		double[][] dataSet = this.getDataRef();
		double[] data = dataSet[0], idata = dataSet[1];
		
		StringBuilder sb = new StringBuilder();
		int maxM = M > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : M;
		int maxN = N > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : N;
		
		System.out.println(((M == 1 || N == 1) ? "vector: " :"matrix: ") + name + (M == 1 ? "^T" : ""));
		if (data != null) {
			for (int i = 0; i < maxM; i++) {
				int iN = i * N;
				if (sb.length() > 0) sb.append("\n");		// don't add newline at start of matrix printout
				sb.append("|");
				for (int j = 0; j < maxN; j++)
					if (nearZero(data[iN + j])) 	sb.append("   -   ");
					else							sb.append(to5chars(data[iN + j], false));

				// if matrix was bigger than allowed printout bounds, indicate the continuation
				if (maxN < N) if (i % 4 == 0) sb.append(" ..."); else sb.append("    ");
				
				// any imaginary data comes as a second line under the real data line
				if (idata != null) {
					sb.append(" |\n|");
					for (int j = 0; j < maxN; j++) sb.append("       ");
					sb.append("     |\n|");
					for (int j = 0; j < maxN; j++)
						if (nearZero(data[iN + j]))	sb.append("       ");
						else						sb.append(to5chars(data[iN + j], true));
					
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
		if (M == N) sb.append("\nDeterminant: " + det + "\n");
		else		sb.append("\n");
		return sb.toString();	
	}
	
	
	public enum Type {
		Null(0), Identity(1), Random(2), Centering(3), Null_Complex(16), Random_Complex(17);
		private int type;
		private Type(int type) { this.type = type; }
	}
}

