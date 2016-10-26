# fractiomatix
Mathematical code repository for matrices and specific combinatorial problems of interest.

Note: this is in no way a finished, streamlined library, but rather a collection of working algorithms with a surrounding structure still needing solid refactoring massage and organisation.

Included are methods for processing FLAT-ARRAY matrices (note: not the ordinary 2D Java arrays), both ordinary row/column, CSR matrices and my sparse matrix implementation called NSPMatrix.

Polymorphic methods getData(), getDataRef(), putData(), putDataRef() allow bringing in data from both ordinary and CSR packing into any algorithm (CSRMatrix class extends Matrix class), the basis data order will naturally be the row/column one.

For CSR matrices, following methods implemented: add/subtract (both with matrix and value), multiply, swap rows, transpose, value insertion, row&column elimination, equal, toString, and conversion to/from CSR ofcourse.

For NSP matrices, adding, multiplication, loading from MatrixMarket files, row/column swapping and sparse LU factorisarion/solving methods are implemented.

Included are a Gaussian elimination with partial pivoting, 2 versions of Gauss-Jordan solver (one with partial pivoting, trivial sparse matrix optimisation and inverse matrix generation, one with full pivoting and destruction of original matrix), a static coefficient solver with LU decomposition.
QR and LU decompositions implemented.
A finder of largest eigenvalue & vector included (the simple Power Method).
A finder of all eigenvalues and vectors using orthogonalisation and Householder reduction implemented. Probably not the fastest on the planet, but we all have to start somewhere.

I've implemented in Java the O(n^2.783) Strassen-Winograd matrix multiplicator (which coincidentally would have to be taking on quite big matrices to be worth the bother, as it is changing stack levels and array context all the time. I've included a little speed test in the main() method). Eventually I will change the operations to a preallocated fixed array for all array operations, but I doubt it will gain any major speed improvement, as I suspect the allocations in Java happen in nearly O(1) time when a lot of free memory is available, as I have in these test cases. Perhaps getting rid of Java's incessant array zeroing will give a little gain. So it's very low priority indeed.

I am currently more interested in sparse-matrix fast-solving methods for physical applications. I'd like to be able to import 3D meshes into Java from ex. FBX format, decompose into tetrahedral volumes and run electromagnetic field simulations and stress/strain/displacement simulations with support for (hopefully) anisotropic materials.
Therefore I have begun implementing the FrontalDAG symbolic factorisation framework class, which utilises Anshul Gupta's theorems for sparse frontal matrix coefficient inclusion, which he implements in the WSMP LU decomposer. A matrix is loaded from MatrixMarket file format into NSP format, then converted into a task-DAG and a data-DAG. Then it is passed to the FrontalMatrix class for postorder traversal and frontal matrix construction, absorption and decomposition.

There is a dependency on [jfreechart] (https://github.com/jfree/jfreechart) for drawing routines.

Current contents (outdated, as FrontalMatrix & FrontalDAG classes have appeared):

		lkr74.mathgenerics:
  			class MiscMath
    				public class HatBasis
      					HatBasis(double span)
      					public double get(double x)
    				public static class RandFill
    					RandFill(int range)
      					public int remainingSlots()
      					public int getRandom()
    				public static double[] solveCubicPolynomial(double A, double B, double C)
    				public static double[] findLinearF(double[] xdata, double[] ydata)
    				public static double[] findLinearPartialF(double[] xdata, double[] ydata, int plength)
    				public static double getInterpolationOfLinPartials(double[] xb0b1array, double x)
    				public static double triangle3DArea(double[] a, double[] b, double[] c, boolean root)
  				 	public static double triangle3DArea2(double[] r1, double[] r2, double[] r3, boolean root)
    				public static class ElevatorMovement
      					public ElevatorMovement(double d, double v, double a1, double a2)
      					public double getAscendingDistance(double t)
      					public double getDescendingDistance(double t)
      					public double getTime(double d)
  
		lkr74.matrixlib:
  			public class Matrix
    				public Matrix(String name, int r, int c)
    				public Matrix(String name, int r, int c, Type type)
    				public Matrix(String name, int M, int N, double[] data, double[] idata)
    				public static void initExecutor(int processors)
    				public static void initTaskList(int processors)
    				public static void finishTaskList() 
    				public void makeThreaded()
    				public void makeNonThreaded()
    				public double[][] getDataRef()
    				public double[][] getData()
    				public void putDataRef(double[] data, double[] idata)
    				public void putData(double[] data, double[] idata)
    				public Matrix clone()
    				public Matrix cloneMatrix()
    				public void zero()
    				public Matrix eliminateRowColumn(int r, int c, boolean makeBitImage)
    				public double valueOf(int r, int c)
    				public double[] valueOfC(int r, int c)
    				public void valueTo(int r, int c, double v)
    				public void valueToC(int r, int c, double[] iv)
    				public Matrix identity(int s, double v)
    				private double norm1()
    				public Matrix center(boolean copy)
    				public static double[][] generateData(int r, int c, Type type, double v)
    				public Matrix rescale(int Mi, int Ni, int Mo, int No, boolean doBitImage)
    				public Matrix orthogonalise(boolean copy)
    				public Matrix decomposeQR(Matrix Q)
    				public Matrix reduceHouseholder(boolean copy)
    				public void analyseRowsColumns()
    				public void swapRows(int r1, int r2)
    				public Matrix transpose(boolean copy)
    				public Matrix normalise(boolean copy)
    				public Matrix add(Matrix B, boolean copy)
    				public Matrix subtract(Matrix B, boolean copy)
    				public Matrix addSub(Matrix B, boolean subtract, boolean copy)
    				public double absLength(int i)
    				public static double multiplyVectors(Matrix A, Matrix B, int ia, int ib, boolean byRows)
    				public Matrix multiply(Matrix B)
    				public Matrix multiply(double v, boolean copy)
    				public static Matrix multiplyStrasWin(Matrix A, Matrix B, int truncP, boolean useDBversion)
    				public static void multiplyStrasWin2(double[] dA, double[] dB, double[] dC, int[] subInfo)
    				public boolean equals(Matrix B)
    				public int halfBandwidth()
    				Matrix[] decomposeLU(boolean copy, boolean generateLandU)
    				public void unpermute(int[] mutator)
    				public boolean convergent()
    				public Matrix solveGaussPP(Matrix rhs)
    				public boolean doGaussEliminationPP()
    				public Matrix solveGaussJordanFP(Matrix B, boolean copy)
    				public Matrix solveGaussJordanPP(Matrix c, Matrix Ai)
    				public Matrix solveGaussJordanPPDO(Matrix c, Matrix Ai)
    				public Matrix solveGaussSeidel(Matrix c, int itersMax, double maxError)
    				public Matrix[] backSubstituteLU(Matrix A, Matrix LU, boolean copy)
    				public boolean hasZeroDeterminant(boolean thoroughTest)
    				public static double determinantLaplace(Matrix M, int version)
    				public static double determinantLaplaceR1(Matrix M, boolean doZeroTest)
    				public boolean hasZeroDeterminant2(int[][] activec, int columns)
    				public static double determinantLaplaceR2(Matrix M, int[][] activec, int columns)
    				public double eigenPowerValue(Matrix y, double thetaError, int iters)
    				public Matrix findEigenQR(Matrix Q, int iters, int vIters, double maxErr)
    				private static String to5chars(double v, boolean complex)
    				public String toString()

  			public class CSRMatrix extends Matrix
    				public CSRMatrix(String name, int M, int N)
    				public CSRMatrix(String name, int M, int N, double[] data, double[] idata)
    				public CSRMatrix(String name, int M, int N, Type type)
    				public Matrix clone()
    				public static CSRMatrix convert(Matrix A)
    				public void zero()
    				public double[][] getDataRef()
    				public double[][] getData()
    				public void putDataRef(double[] data, double[] idata)
    				public void putData(double[] data, double[] idata)
    				public void putDataM(double[] data, double[] idata)
    				public CSRMatrix rescale(int Mi, int Ni, int Mo, int No, boolean doBitImage)
    				public double valueOf(int r, int c)
    				public void valueTo(int r, int c, double v)
    				public CSRMatrix identity(int S, double v)
    				public CSRMatrix eliminateRowColumn(int r, int c, boolean makeBitImage)
    				public Matrix transpose(boolean copy)
    				public static boolean equal(CSRMatrix S, CSRMatrix T)
    				public void swap(int r1, int r2)
    				public CSRMatrix add(double v, boolean copy)
    				public CSRMatrix subtract(double v, boolean copy)
    				public CSRMatrix add(CSRMatrix T, boolean copy)
    				public CSRMatrix subtract(CSRMatrix T, boolean copy)
    				public CSRMatrix addSub(CSRMatrix T, boolean subtract, boolean copy, double scl, boolean addScalar)
    				public CSRMatrix multiply(CSRMatrix T)
    				public String toString()
    
  			public class BinBitImage
    				public BinBitImage()
    				public BinBitImage(Matrix M)
    				public void make()
    				public BinBitImage clone(Matrix M)
    				public void zero()
    				void or(BinBitImage bI)
    				void and(BinBitImage bI)
    				void setBit(int r, int c) 
    				void clearBit(int r, int c)
    				void transposeBit(int r, int c) 
    				void swapRows(int r1, int r2)
    				protected boolean equals(BinBitImage bI)
    				public static boolean compare(double[] c, double[] d, long v)
    				public static int compact(double[] d, long v, int id)
    				public String toString()
    				protected static String binBitToString()
    
  			public class MatrixApp
    				static XYDataset createStatisticSet(double[][] timingLists, String[] testNames, int[] testCases)
    				private static void matrixMultiTest(int[] testCases)
    				public static void main(String[] args)
