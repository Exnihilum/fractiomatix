# fractiomatix
Mathematical code repository for matrices and specific combinatorial problems of interest.<br>

Current contents:

		lkr74.mathgenerics:
  			class MiscMath<br>
    				public class HatBasis
      					HatBasis(double span)
      					public double get(double x)
    				public static class RandFill
    						RandFill(int range)<br>
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
  
		lkr74.matrixlib
  			public class Matrix implements Cloneable
    				public Matrix(String name, int r, int c)
    				public Matrix(String name, int r, int c, Type type)
    				public Matrix(String name, int M, int N, double[] data)
    				public static void initExecutor(int processors)
    				public static void initTaskList(int processors)
    				public static void finishTaskList() 
    				public void makeThreaded()
    				public void makeNonThreaded()
    				public double[] getDataRef()
    				public double[] getData()
    				public void putDataRef(double[] data)
    				public void putData(double[] data) 
    				public Matrix clone()
    				public void zero()
    				public Matrix eliminateRowColumn(int r, int c, boolean makeBitImage)
    				public double valueOf(int r, int c) 
    				public Matrix identity(int s)
    				public static Matrix center(Matrix A)
    				public static double[] generateData(int r, int c, Type type)
    				public Matrix rescale(int Mi, int Ni, int Mo, int No, boolean doBitImage)
    				public void analyseRowsColumns()
    				public void swap(int r1, int r2)
    				public Matrix transpose(boolean copy)
    				public static void normalise(Matrix A)
    				public static Matrix add(Matrix A, Matrix B)
    				public static Matrix subtract(Matrix A, Matrix B)
    				public static double absLength(Matrix v)
    				public static double vmultiply(Matrix a, Matrix b)
    				public static Matrix multiply(Matrix A, Matrix B)
    				public static Matrix multiply(double v, Matrix A)
    				public static Matrix multiplyStrasWin(Matrix A, Matrix B, int truncP, boolean useDBversion)
    				public static void multiplyStrasWin2(double[] dA, double[] dB, double[] dC, int[] subInfo)
    				public boolean equals(Matrix B)
    				public int diagonality()
    				public Matrix[] factorise()
    				public boolean convergent()
    				public int[] conditionDiagonal(Matrix c, boolean doBitImage)
    				public void doGaussElimination()
    				public Matrix solve(Matrix rhs)
    				public Matrix solveGaussJordan(Matrix c, Matrix Ai)
    				public Matrix solveCrout(Matrix U, Matrix V)
    				public boolean hasZeroDeterminant(boolean thoroughTest)
    				public static double determinantLaplace(Matrix M, int version)
    				public static double determinantLaplaceR1(Matrix M, boolean doZeroTest)
    				public boolean hasZeroDeterminant2(int[][] activec, int columns)
    				public static double determinantLaplaceR2(Matrix M, int[][] activec, int columns)
    				public static double eigenPowerMethod(Matrix A, Matrix y, double thetaError, int iters)
    				public String toString()

  			public class CSRMatrix extends Matrix
    				public CSRMatrix(String name, int M, int N)
    				public CSRMatrix(String name, int M, int N, double[] data)
    				public CSRMatrix(String name, int M, int N, Type type)
    				public Matrix clone()
    				public void zero()
    				public double[] getDataRef()
    				public double[] getData()
    				public void putDataRef(double[] data)
    				public void putData(double[] data)
    				public CSRMatrix rescale(int Mi, int Ni, int Mo, int No, boolean doBitImage)
    				public double valueOf(int r, int c)
    				public void valueTo(int r, int c, double v)
    				public CSRMatrix identity(int S)
    				public CSRMatrix eliminateRowColumn(int r, int c, boolean makeBitImage)
    				public Matrix transpose(boolean copy)
    				public static boolean equal(CSRMatrix S, CSRMatrix T)
    				public void swap(int r1, int r2)
    				public static CSRMatrix add(CSRMatrix S, CSRMatrix T)	
  					public static CSRMatrix subtract(CSRMatrix S, CSRMatrix T, boolean subtract)
  					public static CSRMatrix addSub(CSRMatrix S, CSRMatrix T, boolean subtract)
    				public static CSRMatrix multiply(CSRMatrix S, CSRMatrix T)
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
