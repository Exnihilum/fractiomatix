package lkr74.matrixlib;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.InvalidParameterException;
import java.util.Arrays;

import lkr74.sizeof.SizeOf;


public class NSPMatrix extends Matrix {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			DYNAMICALLY RESIZED SPARSE NODE MATRIX AND IT'S OPERATIONS									//
	//			Leonard Krylov 2016																			//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Hsp & Vsp are sparse row/column arrays, preallocated to CSR2_ALLOCBLOCK positions each for fast fill-in
	// at end of list, each array is reallocated on hitting CSR2_ALLOCBLOCK-element bounds with another CSR2_ALLOCBLOCK
	public NSPArray[] Hsp;				// horisontally/row aligned sparse arrays pointing to nodes (resembling JA in typical CSR)
	public NSPArray[] Vsp;				// vertically/column aligned sparse arrays
	public NSPNode[] pivotNsp;			// fast-access to the diagonal nodes
	
	// define how many nodes to allocate on every forced expansion of a Hsp/Vsp array, and how many to deallocate if count falls below so many nodes
	private static final int CSR2_ALLOCBLOCK = 16;
	private static final int CSR2_DEALLOCBLOCK = CSR2_ALLOCBLOCK * 4;
	
	
	private void initNSPMatrix() {						// initialises all NSPMatrix arrays from scratch
		nNZ = 0;
		Hsp = new NSPArray[M];
		for (int i = 0; i < Hsp.length; i++) Hsp[i] = new NSPArray(0, 0, null);
		Vsp = new NSPArray[N];
		for (int i = 0; i < Vsp.length; i++) Vsp[i] = new NSPArray(0, 0, null);
		pivotNsp = new NSPNode[M];
		mutator = new int[2][];
		mutator[0] = new int[M]; mutator[1] = new int[N];
		for (int i = 0; i < M; i++) mutator[0][i] = i;
		for (int j = 0; j < N; j++) mutator[1][j] = j;
		setNull();
	}
	

	public NSPMatrix(String name, int M, int N) {
		super(name, M, N);								// get skeleton Matrix superclass instance
		initNSPMatrix();
	}

	
	public NSPMatrix(String name, int M, int N, double[] data, double[] idata) {
		super(name, M, N);								// get skeleton Matrix superclass instance
		putData(data, idata);							// NSP dual-aspect dynamic list creation
		if (data != null) clearNull();
		if (idata != null) setComplex();
		//bitImage = new BinBitImage(this);
		if (Matrix.DEBUG_LEVEL > 2)  System.out.println(this.toString());		
	}

	public NSPMatrix(String name, int M, int N, Type type, double v) {
		super(name, M, N);								// get skeleton Matrix superclass instance
		initNSPMatrix();
		if (type != Matrix.Type.Null && type != Matrix.Type.Null_Complex) {
			this.generateData(type, v);
			data = idata = null;
			clearNull();
		}
		//bitImage = new BinBitImage(this);
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
	}

	
	
	@Override
	public NSPMatrix clone() {
		NSPMatrix A = (NSPMatrix)super.clone();					// will first call Matrix superclass clone() method

		if (Vsp != null) {										// duplicate Vsp reference arrays first
			A.Vsp = new NSPArray[Vsp.length];					// as we will put references in them when iterating Hsp
			int i = 0;
			for (NSPArray aVsp : Vsp) {
				A.Vsp[i] = aVsp.clone();
				if (aVsp.array != null) A.Vsp[i].array = new NSPNode[A.Vsp[i].size = aVsp.nodes];
				i++;
			}
		}
		A.pivotNsp = new NSPNode[M];							// create pivot fast-access array

		if (Hsp != null) {
			A.Hsp = new NSPArray[Hsp.length];
			int i = 0;
			for (NSPArray aHsp : Hsp) {
				A.Hsp[i] = aHsp.clone();
				NSPArray aHspA = A.Hsp[i];
				if (aHsp.array != null) {
					aHspA.array = new NSPNode[aHspA.size = aHsp.nodes];
					for (int j = 0; j < aHsp.nodes; j++) {
						NSPNode node = aHsp.array[j].clone();
						aHspA.array[j] = node;
						if (node.r == node.c) A.pivotNsp[node.r] = node;
						A.Vsp[node.c].array[node.offV] = node;
					}
				}
				i++;
			}
		}
		if (mutator != null) A.mutator = mutator.clone();
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
		return A;
	}

	
	
	@Override
	public double valueOf(int r, int c) {
		
		if (r < 0 || c < 0 || r >= M || c >= N)
			throw new InvalidParameterException("NSPMatrix.valueOf(): Invalid matrix coordinates.");
		if (Hsp[r] == null) return 0;
		// search in the aspect with fewest contained nodes
		if (Hsp[r].nodes < Vsp[c].nodes) {
			int offset = findHVspNode(Hsp[r].array, 0, Hsp[r].nodes - 1, -1, c);
			if (offset >= 0) return Hsp[r].array[offset].v;					// return the value
		} else {
			int offset = findHVspNode(Vsp[c].array, 0, Vsp[c].nodes - 1, r, -1);
			if (offset >= 0) return Vsp[c].array[offset].v;					// return the value
		}
		return 0;
	}
	
	// looks up a value in a Mx1 vector
	public double valueOf(int r) {
		
		if (r < 0 || r >= M) throw new InvalidParameterException("NSPMatrix.valueOf(): Invalid matrix coordinates.");
		NSPArray aVsp = Vsp[0];
		if (aVsp == null) return 0;
		int offset = findHVspNode(aVsp.array, 0, aVsp.nodes - 1, r, -1);
		if (offset >= 0) return aVsp.array[offset].v;					// return the value
		return 0;
	}

	
	
	@Override
	public void valueTo(int r, int c, double v) {
		if (r < 0 || c < 0 || r >= M || c >= N)
			throw new InvalidParameterException("NSPMatrix.valueTo(): Invalid matrix coordinates.");
		
		if (nearZero(v)) {												// insertion of zero means removal of a value
			int nNodes = Hsp[r].nodes;
			NSPNode[] bHsp = Hsp[r].array;
			int offH = findHVspNode(bHsp, 0, nNodes - 1, -1, c);		// locate offset of node in Hsp
			int offV;
			if (offH >= 0) offV = bHsp[offH].offV;						// get node's offset into Vsp if node existed in Hsp
			else return;											// a node that doesn't exist in Hsp will not exist in corresponding Vsp
			removeLocalHVspNode(Hsp[r], offH, 0);						// dereference node from horisontal sparse array
			removeLocalHVspNode(Vsp[c], offV, 1); nNZ--;				// dereference node from vertical sparse array
			if (r == c) pivotNsp[r] = null;
			return;
		}
		NSPNode node = new NSPNode(r, c, v, 0); nNZ++;
		if (r == c) pivotNsp[r] = node;
		int offset = 		insertHVspNode(r, c, 0, node);				// insert value as new node into Hsp
		if (offset >= 0) 	insertHVspNode(r, c, 1, node);				// if node didn't exist already, insert value as new node into Vsp
	}

	// optimised value insertion method for FrontalDAG library
	public void valueTo2(int p, double v) {
		NSPNode node = new NSPNode(p, p, v, 0); nNZ++;
		pivotNsp[p] = node;
		insertHVspNode(p, p, 0, node);
		insertHVspNode(p, p, 1, node);
	}

		
	
	public NSPMatrix transpose(boolean copy) {
		if (M < 1 || N < 1) throw new InvalidParameterException("NSPMatrix.transpose(): Invalid matrix dimensions.");

		NSPMatrix T = copy ? this.clone() : this;
		T.name = name + "^T"; 

		for (int i = 0; i < M; i++) {
			NSPNode[] bHsp = T.Hsp[i].array;
			int nNodes = T.Hsp[i].nodes;
			// need to swap row & column indexes in every node
			for (int offH = 0; offH < nNodes; offH++) {
				NSPNode node = bHsp[offH];
				int temp = node.r; node.r = node.c; node.c = temp;
				temp = node.offH; node.offH = node.offV; node.offV = temp;
			}
		}
		NSPArray[] aTemp = T.Hsp; T.Hsp = T.Vsp; T.Vsp = aTemp;			// swap horisontal & vertical sparse arrays
		int d = T.M; T.M = T.N; T.N = d;
		return T;
	}

	
	
	public NSPMatrix multiply(NSPMatrix T) {
		if (N != T.M) throw new InvalidParameterException("NSPMatrix.multiply(): Nonmatching matrix dimensions.");

		// results stored in new NSPMatrix U, which is complex if both in-matrices are complex
		NSPMatrix U = new NSPMatrix("M", M, T.N);
		U.clearNull();
		NSPNode node;
		U.name = (Matrix.DEBUG_LEVEL > 1 ? "(" + name + "." + T.name + ")" : "M" + nameCount++);
		int[] offsetV = new int[M];									// array will keep an incremental offset into all Vsp sparse column buffers

		for (int i = 0; i < M; i++) {								// for each row of target matrix
			NSPArray aHspU = U.Hsp[i];
			
			for (int j = 0, j1 = 0; j < T.N; j++) {					// for each column of target matrix		
				double v = multiplyHVsp(Hsp[i], T.Vsp[j], 0, 1);	// multiplies sparse row i with sparse column j
				if (!nearZero(v)) {									// if result != 0, then we need to insert a value into target matrix
					NSPArray aVspU = U.Vsp[j];
					updateArraySize(aHspU.nodes + 1, aHspU);		// increment target horisontal array and check need for reallocation
					updateArraySize(aVspU.nodes + 1, aVspU);		// increment target vertical array and check need for reallocation
					aHspU.array[j1] = node = new NSPNode(i, j, v, 0, j1++, offsetV[j]);
					aVspU.array[offsetV[j]++] = node;
					U.nNZ++;
					if (i == j) U.pivotNsp[i] = node;				// if it's a pivot, insert in fast-access array
				}
			}
		}
		return U;
	}
	
	
	
	// returns Markowitz cost for this element
	public int costMarkowitz(int r, int c) {
		return (Hsp[r].nodes - 1) * (Vsp[c].nodes - 1);
	}
	
	
	// LU decomposer optimised for sparse NSP matrices, note: significant speedups gained only for sparse (<25% NNZ) matrices,
	// as this type of solvers suffer from inner-loop indirect addressing, the typical Achilles heel for sparse non-frontal solvers
	// adapted from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
	public NSPMatrix[] decomposeLU(boolean copy, boolean generateLandU) {
		
		if (M != N)	throw new InvalidParameterException("NSPMatrix.decomposeLU2(): Matrix not square.");
		if (M < 1)	throw new InvalidParameterException("NSPMatrix.decomposeLU2(): Invalid matrix.");
		
		NSPMatrix[] lLU = new NSPMatrix[2];
		int[] mutatorLU = null;
		NSPMatrix LU = this;
		int swaps = 0;
		
		if (copy) {
			LU = this.clone();
			LU.name = "LU" + nameCount++;
			mutatorLU = LU.mutator[0] = new int[M];
		} else {
			name = "LU" + nameCount++;			// this matrix will be modified, change it's name
			mutatorLU = mutator[0] = new int[M];
		}
		
		double[] vv = new double[M];			// vv stores the implicit scaling of each row
		double d = 1;							// every permutation toggles sign of d

		for (int i = 0; i < M; i++) {			// looping over rows getting scaling factor of every row, putting it in vv
			double vMax = 0;
			NSPNode[] bHsp = LU.Hsp[i].array;
			for (int j = 0, jEnd = LU.Hsp[i].nodes; j < jEnd; j++) {
				double v = bHsp[j].v;
				if (v < -vMax || v > vMax) vMax = (v < 0 ? -v : v);
			}
			if (nearZero(vMax)) { status |= SINGULAR_MATRIX; return null; }		// if singular matrix found, return null
			vv[i] = 1.0 / vMax;													// save scaling factor
			mutatorLU[i] = i;
		}
		
		for (int j = 0, jm1 = -1; j < M; j++, jm1++) {
			NSPArray aVsp = LU.Vsp[j];
			NSPNode[] bVsp = aVsp.array;
			int i = 0, offV = 0, nNodes = aVsp.nodes;							// we will continue using i & offH in next loop
			for (int im1 = -1; i < j; i++, im1++) {
				
				double sum = 0;
				if (bVsp[0].r > im1 || LU.Hsp[i].array[0].c > im1) {			// skip row/column mult. if one of the ranges are all zeroes
					if (offV < nNodes && bVsp[offV].r == i) offV++;
					continue;
				}
				if (offV < nNodes && bVsp[offV].r == i) {						// case of operating on a non-zero node
					sum = bVsp[offV].v; 										// sum = LU[i][j]
					if (i == 0) offV++; else {
						sum -= multiplyLUHVsp(LU.Hsp[i], aVsp, 0, 1, im1);		// for(k=1;k<i;k++) sum -= LU[i][k]*LU[k](j] (inner sparse product)
						bVsp[offV++].v = sum;
					}
				} else if (i > 0) {												// node was zero, we'll need to insert result
					sum -= multiplyLUHVsp(LU.Hsp[i], aVsp, 0, 1, im1);			// for(k=1;k<i;k++) sum -= LU[i][k]*LU[k](j] (inner sparse product)
					if (!nearZero(sum)) {
						NSPNode newNsp = new NSPNode(i, j, sum, 0, 0, offV);	// we know the vertical offset for the new node
						LU.nNZ++;
						LU.insertHVspNode(i, j, 0, newNsp);						// horisontal insertion
						LU.insertLocalHVspNode(aVsp, offV++, 1, newNsp);		// LU[i][j] = sum
						bVsp = aVsp.array;
						nNodes = aVsp.nodes;
					}
				}
			}
			
			int iMax = j, costMin = Integer.MAX_VALUE;
			double vMax = 0;													// vMax will store the largest pivot element found
			for (; i < M; i++) {												// i continues iterating from j

				// note, this is same operation as before, only for range of (j <= i < M)
				double sum = 0;
				if (offV < nNodes && bVsp[offV].r == i) {						// case of operating on a non-zero node
					sum = bVsp[offV].v; 										// sum = LU[i][j]
					if (j == 0) offV++; else {
						if (bVsp[0].r <= jm1 && LU.Hsp[i].array[0].c <= jm1)
							sum -= multiplyLUHVsp(LU.Hsp[i], aVsp, 0, 1, j - 1);	// for(k=1;k<j;k++) sum -= LU[i][k]*LU[k][j]
						bVsp[offV++].v = sum;
					}
					
					double v = vv[i] * (sum < 0 ? -sum : sum);
					int costM = LU.costMarkowitz(i, j);
					// have we found a better pivot than the best so far?
					if (costM <= costMin && v >= vMax) {
						vMax = v; iMax = i;
						costMin = costM;
						if (DEBUG_LEVEL > 1) System.out.println("Markowitz cost: " + costM + ", v: " + v);
					}

				// node was zero, we'll need to insert result
				} else if (j > 0 && bVsp[0].r <= jm1 && LU.Hsp[i].array[0].c <= jm1) {
					sum -= multiplyLUHVsp(LU.Hsp[i], aVsp, 0, 1, j - 1);		// for(k=1;k<i;k++) sum -= LU[i][k]*LU[k](j]
					if (!nearZero(sum)) {
						NSPNode newNsp = new NSPNode(i, j, sum, 0, 0, offV);	// we know the vertical offset for the new node
						LU.insertHVspNode(i, j, 0, newNsp);						// horisontal insertion
						LU.insertLocalHVspNode(aVsp, offV++, 1, newNsp);		// LU[i][j] = sum
						LU.nNZ++;
						bVsp = aVsp.array;										// array could have changed allocation, reassign bVsp to it
						nNodes = aVsp.nodes;

						double v = vv[i] * (sum < 0 ? -sum : sum);
						int costM = LU.costMarkowitz(i, j);
						if (costM <= costMin && v >= vMax) {
							vMax = v; iMax = i;
							costMin = costM;
							if (DEBUG_LEVEL > 1) System.out.println("Markowitz cost: " + costM + ", v: " + v);
						}
					}
				}

			}
			if (j != iMax) {									// test if there's a better pivot and we need to swap rows
				swaps++;
				if (DEBUG_LEVEL > 1) System.out.println("Final Markowitz cost: " + costMin + ", vMax: " + vMax);
				LU.swapHVspArrays(iMax, j, 0);					// yes
				if (DEBUG_LEVEL > 1) System.out.println("Swapped " + iMax + " and " + j);
				d = -d;
				vv[iMax] = vv[j];								// also change scale factor
				int temp = mutatorLU[j];						// swap indexes in row mutator array
				mutatorLU[j] = mutatorLU[iMax];
				mutatorLU[iMax] = temp;
			}

			// zero at pivot element implies a singular matrix, but some applications can get by with tweaking away the zero
			if (LU.pivotNsp[j] == null) {
				if (copy) status |= SINGULAR_MATRIX;			// mark original (unchanged) matrix as singular
				double v = ROUNDOFF_ERROR * (LU.norm1 < 0 ? LU.norm1() : LU.norm1);
				NSPNode newNsp = LU.pivotNsp[j] = new NSPNode(j, j, v, 0);
				LU.insertHVspNode(j, j, 0, newNsp);				// horisontal & vertical insertion
				LU.insertHVspNode(j, j, 1, newNsp);
			} else if (nearZero(LU.pivotNsp[j].v))				// did a zero value node end up at pivot position?
				LU.pivotNsp[j].v = ROUNDOFF_ERROR * (LU.norm1 < 0 ? LU.norm1() : LU.norm1);

			// time to divide pivot's submatrix column by the pivot element
			if (j < M - 1) {
				double v = 1.0 / LU.pivotNsp[j].v;
				bVsp = aVsp.array;
				int i2End = aVsp.nodes;
				int i2 = LU.pivotNsp[j].offV + 1;
				for (; i2 < i2End; i2++) bVsp[i2].v *= v;	
			}
		}
		lLU[0] = LU;

		// use d to calculate determinant, put it in LU
		NSPNode[] pivots = LU.pivotNsp;
		for (int i = 0; i < M; i++) d *= pivots[i].v;
		LU.det = d;

		// was generation of individual L & U matrices requested?
		if (generateLandU) {
			// convert LU into L, moving it's upper triangle data into U
			lLU[0].name = "L";									// here we treat dataLU as if containing data of L
			lLU[1] = new NSPMatrix("U", M, N, Matrix.Type.Null, 1);

			for (int i = 0; i < M; i++) {
				NSPArray aHspLU = LU.Hsp[i];
				NSPNode[] bHspLU = aHspLU.array, bHspU = lLU[1].Hsp[i].array;
				for (int j = LU.pivotNsp[i].offH; j < aHspLU.nodes; j++) {	// iterate over the top triangle

					bHspU[j] = bHspLU[j];
					if (bHspLU[j].r == bHspLU[j].c)  bHspLU[j].v = 1;		// L's diagonal has only 1:s
					else bHspLU[j].v = 0;									// L's top triangle is zero
				}
			}
			LU.purgeZeroes();												// purge L's top triangle zeroes
			lLU[0].det = lLU[1].det = d;
			lLU[1].mutator[0] = mutatorLU;
		}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("NSPMatrix.decomposeLU2() result:");
			System.out.println(LU.toString());
			if (generateLandU) System.out.println(lLU[1].toString());
			System.out.println("LU row mutator:");
			System.out.println(Arrays.toString(mutatorLU));
			System.out.println("Row swaps: " + swaps);
		}
		return lLU;
	}
	
	
	
	// unpermute a matrix with mutator values altered
	public void unpermute() {
		if (mutator[0] != null) {
			int[] mutatorR = mutator[0];
			for (int i = M - 1; i >= 0; i--)
				if (mutatorR[i] != i)	swapHVspArrays(i, mutatorR[i], 0);
			halfBandwidth = -1;
		}
		if (mutator[1] != null) {
			int[] mutatorC = mutator[0];
			for (int i = M - 1; i >= 0; i--)
				if (mutatorC[i] != i)		swapHVspArrays(i, mutatorC[i], 1);
			halfBandwidth = -1;
		}
	}

	
	// method takes the combined factorised LU matrix using it to solve a linear system
	// where only the constant vector is changing, the matrix coefficients assumed to be fixed
	// adapted from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
	// the triangular matrices L & U come here in a COMBINED form in matrix LU, split by l(ii)+u(ii) in the diagonal; l(ii) = 1
	// matrix A, LU are applied on this vector, matrix array with solution vector and LU matrix is the return value
	// if LU is null, then the method calls decomposeLU2() with A, if A is null then LU is used
	// method optimisation takes into account that b might contain many leading zero elements
	public NSPMatrix[] backSubstituteLU(NSPMatrix A, NSPMatrix LU, boolean copy) {
		
		if (A == null && LU == null)	throw new RuntimeException("NSPMatrix.backSubstituteLU(): Null matrices/vector.");
		if (N != 1 || M != (A == null ? LU.M : A.M))
										throw new InvalidParameterException("NSPMatrix.backSubstituteLU(): Invalid constant vector.");			
		NSPMatrix[] bLU = new NSPMatrix[2];
		
		// LU = null tells algorithm to create a new LU decomposure, which will be returned in the matrix array
		if (LU == null) {
			if (A.M != A.N)	throw new InvalidParameterException("NSPMatrix.decomposeLU(): Matrix A not square.");			
			if (A.M < 1)	throw new InvalidParameterException("NSPMatrix.decomposeLU(): Invalid matrix.");
			// copy A into a LU matrix and decompose LU, if copy=false it will destroy matrix A returning it as LU
			bLU[1] = A.decomposeLU(copy, false)[0];
			// decompose failed on a singular matrix, return null
			if (bLU[1] == null)	return null;	
		} else
			bLU[1] = LU;
		
		// if copy=true, do not return solution vector in b vector, but in a copy
		if (copy) bLU[0] = this.clone(); else bLU[0] = this;
		if (Matrix.DEBUG_LEVEL > 1) bLU[0].name = "x((" + bLU[1].name + ")^-1*" + name + ")";
		else						bLU[0].name = "x" + nameCount++;

		int ii = -1, N = bLU[1].N;
		int[] mutatorLU = bLU[1].mutator[0];
		double sum = 0;
		
		NSPArray[] lHspLU = bLU[1].Hsp;						// get Hsp array of the LU matrix
		NSPMatrix vectorB = bLU[0];	
		double datab[] = vectorB.getData()[0];				// convert sparse vector b to plain array for quicker depermutation
		
		for (int i = 0; i < N; i++) {
			int ip = mutatorLU[i];
			sum = datab[ip];								// depermute: take a value from it's permuted position
			datab[ip] = datab[i];							// depermute: return former value at this position to it's place
			if (ii >= 0) {									// on turning positive, ii will be the first nonvanishing element of b
				NSPNode[] bHspLU = lHspLU[i].array;
				int nNodesLU = lHspLU[i].nodes;
				int offH = findHVspNode(bHspLU, 0, nNodesLU - 1, -1, ii);
				if (offH < 0) offH = -offH - 1;
				for (; offH < nNodesLU && bHspLU[offH].c <= i - 1; offH++) {
					sum -= bHspLU[offH].v * datab[bHspLU[offH].c];					// sum -= a[i][j] * b[j]
				}
			} else if (!nearZero(sum)) ii = i;
			datab[i] = sum;
		}
						
		NSPNode[] pivotNspA = bLU[1].pivotNsp;
		
		for (int i = N - 1; i >= 0; i--) {											// backsubstitution pass
			sum = datab[i];
			NSPNode[] bHspLU = lHspLU[i].array;
			int nNodesLU = lHspLU[i].nodes;
			int offH = pivotNspA[i].offH + 1;
			for (; offH < nNodesLU; offH++) {
				sum -= bHspLU[offH].v * datab[bHspLU[offH].c];						// sum -= a[i][j] * b[j]
			}
			// store an element of solution vector X
			datab[i] = sum / pivotNspA[i].v;
		}
		vectorB.putData(datab, null);
				
		if (DEBUG_LEVEL > 1) {
			System.out.println("Backsubstitution LU solver:");
			System.out.println(bLU[0].toString());
			if (A != null) System.out.println(bLU[1].toString());
		}
		return bLU;
	}

	

	// calculates 1-norm of a NSPMatrix
	public double norm1() {
		for (int j = 0; j < N; j++) {
			double sum = 0;
			NSPNode[] bVsp = Vsp[j].array;
			int nNodes = Vsp[j].nodes;
			for (int offH = 0; offH < nNodes; offH++) {
				double v = bVsp[offH].v;
				sum += v < 0 ? -v : v; }
			if (sum > norm1) norm1 = sum;
		}
		return norm1;
	}

	// calculates Frobenius norm of a NSPMatrix
	public double normFrobenius() {
		double sum = 0;
		for (int j = 0; j < N; j++) {
			NSPNode[] bVsp = Vsp[j].array;
			int nNodes = Vsp[j].nodes;
			for (int offH = 0; offH < nNodes; offH++) { double v = bVsp[offH].v; sum += v * v; }
		}
		return Math.sqrt(sum);
	}

	

	
	
	// method adds two supplied sparse rows into a new row, the horisontal/vertical aspect chosen individually
	// the new row is returned for integration into a NSP scheme
	static NSPArray addHVsp(NSPArray aHVspA, NSPArray aHVspB, double f, int aspectA, int aspectB, int aspectC) {
		
		if (aHVspA.array == null && aHVspB.array == null) return null;
		int nNa = aHVspA.nodes, nNb = aHVspB.nodes;
		int nNc = trimToAllocBlock(nNa + nNb), a, b;
		NSPNode[] bHVspC = new NSPNode[nNc], bHVspA = aHVspA.array, bHVspB = aHVspB.array;
		int maxN = Integer.MAX_VALUE;
		
		int ic = 0;
		for (int ia = 0, ib = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
			NSPNode nA = bHVspA == null ? null : bHVspA[ia];
			NSPNode nB = bHVspB == null ? null : bHVspB[ib];
			if (aspectA == 0) {
				if (aspectB == 0) 	{ a = (nA == null ? maxN : nA.c); b = (nB == null ? maxN : nB.c); }
				else				{ a = (nA == null ? maxN : nA.c); b = (nB == null ? maxN : nB.r); }
			} else {
				if (aspectB == 0) 	{ a = (nA == null ? maxN : nA.r); b = (nB == null ? maxN : nB.c); }
				else				{ a = (nA == null ? maxN : nA.r); b = (nB == null ? maxN : nB.r); }
			}
			if (a == b) {
				if (aspectC == 0)
						{ bHVspC[ic] = new NSPNode(0, a, nA.v + nB.v * f, nA.iv + nB.iv * f); bHVspC[ic].offH = ic++; }
				else	{ bHVspC[ic] = new NSPNode(a, 0, nA.v + nB.v * f, nA.iv + nB.iv * f); bHVspC[ic].offV = ic++; }
				ia++; ib++;
			} else if (a < b) {
				if (aspectC == 0)
						{ bHVspC[ic] = new NSPNode(0, a, nA.v, nA.iv); bHVspC[ic].offH = ic++; }
				else	{ bHVspC[ic] = new NSPNode(a, 0, nA.v, nA.iv); bHVspC[ic].offV = ic++; }
				ia++;
			} else {
				if (aspectC == 0)
						{ bHVspC[ic] = new NSPNode(0, b, nB.v * f, nB.iv * f); bHVspC[ic].offH = ic++; }
				else	{ bHVspC[ic] = new NSPNode(b, 0, nB.v * f, nB.iv * f); bHVspC[ic].offV = ic++; }
				ib++;	
			}
		}
		return new NSPArray(ic, nNc, bHVspC);
	}

	
	
	// sparse CSR2 inner vector product, with chosen horisontal/vertical (Hsp/Vsp) multiplying aspect
	static double multiplyHVsp(NSPArray aHVspA, NSPArray aHVspB, int aspectA, int aspectB) {
		
		NSPNode[] bHVspA = aHVspA.array, bHVspB = aHVspB.array;
		if (bHVspA == null || bHVspB == null) return 0;				// base case: one of the arrays is all zeroes
		int nNa = aHVspA.nodes, nNb = aHVspB.nodes;
		double v = 0;
	
		for (int iA = 0, iB = 0; iA < nNa || iB < nNb;) {			// march through each node in the two arrays
			NSPNode nA = bHVspA[iA];
			if (nA == null) break;									// if one of nodes is null, then we've finished one of the arrays			
			NSPNode nB = bHVspB[iB];
			if (nB == null) break;		
			if (aspectA == 0) {
				if (aspectB == 0) {
					if (nA.c == nB.c) { v += nA.v * nB.v; iA++; iB++; }
					else if (nA.c < nB.c) iA++; else iB++;
				} else {
					if (nA.c == nB.r) { v += nA.v * nB.v; iA++; iB++; }
					else if (nA.c < nB.r) iA++; else iB++;
				}
			} else {
				if (aspectB == 0) {
					if (nA.r == nB.c) { v += nA.v * nB.v; iA++; iB++; }
					else if (nA.r < nB.c) iA++; else iB++;
				} else {
					if (nA.r == nB.r) { v += nA.v * nB.v; iA++; iB++; }
					else if (nA.r < nB.r) iA++; else iB++;
				}				
			}
		}
		return v;
	}

	
	
	// sparse CSR2 inner vector product for LU decomposition algorithm
	static double multiplyLUHVsp(NSPArray aHVspA, NSPArray aHVspB, int aspectA, int aspectB, int k) {
		
		int nNa = aHVspA.nodes, nNb = aHVspB.nodes;
		//if (nNa == 0 || nNb == 0) return 0;					// base case: one of the arrays is all zeroes
		NSPNode[] bHVspA = aHVspA.array, bHVspB = aHVspB.array;
		double v = 0;
	
		for (int ia = 0, ib = 0; ia < nNa && ib < nNb;) {	// march through each node in the two arrays
			NSPNode nA = bHVspA[ia], nB = bHVspB[ib];	
			if (aspectA == 0) {
				if (aspectB == 0) {
					if (nA.c > k || nB.c > k) return v;
					if (nA.c == nB.c) { v += nA.v * nB.v; ia++; ib++; }
					else if (nA.c < nB.c) ia++; else ib++;
				} else {
					if (nA.c > k || nB.r > k) return v;
					if (nA.c == nB.r) { v += nA.v * nB.v; ia++; ib++; }
					else if (nA.c < nB.r) ia++; else ib++;
				}
			} else {
				if (aspectB == 0) {
					if (nA.r > k || nB.c > k) return v;
					if (nA.r == nB.c) { v += nA.v * nB.v; ia++; ib++; }
					else if (nA.r < nB.c) ia++; else ib++;
				} else {
					if (nA.r > k || nB.r > k) return v;
					if (nA.r == nB.r) { v += nA.v * nB.v; ia++; ib++; }
					else if (nA.r < nB.r) ia++; else ib++;
				}				
			}
		}
		return v;
	}

	
	// sparse CSR2 inner vector product of an index range between start (iS) & end (iE)
	static double multiplyStartEndHVsp(NSPArray aHVspA, NSPArray aHVspB, int aspectA, int aspectB, int iS, int iE) {
		
		int nNa = aHVspA.nodes, nNb = aHVspB.nodes, iA = 0, iB = 0;
		//if (nNa == 0 || nNb == 0) return 0;					// base case: one of the arrays is all zeroes
		NSPNode[] bHVspA = aHVspA.array, bHVspB = aHVspB.array;
		double v = 0;
	
		if (aspectA == 0)	iA = findHVspNode(bHVspA, 0, nNa, -1, iS);
		else				iA = findHVspNode(bHVspA, 0, nNa, iS, -1);
		if (aspectB == 0)	iB = findHVspNode(bHVspB, 0, nNb, -1, iS);
		else				iB = findHVspNode(bHVspB, 0, nNb, iS, -1);

		for (; iA < nNa && iB < nNb;) {	// march through each node in the two arrays
			NSPNode nA = bHVspA[iA], nB = bHVspB[iB];	
			if (aspectA == 0) {
				if (aspectB == 0) {
					if (nA.c > iE || nB.c > iE) return v;
					if (nA.c == nB.c) { v += nA.v * nB.v; iA++; iB++; }
					else if (nA.c < nB.c) iA++; else iB++;
				} else {
					if (nA.c > iE || nB.r > iE) return v;
					if (nA.c == nB.r) { v += nA.v * nB.v; iA++; iB++; }
					else if (nA.c < nB.r) iA++; else iB++;
				}
			} else {
				if (aspectB == 0) {
					if (nA.r > iE || nB.c > iE) return v;
					if (nA.r == nB.c) { v += nA.v * nB.v; iA++; iB++; }
					else if (nA.r < nB.c) iA++; else iB++;
				} else {
					if (nA.r > iE || nB.r > iE) return v;
					if (nA.r == nB.r) { v += nA.v * nB.v; iA++; iB++; }
					else if (nA.r < nB.r) iA++; else iB++;
				}				
			}
		}
		return v;
	}

	
		
	
	
	// swaps two sparse arrays (=row/column) of index a & b, Hsp (2 rows) if aspect = 0, Vsp (2 columns) if aspect = 1
	// method marches through both swapped sparse arrays within OPPOSITE aspect, if there are values in both positions then they're trivially
	// swapped, if a value is swapped with a zero, then the value is shifted along the opposite aspect to it's swap position
	// finally, the two array references within the INTENDED aspect are swapped
	public void swapHVspArrays(int a, int b, int aspect) {
		
		if (a == b) return;
		NSPNode[] bHVspA = null, bHVspB = null;
		NSPNode nA, nB, nTemp = null;
		int nNa = 0, nNb = 0, temp;
		
		switch (aspect) { 
			// switch two rows (Hsp aspect)
			case 0:
				bHVspA = Hsp[a].array; bHVspB = Hsp[b].array; nNa = Hsp[a].nodes; nNb = Hsp[b].nodes;
				
				for (int ia = 0, ib = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
					
					nA = bHVspA == null || ia == nNa ? null : bHVspA[ia];
					nB = bHVspB == null || ib == nNb ? null : bHVspB[ib];
					// provide a fake boundary if one of the marching indexes run out of nodes
					int nAc = (nA == null ? N: nA.c), nBc = (nB == null ? N : nB.c);
					
					if (nAc == nBc) {										// found matching column nodes?
						if (nA.r == nAc) pivotNsp[nAc] = nB;				// if A was pivot, make B pivot
						else if (nB.r == nBc) pivotNsp[nBc] = nA;			// if B was pivot, make A pivot

						nTemp = Vsp[nAc].array[nA.offV];					// exchange their Vsp references
						Vsp[nAc].array[nA.offV] = Vsp[nAc].array[nB.offV];
						Vsp[nAc].array[nB.offV] = nTemp;
						temp = nA.offV; nA.offV = nB.offV; nB.offV = temp;	// exchange their offsets within the arrays	
						nA.r = b; nB.r = a;									// exchange their row indexes
						ia++; ib++;
						
					} else if (nAc < nBc)	{								// array A behind array B (meaning, zero at array B)?
						if (nA.r == nAc) pivotNsp[a] = null;				// see if A was a pivot
						if (nA.c == b) pivotNsp[b] = nA;					// see if A's element ends up in pivot position in row b
						nA.r = b;											// change it's row to the other one
						NSPNode[] bVsp = Vsp[nAc].array;
						if (a < b) {										// do the switch within higher part of array b
							for (int j = nA.offV, k = j + 1, kEnd = Vsp[nAc].nodes; j < kEnd; j++, k++) {
								if (k == kEnd) { nA.offV = kEnd - 1; bVsp[j] = nA; break; }
								if (bVsp[k].r > b) {						// if next Vsp row index is higher than array b's
									nA.offV = bVsp[j].offV;
									bVsp[j] = nA;							// then insert the swapped node at current place in sparse array
									break;
								}
								bVsp[j] = bVsp[k];							// otherwise shift nodes backwards
								bVsp[j].offV--;								// decrease their offsets into Vsp
							}
						} else {
							for (int j = nA.offV, k = j - 1; j >= 0; j--, k--) {
								if (k < 0) { nA.offV = 0; bVsp[j] = nA; break; }
								if (bVsp[k].r < b) {						// if next Vsp row index is lower than array b's
									nA.offV = bVsp[j].offV;
									bVsp[j] = nA;							// then insert the swapped node at current place in sparse array
									break;
								}
								bVsp[j] = bVsp[k];							// otherwise shift nodes forwards
								bVsp[j].offV++;								// increase their offsets into Vsp array
							}
						}
						ia++;
					} else {												// seems like nA.c > nB.c, do the same routine, but for nB
						if (nB.r == nBc) pivotNsp[b] = null;				// see if B was a pivot
						if (nB.c == a) pivotNsp[a] = nB;					// see if B's element ends up in pivot position in row a
						nB.r = a;											// change it's row to the other one
						NSPNode[] bVsp = Vsp[nBc].array;
						if (a < b) {										// do the switch within lower part of array b
							for (int j = nB.offV, k = j - 1; j >= 0; j--, k--) {
								if (k < 0) { nB.offV = 0; bVsp[j] = nB; break; }
								if (bVsp[k].r < a) {						// if next Vsp row index is lower than array b's
									nB.offV = bVsp[j].offV;
									bVsp[j] = nB;							// then insert the swapped node at current place in sparse array
									break;
								}
								bVsp[j] = bVsp[k];							// otherwise shift nodes forwards
								bVsp[j].offV++;								// increase their offsets into Vsp array
							}
						} else {
							for (int j = nB.offV, k = j + 1, kEnd = Vsp[nBc].nodes; j < kEnd; j++, k++) {
								if (k == kEnd) { nB.offV = kEnd - 1; bVsp[j] = nB; break; }
								if (bVsp[k].r > a) {						// if next Vsp row index is higher than array a's
									nB.offV = bVsp[j].offV;
									bVsp[j] = nB;							// then insert the swapped node at current place in sparse array
									break;
								}
								bVsp[j] = bVsp[k];							// otherwise shift nodes backwards
								bVsp[j].offV--;								// decrease their offsets into Vsp
							}
						}
						ib++;
					}
				}
				NSPArray aTemp = Hsp[a]; Hsp[a] = Hsp[b]; Hsp[b] = aTemp;				// finally, exchange the rows in Hsp
				break;	
				
			// switch two columns (Vsp aspect)
			case 1:
				bHVspA = Vsp[a].array; bHVspB = Vsp[b].array; nNa = Vsp[a].nodes; nNb = Vsp[b].nodes;
				for (int ia = 0, ib = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
					
					nA = bHVspA == null ? null : bHVspA[ia];
					nB = bHVspB == null ? null : bHVspB[ib];
					// provide a fake boundary if one of the marching indexes run out of nodes
					int nAr = (nA == null ? M: nA.r), nBr = (nB == null ? M: nB.r);
					
					if (nAr == nBr) {										// found matching row nodes?
						if (nA.c == nAr) pivotNsp[nAr] = nB;				// if A was pivot, make B pivot
						else if (nB.c == nBr) pivotNsp[nBr] = nA;			// if B was pivot, make A pivot

						nTemp = Hsp[nAr].array[nA.offH];					// exchange their Hsp references
						Hsp[nAr].array[nA.offH] = Hsp[nAr].array[nB.offH];
						Hsp[nAr].array[nB.offH] = nTemp;
						temp = nA.offH; nA.offH = nB.offH; nB.offH = temp;	// exchange their offsets within the arrays	
						nA.c = b; nB.c = a;									// exchange their row or column indexes
						ia++; ib++;
						
					} else if (nAr < nBr)	{								// array A behind array B (meaning, zero at array B)?
						if (nA.c == nAr) pivotNsp[a] = null;				// see if A was a pivot
						if (nA.r == b) pivotNsp[b] = nA;					// see if A's element ends up in pivot position in column b
						nA.c = b;											// change it's column to the other one
						NSPNode[] bHsp = Hsp[nAr].array;
						if (a < b) {										// do the switch within right part (relative A) of array B
							for (int j = nA.offH, k = j + 1, kEnd = Hsp[nAr].nodes; j < kEnd; j++, k++) {
								if (k == kEnd) { nA.offH = kEnd - 1; bHsp[j] = nA; break; }
								if (bHsp[k].c > b) {						// if next Hsp column index is higher than array B's
									nA.offH = bHsp[j].offH;
									bHsp[j] = nA;							// then insert the swapped node at current place is sparse array
									break;
								}
								bHsp[j] = bHsp[k];							// otherwise shift nodes backwards
								bHsp[j].offH--;								// decrease their offsets into Hsp array
							}
						} else {
							for (int j = nA.offH, k = j - 1; j >= 0; j--, k--) {
								if (k < 0) { nA.offH = 0; bHsp[j] = nA; break; }
								if (bHsp[k].c < b) {						// if next Hsp column index is lower than array B's
									nA.offH = bHsp[j].offH;
									bHsp[j] = nA;							// then insert the swapped node at current place is sparse array
									break;
								}
								bHsp[j] = bHsp[k];							// otherwise shift nodes forwards
								bHsp[j].offH++;								// increase their offsets into Hsp array
							}
						}
						ia++;
					} else {												// seems like nA.c > nB.c, do the same routine, but for nB
						if (nB.c == nBr) pivotNsp[b] = null;				// see if B was a pivot
						if (nB.r == a) pivotNsp[b] = nB;					// see if B's element ends up in pivot position in column a
						nB.c = a;											// change it's column  to the other one
						NSPNode[] bHsp = Hsp[nB.r].array;
						if (a < b) {										// do the switch within left part (relative A) of array B
							for (int j = nB.offH, k = j - 1; j >= 0; j--, k--) {
								if (k < 0) { nB.offH = 0; bHsp[j] = nB; break; }
								if (bHsp[k].c < a) {				// if next Hsp column index is lower than array B's
									nB.offH = bHsp[j].offH;
									bHsp[j] = nB;							// then insert the swapped node at current place is sparse array
									break;
								}
								bHsp[j] = bHsp[k];							// otherwise shift nodes forwards
								bHsp[j].offH++;								// increase their offsets into Vsp
							}
						} else {
							for (int j = nB.offH, k = j + 1, kEnd = Hsp[nBr].nodes; j < kEnd; j++, k++) {
								if (k == kEnd) { nB.offH = kEnd - 1; bHsp[j] = nB; break; }
								if (bHsp[k].c > a) {						// if next Hsp column index is higher than array B's
									nB.offH = bHsp[j].offH;
									bHsp[j] = nB;							// then insert the swapped node at current place is sparse array
									break;
								}
								bHsp[j] = bHsp[k];							// otherwise shift nodes backwards
								bHsp[j].offH--;								// decrease their offsets into Vsp
							}
						}
						ib++;
					}
				}
				aTemp = Vsp[a]; Vsp[a] = Vsp[b]; Vsp[b] = aTemp;			// finally, exchange the column in Vsp
				break;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			DATA HANDLING METHODS METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	private static int trimToAllocBlock(int v) { return (v & (0xFFFFFFFF - CSR2_ALLOCBLOCK + 1)) + CSR2_ALLOCBLOCK; }
	
	
	// method converts input matrix A into a NSPMatrix
	public static Matrix convert(Matrix A) {
		return new NSPMatrix(A.name, A.M, A.N, A.getData()[0], A.getData()[1]);
	}
	
	
	public void zero() { initNSPMatrix(); setNull(); }


	
	@Override
	public double[][] getDataRef() { return getData(); }

	// method returns NSP-style data as ordinary linear array Matrix data
	@Override
	public double[][] getData() {	
		double[][] dataSet = new double[2][];
		double[] data = dataSet[0] = new double[M * N], idata = null;
		if (isComplex()) idata = dataSet[1] = new double[M * N];
		
		int rN = 0;
		for (NSPArray aHsp : Hsp) {
			NSPNode[] bHsp = aHsp.array;
			if (bHsp != null) {
				for (int i = 0; i < aHsp.nodes; i++)
					if (bHsp[i] != null) data[rN + bHsp[i].c] = bHsp[i].v;
				if (isComplex()) for (int i = 0; i < aHsp.nodes; i++)
					if (bHsp[i] != null) idata[rN + bHsp[i].c] = bHsp[i].iv;
			}
			rN += N;
		}
		return dataSet;
	}

	@Override
	public void putDataRef(double[] data, double[] idata) { putData(data, idata); }
	
	// method converts flat array matrix data into NSP format
	@Override
	public void putData(double[] data, double[] idata) {
				
		initNSPMatrix();		// clear NSP data, we're reinserting from scratch

		// first loop inserts elements in proper order into Hsp
		for (int i = 0; i < M; i++) {

			NSPArray aHsp = Hsp[i];
			NSPNode[] bHsp = aHsp.array;
			int iN = i * N, offsHsp = 0;							// offsHsp = offset into current Hsp buffer
			
			// iterate over values in the matrix row
			for (int j = iN, jEnd = iN + N, col = 0; j < jEnd; j++, col++) {
								
				// we're inserting only nonzero values
				if (!nearZero(data[j]) || (isComplex() && !nearZero(idata[j]))) {
					
					// if current sparse Hsp array is full, add another CSR2_ALLOCBLOCK elements
					if (updateArraySize(aHsp.nodes + 1, aHsp)) bHsp = aHsp.array;

					bHsp[offsHsp] = new NSPNode(i, col, data[j], isComplex() ? idata[j] : 0);	// insert column index
					nNZ++;											// global node count incremented
					if (i == col) pivotNsp[i] = bHsp[offsHsp];		// if it's a pivot, put in fast-access array
					
					bHsp[offsHsp].offH = offsHsp;					// store the node's reference offset in Hsp
					offsHsp++;
					readjustHalfBandwidth(i, col);					// readjust half bandwidth if value coordinates are "sticking out" from diagonal
				}
			}
			if (aHsp.nodes == 0) Hsp[i].array = null;				// dereference empty rows
		}

		crosslinkVsp();
		clearNull();
	}

	
	
	// method utilises a BipartiteDM graph object generated from a matrix to permute nonzeroes to the diagonal in same matrix
	// method alters the mutator reindexer accordingly
	public NSPMatrix permuteMaximumMatching(BipartiteDM bDM, boolean copy, boolean swapColumns) {
		NSPMatrix Mp = new NSPMatrix(this.name + "_permDM", M, N);
		int[] pairing = bDM.pairing;
		int partC = bDM.vertsR + 1, vertsR = bDM.vertsR;
		
		if (swapColumns) {
			Mp.mutator[1] = new int[N];
			for (int i = 0, j = 0; i < M; i++) {
				
				int p = pairing[i];
				if (p == -1) {
					// we 've found an unmatched end C-vertex, backtrack by C->R pairings to the unmatched start C-vertex
					p = pairing[i + partC];
					for (int p2 = p; p2 != vertsR; p2 = pairing[p2 + partC])
						p = p2;
				}
				
				Mp.mutator[1][j] = p;									// store resulting permutation in the mutator's columns array
				Mp.Vsp[j] = copy ? Vsp[p].clone() : Vsp[p];				// use the row's paired column index to place a nz on the diagonal
				
				NSPNode[] bVsp = Mp.Vsp[j].array;						// reindex the column to it's new column indexing
				for (int i2 = 0, i2End = Mp.Vsp[j].nodes; i2 < i2End; i2++) bVsp[i2].c = j;
				j++;
			}
			Mp.crosslinkHsp();											// we've set up the columns, reestablish the rows
		} else {
			Mp.mutator[0] = new int[M];
			for (int j = 0, i = 0; j < N; j++) {
				
				int p = pairing[j + partC];
				if (p == vertsR) {
					// we 've found an unmatched end R-vertex, backtrack by R->C pairings to the unmatched start R-vertex
					p = pairing[i];
					for (int p2 = p; p2 != -1; p2 = pairing[p2]) {
						p = p2;
					}
				}
				
				Mp.mutator[0][i] = p;									// store resulting permutation in the mutator's rows array
				Mp.Hsp[i] = copy ? Hsp[p].clone() : Hsp[p];				// use the column's paired row index to place a nz on the diagonal
				
				NSPNode[] bHsp = Mp.Hsp[i].array;						// reindex the row to it's new row indexing
				for (int j2 = 0, j2End = Mp.Hsp[i].nodes; j2 < j2End; j2++) bHsp[j2].r = i;
				i++;
			}
			Mp.crosslinkVsp();											// we've set up the rows, reestablish the columns
		}
		return Mp;
	}
	
	
	
	// method converts a frontal matrix DAG structure into NSP format
	public void getFrontalData(FrontalDAG tDAG) {
		
		initNSPMatrix();											// clear NSP data, we're reinserting from scratch

		// we're iterating over the frontal matrices in their original master-pivot-index order
		// but recall that there has possibly been internal pivotings AND failed pivots transmigrations to LU-parents
		// thus, the current pivot blocks describe not L & U but the permuted L' & U', and the problem vector
		// will have to be permuted in SAME row permutation order, THEN the resultant solution vector
		// has to be depermuted in REVERSE column permutation order
		FrontalDAGVertex[] vertexL = tDAG.vertex;
		int pivotsFM2, pivotsFM = 0;
		int[] mutatorR = new int[M], mutatorC = new int[N];
		
		for (int i = 0; i < M; i += pivotsFM) {

			// iterate over the pivot block of frontal matrix of current vertex/pivot, except the failed pivots
			// any failed pivots that may exist have been shifted to new positions further on
			FrontalMatrix fm = vertexL[i].fm;
			double[] dataFM = fm.data;
			pivotsFM2 = fm.nExtraP + fm.nPivots;
			pivotsFM = pivotsFM2 - fm.nFailed;
			
			// extract the permuted global indexing pivot order from frontal's idxCR array, excluding the failed pivots
			// this order will be used to permute the problem and solution vectors
			// note: any previous permutation is taken into account by permuting the CURRENT mutator arrays into NEW mutator arrays
			// the internal global indexes of the frontal arrays were arranged SEQUENTIALLY from start, (despite any previous permutations)
			int[] idxCR = fm.idxCR;
			for (int iFM = 0, i1 = i; iFM < pivotsFM; iFM++, i1++) mutatorR[i1] = mutator[0][idxCR[fm.cTot + iFM]];				
			for (int jFM = 0, j1 = i; jFM < pivotsFM; jFM++, j1++) mutatorC[j1] = mutator[1][idxCR[jFM]];			
			
			for (int iFM = 0; iFM < pivotsFM; iFM++) {
				
				int iP = i + iFM;									// iP is the global pivot offset from current frontal master pivot
				NSPArray aHsp = Hsp[iP];
				int offH = aHsp.nodes, jdFM = iFM * fm.cTot;
				// this aHsp row will be extended with current frontal pivot row and it's trailing contribution row
				updateArraySize(aHsp.nodes + pivotsFM + fm.nContribC, aHsp); 
				NSPNode[] bHsp = aHsp.array;
				
				// first we add the proper pivot block part of the current row, according to the running global pivot index
				for (int jFM = 0; jFM < pivotsFM; jFM++, jdFM++) {
					if (!nearZero(dataFM[jdFM])) {
						bHsp[offH] = new NSPNode(iP, i + jFM, dataFM[jdFM]);
						if (iP == i + jFM) pivotNsp[iP] = bHsp[offH];	// store currently added pivot in fast-access array
						offH++;
						nNZ++;
						readjustHalfBandwidth(iP, iP + jFM);		// readjust half bandwidth if value coordinates are "sticking out" from diagonal
					}
				}			

				// next we add the trailing contributions of current row, and here we need to use the global column indexes from idxCR
				for (int jCR = pivotsFM2; jCR < fm.c; jCR++, jdFM++) {
					if (!nearZero(dataFM[jdFM])) {
						bHsp[offH++] = new NSPNode(iP, idxCR[jCR], dataFM[jdFM]);
						nNZ++;
						readjustHalfBandwidth(iP, idxCR[jCR]);		// readjust half bandwidth if value coordinates are "sticking out" from diagonal
					}
				}
				aHsp.nodes = offH;					// store final size (since zeroes might have been skipped) of this aHsp sparse row				
			}
			
			// we're done with the frontal pivot block and it's rows flank, we have pivot block's columns flank left to add
			for (int iFM = pivotsFM2; iFM < fm.r; iFM++) {
				
				int iP = i + iFM;									// iP is the global pivot offset from current frontal master pivot
				NSPArray aHsp = Hsp[iP];
				int offH = aHsp.nodes, jdFM = iFM * fm.cTot;
				// this aHsp row will be extended with column flank contributions equal to number of nonfailed pivots from pivot block
				updateArraySize(aHsp.nodes + pivotsFM, aHsp); 
				NSPNode[] bHsp = aHsp.array;
				
				for (int jFM = 0; jFM < pivotsFM; jFM++, jdFM++) {
					if (!nearZero(dataFM[jdFM])) {
						bHsp[offH++] = new NSPNode(iP, iP + jFM, dataFM[jdFM]);
						nNZ++;
						readjustHalfBandwidth(iP, iP + jFM);		// readjust half bandwidth if value coordinates are "sticking out" from diagonal
					}
				}				
				aHsp.nodes = offH;					// store final size (since zeroes might have been skipped) of this aHsp sparse row				
			}
			
//			int superP = vertexL[i++].superPivot;				// this code block skips past supernodes
//			if (superP == -1) continue;							// if this was just a single-pivot frontal, then i++ is fine
//			else while(vertexL[i].superPivot == superP) i++;	// otherwise skip past all vertices flagged with same superPivot index
		}

		mutator[0] = mutatorR; mutator[1] = mutatorC; 			// update with new potential permutations
		crosslinkVsp();
		clearNull();
	}

	
	

	static int nodeSearch_iterateLim = 40;							// method switches from iterative to linear search if sparse array > 40
	static int nodeSearch_linearLim = Integer.MAX_VALUE;			// DEBUG: binary search deactivated, appears to be slower than linear
	
	// method searches for a NspNode in a bounded array, according to row if r >= 0, otherwise according to column c
	// if not found returns negated offset to the place where the node is supposed to be within a sorted order
	static int findHVspNode(NSPNode[] bHVsp, int start, int end, int r, int c) {

		if (bHVsp == null) return -1;
		// do ordinary linear search for ranges less than 8
		if (end - start <= nodeSearch_iterateLim) {					// do iterative search if no.of nodes is n <= 60
			if (r >= 0) {
				if (r < bHVsp[start].r) return -(start + 1);		// base case: supplied row is before the entire Vsp array
				if (bHVsp[end].r < r) return -(end + 2);			// base case: supplied row is after the entire Vsp array
				
				for (int j = start; j <= end; j++) {
					if (r < bHVsp[j].r) return -(j+1);
					if (bHVsp[j].r == r) return j;
				}
			} else {
				if (c < bHVsp[start].c) return -(start + 1);		// base case: supplied row is before the entire Hsp array
				if (bHVsp[end].c < c) return -(end + 2);			// base case: supplied row is after the entire Hsp array

				for (int j = start; j <= end; j++) {
					if (c < bHVsp[j].c) return -(j+1);
					if (bHVsp[j].c == c) return j;
				}
			}
		} else if (end - start <= nodeSearch_linearLim) {			// do linear approx.search if no. of nodes is 8 < n <= 999999
			// linear approximation search
			if (r >= 0) {
				if (r < bHVsp[start].r) return -(start + 1);		// base case: supplied row is before the entire Vsp array
				if (bHVsp[end].r < r) return -(end + 2);			// base case: supplied row is after the entire Vsp array

				int rStart = bHVsp[start].r;
				// use intermediate (long)predict to avoid integer overrun on larger numbers of nodes
				long predict1 = ((long)(end - start)*(long)(r - rStart)) / (long)(bHVsp[end].r - rStart);
				int predict = (int)predict1;
				
				if (r < bHVsp[predict].r)
						for (int j = predict; j >= start; j--) {	// predict was larger, scan backward from predict towards target
							if (bHVsp[j].r < r) return -(j+1);		// fell below target's rank without finding it, return negative offset
							if (bHVsp[j].r == r) return j; 
						}
				else	for (int j = predict; j <= end; j++) { 		// predict was smaller, scan forward from predict towards target
							if (r < bHVsp[j].r) return -(j+1);		// rose above target's rank without finding it, return negative offset
							if (bHVsp[j].r == r) return j;
						}
			} else {
				if (c < bHVsp[start].c) return -(start + 1);		// base case: supplied row is before the entire Hsp array
				if (bHVsp[end].c < c) return -(end + 2);			// base case: supplied row is after the entire Hsp array

				int cStart = bHVsp[start].c;
				long predict1 = ((long)(end - start)*(long)(c - cStart)) / (long)(bHVsp[end].c - cStart);
				int predict = (int)predict1;

				if (c < bHVsp[predict].c)							// linear prediction higher than target, backtrack & search iteratively
						for (int j = predict; j >= start; j--) { 
							if (bHVsp[j].c < c) return -(j+1);
							if (bHVsp[j].c == c) return j;
						}
				else	for (int j = predict; j <= end; j++) { 		// linear prediction lower than target, advance & search iteratively
							if (c < bHVsp[j].c) return -(j+1);
							if (bHVsp[j].c == c) return j;
						}				
			}			
		} else {
			// make binary search for sought element in Vsp sparse array
			if (r >= 0) {
				if (r < bHVsp[start].r) return -(start + 1);		// base case: supplied row is before the entire Vsp array
				if (bHVsp[end].r < r) return -(end + 2);			// base case: supplied row is after the entire Vsp array

				int seek = (end + start) >> 1, x = bHVsp[seek].r;
				while (start < end) {
					if (x < r) 		{ start = seek + 1; seek = (end + start) >> 1; }
					else if (x > r)	{ end = seek - 1; seek = (end + start) >> 1; }
					else if (x == r) return seek;
					x = bHVsp[seek].r;
				}
				return -(seek+1);
			// make binary search for sought element in Hsp sparse array
			} else {
				if (c < bHVsp[start].c) return -(start + 1);		// base case: supplied row is before the entire Hsp array
				if (bHVsp[end].c < c) return -(end + 2);			// base case: supplied row is after the entire Hsp array

				int seek = (end + start) >> 1, x = bHVsp[seek].c;
				while (start < end) {
					if (x < c) 		{ start = seek + 1; seek = (end + start) >> 1; }
					else if (x > c)	{ end = seek - 1; seek = (end + start) >> 1; }
					else if (x == c) return seek;
					x = bHVsp[seek].c;
				}
				return -(seek+1);
			}
		}
		return -1;
	}

	
	
	// method chooses the smallest aspect (Hsp vs Vsp) to search for an element
	// if aspect = 0, returns offset into Hsp, if aspect = 1 returns offset into Vsp
	// if element wasn't found, returns negative value with 30th bit set if it's from Vsp, otherwise the bit is clear
	private int locateNspNode(int r, int c, int aspect) {
		
		if (r == c && pivotNsp[r] != null)
			return aspect == 0 ? pivotNsp[r].offH : pivotNsp[r].offV;	// if a pivot if sought, pick from pivotNsp[] array
		
		int nodesH = Hsp[r].nodes - 1, nodesV = Vsp[c].nodes - 1;
		if (nodesH < nodesV) {
			int offH = findHVspNode(Hsp[r].array, 0, nodesH, -1, c);
			if (offH < 0) return offH;
			return aspect > 0 ? Hsp[r].array[offH].offV : offH;
		} else {
			int offV = findHVspNode(Hsp[r].array, 0, nodesV, r, -1);
			if (offV < 0) return offV | 0x40000000;
			return aspect > 0 ? offV : Vsp[c].array[offV].offH;
		}
	}
	

	
	// method inserts node into a Hsp/Vsp reference array, reallocating it if necessary
	// method returns the passed buffer, or the new buffer if reallocation happened
	// the aspect parameter decides way to insert: 0 -> row insert, 1 -> column insert
	// method returns the insertion offset or the found offset negated if node existed
	public int insertHVspNode(int r, int c, int aspect, NSPNode node) {

		NSPNode[] bHVsp = null;
		NSPArray aHVsp = null;
		int offset = 0;
		switch (aspect) {
			case 0:	{	bHVsp = Hsp[r].array; aHVsp = Hsp[r];
						offset = findHVspNode(bHVsp, 0, aHVsp.nodes - 1, -1, c); break; }
			case 1:	{	bHVsp = Vsp[c].array; aHVsp = Vsp[c];
						offset = findHVspNode(bHVsp, 0, aHVsp.nodes - 1, r, -1); break; }
		}
		// if we didn't find element, insert it
		if (offset < 0) {
			offset = -offset-1;
			aHVsp.nodes++;
			if (r == c) pivotNsp[r] = node;					// if pivot, insert in fast-access array

			// does Hsp/Vsp reference array need allocation/reallocation ?
			if (aHVsp.nodes >= aHVsp.size) {
				NSPNode[] bHVsp2 = new NSPNode[aHVsp.size = trimToAllocBlock(aHVsp.nodes)];
				
				int i = aHVsp.nodes - 1;
				// choose inner offset update type, depending on whether we're processing Hsp or Vsp
				if (aspect == 0) {
					for (int j = i - 1; j >= offset; i--, j--) { bHVsp2[i] = bHVsp[j]; bHVsp2[i].offH++; }
					node.offH = offset;
					Hsp[r].array = bHVsp2;
				} else {
					for (int j = i - 1; j >= offset; i--, j--) { bHVsp2[i] = bHVsp[j]; bHVsp2[i].offV++; }
					node.offV = offset;
					Vsp[c].array = bHVsp2;
				}
				bHVsp2[offset] = node;										// reference the node
	
				for (i -= 1; i >= 0; i--) bHVsp2[i] = bHVsp[i];				// move the remaining elements
			} else {
				int i = aHVsp.nodes - 1;
				if (aspect == 0) {
					for (int j = i - 1; j >= offset; i--, j--) { bHVsp[i] = bHVsp[j]; bHVsp[i].offH++; }
					node.offH = offset;
				} else {
					for (int j = i - 1; j >= offset; i--, j--) { bHVsp[i] = bHVsp[j]; bHVsp[i].offV++; }
					node.offV = offset;
				}				
				bHVsp[offset] = node;										// reference the node
			}
			return offset;
		}
		bHVsp[offset].v = node.v;											// node existed, overwrite it's values
		bHVsp[offset].iv = node.iv;
		return -offset - 1;														
	}
	
	
	// method inserts node into a Hsp/Vsp reference array, reallocating array if necessary
	// insertion happens at a known local offset
	private NSPArray insertLocalHVspNode(NSPArray aHVsp, int offset, int aspect, NSPNode node) {

		NSPNode[] bHVsp = aHVsp.array;
		int nNodes = ++aHVsp.nodes;
//		if (bHVsp[offset].r == bHVsp[offset].c)
//			pivotNsp[bHVsp[offset].r] = bHVsp[offset];		// if pivot, insert in fast-access array

		// does Hsp/Vsp reference array need allocation/reallocation ?
		if (nNodes >= aHVsp.size) {
			NSPNode[] bHVsp2 = aHVsp.array = new NSPNode[aHVsp.size = trimToAllocBlock(nNodes)];
			
			int i = nNodes - 1;
			// choose inner offset update type, depending on whether we're processing Hsp or Vsp
			if (aspect == 0) {
				for (int j = i - 1; j >= offset; i--, j--) { bHVsp2[i] = bHVsp[j]; bHVsp2[i].offH++; }
				node.offH = offset;
			} else {
				for (int j = i - 1; j >= offset; i--, j--) { bHVsp2[i] = bHVsp[j]; bHVsp2[i].offV++; }
				node.offV = offset;
			}
			bHVsp2[offset] = node;										// reference the node
			for (i -= 1; i >= 0; i--) bHVsp2[i] = bHVsp[i];				// move the remaining elements
		} else {
			int i = nNodes - 1;
			if (aspect ==0 ) {
				for (int j = i - 1; j >= offset; i--, j--) { bHVsp[i] = bHVsp[j]; bHVsp[i].offH++; }
				node.offH = offset;
			} else {
				for (int j = i - 1; j >= offset; i--, j--) { bHVsp[i] = bHVsp[j]; bHVsp[i].offV++; }
				node.offV = offset;
			}				
			bHVsp[offset] = node;										// reference the node
		}
		return aHVsp;												
	}

	
	
	private int removeHVspNode(int r, int c, int aspect) {

		NSPNode[] bHVsp = null;
		NSPArray aHVsp = null;
		int offset = 0;
		
		switch (aspect) {
			case 0:	{	bHVsp = Hsp[r].array; aHVsp = Hsp[r];
						offset = findHVspNode(bHVsp, 0, aHVsp.nodes-1, -1, c); break; }
			case 1:	{	bHVsp = Vsp[c].array; aHVsp = Vsp[c];
						offset = findHVspNode(bHVsp, 0, aHVsp.nodes-1, r, -1); break; }
		}
		
		// if we found element, remove it
		if (offset >= 0) {
			if (r == c) pivotNsp[r] = null;								// if it's a pivot, remove from fast-access array
			
			int nNodes = --aHVsp.nodes;
			if (nNodes <= 0) {
				aHVsp.size = 0;
				if (aspect == 0) Hsp[r].array = null; else Vsp[c].array = null;		// last node in array removed, destroy this array
				return offset;		
			}
			
			// does Hsp/Vsp reference array need deallocation (over 63 nodes been removed)?
			if (aHVsp.size - nNodes >= CSR2_DEALLOCBLOCK) {
				NSPNode[] bHVsp2 = new NSPNode[aHVsp.size = aHVsp.size - CSR2_DEALLOCBLOCK];
				
				int i = 0;
				if ((aspect & 1) == 0)
						for (; i < offset; i++) { bHVsp2[i] = bHVsp[i];	bHVsp2[i].offH--; }		// move the initial elements in Hsp to new Hsp block
				else	for (; i < offset; i++) { bHVsp2[i] = bHVsp[i];	bHVsp2[i].offV--; }		// move the initial elements in Vsp to new Vsp block
				for (int j = i + 1; j < nNodes; i++, j++) bHVsp2[i] = bHVsp[j];					// shift elements left in new Hsp/Vsp block

			} else {
				if ((aspect & 1) == 0)
						for (int i = offset, j = i + 1; j < nNodes; i++, j++) { bHVsp[i] = bHVsp[j]; bHVsp[i].offH--; }
				else	for (int i = offset, j = i + 1; j < nNodes; i++, j++) { bHVsp[i] = bHVsp[j]; bHVsp[i].offV--; }
			}
		}
		return offset;													// if element didn't exist, return negated offset
	}
	
	
	private void removeLocalHVspNode(NSPArray aHVsp, int offset, int aspect) {
		
		NSPNode[] bHVsp = aHVsp.array;
		if (bHVsp[offset].r == bHVsp[offset].c)
			pivotNsp[bHVsp[offset].r] = null;								// if it's a pivot, remove from fast-access array
		int nNodes = --aHVsp.nodes;
		if (nNodes == 0) { aHVsp.size = 0; aHVsp.array = null; return; }	// last node in array removed, destroy array
		
		// does Hsp/Vsp reference array need deallocation (over 63 nodes been removed)?
		if (aHVsp.size - nNodes >= CSR2_DEALLOCBLOCK) {
			NSPNode[] bHVsp2 = new NSPNode[aHVsp.size = aHVsp.size - CSR2_DEALLOCBLOCK];
			
			int i = 0;
			if ((aspect & 1) == 0)
					for (; i < offset; i++) { bHVsp2[i] = bHVsp[i];	bHVsp2[i].offH--; }		// move the initial elements in Hsp to new Hsp block
			else	for (; i < offset; i++) { bHVsp2[i] = bHVsp[i];	bHVsp2[i].offV--; }		// move the initial elements in Vsp to new Vsp block
			for (int j = i + 1; j < nNodes; i++, j++) bHVsp2[i] = bHVsp[j];					// shift elements left in new Hsp/Vsp block
		} else {
			if ((aspect & 1) == 0)
					for (int i = offset, j = i + 1; i < nNodes; i++, j++) { bHVsp[i] = bHVsp[j]; bHVsp[i].offH--; }
			else	for (int i = offset, j = i + 1; i < nNodes; i++, j++) { bHVsp[i] = bHVsp[j]; bHVsp[i].offV--; }
		}
	}
	
	
	
	// method eliminates zeroes from a sparse structure
	public int purgeZeroes() {
		
		nNZ = 0;													// reset non-zeroes counter
		int maxNodes = 0, foundZ = 0;
		// loop eliminates zeroes in the Hsp arrays, crosslinkVsp() updates the Vsp arrays
		for (int i = 0; i < M; i++) {

			NSPArray aHsp = Hsp[i];
			NSPNode[] bHsp = aHsp.array;
			int nodes2 = 0;
			
			// iterate over values in the matrix row
			for (int j = 0; j < aHsp.nodes; j++) {
								
				// we're keeping only nonzero values
				if (!nearZero(bHsp[j].v) || (isComplex() && !nearZero(bHsp[j].iv))) {
					bHsp[nodes2] = bHsp[j];								// shift reference inside buffer
					nNZ++;												// global node count incremented
					if (bHsp[nodes2].r == bHsp[nodes2].c)
						pivotNsp[bHsp[nodes2].r] = bHsp[nodes2];		// if it's a pivot, put in fast-access array					
					bHsp[nodes2].offH = nodes2++;						// change node's offset in Hsp					
				} else {
					// if nonzero found and it was in pivot position, clear that position in fast-access array
					if (bHsp[j].r == bHsp[j].c) pivotNsp[bHsp[j].r] = null;
					foundZ++;
				}
			}
			if (nodes2 > maxNodes) maxNodes = nodes2;
			updateArraySize(nodes2, aHsp);								// check if Hsp array can be reduced
		}
		crosslinkVsp();
		return foundZ;
	}


	
	// method finalises an NSP structure by crossreferencing Hsp's references in vertical Vsp arrays
	// method is called by methods that have preconstructed a NSPMatrix in the Hsp aspect and need finalising Vsp aspect
	public void crosslinkVsp() {
		
		int nodeCnt = 0;											// counts up number of assembled nodes so far
		int[] idxHsp = new int[M];									// stores indexes into Hsp arrays to use in crosslinking
		int[] linkRow = new int[M];
		int[] offsHsp = new int[M];									// holds current offsets from left into every Hsp array

		// iterate over columns, this sets up a left-to-right marching through sparse column data
		for (int c = 0; c < N; c++) {
			
			int toLink = 0;											// counts number of vertical matches found for crosslinking
			// find all Hsp arrays containing current column c
			for (int r = 0; r < M; r++) {
				NSPNode[] bHsp = Hsp[r].array;
				NSPArray aHsp = Hsp[r];
				if (offsHsp[r] >= 0 && bHsp != null && bHsp[offsHsp[r]].c == c) {
					linkRow[toLink] = r;
					idxHsp[toLink++] = offsHsp[r];
					offsHsp[r]++;										// increment the marching offset into Hsp array
					if (offsHsp[r] >= aHsp.nodes) offsHsp[r] = -1;	// deactivate exhausted rows	
				}
			}
			
			// create the column-wise group
			NSPArray aVsp = Vsp[c];
			updateArraySize(toLink, Vsp[c]);						// shrink Vsp array if it's possible

			// reference all found nodes for this column into current Vsp array
			int nodeCnt2 = nodeCnt + toLink;
			for (int r = 0; r < toLink; r++) {
				aVsp.array[r] = Hsp[linkRow[r]].array[idxHsp[r]];
				aVsp.array[r].offV = r;	
			}
			nodeCnt = nodeCnt2;
		}
	}
	
	
	
	// method finalises an NSP structure by crossreferencing Vsp's references in horisontal Hsp arrays
	// method is called by methods that have preconstructed a NSPMatrix in the Vsp aspect and need finalising Hsp aspect
	public void crosslinkHsp() {
		
		int nodeCnt = 0;											// counts up number of assembled nodes so far
		int[] idxVsp = new int[N];									// stores indexes into Vsp arrays to use in crosslinking
		int[] linkCol = new int[N];
		int[] offsVsp = new int[N];									// holds current offsets from left into every Vsp array

		// iterate over columns, this sets up a top-to-bottom marching through sparse column data
		for (int r = 0; r < M; r++) {
			
			int toLink = 0;											// counts number of vertical matches found for crosslinking
			// find all Vsp arrays containing current row r
			for (int c = 0; c < N; c++) {
				NSPArray aVsp = Vsp[c];
				NSPNode[] bVsp = aVsp.array;
				if (offsVsp[c] >= 0 && bVsp != null && bVsp[offsVsp[c]].r == r) {
					linkCol[toLink] = c;
					idxVsp[toLink++] = offsVsp[c];
					offsVsp[c]++;									// increment the marching offset into Vsp array
					if (offsVsp[c] >= aVsp.nodes) offsVsp[c] = -1;	// deactivate exhausted rows	
				}
			}
			
			// create the row-wise group
			if (toLink == 0) {
				Hsp[r].array = null;
				continue;
			}
			updateArraySize(toLink, Hsp[r]);						// set Hsp array to fit the linked nodes
			NSPArray aHsp = Hsp[r];

			// reference all found nodes for this column into current Hsp array
			int nodeCnt2 = nodeCnt + toLink;
			for (int c = 0; c < toLink; c++) {
				aHsp.array[c] = Vsp[linkCol[c]].array[idxVsp[c]];
				aHsp.array[c].offH = c;	
			}
			nodeCnt = nodeCnt2;
		}
	}
	
	void assembleFromArrays(int[] ii, int[] jj, double[] sA, int[] rCounts) {
		
	}

	
	// method returns sparse column or row indexes as an integer index array, from pivot's row if rowAspect = true, otherwise from column
	// range defines what part of indexes to store:
	// -2: store all indexes prior to pivot index
	// -1: store all indexes prior to pivot index including pivot index
	//  0: store all indexes
	//  1: store all indexes after pivot index including pivot index
	//  2: store all indexes after pivot index
	public int[] indexList(int pivot, boolean rowAspect, int range) {
		
		int i=0, iEnd=0;
		if (rowAspect) {
			NSPNode[] bHsp = Hsp[pivot].array;
			switch (range) {
			case -2: i = 0; iEnd = pivotNsp[pivot].offH; break;
			case -1: i = 0; iEnd = pivotNsp[pivot].offH + 1; break;
			case 0: i = 0; iEnd = Hsp[pivot].nodes; break;
			case 1: i = pivotNsp[pivot].offH; iEnd = Hsp[pivot].nodes; break;
			case 2: i = pivotNsp[pivot].offH + 1; iEnd = Hsp[pivot].nodes; break;
			}
			int[] idx = new int[iEnd - i];
			for (int i2 = 0; i < iEnd; i++, i2++) idx[i2] = bHsp[i].c;
			return idx;
		} else {
			NSPNode[] bVsp = Vsp[pivot].array;
			switch (range) {
			case -2: i = 0; iEnd = pivotNsp[pivot].offV; break;
			case -1: i = 0; iEnd = pivotNsp[pivot].offV + 1; break;
			case 0: i = 0; iEnd = Vsp[pivot].nodes; break;
			case 1: i = pivotNsp[pivot].offV; iEnd = Vsp[pivot].nodes; break;
			case 2: i = pivotNsp[pivot].offV + 1; iEnd = Vsp[pivot].nodes; break;
			}
			int[] idx = new int[iEnd - i];
			for (int i2 = 0; i < iEnd; i++, i2++) idx[i2] = bVsp[i].r;
			return idx;
		}
	}
	
	// method checks if a Hsp or Vsp array is full and needs reallocation and duplication
	// method also checks if a Hsp or Vsp array lost at least 64 elements and can be shrunk
	static boolean updateArraySize(int nodes, NSPArray aHVsp) {
		
		// check if array needs increasing
		if (nodes >= aHVsp.size) {
			NSPNode[] newArray = new NSPNode[aHVsp.size = trimToAllocBlock(nodes)], array = aHVsp.array;
			if (array != null) {
				int nodes2 = aHVsp.nodes;
				for (int i = 0; i < nodes2; i++) newArray[i] = array[i];
			}
			aHVsp.array = newArray;
			aHVsp.nodes = nodes;
			return true;
		}
		// check if array needs shrinking
		if (aHVsp.size - nodes >= CSR2_DEALLOCBLOCK) {
			aHVsp.size = trimToAllocBlock(nodes);
			if (aHVsp.size == 0) { aHVsp.array = null; return true; }		// if last element removed, destroy node array
			NSPNode[] newArray = new NSPNode[aHVsp.size], array = aHVsp.array;
			if (array != null)
				for (int i = 0; i < nodes; i++) newArray[i] = array[i];
			aHVsp.array = newArray;
			aHVsp.nodes = nodes;
			return true;			
		}
		aHVsp.nodes = nodes;
		return false;														// array wasn't changed, return false
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	private static int MAX_PRINTEXTENT = 50;
	
	
	@Override
	public String toString() {
		int maxIA = Hsp.length > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : Hsp.length;
		StringBuilder sb = new StringBuilder();
		
		sb.append(super.toString());
		
		if (DEBUG_LEVEL < 3) return sb.toString();
		
		sb.append("NSP data:\nHsp row groups:\n");
		for (int i = 0; i < maxIA; i++) {
			
			NSPNode[] bHsp = Hsp[i].array;
			NSPArray bregHsp = Hsp[i];
			if (bHsp != null) {
				sb.append("Hsp(" + i + "): [");
				int maxHsp = bregHsp.nodes > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : bregHsp.nodes;
				for (int j = 0; j < maxHsp; j++)
					sb.append("(" + 	bHsp[j].r + ", " + 
										bHsp[j].c + ", " + 
										String.format("%.3f", bHsp[j].v) + ", " + 
										bHsp[j].offH + ", " + 
										bHsp[j].offV + (j == maxHsp-1 ? ")" : "), "));
				if (bregHsp.nodes > MAX_PRINTEXTENT)
						sb.append(" ... ]\n");
				else	sb.append("]\n");
				
			} else
				sb.append("Hsp(" + i + "): [empty]\n\n");
		}
		
		sb.append("\n\nVsp column groups:\n");
		maxIA = Vsp.length > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : Vsp.length;
		for (int i = 0; i < maxIA; i++) {
			
			NSPNode[] bVsp = Vsp[i].array;
			NSPArray bregVsp = Vsp[i];
			if (bVsp != null) {
				sb.append("Vsp(" + i + "): [");
				int maxVsp = bregVsp.nodes > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : bregVsp.nodes;
				for (int j = 0; j < maxVsp; j++)
					sb.append("(" +		bVsp[j].r + ", " +
										bVsp[j].c + ", " +
										String.format("%.3f", bVsp[j].v) + ", " +
										bVsp[j].offH + ", " +
										bVsp[j].offV + (j == maxVsp-1 ? ")" : "), "));
				
				if (bregVsp.nodes > MAX_PRINTEXTENT)
					sb.append(" ... ]\n");
			else	sb.append("]\n");
			} else
				sb.append("Vsp(" + i + "): [empty]\n\n");
		}

		sb.append("Matrix size: " + (M*N*4) + "\n");
		sb.append("Total NSP structure size: " + new SizeOf(this, false).size() + "\n");
		return sb.toString();
	}
	
	
	public void toFile(int precision) {
		
		super.toFile(precision);
		
		File file = new File(name + "_NspData.txt");
		if (!file.exists()) {
			try {	file.createNewFile();
			} catch (IOException e) { e.printStackTrace(); }
		}
		BufferedWriter bw = null;
		try {		bw = new BufferedWriter(new FileWriter(file));
		} catch (IOException e) { e.printStackTrace(); }

		StringBuilder sb = new StringBuilder();
		
		sb.append("NSP data:\nHsp row groups:\n");
		for (int i = 0; i < Hsp.length; i++) {
			
			NSPNode[] bHsp = Hsp[i].array;
			NSPArray bregHsp = Hsp[i];
			if (bHsp != null) {
				sb.append("Hsp(" + i + "): [");
				for (int j = 0; j < bregHsp.nodes; j++)
					sb.append("(" + 	bHsp[j].r + ", " + 
										bHsp[j].c + ", " + 
										String.format("%.3f", bHsp[j].v) + ", " + 
										bHsp[j].offH + ", " + 
										bHsp[j].offV + (j == bregHsp.nodes-1 ? ")" : "), "));
				sb.append("]\n");			
			} else
				sb.append("Hsp(" + i + "): [empty]\n\n");
		}
		
		sb.append("\n\nVsp column groups:\n");

		for (int i = 0; i < Vsp.length; i++) {
			
			NSPNode[] bVsp = Vsp[i].array;
			NSPArray bregVsp = Vsp[i];
			if (bVsp != null) {
				sb.append("Vsp(" + i + "): [");
				for (int j = 0; j < bregVsp.nodes; j++)
					sb.append("(" +		bVsp[j].r + ", " +
										bVsp[j].c + ", " +
										String.format("%.3f", bVsp[j].v) + ", " +
										bVsp[j].offH + ", " +
										bVsp[j].offV + (j == bregVsp.nodes-1 ? ")" : "), "));
				
				sb.append("]\n");
			} else
				sb.append("Vsp(" + i + "): [empty]\n\n");
		}

		sb.append("Matrix size: " + (M*N*4) + "\n");
		sb.append("Total NSP structure size: " + new SizeOf(this, false).size() + "\n");

		try {
			bw.write(sb.toString());
			bw.flush();
			bw.close();
		} catch (IOException e) { e.printStackTrace(); }
	}
	
	
	
	public void toGraphViz(boolean showGraph, boolean base1, boolean imageView) {
		
		if (M*N > 400) {
			System.out.println("NSPMatrix.toGraphViz(): 400 node limit struck.");
			return;
		}

		int idxBase = base1 ? 1 : 0;
		boolean colorBlack = true;
		StringBuilder sb = new StringBuilder();
		sb.append("digraph G {\n     mclimit=10\n     rankdir=BT\n");
		
		for (int v = 0; v < M; v++) {
			NSPArray aHsp = Hsp[v];
			NSPNode[] bHsp = aHsp.array;
			for (int eU = 0; eU < aHsp.nodes; eU++) {
				
				// ignore links to oneself
				if (bHsp[eU].c == v) continue;
				// see if for current link L -> U there is a reverse link U -> L
				int offV = findHVspNode(Vsp[v].array, 0, Vsp[v].nodes - 1, bHsp[eU].c, -1);
				if (offV >= 0) {
					if (bHsp[eU].c > v)	{
						if (colorBlack) { sb.append("     edge [color=red];\n"); colorBlack = false; }
						sb.append("     " + (v + idxBase) + " [shape=plaintext];\n");
						sb.append("     " + (bHsp[eU].c + idxBase) + " [shape=plaintext];\n");
						sb.append("     " + (v + idxBase) + " -> " + (bHsp[eU].c + idxBase) + " [label=\"LU\"];\n");
					}
				} else if (bHsp[eU].c > v) {
					if (!colorBlack) { sb.append("     edge [color=black];\n"); colorBlack = true; }
					sb.append("     " + (v + idxBase) + " [shape=plaintext];\n");
					sb.append("     " + (v + idxBase) + " -> " + (bHsp[eU].c + idxBase) + " [label=\"U\"];\n");					
				}
				else {
					if (!colorBlack) { sb.append("     edge [color=black];\n"); colorBlack = true; }
					sb.append("     " + (bHsp[eU].c + idxBase) + " [shape=plaintext];\n");
					sb.append("     " + (bHsp[eU].c + idxBase) + " -> " + (v +idxBase) + " [label=\"L\"];\n");
				}
			}
		}
		sb.append("}\n");

		String gwizName = name + "_gviz.txt";
		File file = new File(gwizName);
		if (!file.exists()) {
			try {	file.createNewFile();
			} catch (IOException e) { e.printStackTrace(); }
		}
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(file));
			bw.write(sb.toString());
			bw.flush();
			bw.close();
		} catch (IOException e) { e.printStackTrace(); }
		
		if (!showGraph) return;
		try {
			ProcessBuilder GVizPB;
			if (imageView)
				GVizPB = new ProcessBuilder("C:\\Program Files\\graphviz\\bin\\dot.exe", "-Tpng", gwizName, "-o", name + ".png");
			else
				GVizPB = new ProcessBuilder("C:\\Program Files\\graphviz\\bin\\dot.exe", "-Tps", gwizName, "-o", name + ".ps");
			try { GVizPB.start().waitFor();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			Process viewP;
			if (imageView)
				viewP = new ProcessBuilder("C:\\Program Files (x86)\\FastStone Image Viewer\\FSViewer.exe", name + ".png").start();
			else
				viewP = new ProcessBuilder("C:\\Program Files\\Ghostgum\\gsview\\gsview64.exe", name + ".ps").start();
					
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}