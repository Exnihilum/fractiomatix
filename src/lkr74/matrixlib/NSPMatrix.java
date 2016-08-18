package lkr74.matrixlib;


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//			HELPER CLASSES FOR SPARSE NODE METHODS
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// the sparse node NspNode is referenced by the NSPMatrix both from a row/horizontal aspect, and
// from a column/vertical aspect, from the relevant referencing NspArrays in Hsp & Vsp
class NspNode {
	int r, c, offH, offV;
	double v, iv;
	NspNode(int r, int c, double v, double iv) { this.r=r; this.c=c; this.v=v; this.iv=iv; }
	NspNode(int r, int c, double v, double iv, int offH, int offV) {
		this.r=r; this.c=c; this.v=v; this.iv=iv; this.offH = offH; this.offV = offV; }
	@Override
	protected NspNode clone() { NspNode n = new NspNode(r, c, v, iv); n.offH = offH; n.offV = offV; return n; }
}

// a NspArray is a freestanding sparse datastructure that can be manipulated by relevant methods,
// multiplied, recombined, added and so on, if the method returns an assembled NspArray, then it
// has to be integrated with the arrays of the opposite aspect: Hsp -> update Vsp, Vsp -> update Hsp
class NspArray {
	int nodes, size;
	NspNode[] array;
	NspArray(int nodes, int size, NspNode[] array) { this.nodes = nodes; this.size = size; this.array = array; }
	@Override
	protected NspArray clone() { return new NspArray(nodes, size, array); }
	
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("Nodes: " + nodes + ", size: " + size + ", data:\n[");
		if (array != null) {
			for (int i = 0; i < nodes; i++) sb.append(Matrix.to5chars(array[i].v, false) + (i == nodes-1 ? "]\n[" : " "));
			for (int i = 0; i < nodes; i++) sb.append(Matrix.to5chars(array[i].iv, true) + (i == nodes-1 ? "]\n" : " "));
		} else sb.append("null]\n");
		return sb.toString();
	}
}

public class NSPMatrix extends Matrix {

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			DYNAMICALLY RESIZED SPARSE NODE METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Hsp & Vsp are sparse row/column arrays, preallocated to CSR2_ALLOCBLOCK positions each for fast fill-in
	// at end of list, each array is reallocated on hitting CSR2_ALLOCBLOCK-element bounds with another CSR2_ALLOCBLOCK
	
	int nNsp;
	NspArray[] Hsp;				// horisontally/row aligned sparse arrays pointing to nodes (resembling JA in typical CSR)
	NspArray[] Vsp;				// vertically/column aligned sparse arrays
	NspNode[] pivotNsp;			// fast-access to the diagonal nodes
	
	// identifiers for offsets into a JA2 buffer
	// tBl = offset to total buffer length param., tEc = offset to total element count param., sEc = offset to sorted element count param.
	private static final int CSR2_ALLOCBLOCK = 16;
	private static final int CSR2_DEALLOCBLOCK = CSR2_ALLOCBLOCK * 4;
	
	
	private void initNSPMatrix() { // initialises all NSPMatric arrays from scratch
		nNsp = 0;
		Hsp = new NspArray[M];
		for (int i = 0; i < Hsp.length; i++) Hsp[i] = new NspArray(0, 0, null);
		Vsp = new NspArray[N];
		for (int i = 0; i < Vsp.length; i++) Vsp[i] = new NspArray(0, 0, null);
		pivotNsp = new NspNode[M];
		setNull();
	}
	

	public NSPMatrix(String name, int M, int N) {
		super(name, M, N);								// get skeleton Matrix superclass instance
		initNSPMatrix();
	}

	
	public NSPMatrix(String name, int M, int N, double[] data, double[] idata) {
		super(name, M, N);								// get skeleton Matrix superclass instance
		initNSPMatrix();
		putData(data, idata);							// NSP dual-aspect dynamic list creation
		if (data != null) clearNull();
		if (idata != null) setComplex();
		//bitImage = new BinBitImage(this);
		if (Matrix.DEBUG_LEVEL > 2)  System.out.println(this.toString());		
	}

	public NSPMatrix(String name, int M, int N, Type type) {
		super(name, M, N);								// get skeleton Matrix superclass instance
		initNSPMatrix();
		this.generateData(type, 1);
		if (type != Matrix.Type.Null && type != Matrix.Type.Null_Complex) clearNull();
		//bitImage = new BinBitImage(this);
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
	}

	
	
	@Override
	public NSPMatrix clone() {
		NSPMatrix A = (NSPMatrix)super.clone();					// will first call Matrix superclass clone() method
		if (Hsp != null) {
			A.Hsp = Hsp.clone();
			int i = 0;
			for (NspArray aHVsp : Hsp) {
				A.Hsp[i] = aHVsp.clone();
				for (int j = 0; j < aHVsp.nodes; j++)
					if (aHVsp.array[j] != null) A.Hsp[i].array[j] = aHVsp.array[j].clone();
				i++;
			}
		}
		if (Vsp != null) {
			A.Vsp = Vsp.clone();
			int i = 0;
			for (NspArray aHVsp : Vsp) {
				A.Vsp[i] = aHVsp.clone();
				for (int j = 0; j < aHVsp.nodes; j++)
					if (aHVsp.array[j] != null) A.Vsp[i].array[j] = aHVsp.array[j].clone();
				i++;
			}
		}
		A.pivotNsp = pivotNsp.clone();
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
		return A;
	}

	
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
		for (NspArray aHsp : Hsp) {
			NspNode[] bHsp = aHsp.array;
			if (bHsp != null) {
				for (int i = 0; i < aHsp.nodes; i++) data[rN + bHsp[i].c] = bHsp[i].v;
				if (isComplex()) for (int i = 0; i < aHsp.nodes; i++) idata[rN + bHsp[i].c] = bHsp[i].iv;
			}
			rN += N;
		}
		return dataSet;
	}

	@Override
	public void putDataRef(double[] data, double[] idata) { putData(data, idata); }
	
	// method converts flat array matrix data into CSR2 format
	@Override
	public void putData(double[] data, double[] idata) {
				
		// first loop  inserts elements in proper order int Hsp
		for (int i = 0; i < M; i++) {

			NspNode[] bHsp = Hsp[i].array;								// set up fast access references for JA2, regJA2, nodeJA2, A2 & iA2
			NspArray bregHsp = Hsp[i];
			int iN = i * N, eLeft = bregHsp.nodes - bregHsp.size;		// eLeft = number of free elements left in current JA2 buffer	
			int offsHsp = 0;											// offset into current A2 & JA2 buffer
			
			// iterate over values in the matrix row
			for (int j = iN, jEnd = iN + N, col = 0; j < jEnd; j++, col++) {
				
				// if current JA2 is full (tBl - tEc = 0), add another 16 elements
				if (eLeft <= 0) {				
					int newSize = bregHsp.size + CSR2_ALLOCBLOCK;
					NspNode[] oldbHsp = bHsp;
					bHsp = Hsp[i].array = new NspNode[newSize];					// allocate row-wise reference array Hsp
					eLeft += CSR2_ALLOCBLOCK;
					
					int eTotal = bregHsp.nodes;
					for (int c = 0; c < eTotal; c++) bHsp[c] = oldbHsp[c];		// copy over data into new JA2 buffer		
					bregHsp.size = newSize;										// update length of buffer
				}
				
				// we're inserting only nonzero values
				if (!nearZero(data[j]) || (isComplex() && !nearZero(idata[j]))) {
					bHsp[offsHsp] = new NspNode(i, col, data[j], idata != null ? idata[j] : 0);	// insert column index
					nNsp++;											// global node count incremented
					if (i == col) pivotNsp[i] = bHsp[offsHsp];		// if it's a pivot, put in fast-access array
					
					bHsp[offsHsp].offH = offsHsp;					// store the node's reference offset in Hsp
					offsHsp++;
					bregHsp.nodes++;								// total element count incremented
					eLeft--; 										// decrease number of free elements left in JA2
				}
			}
			if (bregHsp.nodes == 0) Hsp[i].array = null;					// dereference empty rows
		}

		int nodeCnt = 0;											// counts up number of assembled nodes so far
		int[] idxHsp = new int[M];									// stores indexes into JA2 arrays to use in crosslinking
		int[] linkRow = new int[M];
		int[] offsHsp = new int[M];									// holds current offsets from left into every Hsp sparse row

		// iterate over columns, this sets up a left-to-right marching through sparse column data
		for (int c = 0; c < N; c++) {
			
			int toLink = 0;											// counts number of vertical JA2 matches found for crosslinking
			// find all JA2 arrays containing current column c
			for (int r = 0; r < M; r++) {
				NspNode[] bHsp = Hsp[r].array;
				NspArray bregHsp = Hsp[r];
				if (offsHsp[r] >= 0 && bHsp != null && bHsp[offsHsp[r]].c == c) {
					linkRow[toLink] = r;
					idxHsp[toLink++] = offsHsp[r];
					offsHsp[r]++;									// increment the marching offset into this JA2 row
					if (offsHsp[r] >= bregHsp.nodes)				// if we have exhausted current JA2 row...
						offsHsp[r] = -1;							// ...flag it as deactivated
				}
			}
			
			// create the column-wise group
			if (toLink > 0) {
				int groupSize = trimToAllocBlock(toLink + 1);		// allocates in CSR2_ALLOCBLOCK sizes
				//int groupSize = toLink + 1;						// minimal +1 node allocation option
				Vsp[c].array = new NspNode[groupSize];
				Vsp[c].size = groupSize;							// store total buffer length
			}
			Vsp[c].nodes = toLink;									// store node count in Vsp register

			// reference all found nodes for this column into current Vsp array
			int nodeCnt2 = nodeCnt + toLink;
			for (int node = nodeCnt, r = 0; node < nodeCnt2; node++, r++) {
				Vsp[c].array[r] = Hsp[linkRow[r]].array[idxHsp[r]];
				// for convenience, store the node's reference offset in Vsp, this allows to find a pivot quicker
				Vsp[c].array[r].offV = r;	
			}
			
			nodeCnt = nodeCnt2;
		}

	}
	
	
	
	@Override
	public double valueOf(int r, int c) {
		
		if (r < 0 || c < 0 || r >= M || c >= N)
			throw new RuntimeException("NSPMatrix.valueOf(): Invalid matrix coordinates.");
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
	
	
	@Override
	public void valueTo(int r, int c, double v) {
		if (r < 0 || c < 0 || r >= M || c >= N)
			throw new RuntimeException("NSPMatrix.valueTo(): Invalid matrix coordinates.");
		
		if (nearZero(v)) {												// insertion of zero means removal of a value
			int nNodes = Hsp[r].nodes;
			NspNode[] bHsp = Hsp[r].array;
			int offH = findHVspNode(bHsp, 0, nNodes - 1, -1, c);		// locate offset of node in Hsp
			int offV;
			if (offH >= 0) offV = bHsp[offH].offV;						// get node's offset into Vsp if node existed in Hsp
			else return;												// a node that doesn't exist in Hsp will not exist in corresponding Vsp
			removeLocalHVspNode(Hsp[r], offH, 0);						// dereference node from horisontal sparse array
			removeLocalHVspNode(Vsp[c], offV, 1); nNsp--;				// dereference node from vertical sparse array
			return;
		}
		NspNode node = new NspNode(r, c, v, 0); nNsp++;
		int offset = 		insertHVspNode(r, c, 0, node);				// insert value as new node into Hsp
		if (offset >= 0) 	insertHVspNode(r, c, 1, node);				// if node didn't exist already, insert value as new node into Vsp
	}

		

	
	
	public NSPMatrix multiply(NSPMatrix T) {
		if (N != T.M) throw new RuntimeException("NSPMatrix.multiply(): Nonmatching matrix dimensions.");

		// results stored in new NSPMatrix U, which is complex if both in-matrices are complex
		NSPMatrix U = new NSPMatrix("M", M, T.N);
		U.clearNull();
		NspNode node;
		U.name = (Matrix.DEBUG_LEVEL > 1 ? "(" + name + "." + T.name + ")" : "M" + nameCount++);
		int[] offsetV = new int[M];		// array will keep incremental offsets into all the Vsp sparse column buffers

		for (int i = 0; i < M; i++) {
			NspArray aHspU = U.Hsp[i];
			
			for (int j = 0, j1 = 0; j < T.N; j++) {				
				double v = multiplyHVsp(Hsp[i], T.Vsp[j], 0, 1);
				if (!nearZero(v)) {
					NspArray aVspU = U.Vsp[j];
					aHspU.nodes++;
					aVspU.nodes++;
					updateArraySize(aHspU);				// increment target horisontal array and check need for reallocation
					updateArraySize(aVspU);				// increment target vertical array and check need for reallocation
					aHspU.array[j1] = node = new NspNode(i, j, v, 0, j1++, offsetV[j]);
					aVspU.array[offsetV[j]++] = node;
					U.nNsp++;
					if (i == j) U.pivotNsp[i] = node;		// if it's a pivot, insert in fast-access array
				}
			}
		}
		return U;
	}
	



	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			SPARSE NODE MANIPULATION METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	private static int nodeSearch_iterateLim = 8;
	private static int nodeSearch_linearLim = 256;
	
	// method searches for a NspNode in a bounded array, according to row if r >= 0, otherwise according to column c
	// if not found returns negated offset to the place where the node is supposed to be within a sorted order
	static int findHVspNode(NspNode[] bHVsp, int start, int end, int r, int c) {

		if (bHVsp == null) return -1;
		// do ordinary linear search for ranges less than 8
		if (end - start <= nodeSearch_iterateLim) {					// do iterative search if no.of nodes is n <= 8
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
		} else if (end - start <= nodeSearch_linearLim) {			// do linear approx.search if no. of nodes is 8 < n <= 320
			// linear approximation search
			if (r >= 0) {
				if (r < bHVsp[start].r) return -(start + 1);		// base case: supplied row is before the entire Vsp array
				if (bHVsp[end].r < r) return -(end + 2);			// base case: supplied row is after the entire Vsp array

				int rStart = bHVsp[start].r;
				// use intermediate (long) predict to avoid integer overrun on larger numbers of nodes
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

				if (c < bHVsp[predict].c)
						for (int j = predict; j >= start; j--) { 
							if (bHVsp[j].c < c) return -(j+1);
							if (bHVsp[j].c == c) return j;
						}
				else	for (int j = predict; j <= end; j++) { 
							if (c < bHVsp[j].c) return -(j+1);
							if (bHVsp[j].c == c) return j;
						}				
			}			
		} else {
			// make binary search for sought element in Hsp/Vsp sparse array, naturally only sorted elements apply
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

	
	
	// method adds two supplied sparse rows into a new row, the horisontal/vertical aspect chosen individually
	// the new row is returned for integration into a CSR2 scheme
	static NspArray addHVsp(NspArray aHVspA, NspArray aHVspB, double f, int aspectA, int aspectB, int aspectC) {
		
		if (aHVspA.array == null && aHVspB.array == null) return null;
		int nNa = aHVspA.nodes, nNb = aHVspB.nodes;
		int nNc = trimToAllocBlock(nNa + nNb), a, b;
		NspNode[] bHVspC = new NspNode[nNc], bHVspA = aHVspA.array, bHVspB = aHVspB.array;
		int maxN = Integer.MAX_VALUE;
		
		int ic = 0;
		for (int ia = 0, ib = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
			NspNode nA = bHVspA == null ? null : bHVspA[ia];
			NspNode nB = bHVspB == null ? null : bHVspB[ib];
			if (aspectA == 0) {
				if (aspectB == 0) 	{ a = (nA == null ? maxN : nA.c); b = (nB == null ? maxN : nB.c); }
				else				{ a = (nA == null ? maxN : nA.c); b = (nB == null ? maxN : nB.r); }
			} else {
				if (aspectB == 0) 	{ a = (nA == null ? maxN : nA.r); b = (nB == null ? maxN : nB.c); }
				else				{ a = (nA == null ? maxN : nA.r); b = (nB == null ? maxN : nB.r); }
			}
			if (a == b) {
				if (aspectC == 0)
						{ bHVspC[ic] = new NspNode(0, a, nA.v + nB.v * f, nA.iv + nB.iv * f); bHVspC[ic].offH = ic++; }
				else	{ bHVspC[ic] = new NspNode(a, 0, nA.v + nB.v * f, nA.iv + nB.iv * f); bHVspC[ic].offV = ic++; }
				ia++; ib++;
			} else if (a < b) {
				if (aspectC == 0)
						{ bHVspC[ic] = new NspNode(0, a, nA.v, nA.iv); bHVspC[ic].offH = ic++; }
				else	{ bHVspC[ic] = new NspNode(a, 0, nA.v, nA.iv); bHVspC[ic].offV = ic++; }
				ia++;
			} else {
				if (aspectC == 0)
						{ bHVspC[ic] = new NspNode(0, b, nB.v, nB.iv); bHVspC[ic].offH = ic++; }
				else	{ bHVspC[ic] = new NspNode(b, 0, nB.v, nB.iv); bHVspC[ic].offV = ic++; }
				ib++;	
			}
		}
		return new NspArray(ic, nNc, bHVspC);
	}

	
	
	// sparse CSR2 inner vector product, with chosen horisontal/vertical (Hsp/Vsp) multiplying aspect
	static double multiplyHVsp(NspArray aHVspA, NspArray aHVspB, int aspectA, int aspectB) {
		
		NspNode[] bHVspA = aHVspA.array, bHVspB = aHVspB.array;
		if (bHVspA == null || bHVspB == null) return 0;				// base case: one of the arrays is all zeroes
		int nNa = aHVspA.nodes, nNb = aHVspB.nodes;
		double v = 0;
	
		for (int ia = 0, ib = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
			NspNode nA = bHVspA[ia];
			if (nA == null) break;						// if one of nodes is null, then we've finished one of the arrays			
			NspNode nB = bHVspB[ib];
			if (nB == null) break;		
			if (aspectA == 0) {
				if (aspectB == 0) {
					if (nA.c == nB.c) { v += nA.v * nB.v; ia++; ib++; }
					else if (nA.c < nB.c) ia++; else ib++;
				} else {
					if (nA.c == nB.r) { v += nA.v * nB.v; ia++; ib++; }
					else if (nA.c < nB.r) ia++; else ib++;
				}
			} else {
				if (aspectB == 0) {
					if (nA.r == nB.c) { v += nA.v * nB.v; ia++; ib++; }
					else if (nA.r < nB.c) ia++; else ib++;
				} else {
					if (nA.r == nB.r) { v += nA.v * nB.v; ia++; ib++; }
					else if (nA.r < nB.r) ia++; else ib++;
				}				
			}
		}
		return v;
	}

		
	// method inserts node into a Hsp/Vsp reference array, reallocating it if necessary
	// method returns the passed buffer, or the new buffer if reallocation happened
	// the aspect parameter decides way to insert: 0 -> row insert, 1 -> column insert
	// method returns the insertion offset or the found offset negated if node existed
	private int insertHVspNode(int r, int c, int aspect, NspNode node) {

		NspNode[] bHVsp = null;
		NspArray regHVsp = null;
		int offset = 0;
		switch (aspect) {
			case 0:	{	bHVsp = Hsp[r].array; regHVsp = Hsp[r];
						offset = findHVspNode(bHVsp, 0, regHVsp.nodes - 1, -1, c); break; }
			case 1:	{	bHVsp = Vsp[c].array; regHVsp = Vsp[c];
						offset = findHVspNode(bHVsp, 0, regHVsp.nodes - 1, r, -1); break; }
		}
		// if we didn't find element, insert it
		if (offset < 0) {
			offset = -offset-1;
			regHVsp.nodes++;
			if (r == c) pivotNsp[r] = bHVsp[offset];			// if pivot, insert in fast-access array

			// does Hsp/Vsp reference array need allocation/reallocation ?
			if (regHVsp.nodes > regHVsp.size) {
				NspNode[] bHVsp2 = new NspNode[trimToAllocBlock(regHVsp.nodes)];
				if (aspect == 0) Hsp[r].array = bHVsp2; else Vsp[c].array = bHVsp2;
				
				int i = regHVsp.nodes - 1;
				// choose inner offset update type, depending on whether we're processing Hsp or Vsp
				if ((aspect & 1)==0) {
					for (int j = i - 1; j >= offset; i--, j--) { bHVsp2[i] = bHVsp[j]; bHVsp2[i].offH++; }
					node.offH = offset;
				} else {
					for (int j = i - 1; j >= offset; i--, j--) { bHVsp2[i] = bHVsp[j]; bHVsp2[i].offV++; }
					node.offV = offset;
				}
				bHVsp2[offset] = node;										// reference the node
	
				for (i -= 1; i >= 0; i--) bHVsp2[i] = bHVsp[i];				// move the remaining elements
			} else {
				int i = regHVsp.nodes - 1;
				if ((aspect & 1)==0) {
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
	private NspArray insertLocalHVspNode(NspArray aHVsp, int offset, int aspect, NspNode node) {

		NspNode[] bHVsp = aHVsp.array;
		int nNodes = aHVsp.nodes++;
		if (bHVsp[offset].r == bHVsp[offset].c)
			pivotNsp[bHVsp[offset].r] = bHVsp[offset];		// if pivot, insert in fast-access array

		// does Hsp/Vsp reference array need allocation/reallocation ?
		if (nNodes > aHVsp.size) {
			NspNode[] bHVsp2 = aHVsp.array = new NspNode[trimToAllocBlock(nNodes)];
			
			int i = nNodes - 1;
			// choose inner offset update type, depending on whether we're processing Hsp or Vsp
			if ((aspect & 1)==0) {
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
			if ((aspect & 1)==0) {
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

		NspNode[] bHVsp = null;
		NspArray regHVsp = null;
		int offset = 0;
		
		switch (aspect) {
			case 0:	{	bHVsp = Hsp[r].array; regHVsp = Hsp[r];
						offset = findHVspNode(bHVsp, 0, regHVsp.nodes-1, -1, c); break; }
			case 1:	{	bHVsp = Vsp[c].array; regHVsp = Vsp[c];
						offset = findHVspNode(bHVsp, 0, regHVsp.nodes-1, r, -1); break; }
		}
		
		// if we found element, remove it
		if (offset >= 0) {
			if (r == c) pivotNsp[r] = null;								// if it's a pivot, remove from fast-access array
			
			int nNodes = --regHVsp.nodes;
			if (nNodes == 0) {
				regHVsp.size = 0;
				if (aspect == 0) Hsp[r].array = null; else Vsp[c].array = null;		// last node in array removed, destroy this array
				return offset;		
			}
			
			// does Hsp/Vsp reference array need deallocation (over 63 nodes been removed)?
			if (regHVsp.size - nNodes >= 64) {
				NspNode[] bHVsp2 = new NspNode[regHVsp.size - 64];
				
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
	
	
	private void removeLocalHVspNode(NspArray aHVsp, int offset, int aspect) {
		
		NspNode[] bHVsp = aHVsp.array;
		if (bHVsp[offset].r == bHVsp[offset].c)
			pivotNsp[bHVsp[offset].r] = null;						// if it's a pivot, remove from fast-access array
		int nNodes = --aHVsp.nodes;
		if (nNodes == 0) { aHVsp.size = 0; aHVsp.array = null; return; }	// last node in array removed, destroy array
		
		// does Hsp/Vsp reference array need deallocation (over 63 nodes been removed)?
		if (aHVsp.size - nNodes >= CSR2_DEALLOCBLOCK) {
			NspNode[] bHVsp2 = new NspNode[aHVsp.size - CSR2_DEALLOCBLOCK];
			
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


	
	// swaps two sparse arrays of index a & b, Hsp (rows) if aspect = 0, Vsp (columns) if aspect = 1
	// method marches through both swapped sparse arrays, if there are values in both positions then they're
	// trivially swapped, if a value is swapped with a zero, then the value is moved along the other aspect
	// to it's swapped position, at the end the two array references in the main aspect are simply swapped
	public void swapHVspArrays(int a, int b, int aspect) {
		
		if (a == b) return;
		NspNode[] bHVspA = null, bHVspB = null;
		NspNode nA, nB, nTemp = null;
		int nNa = 0, nNb = 0, temp;
		
		switch (aspect) { 
			// switch two rows (Hsp aspect)
			case 0:
				bHVspA = Hsp[a].array; bHVspB = Hsp[b].array; nNa = Hsp[a].nodes; nNb = Hsp[b].nodes;
				
				for (int ia = 0, ib = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
					
					nA = bHVspA == null ? null : bHVspA[ia];
					nB = bHVspB == null ? null : bHVspB[ib];
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
						nA.r = b;											// change it's row to the other one
						if (nA.r == nAc) pivotNsp[nAc] = nA;				// see if A ended up as a pivot
						NspNode[] bVsp = Vsp[nAc].array;
						if (a < b) {										// do the switch within higher part of array b
							for (int j = nA.offV, k = j + 1, kEnd = Vsp[nAc].nodes; j < kEnd; j++, k++) {
								if (k == kEnd) { nA.offV = bVsp[j].offV; bVsp[j] = nA; break; }
								if (bVsp[k].r > b) {						// if next Vsp row index is higher than array b's
									nA.offV = bVsp[j].offV + 1;
									bVsp[j] = nA;							// then insert the swapped node at current place in sparse array
									break;
								}
								bVsp[j] = bVsp[k];							// otherwise shift nodes backwards
								bVsp[j].offV--;								// decrease their offsets into Vsp
							}
						} else {
							for (int j = nA.offV, k = j - 1; j >= 0; j--, k--) {
								if (k < 0) { nA.offV = bVsp[j].offV; bVsp[j] = nA; break; }
								if (bVsp[k].r < b) {						// if next Vsp row index is lower than array b's
									nA.offV = bVsp[j].offV - 1;
									bVsp[j] = nA;							// then insert the swapped node at current place in sparse array
									break;
								}
								bVsp[j] = bVsp[k];							// otherwise shift nodes forwards
								bVsp[j].offV++;								// increase their offsets into Vsp array
							}
						}
						ia++;
					} else {												// seems like nA.c > nB.c, do the same routine, but for nB
						nB.r = a;											// change it's row to the other one
						if (nB.r == nBc) pivotNsp[nBc] = nB;				// see if B ended up as a pivot
						NspNode[] bVsp = Vsp[nBc].array;
						if (a < b) {										// do the switch within lower part of array b
							for (int j = nB.offV, k = j - 1; j >= 0; j--, k--) {
								if (k < 0) { nB.offV = bVsp[j].offV; bVsp[j] = nB; break; }
								if (bVsp[k].r < a) {						// if next Vsp row index is lower than array b's
									nB.offV = bVsp[j].offV - 1;
									bVsp[j] = nB;							// then insert the swapped node at current place in sparse array
									break;
								}
								bVsp[j] = bVsp[k];							// otherwise shift nodes forwards
								bVsp[j].offV++;								// increase their offsets into Vsp array
							}
						} else {
							for (int j = nB.offV, k = j + 1, kEnd = Vsp[nBc].nodes; j < kEnd; j++, k++) {
								if (k == kEnd) { nB.offV = bVsp[j].offV; bVsp[j] = nB; break; }
								if (bVsp[k].r > a) {						// if next Vsp row index is higher than array b's
									nB.offV = bVsp[j].offV + 1;
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
				NspArray aTemp = Hsp[a]; Hsp[a] = Hsp[b]; Hsp[b] = aTemp;				// finally, exchange the rows in Hsp
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
						nA.c = b;											// change it's column to the other one
						if (nA.c == nAr) pivotNsp[nAr] = nA;				// see if A ended up as a pivot
						NspNode[] bHsp = Hsp[nAr].array;
						if (a < b) {										// do the switch within right part (relative A) of array B
							for (int j = nA.offH, k = j + 1, kEnd = Hsp[nAr].nodes; j < kEnd; j++, k++) {
								if (k == kEnd) { nA.offH = bHsp[j].offH; bHsp[j] = nA; break; }
								if (bHsp[k].c > b) {			// if next Hsp column index is higher than array B's
									nA.offH = bHsp[j].offH + 1;
									bHsp[j] = nA;							// then insert the swapped node at current place is sparse array
									break;
								}
								bHsp[j] = bHsp[k];							// otherwise shift nodes backwards
								bHsp[j].offH--;								// decrease their offsets into Hsp array
							}
						} else {
							for (int j = nA.offH, k = j - 1; j >= 0; j--, k--) {
								if (k < 0) { nA.offH = bHsp[j].offH; bHsp[j] = nA; break; }
								if (bHsp[k].c < b) {						// if next Hsp column index is lower than array B's
									nA.offH = bHsp[j].offH - 1;
									bHsp[j] = nA;							// then insert the swapped node at current place is sparse array
									break;
								}
								bHsp[j] = bHsp[k];							// otherwise shift nodes forwards
								bHsp[j].offH++;								// increase their offsets into Hsp array
							}
						}
						ia++;
					} else {												// seems like nA.c > nB.c, do the same routine, but for nB
						nB.c = a;											// change it's column  to the other one
						if (nB.c == nBr) pivotNsp[nBr] = nB;				// see if B ended up as a pivot
						NspNode[] bHsp = Hsp[nB.r].array;
						if (a < b) {										// do the switch within left part (relative A) of array B
							for (int j = nB.offH, k = j - 1; j >= 0; j--, k--) {
								if (k < 0) { nB.offH = bHsp[j].offH; bHsp[j] = nB; break; }
								if (bHsp[k].c < a) {				// if next Hsp column index is lower than array B's
									nB.offH = bHsp[j].offH - 1;
									bHsp[j] = nB;							// then insert the swapped node at current place is sparse array
									break;
								}
								bHsp[j] = bHsp[k];							// otherwise shift nodes forwards
								bHsp[j].offH++;								// increase their offsets into Vsp
							}
						} else {
							for (int j = nB.offH, k = j + 1, kEnd = Hsp[nBr].nodes; j < kEnd; j++, k++) {
								if (k == kEnd) { nB.offH = bHsp[j].offH; bHsp[j] = nB; break; }
								if (bHsp[k].c > a) {						// if next Hsp column index is higher than array B's
									nB.offH = bHsp[j].offH + 1;
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

	
	
	// method chooses the smallest aspect (Hsp vs Vsp) to search for an element
	// if aspect = 0, returns offset into Hsp, if aspect = 1 returns offset into Vsp
	// if element wasn't found, returns negative value with 30th bit set if it's from Vsp, otherwise the bit is clear
	private int locateNspNode(int r, int c, int aspect) {
		
		if (r == c && pivotNsp[r] != null)
			return aspect == 0 ? pivotNsp[r].offH : pivotNsp[r].offV;	// if a pivot if sought, there is a fast-access buffer
		
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
	
	

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			INLINE/HELPER METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	private static int trimToAllocBlock(int v) {
		if ((v & 0xF) == 0) return (v & 0xFFFFFFF0); else return (v & 0xFFFFFFF0) + CSR2_ALLOCBLOCK;
	}
	
	
	private static boolean updateArraySize(NspArray aHVsp) {
		int nodes = aHVsp.nodes;
		// check if array needs increasing
		if (nodes >= aHVsp.size) {
			NspNode[] newArray = new NspNode[aHVsp.size = trimToAllocBlock(aHVsp.nodes)], array = aHVsp.array;
			if (array != null)
				for (int i = 0; i < nodes; i++) newArray[i] = array[i];
			aHVsp.array = newArray;
			return true;
		}
		// check if array needs shrinking
		if (aHVsp.size - nodes <= CSR2_DEALLOCBLOCK) {
			aHVsp.size = trimToAllocBlock(aHVsp.nodes);
			if (aHVsp.size == 0) { aHVsp.array = null; return true; }		// if last element removed, destroy node array
			NspNode[] newArray = new NspNode[aHVsp.size], array = aHVsp.array;
			if (array != null)
				for (int i = 0; i < nodes; i++) newArray[i] = array[i];
			aHVsp.array = newArray;
			return true;			
		}
		return false;														// array wasn't changed, return false
	}

	// returns total approximate size of the CSR2 structure, a reference counted as 2 ints
	public int sizeOf() {
		final int refSize = 8;
		long nnSize = 16; // ObjectSize.sizeOf(new NspNode(0,0,0,0));
		long naSize = 24; // ObjectSize.sizeOf(new NspArray(0, 0, null));
		
		int totalSize = Hsp.length * refSize + Vsp.length * refSize + pivotNsp.length * refSize;
		for (NspArray l : Hsp)
			if (l != null) {
				totalSize += naSize;											// size of NspArray object
				for (int i = 0; i < l.nodes; i++)
					if (l.array[i] != null) totalSize += nnSize;	// size of NspNode object
			}
		for (NspArray l : Vsp)
			if (l != null) {
				totalSize += naSize;											// size of NspArray object
				for (int i = 0; i < l.nodes; i++)
					if (l.array[i] != null) totalSize += nnSize;	// size of NspNode object
			}

		return totalSize;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	private static int MAX_PRINTEXTENT = 50;
	
	
	@Override
	public String toString() {
		int maxIA = Hsp.length > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : Hsp.length;
		StringBuffer sb = new StringBuffer();
		
		sb.append(super.toString());
		
		sb.append("CSR2 data:\nHsp row groups:\n");
		for (int i = 0; i < maxIA; i++) {
			
			NspNode[] bHsp = Hsp[i].array;
			NspArray bregHsp = Hsp[i];
			if (bHsp != null) {
				sb.append("Hsp(" + i + "): Sorted: [");
				int maxHsp = bregHsp.nodes > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : bregHsp.nodes;
				for (int j = 0; j < maxHsp; j++)
					sb.append("(" + 	bHsp[j].r + ", " + 
										bHsp[j].c + ", " + 
										bHsp[j].v + ", " + 
										bHsp[j].offH + ", " + 
										bHsp[j].offV + (j == maxHsp-1 ? ")]\n" : "), "));
				
			} else
				sb.append("Hsp(" + i + "): [empty]\n\n");
		}
		
		sb.append("\n\nVsp column groups::\n");
		maxIA = Vsp.length > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : Vsp.length;
		for (int i = 0; i < maxIA; i++) {
			
			NspNode[] bVsp = Vsp[i].array;
			NspArray bregVsp = Vsp[i];
			if (bVsp != null) {
				sb.append("Vsp(" + i + "): Sorted: [");
				int maxVsp = bregVsp.nodes > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : bregVsp.nodes;
				for (int j = 0; j < maxVsp; j++)
					sb.append("(" +		bVsp[j].r + ", " +
										bVsp[j].c + ", " +
										bVsp[j].v + ", " +
										bVsp[j].offH + ", " +
										bVsp[j].offV + (j == maxVsp-1 ? ")]\n" : "), "));
				
			} else
				sb.append("Vsp(" + i + "): [empty]\n\n");
		}

		sb.append("Matrix size: " + (M*N*4) + "\n");
		sb.append("Total NSP structure size: " + sizeOf() + "\n");
		return sb.toString();
	}
}