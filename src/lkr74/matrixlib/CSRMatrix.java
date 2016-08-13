package lkr74.matrixlib;

import java.util.Arrays;

public class CSRMatrix extends Matrix {
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			INSTANCE-LEVEL VALUES, CSR SPECIFIC
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// this is the standard flat array style of CSR, simple to analyse but not optimal for searching and inserting
	private int nA, nIA;
	private double[] A, iA;
	private int[] IA, JA;
		
	// IAop & JAop are situational, permamently-allocated/reallocated operational buffers
	private static double[] Aop;
	private static int[] JAop;
	
	{
		Aop = new double[8*8];	// on class load, preallocate situational buffers to default 8x8 matrix size
		JAop = new int[8*8];
	}
	
	public CSRMatrix(String name, int M, int N) { super(name, M, N); nA = nIA = 0; }
	
	public CSRMatrix(String name, int M, int N, double[] data, double[] idata) {
		super(name, M, N);								// get skeleton Matrix superclass instance
		putData(data, idata);							// standard CSR flat list creation
		if (data != null) clearNull();
		bitImage = new BinBitImage(this);
		if (Matrix.DEBUG_LEVEL > 2) {
			System.out.println(this.toString());
			System.out.println(this.toStringCSR2());			
		}
	}

	public CSRMatrix(String name, int M, int N, Type type) {
		super(name, M, N);								// get skeleton Matrix superclass instance
		this.generateData(type, 1);
		if (type != Matrix.Type.Null && type != Matrix.Type.Null_Complex) clearNull();
		bitImage = new BinBitImage(this);
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
	}

	
	
	public CSRMatrix clone() {
		CSRMatrix A;
		A = (CSRMatrix)super.clone();					// will first call Matrix superclass clone() method
		if (this.A != null)		A.A = this.A.clone();
		if (this.iA != null)	A.iA = this.iA.clone();
		if (this.IA != null)	A.IA = IA.clone();
		if (this.JA != null)	A.JA = JA.clone();
		A.nA = nA;
		A.nIA = nIA;
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
		return A;
	}
	
	
	public static CSRMatrix convert(Matrix A) {
		return new CSRMatrix(A.name, A.M, A.N, A.data, A.idata);
	}
	

	public void zero() {
		A = iA = null;
		for (int i = 0; i < nIA; i++) IA[i + 1] = 0;
		JA = null;
		nA = 0;
		bitImage.zero();
		setNull();
	}
	
	
	// superclass operations operate in RCMatrix space, CSR data needs decoding
	@Override
	public double[][] getDataRef() { return getData(); }

	@Override
	public double[][] getData() {	

		double[][] dataSet = new double[2][];
		double[] data = dataSet[0] = new double[M * N], idata = null;
		if (nA == 0) return dataSet;								// there was no data, return empty matrix
		if (iA != null) idata = dataSet[1] = new double[M * N];		// we also have imaginary part, allocate for returning it
		
		for (int r = 0; r < M; r++) {
			int rN = r * N, iar1 = IA[r + 1], nrv = iar1 - IA[r];
			for (int i = iar1 - 1; nrv > 0; nrv--, i--)
				data[rN + JA[i]] = A[i];
			if (iA != null)
				for (int i = iar1 - 1; nrv > 0; nrv--, i--)
					idata[rN + JA[i]] = iA[i];	
		}
		return dataSet;
	}

	// incoming row-column data into CSRMatrix must be converted to native format
	@Override
	public void putDataRef(double[] data, double[] idata) { if (data != null) putDataM(data, idata); }
	@Override
	public void putData(double[] data, double[] idata) {
		if (data != null) { 
			putDataM(data, idata);
			putDataCSR2(data, idata);
		}
	}
	
	// TODO: optimise this algorithm
	public void putDataM(double[] data, double[] idata) {

		if (M < 1 || N < 1)
			throw new RuntimeException("CSRMatrix.putData(): Invalid matrix size.");

		nIA = M + 1;
		IA = new int[nIA];
		IA[0] = 0;
		// iterate through all rows
		for (int j = 0; j < M; j++) {
			int jN = j * N;
			// get sum of all nonzeros in current row
			int c = 0;
			if (idata != null) {
				for (int i = jN, iEnd = jN + N; i < iEnd; i++)
					if (!nearZero(data[i]) || !nearZero(idata[i])) c++;
			} else {
				for (int i = jN, iEnd = jN + N; i < iEnd; i++)
					if (!nearZero(data[i])) c++;
			}
			// assemble compact indexing column IA
			IA[j + 1] = IA[j] + c;
		}

		// gather up all non-zero matrix members in A/iA and column indexes of all non-zero matrix members in JA
		nA = IA[nIA - 1];
		JA = new int[nA];
		if (data != null)	A = new double[nA];
		if (idata != null)	iA = new double[nA];
		
		for (int j = 0, c1 = 0; j < M; j++) {
			// for complex number case, any of real or imaginary = nonzero will store a value doublet
			if (idata != null)
				for (int i = 0, jNi = j * N; i < N; i++, jNi++) {
					if (!nearZero(data[jNi]) || !nearZero(idata[jNi])) { A[c1] = data[jNi]; iA[c1] = idata[jNi]; JA[c1++] = i; }
				}
			else
				for (int i = 0, jNi = j * N; i < N; i++, jNi++) {
					if (!nearZero(data[jNi])) { A[c1] = data[jNi]; JA[c1++] = i; }
				}
		}
	}
	
	
	

	// get any subset of the CSR dataset or an expanded/contracted dataset, converted into RC format
	// Mi = top boundary, Ni = left boundary, Mo = bottom boundary, No = right boundary
	@Override
	public CSRMatrix rescale(int Mi, int Ni, int Mo, int No, boolean doBitImage) {
		
		if (Mi > Mo || Ni > No) throw new RuntimeException("CSRMatrix.getData(): Invalid data subset size.");
		int newM = Mo - Mi, newN = No - Ni;
		
		String newname;
		if (DEBUG_LEVEL > 1)	newname = name + "(r:" + Mi + "," + Ni + "," + "Mo" + "," + No + ")";
		else					newname = name + "(r)";
		CSRMatrix A = new CSRMatrix(newname, newM, newN, Matrix.Type.Null);
		double[] newdata = A.data;

		// nr & nc are indexes into the new dataset
		for (int nr = Mi, r = 0; nr < Mo; nr++, r++) {
			// store a row only if indexing stays inside row bounds of original matrix
			if (nr >= 0 && nr < M ) {
				int nrofs = r * newN, iN = IA[nr + 1];
				for (int nc = Ni, c = 0; nc < No; nc++, c++) {
					// store a value only if indexing stays inside column bounds of original matrix
					if (nc >= 0 && nc < N) {
						// iterate through row's sparse column indices, check if sought column exists -> a sparse value is stored
						for (int i = IA[nr]; i < iN; i++)
							if (JA[i] == nc)	newdata[nrofs + c] = this.A[i];
					}
				}
			}
		}
		A.putData(newdata, null);
		if (Matrix.DEBUG_LEVEL > 1) System.out.println(this.toString());
		if (doBitImage) A.bitImage = new BinBitImage(A);
		return A;
	}
	
	
	// retrieves value of a row & column index
	@Override
	public double valueOf(int r, int c) {
		if (r < 0 || c < 0 || r >= M || c >= N)
			throw new RuntimeException("CSRMatrix.valueOf(): Invalid matrix coordinates.");

		// iterate through row's sparse column indices, check if sought column exists -> a sparse value is stored
		int iN = IA[r + 1]; 
		for (int i = IA[r]; i < iN; i++)
			if (JA[i] == c)	return A[i];
		// no sparse value found, return zero value
		return 0;
	}

	
	
	// inserts value v at row & column index, CSR data insertion
	@Override
	public void valueTo(int r, int c, double v) {
		if (r < 0 || c < 0 || r >= M || c >= N)
			throw new RuntimeException("CSRMatrix.valueTo(): Invalid matrix coordinates.");
		
		double[] A = this.A;
		int[] JA = this.JA;

		// iterate through row's sparse column indices, check if sought column already exists
		int iN = this.IA[r + 1], i_cpos = iN, cmax = c; 
		for (int i = this.IA[r]; i < iN; i++) {
			// it existed, change value and quit
			if (JA[i] == c)	{
				if (nearZero(v)) {				// we're zeroing an existing value, so it must be removed from CSR tables
					int lenA = this.IA[M] - 1;
					for (int j = i; j < lenA ; j++) { A[j] = this.A[j + 1]; JA[j] = this.JA[j + 1]; }
					for (int j = r + 1; j <= M; i++) this.IA[j]--;
				} else {
					A[i] = v;
					clearNull();
				}
				return;
			}
			if (cmax < JA[i]) { cmax = JA[i]; i_cpos = i; }
		}
		
		if (nearZero(v)) return;				// inserting a zero changes nothing for CSR
		clearNull();

		// didn't find column, A & JA must be rescaled, new value(A)-column(JA) pair must be inserted
		// check if A & JA need to be rescaled
		int lenA = this.IA[M];
		int lenAmax = M * M;
		// unexpected catch: CSR table seem to contain all values possible, but we failed finding the sought row & column?
		if (lenA == lenAmax)
			throw new RuntimeException("CSRMatrix.valueTo(): unexpected error.");
		
		if (this.A.length <= lenA) {
			// rescale with some extra elements to optimise for next time A & JA has to be rescaled
			int lenA2 = lenA + (lenAmax - lenA)/3;
			if (lenA2 * 100 / lenAmax < 10) lenA2 = lenAmax;	// if we're only 10% from max size, rescale to max
			if (lenA2 > lenA) {
				A = new double[lenA2];
				JA = new int[lenA2];
			}
		}
		// shift all elements past i_cpos one step forward
		for (int i = lenA - 1; i >= i_cpos ; i--) { A[i + 1] = this.A[i]; JA[i + 1] = this.JA[i]; }
		// insert the new column and value
		A[i_cpos] = v;
		JA[i_cpos] = c;
		// copy remaining elements only if we had to expand A & JA
		if (this.A != A) {
			// copy remaining elements before i_cpos to the rescaled A & JA
			for (int i = i_cpos - 1; i >= 0; i--) { A[i] = this.A[i]; JA[i] = this.JA[i]; }
			this.A = A;
			this.JA = JA;
		}
		// increment IA by 1 between current row & end row, to note increase with 1 element
		for (int i = r + 1; i <= M; i++) this.IA[i]++;

		bitImage.setBit(r, c);
		if (Matrix.DEBUG_LEVEL > 2) System.out.println(this.toString());
	}
	
	
	
	@Override
	public CSRMatrix identity(int S, double v) {
		CSRMatrix I = new CSRMatrix("I", S, S, Type.Identity);
		for (int i = 0; i < I.A.length; i++) I.A[i] *= v;
		return I;
	}
	
	
	
	// returns a CSR matrix with the specified row & column eliminated
	public CSRMatrix eliminateRowColumn(int r, int c, boolean makeBitImage) {

		if (r < 0 || r > M - 1 || c < 0 || c > N - 1)
			throw new RuntimeException("CSRMatrix.eliminateRowColumn(): row or column out of bounds.");
		if (M < 2 || N < 2)
			throw new RuntimeException("CSRMatrix.eliminateRowColumn(): Invalid matrix size.");

		// sA = elimination offset start of A & IA, eA = elimination offset end
		int nrIA = IA.length - 1, sA = IA[r], eA = IA[r + 1], cA = eA - sA;
		// clone a new matrix of same size, don't waste time on recalculating the new reduced size
		// as this operation tends to be heavily used in operations like matrix inverses
		
		String newname;
		if (DEBUG_LEVEL > 1) 	newname = new String(name + "(M-" + r + ",N-" + c);
		else 					newname = name + nameCount++;
		CSRMatrix A = new CSRMatrix(newname, M - 1, N - 1);
		A.A = new double[this.A.length - cA];
		A.IA = new int[nrIA];
		A.JA = new int[this.A.length - cA];
		A.nA = nA - cA;
		A.nIA = nrIA;

		// start by eliminating the row, which will make a column elimination more efficient
		// first copy all offsets before the eliminated row
		r++;
		for (int i = 1; i < r; i++) A.IA[i] = IA[i];
		// copy all rows past eliminated row backwards-shifted and subtracted off the row's constituent elements
		for (int i = r; i < nrIA; i++) A.IA[i] = IA[i + 1] - cA;
		// copy all values coming before the eliminated ones
		int iC = 0;
		for (int i = 0; i < sA; i++) { 
			if (JA[i] == c) {
				// if an eliminated column index found, skip over it and do corresponding decrementing of IA
				for (int j = nrIA - 1; A.IA[j] > i; j--) A.IA[j]--;
				continue;
			}
			A.A[iC] = this.A[i];
			A.JA[iC++] = JA[i] > c ? JA[i] - 1 : JA[i];
		}
		// copy all values after the eliminated ones
		for (int i = eA; i < nA; i++) {
			if (JA[i] == c) {
				for (int j = nrIA - 1; A.IA[j] >= i; j--) A.IA[j]--;
				continue;
			}
			A.A[iC] = this.A[i];
			A.JA[iC++] = JA[i] > c ? JA[i] - 1 : JA[i];
		}
		
		if (Matrix.DEBUG_LEVEL > 1) System.out.println(this.toString());
		// methods that eliminate rows & columns as intermediary states will not need to reconstitute the bitImage
		if (makeBitImage) A.bitImage = new BinBitImage(A);
		return A;
	}
	
	
	
	// transposes a CSR matrix
	@Override
	public Matrix transpose(boolean copy) {
		
		if (M < 1 || N < 1) throw new RuntimeException("CSRMatrix.transpose(): Invalid matrix size.");

		CSRMatrix T = this;
		if (copy) T = new CSRMatrix(name + "^T", N, M);
		int temp = T.M; T.M = T.N; T.N = temp;		// reverse M & N
		
		// transposed matrix will have exactly the same CSR register lengths
		int[] TIA = new int[T.M + 1], offsets = new int[T.M + 1], ncolidx = new int[nA];
		int[] TJA = new int[nA];
		double[] TA = new double[nA];
		nIA = T.M + 1;
		
		// build up the transpose-IA register by summation of counts of elements/column
		// also note down the row for every element to make the rows quickly accessible in next step
		for (int i = 0, n = 0, acc = 0; i < nA; i++, acc++) {
			while (acc == IA[n])  n++;
			TIA[JA[i] + 1]++;
			ncolidx[i] = n - 1;
		}
		// cover up eventual "offset holes" present in the transpose-IA register
		for (int i = 1; i < nIA; i++)
			TIA[i] = offsets[i] = TIA[i] + TIA[i - 1];
		// reorder A & JA into transpose-A & transpose-JA, while updating offsets for sequential value insertions
		for (int i = 0, vJA; i < nA; i++) {
			vJA = JA[i];
			TJA[offsets[vJA]] = ncolidx[i];
			TA[offsets[vJA]] = A[i];
			offsets[vJA]++;
		}
		
		// update registers to transposed registers
		T.A = TA;
		T.IA = TIA;
		T.JA = TJA;
		
		if (Matrix.DEBUG_LEVEL > 1) System.out.println(T.toString());
		T.bitImage.make();					// redo bitImage
		return T;
	}
	
	
	
	// compares two matrices CSR style
	public static boolean equal(CSRMatrix S, CSRMatrix T) {
		if (T.M != S.M || T.N != S.N)
			throw new RuntimeException("CSRMatrix.equal(): Nonmatching matrix dimensions.");

		// fast-compare bitImages for non-matching non-zero fields
		if (!S.bitImage.equals(T.bitImage)) return false;
		
		int nextIA_S = 0, nextIA_T = 0;
		double diff;
		
		for (int i = 0; i < S.M; i++) {
			// set up parallel marching through both rows of S & T
			int jS = nextIA_S, jT = nextIA_T;
			nextIA_S = S.IA[i + 1];
			nextIA_T = T.IA[i + 1];
			while (jS < nextIA_S && jT < nextIA_T) {
				// if one of the S/T column indexes is larger, advance in the opposite row
				int vJAS = S.JA[jS], vJAT = T.JA[jT];
				if (vJAS < vJAT)		jS++;
				else if (vJAS > vJAT)	jT++;
				// matching columns found, compare
				else {
					diff = S.A[jS++] - T.A[jT++];
					if (!nearZero(diff)) return false;
				}
			}
			// we've run out of values either in (S or T), take care of remaining ones in (S or T)
			while(jS < nextIA_S) jS++;
			while(jT < nextIA_T) jT++;
		}
		return true;
	}
	
	
	
	// swaps two rows in a CSR matrix
	public void swap(int r1, int r2) {

		if (r1 < 0 || r2 < 0 || r1 >= M || r2 >= M)
			throw new RuntimeException("CSRMatrix.swap(): Invalid matrix rows.");
		if (r1 == r2) return;
		// change row1 to be a lower index than row2, to simplify algorithmics
		if (r1 > r2) { int t = r1; r1 = r2; r2 = t; }
		
		int i = IA[r2], l1 = IA[r1 + 1] - IA[r1], l2 = IA[r2 + 1] - i, c = 0;
		
		// at least double the situational buffers if they're too short
		if (Aop.length < IA.length) {
			int newsize = Aop.length * 2 < IA.length ? IA.length : Aop.length * 2;
			Aop = new double[newsize];
			JAop = new int[newsize];
		}

		// shift over row2 into correct swapped position in operational buffers JAop & Aop
		for (int j = 0; j < l2; j++) { Aop[c] = A[i]; JAop[c++] = JA[i++]; }
		// shift over everything inbetween row1 & row2 into correct positions within JAop & Aop
		for (i = IA[r1 + 1]; i < IA[r2]; i++) { Aop[c] = A[i]; JAop[c++] = JA[i]; }
		// shift over row1 into correct swapped position in operational buffers JAop & Aop
		i = IA[r1];
		for (int j = 0; j < l1; j++) { Aop[c] = A[i]; JAop[c++] = JA[i++]; }
		// write back the swapped data into correct offsets of A & JA
		for (int j = 0; j < c; j++) { A[IA[r1] + j] = Aop[j]; JA[IA[r1] + j] = JAop[j]; }

		// readjust offsets within IA
		int dl = l2 - l1;
		for (i = r1 + 1; i <= r2; i++) IA[i] += dl;
		
		if (Matrix.DEBUG_LEVEL > 1) System.out.println(this.toString());
		bitImage.swapRows(r1, r2);	// swap bitImage rows
	}

	
	
	public CSRMatrix add(double v, boolean copy)		{ return addSub(this, false, copy, v, true); }
	public CSRMatrix subtract(double v, boolean copy)		{ return addSub(this, true, copy, -v, true); }
	public CSRMatrix add(CSRMatrix T, boolean copy)		{ return addSub(T, false, copy, 0, false); }
	public CSRMatrix subtract(CSRMatrix T, boolean copy)	{ return addSub(T, true, copy, 0, false); }

	// sums two matrices CSR style
	public CSRMatrix addSub(CSRMatrix T, boolean subtract, boolean copy, double scl, boolean addScalar) {
		if (T.M != M || T.N != N)
			throw new RuntimeException("CSRMatrix.addSub(): Nonmatching matrix dimensions.");

		CSRMatrix S = this;
		if (copy) {
			S = this.clone();
			S.name = Matrix.DEBUG_LEVEL > 1 ? "(" + name + (subtract?"-":"+") + T.name + ")" : (subtract?"D":"S") + nameCount++;
		}
		
		if (addScalar && nearZero(scl)) return S;			// in case we're adding a near-zero scalar, nothing more needs to be done
		if (this.isNull() || S.isNull()) return S;			// quick check - are we adding null matrices?
		clearNull();

		// conservatively assume no. of sparse elements in new matrix will be at maximum
		// the sum of sparse elements from the summed matrices
		// in case of adding a scalar, the buffers will effectively be M*N sized
		int lenA = addScalar ? M * N : S.A.length + T.A.length;
		double[] A = new double[lenA];
		int[] IA = new int[S.M + 1];
		int[] JA = new int[lenA];
		
		int accIA = 0;
		int nextIA_S = 0, nextIA_T = 0;
		for (int i = 0; i < S.M; i++) {
			// set up parallel marching through both added rows, add values if column indexes found for both values
			// otherwise, insert the larger value from (S or T)
			int jS = nextIA_S, jT = nextIA_T, vScl = 0;
			nextIA_S = S.IA[i + 1];

			// this code takes care of adding a scalar to a CSR matrix - effectively at this
			// stage one should convert to standard matrix form since the sparsity is destroyed
			if (addScalar)
				while (vScl < N) {
					int vJAS = S.JA[jS];
					if (vJAS != vScl) A[accIA] = scl;
					else A[accIA] = S.A[jS++] + scl;
					JA[accIA++] = vScl++;
				}
			// normal matrix + matrix adding goes on in this segment
			// both marchers keep track of number of values left to read (nextIA are the bounds) for respective marcher
			else {
				nextIA_T = T.IA[i + 1];
				while (jS < nextIA_S && jT < nextIA_T) {		
					int vJAS = S.JA[jS], vJAT = T.JA[jT];
					if (vJAS < vJAT)		{ A[accIA] = S.A[jS++]; JA[accIA++] = vJAS; }
					else if (vJAS > vJAT)	{ A[accIA] = T.A[jT++]; JA[accIA++] = vJAT; }
					// the adding/subtracting happens here
					else {
						A[accIA] = S.A[jS++] + (subtract ? -T.A[jT++] : T.A[jT++]);
						JA[accIA++] = vJAS;
					}
				}
				// we've run out of values either in (S or T), take care of remaining ones in (S or T)
				while(jS < nextIA_S) { A[accIA] = S.A[jS]; JA[accIA++] = S.JA[jS++]; }
				while(jT < nextIA_T) { A[accIA] = T.A[jT]; JA[accIA++] = T.JA[jT++]; }
			}
			IA[i + 1] = accIA;
		}
		
		S.A = A;
		S.IA = IA;
		S.JA = JA;
		if (Matrix.DEBUG_LEVEL > 1) System.out.println(S.toString());
		
		// quick & dirty bitImage copy: assume OR between S and T bitimages, as add/subtract eliminates zeroes from both, into U
		// that is not necessarily the truth, ex: 5 + (-5) will give 0, but zeroing-out of values would happen extremely seldom
		// in real-life float operations
		if (copy) S.bitImage = bitImage.clone(this);
		if (addScalar)	S.bitImage.set();				// adding/subtracting a scalar will set every bit
		else			S.bitImage.or(T.bitImage);
		return S;
	}
	
	
	
	public CSRMatrix multiply(CSRMatrix T) {

		if (N != T.M) throw new RuntimeException("CSRMatrix.multiply(): Nonmatching matrix dimensions.");

		// results stored in new CSRMatrix U
		CSRMatrix U = new CSRMatrix("M", M, T.N);
		U.name = (Matrix.DEBUG_LEVEL > 1 ? "(" + name + "." + T.name + ")" : "M" + nameCount++);

		// create new bitImage, we will fill it up dynamically. Note down what type of bitimage it is
		boolean bitImage8x8 = false;
		int bitSets = 0;
		if (M < 9 && T.N < 9)  bitImage8x8 = true;
		else {
			int sets = T.N / 64, rest = T.N % 64;
			bitSets = sets + (rest == 0 ? 0 : 1);
		}
		U.bitImage = new BinBitImage(U);

		int[] SIA = IA, TIA = T.IA, SJA = JA, TJA = T.JA;
		double[] SA = A, TA = T.A;

		// multiplied result matrix U will never be bigger than either S nor T
		int maxAlen = M * T.N;
		// resize situational buffers if necessary, we'll copy from them later to correctly sized A1 & IA1
		if (Aop.length < maxAlen) {
			int newsize = Aop.length * 2 < maxAlen ? maxAlen : Aop.length * 2;	// at least double the situational buffers
			Aop = new double[newsize];
			JAop = new int[newsize];
		}
		int[] IA1 = new int[M + 1];
		double[] tmp = new double[T.N];
		int l = 0, iIA = 0;

		// sparse multiplication looks a bit as if the currently multiplied row was transposed and stretched out horisontally 
		// over the other matrix, any corresponding values found covered by the stretch of a value are multiplied
		// and added to the result matrix row whose index = index of column of the multiplicator value

		// for every row in first matrix S
		for (int iS = 0; iS < M; iS++) {
			// for every column/value pair in row jS
			int iSIA = SIA[iS + 1];
			
			// copy accumulated value from earlier IA1 slot to next slot for further accumulation
			IA1[++iIA] = IA1[iIA-1];

			// check if entire row is zeroes, if yes then skip temp-row accumulation & gathering
			int jS = SIA[iS];
			if (jS >= iSIA) continue;
			
			for (; jS < iSIA; jS++) {
				// do accumulation iteration over corresponding row in T
				// jS refers to the row of T to be added to temporary summation row
				int iTIA = TIA[SJA[jS]+1];
				for (int jT = TIA[SJA[jS]]; jT < iTIA; jT++)
					// the proper multiplication happens here
					tmp[TJA[jT]] += SA[jS] * TA[jT];
			}
			
			// finished with this row in S, gather up results from temp-row and clear it
			int iTN = iS * T.N, ibitSets = iS * bitSets;
			for (int k=0; k < T.N; k++) {
				if (!nearZero(tmp[k])) {		// only CSR-insert nonzero values, obviously
					Aop[l] = tmp[k];
					JAop[l++] = k;
					IA1[iIA]++;
					tmp[k] = 0;					// don't forget to clear this temp-row value
					// set appropriate bit in bitImage, according to either 8x8 or big-matrix strategy
					if (bitImage8x8)	U.bitImage.data[0] |= 0x1L<<(iTN + k);
					else				U.bitImage.data[ibitSets + (k>>6)] |= 0x1L<<(k%64);
				}
			}
		}
		
		U.A = new double[l];		
		U.JA = new int[l];
		U.bitImage.bitSets = bitSets;
		// move over data from situational buffers to newly allocated CSRmatrix
		for (int i = 0; i < l; i ++) { U.A[i] = Aop[i]; U.JA[i] = JAop[i]; }
		U.IA = IA1;
		U.nA = l;
		U.nIA = U.IA.length;
		
		if (Matrix.DEBUG_LEVEL > 1) System.out.println(U.toString());
		return U;
	}
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			DYNAMICALLY RESIZED CSR METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	// IA2 is a list of sparse row lists JA, preallocated to 16 positions each for fast fill-in at end of list
	// each list is reallocated on hitting 16-element bounds with size -> size + 16
	// each row starts with triplet: (tBl,tEc,sEc) total buffer length, total element count, sorted element count
	// after that, each element of each list is a triplet: (cIdx,uIdx,dIdx)
	// c is column inumber, uIdx & dIdx tell where the equal column number lies in the row above and below
	int nA2;
	private int[][] IA2, regIA2, nodeIA2;
	private double[][] A2, iA2;
	private int toSort;					// counts number of IA2 internal lists on the sorting stack 
	private int[] sortIA2;				// sortIA2 is an index stack for IA2 lists needing sorting
	private int[][] vGraph;				// stores nodes array-wise, for direct indexing
	private int[][] vGraphC;			// arranges nodes for column-wise indexing
	// valueOfColumn() method flushes column data into these temporary buffers
	int[] cCSR2;
	double[] cdCSR2;
	double[] icdCSR2;
	
	// identifiers for offsets into a JA2 buffer
	// tBl = offset to total buffer length param., tEc = offset to total element count param., sEc = offset to sorted element count param.
	private static final int tBl = 0, tEc = 1, sEc = 2, headerJA2 = 3;
	
		
	private void initCSR2(boolean complex) {
		nA2 = 0;
		IA2 = new int[M][];
		regIA2 = new int[M][];
		nodeIA2 = new int[M][];
		A2 = new double[M][];
		if (complex) iA2 = new double[M][];
		for (int i = 0; i < M; i++) {
			IA2[i] = new int[16];
			regIA2[i] = new int[headerJA2];
			nodeIA2[i] = new int[16];
			A2[i] = new double[16];
			if (complex) iA2[i] = new double[16];
			regIA2[i][tBl] = 16;
		}
		//sortIA2 = new int[16];
		cCSR2 = new int[M];
		cdCSR2 = new double[M];
		if (complex) icdCSR2 = new double[M];
	}
	
	
	public void putDataCSR2(double[] data, double[] idata) {
		
		if (IA2 == null) initCSR2(idata == null ? false : true);
		
		// first loop just inserts elements in proper order, the interlinking happens in next loop
		for (int i = 0; i < M; i++) {

			int[] JA2c = IA2[i];									// set up fast access references for JA2, regJA2, nodeJA2, A2 & iA2
			int[] regJA2c = regIA2[i], nodeJA2c = nodeIA2[i];
			double[] A2c = A2[i], iA2c = iA2 == null ? null : iA2[i];
			int iN = i * N, eLeft = regJA2c[tBl] - regJA2c[tEc];	// eLeft = number of free elements left in current JA2 buffer	
			int offsJA = 0;											// offset into current A2 & JA2 buffer
			
			for (int j = iN, jEnd = iN + N, col = 0; j < jEnd; j++, col++) {
				
				// if current JA2 is full (tBl - tEc = 0), add another 16 elements
				if (eLeft <= 0) {				
					int newSize = regJA2c[tBl] + 16;
					int[] oldJA = JA2c;
					JA2c = IA2[i] = new int[newSize];
					
					int eTotal = regJA2c[tEc];
					for (int c = 0; c < eTotal; c++)		// copy over data into new JA2 buffer
						JA2c[c] = oldJA[c];
					regJA2c[tBl] = newSize;					// update length of buffer
					
					double[] oldA = A2c;
					A2c = A2[i] = new double[newSize];		// copy over real values into new A2 buffer
					for (int c = 1; c < eTotal; c++)		A2c[c] = oldA[c];
					
					if (idata != null) {
						oldA = iA2c;						// copy over complex values into new A2 buffer
						iA2c = iA2[i] = new double[newSize];
						for (int c = 1; c < eTotal; c++)	iA2c[c] = oldA[c];						
					}
				}
				
				// we're inserting only nonzero values
				if (!nearZero(data[j]) || (idata != null && !nearZero(idata[j]))) {
					JA2c[offsJA] = col;						// insert column index
					nodeJA2c[offsJA] = -1;					// insert no-link flag
					A2c[offsJA] = data[j];					// insert non-zero value
					if (idata != null)						// insert non-zero complex value
						iA2c[offsJA] = idata[j];
					offsJA++;
					nA2++;									// global element count incremented
					regJA2c[tEc]++;							// total element count incremented
					eLeft--; 								// decrease number of free elements left in JA2
				}
			}
			if (regJA2c[tEc] == 0) {						// dereference empty rows
				IA2[i] = null;
				A2[i] = null;
				if (iA2 != null) iA2[i] = null;
			}
			else regJA2c[sEc] = regJA2c[tEc];
		}
		
		buildVGraph(data, idata);
	}
	
	

	
	private void buildVGraph(double[] data, double[] idata) {
		
		vGraph = new int[nA2][];									// every CSR2 element will be part of the graph, connected or not
		vGraphC = new int[N][];										// indexing of column-wise groups
		int nodeCnt = 0;											// counts up number of assembled nodes so far
		int[] linkJA2 = new int[M];									// stores indexes into JA2 arrays to use in crosslinking
		int[] linkRow = new int[M];
		
		int[] offJA2 = new int[M];									// holds current offsets into every JA2 sparse row

		// iterate over columns, this sets up a left-to-right marching through sparse column data
		for (int c = 0; c < N; c++) {
			
			int toLink = 0;											// counts number of vertical JA2 matches found for crosslinking
			// find all JA2 arrays containing current column c
			for (int r = 0; r < M; r++) {
				int[] JA2c = IA2[r], regJA2c = regIA2[r];
				if (offJA2[r] >= 0 && JA2c != null && JA2c[offJA2[r]] == c) {
					linkRow[toLink] = r;
					linkJA2[toLink++] = offJA2[r];
					offJA2[r]++;									// increment the marching offset into this JA2 row
					if (offJA2[r] >= regJA2c[tEc])					// if we have exhausted current JA2 row...
						offJA2[r] = -1;								// ...flag it as deactivated
				}
			}
			
			// create the column-wise group
			if (toLink > 0) {
				int groupSize = (toLink & 0xFFFFFFF0) + 32;
				vGraphC[c] = new int[groupSize];
			}
			vGraphC[c][0] = toLink;									// first entry of column group is the node count

			// time to allocate nodes for the found JA2 rows and do the vertical crosslinking
			int nodeCnt2 = nodeCnt + toLink;
			for (int node = nodeCnt, r = 0; node < nodeCnt2; node++, r++) {
				
				int nodeSize = 4;
				int[] vGNode =  vGraph[node] = new int[nodeSize];
				vGNode[0] = linkRow[r];
				vGNode[1] = linkJA2[r];
				nodeIA2[linkRow[r]][linkJA2[r]] = node;				// set link index from JA2 to current node
				vGraphC[c][r + 1] = node;
			}
			
			nodeCnt = nodeCnt2;
		}
		return;
	}
	
	
	public double valueOf2(int r, int c) {
		
		if (r < 0 || c < 0 || r >= M || c >= N) throw new RuntimeException("CSRMatrix.valueOf2(): Invalid matrix coordinates.");

		int[] JA2 = IA2[r], regJA2 = regIA2[r];
		if (JA2 == null) return 0;

		// do ordinary linear search for small arrays
		int eCnt = regJA2[sEc];							// search only through all sorted elements (for sake of consistency)
		if (eCnt <= 8)
			for (int j = 0; j < eCnt; j++)
				if (JA2[j] == c) return A2[r][j];
		
		// make binary search for sought element in JA2 sparse array, naturally only sorted elements apply
		int seekS = 0, seekE = eCnt, seekM = (seekE + seekS) >> 1;
//		while (JA2[seekM] != c) {
//			if (JA2[seekM] < c) 		{ seekS = seekM + 1; seekM = (seekE + seekS) >> 1; }
//			else if (JA2[seekM] > c)	{ seekE = seekM - 1; seekM = (seekE + seekS) >> 1; }
//			if (seekS >= seekE) return 0;
//		}
		while (seekS < seekE) {
			if (JA2[seekM] < c) 		{ seekS = seekM + 1; seekM = (seekE + seekS) >> 1; }
			else if (JA2[seekM] > c)	{ seekE = seekM - 1; seekM = (seekE + seekS) >> 1; }
			if (JA2[seekM] == c) return A2[r][seekM];
		}
		return 0;
	}
	
	
	public void valueTo2(int r, int c, double v) {
	
		if (r < 0 || c < 0 || r >= M || c >= N) throw new RuntimeException("CSRMatrix.valueOf2(): Invalid matrix coordinates.");

		int[] JA2 = IA2[r], regJA2 = regIA2[r];
		// if JA2 is a null-element array or if it needs to be rescaled
		if (JA2 == null || regJA2[tEc] >= regJA2[tBl]) {
			IA2[r] = JA2 = new int[16];					// if JA2 was empty, allocate the default 16 elements
			regJA2[tBl] = 16;							// set total buffer length
			regJA2[tEc] = regJA2[sEc] = 1;				// it's the only element of the JA2 array
			A2[r] = new double[16];
			if (isComplex()) iA2[r] = new double[16];
		}
		
		int eCnt = regJA2[sEc];							// search only through all sorted elements (for sake of consistency)
		if (eCnt <= 8)
			for (int j = 0; j < eCnt; j++)
				if (JA2[j] < c) {						// if we found the spot of insertion (next element is of lower column rank)
					
				}
		
	}
	
	
	public void valueOfColumn(int c) {
		if (c < 0 || c >= N) throw new RuntimeException("CSRMatrix.valueOfColumn(): Invalid column.");
		//for (int i = 0; i < )
	}
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	private static int MAX_PRINTEXTENT = 50;
	
	@Override
	public String toString() {
		
		int maxA = A.length > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : A.length;
		int maxIA = IA.length > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : IA.length;
		
		StringBuffer sb = new StringBuffer();
		sb.append(super.toString());						// get matrix data from superclass method
		sb.append("CSR data:\nIA:  [");
		for (int i = 1; i <= maxIA; i++) sb.append(IA[i-1] + (i==maxIA?"":", ") + (i%30==0?"\n      ":""));
		if (maxIA < IA.length)
				sb.append(" ... ]\nA:   [");
		else	sb.append("]\nA:   [");
		for (int i = 1; i <= maxA; i++) sb.append(String.format("%.2f" + (i==maxA-1?" ":", ") + (i%30==0?"\n      ":""), A[i-1]));
		if (maxA < A.length)
				sb.append(" ... ]\nJA:   [");
		else	sb.append("]\nJA:   [");
		if (iA != null) {
			for (int i = 1; i <= maxA; i++) sb.append(String.format("i%.2f" + (i==maxA-1?" ":", ") + (i%30==0?"\n      ":""), iA[i-1]));
			if (maxA < A.length)
					sb.append(" ... ]\nJA:   [");
			else	sb.append("]\nJA:   [");
		}
		for (int i = 1; i <= maxA; i++) sb.append(JA[i-1] + (i==maxA?"":", ") + (i%30==0?"\n      ":""));
		if (maxA < A.length)
				sb.append(" ... ]\n");
		else	sb.append("]\n");
		return sb.toString();
	}
	
	
	
	private String toStringCSR2() {
		int maxIA = IA2.length > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : IA2.length;
		StringBuffer sb = new StringBuffer();
		sb.append("CSR2 data:\n");
		for (int i = 0; i < maxIA; i++) {
			
			int[] JA2 = IA2[i], regJA2 = regIA2[i], nodeJA2 = nodeIA2[i];
			if (JA2 != null) {
				int unsorted = regJA2[tEc] - regJA2[sEc];
				sb.append("JA2_" + String.format("%d", i) + ": Sorted: [");
				for (int j = 0; j < regJA2[sEc]; j++)
					sb.append("(" + JA2[j] + ", node: " + nodeJA2[j] + (j == regJA2[sEc]-1 ? ")] " : "), "));
				
				if (unsorted > 0) sb.append(" unsorted: [");
				for (int j = 0; j < unsorted; j++)
					sb.append("(" + JA2[JA2[sEc] + 1] + ", node: " + JA2[JA2[sEc] + j] + (j == unsorted - 1 ? ")] " : "), "));
					
				sb.append("\nA2: [");
				int maxA = regJA2[tEc] > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : regJA2[tEc];
				for (int j = 0; j < maxA; j++)
					sb.append(String.format("%.2f" + (j==maxA-1?"]":", "), A2[i][j]));
				
				if (iA2 != null) {
					sb.append("\niA2: [");
					for (int j = 0; j < maxA; j++)
						sb.append(String.format("i%.2f" + (j==maxA-1?"]":", "), iA2[i][j]));				
				}
				sb.append("\n\n");
			} else
				sb.append("JA2_" + String.format("%d", i) + ": [empty]\n\n");
		}
		for (int i = 0; i < vGraph.length; i++)
			sb.append("Node " + i + ": [row: " + vGraph[i][0] + ", JA2 offset: " + vGraph[i][1] + "]\n");
		
		for (int i = 0; i < N; i++)
			if (vGraphC[i] != null) {
				sb.append("Column group " + i + ", nodes: [");
				for(int j = 0; j < vGraphC[i][0]; j++) sb.append(vGraphC[i][j+1] + (j == vGraphC[i][0] - 1 ? "]" : ", "));
				sb.append("\n");
			}

		
		return sb.toString();
	}

}
