package lkr74.matrixlib;

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
		if (idata != null) setComplex();
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
		
	int nNsp;
	NspNode[][] Hsp;			// horisontally/row aligned sparse arrays pointing to nodes (resembling JA in typical CSR)
	NspNode[][] Vsp;			// vertically/column aligned sparse arrays
	NspNode[] pivotNsp;			// fast-access to the diagonal nodes
	int[][] regHsp, regVsp;		// registers for Hsp & Vsp arrays
	
	// identifiers for offsets into a JA2 buffer
	// tBl = offset to total buffer length param., tEc = offset to total element count param., sEc = offset to sorted element count param.
	private static final int tBl = 0, tEc = 1, sEc = 2, NVSP_HEADER = 3, CSR2_ALLOCBLOCK = 16;
	
		
	
	// method converts flat array matrix data into CSR2 format
	public void putDataCSR2(double[] data, double[] idata) {
		
		// initialise all arrays
		nNsp = 0;
		Hsp = new NspNode[N][];
		Vsp = new NspNode[M][];
		pivotNsp = new NspNode[M];
		regHsp = new int[N][];
		regVsp = new int[M][];
		for (int i = 0; i < N; i++) regHsp[i] = new int[NVSP_HEADER];
		for (int i = 0; i < M; i++) regVsp[i] = new int[NVSP_HEADER];
		
		// first loop  inserts elements in proper order int Hsp
		for (int i = 0; i < M; i++) {

			NspNode[] bHsp = Hsp[i];								// set up fast access references for JA2, regJA2, nodeJA2, A2 & iA2
			int[] bregHsp = regHsp[i];
			int iN = i * N, eLeft = bregHsp[tBl] - bregHsp[tEc];	// eLeft = number of free elements left in current JA2 buffer	
			int offsHsp = 0;											// offset into current A2 & JA2 buffer
			
			// iterate over values in the matrix row
			for (int j = iN, jEnd = iN + N, col = 0; j < jEnd; j++, col++) {
				
				// if current JA2 is full (tBl - tEc = 0), add another 16 elements
				if (eLeft <= 0) {				
					int newSize = bregHsp[tBl] + CSR2_ALLOCBLOCK;
					NspNode[] oldbHsp = bHsp;
					bHsp = Hsp[i] = new NspNode[newSize];						// allocate row-wise reference array Hsp
					eLeft += CSR2_ALLOCBLOCK;
					
					int eTotal = bregHsp[tEc];
					for (int c = 0; c < eTotal; c++) bHsp[c] = oldbHsp[c];		// copy over data into new JA2 buffer		
					bregHsp[tBl] = newSize;										// update length of buffer
				}
				
				// we're inserting only nonzero values
				if (!nearZero(data[j]) || (isComplex() && !nearZero(idata[j]))) {
					bHsp[offsHsp] = new NspNode(i, col, data[j], idata != null ? idata[j] : 0);	// insert column index
					if (i == col) pivotNsp[i] = bHsp[offsHsp];					// if it's a pivot, put in fast-access array
					// for convenience, store the node's reference offset in Hsp, this allows to find a pivot quicker
					bHsp[offsHsp].offH = offsHsp;
					offsHsp++;
					nNsp++;											// global node count incremented
					bregHsp[tEc]++;									// total element count incremented
					eLeft--; 										// decrease number of free elements left in JA2
				}
			}
			if (bregHsp[tEc] == 0) Hsp[i] = null;					// dereference empty rows+
			else bregHsp[sEc] = bregHsp[tEc];						// TODO: increases algorithm complexity
			bregHsp[sEc] = bregHsp[tEc];
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
				NspNode[] bHsp = Hsp[r];
				int[] bregHsp = regHsp[r];
				if (offsHsp[r] >= 0 && bHsp != null && bHsp[offsHsp[r]].c == c) {
					linkRow[toLink] = r;
					idxHsp[toLink++] = offsHsp[r];
					offsHsp[r]++;									// increment the marching offset into this JA2 row
					if (offsHsp[r] >= bregHsp[tEc])					// if we have exhausted current JA2 row...
						offsHsp[r] = -1;							// ...flag it as deactivated
				}
			}
			
			// create the column-wise group
			if (toLink > 0) {
				int groupSize = trimToAllocBlock(toLink + 1);		// allocates in CSR2_ALLOCBLOCK sizes
				//int groupSize = toLink + 1;						// minimal +1 node allocation option
				Vsp[c] = new NspNode[groupSize];
				regVsp[c][tBl] = groupSize;							// store total buffer length
			}
			regVsp[c][tEc] = regVsp[c][sEc] = toLink;				// store node count in Vsp register

			// reference all found nodes for this column into current Vsp array
			int nodeCnt2 = nodeCnt + toLink;
			for (int node = nodeCnt, r = 0; node < nodeCnt2; node++, r++) {
				Vsp[c][r] = Hsp[linkRow[r]][idxHsp[r]];
				// for convenience, store the node's reference offset in Vsp, this allows to find a pivot quicker
				Vsp[c][r].offV = r;	
			}
			
			nodeCnt = nodeCnt2;
		}

	}
	
	

	private static int trimToAllocBlock(int v) {
		if ((v & 0xF) == 0) return (v & 0xFFFFFFF0); else return (v & 0xFFFFFFF0) + CSR2_ALLOCBLOCK;
	}
	

	// method searches for a NspNode in a bounded array, according to row if r >= 0, otherwise according to column c
	// if not found returns negated offset to the place where the node is supposed to be within a sorted order
	static int findHVspNode(NspNode[] bHVsp, int start, int end, int r, int c) {

		if (bHVsp == null) return -1;
		// do ordinary linear search for ranges less than 8
		if (end - start <= 8) {
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
		} else if (end - start <= 320) {
			// linear approximation finder
			if (r >= 0) {
				if (r < bHVsp[start].r) return -(start + 1);		// base case: supplied row is before the entire Vsp array
				if (bHVsp[end].r < r) return -(end + 2);			// base case: supplied row is after the entire Vsp array

				int rStart = bHVsp[start].r;
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

	// method chooses the smallest aspect (Hsp vs Vsp) to search for an element
	// if aspect = 0, returns offset into Hsp, if aspect = 1 returns offset into Vsp
	// if element wasn't found, returns negative value with 30th bit set if it's from Vsp, otherwise the bit is clear
	private int locateNspNode(int r, int c, int aspect) {
		
		if (r == c && pivotNsp[r] != null)
			return aspect == 0 ? pivotNsp[r].offH : pivotNsp[r].offV;	// if a pivot if sought, there is a fast-access buffer
		
		int nodesH = regHsp[r][sEc] - 1, nodesV = regVsp[c][sEc] - 1;
		if (nodesH < nodesV) {
			int offH = findHVspNode(Hsp[r], 0, nodesH, -1, c);
			if (offH < 0) return offH;
			return aspect > 0 ? Hsp[r][offH].offV : offH;
		} else {
			int offV = findHVspNode(Hsp[r], 0, nodesV, r, -1);
			if (offV < 0) return offV | 0x40000000;
			return aspect > 0 ? offV : Vsp[c][offV].offH;
		}
	}
	
	
	
	static NspNode[] addHVsp(
			NspNode[] bHVspA, NspNode[] bHVspB,
			int[] regHVspA, int[] regHVspB,
			double f, int aspectA, int aspectB, int aspectC) {
		
		int nNa = regHVspA[sEc], nNb = regHVspB[sEc];
		int nNc = trimToAllocBlock(nNa + nNb), a, b;
		NspNode[] bHVspC = new NspNode[nNc];
		
		for (int ia = 0, ib = 0, ic = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
			NspNode nA = bHVspA[ia], nB = bHVspB[ib];
			if (aspectA == 0) {
				if (aspectB == 0) 	{ a = nA.r; b = nB.r; }
				else				{ a = nA.r; b = nB.c; }
			} else {
				if (aspectB == 0) 	{ a = nA.c; b = nB.r; }
				else				{ a = nA.c; b = nB.c; }
			}
			if (a == b) {
				if (aspectC == 0)
						{ bHVspC[ic] = new NspNode(nA.r, nA.c, nA.v + nB.v * f, nA.iv + nB.iv * f); bHVspC[ic].offH = ic++; }
				else	{ bHVspC[ic] = new NspNode(nA.c, nA.r, nA.v + nB.v * f, nA.iv + nB.iv * f); bHVspC[ic].offV = ic++; }
				ia++; ib++;
			} else if (a < b) {
				if (aspectC == 0)
						{ bHVspC[ic] = new NspNode(nA.r, nA.c, nA.v, nA.iv); bHVspC[ic].offH = ic++; }
				else	{ bHVspC[ic] = new NspNode(nA.c, nA.r, nA.v, nA.iv); bHVspC[ic].offV = ic++; }
				ia++;
			} else {
				if (aspectC == 0)
						{ bHVspC[ic] = new NspNode(nB.r, nB.c, nA.v, nA.iv); bHVspC[ic].offH = ic++; }
				else	{ bHVspC[ic] = new NspNode(nB.c, nB.r, nA.v, nA.iv); bHVspC[ic].offV = ic++; }
				ib++;	
			}
		}
		return bHVspC;
	}


	// sparse CSR2 inner vector product, with variable Hsp/Vsp multiplying aspect
	static double multiplyHVsp(NspNode[] bHVspA, NspNode[] bHVspB, int[] regHVspA, int[] regHVspB, int aspectA, int aspectB) {
		
		int nNa = regHVspA[sEc], nNb = regHVspB[sEc];
		double v = 0;
	
		for (int ia = 0, ib = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
			NspNode nA = (bHVspA == null ? null : bHVspA[ia]);
			NspNode nB = (bHVspB == null ? null : bHVspB[ib]);
			if (nA == null || nB == null) break;
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
	
	
	// swaps two sparse arrays of index a & b, Hsp (rows) if aspect = 0, Vsp (columns) if aspect = 1
	public void swapHVspArrays(int a, int b, int aspect) {
		
		if (a == b) return;
		NspNode[] bHVspA = null, bHVspB = null, bTemp = null;
		NspNode nA, nB, nTemp = null;
		int nNa = 0, nNb = 0, temp;
		int[] regTemp;
		
		switch (aspect) { 
			case 0:
				bHVspA = Hsp[a]; bHVspB = Hsp[b]; nNa = regHsp[a][sEc]; nNb = regHsp[b][sEc];
				
				for (int ia = 0, ib = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
					
					nA = bHVspA == null ? null : bHVspA[ia];
					nB = bHVspB == null ? null : bHVspB[ib];
					// provide a fake boundary if one of the marching indexes run out of nodes
					int nAc = (nA == null ? N: nA.c), nBc = (nB == null ? N : nB.c);
					
					if (nAc == nBc) {										// found matching column nodes?
						if (nA.r == nAc) pivotNsp[nAc] = nB;				// if A was pivot, make B pivot
						else if (nB.r == nBc) pivotNsp[nBc] = nA;			// if B was pivot, make A pivot

						nTemp = Vsp[nAc][nA.offV];							// exchange their Vsp references
						Vsp[nAc][nA.offV] = Vsp[nAc][nB.offV];
						Vsp[nAc][nB.offV] = nTemp;
						temp = nA.offV; nA.offV = nB.offV; nB.offV = temp;	// exchange their offsets within the arrays	
						nA.r = b; nB.r = a;									// exchange their row indexes
						ia++; ib++;
						
					} else if (nAc < nBc)	{								// array A behind array B (meaning, zero at array B)?
						nA.r = b;											// change it's row to the other one
						if (nA.r == nAc) pivotNsp[nAc] = nA;				// see if A ended up as a pivot
						NspNode[] bVsp = Vsp[nAc];
						if (a < b) {										// do the switch within higher part of array b
							for (int j = nA.offV, k = j + 1, kEnd = regVsp[nAc][sEc]; j < kEnd; j++, k++) {
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
						NspNode[] bVsp = Vsp[nBc];
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
							for (int j = nB.offV, k = j + 1, kEnd = regVsp[nBc][sEc]; j < kEnd; j++, k++) {
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
				bTemp = bHVspA; Hsp[a] = bHVspB; Hsp[b] = bTemp;						// finally, exchange the rows in Hsp
				regTemp = regHsp[a]; regHsp[a] = regHsp[b]; regHsp[b] = regTemp;		// and swap the array registers
				break;	
				
			case 1:
				bHVspA = Vsp[a]; bHVspB = Vsp[b]; nNa = regVsp[a][sEc]; nNb = regVsp[b][sEc];
				for (int ia = 0, ib = 0; ia < nNa || ib < nNb;) {			// march through each node in the two arrays
					
					nA = bHVspA == null ? null : bHVspA[ia];
					nB = bHVspB == null ? null : bHVspB[ib];
					// provide a fake boundary if one of the marching indexes run out of nodes
					int nAr = (nA == null ? M: nA.r), nBr = (nB == null ? M: nB.r);
					
					if (nAr == nBr) {										// found matching row nodes?
						if (nA.c == nAr) pivotNsp[nAr] = nB;				// if A was pivot, make B pivot
						else if (nB.c == nBr) pivotNsp[nBr] = nA;			// if B was pivot, make A pivot

						nTemp = Hsp[nAr][nA.offH];							// exchange their Hsp references
						Hsp[nAr][nA.offH] = Hsp[nAr][nB.offH];
						Hsp[nAr][nB.offH] = nTemp;
						temp = nA.offH; nA.offH = nB.offH; nB.offH = temp;	// exchange their offsets within the arrays	
						nA.c = b; nB.c = a;									// exchange their row or column indexes
						ia++; ib++;
						
					} else if (nAr < nBr)	{								// array A behind array B (meaning, zero at array B)?
						nA.c = b;											// change it's column to the other one
						if (nA.c == nAr) pivotNsp[nAr] = nA;				// see if A ended up as a pivot
						NspNode[] bHsp = Hsp[nAr];
						if (a < b) {										// do the switch within right part (relative A) of array B
							for (int j = nA.offH, k = j + 1, kEnd = regHsp[nAr][sEc]; j < kEnd; j++, k++) {
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
						NspNode[] bHsp = Hsp[nB.r];
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
							for (int j = nB.offH, k = j + 1, kEnd = regHsp[nBr][sEc]; j < kEnd; j++, k++) {
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
				bTemp = bHVspA; Vsp[a] = bHVspB; Vsp[b] = bTemp;					// finally, exchange the column in Vsp
				regTemp = regVsp[a]; regVsp[a] = regVsp[b]; regVsp[b] = regTemp;	// and swap the array registers
				break;
		}
	}

	
	
	
	public double valueOf2(int r, int c) {
		
		if (r < 0 || c < 0 || r >= M || c >= N)
			throw new RuntimeException("CSRMatrix.valueOf2(): Invalid matrix coordinates.");
		if (Hsp[r] == null) return 0;
		// search in the aspect with fewest contained nodes
		if (regHsp[r][sEc] < regVsp[c][sEc]) {
			int offset = findHVspNode(Hsp[r], 0, regHsp[r][sEc] - 1, -1, c);
			if (offset >= 0) return Hsp[r][offset].v;						// return the value
		} else {
			int offset = findHVspNode(Vsp[c], 0, regVsp[c][sEc] - 1, r, -1);
			if (offset >= 0) return Vsp[c][offset].v;						// return the value
		}
		return 0;
	}
	
	

		
	// method inserts node into a Hsp/Vsp reference array, reallocating it if necessary
	// method returns the passed buffer, or the new buffer if reallocation happened
	// the aspect parameter decides way to insert: 0 -> row insert, 1 -> column insert
	// method returns the insertion offset or the found offset negated if node existed
	private int insertHVspNode(int r, int c, int aspect, NspNode node) {

		NspNode[] bHVsp = null;
		int[] regHVsp = null;
		int offset = 0;
		switch (aspect) {
			case 0:	{	bHVsp = Hsp[r]; regHVsp = regHsp[r];
						offset = findHVspNode(bHVsp, 0, regHVsp[sEc] - 1, -1, c); break; }
			case 1:	{	bHVsp = Vsp[c]; regHVsp = regVsp[c];
						offset = findHVspNode(bHVsp, 0, regHVsp[sEc] - 1, r, -1); break; }
		}
		// if we didn't find element, insert it
		if (offset < 0) {
			offset = -offset-1;
			regHVsp[sEc]++;
			regHVsp[tEc]++;
			if (r == c) pivotNsp[r] = bHVsp[offset];			// if pivot, insert in fast-access array

			// does Hsp/Vsp reference array need alocation/reallocation ?
			if (regHVsp[tEc] > regHVsp[tBl]) {
				NspNode[] bHVsp2 = new NspNode[trimToAllocBlock(regHVsp[tEc])];
				if (aspect == 0) Hsp[r] = bHVsp2; else Vsp[c] = bHVsp2;
				
				int i = regHVsp[tEc] - 1;
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
				int i = regHVsp[tEc] - 1;
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
	private NspNode[] insertLocalHVspNode(NspNode[] bHVsp, int offset, int aspect, int[] regHVsp, NspNode node) {

		regHVsp[sEc]++;
		int nNodes = regHVsp[tEc]++;
		if (bHVsp[offset].r == bHVsp[offset].c)
			pivotNsp[bHVsp[offset].r] = bHVsp[offset];		// if pivot, insert in fast-access array

		// does Hsp/Vsp reference array need allocation/reallocation ?
		if (nNodes > regHVsp[tBl]) {
			NspNode[] bHVsp2 = new NspNode[trimToAllocBlock(nNodes)];
			
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
			return bHVsp2;
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
		return bHVsp;												
	}

	
	
	private int removeHVspNode(int r, int c, int aspect) {

		NspNode[] bHVsp = null;
		int[] regHVsp = null;
		int offset = 0;
		
		switch (aspect) {
			case 0:	{	bHVsp = Hsp[r]; regHVsp = regHsp[r];
						offset = findHVspNode(bHVsp, 0, regHVsp[sEc]-1, -1, c); break; }
			case 1:	{	bHVsp = Vsp[c]; regHVsp = regVsp[c];
						offset = findHVspNode(bHVsp, 0, regHVsp[sEc]-1, r, -1); break; }
		}
		
		// if we found element, remove it
		if (offset >= 0) {
			if (r == c) pivotNsp[r] = null;								// if it's a pivot, remove from fast-access array
			
			if (regHVsp[sEc] > offset) regHVsp[sEc]--;					// if we're removing within bounds of sorted elements, decrement sEc
			int nNodes = --regHVsp[tEc];
			if (nNodes == 0) {
				regHVsp[tBl] = 0;
				if (aspect == 0) Hsp[r] = null;
				else Vsp[c] = null;										// last node in array removed, destroy this array
				return offset;		
			}
			
			// does Hsp/Vsp reference array need deallocation (over 63 nodes been removed)?
			if (regHVsp[tBl] - nNodes >= 64) {
				NspNode[] bHVsp2 = new NspNode[regHVsp[tBl] - 64];
				
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
	
	
	private NspNode[] removeLocalHVspNode(NspNode[] bHVsp, int offset, int aspect, int[] regHVsp) {
		
		if (bHVsp[offset].r == bHVsp[offset].c)
			pivotNsp[bHVsp[offset].r] = null;						// if it's a pivot, remove from fast-access array
		if (regHVsp[sEc] > offset) regHVsp[sEc]--;					// if we're removing within bounds of sorted elements, decrement sEc
		int nNodes = --regHVsp[tEc];
		if (nNodes == 0) { regHVsp[tBl] = 0; return null; }			// last node in array removed, signal array destruction
		
		// does Hsp/Vsp reference array need deallocation (over 63 nodes been removed)?
		if (regHVsp[tBl] - nNodes >= 64) {
			NspNode[] bHVsp2 = new NspNode[regHVsp[tBl] - 64];
			
			int i = 0;
			if ((aspect & 1) == 0)
					for (; i < offset; i++) { bHVsp2[i] = bHVsp[i];	bHVsp2[i].offH--; }		// move the initial elements in Hsp to new Hsp block
			else	for (; i < offset; i++) { bHVsp2[i] = bHVsp[i];	bHVsp2[i].offV--; }		// move the initial elements in Vsp to new Vsp block
			for (int j = i + 1; j < nNodes; i++, j++) bHVsp2[i] = bHVsp[j];					// shift elements left in new Hsp/Vsp block
			return bHVsp2;
		} else {
			if ((aspect & 1) == 0)
					for (int i = offset, j = i + 1; i < nNodes; i++, j++) { bHVsp[i] = bHVsp[j]; bHVsp[i].offH--; }
			else	for (int i = offset, j = i + 1; i < nNodes; i++, j++) { bHVsp[i] = bHVsp[j]; bHVsp[i].offV--; }
		}
		return bHVsp;													// if element didn't exist, return negated offset
	}

	
	
	
	public void valueTo2(int r, int c, double v) {
		if (r < 0 || c < 0 || r >= M || c >= N)
			throw new RuntimeException("CSRMatrix.valueOf2(): Invalid matrix coordinates.");
		
		if (nearZero(v)) {												// insertion of zero means removal of a value
			int nNodes = regHsp[r][sEc];
			NspNode[] bHsp = Hsp[r];
			int offH = findHVspNode(bHsp, 0, nNodes - 1, -1, c);		// locate offset of node in Hsp
			int offV;
			if (offH >= 0) offV = bHsp[offH].offV;						// get node's offset into Vsp if node existed in Hsp
			else return;												// a node that doesn't exist in Hsp will not exist in corresponding Vsp
			Hsp[r] = removeLocalHVspNode(Hsp[r], offH, 0, regHsp[r]);	// dereference node from horisontal sparse array
			Vsp[c] = removeLocalHVspNode(Vsp[c], offV, 1, regVsp[c]);	// dereference node from vertical sparse array
			return;
		}
		NspNode node = new NspNode(r, c, v, 0);
		int offset = 		insertHVspNode(r, c, 0, node);				// insert value as new node into Hsp
		if (offset >= 0) 	insertHVspNode(r, c, 1, node);				// if node didn't exist already, insert value as new node into Vsp
	}



	
	// returns total approximate size of the CSR2 structure, a reference counted as 2 ints
	public int sizeOfCSR2() {
		final int refSize = 2;
		int totalSize = Hsp.length * refSize + regHsp.length * (refSize + 3) + Vsp.length * refSize;
		for (NspNode[] l : Hsp)
			if (l != null) {
				totalSize += l.length;
				for (NspNode n: l) totalSize += 16;
			}
		for (NspNode[] l : Vsp) if (l != null) totalSize += l.length;
		return totalSize;
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
	
	
	
	public String toStringCSR2() {
		int maxIA = Hsp.length > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : Hsp.length;
		StringBuffer sb = new StringBuffer();
		sb.append("CSR2 data:\nHsp row groups:\n");
		for (int i = 0; i < maxIA; i++) {
			
			NspNode[] bHsp = Hsp[i];
			int[] bregHsp = regHsp[i];
			if (bHsp != null) {
				sb.append("Hsp(" + i + "): Sorted: [");
				int maxHsp = bregHsp[sEc] > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : bregHsp[sEc];
				for (int j = 0; j < maxHsp; j++)
					sb.append("(" + 	bHsp[j].r + ", " + 
										bHsp[j].c + ", " + 
										bHsp[j].v + ", " + 
										bHsp[j].offH + ", " + 
										bHsp[j].offV + (j == maxHsp-1 ? ")]\n" : "), "));
				
				int unsorted = (MAX_PRINTEXTENT - bregHsp[sEc]);
				if (unsorted > 0) {
					int unsorted2 = bregHsp[tEc] - bregHsp[sEc];
					if (unsorted2 < unsorted) unsorted = unsorted2;
					if (unsorted > 0) sb.append(" Unsorted: [");
					for (int j = 0; j < unsorted; j++)
						sb.append("(" + bHsp[bregHsp[sEc] + j].r + ", " +
										bHsp[bregHsp[sEc] + j].c + ", " +
										bHsp[bregHsp[sEc] + j].v + ", " + 
										bHsp[bregHsp[sEc] + j].offH + ", " +
										bHsp[bregHsp[sEc] + j].offV + (j == unsorted - 1 ? ")]\n" : "), "));
				}
			} else
				sb.append("Hsp(" + i + "): [empty]\n\n");
		}
		
		sb.append("\n\nVsp column groups::\n");
		maxIA = Vsp.length > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : Vsp.length;
		for (int i = 0; i < maxIA; i++) {
			
			NspNode[] bVsp = Vsp[i];
			int[] bregVsp = regVsp[i];
			if (bVsp != null) {
				sb.append("Vsp(" + i + "): Sorted: [");
				int maxVsp = bregVsp[sEc] > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : bregVsp[sEc];
				for (int j = 0; j < maxVsp; j++)
					sb.append("(" +		bVsp[j].r + ", " +
										bVsp[j].c + ", " +
										bVsp[j].v + ", " +
										bVsp[j].offH + ", " +
										bVsp[j].offV + (j == maxVsp-1 ? ")]\n" : "), "));
				
				int unsorted = (MAX_PRINTEXTENT - bregVsp[sEc]);
				if (unsorted > 0) {
					int unsorted2 = bregVsp[tEc] - bregVsp[sEc];
					if (unsorted2 < unsorted) unsorted = unsorted2;
					if (unsorted > 0) sb.append(" Unsorted: [");
					for (int j = 0; j < unsorted; j++)
						sb.append("(" + bVsp[bregVsp[sEc] + j].r + ", " +
										bVsp[bregVsp[sEc] + j].c + ", " +
										bVsp[bregVsp[sEc] + j].v + (j == unsorted - 1 ? ")] " : "), "));
				}
			} else
				sb.append("Vsp(" + i + "): [empty]\n\n");
		}

		sb.append("Matrix size: " + (M*N*4) + "\n");
		sb.append("Total CSR2 structure size: " + sizeOfCSR2() + "\n");
		return sb.toString();
	}
}


class NspNode {
	int r, c, offH, offV;
	double v, iv;
	NspNode(int r, int c, double v, double iv) { this.r=r; this.c=c; this.v=v; this.iv=iv; }
}

