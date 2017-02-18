package lkr74.matrixlib;

import java.security.InvalidParameterException;

public class FrontalMatrix {

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			FRONTAL MATRIX DATASTRUCTURE FOR MULTIFRONTAL LINEAR SOLVER
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	final static int FM_ALLOCBLOCK = 16, FAILED_VECSIZE = 3;
	final static boolean EXTRA_FM_SPACE = false;
	
	int pivot = -1;												// pivot (uninitialised = -1)
	int r, c, rTot, cTot;										// data height & width and total height & width of frontal matrix
	// nFailed is number of failed pivot row/columns in top left position, nPivots is number of pivot row/columns
	// m/nContrib is MxN rows/columns of contribution part, m/nExtra is MxN extra rows/columns at bottom-left:
	//
	// ffffffffff	<- f:ailed pivots row/cols
	// fPPPPPPPPP	<- P:ivoting row/cols
	// fPPPPPPPPP
	// fPPCCCCCee   <- C:ontribution matrix
	// fPPCCCCCee
	// fPPeeeeeee	<- e:xtra col/rows
	// no. of failed pivot row/cols, no. of pivot row/cols, size of contrib.matrix
	int nFailed = 0, nPivots = 0, nContribR = 0, nContribC = 0;
	// no. of extra rows, no. of extra columns, no. of extra pivots
	int nExtraR = 0, nExtraC = 0, nExtraP = 0, nExtraP2 = 0;
	int[] idxCR;												// global row & column indexes of constituent frontal rows & columns
	int[][] extraR = null, extraC = null;						// holds sorted list of extra appended rows & columns & pivot row/columns from children
	int[] failedP = null;										// holds sorted list of failed pivot indexes
	double[] data;												// data array, linear layout
	byte[] mixCode = null;										// the optimisation of mixcodes allows fast sparse joining of child & parent values
	int iTB;													// offset to temp.row/column assemblies
	int childFlag = 0;											// child type flag for add() method
	
	// creates a skeleton frontal matrix
	public FrontalMatrix() { super(); }

	@Override
	public FrontalMatrix clone() {
		try {
			FrontalMatrix F = (FrontalMatrix)super.clone();		// will first call superclass clone() method
			F.idxCR = idxCR == null ? null : idxCR.clone();
			// the following three are really only transitive index holders for parent nodes during symbolic factorising
//			F.extraR = extraR == null ? null : extraR.clone();
//			F.extraC = extraC == null ? null : extraC.clone();
//			F.failedP = failedP == null ? null : failedP.clone();
			F.data = data == null ? null : data.clone();
			return F;
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	
	// creates a skeleton frontal matrix with parameter setting and allocated datafield
	public FrontalMatrix(int p, int r, int c, int nContribR, int nContribC, int nExtraR, int nExtraC) {
		pivot = p;
		this.r = r; this.c = c;
		if (EXTRA_FM_SPACE) {					// optionally allocate for at least one failed pivot, plus 1.25x extra space
			rTot = r + r/4 + 1;	 cTot = c + c/4 + 1;
		} else { rTot = r; cTot = c; }
		this.nContribR = nContribR; this.nContribC = nContribC;
		this.nExtraR = nExtraR; this.nExtraC = nExtraC;
		iTB = rTot * cTot;
		// additionally allocate space for global indexation of rows & columns and a temporary row/column assembly buffer
		data = new double[iTB + (rTot > cTot ? rTot : cTot)*2];	
		idxCR = new int[rTot + cTot + 1];
		mixCode = new byte[rTot > cTot ? rTot + 1: cTot + 1];
		nFailed = 0;
	}

	
	
	// creates a FrontalMatrix for a supplied array of coefficients, so the frontal matrix can be used directly
	// as a full pivoting factoriser with with decomposeLU() method
	public FrontalMatrix(int r, int c, double[] dataIn) {

		if (r < 1 || c < 1 || r * c < (dataIn == null ? 0 : dataIn.length))
			throw new InvalidParameterException("FrontalMatrix.FrontalMatrix(): Illegal matrix dimensions.");
		this.pivot = 0;
		nPivots = r < c ? r : c;					// for non-square matrix, choose the smallest dimension for no. of pivots
		this.r = rTot = r; this.c = cTot = c;
		iTB = rTot * cTot;
		// additionally allocate space for global indexation of rows & columns and a temporary row/column assembly buffer
		data = new double[iTB + (rTot > cTot ? rTot : cTot)*2];	
		idxCR = new int[rTot + cTot + 1];
		mixCode = new byte[rTot > cTot ? rTot + 1: cTot + 1];
		nFailed = 0;

		// initialise indexing of the rows & columns (will be swapped dutring eventual pivoting)
		for (int j = 0; j < c; j++) idxCR[j] = j;
		for (int i = c; i < c + r; i++) idxCR[i] = i - c;		
		if (dataIn != null) for (int i = 0, iEnd = r * c; i < iEnd; i++) data[i] = dataIn[i];
	}

	
	
	// creates a frontal matrix with pivot row & column interspersed with expected contribution indexes
	// note: gContribC/gContribR,gExtraC are arrays with a special structure, generated in FrontalDAG class
	public FrontalMatrix(NSPMatrix A, int pivot, int nExtraP, int[][] gContribC, int[][] gContribR, int[][] gExtraC, int[][] gExtraR) {
		
		NspNode[] bHsp = A.Hsp[pivot].array, bVsp = A.Vsp[pivot].array;
		// if there are no values in the row or column, return empty frontal matrix, let caller handle problemn
		if (bHsp == null || bVsp == null) return;

		this.pivot = pivot;
		this.nExtraP = nExtraP;
		if (gExtraR[1] != null) { this.nExtraR = gExtraR[2][0]; this.extraR = gExtraR; }
		if (gExtraC[1] != null) { this.nExtraC = gExtraC[2][0]; this.extraC = gExtraC; }
		this.nPivots = 1;

		int offH = A.pivotNsp[pivot].offH + 1, offV = A.pivotNsp[pivot].offV + 1;
		nContribR = gContribR[2][0] - 1;									// the total final dimensions of contribution block (minus pivot)
		nContribC = gContribC[2][0] - 1;
		r = nExtraP + 1 + nContribR + nExtraR;
		c = nExtraP + 1 + nContribC + nExtraC;
		if (EXTRA_FM_SPACE) { rTot = r + r/4; cTot = c + c/4; }				// allocate with at least 1.25x extra space (rTot=r*1.25, cTot=c*1.25)
		else { rTot = r; cTot = c; }
		data = new double[rTot * cTot + (rTot > cTot ? rTot : cTot)*2];
		idxCR = new int[rTot + cTot + 1];
		mixCode = new byte[rTot > cTot ? rTot + 1: cTot + 1];				// mixCode will store codes for consecutive unrolled-loop value insertion
		
		// initialise frontal row & column from pivot's sparse row & column
		// the offset array, gotten from analysis of children contributions, allows placing the values from pivot's row & column
		// directly at their correct positions in frontal matrix, leaving gaps for children contributions
		int offP = nExtraP * cTot + nExtraP;								// offset beyond extra pivot rows & columns
		data[offP] = A.pivotNsp[pivot].v;									// insert pivot value
		int[] contribC = gContribC[0], offset = gContribC[2];
		int j = offH, jO = FrontalDAG.HELPARRAY_SIZE + 1, jIdx = nExtraP;
		int nContribCP = gContribC[2][2] + offH - 1, nContribRP = gContribR[2][2] + offV - 1;
		//idxCR[jIdx++] = pivot;
		//while (contribC[j2] < p) j2++;									// skip past contrib.column indexes lower or equal than current pivot
			
		while (j < nContribCP)	data[offP + offset[jO++]] = bHsp[j++].v;	// insert sparse row values from NSPMatrix
		for (j = 0; j <= nContribC; j++) idxCR[jIdx++] = contribC[j];		// set column indexes for pivots + total contribution
			
		int[] extras = gExtraC[0];
		for (int jExC = 0; jExC < nExtraC; jExC++)							// set the incoming extras global column indexes
			idxCR[jIdx++] = extras[jExC];

		int i = offV, iO = FrontalDAG.HELPARRAY_SIZE + 1,  iIdx = cTot + nExtraP;
		idxCR[iIdx++] = pivot;												//
		int[] contribR = gContribR[0];
		offset = gContribR[2];
		//while (contribR[i2] < p) i2++;									// skip past contrib.column indexes lower or equal than current pivot
			
		while (i < nContribRP)												// insert sparse column values from NSPMatrix
			data[offP + offset[iO++] * cTot] = bVsp[i++].v;
		for (i = 1; i <= nContribR; i++) idxCR[iIdx++] = contribR[i];		// set column indexes for pivots + total contribution
		
		extras = gExtraR[0];
		for (int iExC = 0; iExC < nExtraR; iExC++)							// set the incoming extras global column indexes
			idxCR[iIdx++] = extras[iExC];
		
	}
	
	
	// creates a  frontal matrix with pivot row & column initialised, space in front for failed pivots,
	// space in back for extra rows & columns
	public FrontalMatrix(NSPMatrix M, int pivot, int nExtraP, int nExtraR, int nExtraC) {
				
		NspNode[] bHsp = M.Hsp[pivot].array, bVsp = M.Vsp[pivot].array;
		// if there are no values in the row or column, return empty frontal matrix, let caller handle problemn
		if (bHsp == null || bVsp == null) return;

		this.pivot = pivot;
		nPivots = 1;
		this.nExtraP = nExtraP;
		this.nExtraR = nExtraR;
		this.nExtraC = nExtraC;

		int offH = M.pivotNsp[pivot].offH, offV = M.pivotNsp[pivot].offV;
		nContribR = M.Vsp[pivot].nodes - offV - 1;
		nContribC = M.Hsp[pivot].nodes - offH - 1;
		r = nExtraP + 1 + nContribR + nExtraR;
		c = nExtraP + 1 + nContribC + nExtraC;
		if (EXTRA_FM_SPACE) { rTot = r + r/4; cTot = c + c/4; }				// allocate with space for at least 1.25x extra space
		else { rTot = r; cTot = c; }										// allocate the exact space needed
		data = new double[rTot * cTot + (rTot > cTot ? rTot : cTot)*2];
		idxCR = new int[rTot + cTot + 1];
		//mixCode = new byte[rTot > cTot ? rTot + 1 : cTot + 1];
		
		// initialise frontal row & column from pivot's sparse row & column
		for (int j = offH, jFM = nExtraP; jFM < c; j++, jFM++) {
			data[jFM] = bHsp[j].v;
			idxCR[jFM] = bHsp[j].c;
		}
		int iEnd = M.Vsp[pivot].nodes;
		for (int i = offV, iFM = nExtraP, iFM2 = cTot + iFM; i < iEnd; i++, iFM += cTot, iFM2++) {
			data[iFM] = bVsp[i].v;
			idxCR[iFM2] = bVsp[i].r;
		}
		
		this.setLeaf();														// this is flagged as leaf frontal
	}	
	
		
	// method combinatively checks for pattern match between pivot1 and pivot2 and assembles a supernode frontal matrix
	// method demands that pivot1 < pivot2, if the caller is continuing assembly of a previous matrix, then it's supplied in "fm"
	// pivot1 must be the supernode pivot of the supplied fm
	// TODO: method currently tries matching both horisontally & vertically (which is less prone to make large supernodes), make switch to activate only horisontal matching
	// zAllowed is the base number of zeroes allowed for fill-in to make a looser criterion for a supernode match
	// zAllowedF is a factor of "zero fill-ins allowed"/element (for large supernodes)
	// after initial generation, all consecutive member candidates must fit element-count-wise into the supernode pivot's row & column
	public static FrontalMatrix assembleSupernode(NSPMatrix M, int pivot1, int pivot2, int zAllowed, float zAllowedF, FrontalMatrix fm) {
		
		boolean newFM = false;
		int offPH = M.pivotNsp[pivot1].offH, offPV = M.pivotNsp[pivot1].offV, iMatchC = offPH, iMatchR = offPV;
		NspNode[] bHsp1 = M.Hsp[pivot1].array, bVsp1 = M.Vsp[pivot1].array;
		int p1nodesH = M.Hsp[pivot1].nodes, p1nodesV = M.Vsp[pivot1].nodes; 
		
		// advance index from pivot1 to start of pivot2 along row & column, as comparisons will start there
		while (iMatchC < p1nodesH - 1 && bHsp1[iMatchC + 1].c <= pivot2) iMatchC++;
		while (iMatchR < p1nodesV - 1 && bVsp1[iMatchR + 1].r <= pivot2) iMatchR++;
		
		// do quick check that element counts match (within -/+ zAllowed+zAllowedF*nElems)
		int nElemsR1 = p1nodesV - bVsp1[iMatchR].offV, nElemsC1 = p1nodesH - bHsp1[iMatchC].offH;
		int nElemsR2 = M.Vsp[pivot2].nodes - M.pivotNsp[pivot2].offV, nElemsC2 = M.Hsp[pivot2].nodes - M.pivotNsp[pivot2].offH;
		int zAllowedR = zAllowed + (int)(zAllowedF * (float)nElemsR1);
		int zAllowedC = zAllowed + (int)(zAllowedF * (float)nElemsC1);
		int nDeltaR= nElemsR1 - nElemsR2, nDeltaC = nElemsC1 - nElemsC2;
	
		// fail if element difference over limit of allowed fill-in elements in row or column aspect
		if (nDeltaR > zAllowedR || nDeltaC > zAllowedC) return null;

		if (fm == null) {								// on creation of a frontal matrix, allow zero fill-in flexibility in first pivot
			int r = (nElemsR1 > nElemsR2 ? nElemsR1 : nElemsR2) + 1; 			// +1 includes the pivot
			int c = (nElemsC1 > nElemsC2 ? nElemsC1 : nElemsC2) + 1; 
			fm = new FrontalMatrix(pivot1, r, c, r - 2, c - 2, 0, 0);			// make skeleton frontal matrix
			newFM = true;
			if (nDeltaR < 0) nDeltaR = -nDeltaR;
			if (nDeltaC < 0) nDeltaC = -nDeltaC;
		} else {										// on continuing assembly, all sequential pivots must fit into first pivot's submatrix
			if (nElemsR2 > fm.r - fm.nPivots || nElemsC2 > fm.c - fm.nPivots) return null;
		}
		
			
		double[] data = fm.data;
		int fmPivOff = fm.nExtraP + fm.nPivots - (newFM?0:1);	// since nPivots goes from 0->2 directly on new assembly, subtract 1 in continuation phase
		int n = fm.cTot, iMH1 = fmPivOff * n + fmPivOff, iMH2 = iMH1 + n + 1, iMV1 = iMH1 + n, iMV2 = iMV1 + 1;
		int iTBR = fm.iTB;										// iTB indexes the temporarily built rows
		
		// check for matching pattern by row aspect
		NspNode[] bHsp2 = M.Hsp[pivot2].array;
		// i1 iterates over the outer/lower pivot row, i2 iterates over the inner/higher (which has one element less at the beginning)
		int i1 = iMatchC, i2 = M.pivotNsp[pivot2].offH;
		int iEnd1 = p1nodesH, iEnd2 = M.Hsp[pivot2].nodes;
		int[] idxRC = fm.idxCR;
		
		// the case of a new frontal matrix and adding two rows simultaneously while comparing
		if (newFM) {
			int ixC = fmPivOff;
			if (offPH != iMatchC) {						// store down element of pivot1 (ignore entire pivot1 column/row if continuing assembly)		
						data[iMH1++] = bHsp1[i1 - 1].v;	idxRC[ixC++] = bHsp1[i1 - 1].c;
			} else {	data[iMH1++] = bHsp1[i1].v;		idxRC[ixC++] = bHsp1[i1].c; }
			
			while (i1 < iEnd1 && i2 < iEnd2) {
				int c1 = bHsp1[i1].c, c2 = bHsp2[i2].c;
				if (c1 < c2) {									// was pivot1's column element before pivot2's column element (mismatch)?
					if (--zAllowedR < 0) return null;			// if we ran out of allowed row fill-in elements, fail the supernode matching
					data[iMH1++] = bHsp1[i1++].v;				// otherwise, advance in pivot1's sparse row, copying over to frontal matrix
					iMH2++;
					idxRC[ixC++] = c1;
				} else if (c1 > c2) {							// was pivot1's column element after pivot2's column element (mismatch)?
					if (--zAllowedR < 0) return null;			// if we ran out of allowed row fill-in elements, fail the supernode matching
					data[iMH2++] = bHsp2[i2++].v;				// otherwise, advance in pivot2's sparse row, copying over to frontal matrix
					iMH1++;
					idxRC[ixC++] = c2;
				} else {
					data[iMH1++] = bHsp1[i1++].v;				// column elements matched, advance in both rows
					data[iMH2++] = bHsp2[i2++].v;
					idxRC[ixC++] = c1;
				}
			}
			// if one of the rows still has elements, test if remaining elements do not exceed remaining allowed fill-ins
			nDeltaR = iEnd1 - i1 > 0 ? iEnd1 - i1 : iEnd2 - i2;
			if (nDeltaR > zAllowedR) return null;				// if remaining elements exceed remaining allowed fill-ins, fail
			while (i1 < iEnd1) { data[iMH1++] = bHsp1[i1].v; idxRC[ixC++] = bHsp1[i1++].c; }
			while (i2 < iEnd2) { data[iMH2++] = bHsp2[i2].v; idxRC[ixC++] = bHsp2[i2++].c; }
			
		// the case of adding another row & column to an assembled frontal matrix, pivot1's row is now only used for comparison
		} else {
			while (i1 < iEnd1 && i2 < iEnd2) {
				int c1 = bHsp1[i1].c, c2 = bHsp2[i2].c;
				if (c1 < c2) {									// was pivot1's column element before pivot2's column element (mismatch)?
					if (--zAllowedR < 0) return null;			// if we ran out of allowed row fill-in elements, fail the supernode matching
					data[iTBR++] = 0;
					i1++;
				} else if (c1 > c2) {							// was pivot1's column element after pivot2's column element (mismatch)?
					if (--zAllowedR < 0) return null;			// if we ran out of allowed row fill-in elements, fail the supernode matching
					data[iTBR++] = bHsp2[i2++].v;				// otherwise, advance in pivot2's sparse row, copying over to frontal matrix
				} else {
					data[iTBR++] = bHsp2[i2++].v;
					i1++;
				}
			}
			if (iEnd2 - i2 > zAllowedR) return null;			// if remaining elements of pivot2 exceed remaining allowed fill-ins, fail
			while (i2 < iEnd2) data[iTBR++] = bHsp1[i2++].v;	// otherwise copy them over
		}
		
		// check for matching pattern by column aspect
		NspNode[] bVsp2 = M.Vsp[pivot2].array;
		i1 = iMatchR; i2 = M.pivotNsp[pivot2].offV;
		iEnd1 = p1nodesV; iEnd2 = M.Vsp[pivot2].nodes;
		int iTBC = iTBR;
		
		// the case of a new frontal matrix and adding two columns simultaneously while comparing
		if (newFM) {
			int ixR = fm.cTot + fmPivOff;
			if (offPV != iMatchR)	idxRC[ixR++] = bVsp1[i1 - 1].r;
			else					idxRC[ixR++] = bVsp1[i1].r;
			
			while (i1 < iEnd1 && i2 < iEnd2) {
				int r1 = bVsp1[i1].r, r2 = bVsp2[i2].r;
				if (r1 < r2) {									// was pivot1's row element before pivot2's row element (mismatch)?
					if (--zAllowedC < 0) return null;			// if we ran out of allowed column fill-in elements, fail the supernode matching
					data[iMV1] = bVsp1[i1++].v;					// otherwise, advance in pivot1's sparse column, copying over to frontal matrix
					iMV1 += n; iMV2 += n;
					idxRC[ixR++] = r1;
				} else if (r1 > r2) {							// was pivot1's row element after pivot2's row element (mismatch)?
					if (--zAllowedC < 0) return null;			// if we ran out of allowed column fill-in elements, fail the supernode matching
					data[iMV2] = bVsp2[i2++].v;					// otherwise, advance in pivot2's sparse column, copying over to frontal matrix
					iMV1 += n; iMV2 += n;
					idxRC[ixR++] = r2;
				} else {
					data[iMV1] = bVsp1[i1++].v;					// row elements matched, advance in both columns
					data[iMV2] = bVsp2[i2++].v;
					iMV1 += n; iMV2 += n;
					idxRC[ixR++] = r1;
				}
			}
			// if one of the columns still has elements, test if remaining elements do not exceed remaining allowed fill-ins
			nDeltaC = iEnd1 - i1 > 0 ? iEnd1 - i1 : iEnd2 - i2;
			if (nDeltaC > zAllowedC) return null;				// remaining elements exceed remaining allowed fill-ins, fail
			while (i1 < iEnd1) { data[iMV1] = bVsp1[i1].v; idxRC[ixR++] = bVsp1[i1++].r; iMV1 += n;  }
			while (i2 < iEnd2) { data[iMV2] = bVsp2[i2].v; idxRC[ixR++] = bVsp2[i2++].r; iMV2 += n; }	
			fm.nPivots = 2;
			
		// the case of adding another row & column to an assembled frontal matrix, pivot1's column is now only used for comparison
		} else {
			while (i1 < iEnd1 && i2 < iEnd2) {
				int r1 = bVsp1[i1].r, r2 = bVsp2[i2].r;
				if (r1 < r2) {									// was pivot1's column element before pivot2's column element (mismatch)?
					if (--zAllowedC < 0) return null;			// if we ran out of allowed column fill-in elements, fail the supernode matching
					data[iTBC++] = 0;
					i1++;
				} else if (r1 > r2) {							// was pivot1's column element after pivot2's column element (mismatch)?
					if (--zAllowedC < 0) return null;			// if we ran out of allowed column fill-in elements, fail the supernode matching
					data[iTBC++] = bVsp2[i2++].v;				// otherwise, advance in pivot2's sparse column, copying over to frontal matrix
				} else {
					data[iTBC++] = bVsp2[i2++].v;
					i1++;
				}
			}
			if (iEnd2 - i2 > zAllowedC) return null;			// if remaining elements of pivot2 exceed remaining allowed fill-ins, fail
			while (i2 < iEnd2) data[iTBC++] = bVsp2[i2++].v;	// otherwise copy them over
			fm.nPivots++;
			fm.nContribC--;										// the new pivot row & column has taken away a contributive row & column
			fm.nContribR--;

			// pattern matching tests succeeded, move over data into submatrix from temporary buffer
			for (int iTBR2 = fm.iTB; iTBR2 < iTBR; iTBR2++)	data[iMH2++] = data[iTBR2];
			if (iTBR < iTBC - 1)								// if a single value only has been copied (end of frontal matrix), this loop is redundant
				for (int iTBC2 = iTBR; iTBC2 < iTBC; iTBC2++)	{ data[iMV2] = data[iTBC2]; iMV2 += n; }
		}
		
		return fm;											// return frontal matrix with added elements
	}
	
	
	// method takes a list of expected contribution row & column indexes and index offsets and creates gaps for the
	// contribution indexes that are not present in the absorbing frontal matrix
	public void extend(int[][] gExtraP, int[][] gContribC, int[][] gContribR, int[][] gExtraC, int[][] gExtraR) {

		// get matrix's new extended size
		double[] data2 = null;
		int[] idxCR2 = null;
		int nExtraP3 = gExtraP[1]==null ? 0 : gExtraP[1][0], nExtraP2 = nExtraP3 + nExtraP;
		int nContribC2 = gContribC[1][0], nContribR2 = gContribR[1][0];
		int nExtraC2 = gExtraC[1]==null ? 0 : gExtraC[1][0], nExtraR2 = gExtraR[1]==null ? 0 : gExtraR[1][0];
		int c2 = nExtraP2 + nContribC2 + nExtraC2;							// find out expected size of frontal
		int r2 = nExtraP2 + nContribR2 + nExtraR2;							// note: the pivots are included in nContribC2 & nContribR2
		//if (cTot < c2 || rTot < r2) {										// see if it fits in the original data array
			int cTot2, rTot2;
			if (EXTRA_FM_SPACE) { cTot2 = c2 + c2/4; rTot2 = r2 + r2/4; }	// cTot2 & rTot2 is the new total rank of this frontal
			else { cTot2 = c2; rTot2 = r2; }
			data2 = new double[cTot2 * rTot2 + (rTot2 > cTot2 ? rTot2 : cTot2)*2];
			idxCR2 = new int[cTot2 + rTot2 + 1];
			mixCode = new byte[rTot2 > cTot2 ? rTot2 + 1: cTot2 + 1];
			//------ DEBUG matrix viewer ------//
			FrontalMatrix fm2 = new FrontalMatrix(this.pivot, r2, c2, nContribR2, nContribC2, nExtraR2, nExtraC2);
			fm2.data = data2; fm2.idxCR = idxCR2; fm2.nPivots = nPivots; fm2.nExtraP = nExtraP2;
			//------ DEBUG matrix viewer ------//

			// initiate row-wise offsets for writing data into original and resized frontal
			int jMH1 = nExtraP * (cTot + 1), jMH2 = nExtraP2 * (cTot2 + 1);
			int[] offsetC = gContribC[1];
			// copy from original rows: pivot block + contribution rows
			for (int iP = 0; iP < nPivots; iP++) {									// for every pivot row in original frontal
				int j1 = 0, jMH1b = jMH1 + cTot * iP, jMH2b = jMH2 + cTot2 * iP;
				while (j1 < nPivots) { data2[jMH2b++] = data[jMH1b++]; j1++; }		// copy over pivot block row part

				if (offsetC.length > FrontalDAG.HELPARRAY_SIZE) {
					int jEnd = nPivots + nContribC, jO = FrontalDAG.HELPARRAY_SIZE;	// copy over contribution part with offsets from analysis	
					while (j1 < jEnd) { data2[jMH2b + offsetC[jO++]] = data[jMH1b++]; j1++; }	
				} else {
					int jEnd = nPivots + nContribC;									// there was no contribution, just straight copy	
					while (j1 < jEnd) { data2[jMH2b++] = data[jMH1b++]; j1++; }						
				}
			}
			
			// configure the new idxCR fields: failed pivots+extra pivots, pivots, contributions, extras
			
			// adding unsorted pivot indexes that possibly have disassociated row/column index because of pivoting across frontal matrices
			int ix2 = 0, ix2R = cTot;
			for (int ix1 = 0, ix2End = nExtraP3; ix2 < ix2End; ix1 += 2, ix2++, ix2R++) {
				idxCR2[ix2R] = gExtraP[0][ix1++];			// this is the possibly disassociated row index
				idxCR2[ix2] = gExtraP[0][ix1];				// this is the possibly disassociated column index
			}

			// copy over original extrapivot+pivot indexes into extended idxCR
			for (int ix1 = nExtraP, ix1End = nExtraP + nPivots; ix1 < ix1End; ix1++, ix2++)
				idxCR2[ix2] = idxCR[ix1];
			// copy over the intermixed contribution column indexes into extended idxCR
			for (int ix1 = nPivots, ix2End = gContribC[1][0]; ix2 < ix2End; ix1++, ix2++)
				idxCR2[ix2] = gContribC[0][ix1];
			// copy over the intermixed extras column indexes into extended idxCR
			for (int ix1 = 0, ix2End = nExtraC2; ix2 < ix2End; ix1++, ix2++)
				idxCR2[ix2] = gExtraC[0][ix1];
			
			// initiate column-wise offsets for writing data into original and resized frontal
			int iOffset = nExtraP + nPivots, cTotM1 = cTot-1, cTotM2 = cTot-2, cTotM3 = cTot-3;
			int iMH1 = iOffset * cTot + nExtraP, iMH2 = 0;
			int[] offsetR = gContribR[1];
			// copy from original columns: contribution columns
			if (offsetR.length > FrontalDAG.HELPARRAY_SIZE)
				for (int jP = nPivots; jP > 0;) {											// for every pivot column in original frontal
					if (jP >= 4) {															// unroll inner loop for up to 4 concurrent pivots
						for (int i1 = 0, iO = FrontalDAG.HELPARRAY_SIZE + nPivots; i1 < nContribR; i1++) {
							iMH2 = offsetR[iO++] * cTot2 + nExtraP2;
							data2[iMH2++] = data[iMH1++]; data2[iMH2++] = data[iMH1++];
							data2[iMH2++] = data[iMH1++]; data2[iMH2] = data[iMH1]; iMH1 += cTotM3; }
						jP -= 4;
					} else if (jP == 3) {
						for (int i1 = 0, iO = FrontalDAG.HELPARRAY_SIZE + nPivots; i1 < nContribR; i1++) {
							iMH2 = offsetR[iO++] * cTot2 + nExtraP2;
							data2[iMH2++] = data[iMH1++]; data2[iMH2++] = data[iMH1++];
							data2[iMH2] = data[iMH1]; iMH1 += cTotM2; }
						jP -= 3;
					} else if (jP == 2) {				
						for (int i1 = 0, iO = FrontalDAG.HELPARRAY_SIZE + nPivots; i1 < nContribR; i1++) {
							iMH2 = offsetR[iO++] * cTot2 + nExtraP2;
							data2[iMH2++] = data[iMH1++]; data2[iMH2] = data[iMH1]; iMH1 += cTotM1; }
						jP -= 2;
					} else {
						for (int i1 = 0, iO = FrontalDAG.HELPARRAY_SIZE + nPivots; i1 < nContribR; i1++) {
							iMH2 = offsetR[iO++] * cTot2 + nExtraP2; data2[iMH2] = data[iMH1]; iMH1 += cTot; }
						jP--;
					}					
				}
			// if no horisontal extension happened, do a straight row copy
			else
				for (int jP = nPivots; jP > 0;) {											// for every pivot column in original frontal
					if (jP >= 4) {															// unroll inner loop for up to 4 concurrent pivots
						for (int i1 = 0; i1 < nContribR; i1++) { iMH2 = i1 * cTot2 + nExtraP2;
							data2[iMH2++] = data[iMH1++]; data2[iMH2++] = data[iMH1++];
							data2[iMH2++] = data[iMH1++]; data2[iMH2] = data[iMH1]; iMH1 += cTotM3; }
						jP -= 4;
					} else if (jP == 3) {
						for (int i1 = 0; i1 < nContribR; i1++) { iMH2 = i1 * cTot2 + nExtraP2;
							data2[iMH2++] = data[iMH1++]; data2[iMH2++] = data[iMH1++];
							data2[iMH2] = data[iMH1]; iMH1 += cTotM2; }
						jP -= 3;
					} else if (jP == 2) {				
						for (int i1 = 0; i1 < nContribR; i1++) { iMH2 = i1 * cTot2 + nExtraP2;
							data2[iMH2++] = data[iMH1++]; data2[iMH2] = data[iMH1]; iMH1 += cTotM1; }
						jP -= 2;
					} else {
						for (int i1 = 0; i1 < nContribR; i1++) {
							iMH2 = i1 * cTot2 + nExtraP2; data2[iMH2] = data[iMH1]; iMH1 += cTot; }
						jP--;
					}					
				}
							
			// configure the new idxCR fields: failed pivots+extra pivots, pivots, contributions, extras
			
			// copy over original pivot indexes into extended idxCR
			ix2 = cTot2 + nExtraP3;
			for (int ix1 = cTot + nExtraP, ix1End = ix1 + nPivots; ix1 < ix1End; ix1++, ix2++)
				idxCR2[ix2] = idxCR[ix1];
			// copy over the intermixed contribution row indexes into extended idxCR
			for (int ix1 = nPivots, ix2End = ix2 + nContribR2 - nPivots; ix2 < ix2End; ix1++, ix2++)
				idxCR2[ix2] = gContribR[0][ix1];
			// copy over the intermixed extras row indexes into extended idxCR
			for (int ix1 = 0, ix2End = ix2 + nExtraR2; ix2 < ix2End; ix1++, ix2++)
				idxCR2[ix2] = gExtraR[0][ix1];

			cTot = cTot2; rTot = rTot2;
//		} else {
//			data2 = data;
//			idxCR2 = idxCR;
//		}
		
		// alter parameters of this frontal
		data = data2; idxCR = idxCR2;
		extraC = gExtraC; extraR = gExtraR;
		//extraP = gExtraP[0];
		c = c2; r = r2;
		nExtraP = nExtraP2;			// nExtraP2 is parent's extra + childrens failed pivots, nExtraP3 is just the children's failed pivots
		nContribC = nContribC2 - nPivots; nContribR = nContribR2 - nPivots;
		nExtraC = nExtraC2; nExtraR = nExtraR2;
	}
	
	
	
	
	// method adds in child's contribution to current frontal
	// DEBUG: add() does not use mixCodes
	public void add(FrontalMatrix fmC) {
		
		double[] dataC = fmC.data;
		int[] idxCRC = fmC.idxCR;
		switch (fmC.childFlag) {
		case 1: {	// (1) L-child case, child's contrib.columns added to parent's pivot+contrib.columns

			int cTotC = fmC.cTot, skipExPiC = fmC.nExtraP + fmC.nPivots, skipExPiP = nExtraP + nPivots;
			int jC = skipExPiC, jP = nExtraP;
			int iCEnd = cTotC + skipExPiC + fmC.nContribR, iPEnd = cTot + skipExPiP + nContribR;
			int jCEnd = jC + fmC.nContribC, jPEnd = skipExPiP + nContribC;

			int jdxC = idxCRC[jC], jdxP = idxCR[jP];
			while (jC < jCEnd && jP < jPEnd) {
				if (jdxC < jdxP) jdxC = idxCRC[++jC];								// child's column index lower than parents
				else if (jdxC > jdxP) jdxP = idxCR[++jP];							// child's column index higher than parents
				else {																// columns matching
					int idC = cTotC * skipExPiC + jC, idP = cTot * nExtraP + jP;
					int iC = cTotC + skipExPiC, iP = cTot + nExtraP;

					int idxC = idxCRC[iC], idxP = idxCR[iP];						// get the first row indexes of iteration	
					while (iC < iCEnd && iP < iPEnd) {								// add in those elements of row that match parent's column indexes
						if (idxC < idxP) { idxC = idxCRC[++iC]; idC += cTotC; }		// child's column index lower than parents
						else if (idxC > idxP) { idxP = idxCR[++iP]; idP += cTot; }	// child's column index higher than parents
						else {	
							data[idP] += dataC[idC];
							// a child can have more than one parent, and since every contribution must only happen once, we're zeroing this one
							dataC[idC] = 0;
							idC += cTotC; idP += cTot;
							idxC = idxCRC[++iC]; idxP = idxCR[++iP];				// advance in both row index arrays
						}
					}
					jdxC = idxCRC[++jC]; jdxP = idxCR[++jP];						// advance in both column index arrays
				}
			}
			
			// (2) for L-children, the columns of failed pivots are added to parent's extra columns
			if (fmC.nFailed > 0) {
				
				int[] failedP = fmC.failedP;
				int skipExPiCoP = skipExPiP + nContribC;
				iCEnd = cTotC + skipExPiC + fmC.nContribR; //iPEnd = cTot + skipExPiP + nContribR;
				int jExC = 0, jExP = skipExPiCoP, jExCEnd = fmC.nFailed, jExPEnd = jExP + nExtraC;
								
				// adds columns of failed pivots from an unsorted order triplet list
				int jExC3 = 1;													// we're comparing with column index in child's failedP
				while (jExC < jExCEnd && jExP < jExPEnd) {
					if (failedP[jExC3] == idxCR[jExP]) {						// check if child's failed global column == parent's extras global column
						int idC = failedP[jExC3 + 1], idP = jExP;				// offsets to child's failed column and parent's matching extra column
						int iC = cTotC + idC, iP = cTot;						// offsets to idxCR's frontal row indices of child & parent
						idC = idC * (cTotC + 1);
						
						// start adding contribution along column until end of contribution field
						int idxC = idxCRC[iC], idxP = idxCR[iP];
						while (iC < iCEnd /*&& iP < iPEnd*/) {
							if (idxC < idxP) { idxC = idxCRC[++iC]; idC += cTotC; }			// child's column index lower than parents
							else if (idxC > idxP) { idxP = idxCR[++iP]; idP += cTot; }		// child's column index higher than parents
							else {	
								data[idP] += dataC[idC];
								dataC[idC] = 0;									// every contribution must only happen once, so zero it out
								idC += cTotC; idP += cTot;						// step down both the matching columns
								idxC = idxCRC[++iC]; idxP = idxCR[++iP];		// advance in both index arrays
							}
						}										
						jExC++; jExP++; jExC3 += FAILED_VECSIZE;
						
					} else {
						jExC3 += FAILED_VECSIZE;
						if (++jExC >= jExCEnd) { jExC = 0; jExC3 = 1; }			// if we walked through entire failedP of child, restart
						else jExC3 += FAILED_VECSIZE;							// and continue doing this until parent's extraP is exhausted
					}
				}
				
			}
		} break;
		
		case 2: {	// (1) U-child case, child's contrib.rows added to parent's pivot+contrib.rows
			
			int cTotC = fmC.cTot, skipExPiC = fmC.nExtraP + fmC.nPivots, skipExPiP = nExtraP + nPivots;
			int iC = cTotC + skipExPiC, iP = cTot + nExtraP;
			int jCEnd = skipExPiC + fmC.nContribC, jPEnd = skipExPiP + nContribC;
			int iCEnd = iC + fmC.nContribR, iPEnd = cTot + skipExPiP + nContribR;

			// search for row indexes that match for child and parent
			int idxC = idxCRC[iC], idxP = idxCR[iP];
			while (iC < iCEnd && iP < iPEnd) {
				if (idxC < idxP) idxC = idxCRC[++iC];						// child's row index lower than parents
				else if (idxC > idxP) idxP = idxCR[++iP];					// child's row index higher than parents
				else {														// rows matching
					int jdC = cTotC * (iC - cTotC) + skipExPiC, jdP = cTot * (iP - cTot) + nExtraP;
					int jC = skipExPiC, jP = nExtraP;

					int jdxC = idxCRC[jC], jdxP = idxCR[jP];
					while (jC < jCEnd && jP < jPEnd) {						// add in those elements of row that match parent's column indexes
						if (jdxC < jdxP) { jdxC = idxCRC[++jC]; jdC++; }	// child's column index lower than parents
						else if (jdxC > jdxP) { jdxP = idxCR[++jP]; jdP++; }// child's column index higher than parents
						else {	
							data[jdP++] += dataC[jdC];
							dataC[jdC++] = 0;								// every contribution must only happen once, so zero it out
							jdxC = idxCRC[++jC]; jdxP = idxCR[++jP];		// advance in both index arrays
						}
					}
					idxC = idxCRC[++iC]; idxP = idxCR[++iP];
				}
			}
		
			// (2) for U-children, the rows of failed pivots are added to parent's extra rows		
			if (fmC.nFailed > 0) {
				
				int[] failedP = fmC.failedP;
				jCEnd = skipExPiC + fmC.nContribC; jPEnd = skipExPiP + nContribC; 
				int skipExPiCoP = skipExPiP + nContribC;
				int iExC = 0, iExP = cTot + skipExPiCoP, iExCEnd = fmC.nFailed, iExPEnd = iExP + nExtraR;
	
				// adds columns of failed rows from an unsorted order triplet list
				int iExC3 = 0;													// we're comparing with row index in child's failedP
				while (iExC < iExCEnd && iExP < iExPEnd) {
					if (failedP[iExC3] == idxCR[iExP]) {
						int jC = failedP[iExC3 + 2], jP = 0;					// offsets to idxCR's frontal row indices of child & parent
						int jdC = jC * (cTotC + 1), jdP = cTot * (iExP - cTot);	// offsets to child's failed column and parent's matching extra column
						
						// start adding contribution along column until end of contribution field
						int jdxC = idxCRC[jC], jdxP = idxCR[jP];
						while (jC < jCEnd && jP < jPEnd) {
							if (jdxC < jdxP) { jdxC = idxCRC[++jC]; jdC++; }	// child's column index lower than parents
							else if (jdxC > jdxP) { jdxP = idxCR[++jP]; jdP++; }// child's column index higher than parents
							else {	
								data[jdP++] += dataC[jdC];
								data[jdC++] = 0;								// every contribution must only happen once, so zero it out
								jdxC = idxCRC[++jC]; jdxP = idxCR[++jP];		// advance in both index arrays
							}
						}
											
						iExP++; iExC++; iExC3 += FAILED_VECSIZE;
					} else {
						iExC3 += FAILED_VECSIZE;								// unsorted array scanning case:
						if (++iExC >= iExCEnd) { iExC = 0; iExC3 = 0; }			// if we walked through entire failedP of child, restart
						else iExC3 += FAILED_VECSIZE;							// and continue doing this until parent's extraP is exhausted
					}
				}

			}
		} break;
		
		case 3:	{	// LU-child case, we must add failed pivots to parent's extra pivots
			
			if (fmC.nFailed > 0) {
				// start by adding the extra pivots, iterate trough every extra pivot index
				int[] extraPC = fmC.failedP;
	
				// add failed pivots from an unsorted order triplet list
				for (int jP = nExtraP2, jC = 0, jC3 = 0, jCEnd = fmC.nFailed; jC < jCEnd; jC++, jC3 += FAILED_VECSIZE, jP++) {
					// copy over the row into parent's extra pivots
					int cTotC = fmC.cTot, jC2 = extraPC[jC3 + 2], jP2 = jP;
					int jdC = jC2 * (cTotC + 1);						// get offset into child frontal to failed pivot position
					int jdP = jP2 * (cTot + 1), idC = jdC, idP = jdP;	// get offset into parent frontal to expected extra pivot position
					int jC2End = fmC.c; jC2++;
					data[jdP++] = fmC.data[jdC++];						// copy over failed pivot value
					idxCR[cTot + jP2] = extraPC[jC3];					// copy over failed pivot's row index
					idxCR[jP2++] = extraPC[jC3 + 1];					// copy over failed pivot's column index
	
					// copy over the row into parent's extra pivots
					int idxP = idxCR[jP2], idxC = idxCRC[jC2];
					while (jP2 < c && jC2 < jC2End) {					// copy over failed pivot row only into matching column indexes of parent
						if (idxP < idxC)			{ jdP++; idxP = idxCR[++jP2]; }
						else if (idxP > idxC)	{ jdC++; idxC = idxCRC[++jC2]; }
						else { 	data[jdP++] = fmC.data[jdC++];
								idxP = idxCR[++jP2]; idxC = idxCRC[++jC2]; }
					}
					
					// copy over the column into parent's extra pivots
					int iP = cTot + jP + 1, iC = cTotC + extraPC[jC3 + 2] + 1, iPEnd = cTot + r, iCEnd = cTotC + fmC.r;
					idP += cTot; idC += cTotC;						// skip over pivot value, it's already copied
	
					idxP = idxCR[iP]; idxC = idxCRC[iC];
					while (iP < iPEnd && iC < iCEnd) {					// copy over failed pivot column only into matching row indexes of parent
						if (idxP < idxC)		{ idxP = idxCR[++iP]; idP += cTot; }
						else if (idxP > idxC)	{ idxC = idxCRC[++iC]; idC += cTotC; }
						else {
							data[idP] = fmC.data[idC];
							idP += cTot; idC += cTotC;
							idxP = idxCR[++iP]; idxC = idxCRC[++iC];
						}
					}
				}
				nExtraP2 += fmC.nFailed;								// update temporary nFailed count for consecutive unsorted additions
			}

			// for LU-child, iterate through child's entire contrib.part, index-matching and moving values
			int cTotC = fmC.cTot,  skipExPiC = fmC.nExtraP + fmC.nPivots;
			int iCEnd = cTotC + fmC.r, iPEnd = cTot + c;
			int iC = cTotC + skipExPiC, iP = nExtraP + cTot;
			int idxP = idxCR[iP], idxC = idxCRC[iC];
			for (; iC < iCEnd && iP < iPEnd;) {
				if (idxP < idxC)		idxP = idxCR[++iP];
				else if (idxP > idxC)	idxC = idxCRC[++iC];
				else {
					int jCEnd = fmC.c, idC = (iC - cTotC) * cTotC, idP = (iP - cTot) * cTot;
					int jdxP = idxCR[nExtraP], jdxC = idxCRC[skipExPiC];
					for (int jC = skipExPiC, jP = nExtraP; jC < jCEnd && jP < c;) {
						if (jdxP < jdxC)		jdxP = idxCR[++jP];
						else if (jdxP > jdxC)	jdxC = idxCRC[++jC];
						else { 
							data[idP + jP] += dataC[idC + jC];
							dataC[idC + jC] = 0;
							jdxP = idxCR[++jP]; jdxC = idxCRC[++jC];
						}
					}
					idxP = idxCR[++iP]; idxC = idxCRC[++iC];
				}
			}			
		}
		}
	}

	
	// DEBUG: add2() uses the mixCodes optimisation
	public void add2(FrontalMatrix fmC) {
		
		double[] dataC = fmC.data;
		int[] idxCRC = fmC.idxCR;
		
		int cTotC = fmC.cTot,  skipExPiC = fmC.nExtraP + fmC.nPivots, skipExPiP = nExtraP + nPivots;
		int iCEnd = cTotC + skipExPiC + fmC.nContribR, iPEnd = cTot + r;
		int iC = cTotC + skipExPiC, iP = nExtraP + cTot;
		int idxP = idxCR[iP], idxC = idxCRC[iC];
	
		// the contribution adding part is the same for L/U/LU-children, since the indexes have been joined correctly by factoriseDAG()

		// optimisation: if already 1st pivot of child failed, then no factorisation was made, and if it's a leaf frontal, then contrib.block = 0
		if (!(fmC.firstPivotFailed() == true && fmC.isLeaf())) {

			boolean execMixCodes = false;						// the mixCodes need to be accumulated during first row addition, thus still false

			while (iC < iCEnd && iP < iPEnd) {
				if (idxP < idxC)		idxP = idxCR[++iP]; 	// mismatch: parent's row lower than child's
				else if (idxP > idxC)	idxC = idxCRC[++iC];	// mismatch: child's row lower than parents
				else {											// rows matched and can be added
					if (execMixCodes) {
						mixRowData(dataC, data, skipExPiC + (iC - cTotC) * cTotC, nExtraP + (iP - cTot) * cTot);
					} else {
						int jCEnd = skipExPiC + fmC.nContribC, idC = (iC - cTotC) * cTotC, idP = (iP - cTot) * cTot;
						int jdxP = idxCR[nExtraP], jdxC = idxCRC[skipExPiC];
						int mc = 0, mc3 = 3, code = 0;
						for (int jC = skipExPiC, jP = nExtraP; jC < jCEnd && jP < c;) {
							if (jdxP < jdxC)		{ jdxP = idxCR[++jP]; code = (code<<2) | 1; mc3--; }
							else if (jdxP > jdxC)	{ jdxC = idxCRC[++jC]; code <<= 2; mc3--; }
							else { 
								data[idP + jP] += dataC[idC + jC];
								dataC[idC + jC] = 0;
								jdxP = idxCR[++jP]; jdxC = idxCRC[++jC];
								code = (code<<2) | 2; mc3--;
							}
							if (mc3 == 0) { mixCode[mc++] = (byte)code; mc3 = 3; code = 0; }
						}
						if (mc3 > 0) mixCode[mc++] = (byte)(code << (mc3 * 2));
						mixCode[mc++] = -1;
						execMixCodes = true;
					}
					idxP = idxCR[++iP]; idxC = idxCRC[++iC];
				}
			}
		}

		// take care of the transfers specific for the L/U/LU-children
		switch (fmC.childFlag & 0xF) {
		case 1: {
			
			// (1) for L-children, the columns of failed pivots are added to parent's extra columns
			// note: the footwork of arranging the correct indexing has already been done by either extend() or
			// FrontalMatrix() methods, so these indexes are used to find what columns should be moved
			if (fmC.nFailed > 0) {
				
				int[] failedP = fmC.failedP;
				int skipExPiCoP = skipExPiP + nContribC;
				iCEnd = cTotC + skipExPiC + fmC.nContribR; //iPEnd = cTot + skipExPiP + nContribR;
				int jExC = 0, jExP = skipExPiCoP, jExCEnd = fmC.nFailed;//, jExPEnd = jExP + nExtraC;
								
				// adds columns of failed pivots from an unsorted order duplet list
				int jExC3 = 1;													// we're comparing with column index in child's failedP
				while (jExC < jExCEnd /*&& jExP < jExPEnd*/) {
					if (failedP[jExC3] == idxCR[jExP]) {						// check if child's failed global column == parent's extras global column
						iC = cTotC + failedP[jExC3 + 1]; iP = cTot;				// offsets to idxCR's frontal row indices of child & parent
						int idC = failedP[jExC3 + 1], idP = jExP;				// offsets to child's failed column and parent's matching extra column
						idC = idC * (cTotC + 1);
						
						// start adding contribution along column until end of contribution field
						idxC = idxCRC[iC]; idxP = idxCR[iP];
						while (iC < iCEnd /*&& iP < iPEnd*/) {
							if (idxC < idxP) { idxC = idxCRC[++iC]; idC += cTotC; }			// child's column index lower than parents
							else if (idxC > idxP) { idxP = idxCR[++iP]; idP += cTot; }		// child's column index higher than parents
							else {	
								data[idP] += dataC[idC];
								dataC[idC] = 0;									// every contribution must only happen once, so zero it out
								idC += cTotC; idP += cTot;						// step down both the matching columns
								idxC = idxCRC[++iC]; idxP = idxCR[++iP];		// advance in both index arrays
							}
						}										
						failedP[jExC3] = -1;									// "deactivate" child's column for all other potential parents of child
						jExC++; jExP++; jExC3 += FAILED_VECSIZE;
						
					} else {
						jExC3 += FAILED_VECSIZE;
						if (++jExC >= jExCEnd) { jExC = 0; jExC3 = 1; }			// if we walked through entire failedP of child, restart
						else jExC3 += FAILED_VECSIZE;							// and continue doing this until parent's extraP is exhausted
					}
				}
				
			}
		} break;
		
		case 2: {
			
			// (2) for U-children, the rows of failed pivots are added to parent's extra rows		
			if (fmC.nFailed > 0) {
				
				int[] failedP = fmC.failedP;
				int jCEnd = skipExPiC + fmC.nContribC; //jPEnd = skipExPiP + nContribC; 
				int skipExPiCoP = skipExPiP + nContribC;
				int iExC = 0, iExP = cTot + skipExPiCoP, iExCEnd = fmC.nFailed; //, iExPEnd = iExP + nExtraR;
	
				// adds columns of failed rows from an unsorted order triplet list
				int iExC3 = 0;													// we're comparing with row index in child's failedP
				while (iExC < iExCEnd /*&& iExP < iExPEnd*/) {
					if (failedP[iExC3] == idxCR[iExP]) {
						int jC = failedP[iExC3 + 2], jP = 0;					// offsets to idxCR's frontal row indices of child & parent
						int jdC = jC * (cTotC + 1), jdP = cTot * (iExP - cTot);	// offsets to child's failed column and parent's matching extra column
						
						// start adding contribution along column until end of contribution field
						int jdxC = idxCRC[jC], jdxP = idxCR[jP];
						while (jC < jCEnd /*&& jP < jPEnd*/) {
							if (jdxC < jdxP) { jdxC = idxCRC[++jC]; jdC++; }	// child's column index lower than parents
							else if (jdxC > jdxP) { jdxP = idxCR[++jP]; jdP++; }// child's column index higher than parents
							else {	
								data[jdP++] += dataC[jdC];
								data[jdC++] = 0;								// every contribution must only happen once, so zero it out
								jdxC = idxCRC[++jC]; jdxP = idxCR[++jP];		// advance in both index arrays
							}
						}
						failedP[iExC3] = -1;									// "deactivate" child's row for all other potential parents of child											
						iExP++; iExC++; iExC3 += FAILED_VECSIZE;
						
					} else {
						iExC3 += FAILED_VECSIZE;								// unsorted array scanning case:
						if (++iExC >= iExCEnd) { iExC = 0; iExC3 = 0; }			// if we walked through entire failedP of child, restart
						else iExC3 += FAILED_VECSIZE;							// and continue doing this until parent's extraP is exhausted
					}
				}

			}
		} break;
		
		case 3:	{	// LU-child case, we must add failed pivots to parent's extra pivots
			
			if (fmC.nFailed > 0) {
				// start by adding the extra pivots, iterate trough every extra pivot index
				int[] extraPC = fmC.failedP;
	
				// add failed pivots from an unsorted order triplet list
				for (int jP = nExtraP2, jC = 0, jC3 = 0, jCEnd = fmC.nFailed; jC < jCEnd; jC++, jC3 += FAILED_VECSIZE, jP++) {
					// copy over the row into parent's extra pivots
					int jC2 = extraPC[jC3 + 2], jP2 = jP;
					int jdC = jC2 * (cTotC + 1);						// get offset into child frontal to failed pivot position
					int jdP = jP2 * (cTot + 1), idC = jdC, idP = jdP;	// get offset into parent frontal to expected extra pivot position
					int jC2End = fmC.c; jC2++;
					data[jdP++] = fmC.data[jdC++];						// copy over failed pivot value
					idxCR[cTot + jP2] = extraPC[jC3];					// copy over failed pivot's row index
					idxCR[jP2++] = extraPC[jC3 + 1];					// copy over failed pivot's column index
	
					// copy over the row into parent's extra pivots
					idxP = idxCR[jP2]; idxC = idxCRC[jC2];
					while (/*jP2 < c &&*/ jC2 < jC2End) {				// copy over failed pivot row only into matching column indexes of parent
						if (idxP < idxC)			{ jdP++; idxP = idxCR[++jP2]; }
						else if (idxP > idxC)	{ jdC++; idxC = idxCRC[++jC2]; }
						else { 	data[jdP++] = fmC.data[jdC++];
								idxP = idxCR[++jP2]; idxC = idxCRC[++jC2]; }
					}
					
					// copy over the column into parent's extra pivots
					iP = cTot + jP + 1; iC = cTotC + extraPC[jC3 + 2] + 1;
					idP += cTot; idC += cTotC;						// skip over pivot value, it's already copied
	
					idxP = idxCR[iP]; idxC = idxCRC[iC];
					while (iC < iCEnd /*&& iP < iPEnd*/) {				// copy over failed pivot column only into matching row indexes of parent
						if (idxP < idxC)		{ idxP = idxCR[++iP]; idP += cTot; }
						else if (idxP > idxC)	{ idxC = idxCRC[++iC]; idC += cTotC; }
						else {
							data[idP] = fmC.data[idC];
							idP += cTot; idC += cTotC;
							idxP = idxCR[++iP]; idxC = idxCRC[++iC];
						}
					}
				}
				nExtraP2 += fmC.nFailed;								// update temporary nFailed count for consecutive unsorted additions
			}	
		} break;
		}
	}

	
	
	
	public void decomposeLU(boolean finalParent) {
		double vPivot;
		for (int i = 0, iTot = nExtraP + nPivots; i < iTot; i++) {		// iterate over every pivot (including extra pivots)
			//if (vPivot <= pivotFailFactor) {							// if current pivot is below fail factor		
			if (intraPivoting(i, false)) {								// check numerical stability of current pivot, do pivoting if necessary
				vPivot = data[i * cTot + i];
			} else if (finalParent) {									// for the subproblem's top vertex, we cannot move on failed pivots to parents
				vPivot = Matrix.pivotValid;								// insert a minimal valid pivot value
			} else {
				// intrapivoting failed: store down all consecutive pivots global & local indexes
				failedP = new int[(nFailed = iTot - i) * FAILED_VECSIZE];
				for (int i2 = i, i3 = 0; i2 < iTot; i2++) {
					failedP[i3++] = idxCR[cTot + i2];					// store global row index of failed pivot
					failedP[i3++] = idxCR[i2];							// store global column index of failed pivot
					failedP[i3++] = i2;									// store offset to failed index
				}
				
				// this case can't happen if everything was correctly initialised
				//if (i == 0) firstPivotFailure();						// if failure was on already the first pivot, flag this
				return;													// all consecutive candidates failed, we're done
			}

			// pivot okay, do the factorising
			vPivot = 1.0 / vPivot;
			int jEnd = cTot * r, i2 = i * cTot;
			for (int j = i * (cTot + 1) + cTot; j < jEnd; j += cTot)	// divides column i of L
				data[j] *= vPivot;
			for (int k = i + 1; k < c; k++)	{							// updates rest of matrix
				double di2k = data[i2 + k];
				if (!Matrix.nearZero(di2k))
					for (int j = (i + 1) * cTot; j < jEnd; j += cTot)
						data[j + k] -= data[j + i] * di2k;
			}
		}
	}

	
	
	// method finds a better candidate for pivot p in frontal matrix, swapping in the row & column of the candidate, otherwise failing
	// method also readjusts the indexation of global indexes in idxRC
	// DEBUG: compareWithColumn will flag for comparing current pivot with it's ENTIRE column, and not just within the pivot block
	public boolean intraPivoting(int p, boolean compareWithColumn) {
		
		// iterate over frontal matrix pivot block in search of a better pivot candidate
		double pivot = data[p * (cTot + 1)];
		if (pivot < 0) pivot = - pivot;
		double vMax = pivot / Matrix.pivotCritFactor;
		int p2 = p + 1, iPivots = nExtraP + nPivots;
		int idT = p * (cTot + 1) + cTot, iTEnd = (compareWithColumn ? r : nPivots) * cTot, iMax = -1;
		
		// investigate pivot's column for threshold pivoting criterion
		for (int idT2 = idT, iT = p2; idT2 < iTEnd; idT2 += cTot, iT++)
			if (data[idT2] < -vMax || data[idT2] > vMax ) {					// was a[i,i] < u*max(a[i,k], k>=i), a higher value found?
				vMax = data[idT2] < 0 ? -data[idT2] : data[idT2];
				iMax = iT;													// store down local row index of max value
			}
		
		// pivot failed column threshold test, or pivot is below fail factor
		if (iMax >= 0 || pivot <= Matrix.pivotFailFactor) {
			// if a superior column value was found, check if it's located in the pivoting block and can be used
			if (iMax >= 0 && iMax < iPivots) {
				for (int i1 = p * cTot, i2 = iMax * cTot, i1End = i1 + c; i1 < i1End; i1++, i2++) {
					double v = data[i1]; data[i1] = data[i2]; data[i2] = v;	// swap rows
				}
				p += cTot; iMax += cTot;					
//				mutator[0][mutator[0][0] + 1] = idxCR[p];							// the pivot index that was swapped away
//				mutator[0][mutator[0][0] + 2] = idxCR[iMax];						// the pivot index that was swapped in
//				mutator[0][0] += 2;
//				mutator[0][idxCR[p]] = idxCR[iMax];
				int idx = idxCR[p]; idxCR[p] = idxCR[iMax]; idxCR[iMax] = idx;		// swap global row indexes
				return true;
			}		

			// can't swap in column value OR pivot below fail factor, we have to search this pivot's submatrix for best candidate
			vMax = Matrix.pivotFailFactor;
			// DEBUG: this loop checks pivot's submatrix INCLUSIVE the column of the pivot: (i>p & j>=p)
			int jMax = p;
			for (int i = p2, jP = (p + 1) * cTot, jAdd = cTot - iPivots; i < iPivots; i++, jP += jAdd)
				for (int j = p; j < iPivots; j++, jP++)								// compare to all values except the pivot itself
					if (data[jP] < -vMax || data[jP] > vMax ) {						// found better pivot candidate?
						vMax = data[jP] < 0 ? -data[jP] : data[jP];
						iMax = i; jMax = j;											// store down local row & column indexes of max value
					}
			// DEBUG: this loop checks STRICTLY the submatrix of the pivot: (i>p & j>p)
//			int jMax = p;
//			for (int i = p2, jP = p2 * (cTot + 1), jAdd = cTot - iPivots; i < iPivots; i++, jP += jAdd)
//				for (int j = p2; j < iPivots; j++, jP++)							// compare to all values on subrow
//					if (data[jP] < -vMax || data[jP] > vMax ) {						// found better pivot candidate?
//						vMax = data[jP] < 0 ? -data[jP] : data[jP];
//						iMax = i; jMax = j;											// store down local row & column indexes of max value
//					}
			
			if (vMax <= Matrix.pivotFailFactor) return false;						// if nothing above pivot fail factor found, fail at this pivot

			if (p != jMax) {
				for (int j1 = p, j2 = jMax, j1End = r * cTot; j1 < j1End; j1 += cTot, j2 += cTot) {
					double v = data[j1]; data[j1] = data[j2]; data[j2] = v;			// swap the two columns
				}
	
//				mutator[1][mutator[1][0] + 1] = idxCR[p];							// the column swapped away at this position
//				mutator[1][mutator[1][0] + 2] = idxCR[jMax];						// the column swapped in at this position
//				mutator[1][0] += 2;
//				mutator[0][idxCR[p]] = idxCR[jMax];
				int idx = idxCR[p]; idxCR[p] = idxCR[jMax]; idxCR[jMax] = idx;		// swap global column indexes
			}
			if (p != iMax) {
				for (int i1 = p * cTot, i2 = iMax * cTot, i1End = i1 + c; i1 < i1End; i1++, i2++) {
					double v = data[i1]; data[i1] = data[i2]; data[i2] = v;			// swap the two rows
				}
				p += cTot; iMax += cTot;
//				mutator[0][mutator[0][0] + 1] = idxCR[p];							// the column swapped away at this position
//				mutator[0][mutator[0][0] + 2] = idxCR[iMax];						// the column swapped in at this position
//				mutator[0][0] += 2;
//				mutator[0][idxCR[p]] = idxCR[iMax];
				int idx = idxCR[p]; idxCR[p] = idxCR[iMax]; idxCR[iMax] = idx;		// swap global row indexes
			}
			return true;
		}
		
		// pivot didn't fail the column threshold criterion AND wasn't below pivot fail factor
		return true;
	}
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			DATA HANDLING METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	private void extendRowData(NspNode[] dS, double[] dD, int jS, int jD) {
		int m = 0;
		while (true) {
			switch (mixCode[m++]) {
			case 0:		jD+=4;																				// 0 0 0 0
			case 1:		jD+=3; dD[jD++]=dS[jS++].v;															// 0 0 0 1
			case 2:		jD+=2; dD[jD]=dS[jS++].v; jD+=2;													// 0 0 1 0
			case 3:		jD+=2; dD[jD++]=dS[jS++].v; dD[jD++]=dS[jS++].v;									// 0 0 1 1
			case 4:		jD++; dD[jD]=dS[jS++].v; jD+=3;														// 0 1 0 0
			case 5:		jD++; dD[jD]=dS[jS++].v; jD+=2; dD[jD++]=dS[jS++].v;								// 0 1 0 1
			case 6:		jD++; dD[jD++]=dS[jS++].v; dD[jD]=dS[jS++].v; jD+=2;								// 0 1 1 0
			case 7:		jD++; dD[jD++]=dS[jS++].v; dD[jD++]=dS[jS++].v; dD[jD++]=dS[jS++].v;				// 0 1 1 1
			case 8:		dD[jD]=dS[jS++].v; jD+=4;															// 1 0 0 0
			case 9:		dD[jD]=dS[jS++].v; jD+=3; dD[jD++]=dS[jS++].v;										// 1 0 0 1
			case 10:	dD[jD]=dS[jS++].v; jD+=2; dD[jD]=dS[jS++].v; jD+=2;									// 1 0 1 0
			case 11:	dD[jD]=dS[jS++].v; jD+=2; dD[jD++]=dS[jS++].v; dD[jD++]=dS[jS++].v;					// 1 0 1 1
			case 12:	dD[jD++]=dS[jS++].v; dD[jD]=dS[jS++].v; jD+=3;										// 1 1 0 0
			case 13:	dD[jD++]=dS[jS++].v; dD[jD]=dS[jS++].v; jD+=2; dD[jD++]=dS[jS++].v;					// 1 1 0 1
			case 14:	dD[jD++]=dS[jS++].v; dD[jD++]=dS[jS++].v; dD[jD]=dS[jS++].v; jD+=2;					// 1 1 1 0
			case 15:	dD[jD++]=dS[jS++].v; dD[jD++]=dS[jS++].v; dD[jD++]=dS[jS++].v; dD[jD++]=dS[jS++].v;	// 1 1 1 1
			case -1:	return;
			default:	throw new InvalidParameterException("FrontalMatrix.extendRowData(): Invalid mixCode encountered.");
			}
		}
	}

	private void extendColumnData(NspNode[] dS, double[] dD, int jS, int jD, int rTot) {
		int m = 0, rTot2 = rTot * 2;
		while (true) {
			switch (mixCode[m++]) {
			case 0:		jD+=4*rTot;																				// 0 0 0 0
			case 1:		jD+=3*rTot; dD[jD]=dS[jS++].v; jD+=rTot;												// 0 0 0 1
			case 2:		jD+=rTot2; dD[jD]=dS[jS++].v; jD+=rTot2;												// 0 0 1 0
			case 3:		jD+=rTot2; dD[jD]=dS[jS++].v; jD+=rTot; dD[jD]=dS[jS++].v;	jD+=rTot;					// 0 0 1 1
			case 4:		jD+=rTot; dD[jD]=dS[jS++].v; jD+=3*rTot;												// 0 1 0 0
			case 5:		jD+=rTot; dD[jD]=dS[jS++].v; jD+=rTot2; dD[jD]=dS[jS++].v; jD+=rTot;					// 0 1 0 1
			case 6:		jD+=rTot; dD[jD]=dS[jS++].v; jD+=rTot; dD[jD]=dS[jS++].v; jD+=rTot2;					// 0 1 1 0
			case 7:		jD+=rTot; dD[jD]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; jD+=rTot;// 0 1 1 1
			case 8:		dD[jD]=dS[jS++].v; jD+=4*rTot;															// 1 0 0 0
			case 9:		dD[jD]=dS[jS++].v; jD+=3*rTot; dD[jD]=dS[jS++].v; jD+=rTot;								// 1 0 0 1
			case 10:	dD[jD]=dS[jS++].v; jD+=rTot2; dD[jD]=dS[jS++].v; jD+=rTot2;								// 1 0 1 0
			case 11:	dD[jD]=dS[jS++].v; jD+=rTot2; dD[jD]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; jD+=rTot;		// 1 0 1 1
			case 12:	dD[jD]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; jD+=3*rTot;									// 1 1 0 0
			case 13:	dD[jD]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; jD+=rTot2; dD[jD]=dS[jS++].v; jD+=rTot;		// 1 1 0 1
			case 14:	dD[jD]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; jD+=rTot2;			// 1 1 1 0
			case 15:	dD[jD]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; dD[jD+=rTot]=dS[jS++].v; jD+=rTot;	// 1 1 1 1
			case -1:	return;
			default:	throw new InvalidParameterException("FrontalMatrix.extendRowData(): Invalid mixCode encountered.");
			}
		}
	}

	// method is optimising the sparse mixing of child & parent row data by accepting an array of stepcodes from analysing one row
	// and using these codes to execute one of a set of 3-step mixing access patterns, thus avoiding indirect comparisons between
	// child's & parent's indexation, which leads to a speed-up of execution (affirmed by microbanchmarking)
	// method expects that every row has same access pattern as the predicted analysis given in the array of stepcodes
	// codes: 0 = increment child, 1 = increment father, 2 = copy data from child to father	
	// the mixCode is calculated by: code1<<4 + code2<<2 + code3
	private void mixRowData(double[] dC, double[] dP, int jC, int jP) {
		int m = 0;
		while (true) {
			switch (mixCode[m++]) {
			case 42:		dP[jP++]+=dC[jC]; dC[jC++]=0; dP[jP++]+=dC[jC]; dC[jC++]=0; dP[jP++]+=dC[jC]; dC[jC++]=0; break;// 2 2 2
			case 10:		jC++; dP[jP++]+=dC[jC]; dC[jC++]=0; dP[jP++]+=dC[jC]; dC[jC++]=0; break;						// 0 2 2
			case 26:		jP++; dP[jP++]+=dC[jC]; dC[jC++]=0; dP[jP++]+=dC[jC]; dC[jC++]=0; break;						// 1 2 2
			case 34:		dP[jP++]+=dC[jC]; dC[jC]=0; jC+=2; dP[jP++]+=dC[jC]; dC[jC++]=0; break;							// 2 0 2
			case 38:		dP[jP]+=dC[jC]; dC[jC++]=0; jP+=2; dP[jP++]+=dC[jC]; dC[jC++]=0; break;							// 2 1 2
			case 2:			jC+=2; dP[jP++]+=dC[jC]; dC[jC++]=0; break;														// 0 0 2
			case 18:case 6:	jC++; jP++; dP[jP++]+=dC[jC]; dC[jC++]=0; break;												// 1 0 2 & 0 1 2
			case 22:		jP+=2; dP[jP++]+=dC[jC]; dC[jC++]=0; break;														// 1 1 2
			case 40:		dP[jP++]+=dC[jC]; dC[jC++]=0; dP[jP++]+=dC[jC]; dC[jC]=0; jC+=2; break;							// 2 2 0
			case 41:		dP[jP++]+=dC[jC]; dC[jC++]=0; dP[jP]+=dC[jC]; dC[jC++]=0; jP+=2; break;							// 2 2 1
			case 8:			jC++; dP[jP++]+=dC[jC]; dC[jC]=0; jC+=2; break;													// 0 2 0
			case 24:		jP++; dP[jP++]+=dC[jC]; dC[jC]=0; jC+=2; break;													// 1 2 0
			case 9:			jC++; dP[jP]+=dC[jC]; dC[jC++]=0; jP+=2; break;													// 0 2 1
			case 25:		jP++; dP[jP]+=dC[jC]; dC[jC++]=0; jP+=2; break;													// 1 2 1
			case 32:		dP[jP++]+=dC[jC]; dC[jC]=0; jC+=3; break;														// 2 0 0
			case 36:case 33:dP[jP]+=dC[jC]; dC[jC]=0; jP+=2; jC+=2; break;													// 2 1 0 & 2 0 1
			case 37:		dP[jP]+=dC[jC]; dC[jC++]=0; jP+=3; break;														// 2 1 1
			case 0:			jC+=3; break;																					// 0 0 0
			case 16:case 4:case 1:	jP++; jC+=2; break;																		// 1 0 0 & 0 1 0 & 0 0 1
			case 17:case 5:case 20:	jC++; jP+=2; break;																		// 1 0 1 & 0 1 1 & 1 1 0
			case 21:		jP+=3; break;																					// 1 1 1
			case -1:		return;
			default:		throw new InvalidParameterException("FrontalMatrix.mixRowData(): Invalid mixCode encountered.");
			}
		}
	}

	private void mixColumnData(double[] dC, double[] dP, int iC, int iP, int cTotC) {
		int m = 0;
		while (true) {
			switch (mixCode[m++]) {
			case 42:		dP[iP]+=dC[iC]; iP+=cTot; dC[iC]=0; iC+=cTotC; 
							dP[iP]+=dC[iC]; iP+=cTot; dC[iC]=0; iC+=cTotC;
							dP[iP]+=dC[iC]; iP+=cTot; dC[iC]=0; iC+=cTotC; break;											// 2 2 2
			case 10:		iC+= cTotC; dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC;
							dP[iP]+=dC[iC]; dC[iC]=0;  iP+=cTot; iC+=cTotC; break;											// 0 2 2
			case 26:		iP+=cTot; dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC;
							dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC; break;											// 1 2 2
			case 34:		dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC*2;
							dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC; break;											// 2 0 2
			case 38:		dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot*2; iC+=cTotC;
							dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC; break;											// 2 1 2
			case 2:			iC+=cTotC*2; dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC; break;								// 0 0 2
			case 18:case 6:	iP+=cTot; iC+=cTotC; dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC; break;						// 1 0 2 & 0 1 2
			case 22:		iP+=cTot*2; dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC; break;								// 1 1 2
			case 40:		dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC;
							dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC*2; break;											// 2 2 0
			case 41:		dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC;
							dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot*2; iC+=cTotC; break;											// 2 2 1
			case 8:			iC+=cTotC; dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC*2; break;								// 0 2 0
			case 24:		iP+=cTot; dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC*2; break;								// 1 2 0
			case 9:			iC+=cTotC; dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot*2; iC+=cTotC; break;								// 0 2 1
			case 25:		iP+=cTot; dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot*2; iC+=cTotC; break;								// 1 2 1
			case 32:		dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot; iC+=cTotC*3; break;											// 2 0 0
			case 36:case 33:dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot*2; iC+=cTotC*2; break;										// 2 1 0 & 2 0 1
			case 37:		dP[iP]+=dC[iC]; dC[iC]=0; iP+=cTot*3; iC+=cTotC; break;											// 2 1 1
			case 0:			iC+=cTotC*3; break;																				// 0 0 0
			case 16:case 4:case 1:	iP+=cTot; iC+=cTotC*2; break;															// 1 0 0 & 0 1 0 & 0 0 1
			case 17:case 5:case 20:	iP+=cTot*2; iC+=cTotC; break;															// 1 0 1 & 0 1 1 & 1 1 0
			case 21:		iP+=cTot*3; break;																				// 1 1 1
			case -1:		return;
			default:
				throw new InvalidParameterException("FrontalMatrix.mixColumnData(): Invalid mixCode encountered.");
			}
		}
	}

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			HELPER METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	static boolean notProcessed(FrontalMatrix fm) { return ((fm != null && fm.pivot != -2) ? true : false); }
	void flagProcessed() { pivot = -2; }
	public boolean firstPivotFailed() { return (childFlag & (1<<30)) != 0; }
	public void firstPivotFailure() { childFlag |= (1<<30); }
	public boolean isLeaf() { return (childFlag & (1<<29)) != 0; }
	public void setLeaf() { childFlag |= (1<<29); }
	public boolean isFinal() { return (childFlag & (1<<28)) != 0; }
	public void setFinal() { childFlag |= (1<<28); }
	public void setLcontributor() { childFlag |= 1; }
	public void setUcontributor() { childFlag |= 2; }
	public void setLUcontributor() { childFlag |= 3; }
	public void setLacksContribution() { childFlag |= 4; }
	public void clearLacksContribution() { childFlag &= (0xFFFFFFFF - 4); }
	public boolean lacksContribution() { return (childFlag & 4) != 0; }
	
	private static int trimToAllocBlock(int v) { return (v & (0xFFFFFFFF - FM_ALLOCBLOCK + 1)) + FM_ALLOCBLOCK; }
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	String indexToString(int i) {
		if (i < 10) return String.format("%4d   ", i);
		if (i < 100) return String.format("%5d  ", i);
		return String.format("%6d ", i);
	}
	
	@Override
	public String toString() {		
		StringBuilder sb = new StringBuilder();		
		sb.append("Frontal " + (isLeaf() ? "leaf " : "") + "matrix " + pivot + ", " + r + "x" + c + " (" + rTot + "x" + cTot + "):\n");
		sb.append(	"failed: " + nFailed + ", extra pivots: " + nExtraP + ", pivots: " + nPivots + 
					", contribution RxC: " + nContribR + "x" + nContribC + ", extra RxC: " + nExtraR + "x" + nExtraC + "\n");

		sb.append("        ");
		for (int i = 0; i < cTot; i++) sb.append(indexToString(idxCR[i]));
		
		if (data != null) {
			for (int i = 0; i < rTot; i++) {
				
				int iN = i * cTot;
				if (sb.length() > 0) sb.append("\n");		// don't add newline at start of matrix printout
				sb.append(indexToString(idxCR[cTot + i]) + "|");
				
				for (int j = 0; j < cTot; j++)
					if (Matrix.nearZero(data[iN + j])) 
							sb.append("   -   ");
					else	sb.append(Matrix.to5chars(data[iN + j], false));
				
				// any imaginary data comes as a second line under the real data line
//				if (idata != null) {
//					sb.append("     |\n       |");
//					for (int j = 0; j < cTot; j++)
//						if (Matrix.nearZero(data[iN + j]))
//								sb.append("       ");
//						else	sb.append(Matrix.to5chars(idata[iN + j], true));
//					
//					sb.append(" |\n       |");
//					for (int j = 0; j < cTot; j++) sb.append("       ");
//				}
				
				sb.append(" |");
			}
		} else sb.append("[null]");
		sb.append("\n");
		if (firstPivotFailed()) sb.append("(Failed on 1st pivot)\n");
		return sb.toString();	
	}

}
