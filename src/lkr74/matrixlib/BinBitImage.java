package lkr74.matrixlib;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.Arrays;

import jogamp.opengl.glu.nurbs.Bin;

public class BinBitImage {

	static final double ROUNDOFF_ERROR = Matrix.ROUNDOFF_ERROR;

	private Matrix M;			// each BinBitImage instance is referencing a specific matrix, supplied in the constructor
	
	protected int bitSets;
	protected long[] data;
	
	// variables for the static fast-access binBit buffers
	private static int binBitLength;
	private static DoubleBuffer binBitD;
	private static IntBuffer binBitI;

	static int DEBUG_LEVEL = 2;

	// everything that needs to be statically set up at load time for the Matrix class is done here
	{
		binBitD = ByteBuffer.allocateDirect(128 * 8).asDoubleBuffer();
		binBitI = ByteBuffer.allocateDirect(128 * 4).asIntBuffer();
		// for fast access, binBit buffers should be direct and have no backing arrays
//        System.out.println("Is a direct buffer: " + binBitD.isDirect());
//        System.out.println("Buffer has a backing array: " + binBitD.hasArray());
        binBitLength = 0;
	}

	
	// skeleton BinBitImage constructor method
	public BinBitImage() { M = null; }
	
	// instantiates a BinBitImage and generates a bitwise nonzero-value identifying image of a matrix
	public BinBitImage(Matrix M) {this.M = M; make(); }

	
	// generates a bitwise nonzero-value identifying image of a matrix
	// call this method when the dimensions of bitImage haven't changed and data can be directly overwritten
	public void make() {
		
		// get data in Matrix row-column format from matrix subclass's native format (unless it's an empty matrix
		double[][] dataSet = M.getData();		
		double[] dataM = dataSet[0], idataM = dataSet[1];
		
		// 8x8 matrix fits perfectly into a long
		long bitset = 0;
		if (M.M < 9 && M.N < 9) {
			data = new long[1];
			if (M.isNull()) return;		// matrix is null, no need to do anything else but allocation
			
			// sets bits linearly along the bitset, first row-col in lowest/rightmost position
			for (int i = 0; i < M.M; i++) {
				int irow = i * M.N;
				for (int j=0; j < M.N; j++) {
					if (dataM != null && !Matrix.nearZero(dataM[irow + j]))
						bitset |= (0x1L<<(irow + j));
					else
					if (idataM != null && !Matrix.nearZero(idataM[irow + j]))
						bitset |= (0x1L<<(irow + j));					
				}
			}
			data[0] = bitset;
			bitSets = 0;
			
			if (debugTrigger(3)) System.out.println(this.toString());
			return;
		}
		
		int sets = M.N >> 6, rest = M.N % 64;
		bitSets = sets + (rest == 0 ? 0 : 1);
		data = new long[bitSets * M.M];
		if (M.isNull()) return;		// matrix is null, no need to do anything else but allocation
		
		// for every row of matrix
		for (int i = 0; i < M.M; i++) {
			int iN = i * M.N;
			// go through this row processing 64 elements at a time, find all non-zeroes and set their bit
			for (int j = 0; j < (sets<<6); j += 64) {
				bitset = 0;
				if (dataM != null) {
					for (int k = 0, ij1 = iN + j; k < 64; k++, ij1++)
						if (!Matrix.nearZero(dataM[ij1])) bitset |= (0x1L<<k);
				}
				if (idataM != null) {
					for (int k = 0, ij1 = iN + j; k < 64; k++, ij1++)
						if (!Matrix.nearZero(idataM[ij1])) bitset |= (0x1L<<k);
				}
				data[i * bitSets + (j>>6)] = bitset;
			}
			bitset = 0;
			int j = sets << 6;
			// this takes care of a remainder of this row's elements that didn't fit into previous 64-element parcels
			if (dataM != null) 
				for (int k = 0, ij1 = iN + j; k < rest; k++, ij1++)
					if (!Matrix.nearZero(dataM[ij1])) bitset |= (0x1<<k);
			if (idataM != null) 
				for (int k = 0, ij1 = iN + j; k < rest; k++, ij1++)
					if (!Matrix.nearZero(idataM[ij1])) bitset |= (0x1<<k);
			
			// also need checking if there WAS any remainder, to avoid out of bounds exception!
			if (rest > 0) data[i * bitSets + sets] = bitset;
		}
		
		if (debugTrigger(3)) System.out.println(this.toString());
	}
	
	public BinBitImage clone(Matrix M) {
		BinBitImage bitimage = new BinBitImage();
		bitimage.M = M;
		bitimage.data = data.clone();
		bitimage.bitSets = bitSets;
		
		if (debugTrigger(2)) System.out.println(this.toString());
		return bitimage;
	}
	
	
	public void zero() { for (int i = 0; i < data.length; i++) data[i] = 0; }
	public void set() { for (int i = 0; i < data.length; i++) data[i] = 0xffffffffffffffffL; }
	
	// OR-s this another bitImage into this bitImage
	void or(BinBitImage bI) {
		int l = data.length > bI.data.length ? bI.data.length : data.length;
		for (int i = 0; i < l; i++) this.data[i] |= bI.data[i];
		
		if (debugTrigger(2)) System.out.println(this.toString());
	}

	// AND-s this another bitImage into this bitImage
	void and(BinBitImage bI) {
		int l = data.length > bI.data.length ? bI.data.length : data.length;
		for (int i = 0; i < l; i++) this.data[i] &= bI.data[i];
		
		if (debugTrigger(2)) System.out.println(this.toString());
	}

	
	// set a bit in bitImage
	void setBit(int r, int c) {
		// 8x8 single-long matrix case
		if (bitSets == 0)	data[0] |= 0x1L<<(r * M.N + c);
		// large matrix basecase
		else { int a = c%64; data[r * bitSets + c/64] |= (0x1L << a); }
		
		if (debugTrigger(3)) System.out.println(this.toString());
	}

	// clear a bit in bitImage
	void clearBit(int r, int c) {
		// 8x8 single-long matrix case
		if (bitSets == 0)	data[0] &= 0xFFFFFFFFFFFFFFFFL ^ (0x1L<<(r * M.N + c));
		// large matrix basecase
		else { int a = c%64; data[r * bitSets + c/64] &= 0xFFFFFFFFFFFFFFFFL ^ (0x1L<<a); }
		
		if (debugTrigger(3)) System.out.println(this.toString());
	}
	
	// transposes two bits across a square matrix
	void transposeBit(int r, int c) {
		// 8x8 single-long matrix case
		if (bitSets == 0) {
			long rcmask = 0x1L<<(r * M.N + c), crmask = 0x1L<<(c * M.N + r);
			if ((data[0] & rcmask) != 0) {
				if ((data[0] & crmask) == 0)
						// swap bits only if first is 1 and second is 0
					data[0] = (data[0] - rcmask) | crmask;
			} else if ((data[0] & crmask) != 0)
				// we know first bit is 0, since second is 1, swapping is obligatory
				data[0] = (data[0] - crmask) | rcmask;
			return;
		}
		
		// large matrix basecase
		long bitcr = (0x1L << (c%64)), bitrc = (0x1L << (r%64));
		int rofs = r * bitSets + c/64, cofs = c * bitSets + r/64;
		if ((data[rofs] & bitcr) != 0) {
			if ((data[cofs] & bitrc) == 0)
				data[rofs] -= bitcr; data[cofs] |= bitrc;
		} else if ((data[cofs] & bitrc) != 0) {
			data[rofs] |= bitcr; data[cofs] -= bitrc;
		}
		
		if (debugTrigger(3)) System.out.println(this.toString());
	}
	
	void swapRows(int r1, int r2) {
		// swap bitImage rows, special case for 8x8 matrices and less
		if (bitSets == 0) {
			int r1N = r1 * M.N, r2N = r2 * M.N;
			long row1mask = (0xFFL >>> (8 - M.N)) << r1N;
			long row2mask = (0xFFL >>> (8 - M.N)) << r2N;
			long row1bits = data[0] & row1mask;
			long row2bits = data[0] & row2mask;
			data[0] = (data[0] ^ (row1bits | row2bits)) | ((row1bits>>r1N)<<r2N) | ((row2bits>>r2N)<<r1N);
			return;
		}
		// base case for larger matrices
		for (int i = 0; i < bitSets; i++) {
			long bitset = data[r1*bitSets + i];
			data[r1*bitSets + i] = data[r2*bitSets + i];
			data[r2*bitSets + i] = bitset;
		}
		
		if (debugTrigger(3)) System.out.println(this.toString());
	}
	
	protected boolean equals(BinBitImage bI) {
		long[] data1 = this.data;
		long[] data2 = bI.data;
		int l = data1.length < data2.length ? data1.length : data2.length;
		for (int i = 0; i < l; i++)
			if (data1[i] != data2[i]) return false;
		return true;
	}

	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			BINARY TREE BITSEARCH METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	// TODO: use ROUNDOFF_ERROR in here
	// does fast comparison of only those value pairs that are signed as nonzero in bitvector v
	public static boolean compare(double[] c, double[] d, long v) {	
		for(int i = 0; i < 2; i++, v >>= 32) {
			if ((v & 0xffffffffL)!=0) { if ((v & 0xffff0000L)!=0) { if ((v & 0xff000000L)!=0) {
						if ((v & 0xf0000000L)!=0) { if ((v & 0xC0000000L)!=0) {
								if ((v & 0x80000000L)!=0) 	if (c[31] != d[31]) return false;
								if ((v & 0x40000000L)!=0)	if (c[30] != d[30]) return false; }		
							if ((v & 0x30000000L)!=0) {
								if ((v & 0x20000000L)!=0)	if (c[29] != d[29]) return false;
								if ((v & 0x10000000L)!=0)	if (c[28] != d[28]) return false; }								
						} if ((v & 0x0f000000L)!=0) { if ((v & 0x0C000000L)!=0) {
								if ((v & 0x08000000L)!=0)	if (c[27] != d[27]) return false;
								if ((v & 0x04000000L)!=0)	if (c[26] != d[26]) return false; }		
							if ((v & 0x03000000L)!=0) {
								if ((v & 0x02000000L)!=0)	if (c[25] != d[25]) return false;
								if ((v & 0x01000000L)!=0)	if (c[24] != d[24]) return false; }}}		
					if ((v & 0x00ff0000L)!=0) { if ((v & 0x00f00000L)!=0) { if ((v & 0x00C00000L)!=0) {
								if ((v & 0x00800000L)!=0)	if (c[23] != d[23]) return false;
								if ((v & 0x00400000L)!=0)	if (c[22] != d[22]) return false; }								
							if ((v & 0x00300000L)!=0) {
								if ((v & 0x00200000L)!=0)	if (c[21] != d[21]) return false;
								if ((v & 0x00100000L)!=0)	if (c[20] != d[20]) return false; }}		
						if ((v & 0x000f0000L)!=0) { if ((v & 0x000C0000L)!=0) {
								if ((v & 0x00080000L)!=0)	if (c[19] != d[19]) return false;
								if ((v & 0x00040000L)!=0)	if (c[18] != d[18]) return false; }			
							if ((v & 0x00030000L)!=0) {
								if ((v & 0x00020000L)!=0)	if (c[17] != d[17]) return false;
								if ((v & 0x00010000L)!=0)	if (c[16] != d[16]) return false; }}}}		
				if ((v & 0x0000ffffL)!=0) { if ((v & 0x0000ff00L)!=0) { if ((v & 0x0000f000L)!=0) { if ((v & 0x0000C000L)!=0) {
								if ((v & 0x00008000L)!=0) 	if (c[15] != d[15]) return false;
								if ((v & 0x00004000L)!=0) 	if (c[14] != d[14]) return false; }								
							if ((v & 0x00003000L)!=0) {
								if ((v & 0x00002000L)!=0) 	if (c[13] != d[13]) return false;
								if ((v & 0x00001000L)!=0) 	if (c[12] != d[12]) return false; }}		
						if ((v & 0x00000f00L)!=0) { if ((v & 0x00000C00L)!=0) {
								if ((v & 0x00000800L)!=0)	if (c[11] != d[11]) return false;
								if ((v & 0x00000400L)!=0)	if (c[10] != d[10]) return false; }							
							if ((v & 0x00000300L)!=0) {
								if ((v & 0x00000200L)!=0)	if (c[9] != d[9]) return false;
								if ((v & 0x00000100L)!=0)	if (c[8] != d[8]) return false; }}}		
					if ((v & 0x000000ffL)!=0) { if ((v & 0x000000f0L)!=0) { if ((v & 0x000000C0L)!=0) {
								if ((v & 0x00000080L)!=0)	if (c[7] != d[7]) return false;
								if ((v & 0x00000040L)!=0)	if (c[6] != d[6]) return false; }								
							if ((v & 0x00000030L)!=0) {
								if ((v & 0x00000020L)!=0)	if (c[5] != d[5]) return false;
								if ((v & 0x00000010L)!=0)	if (c[4] != d[4]) return false;	}}		
						if ((v & 0x0000000fL)!=0) { if ((v & 0x0000000CL)!=0) {
								if ((v & 0x00000008L)!=0)	if (c[3] != d[3]) return false;
								if ((v & 0x00000004L)!=0)	if (c[2] != d[2]) return false; }							
							if ((v & 0x00000003L)!=0) {
								if ((v & 0x00000002L)!=0)	if (c[1] != d[1]) return false;
								if ((v & 0x00000001L)!=0)	if (c[0] != d[0]) return false;						
					}}}}}
		}
		return true;
	}
	

	// writes out value-index pair of 64-long input doubles array by checking for nonzero values in bitvector v
	// the write goes to 2 preallocated DirectBuffers of Matrix class. The doubles are read from offset "id"
	public static int compact(double[] d, long v, int id) {
		binBitD.position(0);
		binBitI.position(0);
		for(; v != 0; v = (v>>32) & 0xffffffffL, id += 32) {
			if ((v & 0xffffffffL)!=0) { if ((v & 0x0000ffffL)!=0) {
					if ((v & 0x000000ffL)!=0) { if ((v & 0x0000000fL)!=0) {
							if ((v & 0x00000003L)!=0) {
								if ((v & 0x00000001L)!=0) 	{binBitI.put(id); binBitD.put(d[id]);}
								if ((v & 0x00000002L)!=0)	{binBitI.put(id+1); binBitD.put(d[id+1]);} }		
							if ((v & 0x0000000cL)!=0) {
								if ((v & 0x00000004L)!=0)	{binBitI.put(id+2); binBitD.put(d[id+2]);}
								if ((v & 0x00000008L)!=0)	{binBitI.put(id+3); binBitD.put(d[id+3]);} }}		
						if ((v & 0x000000f0L)!=0) { if ((v & 0x00000030L)!=0) {
								if ((v & 0x00000010L)!=0)	{binBitI.put(id+4); binBitD.put(d[id+4]);}
								if ((v & 0x00000020L)!=0)	{binBitI.put(id+5); binBitD.put(d[id+5]);} }		
							if ((v & 0x000000c0L)!=0) {
								if ((v & 0x00000040L)!=0)	{binBitI.put(id+6); binBitD.put(d[id+6]);}
								if ((v & 0x00000080L)!=0)	{binBitI.put(id+7); binBitD.put(d[id+7]);} }}}		
					if ((v & 0x0000ff00L)!=0) { if ((v & 0x00000f00L)!=0) {
							if ((v & 0x00000300L)!=0) {
								if ((v & 0x00000100L)!=0)	{binBitI.put(id+8); binBitD.put(d[id+8]);}
								if ((v & 0x00000200L)!=0)	{binBitI.put(id+9); binBitD.put(d[id+9]);} }								
							if ((v & 0x00000c00L)!=0) {
								if ((v & 0x00000400L)!=0)	{binBitI.put(id+10); binBitD.put(d[id+10]);}
								if ((v & 0x00000800L)!=0)	{binBitI.put(id+11); binBitD.put(d[id+11]);} }}		
						if ((v & 0x0000f000L)!=0) { if ((v & 0x00003000L)!=0) {
								if ((v & 0x00001000L)!=0)	{binBitI.put(id+12); binBitD.put(d[id+12]);}
								if ((v & 0x00002000L)!=0)	{binBitI.put(id+13); binBitD.put(d[id+13]);} }			
							if ((v & 0x0000c000L)!=0) {
								if ((v & 0x00004000L)!=0)	{binBitI.put(id+14); binBitD.put(d[id+14]);}
								if ((v & 0x00008000L)!=0)	{binBitI.put(id+15); binBitD.put(d[id+15]);} }}}}		
				if ((v & 0xffff0000L)!=0) { if ((v & 0x00ff0000L)!=0) {
						if ((v & 0x000f0000L)!=0) { if ((v & 0x00030000L)!=0) {
								if ((v & 0x00010000L)!=0) 	{binBitI.put(id+16); binBitD.put(d[id+16]);}
								if ((v & 0x00020000L)!=0) 	{binBitI.put(id+17); binBitD.put(d[id+17]);} }								
							if ((v & 0x000c0000L)!=0) {
								if ((v & 0x00040000L)!=0) 	{binBitI.put(id+18); binBitD.put(d[id+18]);}
								if ((v & 0x00080000L)!=0) 	{binBitI.put(id+19); binBitD.put(d[id+19]);} }}		
						if ((v & 0x00f00000L)!=0) {
							if ((v & 0x00300000L)!=0) {
								if ((v & 0x00100000L)!=0)	{binBitI.put(id+20); binBitD.put(d[id+20]);}
								if ((v & 0x00200000L)!=0)	{binBitI.put(id+21); binBitD.put(d[id+21]);} }							
							if ((v & 0x00c00000L)!=0) {
								if ((v & 0x00400000L)!=0)	{binBitI.put(id+22); binBitD.put(d[id+22]);}
								if ((v & 0x00800000L)!=0)	{binBitI.put(id+23); binBitD.put(d[id+23]);} }}}		
					if ((v & 0xff000000L)!=0) { if ((v & 0x0f000000L)!=0) {
							if ((v & 0x03000000L)!=0) {
								if ((v & 0x01000000L)!=0)	{binBitI.put(id+24); binBitD.put(d[id+24]);}
								if ((v & 0x02000000L)!=0)	{binBitI.put(id+25); binBitD.put(d[id+25]);} }								
							if ((v & 0x0c000000L)!=0) {
								if ((v & 0x04000000L)!=0)	{binBitI.put(id+26); binBitD.put(d[id+26]);}
								if ((v & 0x08000000L)!=0)	{binBitI.put(id+27); binBitD.put(d[id+27]);} }}		
						if ((v & 0xf0000000L)!=0) { if ((v & 0x30000000L)!=0) {
								if ((v & 0x10000000L)!=0)	{binBitI.put(id+28); binBitD.put(d[id+28]);}
								if ((v & 0x20000000L)!=0)	{binBitI.put(id+29); binBitD.put(d[id+29]);} }							
							if ((v & 0xc0000000L)!=0) {
								if ((v & 0x40000000L)!=0)	{binBitI.put(id+30); binBitD.put(d[id+30]);}
								if ((v & 0x80000000L)!=0)	{binBitI.put(id+31); binBitD.put(d[id+31]);}							
				}}}}}}
		binBitLength = binBitI.position();
		binBitD.position(0);
		binBitI.position(0);
		return binBitLength;
	}
	


	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			HELPER/INLINE MTHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// debugTrigger() limits local debug level by debug level of Matrix class 
	private static boolean debugTrigger(int trigger) { return Matrix.DEBUG_LEVEL - DEBUG_LEVEL >= 0 && DEBUG_LEVEL >= trigger; }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	private static int MAX_PRINTEXTENT = Matrix.MAX_PRINTEXTENT;

	public String toString() {
		
		StringBuilder sb = new StringBuilder();
		System.out.println("matrix " + this.M.name + "'s bitImage:");
		if (bitSets == 0) {
			for (int i = 0; i < 64 && i < M.M * M.N; i++) {
				sb.append(((data[0] & (0x1L<<i)) > 0) ? "1 " : ". ");
				if ((i + 1) % M.N == 0) sb.append("\n");
			}
			return sb.toString();
		}
		
		int maxM = M.M > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : M.M;
		int maxN = M.N > MAX_PRINTEXTENT ? MAX_PRINTEXTENT : M.N;
		int maxBitSets = (maxN>>6) + (maxN%64 == 0 ? 0 : 1);

		for (int i = 0; i < maxM; i++) {
			for (int j = 0; j < maxBitSets; j++) {
				for (int k = 0; k < maxN; k++)
					if (M.N > j * 64 + k)
						sb.append(((data[i * bitSets + j]&(0x1L<<k))!=0)?"1 ":". ");
			}
			sb.append("\n");
		}
		sb.append("\n");
		return sb.toString();
	}


	// returns String version of current binBitI & binBitD contents
	protected static String binBitToString() {
		if (!debugTrigger(DEBUG_LEVEL)) return "";
		double[] ad = new double[binBitLength];
		int[] ai = new int[binBitLength];
		binBitI.position(0);
		binBitD.position(0);
		binBitI.get(ai, 0, binBitLength);
		binBitD.get(ad, 0, binBitLength);
		return "binBitI: " + binBitLength + " " + Arrays.toString(ai) + "\binBitI length: " + binBitLength + " " + Arrays.toString(ad) + "\n";
	}

	
}
