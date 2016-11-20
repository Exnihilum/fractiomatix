
package lkr74.matrixlib;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class MatrixMarketIO {
	
	BufferedReader br = null;
	StringBuilder sb;
	Matrix M = null;
	CSRMatrix M1 = null;
	NSPMatrix M2 = null;
	
//	private static boolean isSparse(String typeCode) { return typeCode.charAt(0) == 'C'; }
//	private static boolean isCoordinate(String typeCode) { return typeCode.charAt(0) == 'C'; }
//	private static boolean isDense(String typeCode) { return typeCode.charAt(0) == 'A'; }
//	private static boolean isArray(String typeCode) { return typeCode.charAt(0) == 'A'; }
//	private static boolean isComplex(String typeCode) { return typeCode.charAt(1) == 'C'; }
//	private static boolean isReal(String typeCode) { return typeCode.charAt(1) == 'R'; }
//	private static boolean isPattern(String typeCode) { return typeCode.charAt(1) == 'P'; }
//	private static boolean isInteger(String typeCode) { return typeCode.charAt(1) == 'I'; }
//	private static boolean isSymmetric(String typeCode) { return typeCode.charAt(2) == 'S'; }
//	private static boolean isGeneral(String typeCode) { return typeCode.charAt(2) == 'G'; }
//	private static boolean isSkew(String typeCode) { return typeCode.charAt(2) == 'K'; }
//	private static boolean isHermitian(String typeCode) { return typeCode.charAt(2) == 'H'; }
	
	// constructor reads a MatrixMarket file into a matrix type of choice:
	// toType = 0 -> read into Matrix
	// toType = 2 -> read into NSPMatrix
	public MatrixMarketIO(String fName, int toType, boolean threshold) {
		
		int state = 1, rows = 0, cols = 0, entries = 0, r = 0, c = 0, cOld = -1, i = -1;
		double[] iv = new double[2];
		boolean keepReading = true;
		String[] dataRow = null, mName = fName.split("[//.]+");
		NspArray aVsp = null;
		NspNode[] bVsp = null;							
		int offV = 0;							// offsVsp = offset into current Vsp buffer
		
		int debug_level = Matrix.DEBUG_LEVEL;
		Matrix.DEBUG_LEVEL = 0;
		
		try {
			br = new BufferedReader(new FileReader(fName));
			String s;
						
			try {
				// load in lines until EOF or a signal to stop
				while (keepReading) {
					
					switch(state) {
					
					// state 1 takes care of the header data
					case 1:
						if ((s = br.readLine()) == null) { state = 5; break; }
						// skip all commenting
						if (s.startsWith("%")) break;
						
						String[] header = s.split("[ ]+");
						// skip empty lines
						if (header == null) break;
						
						rows = Integer.valueOf(header[0]);
						cols = Integer.valueOf(header[1]);
						entries = Integer.valueOf(header[2]);
						if (rows < 1 || cols < 1 || entries < 1) throw new RuntimeException("MatrixMarketIO(): Invalid header.");
						state = 2;
						//break;
						
					// state 2 allocates correct matrix type
					case 2:
						if ((s = br.readLine()) == null) { state = 5; break; }
						if (s.startsWith("%")) break;

						dataRow = s.split("[ ]+");
						// skip empty lines
						if (dataRow == null) break;
						
						if (dataRow.length < 3)
							throw new RuntimeException("MatrixMarketIO(): Invalid data row.");

						// do we have complex numbers?
						if (dataRow.length == 4) {
							// allocate empty complex matrix to fill up
							switch (toType) {
								case 0: M = new Matrix(mName[mName.length-2], rows, cols, Matrix.Type.Null_Complex, 1);
										state = 9; break;
								case 2: M2 = new NSPMatrix(mName[mName.length-2], rows, cols, Matrix.Type.Null_Complex, 1);
										aVsp = M2.Vsp[0];
										bVsp = aVsp.array;							
										state = 29; break;
							}
						// or only real numbers?
						} else {
							// allocate empty real matrix to fill up
							switch (toType) {
								case 0: M = new Matrix(mName[mName.length-2], rows, cols, Matrix.Type.Null, 1);
										state = 8; break;
								case 2: M2 = new NSPMatrix(mName[mName.length-2], rows, cols, Matrix.Type.Null, 1);
										aVsp = M2.Vsp[0];
										bVsp = aVsp.array;							
										state = 28; break;
							}
						}
						break;
						
					// The states for reading into a Matrix type follow
					// state 3 reads in all sparse coordinated values for a real matrix
					case 3:
						if ((s = br.readLine()) == null) { state = 5; break; }
						if (s.startsWith("%")) break;
						dataRow = s.split("[ ]+");
						if (dataRow == null) break;
					case 8:
						state = 3;
						double v = Double.valueOf(dataRow[2]);
						// only read values above the configured threshold for a non-zero value
						if (!threshold || (threshold && !Matrix.nearZero(v))) {
							M.valueTo(Integer.valueOf(dataRow[0]) - 1, Integer.valueOf(dataRow[1]) - 1, v);
							M.nNZ++;
						}
						if (--entries <= 0) state = 5;					// we know the number of nonzeroes we expect and stop after they're read
						break;
						
					// complex matrix loop
					case 4:
						if ((s = br.readLine()) == null) { state = 5; break; }
						if (s.startsWith("%")) break;
						dataRow = s.split("[ ]+");
						if (dataRow == null) break;
					case 9: {
						state = 4;
						iv[0] = Double.valueOf(dataRow[2]);
						iv[1] = Double.valueOf(dataRow[3]);
						// only read values above the configured threshold for a non-zero value
						if (!threshold || (threshold && (!Matrix.nearZero(iv[0]) || !Matrix.nearZero(iv[1])))) {
							M.valueToC(Integer.valueOf(dataRow[0]) - 1, Integer.valueOf(dataRow[1]) - 1, iv);
							M.nNZ++;
						}
						if (--entries <= 0) state = 5;					// we know the number of nonzeroes we expect and stop after they're read
						break;
					}
					
					// The states for reading into a NSPMatrix type follow
					// real matrix loop
					case 23:
						if ((s = br.readLine()) == null) { state = 5; break; }
						if (s.startsWith("%")) break;
						dataRow = s.split("[ ]+");
						if (dataRow == null) break;
					case 28:
						state = 23;
						r = Integer.valueOf(dataRow[0]) - 1;
						c = Integer.valueOf(dataRow[1]) - 1;
						double v1 = Double.valueOf(dataRow[2]), v2 = 0;
						if (dataRow.length == 4) v2 = Double.valueOf(dataRow[3]);

						// update to next Vsp column array if we're on the next column in the file
						if (c > cOld) {
							if (aVsp.nodes == 0) M2.Vsp[i < 0 ? 0 : i].array = null;	// check previous Vsp array, dereference if empty
							cOld = c; i++;
							aVsp = M2.Vsp[i];
							bVsp = aVsp.array;
							offV = 0;
						}
						
						// only read values above the configured threshold for a non-zero value
						if (!threshold || (threshold && (!Matrix.nearZero(v1) || !Matrix.nearZero(v2)))) {
							
							// if current sparse Vsp array is full, add another CSR2_ALLOCBLOCK elements
							if (NSPMatrix.updateArraySize(aVsp.nodes + 1, aVsp))
								bVsp = aVsp.array;
	
							bVsp[offV] = new NspNode(r, c, v1, v2, 0, offV);	// insert row index w. vertical array offset
							M2.nNZ++;											// total node count incremented
							if (r == c) M2.pivotNsp[r] = bVsp[offV];			// if it's a pivot, put in fast-access array
							M2.readjustHalfBandwidth(r, c);
							offV++;
						}

						if (--entries <= 0) state = 25;					// we know the number of nonzeroes we expect and stop after they're read
						break;
						
					// state 5 = end of file case
					case 5:
						keepReading = false;
						break;
						
					// state 25 = end of file for NSPMatrix, fix the horisontal references
					case 25:
						M2.crosslinkHsp();
						keepReading = false;
						break;
						
					default:
					}
				}
			} finally { br.close(); }

		} catch (IOException e) { e.printStackTrace(); }
		
		Matrix.DEBUG_LEVEL = debug_level;

	}
	
	public Matrix getMatrix() { return M; }
	public NSPMatrix getNSPMatrix() { return M2; }

}