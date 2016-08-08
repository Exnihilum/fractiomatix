package lkr74.matrixlib;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import javax.management.RuntimeErrorException;

public class MatrixMarketIO {
	
	BufferedReader br = null;
	StringBuilder sb;
	Matrix M = null;
	
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
	
	// constructor reads in file into stringbuffer in memory
	public MatrixMarketIO(String fName) {
		
		int state = 1, rows = 0, cols = 0, entries = 0;
		double[] iv = new double[2];
		boolean keepReading = true;
		String[] dataRow = null, mName = fName.split("[//.]+");
		
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
							M = new Matrix(mName[mName.length-2], rows, cols, Matrix.Type.Null_Complex);
							state = 8;
							break;
						// or only real numbers?
						} else {
							// allocate empty real matrix to fill up
							M = new Matrix(mName[mName.length-2], rows, cols, Matrix.Type.Null);
							state = 9;
							break;
						}
						
					// state 3 reads in all sparse coordinated values for a real matrix
					case 3:
						if ((s = br.readLine()) == null) { state = 5; break; }
						if (s.startsWith("%")) break;
						dataRow = s.split("[ ]+");
						if (dataRow == null) break;
					case 8:
						state = 3;
						M.valueTo(	Integer.valueOf(dataRow[0]) - 1,
									Integer.valueOf(dataRow[1]) - 1,
									Double.valueOf(dataRow[2]));
						if (--entries <= 0) state = 5;
						break;

					case 4:
						if ((s = br.readLine()) == null) { state = 5; break; }
						if (s.startsWith("%")) break;
						dataRow = s.split("[ ]+");
						if (dataRow == null) break;
					case 9:
						state = 4;
						iv[0] = Double.valueOf(dataRow[2]);
						iv[1] = Double.valueOf(dataRow[3]);
						M.valueToC(Integer.valueOf(dataRow[0]) - 1, Integer.valueOf(dataRow[1]) - 1, iv);
						if (--entries <= 0) state = 5;
						break;

					// state 5 = end of file case
					case 5:
						keepReading = false;
						break;
					default:
					}
				}
			} finally { br.close(); }

		} catch (IOException e) { e.printStackTrace(); }
	}
	
	public Matrix getMatrix() { return M; }
	
}
