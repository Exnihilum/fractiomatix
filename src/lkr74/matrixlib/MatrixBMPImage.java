package lkr74.matrixlib;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

import io.nayuki.bmpio.AbstractRgb888Image;
import io.nayuki.bmpio.BmpImage;
import io.nayuki.bmpio.BmpWriter;

final class MatrixBMPImage extends AbstractRgb888Image {

	// Class uses https://github.com/nayuki/BMP-IO image file in/out class to output a BMP image of a flat array Matrix
	double vMax = Matrix.ROUNDOFF_ERROR, vMedian = 0;
	BmpImage bmp = null;
	Matrix A;

	public MatrixBMPImage(Matrix A) {
		super(A.isComplex() ? A.N * 2 : A.N, A.M * 2);
		bmp = new BmpImage();
		bmp.image = this;
		this.A = A;
		// find largest value of matrix to divide all others by
		for(double v: A.data) {
			if (v < 0) v = -v;
			if (vMax < v) vMax = v;
			vMedian += v;
		}
		vMedian /= (double)(A.data.length);
	}

	public int getRgb888Pixel(int c, int r) {

		int r1 = (r >= A.M ? r - A.M : r);
		double v1 = (c >= A.N ? A.valueOfC(r1, c-A.N)[1] : A.valueOf(r1, c));
		double v = v1 / vMedian;
		int pixel = 0;
		
		if (r < A.M) {
			// the top half of image is color coded matrix values
			// positive values as blues, negative as reds
			pixel = v > 0 ? ((int)(v * 255))&0xff : (((int)(v * 255))&0xff)<<16;
			// make relatively  microscopic but nonzero values visible
			if (!Matrix.nearZero(v) && pixel == 0) pixel = 0x121212;
		} else {
			// the bottom half of image is white if zero value, black otherwise
			if (Matrix.nearZero(v1)) pixel = 0xffffff;
		}
		return pixel;
	}
	
	public void write() {
		
		File file = new File(A.name + ".bmp");
		FileOutputStream fos = null;
		try {	fos = new FileOutputStream(file);
		} catch (FileNotFoundException e) { e.printStackTrace(); }

		try {	BmpWriter.write(fos, bmp);
				fos.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void write(String fName) {
		
		File file = new File(fName.replace(".", "_") + ".bmp");
		FileOutputStream fos = null;
		try {	fos = new FileOutputStream(file);
		} catch (FileNotFoundException e) { e.printStackTrace(); }

		try {	BmpWriter.write(fos, bmp);
				fos.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


}
