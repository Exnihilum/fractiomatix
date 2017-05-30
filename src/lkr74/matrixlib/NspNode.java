package lkr74.matrixlib;

public class NSPNode {
	
	// The sparse node NSPNode is the building block of a NSPMAtrix, constituting every NONZERO entry in the matrix
	// it is referenced by the NSPMatrix both from a row/horizontal aspect, and
	// from a column/vertical aspect, from the relevant referencing NspArrays in Hsp & Vsp
	// Leonard Krylov 2016
	
	int r, c, offH, offV;
	double v, iv;
	public NSPNode(int r, int c, double v) { this.r=r; this.c=c; this.v=v; }
	public NSPNode(int r, int c, double v, double iv) { this.r=r; this.c=c; this.v=v; this.iv=iv; }
	public NSPNode(int r, int c, double v, double iv, int offH, int offV) {
		this.r=r; this.c=c; this.v=v; this.iv=iv; this.offH = offH; this.offV = offV; }
	@Override
	protected NSPNode clone() { NSPNode n = new NSPNode(r, c, v, iv); n.offH = offH; n.offV = offV; return n; }
	
	public int c() { return c; }
	public int r() { return r; }
	public void add(double v) { this.v += v; }
	public void addC(double v, double iv) { this.v += v; this.iv += iv; }
}

