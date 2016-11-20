package lkr74.matrixlib;

//the sparse node NspNode is referenced by the NSPMatrix both from a row/horizontal aspect, and
//from a column/vertical aspect, from the relevant referencing NspArrays in Hsp & Vsp
public class NspNode {
	public int r, c, offH, offV;
	public double v, iv;
	NspNode(int r, int c, double v, double iv) { this.r=r; this.c=c; this.v=v; this.iv=iv; }
	NspNode(int r, int c, double v, double iv, int offH, int offV) {
		this.r=r; this.c=c; this.v=v; this.iv=iv; this.offH = offH; this.offV = offV; }
	@Override
	protected NspNode clone() { NspNode n = new NspNode(r, c, v, iv); n.offH = offH; n.offV = offV; return n; }
}
