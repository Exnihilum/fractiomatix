package lkr74.matrixlib;

public class NSPArray {
	
	// A NspArray is a freestanding sparse datastructure that can be manipulated by relevant methods,
	// multiplied, recombined, added and so on, if the method returns an assembled NspArray, then it
	// has to be integrated with the arrays of the opposite aspect: Hsp -> update Vsp, Vsp -> update Hsp
	// This is a helper class for the NSPMatrix class, since a NSPMAtrix consists of NSPArrays
	// Leonard Krylov 2016
	
	public int nodes, size;
	public NSPNode[] array;
	
	public NSPArray(int nodes, int size, NSPNode[] array) { this.nodes = nodes; this.size = size; this.array = array; }

	@Override
	protected NSPArray clone() { return new NSPArray(nodes, size, array != null ? array.clone() : null); }
	
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
		sb.append("Nodes: " + nodes + ", size: " + size + ", data:\n[");
		if (array != null && array[0] != null) {
			for (int i = 0; i < nodes; i++) sb.append(Matrix.to5chars(array[i].v, false) + (i == nodes-1 ? "]\n[" : " "));
			for (int i = 0; i < nodes; i++) sb.append(Matrix.to5chars(array[i].iv, true) + (i == nodes-1 ? "]\n[" : " "));
			for (int i = 0; i < nodes; i++) sb.append(indexToString(array[i].r) + (i == nodes-1 ? "]\n[" : " "));
			for (int i = 0; i < nodes; i++) sb.append(indexToString(array[i].c) + (i == nodes-1 ? "]\n" : " "));
		} else sb.append("null]\n");
		return sb.toString();
	}
}
