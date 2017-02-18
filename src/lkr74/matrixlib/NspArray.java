package lkr74.matrixlib;

//a NspArray is a freestanding sparse datastructure that can be manipulated by relevant methods,
//multiplied, recombined, added and so on, if the method returns an assembled NspArray, then it
//has to be integrated with the arrays of the opposite aspect: Hsp -> update Vsp, Vsp -> update Hsp
public class NspArray {
	public int nodes, size;
	public NspNode[] array;
	
	public NspArray(int nodes, int size, NspNode[] array) { this.nodes = nodes; this.size = size; this.array = array; }

	@Override
	protected NspArray clone() { return new NspArray(nodes, size, array != null ? array.clone() : null); }
	
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
