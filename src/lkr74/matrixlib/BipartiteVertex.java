package lkr74.matrixlib;

public class BipartiteVertex {
	
	// Definition of a graph vertex in a Bipartite decomposition
	// Leonard Krylov 2016

	protected static final int BPVERT_ALLOCBLOCK = 16;
	public static final int C_VERTEX = 1<<30, SUPERSINK = 1<<29, FLAGS = C_VERTEX + SUPERSINK;
	public int i;						// i = index of this vertex in bipartite graph (on either partition side)
	int distance, idx;					// tells alternating path distance in Hopcroft-Karp algorithm, idx = type and true index of vertex
	public int edges=0;					// edges = edge count
	int cutedges=0;						// cutedges are the number of previously cut edges, they are appended from end of array
	int checked=0;						// checked is an index telling how many edges have been checked by a DAG-walking algorithm
	int[] edge = null;					// the array of edges
	
	// skeleton vertex constructors
	BipartiteVertex(int i, int idx) { this.i = i; edge = new int[BPVERT_ALLOCBLOCK]; this.idx = idx; }
	BipartiteVertex(int i, int[] edge, int idx) { this.i = i; this.edge = edge; this.idx = idx; }

	@Override
	public BipartiteVertex clone() {
		// first clone the static data of this vertex
		Object O = null;
		try { O = super.clone(); } catch (CloneNotSupportedException e) { e.printStackTrace(); }
		BipartiteVertex vtx = (BipartiteVertex) O;
		vtx.edge = vtx.edge == null ? null : edge.clone();
		return vtx;
	}

	// method fills in the holes created by cut edges, progressing from right
	// optimally, this method is called when more than half of the left part of array constitutes holes
	public void purgeFlagged() {
		for (int e = 0, e2 = 0; e < edges + cutedges; e++) if (edge[e] != -1) edge[e2++] = edge[e];
		cutedges = 0;
	}
	
	
	// method test if edge array needs expanding and is called with the requested number of additional elements
	boolean updateEdgeArray(int i) {
		
		int totaledges = edges + i, alength = edge.length;
		if (totaledges >= alength) {
			int[] edge2 = new int[trimToAllocBlock(totaledges)];
			
			if (cutedges > 0) {							// copy over existing edges, skip over cut edge "holes" assigned to -1
				for (int e = 0, e2 = 0, eEnd = edges + cutedges; e < eEnd; e++)
					if (edge[e] != -1) edge2[e2++] = edge[e];
			} else {									// straight copy-over of existing edges
				for (int e = 0; e < edges; e++) edge2[e] = edge[e];
			}
			edge = edge2;
			return true;
		}		
		return false;					// no resize needed, return false
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			HELPER METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	public int edge(int e) { return edge[e] & (0xFFFFFFFF - FLAGS); }
	public boolean is_C_vertex(int e) { return ((edge[e] & FLAGS) == C_VERTEX); }
	public boolean is_supersink(int e) { return ((edge[e] & FLAGS) == SUPERSINK); }
	public void make_C_edge(int e) { edge[e] |= C_VERTEX; }
	public void supersink_edge(int e) { edge[e] |= SUPERSINK; }
	public void make_R_edge(int e) { edge[e] &= (0xFFFFFFFF - C_VERTEX); }

	private static int trimToAllocBlock(int v) { return (v & (0xFFFFFFFF - BPVERT_ALLOCBLOCK + 1)) + BPVERT_ALLOCBLOCK; }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
				
		String dst = distance == Integer.MAX_VALUE ? dst = "inf" : Integer.toString(distance);
		String idx = ((this.idx & SUPERSINK) != 0 ? "ssink" : Integer.toString(this.idx & (0xFFFFFFFF - FLAGS)) + "(" + i + ")");
		
		sb.append(((this.idx & C_VERTEX) != 0 ? "C" : "R") + idx + ((this.idx & C_VERTEX) == 0 ? " (d:" + dst + ") :" : " :") + "[");
		for (int e = 0, e2 = 0; e < edges + cutedges; e++) {
			if (edge[e] != -1) {
				if (is_supersink(e)) sb.append("ssink");
				else sb.append(edge(e));
				if (++e2 < edges) sb.append(", ");
			}
		}
		sb.append("]");
		
		if (cutedges > 0) {
			sb.append(" cut: [");
			for (int ce = edge.length - 1; ce >= edge.length - cutedges; ce--) {
				sb.append(edge[ce]);
				if (ce > edge.length - cutedges) sb.append(", "); }
			sb.append("]");
		}
		sb.append("\n");
		return sb.toString();
	}

}
