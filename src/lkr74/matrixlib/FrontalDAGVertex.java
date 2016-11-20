package lkr74.matrixlib;

public class FrontalDAGVertex implements Cloneable {
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			DIRECTED ACYCLIC GRAPH VERTEX DATASTRUCTURE FOR MULTIFRONTAL LINEAR SOLVER
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	protected static final int FDAGV_ALLOCBLOCK = 16;
	public static final int L_CHILD = 1<<29, U_CHILD = 1<<30, LU_CHILD = L_CHILD | U_CHILD, SUPERNODE = 1<<28, N_CHILD = 0;
	public static final int L_PATH = L_CHILD, U_PATH = U_CHILD, FLAGS = (L_CHILD|U_CHILD|SUPERNODE);
	public int i;						// i = index of this vertex in DAG (=pivot)
	public int edges=0;					// edges = edge count
	int cutedges=0;						// cutedges are the number of previously cut edges, they are appended from end of array
	int checked=0;						// checked is an index telling how many edges have been checked by a DAG-walking algorithm
	int parentLU = Integer.MAX_VALUE;	// fast-access-stores this vertex's neares LU parent
	int[] edge = null, edgeP = null;	// edge = parent->child edges, edgeP = child->parent edges
	public FrontalMatrix fm = null;
	// index to chief pivot of supernode, -1 if not a supernode, superChild is the last pivot of supernode and is referred to by the parents
	int superPivot = -1, superChild = -1;
	int parentsL = 0, parentsU = 0;		// L & U parents count
	
	
	FrontalDAGVertex(int i) {
		this.i = i;
		edge = new int[FrontalDAG.DAG_ALLOCBLOCK];
	}

	FrontalDAGVertex(int i, int[] edge) {
		this.i = i;
		this.edge = edge;
	}

	@Override
	public FrontalDAGVertex clone() {
		// first clone the static data of this vertex
		Object O = null;
		try { O = super.clone(); } catch (CloneNotSupportedException e) { e.printStackTrace(); }
		FrontalDAGVertex vtx = (FrontalDAGVertex) O;
		vtx.edge = vtx.edge == null ? null : edge.clone();
		vtx.edgeP = vtx.edgeP == null ? null : edgeP.clone();
		// frontal matrix is not initialised at DAG reduction stages, assumed null
		return vtx;
	}
	
	
	
	// method tells what type of child this is of vertex vtx
/*	public int whatChildOf(FrontalDAGVertex vtx) {
		for (int e = 0, eEnd = edges + cutedges; e < eEnd; e++) {
			if (edge(e) == vtx.i)
				return edge[e] & LU_CHILD;
		}
		return -1;
	}*/
	
	
	// method takes a vertex index and returns it's offset in a supplied edge array, or -1 if not found
	public int edgeIndex(int v) {
		v &= 0xFFFFFFFF - FLAGS;
		for (int e = 0, eEnd = edges + cutedges; e < eEnd; e++)
			if (edge(e) == v) return e;
		return -1;
	}

	// method takes a vertex index and returns it's offset in edge array, or -1 if not found
	// method works on a sorted array, interrupting search on exceeding/underrunning the current index
	// TODO: make an element-count based decision between iterative and linear search instead
	public int edgeIndexS(int v, int start, boolean ascending) {
		v &= 0xFFFFFFFF - FLAGS;
		if (ascending) {
			for (int e = start, eEnd = edges + cutedges; e < eEnd && edge(e) <= v; e++)
				if (edge(e) == v) return e;
		} else {
			for (int e = start, eEnd = edges + cutedges; e < eEnd && edge(e) >= v; e++)
				if (edge(e) == v) return e;
		}
		return -1;
	}

	
	// method fills in the holes created by cut edges, progressing from right
	// optimally, this method is called when more than half of the left part of array constitutes holes
	public void purgeFlagged() {
		for (int e = 0, e2 = 0; e < edges + cutedges; e++)
			if (edge[e] != -1) edge[e2++] = edge[e];
		cutedges = 0;
	}
	
	
	// method test if edge array needs expanding and is called with the requested number of additional elements
	boolean updateEdgeArray(int i) {
		
		int totaledges = edges + i, alength = edge.length;
		if (totaledges >= alength) {
			int[] edge2 = new int[trimToAllocBlock(totaledges)];
			
			if (cutedges > 0) {							// copy over existing edges, skip over cut edge "holes" assigned to -1
				for (int e = 0, e2 = 0, eEnd = edges + cutedges; e < eEnd; e++)
					if (edge[e] != -1)
						edge2[e2++] = edge[e];
			} else {									// straight copy-over of existing edges
				for (int e = 0; e < edges; e++) edge2[e] = edge[e];
			}
			edge = edge2;
			return true;
		}		
		return false;					// no resize needed, return false
	}
	
//	boolean updateEdgeArrayP(int i) {
//		
//		int totaledges = parentsU + i, alength = edge.length;
//		if (totaledges >= alength) {
//			int[] edgeP = new int[trimToAllocBlock(totaledges)];
//			
//			for (int e = 0; e < edges; e++) edge2[e] = edge[e];
//
//			edge = edge2;
//			return true;
//		}		
//		return false;					// no resize needed, return false
//	}
	

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			HELPER METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	int edge(int e) { return edge[e] & (0xFFFFFFFF - FLAGS); }
	int edgeP(int e) { return edgeP[e] & (0xFFFFFFFF - FLAGS); }
	
	public boolean is_edgeType(int e, int type) { return ((edge[e] & FLAGS) == type); }
	
	// the bits 30 & 29 of edges flag for L/U/LU status of an edge
	public boolean is_L_edge(int e) { return ((edge[e] & FLAGS) == L_CHILD); }
	public boolean is_L_edgeP(int e) { return ((edgeP[e] & FLAGS) == L_CHILD); }
	public void set_L_edge(int e) { edge[e] |= L_CHILD; }
	public boolean is_U_edge(int e) { return ((edge[e] & FLAGS) == U_CHILD); }
	public boolean is_U_edgeP(int e) { return ((edgeP[e] & FLAGS) == U_CHILD); }
	public void set_U_edge(int e) { edge[e] |= U_CHILD; }
	public boolean is_LU_edge(int e) { return ((edge[e] & FLAGS) == LU_CHILD); }
	public boolean is_LU_edgeP(int e) { return ((edgeP[e] & FLAGS) == LU_CHILD); }
	public void set_LU_edge(int e) { edge[e] |= LU_CHILD; }
	// a N-edge has neither L nor U bits set
	public boolean is_N_edge(int e) { return ((edge[e] & FLAGS) == 0); }
	public void set_N_edge(int e) { edge[e] &= 0xFFFFFFFF - LU_CHILD; }
	// supernode test/set
	public boolean is_S_edge(int e) { return ((edge[e] & FLAGS) == SUPERNODE); }
	public void set_S_edge(int e) { edge[e] |= SUPERNODE; }

	private static int trimToAllocBlock(int v) { return (v & (0xFFFFFFFF - FDAGV_ALLOCBLOCK + 1)) + FDAGV_ALLOCBLOCK; }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		
		sb.append(	i + 
					" pLU(" + (parentLU == Integer.MAX_VALUE ? "-" : parentLU) + 
					") ptsLU(" + parentsL + "," + (parentsU - parentsL) +
					"): [");
		for (int e = 0, e2 = 0; e < edges + cutedges; e++) {
			if (edge[e] != -1) {
				if (is_LU_edge(e)) sb.append("LU" + edge(e));
				else if (is_L_edge(e)) sb.append("L" + edge(e));
				else if (is_U_edge(e)) sb.append("U" + edge(e));
				else if (is_S_edge(e)) sb.append("S" + edge(e));
				else sb.append("N" + edge(e));
				if (++e2 < edges) sb.append(", ");
			}
		}
		sb.append("]");
		
		if (cutedges > 0) {
			sb.append(" cut: [");
			for (int ce = edge.length - 1; ce >= edge.length - cutedges; ce--) {
				sb.append(edge[ce]);
				if (ce > edge.length - cutedges) sb.append(", ");
			}
			sb.append("]");
		} 

		if (edgeP != null) {
			sb.append(", edgeP: [");
			for (int eP = 0; eP < parentsU;) {
				if (is_LU_edgeP(eP)) sb.append("LU" + edgeP(eP));
				else if (is_L_edgeP(eP)) sb.append("L" + edgeP(eP));
				else if (is_U_edgeP(eP)) sb.append("U" + edgeP(eP));
				if (++eP < parentsU) sb.append(", ");
			}
			sb.append("]\n");
		} else
			sb.append("\n");
		
		if (fm != null) sb.append("\n" + fm.toString());
		
		return sb.toString();
	}

}
