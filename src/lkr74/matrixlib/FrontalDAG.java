package lkr74.matrixlib;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.InvalidParameterException;

import lkr74.matrixlib.FrontalDAGVertex;

public class FrontalDAG implements Cloneable {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			FRONTAL DAG GRAPH OPERATIONS FOR THE FRONTAL MATRIX LU SOLVER								//
	//			Leonard Krylov 2016																			//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	protected static final int DAG_ALLOCBLOCK = 16;
	
	String name = "DAG";
	int verts = 0, bitSets = 0, edges = 0, subProblems = 0;
	FrontalDAGVertex[] vertex = null;			// the vertices of the DAG
	private long[] visitBits = null;
	// subProblem holds indexes to first vertex of every found subproblem, parentsL & parentsU are accumulative counters for every vertex's L/U-parents
	private int[] subProblem = null, parentsL = null, parentsU = null;
	private boolean distantFirst = true;		// insertion order of vertex edges
	FrontalDAG dataDAG = null;					// keep reference to the dataDAG created from this DAG
	private int allocBlockExtras = 4;
	int[][] mutator;
	
	// skeleton DAG instantiator
	FrontalDAG (String name, int verts) {
		this.name = name;
		this.verts = verts;
		vertex = new FrontalDAGVertex[verts];
		visitBits = visitBitArray(verts);
		subProblem = new int[verts];
		mutator = new int[2][];
		mutator[0] = new int[verts]; mutator[1] = new int[verts];
	}
	
	@Override
	public FrontalDAG clone() {
		// first clone the static data of this DAG
		Object O = null;
		try { O = super.clone(); } catch (CloneNotSupportedException e) { e.printStackTrace(); }
		FrontalDAG dag = (FrontalDAG) O;
		
		dag.vertex = new FrontalDAGVertex[verts];							// new vertexes array
		for (int i = 0; i < verts; i++) dag.vertex[i] = vertex[i].clone();	// clone every vertex in vertexes array
		dag.visitBits = visitBits.clone();
		dag.subProblem = subProblem.clone();
		dag.parentsL = parentsL == null ? null : parentsL.clone();
		dag.parentsU = parentsU == null ? null : parentsU.clone();
		dag.mutator = new int[2][];
		dag.mutator[0] = mutator[0] == null ? null : mutator[0].clone();
		dag.mutator[1] = mutator[1] == null ? null : mutator[1].clone();
		return dag;
	}
	
	
	// constructs a DAG from a NSPMatrix, we are building a parent->child looking DAG, from U part and from L part separately
	// if distantFirst=true, the most distant children will be first in edge arrays, otherwise they will be last
	// note: distantFirst choice will influence how the DAG is transitively reduced, if =true, the long-range edges will prevail
	public static FrontalDAG taskDAG(NSPMatrix M, boolean distantFirst) {
		
		if (M.M != M.N) throw new InvalidParameterException("FrontalDAG.buildDAGfromNSP(): Matrix not square.");
		
		FrontalDAG dagL = new FrontalDAG("TaskDAG_" + M.name, M.M), dagU = new FrontalDAG("", M.M);
		dagL.distantFirst = dagU.distantFirst = distantFirst;

		//////////////// Supernode analysis happens here ////////////////

		// first we do symbolic assembly of all supernodes, testing for membership and assigning the top supernode indexes
		// to save processing time, the frontal supernodes are constructed during analysis and discarded on failure of a pattern match
		dagL.vertex[0] = new FrontalDAGVertex(0);

		int vEnd = dagL.verts - 1;
		for (int v = 0, superCount = 0; v <= vEnd; v++) {
			
			if (v + 1 <= vEnd) {
				FrontalDAGVertex vtx = dagL.vertex[v + 1] = new FrontalDAGVertex(v + 1, new int[1]);		// allocate next vertex
	
				// take care of the case of a zero pivot, salvage by inserting a very small value and leave the problem to consecutive pivoting
				// note: this should have been handled already by preconditioning
				if (M.pivotNsp[v] == null)
					// insert a minimal pivot value that will pass intrapivoting = (1-norm of matrix) * (smallest allowed nonzero)
					M.valueTo2(v, Matrix.ROUNDOFF_ERROR * (M.norm1 < 0 ? M.norm1() : M.norm1));	
	
				// on looking for a candidate pivot for new supernode, skip vertices without L&U-ancestors
				if (superCount == 0 &&
					(M.pivotNsp[v].offH >= M.Hsp[v].nodes - 1 ||
					M.pivotNsp[v].offV >= M.Vsp[v].nodes - 1))
					continue;

				// take care of the case of NEXT pivot also being zero, salvage by inserting a very small value
				// note: this should have been handled already by preconditioning
				int vNext = v + 1;
				if (M.pivotNsp[vNext] == null)
					// insert a minimal pivot value that will pass intrapivoting
					M.valueTo2(vNext, Matrix.ROUNDOFF_ERROR * (M.norm1 < 0 ? M.norm1() : M.norm1));	

				// check if supernode pivot & this sequential pivot have similar row/column patterns
				FrontalMatrix fm = FrontalMatrix.assembleSupernode(M, v - superCount, vNext, 2, 1f/4f, dagL.vertex[v].fm);
				if (fm != null) {
					// next vertex will only have one child, and that child will be in it's supernode group
					vtx.fm = fm;													// reference frontal matrix in next vertex
					vtx.edge[0] = v | FrontalDAGVertex.SUPERNODE;					// pointing to the lower pivot as a lower supernode member
					vtx.edges = 1;
					dagL.edges++;
					superCount++;
					continue;
				}
			}
			
			// we end up here if either a new supernode assembly failed, or an ongoing assembly failed and needs to be finalised
			// superCount>0 means that there is a superPivot assignment needed in all members
			if (superCount > 0) {
				dagL.vertex[v - superCount].fm = dagL.vertex[v].fm;
				int superPivot = v - superCount;
				while (superCount >= 0)	{
					dagL.vertex[v - superCount].superPivot = superPivot;		// assign supernode pivot index to all members of supernode
					dagL.vertex[v - superCount--].superChild = v;				// assign final supernode pivot index to all members
				}
				superCount = 0;
			}
		}

		// take care of the case of a zero last pivot of matrix
		if (M.pivotNsp[vEnd] == null)
			// insert a minimal pivot value that will pass intrapivoting
			M.valueTo2(vEnd, Matrix.ROUNDOFF_ERROR * (M.norm1 < 0 ? M.norm1() : M.norm1));

		//////////////// DAG-L&U Parent->Child L&U-edge addition happens here ////////////////
		
		FrontalDAGVertex vtxFinal = dagL.vertex[dagL.verts - 1];			// memorise the highest vertex, needed for linking disjoint submatrices
		
		// disjointU stores indexes of vertexes that are disjoint from leading indexes in U
		int[] disjointU = new int[dagL.verts];
		int d = 0;
		boolean superNodeDisjointC = true, checkedSuperNode = true;		// flags for supernode offspring disjointness (leaf supernodes)
		boolean superNodeDisjointA = true;								// flag for supernode ancestor disjointness

		// iterate through M's sparse row arrays, adding edges to the U-DAG
		for (int v = dagL.verts - 1; v > 0; v--) {
			
			NSPArray aHsp = M.Hsp[v], aVsp = M.Vsp[v];
			
			int toPivotH = M.pivotNsp[v].offH, toPivotV = M.pivotNsp[v].offV;
			NSPNode[] bHsp = aHsp.array, bVsp = aVsp.array;
			FrontalDAGVertex vtx = dagL.vertex[v];
			int[] edge = null;
			// superNodeDisjoint = false if ANY offspring to a supernode is found, if loop iterates through a superpivot
			// without the flag ever setting to false, the superpivot block will be resembling a block withinin a BTF matrix
			
			if (vtx.superPivot < 0) {									// if this vertex is not part of a supernode
				vtx.edge = edge = new int[trimToAllocBlock(toPivotH)];	// allocate for it an alloc.block that fits at least it's L-edges
			
//				if (aVsp.nodes - 1 <= toPivotV) {						// does row lack a U-parent? (disjoint submatrix in L?)
//					if (	aHsp.nodes - 1 <= toPivotH && 				// does column lack an L-parent? (disjoint submatrix in U?)
//							vtx != vtxFinal) {							// for final vertex, avoid adding disjoint edge to itself
//						vtxFinal.edge[vtxFinal.edges++] = v;			// add N-edge to last pivot of disjoint matrix to top vertex of L
//						dagL.edges++;
//						dagL.subProblem[dagL.subProblems++] = v;		// mark the completely disjointed submatrix as subproblem
//					}
//					disjointU[d++] = v;
//				}
				
				if (toPivotV == 0) {									// does pivot lack U-offspring? (disjoint submatrix in L?)
					if (	toPivotH == 0 && 							// does pivot lack L-offspring? (disjoint submatrix in U?)
							vtx != vtxFinal) {							// for final vertex, avoid adding disjoint edge to itself
						vtxFinal.updateEdgeArray(1);					// we're adding one N-edge into final vertex, make sure it fits
						vtxFinal.edge[vtxFinal.edges++] = v;			// add N-edge to last pivot of disjoint matrix to top vertex of L
						dagL.edges++;
						dagL.subProblem[dagL.subProblems++] = v;		// mark the completely disjointed submatrix as subproblem
					}
					disjointU[d++] = v;
				}
			} else {
				// reset disjoint flag to true if we're at beginning of a new supernode
				if (vtx.superPivot == vtx.i) { checkedSuperNode = false; superNodeDisjointC = superNodeDisjointA = true; }
				// flag for an N-edge if we checked an entire supernode (N-edge is only added if superNodeDisjointC == true or superNodeDisjointA == true)
				else if (vtx.superChild == vtx.i) checkedSuperNode = true;
				
				//int toSPivotH = M.pivotNsp[vtx.superPivot].offH, toSPivotV = M.pivotNsp[vtx.superPivot].offV;
//				int toSPivotH = M.pivotNsp[v].offH, toSPivotV = M.pivotNsp[v].offV;
//				while (aHsp.array[toSPivotH].c <= vtx.superChild && )
//				if (aVsp.nodes - 1 <= toSPivotV) {							// does superpivot row lack a U-parent? (disjoint supernode in L?)
//					if (aHsp.nodes - 1 <= toSPivotH) {						// does superpivot column lack an L-parent? (disjoint supernode in U?)
//						if(superNodeDisjointA && checkedSuperNode &&		// did all other pivots of supernode lack offspring?
//							vtx != vtxFinal) {								// for final vertex, avoid adding disjoint edge to itself
//							vtxFinal.updateEdgeArray(1);					// we're adding one N-edge into final vertex, make sure it fits
//							vtxFinal.edge[vtxFinal.edges++] = v;			// add N-edge to last pivot of disjoint matrix to top vertex of L
//							dagL.edges++;
//						}
//					} else superNodeDisjointA = false;
//					disjointU[d++] = vtx.superPivot;
//				} else superNodeDisjointA = false;
				
				// for a superpivot, check for children ABOVE and BEFORE the FIRST pivot index of superpivot, for EVERY pivot in superpivot block
				int toSPivotH = M.pivotNsp[v].offH, toSPivotV = M.pivotNsp[v].offV;
				// retract to BEFORE the superpivot block of current pivot/vertex
				while (toSPivotH > 0 && aHsp.array[toSPivotH].c >= vtx.superPivot) toSPivotH--;
				while (toSPivotV > 0 && aVsp.array[toSPivotV].r >= vtx.superPivot) toSPivotV--;

				if (toSPivotV == 0) {										// does superpivot lack U-offspring? (disjoint supernode in L?)
					if (toSPivotH == 0) { 									// does superpivot lack L-offspring? (disjoint supernode in U?)
						if (	superNodeDisjointC && checkedSuperNode &&	// did all other pivots of supernode lack offspring?
								v == vtx.superChild) {						// have we reached last pivot of supernode?
							if (vtx != vtxFinal ) {							// for final vertex, avoid adding disjoint edge to itself
								vtxFinal.updateEdgeArray(1);				// we're adding one N-edge into final vertex, make sure it fits
								vtxFinal.edge[vtxFinal.edges++] = v;		// add N-edge to last pivot of disjoint matrix to top vertex of L
								dagL.edges++;
							}
							dagL.subProblem[dagL.subProblems++] = vtx.superPivot;	// mark the completely disjointed supernode + trailings as subproblem
							
						} 
					} else superNodeDisjointC = false;
					disjointU[d++] = vtx.superPivot;
				} else superNodeDisjointC = false;
				
			}

			//////////////// L-DAG edge insertion happens here ////////////////
			
			int e = 0;
			boolean accumEdges = true;
			if (vtx.superPivot >= 0) {								// if we're adding L-edges to a supernode
				while (toPivotH > 0 && bHsp[toPivotH].c > vtx.superPivot)
					toPivotH--;										// retract column index to before start of supernode pivot index
				if (toPivotH <= 0) accumEdges = false;				// there were no L-edges to this supernode vertex, signal it to edge setter
				else {
					vtx = dagL.vertex[vtx.superPivot];
					vtx.updateEdgeArray(toPivotH);					// we're adding edges to pivot vertex of supernode, make sure they fit
	
					e = vtx.edges;
					edge = vtx.edge;
					vtx.edges += toPivotH;							// DEBUG: set preliminary L-edge count so toString() method displays them
				}
			} else
				vtx.edges = toPivotH;								// DEBUG: set preliminary L-edge count so toString() method displays them
			
			if (accumEdges) {
				// choose how edges are laid out in the edge arrays, distant index (to pivot index) or closest index first
				int e1 = 0;
				if (distantFirst) {
					for (int i = 0; i < toPivotH; i++)	{				// iterate from start of row-wise nodes of L-part of M
						int c = bHsp[i].c;
						if (dagL.vertex[c].superPivot >= 0)				// if edge is to a supernode member, then point it to supernode's final vertex
								edge[e] = dagL.vertex[c].superChild | FrontalDAGVertex.L_CHILD;
						else	edge[e] = c | FrontalDAGVertex.L_CHILD;
						if (e > 0 && edge[e - 1] == edge[e]) continue;	// skip duplicates (edges going to the same supernode)
						e++; e1++;
					}
				} else {
					for (int i = toPivotH - 1; i >= 0; i--)	{			// iterate from end of row-wise nodes of L-part of M
						int c = bHsp[i].c;
						if (dagL.vertex[c].superPivot >= 0)				// if edge is to a supernode member, then point it to supernode's final vertex
								edge[e] = dagL.vertex[c].superChild | FrontalDAGVertex.L_CHILD;
						else	edge[e] = c | FrontalDAGVertex.L_CHILD;
						if (e > 0 && edge[e - 1] == edge[e]) continue;	// skip duplicates (edges going to the same supernode)
						e++; e1++;
					}
				}
				vtx.edges = e;
				dagL.edges += e1;
			}
			
			//////////////// U-DAG edge insertion happens here ////////////////
			
			int superPivot = vtx.superPivot;
			if (dagU.vertex[v] == null) dagU.vertex[v] = new FrontalDAGVertex(v, new int[1]);
			vtx = dagU.vertex[v];
			if (superPivot < 0) 									// if this vertex is not part of a supernode, assign regular edge allocblock
				vtx.edge = edge = new int[trimToAllocBlock(toPivotV)];

			e = 0;
			accumEdges = true;
			if (superPivot >= 0) {									// if we're adding U-edges to a supernode
				
//				if (v != superPivot) {								// for correct transitional reduction, make edges between members in supernode
//					vtx.edge[0] = (v - 1) | FrontalDAGVertex.SUPERNODE;	
//					vtx.edges++;
//				}
				while (toPivotV > 0 && bVsp[toPivotV].r > superPivot)
					toPivotV--;										// retract row index to before start of supernode pivot index
				if (toPivotV <= 0) accumEdges = false;				// there were no U-edges to this supernode vertex, signal it to edge setter
				else {
					vtx = dagU.vertex[superPivot];
					if (vtx == null)								// create superpivot vertex also in U-DAG, if it doesn't exist
						vtx = dagU.vertex[superPivot] = new FrontalDAGVertex(superPivot);
					vtx.updateEdgeArray(toPivotV);					// we're adding edges to pivot vertex of supernode, make sure they fit
	
					e = vtx.edges;									// if we already have edges, start adding from last index
					edge = vtx.edge;
					vtx.edges += toPivotV;							// DEBUG: set preliminary U-edge count so toString() method displays them
				}
			} else {
				vtx.edges = toPivotV;								// DEBUG: set preliminary U-edge count so toString() method displays them
			}
			
			if (accumEdges) {
				int e1 = 0;
	 			if (distantFirst) {
					for (int i = 0; i < toPivotV; i++) {				// iterate from start of row-wise nodes of L-part of M
						int r = bVsp[i].r;
						if (dagU.vertex[r].superPivot >= 0)				// if edge is to a supernode member, then point it to supernode's final vertex
								edge[e] = dagU.vertex[r].superChild | FrontalDAGVertex.U_CHILD;
						else	edge[e] = r | FrontalDAGVertex.U_CHILD;
						if (e > 0 && edge[e - 1] == edge[e]) continue;	// skip duplicates (edges going to the same supernode)
						e++; e1++;
					}
				} else {
					for (int i = toPivotV - 1; i >= 0; i--)	{			// iterate from end of row-wise nodes of L-part of M
						int r = bVsp[i].r;
						if (dagL.vertex[r].superPivot >= 0)				// if edge is to a supernode member, then point it to supernode's final vertex
								edge[e] = dagL.vertex[r].superChild | FrontalDAGVertex.U_CHILD;
						else	edge[e] = r | FrontalDAGVertex.U_CHILD;
						if (e > 0 && edge[e - 1] == edge[e]) continue;	// skip duplicates (edges going to the same supernode)
						e++; e1++;
					}
				}
	 			vtx.edges = e;											// update to the final edge count without duplicates
				dagU.edges += e1;
			}
		}
		dagU.vertex[0] = new FrontalDAGVertex(0, new int[1]);
		
		// transitively reduce L-DAG (will reduce everything in a single call as disjoint submatrices are linked through vertexFinal)
		int[] vStack = new int[dagL.verts];
		dagL.transitiveReduction(dagL.verts - 1, vStack, false);
		
		// transitively reduce disjoint submatrices of U-DAG from the highest downwards
		for (int d1 = 0; d1 < d; d1++)	{
			int v = disjointU[d1];
			// transitively reduce this submatrix from every disjoint vertex down, do not clear visitation list
			if (!dagU.visited(v))
				dagU.transitiveReduction(v, vStack, false);
		}
		
		FrontalDAG taskDAG = dagL.uniteEdgeSetLU(dagU);					// unite U-DAG & L-DAG edgesets (also calculates parentsL & parentsU)
		//taskDAG.toGraphViz(true, false, true);
		
		// convert to LU-edges the direct edges i-L->j that also have a U-path i~U~>j, and
		// convert to LU-edges the direct edges i-U->j that also have an L-path i~L~>j
		for (int v = 0; v < taskDAG.verts; v++) {
			FrontalDAGVertex vtx = taskDAG.vertex[v];
			for (int e = 0, v2 = 0, eEnd = vtx.edges + vtx.cutedges; e < eEnd; e++) {
				if (	vtx.is_L_edge(e) &&
						taskDAG.pathTo(v, v2 = vtx.edge(e), FrontalDAGVertex.U_PATH, 0, vStack, true) != null) {
					vtx.set_LU_edge(e);
					if (taskDAG.vertex[v2].parentLU > v)
						taskDAG.vertex[v2].parentLU = v;
					taskDAG.vertex[v2].parentsL--;				// the child stops being L-child becoming a LU-child
				}
				else if(vtx.is_U_edge(e) &&
						taskDAG.pathTo(v, v2 = vtx.edge(e), FrontalDAGVertex.L_PATH, 0, vStack, true) != null) {
					vtx.set_LU_edge(e);
					if (taskDAG.vertex[v2].parentLU > v)
						taskDAG.vertex[v2].parentLU = v;
					taskDAG.vertex[v2].parentsU--;				// the child stops being U-child becoming an LU-child
				}
			}
		}
		
		//taskDAG.toGraphViz(true, false, true);				// output to GraphViz graph visualiser, will also execute GrapViz utility
		// get the approximate block allocation size for extras index contributions, calculated from the factor of nonzeroes in the matrix
		taskDAG.allocBlockExtras = trimToPowerOfTwo((int)((long)(M.M * M.nNZ) / (M.M * M.N)));
		return taskDAG;
	}
	
	
	// method applies Anshul Gupta's criterions for all i-U->j and i->L->j edges, i<j, that allows non-failing transfer
	// of contribution data upwards in a DAG called non-pivoting Data-DAG, which is an edge-extended version of the Task-DAG
	public FrontalDAG dataDAG(NSPMatrix M) {

		FrontalDAG dataDAG = this.clone();
		int[] parentsL = dataDAG.parentsL = new int[verts], parentsU = dataDAG.parentsU = new int[verts];
		
		// child->parent edges need to be added (in dataDAG) for more effectively pursuing the following heuristics
		for (int v = 0; v < verts; v++) {
			
			// we can skip over internal supernode pivots directly to the final/superchild pivot
			FrontalDAGVertex vtx = dataDAG.vertex[v = dataDAG.superChild(v)];
			
			vtx.parentsU += vtx.parentsL;					// readjust parentsU to become an offset for U-edges in array edgeP
			parentsL[v] = vtx.parentsL;						// save parentsL & parentsU counters, will be needed in next section
			parentsU[v] = vtx.parentsU;
			// allocate for the expected L&U-parents of current vertex, fill-in will happen when this vertex's parents are iterated
			if (vtx.parentsU > 0)
				vtx.edgeP = new int[vtx.parentsU];
			// a child's edgeP will be a lower/upper halves array with first half holding L-parent edges, second holding U-parent edges
			// the child's parentsL & parentsU counters will be used as steppers during edge addition, and will thus be destroyed
			for (int e = 0, eEnd = vtx.edges + vtx.cutedges; e < eEnd; e++) {
				if (vtx.edge[e] == -1) continue;
				int v2 = vtx.edge(e);
				if (vtx.is_L_edge(e))
					dataDAG.vertex[v2].edgeP[--dataDAG.vertex[v2].parentsL] = v | FrontalDAGVertex.L_CHILD;
				else if (vtx.is_U_edge(e))
					dataDAG.vertex[v2].edgeP[--dataDAG.vertex[v2].parentsU] = v | FrontalDAGVertex.U_CHILD;
			}
		}
		
		// reinstate the vertex-localised parentsL/parentsU counters
		FrontalDAGVertex[] dDAGvertex = dataDAG.vertex;
		for (int v = 0; v < verts; v++) { dDAGvertex[v].parentsL = parentsL[v]; dDAGvertex[v].parentsU = parentsU[v]; }
		int[] vStack = new int[verts];
		
		// TODO: sort out the correct order of all the heuristics that follow below

		// this loop starts converting the nonpivoting data DAG into a pivoting data DAG by following edge additions:
		// for each vertex/supernode g, find the smallest vertex/supernode h such that g~L~>h & g~U~>h exist and if found:
		// add an g-LU->h edge, OR if a g-L->h or g-U->h edge exists, convert it, then, for all k > h delete edges g->k,
		// on this modification, each supernode except root will have an LU-parent to accommodate failed row/column pairs
		for (int v = 0, s = subProblems - 1, subEnd = 0; v < verts; v++) {
			
			if (s <= 1) subEnd = verts;											// find the boundaries of this subproblem to search-iterate over
			else if (v >= subProblem[s])
					subEnd = subProblem[--s];
			else	subEnd = subProblem[s];
			
			// note: supernodes automatically caretaken since only superchild pivot/vertex will have nonzero parentsL & parentsU
			if (parentsL[v] > 0 && parentsU[v] != parentsL[v]) {				// if child has both an U-parent and L-parent
				// use visit flags for the search of smallest LU-parent of vertex v, since we don't want to revisit unattainable paths to v
				// at the same time, the pathTo() method will be called with clearVisits = false, to mark those unattainable paths
				clearVisits();
				boolean clearAboveH = false;
				int vH = 0;
				for (int v2 = v + 1; v2 < subEnd; v2++) {						// then check all vertices a step above v within current subproblem
					
					// note: supernodes automatically caretaken since only superparent pivot/vertex will have edges towards children
					if (vertex[v2].edges < 2) continue;							// vertex v2 must have both an U-child and L-child
					FrontalDAGVertex vtx2 = vertex[v2];
					int eLdirect = -1, eUdirect = -1, eLfound = -1, eUfound = -1;
					
					for (int e = 0, eEnd = vtx2.edges + vtx2.cutedges; e < eEnd; e++) {
						if (vtx2.edge[e] == -1) continue;						// skip erased edges
						if (vtx2.is_L_edge(e)) {
							eLfound = e;
							if (vtx2.edge(e) == v) eLdirect = eLfound;			// flag the case of finding a g-L->h (there can only be one)
							if (eLdirect >= 0 && eUfound >= 0) break;			// on finding a g-L->h and at least one U-child, stop
						}
						else if (vtx2.is_U_edge(e)) {
							eUfound = e;
							if (vtx2.edge(e) == v) eUdirect = eUfound;			// flag the case of finding a g-U->h (there can only be one)
							if (eUdirect >= 0 && eLfound >= 0) break;			// on finding a g-U->h and at least one L-child, stop
						}
					}
					
					if (eLfound >= 0 && eUfound >= 0) {
						// if we found a g-L->h and at least one U-child edge, we only need to check the g~U~>h case
						if (eLdirect >= 0 && pathTo(v2, v, FrontalDAGVertex.U_PATH, eUfound, vStack, false) != null) {
							dataDAG.vertex[v2].set_LU_edge(eLdirect);			// convert g-L->h to g-LU->h
							clearAboveH = true;									// flag for removal of all edges g->k for all k > h
							vH = v2; break;
						} else
						// if we found a g-U->h and at least one L-child edge, we only need to check the g~L~>h case
						if (eUdirect >= 0 && pathTo(v2, v, FrontalDAGVertex.L_PATH, eLfound, vStack, false) != null) {
							dataDAG.vertex[v2].set_LU_edge(eUdirect);			// convert g-L->h to g-LU->h
							clearAboveH = true;									// flag for removal of all edges g->k for all k > h
							vH = v2; break;
						} else
						if (pathTo(v2, v, FrontalDAGVertex.L_PATH, eLfound, vStack, false) != null &&	// if both L-path & U-path between vertex v & v2 exists
							pathTo(v2, v, FrontalDAGVertex.U_PATH, eUfound, vStack, false) != null) {
							clearAboveH = true;
							dataDAG.addEdge(v2, v, FrontalDAGVertex.LU_CHILD, true);
							vH = v2; break;
						}
					}
				}
				
				if (clearAboveH) {
					int[] edgeP = dataDAG.vertex[v].edgeP;
					for (int e = 0, eEnd = parentsU[v]; e < eEnd; e++) {
						int vK = edgeP[e] & (0xFFFFFFFF - FrontalDAGVertex.FLAGS);
						if (vK > vH)														// if some k > h
						dataDAG.cutEdge(vH, v, e);											// cut g->k
					}
				}
			}
		}

		for (int v = 0; v < verts; v++) {
			int vSC = dataDAG.superChild(v), pLU = vertex[vSC].parentLU;
			FrontalDAGVertex vtx = dataDAG.vertex[v], vtxSC = dataDAG.vertex[vSC];
			
			// a heuristic that tests every pivot for three cases, which, if fulfilled, demand adding an U-edge i-U->j that
			// will guarantee ancestor absorption of all contribution matrices in a non-pivoting data-DAG:
			// for element j in Struct(U[i,*]) / Struct(L[*,i]):
			// 1) if LU-parent of i (if existing) is greater than j
			// 2) if none of i's U-parents are in Struct(U[*,j]) / i's L-parents are in Struct(L[j,*])
			// 3a) if a k exists in Struct(L[*,i]) / Struct(U[i,*]) such that k > j OR
			// 3b) if a path k~U~>i / k~L~>i exists in Task-DAG AND
			// 3c) if LU-parent(k) > j

			NSPArray aHsp = M.Hsp[v];
			NSPNode[] bHsp = aHsp.array;
			int i = M.pivotNsp[v].offH + 1, pntsL = parentsL[vSC], pntsU = parentsU[vSC];
			if (vtx.superPivot >= 0)															// if we're dealing with a supernode
				while (i < aHsp.nodes && bHsp[i].c < vtx.superChild) i++;						// go outside supernode pivot block boundary
			
			for (; i < aHsp.nodes; i++) {														// for every element j in Struct(U[i,*])
				int j = bHsp[i].c;
				if (j < pLU) {																	// if LU-parent of i (if existing) is greater than j
					
					// for every U-parent of i
					boolean inStructJ = false;
					NSPArray aVsp = M.Vsp[j];
					NSPNode[] bVsp = aVsp.array;
					int r2 =  M.pivotNsp[j].offV + 1;
					for (int p = pntsL; p < pntsU; p++) {
						
						for (int r = r2; r < aVsp.nodes; r++)									// if none of i's U-parents are in Struct(U[*,j])
							if (bVsp[r].r == vtxSC.edgeP(p)) {
								inStructJ = true; break;										// a U-parent found in Struct(U[*,j])...
							}
						if (inStructJ) break;													// ...try next j in Struct(U[i,*])				
					}
					
					if (!inStructJ) {															// (1) and (2) succeeded, try last test (3abc)
						aVsp = M.Vsp[v];
						bVsp = aVsp.array;
						r2 = M.pivotNsp[i].offV + 1;
						if (vtx.superPivot >= 0)												// if we're dealing with a supernode
							while (r2 < aVsp.nodes && bVsp[r2].r < vtx.superChild) r2++;		// go outside supernode pivot block boundary

						for (; r2 < aVsp.nodes; r2++) {
							int k = bVsp[r2].r;
							if (k > j ||														// if a k exists in Struct(L[*,i]) such that k > j
								(dataDAG.vertex[k].parentLU > j && i > k &&						// or LU-parent(k) > j
								dataDAG.pathTo(i,k,FrontalDAGVertex.U_PATH,0,vStack,true) != null)) {	// and a path k~U~>i exists
								dataDAG.addEdge(												// then add i-U->j (if it didn't exist)
										dataDAG.superParent(j),									// and pay attention to existing supernodes when adding edge
										dataDAG.superChild(v), FrontalDAGVertex.U_CHILD, true);
								break;
							}
						}
					}
				}
			}
			
			// do the same heuristic for adding L-edge i-L->j from a transposed aspect			
			NSPArray aVsp = M.Vsp[v];
			NSPNode[] bVsp = aVsp.array;
			i = M.pivotNsp[v].offV + 1;
			if (vtx.superPivot >= 0)													// if we're dealing with a supernode
				while (i < aVsp.nodes && bVsp[i].r < vtx.superChild) i++;				// go outside supernode pivot block boundary

			for (; i < aVsp.nodes; i++) {												// for every element j in Struct(L[*,i])
				int j = bVsp[i].r;
				if (j < pLU) {															// if LU-parent of i (if existing) is greater than j
					
					// for every L-parent of i
					boolean inStructJ = false;
					aHsp = M.Hsp[j];
					bHsp = aHsp.array;
					int c2 = M.pivotNsp[j].offH + 1;
					for (int p = 0; p < pntsL; p++) {
						
						for (int c = c2; c < aHsp.nodes; c++)							// if none of i's L-parents are in Struct(L[j,*])
							if (bHsp[c].c == vtxSC.edgeP(p)) {
								inStructJ = true; break;								// L-parent was found in Struct(L[j,*])...
							}
						if (inStructJ) break;											// ...try next j in Struct(L[*,i])					
					}
					
					if (!inStructJ) {
						// (1) and (2) succeeded, try last test (3)
						aHsp = M.Hsp[v];
						bHsp = aHsp.array;
						c2 = M.pivotNsp[i].offH + 1;
						if (vtx.superPivot >= 0)												// if we're dealing with a supernode
							while (c2 < aHsp.nodes && bHsp[c2].c < vtx.superChild) c2++;		// go outside supernode pivot block boundary

						for (; c2 < aHsp.nodes; c2++) {											// if a k exists in Struct(U[i,*]) such that k > j
							int k = bHsp[c2].c;
							if (k > j ||														// if a k exists in Struct(U[i,*]) such that k > j
								(dataDAG.vertex[k].parentLU > j && i > k &&								// or LU-parent(k) > j
								dataDAG.pathTo(i,k,FrontalDAGVertex.L_PATH,0,vStack,true) != null)){	// and a path k~L~>i exists in Task-DAG
								dataDAG.addEdge(												// then add i-L->j (if it didn't exist)
										dataDAG.superParent(j),									// and pay attention to existing supernodes when adding edge
										dataDAG.superChild(v), FrontalDAGVertex.L_CHILD, true);
								break;
							}
						}						
					}
				}
			}
		}

		for (int v = 0; v < verts; v++) {
			int vSC = dataDAG.superChild(v), pLU = dataDAG.vertex[vSC].parentLU;
			//FrontalDAGVertex vtxSP = dataDAG.vertex[dataDAG.superParent(v)];
			
			// a heuristic that tests for cases when intersupernodal pivoting creates fill-ins
			// that must be taken care of with extra contribution inheritance edges
			// if h is the LU-parent of j and all of following conditions hold, then i-U->h / i-L->h
			// is necessary for complete assembly of contribution C(i) into it's parental frontal matrices,
			// in the event that j fails meeting the pivot criterion in it's original position:
			// 1) there exists an L-path j~L~>i / U-path j~U~>i such that i < h and LU-parent(i) > h
			// 2) none of i's U-parents are in Struct(L[*,j]) / L-parents are in Struct(U[j,*])
			// 3) either there's a k in Struct(L[*,i]) / Struct(U[i,*]) where k > h or a path k~U~>i / k~L~>i and LU-parent(k) > h
			// ((3) is not really tested according to [A.Gupta], an edge is added anyway, a few extra edges will not add too much computation)
			// (and the same heuristic done for the transposed case)
			
			if (pLU < Integer.MAX_VALUE) {												// if h is the LU-parent of j (current vertex) and...
				NSPNode pivNj = M.pivotNsp[v];
				int h = pLU, j = vSC, r1 = pivNj.offV + 1, c1 = pivNj.offH + 1;
				for (int i = j + 1; i < h; i++)	{										// for j < i < h ...
					
					FrontalDAGVertex vtxSCi = dataDAG.vertex[superChild(i)];			// (get superchild vertex of a supernode)
					i = dataDAG.superParent(i);											// (get superparent vertex index of a supernode)
					if (vtxSCi.edges <= 0) continue;									// (skip i if it has no children)
					int pLUi = vtxSCi.parentLU;
					
					if (pLUi <= h) continue;											// 1) if LU-parent(i) > h
					if (i > j &&														// and there exists an L-path j~L~>i such that i < h (for-loop fulfills i < h)
						dataDAG.pathTo(i, j, FrontalDAGVertex.L_PATH, 0, vStack, true) != null) {

						NSPArray aVsp = M.Vsp[v];
						NSPNode[] bVsp = aVsp.array;
						boolean inStructJ = false;
						for (int p = vtxSCi.parentsL; p < vtxSCi.parentsU; p++) {		// 2) none of i's U-parents...
							for (int r = r1; r < aVsp.nodes; r++)						// are in Struct(U[*,j])
								if (bVsp[r].r == vtxSCi.edgeP(p)) {
									inStructJ = true; break;							// (U-parent(i) found in Struct(U[*,j]), abort)
								}
							if (inStructJ) break;							
						}
						if (!inStructJ)
							dataDAG.addEdge(h, dataDAG.superChild(i), FrontalDAGVertex.U_CHILD, true);		// add i-U->h (if it didn't exist)	
					}
					
					// (do the same heuristic for transpose case)
					if (i > j &&														// if there exists an U-path j~U~>i such that i < h ...
						dataDAG.pathTo(i, j, FrontalDAGVertex.U_PATH, 0, vStack, true) != null) {

						NSPArray aHsp = M.Hsp[v];
						NSPNode[]  bHsp = aHsp.array;
						boolean inStructJ = false;
						for (int p = 0; p < vtxSCi.parentsL; p++) {						// 2) none of i's L-parents...
							for (int c = c1; c < aHsp.nodes; c++)						// are in Struct(L[j,*])
								if (bHsp[c].c == vtxSCi.edgeP(p)) {
									inStructJ = true; break;							// (L-parent(i) found in Struct(L[j,*]), abort)
								}
							if (inStructJ) break;							
						}
						if (!inStructJ)
							dataDAG.addEdge(h, dataDAG.superChild(i), FrontalDAGVertex.L_CHILD, true);		// add i-L->h (if it didn't exist)	
					}
				}
			}
			
		}
		
		this.dataDAG = dataDAG;
		dataDAG.clearVisits();
		clearVisits();
		return dataDAG;
	}
	
	
	NSPMatrix decomposeLU(NSPMatrix M) {
		if (M.M != M.N)	throw new InvalidParameterException("FrontalDAG.decomposeLU(): Matrix not square.");
		
		for (int i = 0; i < M.M; i++) mutator[0][i] = mutator[1][i] = i;	// initialise mutators to initial pivot order
		
		// option1: iterate through subproblems at this level
		// if we found separate subproblems in the matrix, process them sequentially (or perhaps entire matrix is one single subproblem)
		// this is the prime place for parallel processing of a factorisation
//		int spStart = verts - 1;
//		for (int sp = 0, spEnd; sp < subProblems; sp++) {
//			spEnd = subProblem[sp];
//			factoriseDAG(vertex[spStart], M, spEnd);						// tell factoriser to stop at depth of "spEnd" (final subproblem index)
//			spStart = spEnd - 1;
//		}
		// option2: we already have stored N-edges to all subproblems in topmost vertex, factoriser takes care of processing thrm
		factoriseDAG(vertex[verts - 1], M, 0);
		
		NSPMatrix LU = new NSPMatrix("fmLU", M.M, M.N);
		LU.getFrontalData(this);
		return LU;
	}

	static int HELPARRAY_SIZE = 6;
	
	// top-level code for multifrontal assembly & factorisation, structured according to Anshul Gupta
	// subPi is an index into the subproblems delineation array
	private void factoriseDAG(FrontalDAGVertex vtx, NSPMatrix M, int subPi) {
		
		while (vtx.i < subProblem[subPi]) subPi++;			// update vertex location within the subproblem indexing
		
		// recursively factor children of vertex in Task-DAG if they weren't factored already (if they haven't got a built-up frontal matrix)
		// note that only the first parent of a child reaching that child performs it's recursive computation
		// the other parents are only picking up it's constituent data later in the algorithm
		// also note that only the supernode's final pivot vertex will keep edges to children, thus internal supernode edges are bypassed
		for (int v, e = 0; e < vtx.edges; e++) {
			v = vtx.edge(e);
			if (!visited(v)) {								// only first parent to reach this child will factor it
				visit(v);									// signal to other parents that this vertex has been processed
				if (vtx.is_S_edge(e)) {
					// skip supernodal structures, if this is top vertex's supernode, signal that this is final vertex of subproblem
					factoriseDAG(vertex[vtx.superPivot], M, subPi);
					// if this was a superpivot branch and not the topmost vertex, then this vertex is factorised now, we're done
					// the topmost vertex CAN be a superpivot branch, but can ALSO hold subproblems in N-edges to be solved
					if (vtx.superPivot != -1 && vtx.i != verts - 1) return;	
				}
				else factoriseDAG(vertex[v], M, subPi);		// this an N-edge or an L/U/LU-edge
			}
		}
		if (vtx.i == verts - 1) return;						// for top vertex, all children done means all matrix subproblems are done
		
		// we reach this spot EITHER if this is not the top vertex, and is a leaf vertex OR all children have been processed
		// find out whether we're processing topmost vertex/supernode of current subproblem
		int subPfinal = subPi - 1 < 0 ? verts - 1 : subProblem[subPi - 1] - 1;
		boolean isFinal = vtx.superChild == subPfinal || vtx.i == subPfinal;
		
		if (vtx.edges == 0) {											// if this is a leaf vertex, construct and factorise immediately
			if (vtx.fm == null)											// construct only if it isn't a preconstructed supernode
				vtx.fm = new FrontalMatrix(M, vtx.i, 0, 0, 0);			// leaf basecase - new root with pivot row&column initialised
			vtx.fm.decomposeLU(isFinal);								// factorise root frontal matrix (flag if this is final pivot of a subproblem)
			return;
		}
		
		FrontalDAGVertex vtxD = dataDAG.vertex[vtx.i];
		// data structures for intermixing extras indexes from children for addition to root and for intermixing contribution indexes from children
		// the first two arrays are: No.0 = the global row or column indexes, No.1 = the index pairings of No.0
		// third array is a help-array holding data counts and some optional flags/counts/parameters needed/produced by joinIndexes() method
		int[][] extraJoinDataC = new int[3][], extraJoinDataR = new int[3][];
		int[][] contribJoinDataC = new int[3][], contribJoinDataR = new int[3][];
		int[][] extraJoinDataP = new int[3][], joinData = new int[3][];
		int[] eColsC = null, eColsR = null;
		boolean gotContribution = false, gotExtras = false;

		// collect contributions/extras/failed pivots from root's children, join into appropriate arrays
		for (int e = 0; e < vtxD.edges; e++) {
			
			int childType = vtxD.edge[e] & FrontalDAGVertex.LU_CHILD;
			if (childType == FrontalDAGVertex.N_CHILD) continue;		// N-edges (subproblems) have no parents to contribute to
			
			FrontalMatrix fmC = vertex[vtxD.edge(e)].fm;
			int[] idxCR = fmC.idxCR;
			boolean addExtraColumns = false, addExtraRows = false;
			
			// if we're L-Parent to this child
			if (childType == FrontalDAGVertex.L_CHILD) {
				
				// if L-child has failed pivots
				if (fmC.nFailed > 0) {
					fmC.setLcontributor();													// signal that child has contribution to make
					// add them in sorted order to the sorted list of extra columns
					if (extraJoinDataC[2] == null) extraJoinDataC[2] = new int[HELPARRAY_SIZE];
					extraJoinDataC[2][2] = fmC.nFailed;
					// flag for adding unsorted indexes into sorted order, them being unsorted because of possible intrapivoting in child
					extraJoinDataC[2][4] = 4;
					joinData[0] = fmC.failedP;
					int oldCount = extraJoinDataC[2][0];
					joinIndexes(extraJoinDataC, joinData);
					if (!gotExtras && extraJoinDataC[2][0] > oldCount) gotExtras = true;	// flag if any contribution was made
				}

				if (fmC.nContribR > 0) {
					fmC.setLcontributor();													// signal that child has contribution to make
					// calculate in expected contrib.indexes from rows of child, skipping indexes before pivot's of root
					int toC = fmC.cTot + fmC.nExtraP + fmC.nPivots, toCend = toC + fmC.nContribR;
					while (idxCR[toC] < vtx.i) toC++;										// skip child's row indexes < parent's pivot
					if (contribJoinDataR[2] == null) contribJoinDataR[2] = new int[HELPARRAY_SIZE];
					contribJoinDataR[2][2] = toCend - toC;									// how many contribution indexes from child
					contribJoinDataR[2][3] = toC;											// the offset to contribution indexes in child
					joinData[0] = idxCR;
					int oldCount = contribJoinDataR[2][0];
					joinIndexes(contribJoinDataR, joinData);
					if (!gotContribution && contribJoinDataR[2][0] > oldCount) gotContribution = true;	// flag if any contribution was made
				}
				addExtraRows = true;
				
			} else if (childType == FrontalDAGVertex.U_CHILD) {

				//fmC.setLackingContribution();
				// if U-child has failed pivots
				if (fmC.nFailed > 0) {
					fmC.setUcontributor();													// signal that child has contribution to make
					// add them in sorted order to the sorted list of extra rows
					if (extraJoinDataR[2] == null) extraJoinDataR[2] = new int[HELPARRAY_SIZE];
					extraJoinDataR[2][2] = fmC.nFailed;
					// flag for adding unsorted indexes into sorted order, them being unsorted because of possible intrapivoting in child
					extraJoinDataR[2][4] = 4;
					joinData[0] = fmC.failedP;
					int oldCount = extraJoinDataR[2][0];
					joinIndexes(extraJoinDataR, joinData);
					if (!gotExtras && extraJoinDataR[2][0] > oldCount) gotExtras = true;	// flag if any contribution was made
				}
				
				// calculate in expected contrib.indexes from columns of child, skipping indexes before pivot's of root
				if (fmC.nContribC > 0) {
					fmC.setUcontributor();													// signal that child has contribution to make
					int toC = fmC.nExtraP + fmC.nPivots, toCend = toC + fmC.nContribC;
					while (idxCR[toC] < vtx.i) toC++;										// skip child's column indexes < parent's pivot
					if (contribJoinDataC[2] == null) contribJoinDataC[2] = new int[HELPARRAY_SIZE];
					contribJoinDataC[2][2] = toCend - toC;									// how many contribution indexes from child
					contribJoinDataC[2][3] = toC;											// the offset to contribution indexes in child
					joinData[0] = idxCR;
					int oldCount = contribJoinDataC[2][0];
					joinIndexes(contribJoinDataC, joinData);
					if (!gotContribution && contribJoinDataC[2][0] > oldCount) gotContribution = true;	// flag if any contribution was made
				}
				addExtraColumns = true;
				
			} else if (childType == FrontalDAGVertex.LU_CHILD) {

				//fmC.setLackingContribution();
				// if LU-child has failed pivots
				if (fmC.nFailed > 0) {
					fmC.setLUcontributor();													// signal that child has contribution to make
					// add them to the sorted list of extra pivots
					if (extraJoinDataP[2] == null) extraJoinDataP[2] = new int[HELPARRAY_SIZE];
					extraJoinDataP[2][2] = fmC.nFailed;
					// flag unsorted addition of possibly pivoted index duplets, as assymetric pivoting will disassociate
					// row & column indexes, we need to store global row & column index and datafield offset for every extra pivot
					extraJoinDataP[2][4] = 3;
					joinData[0] = fmC.failedP;
					int oldCount =  extraJoinDataP[2][0];
					joinIndexes(extraJoinDataP, joinData);
					if (!gotExtras && extraJoinDataP[2][0] > oldCount) gotExtras = true;	// flag if any contribution was made
				}
				
				// calculate in expected contrib.indexes from rows and columns of child
				int toC, toCend, p;
				if (fmC.nContribR > 0) {
					fmC.setLUcontributor();													// signal that child has contribution to make
					toC = fmC.cTot + fmC.nExtraP + fmC.nPivots; toCend = toC + fmC.nContribR; p = vtx.i;
					while (idxCR[toC] < p) toC++;											// skip child's row indexes < parent's pivot
					if (contribJoinDataR[2] == null) contribJoinDataR[2] = new int[HELPARRAY_SIZE];
					contribJoinDataR[2][2] = toCend - toC;									// how many contribution indexes from child
					contribJoinDataR[2][3] = toC;											// the offset to contribution indexes in child
					joinData[0] = idxCR;
					int oldCount = contribJoinDataR[2][0];
					joinIndexes(contribJoinDataR, joinData);
					if (!gotContribution && contribJoinDataR[2][0] > oldCount) gotContribution = true;	// flag if any contribution was made
				}
				if (fmC.nContribC > 0) {
					fmC.setLUcontributor();													// signal that child has contribution to make
					toC = fmC.nExtraP + fmC.nPivots; toCend = toC + fmC.nContribC; p = vtx.i;
					while (idxCR[toC] < p) toC++;											// skip child's column indexes < parent's pivot
					if (contribJoinDataC[2] == null) contribJoinDataC[2] = new int[HELPARRAY_SIZE];
					contribJoinDataC[2][2] = toCend - toC;									// how many contribution indexes from child
					contribJoinDataC[2][3] = toC;											// the offset to contribution indexes in child
					joinData[0] = idxCR;
					int oldCount = contribJoinDataC[2][0];
					joinIndexes(contribJoinDataC, joinData);
					if (!gotContribution && contribJoinDataC[2][0] > oldCount) gotContribution = true;	// flag if any contribution was made
				}
				addExtraColumns = addExtraRows = true;
			}
			
			if (addExtraColumns && fmC.nExtraC > 0) {
				
				// TODO: investigate case when child produces no contribution (potential subproblem?) but still has extra rows/columns?
				//if (fmC.lacksContribution()) fmC.clearLacksContribution();
				
				// if child has extra columns, add those whose L/LU-parent is above current root (frontal pivot)
				// worst-case allocate temp.arrays for the extras columns that might be joined, and also for their pairings
				if (eColsC == null || fmC.nExtraC > eColsC.length)		// see if this child's contrib.indexes will fit in temp.array
					eColsC = joinData[0] = new int[fmC.nExtraC];
				if (eColsR == null || fmC.nExtraC > eColsR.length)		// see if this child's pairing indexes will fit in temp.array
					eColsR = joinData[1] = new int[fmC.nExtraC];
				
				int j1 = 0, j2 = 0;	
				int[] extraC1 = fmC.extraC[0], extraC2 = fmC.extraC[1];
				for (int j = 0, vP = vtx.i; j < fmC.nExtraC; j++, j2++) {	// first gather all valid edges into temporary array eColsC
					int vExt = extraC1[j], vPLU = dataDAG.vertex[vExt].parentLU;
					if (/*vPLU != Integer.MAX_VALUE &&*/ vPLU > vP) {
						eColsC[j1] = vExt;								// copy column index that passed criterion
						eColsR[j1++] = extraC2[j2];						// copy it's pairing
					}
				}
				// join in valid children's extra columns en-bloque
				if (j1 > 0) {
					if (extraJoinDataC[2] == null) extraJoinDataC[2] = new int[HELPARRAY_SIZE];
					extraJoinDataC[2][2] = j1;
					int oldCount = extraJoinDataC[2][0];
					joinIndexes(extraJoinDataC, joinData);
					if (!gotExtras && extraJoinDataC[2][0] > oldCount) gotExtras = true;	// flag if any contribution was made
				}
			}

			if (addExtraRows && fmC.nExtraR > 0) {
				// if child has extra rows, add those whose U/LU-parent is above current root (frontal pivot)
				if (eColsR == null || fmC.nExtraR > eColsR.length) 		// see if this child's contrib.indexes will fit in temp.array
					eColsR = joinData[0] = new int[fmC.nExtraR];
				if (eColsC == null || fmC.nExtraR > eColsC.length)		// see if this child's pairing indexes will fit in temp.array
					eColsC = joinData[1] = new int[fmC.nExtraR];
				
				int i1 = 0, i2 = 0;	
				int[] extraR1 = fmC.extraR[0], extraR2 = fmC.extraR[1];
				for (int i = 0, vP = vtx.i; i < fmC.nExtraR; i++, i2++) {	// first gather all valid edges into temporary array eColsC
					int vExt = extraR1[i], vPLU = dataDAG.vertex[vExt].parentLU;
					if (/*vPLU != Integer.MAX_VALUE &&*/ vPLU > vP) {
						eColsR[i1] = vExt;
						eColsC[i1++] = extraR2[i2];
					}
				}
				// join in valid children's extra columns en-bloque
				if (i1 > 0) {
					if (extraJoinDataR[2] == null) extraJoinDataR[2] = new int[HELPARRAY_SIZE];
					extraJoinDataR[2][2] = i1;
					int oldCount = extraJoinDataR[2][0];
					joinIndexes(extraJoinDataR, joinData);
					if (!gotExtras && extraJoinDataR[2][0] > oldCount) gotExtras = true;	// flag if any contribution was made
				}
			}			
		}
		
		// assemble the frontal matrix from what we have gathered			
		// set flag for calculation of offsets for root's indexes into collective indexation
		int[] contribC = contribJoinDataC[2], contribR = contribJoinDataR[2];
		boolean hasContribC = true, hasContribR = true;
		if (contribR == null) { contribJoinDataR[2] = contribR = new int[HELPARRAY_SIZE]; hasContribR = false; }
		if (contribC == null) { contribJoinDataC[2] = contribC = new int[HELPARRAY_SIZE]; hasContribC = false; }
		
		// if this is a LU-parent of some child/children with failed pivots, then we need to move the extras indexes for rows/columns
		// of those children that must be moved into extra pivots, the LU-parent will add them to it's pivot block
		if (gotExtras)
			moveLUparentExtraPivots(vtx.i, extraJoinDataC, extraJoinDataR, extraJoinDataP);

		// this block joins in an existing parent's frontal's indexes directly from it's idxCR array
		if (vtx.fm != null) {											// if root non-null, then it's an initialised frontal/supernode
			FrontalMatrix fmP = vtx.fm;
			joinData[0] = fmP.idxCR;
			if (hasContribC) {
				contribC[2] = fmP.nPivots + fmP.nContribC;				// how many pivots + contribution indexes from root
				contribC[3] = fmP.nExtraP;								// the offset to pivots + contribution indexes in root
				contribC[4] = 2;										// we want to get offsets for parent's indexes
				joinIndexes(contribJoinDataC, joinData);
			}
			if (hasContribR) {
				contribR[2] = fmP.nPivots + fmP.nContribR;
				contribR[3] = fmP.nExtraP + fmP.cTot;
				contribR[4] = 2;
				joinData[0] = fmP.idxCR;
				joinIndexes(contribJoinDataR, joinData);
			}
			
			// extend existing root frontal with new contribution indexes (only if any new indexes were added)
			if (gotContribution || gotExtras)
				vtx.fm.extend(extraJoinDataP, contribJoinDataC, contribJoinDataR, extraJoinDataC, extraJoinDataR);
			
			vtx.fm.extraC = extraJoinDataC;								// store down the gathered extras information in root
			vtx.fm.extraR = extraJoinDataR;								// it will be used by the parents of root

		// this block is called when parent's frontal isn't constructed and indexes must be accumulated from the NSPMatrix
		} else {
			int[] rootContribC = joinData[0] = M.indexList(vtx.i, true, 1);		// collect column indexes from NSP matrix, including pivot
			contribC[2] = rootContribC.length;
			contribC[3] = 0;
			contribC[4] = 2;											// we want to get offsets for parent's indexes
			joinIndexes(contribJoinDataC, joinData);
			int[] rootContribR = joinData[0] = M.indexList(vtx.i, false, 1);	// collect row indexes from NSP matrix, including pivot
			contribR[2] = rootContribR.length;
			contribR[3] = 0;
			contribR[4] = 2;											// we want to get offsets for parent's indexes
			joinIndexes(contribJoinDataR, joinData);

			vtx.fm = new FrontalMatrix(
					M, vtx.i, extraJoinDataP[1] == null ? 0 : extraJoinDataP[1][0], 
					contribJoinDataC, contribJoinDataR, extraJoinDataC, extraJoinDataR);
		}
					
		// add in contribution from each child, use offsets in contribJoinData arrays for optimising iteration through parent's indexes
		for (int e = 0; e < vtxD.edges; e++) {
			int v = vtxD.edge(e);
			// process only legitimate children of this parent (that have childFlag properly set, also skips N-edges) 
			if ((vertex[v].fm.childFlag & 3) > 0)
				vtx.fm.add2(vertex[v].fm);
		}
		
		vtx.fm.decomposeLU(isFinal);							// factorise root frontal matrix
	}
	
	
	
	
	// method singles out the extras row & column indexes that come from frontals that are LU-children of currently processed vertex
	// the incoming row/column also has it's row/colum aspect pairing incoming (because of asymmetric pivoting mutating the indexes)
	// therefore we have following cases to take care of:
	// 1) an incoming column that needs to check for it's row half-pairing in the already absorbed failed pivots
	// 2) an incoming column that found no pairing and is added to wait for an incoming row that matches it's pairing
	// 3) an incoming row that needs to check for it's column half-pairing in the already absorbed failed pivots
	void moveLUparentExtraPivots(int parentLU, int[][] eJDC, int[][] eJDR, int[][] eJDP) {

		int[] exPR = eJDP[0], exPC = eJDP[1], helpArrayP = eJDP[2], exPRn = null, exPCn = null;
		int xPSum = helpArrayP == null ? 0 : helpArrayP[0], xPEnd = xPSum;

		if (eJDC[0] != null) {											// do we have any extras columns?
			int[] exR = eJDC[0], exC = eJDC[1], helpArrayC = eJDC[2];
			int xCEnd = helpArrayC[0];

			for (int xC = 0, xPn = 0; xC < xCEnd; xC++) {
				int superC = vertex[exR[xC]].superChild;
				// if the vertex referred to by the column index has root (=currently processed vertex) as LU-parent					
				if ((superC >= 0 && vertex[superC].parentLU == parentLU) || vertex[exR[xC]].parentLU == parentLU) {

					boolean foundPairing = false;
					// loop over and check currently joined failed pivots
					for (int xR = 0, xP = 0; xR < xPEnd; xR++, xP++) {
						if (exR[xR] == exPR[xP]) {						// found matching row index pairing for column index?
							exPC[xP] = exC[xC];							// complete the pairing by inserting the missing column index
							foundPairing = true; break;					// loop is done, signal outer loop to continue
						}
					}
					if (foundPairing) continue;
					// only allocate new failed pivots arrays if new duplets are found
					if (exPRn == null) exPRn = new int[xCEnd + xPEnd];	// worst case allocation
					if (exPCn == null) exPCn = new int[xCEnd + xPEnd];	// worst case allocation
					exPRn[xPn] = -1;									// new failed/extra pivot add:, we have no row index yet
					exPCn[xPn++] = exC[xC];								// but we have the column index

					xPSum++;											// increment resultant extra pivots count
				}
			}
		}
		
		// if we located new extra pivots and made a new extra pivots array, append the existing array
		// same case will apply if extra pivots array was null/empty
		if (xPSum > 0) {
			for (int xP = 0, xP2 = 0, xPn = xPSum; xP < xPEnd; xP++, xPn++, xP2++) {
				exPRn[xPn] = exPR[xP2]; exPCn[xPn] = exPC[xP2];
			}
			eJDP[0] = exPRn; eJDP[1] = exPCn;							// reassign existing extra pivots array to the new joined one
			helpArrayP[0] += xPSum;
			xPEnd += xPSum;
		}

		// since we don't expect the faulty situation of having UNMATCHED index pairs in the extras pivots array,
		// the loop for matching row extras (otherwise analogous to previous loop) will not have the unmatched row addition option
		
		if (eJDR[0] != null) {											// do we have any extras rows?
			int[] exR = eJDR[0], exC = eJDR[1], helpArrayR = eJDR[2];
			int xREnd = helpArrayR[0];

			for (int xR = 0; xR < xREnd; xR++) {
				int superC = vertex[exR[xR]].superChild;
				// if the vertex referred to by the column index has root (=currently processed vertex) as LU-parent					
				if ((superC >= 0 && vertex[superC].parentLU == parentLU) || vertex[exR[xR]].parentLU == parentLU) {

					// loop over and check currently joined failed pivots
					for (int xC = 0, xP = 1; xC < xPEnd; xC++, xP++) {
						if (exC[xC] == exPC[xP]) {						// found matching column index pairing for row index?
							exPR[xP] = exR[xR];							// complete the pairing by inserting the missing row index
							break;
						}
					}
				}
			}
		}
	}

	
	// method joins two sorted index arrays in sorted order, with array lengths & offsets into array ix1 and ix2 provided
	// by elements 0-3 of helpArray: [ix1 length, ix1 offset, ix2 length, ix2 offset]
	// the result arrays size will be summed in element 0 of helpArray
	// if helpArray[4]=0: just do simple sorted joining
	// if helpArray[4]=1: do joining with duplicates added (index duplicates count will be put in helpArray element 5)
	// if helpArray[4]=2: while joining, offsets of indexes from ix1 into the joined array ix3 will be stored in helpArray
	// if helpArray[4]=3: join duplets by unsorted add to destination duplets array
	// if helpArray[4]=4: join to new arrays by sorted add the column indexes and their pairings out of unsorted duplets
	// if helpArray[4]=5: join to new arrays by sorted add the row indexes and their pairings out of unsorted duplets
	// note: helpArray is assumed to be allocated properly = (6 elements) or (6 elements + element count of ix1)
	private static int FAILED_VECSIZE = 2; //FrontalMatrix.FAILED_VECSIZE;
	
	void joinIndexes(int[][] idxJgroup, int[][] ix2group) {
		
		int[] ix1 = idxJgroup[0], ix3 = null, ix3b = null, helpArray = idxJgroup[2], ix2 = ix2group[0];
		int l1 = helpArray[0], l2 = helpArray[2];

		// treat extra pivots as special case of unsorted addition of duplets: [global row, global column]
		if (helpArray[4] == 3) {
			// this is unsorted addition of duplets to a destination array of duplets
			int lExP = (l1 + l2) * FAILED_VECSIZE, i2 = 0;
			if (ix1 != null && lExP <= ix1.length)
				ix3 = ix1;					// we're adding to ix1, but if original array can fit all elements, add directly to it
			else ix3 = idxJgroup[0] = new int[trimToAllocBlock(lExP)];	// worst case allocation
			l1 *= FAILED_VECSIZE;
			for (int i = 0; i < l2; i++) {
				ix3[l1++] = ix2[i2++]; ix3[l1++] = ix2[i2++];
				i2++;													// skip the offset to index of failed pivot (used by add() method locally)
			}
			helpArray[0] += l2;
			return;
			
		} else if (helpArray[4] >= 4) {
			// this is sorted addition of elements out of unsorted duplets to a destination array: flag=4 stores row, flag=5 stores column
			// note: as asymmetric pivoting decouples pivot row & column indexes, the row or column pairing will also be stored
			ix3 = idxJgroup[0] = new int[trimToAllocExtras(l1 + l2)];	// worst case allocation
			ix3b = idxJgroup[1] = new int[trimToAllocExtras(l1 + l2)];
			
			int x3 = 0, l23 = l2 * FAILED_VECSIZE;
			if (helpArray[4] == 5)										// this is the row taken out of a failed pivot index pair
				for (int x2 = 0, x3b = 0; x2 < l23; x2 += FAILED_VECSIZE) {
					int i_ix2 = ix2[x2];
					if (i_ix2 > -1) {									// if some L-parent haven't already taken this index
						for (int x1 = 0; x1 < l1;)
							if (i_ix2 > ix1[x1]) ix3[x3++] = ix1[x1++];
						ix3[x3++] = i_ix2; //ix2[x2] = -1;			// flag this index as absorbed after moving it
						ix3b[x3b++] = ix2[x2 + 1];						// also store this failed pivot index's column pairing
					}
				}
			else														// this is the column taken out of a failed pivot index pair
				for (int x2 = 1, x3b = 0; x2 < l23; x2 += FAILED_VECSIZE) {
					int i_ix2 = ix2[x2];
					if (i_ix2 > -1) {									// if some L-parent haven't already taken this index
						for (int x1 = 0; x1 < l1;)						// copy over any indexes that were lower than the joined one
							if (i_ix2 > ix1[x1]) ix3[x3++] = ix1[x1++];
						ix3[x3++] = i_ix2; //ix2[x2] = -1;			// flag this index as absorbed after moving it
						ix3b[x3b++] = ix2[x2 - 1];						// also store this failed pivot index's row pairing
					}
				}
			// now we need to eliminate those row & column pairings that have been completely absorbed by an L-parent and an U-parent
			// the remaining full-pairings & half-pairings are waiting to be absorbed by either some L/U-parent or a LU-parent
			for (int x2s = 0, x2d = 0; x2s < l23;) {
				if (ix2[x2s] < 0 && ix2[x2s + 1] < 0) {				// if both row & column part have been absorbed, skip over duplet
					x2s += FAILED_VECSIZE; continue; }
				ix2[x2d++] = ix2[x2s++]; ix2[x2d++] = ix2[x2s++];
			}
			bubbleSort(ix3, ix3b, x3);									// as the array is not guaranteedly sorted, a do final bubble sort
			helpArray[0] = x3;
			return;
		}
		
		// following cases are regular index joining cases
		ix3 = idxJgroup[0] = new int[l1 + l2];							// worst case allocation
		ix3b = idxJgroup[1] = new int[l1 + l2];	

		int[] ix1b = idxJgroup[1], ix2b = ix2group[1];
		int x1 = helpArray[1], x2 = helpArray[3], x3 = 0;				// get the offsets into source index arrays
		l1 += x1; l2 += x2;
		if (helpArray[4] < 2) {											
			boolean joinDuplicates = helpArray[4] == 1 ? true : false;
			if (ix2b == null) {											// this block sort-copies indexes without pairings
				while (x1 < l1 && x2 < l2) {
					if (ix1[x1] < ix2[x2]) ix3[x3++] = ix1[x1++];
					else if (ix1[x1] > ix2[x2]) ix3[x3++] = ix2[x2++];
					else {
						ix3[x3++] = ix1[x1++];
						if (joinDuplicates) { ix3[x3++] = ix2[x2++]; helpArray[5]++; }
						else x2++;
					}
				}
				// we also end up here if either ex1 or ex2 is null, effectively just copying the arrays
				while (x1 < l1) ix3[x3++] = ix1[x1++];
				while (x2 < l2) ix3[x3++] = ix2[x2++];
				
			} else  {													// this block sort-copies indexes with pairings in the second array
				while (x1 < l1 && x2 < l2) {
					if (ix1[x1] < ix2[x2]) { ix3[x3] = ix1[x1]; ix3b[x3++] = ix1b[x1++]; }
					else if (ix1[x1] > ix2[x2]) { ix3[x3] = ix2[x2]; ix3b[x3++] = ix2b[x2++]; }
					else {
						ix3[x3] = ix1[x1];
						ix3b[x3++] = ix1b[x1++];
						if (joinDuplicates) { ix3[x3] = ix2[x2]; ix3b[x3++] = ix2b[x2++]; helpArray[5]++; }
						else x2++;
					}
				}
				// we also end up here if either ex1 or ex2 is null, effectively just copying the arrays
				while (x1 < l1) { ix3[x3] = ix1[x1]; ix3b[x3++] = ix1b[x1++]; }
				while (x2 < l2) { ix3[x3] = ix2[x2]; ix3b[x3++] = ix2b[x2++]; }
			}
			helpArray[0] = x3;
			
		} else {
			if (helpArray.length < HELPARRAY_SIZE + l2) {				// helpArray must at least fit it's own variables + length of ix2
				int[] helpArray2 = helpArray;
				helpArray = idxJgroup[2] = new int[HELPARRAY_SIZE + l2];
				helpArray[0] = helpArray2[0]; helpArray[1] = helpArray2[1]; helpArray[2] = helpArray2[2];
				helpArray[3] = helpArray2[3]; helpArray[4] = helpArray2[4];
			}
			
			int h = HELPARRAY_SIZE;
			if (ix2b == null) {											// this block copies an array without a twin pairings array
				while (x1 < l1 && x2 < l2) {
					if (ix1[x1] < ix2[x2]) ix3[x3++] = ix1[x1++];
					else if (ix1[x1] > ix2[x2]) {
						helpArray[h++] = x3;							// store current offset of index from ix2 into the joined array ix3
						ix3[x3++] = ix2[x2++];
					}
					else {
						helpArray[h++] = x3;							// store current offset of index from ix1 into the joined array ix3
						ix3[x3++] = ix1[x1++]; x2++;
					}
				}
				// we also end up here if either ex1 or ex2 is null, effectively just copying the array
				while (x1 < l1) ix3[x3++] = ix1[x1++];
				while (x2 < l2) { helpArray[h++] = x3; ix3[x3++] = ix2[x2++]; }
				
			} else {													// this block copies an array AND it's twin pairings array
				while (x1 < l1 && x2 < l2) {
					if (ix1[x1] < ix2[x2]) { ix3[x3] = ix1[x1]; ix3b[x3++] = ix1b[x1++]; }
					else if (ix1[x1] > ix2[x2]) {
						helpArray[h++] = x3;							// store current offset of index from ix2 into the joined array ix3
						ix3[x3] = ix2[x2]; ix3b[x3++] = ix2b[x2++];
					}
					else {
						helpArray[h++] = x3;							// store current offset of index from ix1 into the joined array ix3
						ix3[x3] = ix1[x1]; ix3b[x3++] = ix1b[x1++];
						x2++;
					}
				}
				// we also end up here if either ex1 or ex2 is null, effectively just copying the array
				while (x1 < l1) { ix3[x3] = ix1[x1]; ix3b[x3++] = ix1b[x1++]; }
				while (x2 < l2) { helpArray[h++] = x3; ix3[x3] = ix2[x2]; ix3b[x3++] = ix2b[x2++]; }
			}
			helpArray[0] = x3;
		}
	}
	
	// method bubble sorts failed row or column pivot index together with it's aspect pairing index in help array
	static void bubbleSort(int[] a1, int[] a2, int n) {
		n--;
		for (boolean sorted = false; !sorted;) {
			sorted = true;
			for (int i = 0, i1 = 1; i < n; i++, i1++)
				if (a1[i] > a1[i1]) {
					int temp = a1[i]; a1[i] = a1[i1]; a1[i1] = temp;
					temp = a2[i]; a2[i] = a2[i1]; a2[i1] = temp;
					sorted = false;
				}
		}		
	}

	
	
	// adds edge between vtx1 and vtx2, doing necessary reallocations of edge list
	// method is generic, doing fill-in of cut-edge holes and allowing testing for existing edge
	public boolean addEdge(int v1, int v2, int flags, boolean test) {
		
		FrontalDAGVertex vtx = vertex[v1];
		// if test=true, test if edge is already existing
		if (test)
			for (int e = 0, eEnd = vtx.edges + vtx.cutedges; e < eEnd; e++) 	
				if (vtx.edge(e) == v2) return false;
		
		// check if some previously cut edge can be reoccupied
		if (vtx.cutedges > 0) {
			vtx.edge[vtx.edge.length - vtx.cutedges--] = v2 | flags;
			vtx.edges++;
			edges++;
			return true;
		}
		
		// pushing edge onto array, check if array size needs expanding
		vtx.updateEdgeArray(1);		
		vtx.edge[vtx.edges++] = v2 | flags;
		edges++;
				
		int[] edgeP2;
		switch (flags) {
		case FrontalDAGVertex.L_CHILD:	vertex[v2].parentsL++; vertex[v2].parentsU++; break;
		case FrontalDAGVertex.U_CHILD:	vertex[v2].parentsU++; break;
		case FrontalDAGVertex.LU_CHILD:
			// readjust child's LU-ancestor if this was a LU-edge and this LU-ancestor is a lower/closer LU-ancestor
			// note: since there is only one LU-parent, it is indexed directly by internal variable parentLU
			if (vertex[v2].parentLU > v1) vertex[v2].parentLU = v1;
		}
		
		// add return L/U edge from v2 to v1, note: the L & U edges have been laid in order, first L-edges, then U-edges
		// this creates extra hassle of pushing entries +1 to insert an L-edge, besides the case of a necessary reallocation
		// parentsU[] counts have on edgeP initialisation been redefined as a L + U edge sum counts, while parentsL[] are counting L-edges
		int[] edgeP2new = edgeP2 = vertex[v2].edgeP;
		if (edgeP2 == null || edgeP2.length < parentsU[v2] + 1) {
			edgeP2new = vertex[v2].edgeP = new int[trimToAllocBlock(parentsU[v2] + 1)];
			// case of reinitialising the edgeP array of child vertex (reallocation was needed)
			switch (flags) {
			case FrontalDAGVertex.L_CHILD: {		// to insert L-edge, copy all L-edges, insert L-edge into the space, append all U-edges
				int i2 = 0, i3 = 0; 			for (int i2End = parentsL[v2]; i2 < i2End; i2++, i3++) edgeP2new[i3] = edgeP2[i2];
				edgeP2new[i3++] = v1 | flags;	for (int i2End = parentsU[v2]; i2 < i2End; i2++, i3++) edgeP2new[i3] = edgeP2[i2];
				parentsL[v2]++; parentsU[v2]++; break; }
			case FrontalDAGVertex.U_CHILD:			// to insert U-edge, copy L & U-edges, append U-edge on top
				int i2 = 0; for (int i2End = parentsU[v2]; i2 < i2End; i2++) edgeP2new[i2] = edgeP2[i2];	
				edgeP2new[i2] = v1 | flags;
				parentsU[v2]++;
			}
		} else {
			// case of reusing the edgeP array of child vertex (insertion into same array):
			switch (flags) {
			case FrontalDAGVertex.L_CHILD:		// for L-edge, shift U-edges one step up, insert L-edge in the space
				int i2End = parentsL[v2];
				for (int i3 = parentsU[v2], i2 = i3 - 1; i2 >= i2End; i2--, i3--) edgeP2new[i3] = edgeP2[i2];
				edgeP2new[i2End] = v1 | flags;
				parentsL[v2]++; parentsU[v2]++; break;
			case FrontalDAGVertex.U_CHILD:		// for U-edge, simply append it
				edgeP2new[parentsU[v2]] = v1 | flags; parentsU[v2]++; break;
			}
		}

		return true;
	}
	

	// edge cut method takes an optional edge index if it has been found beforehand
	public boolean cutEdge(int v1, int v2, int e) {
		
		FrontalDAGVertex vtx = vertex[v1];
		if (vtx.cutedges > vtx.edges) vtx.purgeFlagged();			// clear exsessive erased edges
		int e1 = 0, e1End = vtx.edges + vtx.cutedges;
		boolean edgeFound = false;
		
		if (e == -1) {												// edge index not supplied, search the edge
			for (; e1 < e1End; e1++)
				if (vtx.edge(e1) == v2) { edgeFound = true; break; }		// if sought edge was found	
			if (!edgeFound) return false;
		} else e1 = e;

		int flags = vtx.edge[e1] & FrontalDAGVertex.FLAGS;
		// save the cut-flagged edge index at back of the array
		if (--vtx.edges == 0) vtx.cutedges = 0;			// if last edge was cut, cut-flag count becomes irrelevant
		else {
			vtx.edge[vtx.edge.length - ++vtx.cutedges] = e1;
			vtx.edge[e1] = -1;							// -1 cutflags an edge, making impossible to match it with another vertex
		}
		edges--;

		FrontalDAGVertex vtx2 = vertex[v2];
		switch (flags) {
		case FrontalDAGVertex.L_CHILD:	vtx2.parentsL--; vtx2.parentsU--; break;
		case FrontalDAGVertex.U_CHILD:	vtx2.parentsU--; break;
		case FrontalDAGVertex.LU_CHILD:
			// take care of case of initialised parent edges and the case of LU-parent edge being cut
			if (vtx2.parentLU == v1) vtx2.parentLU = Integer.MAX_VALUE;
		}

		if (vtx2.edgeP != null) {									// parent edges can be uninitialised, execute this block only if they're initialised
			int[] edgeP = vtx2.edgeP;
			int e2 = 0;
			switch (vtx.edge[e1] & FrontalDAGVertex.FLAGS) {
			case FrontalDAGVertex.L_CHILD: 							// for L-child, find edge in the 0 to parentsL region
				for (int e2End = parentsL[v2]; e2 < e2End; e2++)
					if ((edgeP[e2]  & (0xFFFFFFFF - FrontalDAGVertex.FLAGS)) == v1) break;
				parentsL[v2]--;  break; 
			case FrontalDAGVertex.U_CHILD:							// for U-child,  find edge in the parentsL to parentsU region
				e2 = parentsL[v2];
				for (int e2End = parentsU[v2]; e2 < e2End; e2++)
					if ((edgeP[e2]  & (0xFFFFFFFF - FrontalDAGVertex.FLAGS)) == v1) break;
			}
			for (int e3 = e2 + 1, e3End = parentsU[v2]; e3 < e3End; e2++, e3++)		// shift all edges left over the cut parent edge
				edgeP[e2] = edgeP[e3];
				parentsU[v2]--;
		}

		return true;
	}
	

	
	// method makes union of edges between this DAG and DAG2, method expects vertex sets to be identical
	// the edges are added in a sorted fashion, but the sorting is NOT a proper sort, only a size-ordered additon
	// two properly sorted arrays WILL be joined in sorted fashion
	public FrontalDAG uniteEdgeSetLU(FrontalDAG dag2) {
		
		boolean order = distantFirst;
		
		for (FrontalDAGVertex vtx : vertex) {
			FrontalDAGVertex vtx2 = dag2.vertex[vtx.i];			// find equivalent vertex in DAG2
			int e1End = vtx.edges + vtx.cutedges, e2End = vtx2.edges + vtx2.cutedges;
			vtx.parentsU = vtx2.parentsU;
			vtx.parentsU = vtx.parentsL = 0;

			int[] edgeN = new int[trimToAllocBlock(vtx.edges + vtx2.edges)];
				
			int e1 = 0, e2 = 0, e1c = 0;
			while (e1 < e1End && e2 < e2End) {
				
				int edge1 = vtx.edge[e1], edge2 = vtx2.edge[e2];				// skip edge holes in both arrays
				if (edge1 == -1) { e1++; continue; }
				if (edge2 == -1) { e2++; continue; }
				
				int edgeD = vtx.edgeIndexS(edge2, e1, order);					// check if we're adding a duplicate edge
				if (edgeD >= 0) {												// if yes, we're assuming we're adding an equal U-edge
					if (vtx.is_S_edge(edgeD)) {									// unless it is supernodal edge chains, which we skip over
						edgeN[e1c++] = vtx.edge[e1++];
						continue;
					}
					edgeN[e1c++] = vtx.edge[e1++] | FrontalDAGVertex.LU_CHILD;	// so instead we convert the L-edge to an LU-edge
					edge2 &= (0xFFFFFFFF - FrontalDAGVertex.FLAGS);

					if (vertex[edge2].parentLU > vtx.i)							// set child's parent-LU index if a closer parent is found
						vertex[edge2].parentLU = vtx.i;
					e2++; continue;
				}
				
				// edges from DAG1 & DAG2 were not duplicates, now add them in proper order
				
				int edge1T = edge1, edge2T = edge2;
				edge1 &= (0xFFFFFFFF - FrontalDAGVertex.FLAGS);					// strip flags from vertex indexers
				edge2 &= (0xFFFFFFFF - FrontalDAGVertex.FLAGS);
				if (distantFirst) {
					if (edge1 < edge2) {										// edge from vertex1 was lower and is added first
						edgeN[e1c++] = vtx.edge[e1++];
						if ((edge1T & FrontalDAGVertex.L_CHILD) != 0)		vertex[edge1].parentsL++;
						else if ((edge1T & FrontalDAGVertex.U_CHILD) != 0)	vertex[edge1].parentsU++;
					} else  {													// edge from vertex2 was lower and is added first
						edgeN[e1c++] = vtx2.edge[e2++];
						if ((edge2T & FrontalDAGVertex.L_CHILD) != 0)		vertex[edge2].parentsL++;
						else if ((edge2T & FrontalDAGVertex.U_CHILD) != 0)	vertex[edge2].parentsU++;
						edges++;
					}
				} else {
					if (edge1 > edge2) {										// edge from vertex1 was higher and is added first
						edgeN[e1c++] = vtx.edge[e1++];
						if ((edge1T & FrontalDAGVertex.L_CHILD) != 0)		vertex[edge1].parentsL++;
						else if ((edge1T & FrontalDAGVertex.U_CHILD) != 0)	vertex[edge1].parentsU++;
					} else  {													// edge from vertex2 was higher and is added first
						edgeN[e1c++] = vtx2.edge[e2++];
						if ((edge2T & FrontalDAGVertex.L_CHILD) != 0)		vertex[edge2].parentsL++;
						else if ((edge2T & FrontalDAGVertex.U_CHILD) != 0)	vertex[edge2].parentsU++;
						edges++;
					}
				}
			}
			
			// if edges remained in whichever array, copy them over
			
			while (e1 < e1End) {											
				int edge1 = vtx.edge[e1++];
				if (edge1 != -1) {
					edgeN[e1c++] = edge1;
					if ((edge1 & FrontalDAGVertex.L_CHILD) != 0)		vertex[edge1 & 0x0FFFFFFF].parentsL++;
					else if ((edge1 & FrontalDAGVertex.U_CHILD) != 0)	vertex[edge1 & 0x0FFFFFFF].parentsU++;
				}
			}
			while (e2 < e2End) {
				int edge2 = vtx2.edge[e2++];
				if (edge2 != -1) {
					edgeN[e1c++] = edge2; edges++;
					if ((edge2 & FrontalDAGVertex.L_CHILD) != 0)		vertex[edge2 & 0x0FFFFFFF].parentsL++;
					else if ((edge2 & FrontalDAGVertex.U_CHILD) != 0)	vertex[edge2 & 0x0FFFFFFF].parentsU++;
				}
			}
			
			vtx.edges = e1c;
			vtx.cutedges = 0;
			vtx.edge = edgeN;
		}
		return this;
	}
	
	
		
	
	// method walks down the DAG from vertex v1 looking for a path to vertex v2 along edges of pathType
	// note: method uses a non-recursive stack-based depth-search strategy
	public int[] pathTo(int v1, int v2, int pathType, int startEdge, int[] vStack, boolean clearVisits) {
		
		// allocate for maximal stacking length, potentially tracing through every vertex top-down
		if (vStack == null) vStack = new int[verts]; 				
		int vCnt = 0, v = v1;								// vCnt keeps track of position within the stack
		FrontalDAGVertex vtx = null;
		if (clearVisits) clearVisits();
		vStack[0] = v;
		visit(v);
		vertex[v].checked = startEdge;						// caller decides what edge to start checking from, normally from edge 0

		// iterate looking for path to v2, until every vertex has been visited by every possible path
		while (vCnt >= 0) {
					
			vtx = vertex[vStack[vCnt]];
			int eCnt = vtx.edges + vtx.cutedges;
			boolean noAdvance = true;
			// iterate while a non-visited vertex isn't found and edges of current vertex aren't exhausted
			do {
				while (vtx.checked < eCnt && vtx.edge[vtx.checked] == -1)			// skip over delete-flagged edges
					vtx.checked++;
				if (vtx.checked < eCnt) {											// while current vertex has edges left to check
					int edge = vtx.edge[vtx.checked];
					// go by either an L/LU path or an U/LU path and pass unhindered through supernodes
					if ((edge & pathType) != 0 || (edge & FrontalDAGVertex.SUPERNODE) != 0) {
						v = edge & (0xFFFFFFFF - FrontalDAGVertex.FLAGS);										
						if (!visited(v) && v >= v2) {
							if (v == v2) return vStack;								// if destination vertex found, return path to it
							visit(v);
							vStack[++vCnt] = v;										// push non-visited child vertex onto stack
							vertex[v].checked = 0;									// first time visit of vertex, clear it's checked edge count
							noAdvance = false;										// we found a non-visited vertex, advance to it
						} else vtx.checked++;
					} else {
						vtx.checked++;
					}
				} else {															// current vertex's edges exhausted...
					vtx.checked++;													// DEBUG: just a way to flag a checked vertex
					vCnt--;															// backtrack to previous vertex
					break;															// get out to main loop
				}
			} while (noAdvance);
		}
		return null;
	}	
	
		
	
	// method traces paths through the DAG, putting edges on a stack, thus using a non-recursive depth-first-search strategy
	// if an earlier visited vertex is reached, the edge to that vertex can be cut
	// if that is the last edge of that vertex, the vertex effectively becomes a leaf of a branch
	// method can take a preallocated vertex stack (saves transient memory), or allocates a temporary one if vStack=null
	public void transitiveReduction(int v, int[] vStack, boolean clearVisits) {
		
		if (visited(v)) return;							// if this entire subDAG has been visited, we're done
		// allocate for maximal stacking length, potentially tracing through every vertex top-down
		if (vStack == null) vStack = new int[v + 1]; 				
		int vCnt = 0;									// vCnt keeps track of position within the stack
		FrontalDAGVertex vtx = null;
		if (clearVisits) clearVisits();
		vStack[0] = v;
		visit(v);
		vertex[v].checked = 0;							// clear checked edge index from previous reductions

		// iterate until every vertex has been visited by every possible trail, cutting secondary trails
		while (vCnt >= 0) {
					
			vtx = vertex[vStack[vCnt]];
			int eCnt = vtx.edges + vtx.cutedges, eArrayL = vtx.edge.length, pastNonVisited = vtx.checked;
			boolean noAdvance = true;
			// iterate while a non-visited vertex isn't found and edges of current vertex aren't exhausted
			do {
				while (vtx.checked < eCnt && vtx.edge[vtx.checked] == -1)			// skip over delete-flagged edges
					vtx.checked++;
				if (vtx.checked < eCnt) {											// while current vertex has edges left to check
					v = vtx.edge(vtx.checked);										
					if (visited(v)) {												// if it's been visited...
						edges--;
						
						if (vtx.edges + vtx.cutedges * 2 >= eArrayL) {				// is array full disallowing edge flagging?
							eCnt -= vtx.cutedges + 1;								// take the cut-flag count off the edge count
							vtx.edge[vtx.checked] = -1;								// flag the edge
							vtx.purgeFlagged();										// purge flagged
							vtx.checked = pastNonVisited;
							vtx.edges--;
						} 
						else if (--vtx.edges == 0) {
							vtx.checked++;											// DEBUG: just a way to flag a checked vertex
							vtx.cutedges = 0;										// if last edge was cut, hole count becomes irrelevant
							vCnt--;													// nothing more to investigate, step to previous vertex
							break;
						} else { 													// edge cutting happens here
							vtx.edge[eArrayL - ++vtx.cutedges] = vtx.checked;		// cut the edge to it
							vtx.edge[vtx.checked] = -1;
							
							if (vtx.cutedges > vtx.edges) {							// if flagged holes increase over the edges...
								eCnt -= vtx.cutedges;								// take the cut-flagged off the edge count
								vtx.purgeFlagged();									// purge flagged
								vtx.checked = pastNonVisited;
							} else
								vtx.checked++;
						}
					} else {
						pastNonVisited = ++vtx.checked;
						vStack[++vCnt] = v;											// push non-visited child vertex onto stack
						visit(v);
						vertex[v].checked = 0;										// first time visit if vertex, clear it's checked edge count
						noAdvance = false;											// we found a non-visited vertex, advance to it
					}
				} else {															// current vertex's edges exhausted...
					vtx.checked++;	// DEBUG: just a way to flag a checked vertex
					vCnt--;															// backtrack to previous vertex
					break;															// get out to main loop
				}
			} while (noAdvance);
		}
	}


	
	// method checks for existing direct edge from vtx1 to vtx2
	public boolean edgeTo(int v1, int v2) {
		FrontalDAGVertex vtxB = vertex[v1];
		for (int e = 0, eEnd = vtxB.edges + vtxB.cutedges; e < eEnd; e++)
			if (vtxB.edge(e) == v2) return true;
		return false;
	}

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			HELPER METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	public int superChild(int v) { return vertex[v].superPivot >= 0 ? vertex[v].superChild : v; }
	public int superParent(int v) { return vertex[v].superPivot >= 0 ? vertex[v].superPivot : v; }
	
	long[] visitBitArray(int verts) { return new long[bitSets = (verts >> 6) + 1]; }
	public void clearVisits() { visitBits = new long[bitSets]; }
	public boolean visited(int vtx) { return (visitBits[vtx >> 6] & (0x1L << (vtx & 63))) != 0; }
	public void visit(int vtx) { visitBits[vtx >> 6] |= (0x1L << (vtx & 63)); }

	public boolean isProcessed(FrontalDAGVertex vtx) { return vtx.fm != null; }
	private static int trimToAllocBlock(int v) { return (v & (0xFFFFFFFF - DAG_ALLOCBLOCK + 1)) + DAG_ALLOCBLOCK; }
	private int trimToAllocExtras(int v) { return (v & (0xFFFFFFFF - allocBlockExtras + 1)) + allocBlockExtras; }
	public static int trimToPowerOfTwo(int i) {
		if (i > 65536) i <<= 16;
		if (i < 65536) {
			if (i < 256) {
				if (i < 16) {
					if (i < 4) return 4; else return 8;
				} else { if (i < 64) return i < 32 ? 32 : 64; else return i < 128 ? 128 : 256; }
			} else {
				if (i < 4096) {
					if (i < 1024) return i < 512 ? 512 : 1024; else return i < 2048 ? 2048 : 4096;
				} else { if (i < 16384) return i < 8192 ? 8192 : 16384; else return i < 32768 ? 32768 : 65536; }
			}
		} else {
			if (i < 16777216) {
				if (i < 1048576) {
					if (i < 262144) return i < 131072 ? 131072 : 262144; else return i < 524288 ? 524288 : 1048576;
				} else { if (i < 64) return i < 32 ? 32 : 64; else return i < 128 ? 128 : 256; }
			} else {
				if (i < 268435456) {
					if (i < 67108864) return i < 33554432 ? 33554432 : 67108864; else return i < 134217728 ? 134217728 : 268435456;
				} else { if (i < 1073741824) return i < 536870912 ? 536870912 : 1073741824; else return i < (1<<31) ? (1<<31) : 0xFFFFFFFF; }
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	@Override
	public String toString() {		
		StringBuilder sb = new StringBuilder();
		FrontalMatrix fm;
		
		sb.append("Frontal DAG: " + name + ", edges: " + edges + ", subproblems: " + subProblems + "\n");
		for (int i = 0; i < verts; i++)
			if (vertex[i] != null) {
				if (visited(i)) sb.append("#"); else sb.append(" ");
				if (vertex[i].fm != null) {
					fm = vertex[i].fm;
					vertex[i].fm = null;
					sb.append("f" + vertex[i].toString());
					vertex[i].fm = fm;
				} else sb.append(" " + vertex[i].toString());
			} 
		return sb.toString();
	}


	public void toGraphViz(boolean showGraph, boolean base1, boolean imageView) {
		
		if (verts > 400) {
			System.out.println("FrontalDAG.toGraphViz(): 400 node limit struck.");
			return;
		}
		int idxBase = base1 ? 1 : 0;
				StringBuilder sb = new StringBuilder();
		sb.append("digraph G {\n     size=\"75,100\"\n     mclimit=10\n     rankdir=BT\n");
		
		for (int v = 0; v < verts; v++) {
			FrontalDAGVertex vtx = vertex[v];
			for (int e = 0, eEnd = vtx.edges + vtx.cutedges; e < eEnd; e++) {
				if (vtx.edge[e] != -1) {
					if (vtx.is_L_edge(e)) {
						sb.append("     edge [color=black];\n");
						sb.append("     " + (v + idxBase) + " [shape=ellipse, width=.2, height=.2];\n");
						sb.append("     " + (v + idxBase) + " -> " + (vtx.edge(e) + idxBase) + " [label=\"L\"];\n");										
					}
					else if (vtx.is_U_edge(e)) {
						sb.append("     edge [color=black];\n");
						sb.append("     " + (v + idxBase) + " [shape=ellipse, width=.2, height=.2];\n");
						sb.append("     " + (v + idxBase) + " -> " + (vtx.edge(e) + idxBase) + " [label=\"U\"];\n");										
					}
					else if (vtx.is_LU_edge(e)) {
						sb.append("     edge [color=red];\n");
						sb.append("     " + (v + idxBase) + " [shape=ellipse, width=.2, height=.2];\n");
						sb.append("     " + (v + idxBase) + " -> " + (vtx.edge(e) + idxBase) + " [label=\"LU\"];\n");										
					}
					else if (vtx.is_N_edge(e)) {
						sb.append("     edge [color=green];\n");
						sb.append("     " + (v + idxBase) + " [shape=ellipse, width=.2, height=.2];\n");
						sb.append("     " + (v + idxBase) + " -> " + (vtx.edge(e) + idxBase) + " [label=\"N\"];\n");										
					}
					else if (vtx.is_S_edge(e)) {
						sb.append("     edge [color=blue];\n");
						sb.append("     " + (v + idxBase) + " [shape=ellipse, width=.2, height=.2];\n");
						sb.append("     " + (v + idxBase) + " -> " + (vtx.edge(e) + idxBase) + " [label=\"S\"];\n");										
					}
				}
			}
		}
		sb.append("}\n");

		String gwizName = name + "_gviz.txt";
		File file = new File(gwizName);
		if (!file.exists()) {
			try {	file.createNewFile();
			} catch (IOException e) { e.printStackTrace(); }
		}
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(file));
			bw.write(sb.toString());
			bw.flush();
			bw.close();
		} catch (IOException e) { e.printStackTrace(); }
		
		if (!showGraph) return;
		try {
			ProcessBuilder GVizPB;
			if (imageView)
				GVizPB = new ProcessBuilder("C:\\Program Files\\graphviz\\bin\\dot.exe", "-Tpng", gwizName, "-o", name + ".png");
			else
				GVizPB = new ProcessBuilder("C:\\Program Files\\graphviz\\bin\\dot.exe", "-Tps", gwizName, "-o", name + ".ps");
			try { GVizPB.start().waitFor();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			Process viewP;
			if (imageView)
				viewP = new ProcessBuilder("C:\\Program Files (x86)\\FastStone Image Viewer\\FSViewer.exe", name + ".png").start();
			else
				viewP = new ProcessBuilder("C:\\Program Files\\Ghostgum\\gsview\\gsview64.exe", name + ".ps").start();
					
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
}
