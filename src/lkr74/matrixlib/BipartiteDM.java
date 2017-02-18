package lkr74.matrixlib;


// The class of the bipartite graph for Dulmage-Mendelsohn decomposition will find the maximum matching of a matrix
// facilitating the tranformation of matrix into block triangular form (BTF)
// the class utilises Hopcroft-Karp algorithm for the maximum matching finding, which is stored in "pairing" array
// the array consists of [ pairings-of-rows x m, supersink x 1, pairings-of-columns x n ],
// so one adds an offset into this array to access the pairings of the columns (although the pairings are symmetrical anyway)
// the pairing indexes themselves are ranged as they are in the matrix: (0 to m) and (0 to n)
public class BipartiteDM {

	String name = "DAG";
	int verts = 0, vertsR = 0, vertsC = 0;
	int bitSets = 0, edgesRC = 0, edgesCR = 0;
	BipartiteVertex[] vertex = null;
	// pairing consists of references to what row is paired to what column in the max.matching, the first vertsR elements
	// are column pairings of the rows, then the sink element, the last elements are the row pairings of the columns: [ R=C, S, C=R ]
	int[] pairing = null;
	long[] visitBits = null;
	boolean matched = false;

	// bipartite graph constructor
	BipartiteDM (String name, int vertsR, int vertsC) {
		this.name = name;
		this.vertsR = vertsR; this.vertsC = vertsC; 
		this.verts = vertsR + 1 + vertsC;
		vertex = new BipartiteVertex[verts];				// allocate left (rows) + right (columns) side of bipartite, plus a supersink vertex
		for (int i = 0; i < verts; i++) {
			if (i == vertsR) vertex[i] = new BipartiteVertex(i, i | BipartiteVertex.SUPERSINK);
			else vertex[i] = new BipartiteVertex(i, i <= vertsR ? i : (i - vertsR - 1) | BipartiteVertex.C_VERTEX);
		}
		pairing = new int[verts];
		visitBits = new long[bitSets = (verts >> 6) + 1];
	}


	// constructor initialises the bipartite graph fom an NSP matrix and it's nz entries
	public BipartiteDM(NSPMatrix M) {

		// a bipartite with left side being the row vertexes and right side being the column vertexes, thus M + N vertexes
		this.name = "G(E,R,C)_" + M.name;
		this.vertsR = M.M; this.vertsC = M.N; 
		this.verts = vertsR + vertsC + 1;
		vertex = new BipartiteVertex[verts];				// allocate left (rows) + right (columns) side of bipartite, plus a supersink vertex
		for (int i = 0; i < verts; i++) {
			if (i == vertsR)	vertex[i] = new BipartiteVertex(i, i | BipartiteVertex.SUPERSINK);		// supersink case
			else				vertex[i] = new BipartiteVertex(i, i <= vertsR ? i : (i - vertsR - 1) | BipartiteVertex.C_VERTEX);
		}
		pairing = new int[verts];
		visitBits = new long[bitSets = (verts >> 6) + 1];	// visit-flag bits are not used by Copcroft-Karp which uses another tracking mechanism

		int partC = vertsR + 1;
		
		// build up the bipartite edges from existing nz values
		// DEBUG: two iteration options: only add edges from current (not including it) PIVOT, or add ALL nz edges within current row or column
		for (int i = 0; i < M.M; i++) {
			NspArray aHsp = M.Hsp[i];
			NspNode[] bHsp = aHsp.array;
			// opt1: look left-in-row from pivot, adding found nz as R->C edges
			//for (int i1 = M.pivotNsp[i].offH + 1; i1 < aHsp.nodes; i++)
			// opt2: add all R->C edges that exist in current row of NSP matrix
			for (int j1 = 0; j1 < aHsp.nodes; j1++)
				addEdge(i, bHsp[j1].c + partC, false);
				// opt3: do not add edge to oneself
				//if (i != bHsp[j1].c) addEdge(i, bHsp[j1].c + partC, false);
		}
				
		for (int j = 0; j < M.N; j++) {
			NspArray aVsp = M.Vsp[j];
			NspNode[] bVsp = aVsp.array;
			// opt1: look down-in-column from pivot, adding found nz as C->R edges
			//for (int i1 = M.pivotNsp[j].offV + 1; i1 < aVsp.nodes; i1++)
			// opt2: add all C->R edges that exist in current column of NSP matrix
			for (int i1 = 0; i1 < aVsp.nodes; i1++)
				// add them within the index range of the R->C vertexes, which begin at index vertsR
				addEdge(j + partC, bVsp[i1].r, false);
				// opt3: do not add edge to oneself
				//if (j != bVsp[i1].r) addEdge(j + partC, bVsp[i1].r, false);
		}
		
		// add edges from all column-partition vertexes to supersink vertex, which lies as last vertex of row partition
		for (int v = vertsR + 1, vEnd = v + vertsC; v < vEnd; v++)
			addEdge(v, vertsR, false);
		
		// the bipartite graph is ready for matching and transformations
	}
	
	
	@Override
	public BipartiteDM clone() {
		// first clone the static data of this DAG
		Object O = null;
		try { O = super.clone(); } catch (CloneNotSupportedException e) { e.printStackTrace(); }
		BipartiteDM bDM = (BipartiteDM) O;
		bDM.vertex = new BipartiteVertex[verts];
		for (int i = 0; i < verts; i++) bDM.vertex[i] = vertex[i].clone();
		bDM.pairing = pairing.clone();
		bDM.visitBits = visitBits.clone();
		return bDM;
	}

	
	
	// method finds the maximum matching for a supplied bipartite graph over a matrix, utilising Hopcroft-Karp algorithm
	public int maximumMatchingHK() {

		int maxAspect = vertsR > vertsC ? vertsR : vertsC;
		int[] vStack = new int[maxAspect * 2];
	
		for (int i = 0; i <= vertsR; i++) pairing[i] = -1; 				// flag all vertexes in pairings R=C as free vertexes
		for (int i = vertsR + 1; i < verts; i++) pairing[i] = vertsR; 	// set all vertexes in pairings C=R to supersink vertex

		// queue will keep vertex indexes for breadth first search and is allocated for worst case scenario, which is
		// that every vertex will point to adjacent vertexes of a nz FULL density row/column, but we also need cheap allocation
		// so with a buffer of buffers one has possibility to dynamically allocate consecutive buffers to keep ALL adjacencies
		int[][] queue = new int[maxAspect][];
		queue[0] = new int[vertsR];					// prepare for worst case element count (=rows of the matrix) of first layer of BFS
		int[] queueCount = new int[maxAspect];
		int matching = 0;
		
		while (breadthFirstSearchHK(queue, queueCount)) {	// BFS finds shortest paths from all free vertexes to supersink + marks out distances along them
			for (int r = 0; r < vertsR; r++)
				if (pairing[r] == -1 && depthFirstSearchHK(r, vStack, true)) matching++;
		}
		
		// erase the supersink edges, we don't need them anymore (eliminates unnecessary computation)
		for (int i = verts - 1; i > vertsR; i--) vertex[i].edges--;
		
		matched = true;													// allow toString() of the matched chains		
		if (Matrix.DEBUG_LEVEL > 1) {
			System.out.println("Hopcroft-Karp maximum matching:\n");
			System.out.println(this.toString());
		}
		return matching;
	}
	
	
	
	// method partitions into layers the bipartite graph, starting from the remaining free vertices at current max.matching iteration
	// the traversed edges must alternate between matched & unmatched (augmentation path)
	// the BFS search terminates at first layer where one or more free C-vertices are reached, these vertices are placed into a set
	// this set is defined by setting equal distances at these terminal C-vertices
	public boolean breadthFirstSearchHK(int[][] queue, int[] queueCount) {
		
		int layerBFS = 0, qCnt = 0;
		BipartiteVertex superSink = vertex[vertsR];
		
		for (int r = 0; r < vertsR; r++)									// collect the yet free vertexes (that are not in a matching)
			if (pairing[r] == -1) {											// is this a free vertex?
				vertex[r].distance = 0;
				queue[layerBFS][qCnt++] = r;								// enqueue it in current layer
			} else vertex[r].distance = Integer.MAX_VALUE;
		superSink.distance = Integer.MAX_VALUE;								// set supersink vertex to "infinity"
		queueCount[layerBFS] = qCnt - 1;
		
		while (layerBFS >= 0) {												// while layers of element queue of BFS are not exhausted

			if (--queueCount[layerBFS] < 0) { layerBFS--; continue; }		// if last element in layer was dequeued, retract to previous layer

			int r = queue[layerBFS][queueCount[layerBFS]];					// dequeue an element
			if (vertex[r].distance < superSink.distance) {
				
				BipartiteVertex vtx = vertex[r];
				int eEnd = vtx.edges;
				// allocate next BFS layer for worst case = no. of edges from previous layer's vertex
				if (queue[++layerBFS] == null || queue[layerBFS].length < eEnd)
					queue[layerBFS] = new int[eEnd];
				queueCount[layerBFS] = 0;									// we start with zero elements in new layer
				
				for (int e = 0; e < eEnd; e++) {							// process next layer in BFS shortest augm.path search
					int cPair = pairing[vtx.edge(e)];
					if (vertex[cPair].distance == Integer.MAX_VALUE) {
						vertex[cPair].distance = vtx.distance + 1;
						queue[layerBFS][queueCount[layerBFS]++] = cPair;	// enqueue the found vertex
					}
				}
			}
		}
		
		// since supersink vertex is last BFS layer, if it's set to "infinity", we're done
		return (superSink.distance != Integer.MAX_VALUE);	
	}
	
	
	
	
	// method does depth first search for either the supersink or for the currently shortest augmenting path
	// method is used by Hopcroft-Karp algorithm
	public boolean depthFirstSearchHK(int v1, int[] vStack, boolean clearVisits) {
		
		// allocate for maximal stacking length, potentially tracing through every vertex top-down
		int vCnt = 0, v = v1, partC = vertsR + 1;									// vCnt keeps track of position within the stack
		BipartiteVertex vtx = null;
		boolean foundSupersink = false;
		vStack[0] = v;
		vertex[v].checked = 0;								// start checking from edge 0

		// iterate looking for path to v2, until every vertex has been visited by every possible path
		while (vCnt >= 0 && !foundSupersink) {
					
			// we're doing the depth first search through vertex PAIRINGS, thus always from same R->C edge aspect
			vtx = vertex[vStack[vCnt]];
			int eCnt = vtx.edges;
			boolean noAdvance = true;
			// iterate while a non-visited vertex isn't found and edges of current vertex aren't exhausted
			do {
				if (!foundSupersink && vtx.checked < eCnt) {						// while supersink not found and current vertex has edges left to check
					v = vtx.edge(vtx.checked++);									// strip edge flags	
					int vPair = pairing[v];
					
					// if vertex belongs to next length for augmentation paths
					if (vertex[vPair].distance == vtx.distance + 1) {				
						if (vPair == vertsR) {										// did we reach supersink?
							foundSupersink = true;									// yes, flag to quit loop and do backtracking with adding of pairings
							vStack[++vCnt] = v;										// we need this R-vertex's paired supersink vertex
							break;
						}
						vStack[++vCnt] = v;											// push C-vertex onto stack, needed for backtrack pairings creation
						vStack[++vCnt] = vPair;										// push the R-pairing of C-vertex, we'll scan it's adjacents next
						vertex[vPair].checked = 0;									// first time visit of pairing, clear it's checked edge count
						noAdvance = false;											// we found next vertex in the augm.path, advance to it
					}
					
				} else {															// current vertex's edges exhausted w.out finding shortest augm.path
					vtx.distance = Integer.MAX_VALUE;								// set this free vertex's distance to "infinity"
					vCnt -=2 ;														// backtrack to previous vertex and it's pairing
					break;															// get out to main loop
				}
			} while (noAdvance);
		}
		
		if (foundSupersink) {
			while (vCnt >= 0) {
				int u = vStack[vCnt - 1];
				pairing[vStack[vCnt]] = u;											// do the pairing: C = R
				pairing[u] = vStack[vCnt] - partC;									// do the pairing: R = C
				vCnt -= 2;															// backtrack the found path to beginning
			}
			return true;
		}
		return false;
	}	
	
	
	// method rehashes row & column indexes according to found Hr/Hc, Sr/Sc, Vr/Vc matrix division of rows & columns
	// method uses the maximum matching of this bipartite graph, so a matching algorithm must be called beforehand
	public int[][] findBTF() {
		if (!matched) return null;
		
		int[][] idxBTF = new int[2][];
		idxBTF[0] = new int[vertsR];											// contains Hr, Sr, Vr
		idxBTF[1] = new int[vertsC];											// contains Hc, Sc, Vc
		int[] U = new int[vertsR > vertsC ? vertsR : vertsC];					// unmatched vertexes stored here
		boolean[] flagR = new boolean[vertsR], flagC = new boolean[vertsC];		// the vertexes that were moved to Hr/Hc or Vr/Vc are flagged off here
		int partC = vertsR + 1, iU = 0, iEndU = 0;								// iU & iEndU define boundaries for currently processed unmatched vertexes
		
		for (int c = partC, cEnd = c + vertsC; c < cEnd; c++)					// all unmatched C-vertexes -> U
			if (pairing[c] == vertsR) {
				U[iEndU++] = c - partC;
				flagC[c - partC] = true; }
		
		// note: Hr = rows reachable via altern.path from some unmatched column
		// note: Hc = columns reachable by altern.path from some unmatched column
		int[] Hc = idxBTF[1], Hr = idxBTF[0];
		int iHr = 0, iHc = 0;

		while (iU != iEndU) {
			int c = Hc[iHc++] = U[iU++];
			int[] edge = vertex[c + partC].edge;
			for (int e = 0, eEnd = vertex[c + partC].edges; e < eEnd; e++) {
				int r = edge[e] & (0xFFFFFFFF - BipartiteVertex.FLAGS);
				if (!flagR[r]) {												// if r of adj(c) is not in Hr (if it's not used up, basically)
					Hr[iHr++] = r;												// insert r into Hr[]
					flagR[r] = true;
					int cMate = pairing[r];
					if (cMate == -1) { continue; }								// unmatched vertex case, should flag an augmenting path
					if (!flagC[cMate]) {										// if c pairing of r is not in U or Hc (if it's not used up, basically)
						U[iEndU++] = cMate;										// then add it to U
						flagC[cMate] = true;
					}
				}
			}
		}
				
		iU = 0; iEndU = 0;
		for (int r = 0; r < vertsR; r++)										// all unmatched R-vertexes -> U
			if (pairing[r] == -1) {
				U[iEndU++] = r;
				flagR[r] = true; }
		
		// note: Vr = rows reachable via altern.path from some unmatched row
		// note: Vc = columns reachable by altern.path from some unmatched row
		flagR = new boolean[vertsR]; flagC = new boolean[vertsC];
		int[] Vc = idxBTF[1], Vr = idxBTF[0];									// give clarifying names, but we continue putting into same arrays
		int iVr = iHr, iVc = iHc;												// iHr & iHc mark end of H-part of BTF, iVr & iVc continue onward
		iVr = 0; iVc = 0;	
		
		while (iU != iEndU) {
			int r = Vr[iVr++] = U[iU++];
			int[] edge = vertex[r].edge;
			for (int e = 0, eEnd = vertex[r].edges; e < eEnd; e++) {
				int c = (edge[e] & (0xFFFFFFFF - BipartiteVertex.FLAGS)) - partC;
				if (!flagC[c]) {												// if c of adj(r) is not in Vc (if it's not used up, basically)
					Vc[iVc++] = c;												// insert c into Vc[]
					flagC[c] = true;
					int rMate = pairing[c + partC];
					if (rMate == vertsR) { continue; }							// unmatched vertex case, should flag an augmenting path
					if (!flagR[rMate]) {										// if r pairing of c is not in U or Vr (if it's not used up, basically)
						U[iEndU++] = rMate;										// then add it to U
						flagR[rMate] = true;
					}
				}
			}
		}
		
		return idxBTF;
	}
	
	
	
	// adds edge between vtx1 and vtx2, doing necessary reallocations of edge list
	// method is generic, doing fill-in of cut-edge holes and allowing testing for existing edge
	public boolean addEdge(int v1, int v2, boolean test) {
		
		BipartiteVertex vtx = vertex[v1];
		// if test=true, test if edge is already existing
		if (test) for (int e = 0, eEnd = vtx.edges + vtx.cutedges; e < eEnd; e++) if (vtx.edge(e) == v2) return false;
		
		// check if edge is to the supersink, set flag accordingly
		boolean superSink = false;
		if ((vertex[v2].idx & BipartiteVertex.SUPERSINK) != 0) { v2 |= BipartiteVertex.SUPERSINK; superSink = true; }
		
		// check if some previously cut edge can be reoccupied
		if (vtx.cutedges > 0) {
			vtx.edge[vtx.edge.length - vtx.cutedges--] = v2;
			vtx.edges++;
			if (superSink) return true;							// don't count edges to supersink
			if (v1 >= vertsR) edgesRC++; else edgesCR++;		// increment R->C or C->R edge count, depending on what partition v1 belongs to
			return true;
		}
		
		// pushing edge onto array, check if array size needs expanding
		vtx.updateEdgeArray(1);		
		vtx.edge[vtx.edges++] = v2;
		if (superSink) return true;								// don't count edges to supersink
		if (v1 >= vertsR) edgesRC++; else edgesCR++;
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			HELPER/INLINE METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	public void clearVisits() { visitBits = new long[bitSets]; }
	public boolean visited(int vtx) { return (visitBits[vtx >> 6] & (0x1L << (vtx & 63))) != 0; }
	public void visit(int vtx) { visitBits[vtx >> 6] |= (0x1L << (vtx & 63)); }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	@Override
	public String toString() {		
		StringBuilder sb = new StringBuilder();
		boolean typeGap = false;
		
		sb.append("Bipartite Graph: " + name + ", R->C: " + edgesRC + " C->R: " + edgesCR + "\n\n===== R->C: =====\n");
		for (int i = 0; i < verts; i++)
			if (vertex[i] != null) {

				if (i > vertsR && !typeGap) { sb.append("\n===== C->R: =====\n"); typeGap = true; }

				if (visited(i)) sb.append("#"); else sb.append(" ");
				sb.append(vertex[i].toString());
				if (i < vertsR) {
					if (pairing[i] != -1)
						sb.append("     p:(R" + i + " = C" + pairing[i] + ")\n");
				} else {
					if (i == vertsR) continue;			// skip supersink vertex pairing writeout
					if (pairing[i] != vertsR)
						sb.append("     p:(C" + (i - vertsR - 1) + " = R" + pairing[i] + ")\n");
				}
			} 
		
		sb.append("\n");
		
		if (!matched) return sb.toString();
		
		// sort the index chains by length
		int[] counts = new int[vertsR], order = new int[vertsR];
		for (int i = 0; i < vertsR; i++) {
			order[i] = i;
			if (pairing[i] != -1)
				for (int vR = i, vC = pairing[vR]; ; vR = vC, vC = pairing[vR]) {
					if (vC == i || vC == -1) break;
					counts[i]++;
				}
		}
		
		int elems = vertsR - 1;
		for (boolean sorted = false; !sorted;) {
			sorted = true;
			for (int i = 0, i1 = 1; i < elems; i++, i1++)
				if (counts[i] > counts[i1]) {
					int temp = counts[i]; counts[i] = counts[i1]; counts[i1] = temp;
					temp = order[i]; order[i] = order[i1]; order[i1] = temp;
					sorted = false;
				}
		}		

		// print the index chains in length-sorted order
		for (int i = 0; i < vertsR; i++) {
			int o = order[i];
			sb.append("Chain of R" + o + " (x" + (counts[i] + 1) + "): ");
			if (pairing[o] != -1)
				for (int vR = o, vC = pairing[vR]; ; vR = vC, vC = pairing[vR]) {
					if (vC == o || vC == -1) { sb.append(vR); break; }
					else sb.append(vR + "->");
				}
			sb.append("\n");
		}

		return sb.toString();
	}

}
