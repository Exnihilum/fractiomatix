package lkr74.fem1;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class FEM1Octree {
	
	// the class helps localising incoming FEM nodes and do a neighbourhood search on closest node 

	double xM, yM, zM, xP, yP, zP;			// boundaries of this octree octant (xM,yM,zM = x-,y-,z- and xP,yP,zP = x+,y+,z+)
	double xC, yC, zC;						// subdivision centroid of this octree octant
	short status = 0, level = 0;			// status tells what octant it is, level tells how high up we are in the octree

	public int nodes, nodesX, branches = 0;	// nodesX = extraneous nodes past the true nodes, inserted/deleted for in-tree processing (extra space allocated for those)
	int[] node;								// node references of this octant container
	public int edges;
	int[] edge;								// edge references of this octant container
	
	FEM1Octree parent =  null;				// backlink to parent
	public FEM1Octree[] octant = null;		// array of 8 child octant references
	int octantIterator = 0;					// bitwise "array" holding the created octant indexes (spees up iteration)
	int octantIteratorMN = 0;				// bitwise "array" holding the created octant indexes that exceed maxNodes criterion (speeds up iteration)
	
	static float disbalanceCap = 0.85f;		// the max. level of disbalance allowed before stopping further branching
	static int vOBJcnt = 0, pOBJidx = 1;	// running vertex count & poly.index during OBJ export of octree visualisation
	// codes for the 8 subdivisions of an octree node (M = minus, P = plus, ex: OCT_PMP = (x+,y-,z+))
	final static byte OCT_MMM=0, OCT_PMM=1, OCT_MPM=2, OCT_PPM=3, OCT_MMP=4, OCT_PMP=5, OCT_MPP=6, OCT_PPP=7;
	public final static short MAX_LEVEL = 128; //Short.MAX_VALUE;
	public final static short FACTOR_NODESX = 128; //Short.MAX_VALUE;
	
	// instantiates first octant, finding bounding box and setting up node reference array
	public FEM1Octree(FEM1 fem) {
		
		if (fem.bBox == null) fem.bBox = new double[6];	
		double[] nodeCoord = fem.node;
		node = new int[extra_space(fem.nodes)];
		int nodes3 = fem.nodes * 3;
		xM = nodeCoord[0]; yM = nodeCoord[1]; zM = nodeCoord[2];					// arbitrarily choose first incoming vertex for initial bounding box
		xP = xM; yP = yM; zP = zM;

		for (int n = 3; n < nodes3;) {
			if (nodeCoord[n] < xM) xM = nodeCoord[n]; else if (nodeCoord[n] > xP) xP = nodeCoord[n]; n++;
			if (nodeCoord[n] < yM) yM = nodeCoord[n]; else if (nodeCoord[n] > yP) yP = nodeCoord[n]; n++;
			if (nodeCoord[n] < zM) zM = nodeCoord[n]; else if (nodeCoord[n] > zP) zP = nodeCoord[n]; n++;
		}
		double xDlt = xP - xM, yDlt = yP - yM, zDlt = zP - zM, maxBB = xDlt > yDlt ? xDlt : yDlt;
		if (maxBB < zDlt) maxBB = zDlt;
		double xP2 = (xP + xM + maxBB) * .5; xM = fem.bBox[0] = (xP + xM - maxBB) * .5; xP = fem.bBox[3] = xP2;	// equalise bounding box sides to the one maximal side
		double yP2 = (yP + yM + maxBB) * .5; yM = fem.bBox[1] = (yP + yM - maxBB) * .5; yP = fem.bBox[4] = yP2;
		double zP2 = (zP + zM + maxBB) * .5; zM = fem.bBox[2] = (zP + zM - maxBB) * .5; zP = fem.bBox[5] = zP2;
		xC = (xP + xM) * .5; yC = (yP + yM) * .5; zC = (zP + zM) * .5;
		
		int n = 0, nEnd = fem.nodes & 0xFFFFFFFC;
		while (n < nEnd) {						// set up initial node reference array
			node[n] = n++; node[n] = n++; node[n] = n++; node[n] = n++; }			// note: enforced simple loop unrolling
		for (int n2 = fem.nodes & 3; n2 > 0; n2--) node[nEnd++] = n++;
		nodes = fem.nodes;
		branches = 1;
		
	}
	
	// method instantiates a skeleton child
	public FEM1Octree(FEM1Octree parent, short status, int nodes, int[] node) {
		this.parent = parent; this.status = status;
		this.level = (short)(parent.level + 1);
		this.nodes = nodes; this.node = node;
		switch (status) {
			case OCT_MMM: xM=parent.xM; yM=parent.yM; zM=parent.zM; xP=parent.xC; yP=parent.yC; zP=parent.zC; break;
			case OCT_PMM: xM=parent.xC; yM=parent.yM; zM=parent.zM; xP=parent.xP; yP=parent.yC; zP=parent.zC; break;
			case OCT_MPM: xM=parent.xM; yM=parent.yC; zM=parent.zM; xP=parent.xC; yP=parent.yP; zP=parent.zC; break;
			case OCT_PPM: xM=parent.xC; yM=parent.yC; zM=parent.zM; xP=parent.xP; yP=parent.yP; zP=parent.zC; break;
			case OCT_MMP: xM=parent.xM; yM=parent.yM; zM=parent.zC; xP=parent.xC; yP=parent.yC; zP=parent.zP; break;
			case OCT_PMP: xM=parent.xC; yM=parent.yM; zM=parent.zC; xP=parent.xP; yP=parent.yC; zP=parent.zP; break;
			case OCT_MPP: xM=parent.xM; yM=parent.yC; zM=parent.zC; xP=parent.xC; yP=parent.yP; zP=parent.zP; break;
			case OCT_PPP: xM=parent.xC; yM=parent.yC; zM=parent.zC; xP=parent.xP; yP=parent.yP; zP=parent.zP; break; }
		xC = (xP + xM) * .5; yC = (yP + yM) * .5; zC = (zP + zM) * .5;
	}

	public FEM1Octree(FEM1Octree parent, short status, int nodes, int[] node, int edges, int[] edge) {
		this.parent = parent; this.status = status;
		this.level = (short)(parent.level + 1);
		this.nodes = nodes; this.node = node;
		this.edges = edges; this.edge = edge;
		switch (status) {
			case OCT_MMM: xM=parent.xM; yM=parent.yM; zM=parent.zM; xP=parent.xC; yP=parent.yC; zP=parent.zC; break;
			case OCT_PMM: xM=parent.xC; yM=parent.yM; zM=parent.zM; xP=parent.xP; yP=parent.yC; zP=parent.zC; break;
			case OCT_MPM: xM=parent.xM; yM=parent.yC; zM=parent.zM; xP=parent.xC; yP=parent.yP; zP=parent.zC; break;
			case OCT_PPM: xM=parent.xC; yM=parent.yC; zM=parent.zM; xP=parent.xP; yP=parent.yP; zP=parent.zC; break;
			case OCT_MMP: xM=parent.xM; yM=parent.yM; zM=parent.zC; xP=parent.xC; yP=parent.yC; zP=parent.zP; break;
			case OCT_PMP: xM=parent.xC; yM=parent.yM; zM=parent.zC; xP=parent.xP; yP=parent.yC; zP=parent.zP; break;
			case OCT_MPP: xM=parent.xM; yM=parent.yC; zM=parent.zC; xP=parent.xC; yP=parent.yP; zP=parent.zP; break;
			case OCT_PPP: xM=parent.xC; yM=parent.yC; zM=parent.zC; xP=parent.xP; yP=parent.yP; zP=parent.zP; break; }
		xC = (xP + xM) * .5; yC = (yP + yM) * .5; zC = (zP + zM) * .5;
	}

	
	
	
	// splits current octree node into 8 partitions, divides up the nodes, returns rate of disbalance within the current branch
	public double split(FEM1 fem, int maxItems) {
		
		double[] nodeCoord = fem.node;
		// worst case allocation of octosplitted node arrays
		int nMMM=0, nPMM=nodes, nMPM=nPMM+nodes, nPPM=nMPM+nodes, nMMP=nPPM+nodes, nPMP=nMMP+nodes, nMPP=nPMP+nodes, nPPP=nMPP+nodes;
		int eMMM=0, ePMM=edges, eMPM=ePMM+edges, ePPM=eMPM+edges, eMMP=ePPM+edges, ePMP=eMMP+edges, eMPP=ePMP+edges, ePPP=eMPP+edges;
		int nodes2 = nodes, edges2 = edges;
		int[] nodeOct = nodes > 0 ? new int[nodes * 8] : null, edgeOct = edges > 0 ? new int[edges * 8] : null;
		
		for (int n = 0; n < nodes; n++) {											// distribute node references into the 8 octants
			int n3 = node[n] * 3;
			double xN = nodeCoord[n3++], yN = nodeCoord[n3++], zN = nodeCoord[n3];
			if (xN <= xC) {
				if (yN <= yC) {	if (zN <= zC)	{ nodeOct[nMMM++] = node[n]; } else { nodeOct[nMMP++] = node[n]; }}
				else {			if (zN <= zC)	{ nodeOct[nMPM++] = node[n]; } else { nodeOct[nMPP++] = node[n]; }}
			} else {
				if (yN <= yC) {	if (zN <= zC)	{ nodeOct[nPMM++] = node[n]; } else { nodeOct[nPMP++] = node[n]; }}
				else {			if (zN <= zC)	{ nodeOct[nPPM++] = node[n]; } else { nodeOct[nPPP++] = node[n]; }}
			}
		}
		for (int e = 0; e < edges; e++) {											// distribute edge references into the 8 octants
			int eN = edge[e] * 2;
			int[] pEdgeN = fem.edgeNode;
			int bits = edgeMembership(pEdgeN[eN++], pEdgeN[eN++], pEdgeN[eN++], pEdgeN[eN++], pEdgeN[eN++], pEdgeN[eN++], 0, 0xFF);
			if ((bits & 1) != 0) nodeOct[eMMM++] = edge[e]; bits >>= 1; if ((bits & 1) != 0) nodeOct[ePMM++] = edge[e]; bits >>= 1;
			if ((bits & 1) != 0) nodeOct[eMPM++] = edge[e]; bits >>= 1; if ((bits & 1) != 0) nodeOct[ePPM++] = edge[e]; bits >>= 1;
			if ((bits & 1) != 0) nodeOct[eMMP++] = edge[e]; bits >>= 1; if ((bits & 1) != 0) nodeOct[ePMP++] = edge[e]; bits >>= 1;
			if ((bits & 1) != 0) nodeOct[eMPP++] = edge[e]; bits >>= 1; if ((bits & 1) != 0) nodeOct[ePPP++] = edge[e]; bits >>= 1;		
		}

		octant = new FEM1Octree[8];
		int[] nodeChild = null, edgeChild = null;
		
		int nPMM2=0, nMPM2=0, nPPM2=0, nMMP2=0, nPMP2=0, nMPP2=0, nPPP2=0, ePMM2=0, eMPM2=0, ePPM2=0, eMMP2=0, ePMP2=0, eMPP2=0, ePPP2=0;
		if (nodes > 0) {
			nPMM2 = nPMM - nodes2; nMPM2 = nMPM - (nodes2+=nodes); nPPM2 = nPPM - (nodes2+=nodes);
			nMMP2 = nMMP - (nodes2+=nodes); nPMP2 = nPMP - (nodes2+=nodes); nMPP2 = nMPP - (nodes2+=nodes); nPPP2 = nPPP - (nodes2+=nodes); }
		if (edges > 0) {
			ePMM2 = ePMM - edges2; eMPM2 = eMPM - (edges2+=edges); ePPM2 = ePPM - (edges2+=edges);
			eMMP2 = eMMP - (edges2+=edges); ePMP2 = ePMP - (edges2+=edges); eMPP2 = eMPP - (edges2+=edges); ePPP2 = ePPP - (edges2+=edges); }
		int octants = 0, octantsMN = 0;
		
		if (nMMM > 0 || eMMM > 0) {
			nodeChild = nMMM > 0 ? new int[extra_space(nMMM)] : null;
			edgeChild = eMMM > 0 ? new int[extra_space(eMMM)] : null;
			for (int n = 0; n < nMMM; n++) nodeChild[n] = nodeOct[n];
			for (int e = 0; e < eMMM; e++) edgeChild[e] = edgeOct[e];
			octant[OCT_MMM] = new FEM1Octree(this, OCT_MMM, nMMM, nodeChild, eMMM, edgeChild);		// initialise octant MMM
			octantIterator = octantIterator | OCT_MMM; octants++;
			if (nMMM > maxItems || eMMM > maxItems) { octantIteratorMN = octantIteratorMN | OCT_MMM; octantsMN++; }
		}
		nodes2 = nodes; edges2 = edges;
		if (nPMM2 > 0 || ePMM2 > 0) {
			nodeChild = nPMM2 > 0 ? new int[extra_space(nPMM2)] : null;
			edgeChild = ePMM2 > 0 ? new int[extra_space(ePMM2)] : null;
			for (int n = nodes, n2 = 0; n < nPMM; n++, n2++) nodeChild[n2] = nodeOct[n];
			for (int e = edges, e2 = 0; e < ePMM; e++, e2++) edgeChild[e2] = edgeOct[e];
			octant[OCT_PMM] = new FEM1Octree(this, OCT_PMM, nPMM2, nodeChild, ePMM2, edgeChild);	// initialise octant PMM
			octantIterator = (octantIterator << 3) | OCT_PMM; octants++;
			if (nPMM2 > maxItems || ePMM2 > maxItems) {  octantIteratorMN = (octantIteratorMN << 3) | OCT_PMM; octantsMN++; }
		}
		nodes2 += nodes; edges2 += edges;
		if (nMPM2 > 0 || eMPM2 > 0) {
			nodeChild = nMPM2 > 0 ? new int[extra_space(nMPM2)] : null;
			edgeChild = eMPM2 > 0 ? new int[extra_space(eMPM2)] : null;
			for (int n = nodes2, n2 = 0; n < nMPM; n++, n2++) nodeChild[n2] = nodeOct[n];
			for (int e = edges2, e2 = 0; e < eMPM; e++, e2++) edgeChild[e2] = edgeOct[e];
			octant[OCT_MPM] = new FEM1Octree(this, OCT_MPM, nMPM2, nodeChild, eMPM2, edgeChild);		// initialise octant MPM
			octantIterator = (octantIterator << 3) | OCT_MPM; octants++;
			if (nMPM2 > maxItems || eMPM2 > maxItems) { octantIteratorMN = (octantIteratorMN << 3) | OCT_MPM; octantsMN++; }
		}
		nodes2 += nodes; edges2 += edges;
		if (nPPM2 > 0 || ePPM2 > 0) {
			nodeChild = nPPM2 > 0 ? new int[extra_space(nPPM2)] : null;
			edgeChild = ePPM2 > 0 ? new int[extra_space(ePPM2)] : null;
			for (int n = nodes2, n2 = 0; n < nPPM; n++, n2++) nodeChild[n2] = nodeOct[n];
			for (int e = edges2, e2 = 0; e < ePPM; e++, e2++) edgeChild[e2] = edgeOct[e];
			octant[OCT_PPM] = new FEM1Octree(this, OCT_PPM, nPPM2, nodeChild, ePPM2, edgeChild);		// initialise octant MMM
			octantIterator = (octantIterator << 3) | OCT_PPM; octants++;
			if (nPPM2 > maxItems || ePPM2 > maxItems) { octantIteratorMN = (octantIteratorMN << 3) | OCT_PPM; octantsMN++; }
		}
		nodes2 += nodes; edges2 += edges;
		if (nMMP2 > 0 || eMMP2 > 0) {
			nodeChild = nMMP2 > 0 ? new int[extra_space(nMMP2)] : null;
			edgeChild = eMMP2 > 0 ? new int[extra_space(eMMP2)] : null;
			for (int n = nodes2, n2 = 0; n < nMMP; n++, n2++) nodeChild[n2] = nodeOct[n];
			for (int e = edges2, e2 = 0; e < eMMP; e++, e2++) edgeChild[e2] = edgeOct[e];
			octant[OCT_MMP] = new FEM1Octree(this, OCT_MMP, nMMP2, nodeChild, eMMP2, edgeChild);		// initialise octant MMP
			octantIterator = (octantIterator << 3) | OCT_MMP; octants++;
			if (nMMP2 > maxItems || eMMP2 > maxItems) { octantIteratorMN = (octantIteratorMN << 3) | OCT_MMP; octantsMN++; }
		}
		nodes2 += nodes; edges2 += edges;
		if (nPMP2 > 0 || ePMP2 > 0) {
			nodeChild = nPMP2 > 0 ? new int[extra_space(nPMP2)] : null;
			edgeChild = ePMP2 > 0 ? new int[extra_space(ePMP2)] : null;
			for (int n = nodes2, n2 = 0; n < nPMP; n++, n2++) nodeChild[n2] = nodeOct[n];
			for (int e = edges2, e2 = 0; e < ePMP; e++, e2++) edgeChild[e2] = edgeOct[e];
			octant[OCT_PMP] = new FEM1Octree(this, OCT_PMP, nPMP2, nodeChild, ePMP2, edgeChild);		// initialise octant PMP
			octantIterator = (octantIterator << 3) | OCT_PMP; octants++;
			if (nPMP2 > maxItems || ePMP2 > maxItems) { octantIteratorMN = (octantIteratorMN << 3) | OCT_PMP; octantsMN++; }
		}
		nodes2 += nodes; edges2 += edges;
		if (nMPP2 > 0 || eMPP2 > 0) {
			nodeChild = nMPP2 > 0 ? new int[extra_space(nMPP2)] : null;
			edgeChild = eMPP2 > 0 ? new int[extra_space(eMPP2)] : null;
			for (int n = nodes2, n2 = 0; n < nMPP; n++, n2++) nodeChild[n2] = nodeOct[n];
			for (int e = edges2, e2 = 0; e < eMPP; e++, e2++) edgeChild[e2] = edgeOct[e];
			octant[OCT_MPP] = new FEM1Octree(this, OCT_MPP, nMPP2, nodeChild, eMPP2, edgeChild);		// initialise octant MPP
			octantIterator = (octantIterator << 3) | OCT_MPP; octants++;
			if (nMPP2 > maxItems || eMPP2 > maxItems) { octantIteratorMN = (octantIteratorMN << 3) | OCT_MPP; octantsMN++; }
		}
		nodes2 += nodes; edges2 += edges;
		if (nPPP2 > 0 || ePPP2 > 0) {
			nodeChild = nPPP2 > 0 ? new int[extra_space(nPPP2)] : null;
			edgeChild = ePPP2 > 0 ? new int[extra_space(ePPP2)] : null;
			for (int n = nodes2, n2 = 0; n < nPPP; n++, n2++) nodeChild[n2] = nodeOct[n];
			for (int e = edges2, e2 = 0; e < ePPP; e++, e2++) edgeChild[e2] = edgeOct[e];
			octant[OCT_PPP] = new FEM1Octree(this, OCT_PPP, nPPP2, nodeChild, ePPP2, edgeChild);		// initialise octant PPP
			octantIterator = (octantIterator << 3) | OCT_PPP; octants++;
			if (nPPP2 > maxItems || ePPP2 > maxItems) { octantIteratorMN = (octantIteratorMN << 3) | OCT_PPP; octantsMN++; }
		}
		octantIterator |= (octants << 24);
		octantIteratorMN |= (octantsMN << 24);
		node = null; edge = null;					// a branch doesnt need item references (but keeping nodes count of entire branch can be useful)
		// branch disbalance deducted from ALL item types clustered across octants
		double disbalance = octantNodeDisbalance(nodes, nMMM+eMMM, nPMM2+ePMM2, nMPM2+eMPM2, nPPM2+ePPM2, nMMP2+eMMP2, nPMP2+nPMP2, nMPP2+eMPP2, nPPP2+ePPP2);
		return disbalance;
	}
	
	
	
	// method recursively constructs octree with at most maxItems items per container and splitting limited by maximal branch disbalance
	public int buildOctree(FEM1 fem, int maxItems, double disbalance) {
		
		// if current container has over maxNodes nodes & branch levels are not overrun & terminate limiter hasn't been assigned
		if (level < MAX_LEVEL && disbalance < disbalanceCap) {
			double disbalanceB = split(fem, maxItems);							// split into octants (if possible), becoming a branch
			branches = octantIterator >> 24;
			int octantIt = octantIteratorMN, oEnd = octantIteratorMN >> 24;
			for (int o = 0; o < oEnd; o++, octantIt >>= 3) {					// we are iterating over octants that contain more than "maxItems" items
				FEM1Octree octant = this.octant[octantIt & 7];
				branches += octant.buildOctree(fem, maxItems, (double)((disbalance + disbalanceB) * 0.5));				
			}
		}
		return branches;
	}
	
	
	// method locates the best fitting octant container of the supplied coordinate, returning EITHER a branch OR a leaf octant
	// method will stop if level is reached, setting level to MAX_LEVEL is a way to "ignore" the level limit
	public FEM1Octree locateCoordinate(double x, double y, double z, short level) {

		FEM1Octree oct = this, oct2 = this;
		while (oct2.octant != null && oct2.level < level) {
			if (x <= oct2.xC) {	if (y <= oct2.yC) {	if (z <= oct2.zC)	{ oct2 = oct2.octant[OCT_MMM]; } else { oct2 = oct2.octant[OCT_MMP]; }}
								else {				if (z <= oct2.zC)	{ oct2 = oct2.octant[OCT_MPM]; } else { oct2 = oct2.octant[OCT_MPP]; }}
			} else {			if (y <= oct2.yC) {	if (z <= oct2.zC)	{ oct2 = oct2.octant[OCT_PMM]; } else { oct2 = oct2.octant[OCT_PMP]; }}
								else {				if (z <= oct2.zC)	{ oct2 = oct2.octant[OCT_PPM]; } else { oct2 = oct2.octant[OCT_PPP]; }}}
			if (oct2 == null) break;
			oct = oct2;
		}
		return oct;
	}
	
	
	// method locates the best fitting octant container of the supplied coordinate, then adds the coordinate as the supplied node index
	// if the found topmost container was a branch, method creates a suboctant, inserting node as the only member
	// method returns the octant that became final container of the inserted node
	// the nodeX flag indicates whether to insert a real or extraneous node
	public FEM1Octree addNode(FEM1 fem, double x, double y, double z, int n, int maxItems, boolean nodeX) {
		int topmost = -1;
		FEM1Octree oct = this, oct2 = this;
		while (oct2.octant != null && oct2.level < level) {
			if (x <= oct2.xC) {
				if (y <= oct2.yC) {	if (z <= oct2.zC)	{ oct2 = oct2.octant[topmost = OCT_MMM]; } else { oct2 = oct2.octant[topmost = OCT_MMP]; }}
				else {				if (z <= oct2.zC)	{ oct2 = oct2.octant[topmost = OCT_MPM]; } else { oct2 = oct2.octant[topmost = OCT_MPP]; }}
			} else {
				if (y <= oct2.yC) {	if (z <= oct2.zC)	{ oct2 = oct2.octant[topmost = OCT_PMM]; } else { oct2 = oct2.octant[topmost = OCT_PMP]; }}
				else {				if (z <= oct2.zC)	{ oct2 = oct2.octant[topmost = OCT_PPM]; } else { oct2 = oct2.octant[topmost = OCT_PPP]; }}}
			if (oct2 == null) break;
			oct = oct2;
		}
		if (oct2 == null) {															// if containing child octant within branch octant isn't initialised
			oct2 = oct.octant[topmost] = new FEM1Octree(oct, (short)topmost, 0, new int[extra_space(1)]);
			if (nodeX) oct2.nodesX++; else oct2.nodes++;	
			oct2.node[0] = n;
		} else {																	// if target leaf container octant was found (or this is octree was a leaf)
			oct2.octantAddNode(x, y, z, n, nodeX);
			if (oct2.nodes > maxItems || oct2.edges > maxItems)
				oct2.split(fem, maxItems);											// if node+edge count exceeds maxItems, split up the octant
		}
		return oct2;
	}
	
	// method directly inserts node into a known octant, whether as extraneous or as a real node
	public void octantAddNode(double x, double y, double z, int n, boolean nodeX) {
		int nodesTot = nodes + nodesX;
		if (nodesTot + 1 >= node.length) {														// need to expand octant's node array?
			int[] nodeNew = new int[FEM1Octree.extra_space(nodesTot + 1)];
			
			if (nodeX) {for (int n1 = 0; n1 < nodesTot; n1++) nodeNew[n1] = node[n1];			// if adding an extraneous node, copy everything straight-on
						nodeNew[nodesTot] = n; 													// add extraneous node at end
						nodesX++; }
			else {		for (int n1 = 0; n1 < nodes; n1++) nodeNew[n1] = node[n1];				// if adding real node, leave gap between nodes & nodesX
						for (int n1 = nodes, n2 = n1 + 1; n1 < nodesTot; n1++, n2++) nodeNew[n2] = node[n1];
						nodeNew[nodes++] = n; }													// add real node in gap
			node = nodeNew;
		} else {
			if (nodeX)	{	node[nodesTot] = n; nodesX++; } 									// regular append at end of array, in nodesX interval
			else {			node[nodesX] = node[nodes]; node[nodes++] = n; }					// move out one node from nodes interval, insert new node
		}
	}

	
	
	// method divides down an edge into the 8 octants, bitflagging only the visited octants, using recursive calls
	// octant bitflags (high to low): PPP,MPP,PMP,MMP,PPM,MPM,PMM,MMM
	int edgeMembership(double xa, double ya, double za, double xb, double yb, double zb, int dim, int bits) {
		
		switch (dim) {
		case 0:
			if (xa <= xC) {
				if (xC < xb) {																		// point a in lower x-octants, b in higher
						double xbaD = 1. / (xa - xb), f1 = (xC - xa) * xbaD, f2 = (xb - xC) * xbaD;
						double ys = f1 > f2 ? ya+(yb-ya)*f1 : yb+(ya-yb)*f2, zs = f1 > f2 ? za+(zb-za)*f1 : zb+(za-zb)*f2;
						return 	bits & (edgeMembership(xa, ya, za, xC, ys, zs, dim + 1, 0x55) | edgeMembership(xC, ys, zs, xb, yb, zb, dim + 1, 0xAA));
				} else	return 	bits & edgeMembership(xa, ya, za, xb, yb, zb, dim + 1, 0x55);		// edge completely in OCT_M** domain
			} else if (xb <= xC) {
				if (xC < xb) {																		// reverse: point b in lower x-octants, a in higher
						double xbaD = 1. / (xb - xa), f1 = (xC - xb) * xbaD, f2 = (xa - xC) * xbaD;
						double ys = f1 > f2 ? yb+(ya-yb)*f1 : ya+(yb-ya)*f2, zs = f1 > f2 ? zb+(za-zb)*f1 : za+(zb-za)*f2;
						return 	bits & (edgeMembership(xb, yb, zb, xC, ys, zs, dim + 1, 0x55) | edgeMembership(xC, ys, zs, xa, ya, za, dim + 1, 0xAA));
				} else	return 	bits & edgeMembership(xb, yb, zb, xa, ya, za, dim + 1, 0x55);		// edge completely in OCT_M** domain
			} else		return	bits & edgeMembership(xa, ya, za, xb, yb, zb, dim + 1, 0xAA);		// edge completely in OCT_P** domain	
		case 1:
			if (ya <= yC) {
				if (yC < yb) {																		// point a in lower y-octants, b in higher
						double ybaD = 1. / (ya - yb), f1 = (yC - ya) * ybaD, f2 = (yb - yC) * ybaD;
						double xs = f1 > f2 ? xa+(xb-xa)*f1 : xb+(xa-xb)*f2, zs = f1 > f2 ? za+(zb-za)*f1 : zb+(za-zb)*f2;
						return 	bits & (edgeMembership(xa, ya, za, xs, yC, zs, dim + 1, 0x33) | edgeMembership(xs, yC, zs, xb, yb, zb, dim + 1, 0xCC));
				} else	return 	bits & edgeMembership(xa, ya, za, xb, yb, zb, dim + 1, 0x33);		// edge completely in OCT_*M* domain
			} else if (yb <= yC) {
				if (yC < yb) {																		// reverse: point b in lower y-octants, a in higher
						double ybaD = 1. / (yb - ya), f1 = (yC - yb) * ybaD, f2 = (ya - yC) * ybaD;
						double xs = f1 > f2 ? xb+(xa-xb)*f1 : xa+(xb-xa)*f2, zs = f1 > f2 ? zb+(za-zb)*f1 : za+(zb-za)*f2;
						return 	bits & (edgeMembership(xb, yb, zb, xs, yC, zs, dim + 1, 0x33) | edgeMembership(xs, yC, zs, xa, ya, za, dim + 1, 0xCC));
				} else	return 	bits & edgeMembership(xb, yb, zb, xa, ya, za, dim + 1, 0x33);		// edge completely in OCT_*M* domain
			} else		return	bits & edgeMembership(xa, ya, za, xb, yb, zb, dim + 1, 0xCC);		// edge completely in OCT_*P* domain	
		case 2:
			if (za <= zC) {
				if (zC < zb) {	return 	bits;														// point a in lower z-octants, b in higher	
				} else	return 	bits & 0x0F;														// edge completely in OCT_**M domain
			} else if (zb <= zC) {
				if (zC < zb) {	return 	bits;														// reverse: point b in lower x-octants, a in higher	
				} else	return 	bits & 0x0F;														// edge completely in OCT_**M domain
			} else		return	bits & 0xF0;														// edge completely in OCT_**P domain	
		}
		return 0;
	}
	
	// method inserts an edge segment into octree by recursive splitting & membership testing
	// TODO: debug method
	public int addEdge(FEM1 fem, double xa, double ya, double za, double xb, double yb, double zb, int edgeIndex, int maxItems, double disbalance) {
		
		// max depth reached or disbalanced exceeded at leaf level (thus overriding maxEdges criterion)
		if (level >= MAX_LEVEL || (edges <= maxItems && octant == null)) {
			octantAddEdge(edgeIndex);
			if (nodes > maxItems || edges > maxItems) split(fem, maxItems);				// split this octant if maxItems exceeded
			return 0;
		}	
		int bits = edgeMembership(xa, ya, za, xb, yb, zb, 0, 0xFF);
		
		int[] eDisb = { 0,0,0,0,0,0,0,0 };
		for (int o = 0, bits2 = bits; o < 8; o++, bits2 >>= 1)							// precalculate expected branch disbalance with the added edges
			if ((bits2 & 1) != 0 && octant[o] != null)
				eDisb[o] = octant[o].nodes + octant[o].edges + 1;
		double disbalanceB = octantNodeDisbalance(edges, eDisb[0], eDisb[1], eDisb[2], eDisb[3], eDisb[4], eDisb[5], eDisb[6], eDisb[7]);

		for (int o = 0; o < 8; o++, bits >>= 1) {										// start adding edges
			if ((bits & 1) != 0) {
				if (octant[o] == null) {												// if containing octant doesn't exist
					if (disbalance < disbalanceCap) {									// if disbalanceCap allows dividing edge membership further
						octant[o] = new FEM1Octree(this, (short) o, 0, null, 1, null);	// create branch containing 1 edge, account for it in disbalance
						branches++;														// which adds one branch
						branches += octant[o].addEdge(fem, xa, ya, za, xb, yb, zb, edgeIndex, maxItems, (disbalanceB + disbalance) * 0.5);
						
					} else {															// case of disbalanceCap reached
						int[] edge1 = {edgeIndex};
						octant[o] = new FEM1Octree(this, (short) o, 0, null, 1, edge1);	// make leaf octant with one edge
						branches++;														// which adds one branch
					}
				} else	octant[o].addEdge(fem, xa, ya, za, xb, yb, zb, edgeIndex, maxItems, disbalance);	// we're not at leaf level yet, continue
			}
		}
		return branches;		
	}
	
	
	// method directly inserts edge into a known octant
	public void octantAddEdge(int e) {
		if (edges + 1 >= edge.length) {															// need to expand octant's edge array?
			int[] edgeNew = new int[FEM1Octree.extra_space(edges + 1)];		
			for (int n1 = 0; n1 < edges; n1++) edgeNew[n1] = edge[n1];							// copy everything straight-on
			edgeNew[edges++] = e; edge = edgeNew;												// add extraneous edge at end	
		} else edge[edges++] = e;
	}

	
	
	// method returns array of all leaves within supplied octant, or supplied octant as single member if it's a leaf itself
	// method utilises nonrecursive tree walk with expansion into a stack
	public FEM1Octree[] leafOctantArray() {
		
		if (octant == null) { FEM1Octree[] leafArray = { this }; return leafArray; }	// octant was a leaf itself
		
		FEM1Octree[] leafArray = new FEM1Octree[branches];					// worst-case-allocate for octant leaves
		FEM1Octree[] stackO = new FEM1Octree[MAX_LEVEL - level];			// worst-case-allocate for octant stack
		stackO[0] = this;
		int[] stack = new int[MAX_LEVEL - level];							// worst-case stack allocation
		stack[0] = octantIterator;
		int s = 0, leaves = 0;
		while (s >= 0) {
			int c = stack[s] & 0xFF000000;									// get remaining children of current octant on stack
			if (c > 0) {						
				FEM1Octree octantChild = stackO[s].octant[stack[s] & 7];	// get outstanding child
				stack[s] = ((stack[s]&0xFFFFFF)>>3) | (c-0x1000000);		// recompose octant stack entry: pop child, make next outstanding
				stack[++s] = octantChild.octantIterator;
				stackO[s] = octantChild;
			} else {
				if (stackO[s].octant == null)								// if it was a leaf on top of octant stack
					leafArray[leaves++] = stackO[s--];						// confirm latest leaf octant, backtrack to parent
				else s--;
			}
		}
		return leafArray;
	}


	// method returns array of all leaves close enough to a coordinate within supplied octant, or supplied octant as single member if it's a leaf itself
	// method utilises nonrecursive tree walk with expansion into a stack
	public FEM1Octree[] leafOctantArrayByDistance(double x, double y, double z, double dist) {
		
		if (octant == null) { FEM1Octree[] leafArray = { this }; return leafArray; }	// octant was a leaf itself
		
		FEM1Octree[] leafArray = new FEM1Octree[branches];					// worst-case-allocate for octant leaves
		FEM1Octree[] stackO = new FEM1Octree[MAX_LEVEL - level];			// worst-case-allocate for octant stack
		stackO[0] = this;
		int[] stack = new int[MAX_LEVEL - level];							// worst-case stack allocation
		stack[0] = octantIterator;
		int s = 0, leaves = 0;
		while (s >= 0) {
			int c = stack[s] & 0xFF000000;									// get remaining children of current octant on stack
			if (c > 0) {						
				FEM1Octree octantChild = stackO[s].octant[stack[s] & 7];	// get outstanding child
				stack[s] = ((stack[s]&0xFFFFFF)>>3) | (c-0x1000000);		// recompose octant stack entry: pop child, make next outstanding
				if (octantChild.octantDistance(x, y, z, false) < dist) {	// stack child octant only if it fulfills closeness criterion
					stack[++s] = octantChild.octantIterator;
					stackO[s] = octantChild;
				}
			} else {														// if no children remain/exist
				if (stackO[s].octant == null)								// if it was a leaf on top of octant stack
					leafArray[leaves++] = stackO[s]; 						// add to leaf collection										
				s--;														// backtrack
			}
		}
		return leafArray;
	}

	
	
	// method finds the closest node of node n within supplied octree by neighbourhood backtracking/collection heuristics:
	// 1) if there are other local nodes, first the closest node distance is found, then the scope is expanded to the ancestor octant that fully
	// encompasses a bounding box defined by that distance, the leaf octants that fulfill the distance criterion are collected from that ancestor,
	// and their local nodes are recursively compared against the shortest distance attained so far
	// 2) if octant only holds one node, backtracking to first node-holding parent is made, the leaf octant closest to node is found,
	// the closest node distance is found from it's local nodes, then one proceeds as in point 1)
	// octants[0] can supply the node's container, if it's null then the node's container is found first
	// if distance[] supplied, closest distance will be returned
	// if octants[] supplied, method will return the octants that held the closest nodes
	final static double OCT_MARGIN = 1e-25;
	public int closestNode(FEM1 fem, int n, double[] distance, FEM1Octree[] octants) {
		
		double[] distC = {0}, nodeG = fem.node;
		FEM1Octree octant1;
		if (octants == null || octants[0] == null) {
			int n3 = n * 3;
			if (n3 > nodeG.length) {
					n3 -= nodeG.length;
					octant1 = locateCoordinate(fem.nodeWork[n3++], fem.nodeWork[n3++], fem.nodeWork[n3], MAX_LEVEL);
			} else	octant1 = locateCoordinate(nodeG[n3++], nodeG[n3++], nodeG[n3], MAX_LEVEL);
		} else
			octant1 = octants[0];
		
		int n3 = n * 3, nClosest;
		double x, y, z;
		if (n3 > nodeG.length) {	n3 -= nodeG.length; x = fem.nodeWork[n3++]; y = fem.nodeWork[n3++]; z = fem.nodeWork[n3]; }
		else				 	  {	x = nodeG[n3++]; y = nodeG[n3++]; z = nodeG[n3]; }		
		FEM1Octree octantC = null;
		
		// this is the case of only one node in container, start by finding the first ancestral relation branch besides direct parentage
		if (octant1.nodes == 1) {
			FEM1Octree octantRB = octant1;
			while ((octantRB.octantIterator >> 24) <= 1)						// backtrack until some other relation branch found
				if (octantRB == this) {											// special case: entire octree has only one node, so return itself
					if (distance != null) distance[0] = 0;
					if (octants != null) octants[1] = null;
					return n;
				}
				else octantRB = octantRB.parent;
			FEM1Octree[] octantsN = octantRB.leafOctantArray();					// gather all leaf octants of ancestor with relation branches

			double dist = Double.MAX_VALUE;
			for (FEM1Octree octantN : octantsN) {								// find the closest leaf octant in array
				if (octantN == null) break;
				if (octantN == octant1) continue;								// skip node's octant
				double distN = octantN.octantDistance(x, y, z, false);
				if (distN < dist) { dist = distN; octantC = octantN; }
			}
			nClosest = fem.closestNode(x, y, z, octantC.node, distC);			// find the closest node in closest leaf octant
		} else {
			// this is the base case of there being more than one node within container
			nClosest = fem.closestNode(x, y, z, octant1.node, distC);			// find out what local node that is closest
			if (octant1 == this) {												// if this is "archparent" octant (has no neighbours), we're done 
				if (distance != null) distance[0] = Math.sqrt(distC[0]);
				if (octants != null) octants[1] = null;
				return nClosest;
			}
		}
			
		// in both cases we have found a closest node candidate and the distance that will define a testing bounding box
		double dist = distC[0];
		// find the smallest ancestor octant that encompasses the bounding box
		double distRt = Math.sqrt(dist);
		double xMn = x - distRt, xPn = x + distRt, yMn = y - distRt, yPn = y + distRt, zMn = z - distRt, zPn = z + distRt;
		FEM1Octree octantE = octant1;
		// expand scope in the direction that isn't overshooting global bounding box and that still is farther than matching octree wall
		if (xMn > fem.bBox[0]) while (octantE.xM >= xMn && octantE != this) octantE = octantE.parent;
		if (xPn < fem.bBox[3]) while (octantE.xP <= xPn && octantE != this) octantE = octantE.parent;
		if (yMn > fem.bBox[1]) while (octantE.yM >= yMn && octantE != this) octantE = octantE.parent;
		if (yPn < fem.bBox[4]) while (octantE.yP <= yPn && octantE != this) octantE = octantE.parent;
		if (zMn > fem.bBox[2]) while (octantE.zM >= zMn && octantE != this) octantE = octantE.parent;
		if (zPn < fem.bBox[5]) while (octantE.zP <= zPn && octantE != this) octantE = octantE.parent;
		if (octantE == octant1)	{												// if smallest encompasser is node's octant itself, then we're done
			if (distance != null) distance[0] = distRt;
			if (octants != null) octants[1] = null;
			return nClosest;
		}
		// call method that collects only the octant's leaves that lie close enough to distance sphere of node's coordinate (criterion neighbours)
		FEM1Octree[] octantsN = octantE.leafOctantArrayByDistance(x, y, z, dist);
		FEM1Octree octantC2 = null;
		
		for (FEM1Octree octantN : octantsN) {									// time to check closeness of nodes of each close-enough leaf octant
			if (octantN == null) break;
			if (octantN == octant1 || octantN == octantC) continue;				// skip the already tested octant of the node
			int nClosestN = fem.closestNode(x, y, z, octantN.node, distC);		// find out closest node in criterion-neighbour octant
			if (distC[0] < dist) {												// if it was closer than closest node in this octant, assign to it
				nClosest = nClosestN; dist = distC[0]; octantC2 = octantN; }
		}
		if (distance != null) distance[0] = distRt;
		if (octants != null) octants[1] = octantC2;
		return nClosest;
	}

	
	
	// method applies closeness criterions of a coordinate to an octant, returning a distance
	double octantDistance(double x, double y, double z, boolean root) {
		double dx, dy, dz;
		if (x < xM) {
			if (y < yM) {		if (z < zM) { dx=xM-x; dy=yM-y; dz=zM-z; } else if (z > zP) { dx=xM-x; dy=yM-y; dz=z-zP; } else { dx=xM-x; dy=yM-y; dz=0; }
			} else if (y > yP) {if (z < zM) { dx=xM-x; dy=y-yP; dz=zM-z; } else if (z > zP) { dx=xM-x; dy=y-yP; dz=z-zP; } else { dx=xM-x; dy=y-yP; dz=0; }
			} else {			if (z < zM) { dx=xM-x; dy=0; dz=zM-z; }    else if (z > zP) { dx=xM-x; dy=0; dz=z-zP; }    else { dx=xM-x; dy=0; dz=0; }}
		} else if (x > xP) {
			if (y < yM) {		if (z < zM) { dx=x-xP; dy=yM-y; dz=zM-z; } else if (z > zP) { dx=x-xP; dy=yM-y; dz=z-zP; } else { dx=x-xP; dy=yM-y; dz=0; }
			} else if (y > yP) {if (z < zM) { dx=x-xP; dy=y-yP; dz=zM-z; } else if (z > zP) { dx=x-xP; dy=y-yP; dz=z-zP; } else { dx=x-xP; dy=y-yP; dz=0; }
			} else {			if (z < zM) { dx=x-xP; dy=0; dz=zM-z; }    else if (z > zP) { dx=x-xP; dy=0; dz=z-zP; }    else { dx=x-xP; dy=0; dz=0; }}
		} else {
			if (y < yM) {		if (z < zM) { dx=0; dy=yM-y; dz=zM-z; } else if (z > zP) { dx=0; dy=yM-y; dz=z-zP; } else { dx=0; dy=yM-y; dz=0; }
			} else if (y > yP) {if (z < zM) { dx=0; dy=y-yP; dz=zM-z; } else if (z > zP) { dx=0; dy=y-yP; dz=z-zP; } else { dx=0; dy=y-yP; dz=0; }
			} else {			if (z < zM) { dx=0; dy=0; dz=zM-z; }    else if (z > zP) { dx=0; dy=0; dz=z-zP; }
			  else { dx=0; dy=0; dz=0; }}	// this last case signifies that the coordinate lies inside octant
		}
		return root ? Math.sqrt(dx*dx + dy*dy + dz*dz) : dx*dx + dy*dy + dz*dz;
	}
	
	
	

	// method vectorises by centroid deviation into a scalar the level of nodal distribution disbalance between octants,
	// where 1 = max. disbalance (all nodes ending up in a single octant)
	double octantNodeDisbalance(int total, int nMMM, int nPMM, int nMPM, int nPPM, int nMMP, int nPMP, int nMPP, int nPPP) {
		double totalD = 0.5773 / total, vDisbX = 0, vDisbY = 0, vDisbZ = 0;
		if (nMMM > 0) { double nMMMf = (double)nMMM * totalD; vDisbX -= nMMMf; vDisbY -= nMMMf; vDisbZ -= nMMMf; }
		if (nPMM > 0) { double nPMMf = (double)nPMM * totalD; vDisbX += nPMMf; vDisbY -= nPMMf; vDisbZ -= nPMMf; }
		if (nMPM > 0) { double nMPMf = (double)nMPM * totalD; vDisbX -= nMPMf; vDisbY += nMPMf; vDisbZ -= nMPMf; }
		if (nPPM > 0) { double nPPMf = (double)nPPM * totalD; vDisbX += nPPMf; vDisbY += nPPMf; vDisbZ -= nPPMf; }
		if (nMMP > 0) { double nMMPf = (double)nMMP * totalD; vDisbX -= nMMPf; vDisbY -= nMMPf; vDisbZ += nMMPf; }
		if (nPMP > 0) { double nPMPf = (double)nPMP * totalD; vDisbX += nPMPf; vDisbY -= nPMPf; vDisbZ += nPMPf; }
		if (nMPP > 0) { double nMPPf = (double)nMPP * totalD; vDisbX -= nMPPf; vDisbY += nMPPf; vDisbZ += nMPPf; }
		if (nPPP > 0) { double nPPPf = (double)nPPP * totalD; vDisbX += nPPPf; vDisbY += nPPPf; vDisbZ += nPPPf; }
		return Math.sqrt(vDisbX * vDisbX + vDisbY * vDisbY + vDisbZ * vDisbZ);
	}
	
	
	// method tells how many extraneous nodes to allocate additional space for
	public static int extra_space(int n) { return (n >> 1) == 0 ? n + 1 : n + (n >> 1); }
	
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Octree, total branches: " + branches + "\n");
		sb.append(appendToString());
		sb.append("\n");
		return sb.toString();
	}
	
	private int maxItemPrint = 10;
	StringBuilder appendToString() {
		StringBuilder sb2 = new StringBuilder();
		for (int t = 0; t < level; t++) sb2.append("   ");
		switch(status) {
		case OCT_MMM: sb2.append("MMM: "); break; case OCT_MPM: sb2.append("MPM: "); break;
		case OCT_PMM: sb2.append("PMM: "); break; case OCT_PPM: sb2.append("PPM: "); break;
		case OCT_MMP: sb2.append("MMP: "); break; case OCT_MPP: sb2.append("MPP: "); break;
		case OCT_PMP: sb2.append("PMP: "); break; case OCT_PPP: sb2.append("PPP: "); break;
		}
		sb2.append(String.format(" x(%.3f,%.3f)y(%.3f,%.3f)z(%.3f,%.3f)\n", xM,xP,yM,yP,zM,zP));
		for (int t = 0; t < level; t++) sb2.append("   ");
		sb2.append("Level: " + level + ", nodes: " + nodes + ", edges: " + edges + ", branches: " + branches + "\n");
		if (node != null) {
			for (int t = 0; t < level; t++) sb2.append("   ");
			sb2.append("n:[");
			if (maxItemPrint < nodes)
					for (int n = 0; n < maxItemPrint; n++) sb2.append(node[n] + (n == maxItemPrint - 1 ? "...]\n" : ","));
			else	for (int n = 0; n < nodes; n++) sb2.append(node[n] + (n == nodes - 1 ? "]" : ","));
			sb2.append(edges > 0 ? " e:[" : "\n");
			if (maxItemPrint < edges)
					for (int e = 0; e < maxItemPrint; e++) sb2.append(edge[e] + (e == maxItemPrint - 1 ? "...]\n" : ","));
			else	for (int e = 0; e < edges; e++) sb2.append(edge[e] + (e == edges - 1 ? "]\n" : ","));
		}
		if (octant != null) {
			if (octant[0] != null) sb2.append(octant[0].appendToString());
			if (octant[1] != null) sb2.append(octant[1].appendToString());
			if (octant[2] != null) sb2.append(octant[2].appendToString());
			if (octant[3] != null) sb2.append(octant[3].appendToString());
			if (octant[4] != null) sb2.append(octant[4].appendToString());
			if (octant[5] != null) sb2.append(octant[5].appendToString());
			if (octant[6] != null) sb2.append(octant[6].appendToString());
			if (octant[7] != null) sb2.append(octant[7].appendToString());
		}
		return sb2;
	}
	
	public void toOBJ(FEM1 fem, boolean leafsOnly) {
		File file = new File(fem.name + "_octree.obj");
		if (!file.exists()) { try {	file.createNewFile(); } catch (IOException e) { e.printStackTrace(); }}
		BufferedWriter bw = null;
		try {		bw = new BufferedWriter(new FileWriter(file));
		} catch (IOException e) { e.printStackTrace(); }
		StringBuilder sb = new StringBuilder();
		String precFormat = "%." + fem.precision_OBJ + "f";
		vOBJcnt = 0; pOBJidx = 1;
		
		sb.append("# Visualisation of octree\n#\n");
		sb.append(appendVerticesToOBJ(fem, precFormat, leafsOnly));
		sb.append("# " + vOBJcnt + " vertices\n\ng " + fem.name + "_octree\n");
		sb.append(appendPolygonsToOBJ(fem, leafsOnly));
		sb.append("# " + (pOBJidx - 1) / 6 + " faces\n\ng\n");
		
		try {	bw.write(sb.toString());									// write out & close file
				bw.flush(); bw.close();
		} catch (IOException e) { e.printStackTrace(); }

	}
	
	StringBuilder appendVerticesToOBJ(FEM1 fem, String precFormat, boolean leafsOnly) {
		StringBuilder sb2 = new StringBuilder();
		if (!leafsOnly || (leafsOnly && octant == null)) {
			double[] pCoord = {xM,yM,zM, xP,yM,zM, xM,yP,zM, xP,yP,zM, xM,yM,zP, xP,yM,zP, xM,yP,zP, xP,yP,zP};
			for (int p = 0, c = 0; p < 8; p++)
				sb2.append("v  "+String.format(precFormat,pCoord[c++])+" "+String.format(precFormat,pCoord[c++])+" "+String.format(precFormat,pCoord[c++])+"\n");
			vOBJcnt += 8;
		}
		if (octant != null) {
			for (int o = 0; o < 8; o++) if (octant[o] != null) sb2.append(octant[o].appendVerticesToOBJ(fem, precFormat, leafsOnly)); }
		return sb2;
	}

	StringBuilder appendPolygonsToOBJ(FEM1 fem, boolean leafsOnly) {
		StringBuilder sb2 = new StringBuilder();
		if (!leafsOnly || (leafsOnly && octant == null)) {
			sb2.append("f " + pOBJidx + " " + (pOBJidx+1) + " " + (pOBJidx+5) + " " + (pOBJidx+4) + "\n");
			sb2.append("f " + (pOBJidx+1) + " " + (pOBJidx+3) + " " + (pOBJidx+7) + " " + (pOBJidx+5) + "\n");
			sb2.append("f " + (pOBJidx+3) + " " + (pOBJidx+2) + " " + (pOBJidx+6) + " " + (pOBJidx+7) + "\n");
			sb2.append("f " + pOBJidx + " " + (pOBJidx+2) + " " + (pOBJidx+6) + " " + (pOBJidx+4) + "\n");
			sb2.append("f " + pOBJidx + " " + (pOBJidx+1) + " " + (pOBJidx+3) + " " + (pOBJidx+2) + "\n");
			sb2.append("f " + (pOBJidx+4) + " " + (pOBJidx+5) + " " + (pOBJidx+7) + " " + (pOBJidx+6) + "\n");
			pOBJidx += 8;
		}
		if (octant != null) {
			for (int o = 0; o < 8; o++) if (octant[o] != null) sb2.append(octant[o].appendPolygonsToOBJ(fem, leafsOnly)); }
		return sb2;
	}


}
