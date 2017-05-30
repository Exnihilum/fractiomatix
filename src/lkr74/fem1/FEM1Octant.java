package lkr74.fem1;

import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class FEM1Octant {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			FEM1 OCTANT/SUBTREE PROCESSING METHODS														//
	//			Leonard Krylov 2017																			//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	// the class helps localising incoming FEM nodes and do a neighbourhood search on closest node 

	public final double xM,yM,zM,xP,yP,zP;	// boundaries of this octree octant (xM,yM,zM = x-,y-,z- and xP,yP,zP = x+,y+,z+)
	public final double xC, yC, zC;			// subdivision centroid of this octree octant
	volatile short status = IN_PROGRESS;	// status tells what octant it is and it's state
	final short level;						// level tells how high up we are in the octree
	volatile int int_ext = 0;				// bits flagging internal/external status of octant, including mid-edge & mid-face points
	public int branches=0,leaves=0;			// each subtree's branch & leaf count

	int nodes=0,nodesX=0;					// nodesX = extraneous nodes past the true nodes, temporaries for in-tree processing
	int edges=0,facets=0;
	int[] nodeI, edgeI, facetI;				// node, edge & facet references of this octant container
	double[] tmpR = null;					// temporary holder array for recursively processed data
	
	final FEM1Octant parent;				// backlink to parent
	public volatile FEM1Octant[] octant=null;		// array of 8 child octant references
	int iterator = 0;						// bitwise "array" holding the created octant indexes (helpful for nonrecursive tree searches)
	int enumerator = 0;						// primarily for ID numbering of octants, during building it enumerates octants that exceed maxNodes criterion
	
	float disbalance = 0;					// how disbalanced current branch is (=0 for a leaf)
	
	static int DEBUG_LEVEL;
	
	// codes for the 8 subdivisions of an octree node (M = minus, P = plus, ex: OCT_PMP = (x+,y-,z+))
	final static byte OCT_MMM=0, OCT_PMM=1, OCT_MPM=2, OCT_PPM=3, OCT_MMP=4, OCT_PMP=5, OCT_MPP=6, OCT_PPP=7;
	// following flags mask different children combinations to check their existence
	final static byte FLG_MMM=1, FLG_PMM=1<<1, FLG_MPM=1<<2, FLG_PPM=1<<3, FLG_MMP=1<<4, FLG_PMP=1<<5, FLG_MPP=1<<6, FLG_PPP=(byte)(1<<7);
	static final byte FLG_MXX = (byte)0x55, FLG_PXX = (byte)0xAA, FLG_XMX = (byte)0x33, FLG_XPX = (byte)0xCC, FLG_XXM = (byte)0x0F, FLG_XXP = (byte)0xF0;
	final static byte FLG_XMM = FLG_XMX&FLG_XXM, FLG_MXM = FLG_MXX&FLG_XXM, FLG_PXM = FLG_PXX&FLG_XXM, FLG_XPM = FLG_XPX&FLG_XXM;
	final static byte FLG_MMX = FLG_MXX&FLG_XMX, FLG_PMX = FLG_PXX&FLG_XMX, FLG_MPX = FLG_MXX&FLG_XPX, FLG_PPX = FLG_PXX&FLG_XPX;
	final static byte FLG_XMP = FLG_XMX&FLG_XXP, FLG_MXP = FLG_MXX&FLG_XXP, FLG_PXP = FLG_PXX&FLG_XXP, FLG_XPP = FLG_XPX&FLG_XXP;	

	public static final int DO_NODES = 1, DO_EDGES = 2, DO_FACETS = 4;	// flagging for processing of particular container data


	// special instantiator of a null-octant, which acts purely as an arrayspace occupier
	public FEM1Octant() { level = 0; xM = yM = zM = xP = yP = zP = xC = yC = zC = 0; parent = null; }
	
	// instantiates first octant, finding bounding box and setting up node reference array
	public FEM1Octant(FEM1Octree octree, int doWhat) {
		
		FEM1 fem = octree.fem;
		DEBUG_LEVEL = FEM1Octree.DEBUG_LEVEL;
		level = 0; parent = null;
		if (fem.bBox == null) fem.bBox = new double[6];	
		double[] nodeCoord = fem.node;
		int nodes3 = fem.nodes * FEM1.NCOORD;
		double xM1, yM1, zM1, xP1, yP1, zP1;
		xM1 = nodeCoord[0]; yM1 = nodeCoord[1]; zM1 = nodeCoord[2];				// arbitrarily choose first incoming vertex for initial bounding box
		xP1 = xM1; yP1 = yM1; zP1 = zM1;

		for (int n = FEM1.NCOORD; n < nodes3;) {
			if (nodeCoord[n] < xM1) xM1 = nodeCoord[n]; else if (nodeCoord[n] > xP1) xP1 = nodeCoord[n]; n++;
			if (nodeCoord[n] < yM1) yM1 = nodeCoord[n]; else if (nodeCoord[n] > yP1) yP1 = nodeCoord[n]; n++;
			if (nodeCoord[n] < zM1) zM1 = nodeCoord[n]; else if (nodeCoord[n] > zP1) zP1 = nodeCoord[n]; n++;
		}
		xM1-=OCT_MARGIN; yM1-=OCT_MARGIN; zM1-=OCT_MARGIN; xP1+=OCT_MARGIN; yP1+=OCT_MARGIN; zP1+=OCT_MARGIN;
		double xDlt = xP1 - xM1, yDlt = yP1 - yM1, zDlt = zP1 - zM1, maxBB = xDlt > yDlt ? xDlt : yDlt;
		if (maxBB < zDlt) maxBB = zDlt;
		// equalise bounding box sides to the one maximal side
		double xP2 = (xP1 + xM1 + maxBB) * .5; xM = fem.bBox[0] = (xP1 + xM1 - maxBB) * .5; xP = fem.bBox[3] = xP2;	
		double yP2 = (yP1 + yM1 + maxBB) * .5; yM = fem.bBox[1] = (yP1 + yM1 - maxBB) * .5; yP = fem.bBox[4] = yP2;
		double zP2 = (zP1 + zM1 + maxBB) * .5; zM = fem.bBox[2] = (zP1 + zM1 - maxBB) * .5; zP = fem.bBox[5] = zP2;
		xC = (xP + xM) * .5; yC = (yP + yM) * .5; zC = (zP + zM) * .5;
		
		if ((doWhat & DO_NODES) != 0) {
			nodeI = new int[fem.nodes <= octree.maxNodes ? extra_space(fem.nodes) : fem.nodes];
			int n = 0;
			while (n < fem.nodes) nodeI[n] = n++;									// set up initial node reference array
			nodes = fem.nodes;
		}
		if ((doWhat & DO_EDGES) != 0) edges = fem.edges;							// the counts signal to tree builder that these items must be included
		if ((doWhat & DO_FACETS) != 0) facets = fem.polygons;
	}
	
	
	// instantiates octree from boundary-intersecting lattice cubes, the final octree size supplied as boundbox zorners
	public FEM1Octant(FEM1Octree sourceTree, double xM0, double yM0, double zM0, double xP0, double yP0, double zP0) {	
		level = 0; parent = null;
		DEBUG_LEVEL = FEM1Octree.DEBUG_LEVEL;
		// note: since isosurface stuffing algorithm uses a (graded) BCC grid, we might need to include border lattice cubes and need
		// an extra buffer of one lattice step in the global bounding box
		xM = xM0; yM = yM0; zM = zM0; xP = xP0; yP = yP0; zP = zP0;
		xC = (xP + xM) * .5; yC = (yP + yM) * .5; zC = (zP + zM) * .5;
		nodeI = new int[nodes = sourceTree.fem.bLatticeNode.length / 3];
		int n = 0; while (n < nodes) nodeI[n] = n++;								// set up initial node reference array
	}

	
	// method instantiates skeleton child, if requested it will keep the IN_PROGRESS flag set
	public FEM1Octant(FEM1Octant parent, short status, boolean inProgress, boolean updateIterator) {
		this.parent = parent; this.status |= status;								// preserve IN_PROGRESS flag
		//parent.octant[status & 7] = this;											// TODO: pointing from parent to this octant done by client?
		if (updateIterator) parent.iteratorAdd(status & 7);							// updates both iterator and status bits
		level = (short)(parent.level + 1);
		switch (status) {
			case OCT_MMM: xM=parent.xM; yM=parent.yM; zM=parent.zM; xP=parent.xC; yP=parent.yC; zP=parent.zC; break;
			case OCT_PMM: xM=parent.xC; yM=parent.yM; zM=parent.zM; xP=parent.xP; yP=parent.yC; zP=parent.zC; break;
			case OCT_MPM: xM=parent.xM; yM=parent.yC; zM=parent.zM; xP=parent.xC; yP=parent.yP; zP=parent.zC; break;
			case OCT_PPM: xM=parent.xC; yM=parent.yC; zM=parent.zM; xP=parent.xP; yP=parent.yP; zP=parent.zC; break;
			case OCT_MMP: xM=parent.xM; yM=parent.yM; zM=parent.zC; xP=parent.xC; yP=parent.yC; zP=parent.zP; break;
			case OCT_PMP: xM=parent.xC; yM=parent.yM; zM=parent.zC; xP=parent.xP; yP=parent.yC; zP=parent.zP; break;
			case OCT_MPP: xM=parent.xM; yM=parent.yC; zM=parent.zC; xP=parent.xC; yP=parent.yP; zP=parent.zP; break;
			case OCT_PPP: xM=parent.xC; yM=parent.yC; zM=parent.zC; xP=parent.xP; yP=parent.yP; zP=parent.zP; break;
			default: xM=parent.xM; yM=parent.yM; zM=parent.zM; xP=parent.xP; yP=parent.yP; zP=parent.zP;
		}
		xC = (xP + xM) * .5; yC = (yP + yM) * .5; zC = (zP + zM) * .5;
		if (!inProgress) finish();
	}

	// full instantiator
	public FEM1Octant(FEM1Octant parent,
			short status, int nodes, int[] node, int edges, int[] edge, int facets, int[] facet, boolean inProgress, boolean updateIterator) {
		this.parent = parent; this.status |= status;								// preserve IN_PROGRESS flag
		//parent.octant[status & 7] = this;											// TODO: pointing from parent to this octant done by client?
		if (updateIterator) parent.iteratorAdd(status & 7);							// updates both iterator and status bits
		this.level = (short)(parent.level + 1);
		this.nodes = nodes; this.nodeI = node;
		this.edges = edges; this.edgeI = edge;
		this.facets = facets; this.facetI = facet;
		switch (status) {
			case OCT_MMM: xM=parent.xM; yM=parent.yM; zM=parent.zM; xP=parent.xC; yP=parent.yC; zP=parent.zC; break;
			case OCT_PMM: xM=parent.xC; yM=parent.yM; zM=parent.zM; xP=parent.xP; yP=parent.yC; zP=parent.zC; break;
			case OCT_MPM: xM=parent.xM; yM=parent.yC; zM=parent.zM; xP=parent.xC; yP=parent.yP; zP=parent.zC; break;
			case OCT_PPM: xM=parent.xC; yM=parent.yC; zM=parent.zM; xP=parent.xP; yP=parent.yP; zP=parent.zC; break;
			case OCT_MMP: xM=parent.xM; yM=parent.yM; zM=parent.zC; xP=parent.xC; yP=parent.yC; zP=parent.zP; break;
			case OCT_PMP: xM=parent.xC; yM=parent.yM; zM=parent.zC; xP=parent.xP; yP=parent.yC; zP=parent.zP; break;
			case OCT_MPP: xM=parent.xM; yM=parent.yC; zM=parent.zC; xP=parent.xC; yP=parent.yP; zP=parent.zP; break;
			case OCT_PPP: xM=parent.xC; yM=parent.yC; zM=parent.zC; xP=parent.xP; yP=parent.yP; zP=parent.zP; break;
			default: xM=parent.xM; yM=parent.yM; zM=parent.zM; xP=parent.xP; yP=parent.yP; zP=parent.zP;			
		}
		xC = (xP + xM) * .5; yC = (yP + yM) * .5; zC = (zP + zM) * .5;
		if (!inProgress) finish();
	}

	
	
	
	// method distributes nodes
	private void distributeNodes(FEM1Octree octree, int[] nNcnt, int maxItems, boolean multitask) {

		int nMMM=0, nPMM=nodes, nMPM=nPMM+nodes, nPPM=nMPM+nodes, nMMP=nPPM+nodes, nPMP=nMMP+nodes, nMPP=nPMP+nodes, nPPP=nMPP+nodes;
		nNcnt[1] = nPMM; nNcnt[2] = nMPM; nNcnt[3] = nPPM; nNcnt[4] = nMMP; nNcnt[5] = nPMP; nNcnt[6] = nMPP; nNcnt[7] = nPPP;
		int[] nodeOct = new int[nodes * 8];											// worst case allocation of octosplitted node arrays
		double[] nodeCoord = octree.node;

		// ********* NODE SEPARATION PART ********* //
		// this block eliminates extraneous internal nodes if a lattice tree gradation was requested
		// note: right now it eliminates the top gradation level's internal octants, leaving only boundary octants
		// the boundary flags are configured to remind of internal/external status of the removed suboctants
		if (octree.latticeTree && octree.gradations > 1 && level == octree.maxLevel - 1) {	
			short[] bLstatus = octree.fem.bLatticeStatus;
			for (int n = 0; n < nodes; n++) {											// distribute last level's node references into 8 octants
				int n3 = nodeI[n], inEx = bLstatus[n3] & 0xFF;
				n3 *= 3;
				double xN = nodeCoord[n3++], yN = nodeCoord[n3++], zN = nodeCoord[n3];	// access supplied global coordinates array
				if (xN<=xC) {
					if (yN<=yC) {	if (zN<=zC)	{	if(inEx!=0x1FF) nodeOct[nMMM++]=nodeI[n]; int_ext |= fMap_bFC[inEx&15] |      fMap_bFC[64+(inEx>>4)]; }
									else {			if(inEx!=0x1FF) nodeOct[nMMP++]=nodeI[n]; int_ext |= fMap_bFC[64+(inEx&15)] | fMap_bFC[128+(inEx>>4)]; }}
					else {			if (zN<=zC)	{	if(inEx!=0x1FF) nodeOct[nMPM++]=nodeI[n]; int_ext |= fMap_bFC[32+(inEx&15)] | fMap_bFC[96+(inEx>>4)]; }
									else 		{	if(inEx!=0x1FF) nodeOct[nMPP++]=nodeI[n]; int_ext |= fMap_bFC[96+(inEx&15)] | fMap_bFC[160+(inEx>>4)]; }}
				} else {
					if (yN<=yC) {	if (zN<=zC)	{	if(inEx!=0x1FF) nodeOct[nPMM++]=nodeI[n]; int_ext |= fMap_bFC[16+(inEx&15)] | fMap_bFC[80+(inEx>>4)]; }
									else {			if(inEx!=0x1FF) nodeOct[nPMP++]=nodeI[n]; int_ext |= fMap_bFC[80+(inEx&15)] | fMap_bFC[144+(inEx>>4)]; }}
					else {			if (zN<=zC)	{	if(inEx!=0x1FF) nodeOct[nPPM++]=nodeI[n]; int_ext |= fMap_bFC[48+(inEx&15)] | fMap_bFC[112+(inEx>>4)]; }
									else {			if(inEx!=0x1FF) nodeOct[nPPP++]=nodeI[n]; int_ext |= fMap_bFC[112+(inEx&15)]| fMap_bFC[176+(inEx>>4)]; }}
				}
			}
		} else {
			// this block does ordinary octree distribution of an octant's nodeset
			for (int n = 0; n < nodes; n++) {											// distribute node references into the 8 octants
				int n3 = nodeI[n] * 3;
				double xN = nodeCoord[n3++], yN = nodeCoord[n3++], zN = nodeCoord[n3];	// access supplied global coordinates array
				if (xN <= xC) {
					if (yN <= yC) {	if (zN <= zC)	{ nodeOct[nMMM++] = nodeI[n]; } else { nodeOct[nMMP++] = nodeI[n]; }}
					else {			if (zN <= zC)	{ nodeOct[nMPM++] = nodeI[n]; } else { nodeOct[nMPP++] = nodeI[n]; }}
				} else {
					if (yN <= yC) {	if (zN <= zC)	{ nodeOct[nPMM++] = nodeI[n]; } else { nodeOct[nPMP++] = nodeI[n]; }}
					else {			if (zN <= zC)	{ nodeOct[nPPM++] = nodeI[n]; } else { nodeOct[nPPP++] = nodeI[n]; }}
				}
			}
		}

		// ********* NODE DISTRIBUTION PART ********* //
		nNcnt[0] = nMMM; nNcnt[1] = nPMM; nNcnt[2] = nMPM; nNcnt[3] = nPPM; nNcnt[4] = nMMP; nNcnt[5] = nPMP; nNcnt[6] = nMPP; nNcnt[7] = nPPP;

		for (int o = 0; o < 8; o++) {												// create the 8 octants (if they contain items)
			int[] nodeChild;
			int nN = nNcnt[o] - nodes * o;
			if (nN > 0) {															// create this octant only if it received any items
				nodeChild = new int[nN <= maxItems ? extra_space(nN): nN];			// provide extra space only to probable leaf nodes
				for (int n = o*nodes, n2 = 0, nEnd = nNcnt[o]; n < nEnd; n2++, n++)
					nodeChild[n2] = nodeOct[n];
				if (multitask) {
					// need synchronisation since node/edge/facet splitters can work in parallel
					synchronized (this) {
						if (octant[o] == null)										// need to check for octant existence if doing concurrent processing
							octant[o] = new FEM1Octant(this,(short)o,nN,nodeChild,0,null,0,null,true,true);	// initialise octant with nodes (multitasked mode)
						else { octant[o].nodeI = nodeChild; octant[o].nodes = nN; }
					}
				} else {
					// note: flag inProgress=true because method can be multitasked in two ways, either the nodes/edges/facets get separated
					// in parallel, or an entire branch is split by a dedicated thread, so any processor must wait until flags are cleared backwards
					octant[o] = new FEM1Octant(this,(short)o,nN,nodeChild,0,null,0,null,true,false);
					iteratorAdd(o);
				}
			} else nodeChild = null;
			nNcnt[o] = nN;															// remove offset from node count
		}
	}
	
	
	private boolean distributeEdges(FEM1 fem, int[] nEcnt, int maxItems, boolean multitask) {

		nEcnt[1] = edges; nEcnt[2] = edges*2; nEcnt[3] = edges*3; nEcnt[4] = edges*4; nEcnt[5] = edges*5; nEcnt[6] = edges*6; nEcnt[7] = edges*7;
		
		int e1, pageSz = ((edges+8)*6)/8;											// pageSz = edge segment coords page allocation size
		double[] nodea, nodeb;
		int distribFailureTest = 0xFF;

		int[] edgeOct = new int[edges*8];											// temporary distribution array of edge indexes
		double[][] edgeCrd = new double[8*8][];										// prematurely assume we will have a perfect even split-up (unlikely)
		for (int c = 0; c < 8*8; c+= 8) edgeCrd[c] = new double[pageSz];			// the trick is to allocate further "pages" if necessary
		int[] edgeCrdP = new int[8], edgeCnt = new int[8*8];						// holds counters of the pages
		double[] split = new double[6*8];

		int[] pEdgeN = fem.edgeNode;
		// ********* EDGE SEPARATION PART ********* //
		for (int e=0, e3b=0, oC=0, bits=0; e < edges; e++) {						// then distribute edge references into the 8 octants

			if (tmpR == null) {							// if no local segments found, gather edge coordinates from FEM system's edgeNode[] referencer
				e1 = e;
				int ee = e * 2, ea3 = pEdgeN[ee++] * 3, eb3 = pEdgeN[ee] * 3;
				if (fem.nodeWork == null) {
					nodea = nodeb = fem.node;
				} else {
					if (ea3 >= fem.node.length) { ea3 -= fem.node.length; nodea = fem.nodeWork; } else nodea = fem.node;
					if (eb3 >= fem.node.length) { eb3 -= fem.node.length; nodeb = fem.nodeWork; } else nodeb = fem.node;
				}
				bits = edgeMembership(nodea[ea3++], nodea[ea3++], nodea[ea3++], nodeb[eb3++], nodeb[eb3++], nodeb[eb3++], split, 0, 0, 0xFF);
			} else {															// for rest of tree, read from locally temporary stored edge segments
				e1 = edgeI[e];
				bits = edgeMembership(tmpR[e3b++], tmpR[e3b++], tmpR[e3b++], tmpR[e3b++], tmpR[e3b++], tmpR[e3b++], split, 0, 0, 0xFF);
			}
			distribFailureTest &= bits;
			
			int o = bits < 0x0F ? 0 : ((bits & 0x0F)==0 ? 4 : 0); if (o == 4) bits >>= 4;
			while (bits != 0) {														// pick out the octants that the edge was segmented into
				if ((bits & 1) == 1) {
					edgeOct[nEcnt[o]++] = e1; oC = o * 6;
					int eCofs = o * 8 + edgeCrdP[o];
					int e3a = edgeCnt[eCofs];										// edgeCnt[] holds the indexes into current pages
					if (pageSz < e3a+6) {											// allocate another page of coordinates if necessary
						edgeCrd[++eCofs] = new double[pageSz];
						edgeCrdP[o]++; e3a = 0;
					}
					double[] edgeCrd1 = edgeCrd[eCofs];
					edgeCrd1[e3a++]=split[oC++]; edgeCrd1[e3a++]=split[oC++]; edgeCrd1[e3a++]=split[oC++];	// add edge segments to temporary edgeCrd[] arrays
					edgeCrd1[e3a++]=split[oC++]; edgeCrd1[e3a++]=split[oC++]; edgeCrd1[e3a++]=split[oC++];
					edgeCnt[eCofs] += 6;
				}
				bits >>= 1; o++;
			}
		}
		tmpR = null;																// enforce deallocation of current tmpR[] array
		if (distribFailureTest == 0xFF) {
			System.out.println("FEM1.distributeEdges(): distribution criterion unattainable in:\n" + this.toString());
			return false;															// if all edges ended up in ALL octants, distrib.criterion unattainable
		}
		
		// ********* EDGE DISTRIBUTION PART ********* //
		int[] edgeChild = null;	
		for (int o = 0; o < 8; o++) {												// create remaining octants necessary for holding edges
			int eN = nEcnt[o] - o * edges;
			double[] tmpR1 = null;
			if (eN > 0) {
				edgeChild = new int[eN <= maxItems ? extra_space(eN) : eN];			// provide extra space only to probable leaf nodes
				for (int e = o*edges, e2=0, eEnd = nEcnt[o]; e < eEnd; e2++, e++) edgeChild[e2] = edgeOct[e];
				int eCP = edgeCrdP[o];
				tmpR1 = new double[eCP * pageSz + edgeCnt[o*8 + eCP]];
				for (int p = 0, et3 = 0; p <= eCP; p++) {							// join edge segment coordinates from pages into a single tmpR[] array
					double[] eCrd = edgeCrd[o*8 + p];
					for (int e3=0, e3end = edgeCnt[o*8 + p]; e3 < e3end; e3++) tmpR1[et3++] = eCrd[e3]; }
				
				if (multitask) {
					synchronized (this) {											// need synchronisation since node & facet splitter can work in parallel
						if (octant[o] == null) {									// if octant holding edges was missing
								octant[o] = new FEM1Octant(this, (short)o, 0, null, eN, edgeChild, 0, null,true,true);	// initialise it
						} else {octant[o].edges = eN; octant[o].edgeI = edgeChild; }
					}
				} else {
					if (octant[o] == null) {										// if octant holding edges was missing
							octant[o] = new FEM1Octant(this, (short)o, 0, null, eN, edgeChild, 0, null,true,false);		// initialise it
							iteratorAdd(o);
					} else {octant[o].edges = eN; octant[o].edgeI = edgeChild; }
				}
				octant[o].tmpR = tmpR1;
			}
			nEcnt[o] = eN;															// remove offset from edge count
		}
		return true;
	}
	
	// facet distributor is able to take care of degenerate case when same facets end up in all suboctants, indicating
	// a failure of the facet count criterion: count per octant cannot be achieved since constituent facet boundboxes always cover that octant's space
	// such a case is impossible for nondimensional objects like points
	private boolean distributeFacets(FEM1 fem, int[] nFcnt, int maxItems, boolean multitask) {

		nFcnt[1] = facets; nFcnt[2] = facets*2; nFcnt[3] = facets*3; nFcnt[4] = facets*4; nFcnt[5] = facets*5; nFcnt[6] = facets*6; nFcnt[7] = facets*7;
		int distribFailureTest = 0xFF;

		int[] facetOct = new int[facets*8];											// temporary distribution array of facet indexes
		// ********* FACET SEPARATION PART ********* //
		for (int f = 0; f < facets; f++) {
			int f1 = facetI == null ? f : facetI[f];								// if this is root branch, use iteration index as facet index
			double[] fBBox = fem.facetBBox(f1);
			if (fBBox != null) {													// accept only facets from the polygon array
				int bits = octantBBoxOverlap(fBBox[0], fBBox[1], fBBox[2], fBBox[3], fBBox[4], fBBox[5]) & 0xFF;	// get facet's distribution in this branch
				distribFailureTest &= bits;
				int o = bits < 0x0F ? 0 : ((bits & 0x0F)==0 ? 4 : 0); if (o == 4) bits >>= 4;
				while (bits != 0) {													// pick out the octants that the facet bounding box was located in
					if ((bits & 1) == 1) facetOct[nFcnt[o]++] = f1;
					bits >>= 1; o++;
				}
			}
		}
		if (distribFailureTest == 0xFF) {
			if (DEBUG_LEVEL > 1) System.out.println("FEM1.distributeFacets(): distribution criterion unattainable in:\n" + this.toString());
			return false;															// all facets ended up in ALL octants, distrib.criterion unattainable
		}

		// ********* FACET DISTRIBUTION PART ********* //
		int[] facetChild = null;
		for (int o = 0; o < 8; o++) {												// create remaining octants necessary for holding facets
			int fN = nFcnt[o] - o * facets;
			if (fN > 0) {
				facetChild = new int[fN <= maxItems ? extra_space(fN) : fN];      	// make extra space only if this is an expected leaf node
				for (int f = o*facets, f2=0, fEnd = nFcnt[o]; f < fEnd; f2++, f++) facetChild[f2] = facetOct[f];
				if (multitask) {
					synchronized (this) {											// need synchronisation since node & edge splitter can work in parallel
						if (octant[o] == null) {									// if octant holding facets was missing
								octant[o] = new FEM1Octant(this, (short)o, 0, null, 0, null, fN, facetChild, true, true);	// initialise it
						} else {octant[o].facets = fN; octant[o].facetI = facetChild; }
					}
				} else {
					if (octant[o] == null) {										// if octant holding facets was missing
							octant[o] = new FEM1Octant(this, (short)o, 0, null, 0, null, fN, facetChild, true, false);		// initialise it
							iteratorAdd(o);
					} else {octant[o].facets = fN; octant[o].facetI = facetChild; }
				}
			}
			nFcnt[o] = fN;
		}
		return true;
	}


	
	// splits current octree node into 8 partitions, divides up the items, returns rate of disbalance within the current branch
	// edges are recursively split into segments which traverse as temporary coordinates up the tree
	// facets are for now expected to be small enough to fulfill only bounding box overlap criterion
	public int split(FEM1Octree octree, boolean multitask, short enforcedLevel) {
		
		if (octree.latticeTree && octree.gradations > 1 && level == octree.maxLevel - 1) {
			// if gradation of a lattice tree volume's level of detail is requested, we can create first level already here
			// by not doing the final subdivision of the octants whose children are all flagged as fully internal
			short[] bLstatus = octree.fem.bLatticeStatus;
			if ( nodes == 8 && (bLstatus[nodeI[0]] & bLstatus[nodeI[1]] & bLstatus[nodeI[2]] & bLstatus[nodeI[3]] & 
								bLstatus[nodeI[4]] & bLstatus[nodeI[5]] & bLstatus[nodeI[6]] & bLstatus[nodeI[7]] & 0x1FF) == 0x1FF) {
				int_ext = FULL_INTERNAL;
				nodeI = null; nodes = 0; nodesX = level | 0xFF<<16;		// indicate gradation level in nodesX & set IST-internation flags
				return 0;
			}
		}
		
		int[] nNcnt = new int[8], nEcnt = new int[8], nFcnt = new int[8], nIcnt = new int[8];	// nIcnt[] holds total item count per octant
		octant = new FEM1Octant[8];

		// do concurrent node & edge & facet distribution on root node, the criterions are either that this is initial distributed work on root,
		// or the tree is completed, but added elements motivate a branch split, which can be safely done concurrently
		final AtomicInteger pCtr = new AtomicInteger(0);						// make every separate mesh op. decrement a local count
		if (nodes > 0) {
			// DEBUG: relative element counts taken from measuring relative speeds of node/edge/facet splits, which ofcourse
			// will only balance well if all three operations below are enacted at once
			if (multitask && nodes > octree.threadSplitF*octree.splitNf && FEM1.freeTasks.get() > 0) {
				pCtr.incrementAndGet(); FEM1.freeTasks.decrementAndGet();
				FEM1.getExecutor(0).execute(new Runnable() {@Override public void run() {
					distributeNodes(octree, nNcnt, octree.maxNodes, true);
					pCtr.decrementAndGet(); FEM1.freeTasks.incrementAndGet(); } });
			} else {
				//long timer = level == 0 ? System.nanoTime() : 0;				// DEBUG: gets approx.timing of operation
				distributeNodes(octree, nNcnt,  octree.maxNodes, false);
				//if (level == 0) octree.splitNtiming = (octree.splitNtiming + (System.nanoTime() - timer)) / 2;
			}
		}
		if (edges > 0) {
			if (multitask && edges > octree.threadSplitF*octree.splitEf && FEM1.freeTasks.get() > 0) {
				pCtr.incrementAndGet(); FEM1.freeTasks.decrementAndGet();
				FEM1.getExecutor(0).execute(new Runnable() {@Override public void run() {
					distributeEdges(octree.fem, nEcnt, octree.maxEdges, true); 
					pCtr.decrementAndGet(); FEM1.freeTasks.incrementAndGet(); } });
			} else {
				//long timer = level == 0 ? System.nanoTime() : 0;				// DEBUG: gets approx.timing of operation
				distributeEdges(octree.fem, nEcnt, octree.maxEdges, false);
				//if (level == 0) octree.splitEtiming = (octree.splitEtiming + (System.nanoTime() - timer)) / 2;
			}
		}
		if (facets > 0) {
			if (multitask && facets > octree.threadSplitF*octree.splitFf &&  FEM1.freeTasks.get() > 0) {
				pCtr.incrementAndGet(); FEM1.freeTasks.decrementAndGet();
				FEM1.getExecutor(0).execute(new Runnable() {@Override public void run() {
					distributeFacets(octree.fem, nFcnt, octree.maxFacets, true);
					pCtr.decrementAndGet(); FEM1.freeTasks.incrementAndGet(); } });
			} else {
				//long timer = level == 0 ? System.nanoTime() : 0;				// DEBUG: gets approx.timing of operation
				distributeFacets(octree.fem, nFcnt, octree.maxFacets, false);
				//if (level == 0) octree.splitFtiming = (octree.splitFtiming + (System.nanoTime() - timer)) / 2;
			}
		}
		long timer = System.currentTimeMillis();
		while (pCtr.get() > 0)											// wait for the three steps to finish, we need results for further splitting
			if (System.currentTimeMillis() - timer > 10000) {			// certify that no unexpected lockups happen
				FEM1.getExecutor(0).shutdown();							// TODO: adapt approximate waiting time to datasize
				try { FEM1.getExecutor(0).awaitTermination(8, TimeUnit.SECONDS);
				} catch (InterruptedException e) { e.printStackTrace(); }				
				throw new RuntimeException("FEM1Octant.split(): Completion wait timeout.");
			}		
				
		// enumerate the octants eligible for further splitting in "enumerator"
		boolean levelIncreased = false;
		int facetSum = 0, edgeSum = 0, octantsMN = 0;
		for (int o = 0; o < 8; o++) {
			FEM1Octant octant1 = octant[o];
			if (octant1 != null) {
				levelIncreased = true;
				facetSum += nFcnt[o]; edgeSum += nEcnt[o];
				// any group, nodes/edges/facets, can be large enough to provoke a split, or enforcedLevel can force attainment a specified level
				// note: set enforcedLevel to 0 to deactivate level enforcement criterion
				if (octant1.level < enforcedLevel || nNcnt[o] > octree.maxNodes || nEcnt[o] > octree.maxEdges || nFcnt[o] > octree.maxFacets) {
					enumerator = (enumerator<<3) | o; octantsMN++; }
				else {
					// for lattice tree, assign internal/external status from the lattice cube array's status words
					if (octree.latticeTree) {
						int bnd = octree.fem.bLatticeStatus[octant1.nodeI[0]];
						octant1.int_ext = bnd<<9 & CENTROID_INTERNAL | (bnd&0xFF);
						if (octant1.internal()) octant1.nodesX = octant1.level;	// if leaf was internal, assign it as an IST-internal octant
						octant1.nodeI = null; octant1.nodes = 0;				// node index not needed anymore
					}
					octant1.finish();											// nodes that are definitive leaves are no longer IN_PROGRESS on creation
					leaves++;
				}
				if (!multitask) octree.levelSums[level+1]++;
				nIcnt[o] = nNcnt[o] + nEcnt[o] + nFcnt[o];
			}
		}

		if (levelIncreased) {
			if (level + 1 > octree.topLevel) octree.topLevel = (short)(level + 1);	// keep track of maximal level attained in tree
			enumerator |= (octantsMN << 24);
			nodeI = edgeI = facetI = null;			// a branch doesn't need item references (but keeping nodes count of entire branch can be useful)
			branches = 1;
		} else octant = null;													// no children could be created
		// branch disbalance calculated from items clustered across octants, set disbalanceCap=1 to stop disbalance heuristic from deciding splits
		// note: since facet count of parent != facet sum of childre, facetSum & edgeSum are used to stop disbalance from exceeding 1.0
		if (octree.disbalanceCap < 1)
			disbalance = octantNodeDisbalance(nodes + edgeSum + facetSum, nIcnt[0], nIcnt[1], nIcnt[2], nIcnt[3], nIcnt[4], nIcnt[5], nIcnt[6], nIcnt[7]);
		return octantsMN;
	}
	
	
	// method recursively constructs octree with at most maxItems items per container and splitting limited by maximal branch disbalance
	// method will send off individual subtree constructions to threads (if current thread can't do it itself) if multitask=true
	public int build(FEM1Octree octree, double disbalanceGrowth, boolean multitask) {
		
		boolean backTrack = true;
		// if branch levels are not overrun & branch is not disbalanced
		if (level < octree.maxLevel && disbalanceGrowth < octree.disbalanceCap) {
			int splittable = split(octree, multitask, octree.enforcedLevel);		// split into octants (if possible), becoming a branch

			int octantIt = enumerator;
			if (splittable > 0) backTrack = false;									// make sure to backtrack on branches that didn't generate splittable leaves
			
			for (int o = 0; o < splittable; o++, octantIt >>= 3) {					// note: we recurse only octants that contain more than "maxItems" items
				FEM1Octant octantC = this.octant[octantIt & 7];
				// if multitask requested and this branch contains enough items to bother with threading overhead
				// and this isn't the last call (which this task can complete itself)
				if (multitask && FEM1.freeTasks.get() > 0 && nodes + edges + facets >= octree.threadSpawnItemsCap && o < splittable-1) {
					 FEM1.freeTasks.decrementAndGet();
					 FEM1.getExecutor(0).execute(new Runnable() {
						@Override public void run() {
							octantC.build(octree, (disbalance + disbalanceGrowth) * 0.5, multitask);
							FEM1.freeTasks.incrementAndGet();
						}
					});
				} else {
					branches += octantC.build(octree, (disbalance + disbalanceGrowth) * 0.5, multitask);
					leaves += octantC.octant == null ? 1 : octantC.leaves;
				}
			}
			if (octree.latticeTree) {											// if doing a lattice tree
				if (!multitask && splittable > 0) {
					int_extFromChildren(true);									// on split of a latticeTree branch, resolve parent's internal status
					internOfISTFromChildren(); }								// to cull octants within IST gradations, get their IST-internality status
				nodes = 0;														// nodes parameter will be used for a different purpose								
			}
			if (parent == null) {
				finish();
				if (octree.latticeTree) nodes = 0;
				// if we're at root level in multitask mode, wait for all threads to finish, then count the branches properly
				// note: branches CANNOT be recursively counted by threads, since each parent exits before subbranches are finished
				// fix the incorrect branch counts (will also block up on unfinished nodes) and for a lattice tree, sort out internal status of branches
				if (multitask) process(octree, OCTREE_SUM_BRANCHES|OCTREE_SUM_LEVELS|(octree.latticeTree ? OCTREE_IST : 0));
//				else if (octree.latticeTree && !octree.limitGradation)			// if >2 gradations wanted, cull unnecessary IST-internal leaves
//					process(octree, OCTREE_IST_CULL);
			}
		}
		
		if (backTrack) {														// if split generated no consecutive recursions or if this is a leaf
			FEM1Octant octant = this;
			while (octant != null && !octant.finished()) {						// backtrack from finished leaf clearing IN_PROGRESS as far as necessary
				octant.finish(); octant = octant.parent;
			}
		}
		return branches;
	}

	
	// method applies processing actions of choice on octree
	public static final int OCTREE_WIPE_EDGES = 1, OCTREE_SUM_BRANCHES = 2, OCTREE_FIX_ITERATORS = 4, OCTREE_ENUMERATE = 8;
	public static final int OCTREE_SUM_LEVELS = 16, OCTREE_SUMMING_LEVELS = 32, OCTREE_IST = 64;
	public int process(FEM1Octree octree, int action) {
		
		long timer = System.currentTimeMillis();
		while (in_progress()) {													// wait for branch to complete (for multitask mode)
			if (System.currentTimeMillis() - timer > 10000) {
				FEM1.getExecutor(0).shutdown();
				try { FEM1.getExecutor(0).awaitTermination(8, TimeUnit.SECONDS);
				} catch (InterruptedException e) { e.printStackTrace(); }				
				throw new RuntimeException("FEM1Octant.process(): Completion wait timeout on level " + level + ".");
			}
		}
		
		if ((action & OCTREE_SUM_LEVELS) != 0) {
			if ((action & OCTREE_SUMMING_LEVELS) == 0) {
				action |= OCTREE_SUMMING_LEVELS;									// on first entry, reset level sums, switch flag
				octree.levelSums = new int[octree.topLevel + 1];					// levels summable from any octant's aspect
			}
			octree.levelSums[level]++;												// sum up distribution across levels
		}

		if ((action & OCTREE_ENUMERATE) != 0) {									// if octant ID enumeration requested
			if (parent == null) octree.enumerator = 0;
			enumerator = octree.enumerator++;
		}
		
		
		if (octant!=null) {														// if this is a branch
			if ((action & OCTREE_FIX_ITERATORS) != 0) {							// if iterator recalculation requested
				iterator = 0; status &= 0xFFFF - (0xFF<<3);
				for (int o = 7; o >= 0; o--) if (octant[o] != null) iteratorAdd(o);
			}
			int branchSum = 1, leafSum = 0;										// not forgetting to sum in the current branch octant
			int octantIt = iterator, oEnd = iterator >> 24;
			for (int o = 0; o < oEnd; o++, octantIt >>= 3) {
				FEM1Octant octantC = octant[octantIt & 7];
				branchSum += octantC.process(octree, action);
				leafSum += octantC.octant == null ? 1 : octantC.leaves;			// sum the leaves
			}
			if ((action & OCTREE_IST) != 0) {
				int_extFromChildren(true);										// we only need to find out internal status for branches
				internOfISTFromChildren();
			}
			if ((action & OCTREE_SUM_BRANCHES) != 0) { branches = branchSum; leaves = leafSum; }
			return branches;
		}
		if ((action & OCTREE_WIPE_EDGES) != 0) tmpR = null;
		return 0;																// leaves don't count as branches
	}

	
	
	// method returns octant's same-level face neighbours (if they exist, otherwise lower level octants),
	// neighbour storage array and level retrace arrays can be supplied
	// octants are enumerated by priority x-,x+,y-,y+,z-,z+
	public FEM1Octant[] getFaceNeighbours(FEM1Octree octree, FEM1Octant[] neighbour, byte[] retrace) {
		//if (parent == null) return null;			// DEBUG: assume caller knows that obviously, method is not appliable on root octant
		if (neighbour == null) neighbour =  new FEM1Octant[7];
		if (retrace == null) retrace = new byte[octree.topLevel];
		short n0=0, n1=0, n2=0, n3=0, n4=0 ,n5=0, o0=0, o1=0, o2=0;
		switch (status & 7) {	
		case OCT_MMM:	n0=0; n1=1; n2=2; n3=3; n4=4; n5=5; o0=OCT_PMM; o1=OCT_MPM; o2=OCT_MMP; break;
		case OCT_PMM:	n0=1; n1=0; n2=2; n3=3; n4=4; n5=5; o0=OCT_MMM; o1=OCT_PPM; o2=OCT_PMP; break;
		case OCT_MPM:	n0=0; n1=1; n2=3; n3=2; n4=4; n5=5; o0=OCT_PPM; o1=OCT_MMM; o2=OCT_MPP; break;
		case OCT_PPM:	n0=1; n1=0; n2=3; n3=2; n4=4; n5=5; o0=OCT_MPM; o1=OCT_PMM; o2=OCT_PPP; break;
		case OCT_MMP:	n0=0; n1=1; n2=2; n3=3; n4=5; n5=4; o0=OCT_PMP; o1=OCT_MPP; o2=OCT_MMM; break;
		case OCT_PMP:	n0=1; n1=0; n2=2; n3=3; n4=5; n5=4; o0=OCT_MMP; o1=OCT_PPP; o2=OCT_PMM; break;
		case OCT_MPP:	n0=0; n1=1; n2=3; n3=2; n4=5; n5=4; o0=OCT_PPP; o1=OCT_MMP; o2=OCT_MPM; break;
		case OCT_PPP:	n0=1; n1=0; n2=3; n3=2; n4=5; n5=4; o0=OCT_MPP; o1=OCT_PMP; o2=OCT_PPM;
		}
		FEM1Octant[] octantP = parent.octant;
		neighbour[n0] = symmetricDistantNeighbour(o0, 1, retrace, false); neighbour[n1] = octantP[o0];
		neighbour[n2] = symmetricDistantNeighbour(o1, 2, retrace, false); neighbour[n3] = octantP[o1];	
		neighbour[n4] = symmetricDistantNeighbour(o2, 4, retrace, false); neighbour[n5] = octantP[o2];
		return neighbour;
	}
	
	// array for scanning neighbour aspects within neighbours1() method
	final static byte[] nbrs1aspect = {
			6,OCT_MPP, 5,OCT_PMP, 4,OCT_MMP, 5,OCT_PMP, 6,OCT_MPP, 3,OCT_PPM, 2,OCT_MPM, 3,OCT_PPM, 1,OCT_PMM,
			-1,OCT_PMM, 3,OCT_PPM, -1,OCT_MPM, -1,OCT_PPM, 6,OCT_MPP, 5,OCT_PMP, -1,OCT_MMP, -1,OCT_PMP, -1,OCT_MPP,
				6,OCT_PPP, 5,OCT_MMP, 4,OCT_PMP, 5,OCT_MMP, 6,OCT_PPP, 3,OCT_MPM, 2,OCT_PPM, 3,OCT_MPM, -1,OCT_MMM,
				1,OCT_MMM, -1,OCT_MPM, -1,OCT_PPM, 3,OCT_MPM, 6,OCT_PPP, -1,OCT_MMP, -1,OCT_PMP, 5,OCT_MMP, -1,OCT_PPP,
			6,OCT_MMP, 5,OCT_PPP, 4,OCT_MPP, 5,OCT_PPP, 6,OCT_MMP, 3,OCT_PMM, -1,OCT_MMM, -1,OCT_PMM, 1,OCT_PPM,
			-1,OCT_PPM, 3,OCT_PMM, 2,OCT_MMM, 3,OCT_PMM, -1,OCT_MMP, 5,OCT_PPP, -1,OCT_MPP, -1,OCT_PPP, 6,OCT_MMP,
				6,OCT_PMP, 5,OCT_MPP, 4,OCT_PPP, 5,OCT_MPP, 6,OCT_PMP, -1,OCT_MMM, -1,OCT_PMM, 3,OCT_MMM, -1,OCT_MPM,
				1,OCT_MPM, 3,OCT_MMM, 2,OCT_PMM, 3,OCT_MMM, -1,OCT_PMP, -1,OCT_MPP, -1,OCT_PPP, 5,OCT_MPP, 6,OCT_PMP,
			6,OCT_MPM, 5,OCT_PMM, -1,OCT_MMM, -1,OCT_PMM, -1,OCT_MPM, 3,OCT_PPP, 2,OCT_MPP, 3,OCT_PPP, 1,OCT_PMP,
			-1,OCT_PMP, 3,OCT_PPP, -1,OCT_MPP, -1,OCT_PPP, 6,OCT_MPM, 5,OCT_PMM, 4,OCT_MMM, 5,OCT_PMM, 6,OCT_MPM,
				6,OCT_PPM, -1,OCT_MMM, -1,OCT_PMM, 5,OCT_MMM, -1,OCT_PPM, 3,OCT_MPP, 2,OCT_PPP, 3,OCT_MPP, -1,OCT_MMP,
				1,OCT_MMP, -1,OCT_MPP, -1,OCT_PPP, 3,OCT_MPP, 6,OCT_PPM, 5,OCT_MMM, 4,OCT_PMM, 5,OCT_MMM, 6,OCT_PPM,
			-1,OCT_MMM, 5,OCT_PPM, -1,OCT_MPM, -1,OCT_PPM, 6,OCT_MMM, 3,OCT_PMP, -1,OCT_MMP, -1,OCT_PMP, 1,OCT_PPP,
			-1,OCT_PPP, 3,OCT_PMP, 2,OCT_MMP, 3,OCT_PMP, 6,OCT_MMM, 5,OCT_PPM, 4,OCT_MPM, 5,OCT_PPM, 6,OCT_MMM,
				-1,OCT_PMM, -1,OCT_MPM, -1,OCT_PPM, 5,OCT_MPM, 6,OCT_PMM, -1,OCT_MMP, -1,OCT_PMP, 3,OCT_MMP, -1,OCT_MPP,
				1,OCT_MPP, 3,OCT_MMP, 2,OCT_PMP, 3,OCT_MMP, 6,OCT_PMM, 5,OCT_MPM, 4,OCT_PPM, 5,OCT_MPM, 6,OCT_PMM };
	// the skiplist checker of facets bits assumes that 32nd bit is ALWAYS zero
	final static int skI = 0x80000000;
	final static int[] skipList_nbrs1 = {skI,skI,1<<6,skI,skI,skI,1<<7,skI,1<<8,1<<9,skI,1<<10,skI,skI,skI,1<<11,skI,skI};

	// method returns same-level face & edge neighbours (if they exist, otherwise lower level octants or nulls)
	// ofsIdx tells which offser of supplied neighbourhood array to start writing to (fitting for batch searches)
	// oFlags tell which neighbour aspect to check for a neighbour, one bit per aspect
	public FEM1Octant[] getFaceEdgeNeighbours(FEM1Octree octree, FEM1Octant[] neighbour, int ofsIdx, int oFlags) {

		if (neighbour == null) neighbour = new FEM1Octant[19];
		double step = xC - xM + OCT_MARGIN;
		
		for (int i = (status&7) * 36, a = 0; a < 18; a++, oFlags>>=1) {
			int a1 = ofsIdx+a;
			neighbour[a1] = null;											// DEBUG: if a neighbour array is reused by caller, clearing it is a must
			if ((oFlags&1)==0) { i += 2; continue; }						// check if this aspect's neighbour was requested
			
			if (nbrs1aspect[i++] >= 0) {									// get the plane-mirror flags for symmetric neighbour location
				i++;														// get the neighbour code we're looking for
				switch (a) {
				case  0: neighbour[a1] = octree.root.locateCoordinate(xC, yC-step, zC-step, level); break;
				case  1: neighbour[a1] = octree.root.locateCoordinate(xC-step, yC, zC-step, level); break;
				case  2: neighbour[a1] = octree.root.locateCoordinate(xC, yC, zC-step, level); break;
				case  3: neighbour[a1] = octree.root.locateCoordinate(xC+step, yC, zC-step, level); break;
				case  4: neighbour[a1] = octree.root.locateCoordinate(xC, yC+step, zC-step, level); break;
				case  5: neighbour[a1] = octree.root.locateCoordinate(xC-step, yC-step, zC, level); break;
				case  6: neighbour[a1] = octree.root.locateCoordinate(xC, yC-step, zC, level); break;
				case  7: neighbour[a1] = octree.root.locateCoordinate(xC+step, yC-step, zC, level); break;
				case  8: neighbour[a1] = octree.root.locateCoordinate(xC-step, yC, zC, level); break;
				case  9: neighbour[a1] = octree.root.locateCoordinate(xC+step, yC, zC, level); break;
				case 10: neighbour[a1] = octree.root.locateCoordinate(xC-step, yC+step, zC, level); break;
				case 11: neighbour[a1] = octree.root.locateCoordinate(xC, yC+step, zC, level); break;
				case 12: neighbour[a1] = octree.root.locateCoordinate(xC+step, yC+step, zC, level); break;
				case 13: neighbour[a1] = octree.root.locateCoordinate(xC, yC-step, zC+step, level); break;
				case 14: neighbour[a1] = octree.root.locateCoordinate(xC-step, yC, zC+step, level); break;
				case 15: neighbour[a1] = octree.root.locateCoordinate(xC, yC, zC+step, level); break;
				case 16: neighbour[a1] = octree.root.locateCoordinate(xC+step, yC, zC+step, level); break;
				case 17: neighbour[a1] = octree.root.locateCoordinate(xC, yC+step, zC+step, level); }
			} else {														// the case of aspect pointing to a sibling
				int nbr = nbrs1aspect[i++];
				neighbour[a1] = parent.octant[nbr];
				if (neighbour[a1]==null) neighbour[a1] = parent;			// if sibling didn't exist, return octant's own parent
			}
		}
		return neighbour;
	}


	// method creates octant's 6 face neighbours them if they don't exist, the booleans decide which ones that are needed
	public void generateFaceNeighbours(FEM1Octree octree, boolean xM, boolean xP, boolean yM, boolean yP, boolean zM, boolean zP) {
		byte[] retrace = new byte[octree.topLevel];
		switch (status & 7) {	
		case OCT_MMM:
			symmetricDistantNeighbour(OCT_PMM, 1, retrace, xM);
			if (xP && parent.octant[OCT_PMM]==null) parent.octant[OCT_PMM] = new FEM1Octant(parent, OCT_PMM, false, true);
			symmetricDistantNeighbour(OCT_MPM, 2, retrace, yM);
			if (yP && parent.octant[OCT_MPM]==null) parent.octant[OCT_MPM] = new FEM1Octant(parent, OCT_MPM, false, true);
			symmetricDistantNeighbour(OCT_MMP, 4, retrace, zM);
			if (zP && parent.octant[OCT_MMP]==null) parent.octant[OCT_MMP] = new FEM1Octant(parent, OCT_MMP, false, true); break;
		case OCT_PMM:
			symmetricDistantNeighbour(OCT_MMM, 1, retrace, xP);
			if (xM && parent.octant[OCT_MMM]==null) parent.octant[OCT_MMM] = new FEM1Octant(parent, OCT_MMM, false, true);
			symmetricDistantNeighbour(OCT_PPM, 2, retrace, yM);
			if (yP && parent.octant[OCT_PPM]==null) parent.octant[OCT_PPM] = new FEM1Octant(parent, OCT_PPM, false, true);
			symmetricDistantNeighbour(OCT_PMP, 4, retrace, zM);
			if (zP && parent.octant[OCT_PMP]==null) parent.octant[OCT_PMP] = new FEM1Octant(parent, OCT_PMP, false, true); break;
		case OCT_MPM:
			symmetricDistantNeighbour(OCT_PPM, 1, retrace, xM);
			if (xP && parent.octant[OCT_PPM]==null) parent.octant[OCT_PPM] = new FEM1Octant(parent, OCT_PPM, false, true);
			symmetricDistantNeighbour(OCT_MMM, 2, retrace, yP);
			if (yM && parent.octant[OCT_MMM]==null) parent.octant[OCT_MMM] = new FEM1Octant(parent, OCT_MMM, false, true);
			symmetricDistantNeighbour(OCT_MPP, 4, retrace, zM);
			if (zP && parent.octant[OCT_MPP]==null) parent.octant[OCT_MPP] = new FEM1Octant(parent, OCT_MPP, false, true); break;
		case OCT_PPM:
			symmetricDistantNeighbour(OCT_MPM, 1, retrace, xP);
			if (xM && parent.octant[OCT_MPM]==null) parent.octant[OCT_MPM] = new FEM1Octant(parent, OCT_MPM, false, true);
			symmetricDistantNeighbour(OCT_PMM, 2, retrace, yP);
			if (yM && parent.octant[OCT_PMM]==null) parent.octant[OCT_PMM] = new FEM1Octant(parent, OCT_PMM, false, true);
			symmetricDistantNeighbour(OCT_PPP, 4, retrace, zM);
			if (zP && parent.octant[OCT_PPP]==null) parent.octant[OCT_PPP] = new FEM1Octant(parent, OCT_PPP, false, true); break;
		case OCT_MMP:
			symmetricDistantNeighbour(OCT_PMP, 1, retrace, xM);
			if (xP && parent.octant[OCT_PMP]==null) parent.octant[OCT_PMP] = new FEM1Octant(parent, OCT_PMP, false, true);
			symmetricDistantNeighbour(OCT_MPP, 2, retrace, yM);
			if (yP && parent.octant[OCT_MPP]==null) parent.octant[OCT_MPP] = new FEM1Octant(parent, OCT_MPP, false, true);
			symmetricDistantNeighbour(OCT_MMM, 4, retrace, zP);
			if (zP && parent.octant[OCT_MMM]==null) parent.octant[OCT_MMM] = new FEM1Octant(parent, OCT_MMM, false, true); break;
		case OCT_PMP:
			symmetricDistantNeighbour(OCT_MMP, 1, retrace, xP);
			if (xM && parent.octant[OCT_MMP]==null) parent.octant[OCT_MMP] = new FEM1Octant(parent, OCT_MMP, false, true);
			symmetricDistantNeighbour(OCT_PPP, 2, retrace, yM);
			if (yP && parent.octant[OCT_PPP]==null) parent.octant[OCT_PPP] = new FEM1Octant(parent, OCT_PPP, false, true);
			symmetricDistantNeighbour(OCT_PMM, 4, retrace, zP);
			if (zP && parent.octant[OCT_PMM]==null) parent.octant[OCT_PMM] = new FEM1Octant(parent, OCT_PMM, false, true); break;
		case OCT_MPP:
			symmetricDistantNeighbour(OCT_PPP, 1, retrace, xM);
			if (xP && parent.octant[OCT_PPP]==null) parent.octant[OCT_PPP] = new FEM1Octant(parent, OCT_PPP, false, true);
			symmetricDistantNeighbour(OCT_MMP, 2, retrace, yP);
			if (yM && parent.octant[OCT_MMP]==null) parent.octant[OCT_MMP] = new FEM1Octant(parent, OCT_MMP, false, true);
			symmetricDistantNeighbour(OCT_MPM, 4, retrace, zP);
			if (zP && parent.octant[OCT_MPM]==null) parent.octant[OCT_MPM] = new FEM1Octant(parent, OCT_MPM, false, true); break;
		case OCT_PPP:
			symmetricDistantNeighbour(OCT_MPP, 1, retrace, xP);
			if (xM && parent.octant[OCT_MPP]==null) parent.octant[OCT_MPP] = new FEM1Octant(parent, OCT_MPP, false, true);
			symmetricDistantNeighbour(OCT_PMP, 2, retrace, yP);
			if (yM && parent.octant[OCT_PMP]==null) parent.octant[OCT_PMP] = new FEM1Octant(parent, OCT_PMP, false, true);
			symmetricDistantNeighbour(OCT_PPM, 4, retrace, zP);
			if (zP && parent.octant[OCT_PPM]==null) parent.octant[OCT_PPM] = new FEM1Octant(parent, OCT_PPM, false, true);
		}
	}
	
	// this method collects the 6 facet neighbours, generating them if generate=true AND they are requested by the booleans and missing
	// if fillIn=true, existent octants will NOT be overwritten
	// nbrOfs will place the 6 neighbours at arbitrary position within a supplied octantN[] array
	public FEM1Octant[] generateGetFaceNeighbours(FEM1Octree octree, FEM1Octant[] octantN, int nFlags, int nbrOfs, boolean fillIn) {

		if (octantN == null) octantN =  new FEM1Octant[7];		// one extra space to include the central octant of this neighbourhood, if caller needs it
		else if (!fillIn) {
			octantN[nbrOfs++] = null; octantN[nbrOfs++] = null; octantN[nbrOfs++] = null;
			octantN[nbrOfs++] = null; octantN[nbrOfs++] = null; octantN[nbrOfs] = null; nbrOfs -= 5;
		}
		byte[] retrace = new byte[this.level];
		int n0=0, n1=0, n2=0, n3=0, n4=0 ,n5=0;
		short o0=0, o1=0, o2=0;
		switch (status & 7) {	
		case OCT_MMM:	n0=nbrOfs+2; n1=nbrOfs+3; n2=nbrOfs+1; n3=nbrOfs+4; n4=nbrOfs; n5=nbrOfs+5; o0=OCT_PMM; o1=OCT_MPM; o2=OCT_MMP; break;
		case OCT_PMM:	n0=nbrOfs+3; n1=nbrOfs+2; n2=nbrOfs+1; n3=nbrOfs+4; n4=nbrOfs; n5=nbrOfs+5; o0=OCT_MMM; o1=OCT_PPM; o2=OCT_PMP; break;
		case OCT_MPM:	n0=nbrOfs+2; n1=nbrOfs+3; n2=nbrOfs+4; n3=nbrOfs+1; n4=nbrOfs; n5=nbrOfs+5; o0=OCT_PPM; o1=OCT_MMM; o2=OCT_MPP; break;
		case OCT_PPM:	n0=nbrOfs+3; n1=nbrOfs+2; n2=nbrOfs+4; n3=nbrOfs+1; n4=nbrOfs; n5=nbrOfs+5; o0=OCT_MPM; o1=OCT_PMM; o2=OCT_PPP; break;
		case OCT_MMP:	n0=nbrOfs+2; n1=nbrOfs+3; n2=nbrOfs+1; n3=nbrOfs+4; n4=nbrOfs+5; n5=nbrOfs; o0=OCT_PMP; o1=OCT_MPP; o2=OCT_MMM; break;
		case OCT_PMP:	n0=nbrOfs+3; n1=nbrOfs+2; n2=nbrOfs+1; n3=nbrOfs+4; n4=nbrOfs+5; n5=nbrOfs; o0=OCT_MMP; o1=OCT_PPP; o2=OCT_PMM; break;
		case OCT_MPP:	n0=nbrOfs+2; n1=nbrOfs+3; n2=nbrOfs+4; n3=nbrOfs+1; n4=nbrOfs+5; n5=nbrOfs; o0=OCT_PPP; o1=OCT_MMP; o2=OCT_MPM; break;
		case OCT_PPP:	n0=nbrOfs+3; n1=nbrOfs+2; n2=nbrOfs+4; n3=nbrOfs+1; n4=nbrOfs+5; n5=nbrOfs; o0=OCT_MPP; o1=OCT_PMP; o2=OCT_PPM;
		}
		FEM1Octant[] octantP = parent.octant;
		if (n0 == nbrOfs+2) {					// handles -x / +x neighbours
			if ((nFlags&16)!=0) octantN[n0] = symmetricDistantNeighbour(o0, 1, retrace, (nFlags&32)!=0);
			if ((nFlags&64)!=0) {
                if (octantP[o0]==null) { if ((nFlags&128)!=0) { octantN[n1] = octantP[o0] = new FEM1Octant(parent, o0, true, true); }}
				else octantN[n1] = octantP[o0]; }
		} else {
			if ((nFlags&64)!=0) octantN[n0] = symmetricDistantNeighbour(o0, 1, retrace, (nFlags&128)!=0);
			if ((nFlags&16)!=0) {
                if (octantP[o0]==null) { if ((nFlags&32)!=0) { octantN[n1] = octantP[o0] = new FEM1Octant(parent, o0, true, true); }}
				else octantN[n1] = octantP[o0]; }
		}
		if (n2 == nbrOfs+1) {					// handles -y / +y neighbours
			if ((nFlags&4)!=0) octantN[n2] = symmetricDistantNeighbour(o1, 2, retrace, (nFlags&8)!=0);
			if ((nFlags&256)!=0) {
                if (octantP[o1]==null) { if ((nFlags&512)!=0) { octantN[n3] = octantP[o1] = new FEM1Octant(parent, o1, true, true); }}
				else octantN[n3] = octantP[o1]; }
		} else {
			if ((nFlags&256)!=0) octantN[n2] = symmetricDistantNeighbour(o1, 2, retrace, (nFlags&512)!=0);
			if ((nFlags&4)!=0) {
                if (octantP[o1]==null) { if ((nFlags&8)!=0) { octantN[n3] = octantP[o1] = new FEM1Octant(parent, o1, true, true); }}
				else octantN[n3] = octantP[o1]; }
		}
		if (n4 == nbrOfs) {					// handles -z / +z neighbours
			if ((nFlags&1)!=0) octantN[n4] = symmetricDistantNeighbour(o2, 4, retrace, (nFlags&2)!=0);
			if ((nFlags&1024)!=0) {
                if (octantP[o2]==null) { if ((nFlags&2048)!=0) { octantN[n5] = octantP[o2] = new FEM1Octant(parent, o2, true, true); }
				else octantN[n5] = octantP[o2]; }}
		} else {
			if ((nFlags&1024)!=0) octantN[n4] = symmetricDistantNeighbour(o2, 4, retrace, (nFlags&2048)!=0);
			if ((nFlags&1)!=0) {
                if (octantP[o2]==null) { if ((nFlags&2)!=0) { octantN[n5] = octantP[o2] = new FEM1Octant(parent, o2, true, true); }}
				else octantN[n5] = octantP[o2]; }
		}
		return octantN;
	}

	
	
	
	// this is the order the child qualification will be tested for every octant, geometrically checking around octant
	final static byte[] s2to1flag = {
		FLG_XMM,FLG_MXM,FLG_XXM,FLG_PXM,FLG_XPM,FLG_MMX,FLG_XMX,FLG_PMX,FLG_MXX,FLG_PXX,FLG_MPX,FLG_XPX,FLG_PPX,FLG_XMP,FLG_MXP,FLG_XXP,FLG_PXP,FLG_XPP };
	// neghbour indexes & aspect vectors for neighbour inclusion for child, 1024+x signify inclusion of child's siblings, positives are
	// properly aligned suboctants taken from neighbours, each child receives it's own neighbour list
	final static byte[] s2to1aspect = {
		0,OCT_MPP, 1,OCT_PMP, 2,OCT_MMP, 2,OCT_PMP, 2,OCT_MPP, 5,OCT_PPM, 6,OCT_MPM, 6,OCT_PPM, 8,OCT_PMM,
		-1,OCT_PMM, 8,OCT_PPM, -1,OCT_MPM, -1,OCT_PPM, 6,OCT_MPP, 8,OCT_PMP, -1,OCT_MMP, -1,OCT_PMP, -1,OCT_MPP,		
/*36*/		0,OCT_PPP, 2,OCT_MMP, 2,OCT_PMP, 3,OCT_MMP, 2,OCT_PPP, 6,OCT_MPM, 6,OCT_PPM, 7,OCT_MPM, -1,OCT_MMM,
			9,OCT_MMM, -1,OCT_MPM, -1,OCT_PPM, 9,OCT_MPM, 6,OCT_PPP, -1,OCT_MMP, -1,OCT_PMP, 9,OCT_MMP, -1,OCT_PPP,	
/*72*/	2,OCT_MMP, 1,OCT_PPP, 2,OCT_MPP, 2,OCT_PPP, 4,OCT_MMP, 8,OCT_PMM, -1,OCT_MMM, -1,OCT_PMM, 8,OCT_PPM,
		-1,OCT_PPM, 10,OCT_PMM, 11,OCT_MMM, 11,OCT_PMM, -1,OCT_MMP, 8,OCT_PPP, -1,OCT_MPP, -1,OCT_PPP, 11,OCT_MMP,
/*108*/		2,OCT_PMP, 2,OCT_MPP, 2,OCT_PPP, 3,OCT_MPP, 4,OCT_PMP, -1,OCT_MMM, -1,OCT_PMM, 9,OCT_MMM, -1,OCT_MPM,
			9,OCT_MPM, 11,OCT_MMM, 11,OCT_PMM, 12,OCT_MMM, -1,OCT_PMP, -1,OCT_MPP, -1,OCT_PPP, 9,OCT_MPP, 11,OCT_PMP,	
/*144*/	6,OCT_MPM, 8,OCT_PMM, -1,OCT_MMM, -1,OCT_PMM, -1,OCT_MPM, 5,OCT_PPP, 6,OCT_MPP, 6,OCT_PPP, 8,OCT_PMP,
		-1,OCT_PMP, 8,OCT_PPP, -1,OCT_MPP, -1,OCT_PPP, 13,OCT_MPM, 14,OCT_PMM, 15,OCT_MMM, 15,OCT_PMM, 15,OCT_MPM,	
/*180*/		6,OCT_PPM, -1,OCT_MMM, -1,OCT_PMM, 9,OCT_MMM, -1,OCT_PPM, 6,OCT_MPP, 6,OCT_PPP, 7,OCT_MPP, -1,OCT_MMP,
			9,OCT_MMP, -1,OCT_MPP, -1,OCT_PPP, 9,OCT_MPP, 13,OCT_PPM, 15,OCT_MMM, 15,OCT_PMM, 16,OCT_MMM, 15,OCT_PPM,	
/*216*/	-1,OCT_MMM, 8,OCT_PPM, -1,OCT_MPM, -1,OCT_PPM, 11,OCT_MMM, 8,OCT_PMP, -1,OCT_MMP, -1,OCT_PMP, 8,OCT_PPP,
		-1,OCT_PPP, 10,OCT_PMP, 11,OCT_MMP, 11,OCT_PMP, 15,OCT_MMM, 14,OCT_PPM, 15,OCT_MPM, 15,OCT_PPM, 17,OCT_MMM,	
/*252*/		-1,OCT_PMM, -1,OCT_MPM, -1,OCT_PPM, 9,OCT_MPM, 11,OCT_PMM, -1,OCT_MMP, -1,OCT_PMP, 9,OCT_MMP, -1,OCT_MPP,
			9,OCT_MPP, 11,OCT_MMP, 11,OCT_PMP, 12,OCT_MMP, 15,OCT_PMM, 15,OCT_MPM, 15,OCT_PPM, 16,OCT_MPM, 17,OCT_PMM };
	
	// this is a dummy array used by split2to1recursor() method to store neighbour values that don't matter anymore
	private static FEM1Octant[] s2to1neighbourLeaf = new FEM1Octant[18];
	
	// method applies recursive Weak Condition 2:1 refinement on an octree, every branch tests against neighbour families and
	// brings up-level the neighbour octants of their children for a recursive child-check
	public boolean split2to1recursor(FEM1Octree octree, FEM1Octant[] neighbour) {
		boolean splitMade = false;
		int octantIt = iterator;
		boolean stopRecursion = (level < octree.topLevel - 2) ? false : true;

		// every child of this octant is first enforcing the neighbours to conform to 2:1 split, then collects the resulting neighbourhood for recursion
		for (int o1 = 0, oEnd = iterator >> 24; o1 < oEnd; o1++, octantIt >>= 3) {
			int o = octantIt & 7;
			FEM1Octant octantO = octant[o];
			if (octantO.octant == null) continue;						// skip leaf child (but it will be included as neighbour by siblings)

			int bitsC = octantO.status >> 3;
			// create neighbourhood if 4:1 relation holds at next level, otherwise let method write to a dummy neighbourhood array, it won't be used
			FEM1Octant[] neighbourC = stopRecursion ? s2to1neighbourLeaf : new FEM1Octant[18];
			
			for (int i = o * 36, a = 0, iEnd = i + 36; i < iEnd; a++) {	// collect or create child level's neighbourhood for this child
				if ((bitsC & s2to1flag[a]) != 0) {						// if octant has children vs tested 2:1 aspect
					short nbr = s2to1aspect[i++];
					if (nbr >= 0) {
						if (neighbour[nbr] == null) { i++; continue; }	// there was no neighbour to test
						short oN = s2to1aspect[i++];					// get the 2:1-aspect child of neighbour family
						FEM1Octant octantN = neighbour[nbr];
						if (octantN.octant == null) {					// case of a 2:1 leaf
							octantN.octant = new FEM1Octant[8];
							neighbourC[a] = octantN.octant[oN] = new FEM1Octant(octantN, oN, false, true);
							if (octree.latticeTree) {	neighbourC[a].int_ext = INTERNAL; 
														octree.levelSums[level + 2]++; }
							splitMade = true;
						} else if (octantN.octant[oN] == null) {		// case of child not existing
							neighbourC[a] = octantN.octant[oN] = new FEM1Octant(octantN, oN, false, true);
							if (octree.latticeTree) {	neighbourC[a].int_ext = INTERNAL; 
														octree.levelSums[level + 2]++; }
							splitMade = true;
						} else {
							neighbourC[a] = octantN = octantN.octant[oN];		// case of child having proper 2:1 neighbour
//							if (((octantN.status>>3) & s2to1flag[17-a]) != 0)	// if it's a true 2:1 relationship (neighbour has no children facing aspect)
//								octantN.int_ext |= 1 << (17 - a + 8);			// flag neighbour's counteraspect as true 2:1 relationship
						}
					} else {
						nbr = s2to1aspect[i++];
						if (octant[nbr] == null) {						// case of a 2:1 parent-sibling not existing
							octant[nbr] = neighbourC[a] = new FEM1Octant(this, nbr, false, true);
							if (octree.latticeTree) {	octant[nbr].int_ext = INTERNAL;
														octree.levelSums[level + 2]++; }
							splitMade = true;
						}
					}
				} else i += 2;
			}
			if (!stopRecursion)													// if a 4:1 relationship still can be tested for next level
				splitMade |= octantO.split2to1recursor(octree, neighbourC);		// adjust child's family & neighbourhood 2:1 recursively
		}
		return splitMade;														// report downwards if 2:1 imbalances were found
	}
	
	
	// the multitasked method for 2:1 Weak Condition refinement calls new threads on every recursion if a thread is available,
	// except the last recursion which is continued by current thread itself, generated offshoots (which WILL have redundant data)
	// are simply added to a concurrent Queue, and later merged back into the tree in serial fashion
	// thus threads are spawning offshoots independently each in it's own octant branch
	public void split2to1recursor(FEM1Octree octree, FEM1Octant[] neighbour, Queue<FEM1Octant> queue2to1) {
		
		boolean stopRecursion = (level < octree.topLevel - 2) ? false : true;
		// every child of this octant is first enforcing the neighbours to conform to 2:1 split, then collects the resulting neighbourhood for recursion
		int octantIt = iterator;
		for (int o1 = 0, oEnd = iterator >> 24; o1 < oEnd; o1++, octantIt >>= 3) {
			int o = octantIt & 7;
			FEM1Octant octantO = octant[o];
			if (octantO.octant == null) continue;						// skip leaf or offshoot child (will be included as neighbour by siblings)						

			int bitsC = octantO.status >> 3;
			// create neighbourhood if 4:1 relation holds at next level, otherwise let method write to a dummy neighbourhood array, it won't be used
			FEM1Octant[] neighbourC = stopRecursion ? s2to1neighbourLeaf : new FEM1Octant[18];
			
			for (int i = o * 36, a = 0, iEnd = i + 36; i < iEnd; a++) {	// collect or create child level's neighbourhood for this child
				if ((bitsC & s2to1flag[a]) != 0) {						// if test octant has children versus tested 2:1 aspect
					short nbr = s2to1aspect[i++];
					if (nbr >= 0) {
						if (neighbour[nbr] == null) { i++; continue; }	// there was no neighbour to test
						short oN = s2to1aspect[i++];					// get the 2:1-aspect child of neighbour family
						FEM1Octant octantN = neighbour[nbr];
						
						if (octantN.octant == null) 					// case of a leaf neighbour
							octantN.octant = new FEM1Octant[8];
						if (octantN.octant[oN] == null) {
							// every threaded construction is "in progress" and must not touch parent data, such as the iterator
							FEM1Octant octant1 = neighbourC[a] = new FEM1Octant(octantN, oN, false, false);	
							if (octree.latticeTree) {
								octant1.int_ext = INTERNAL;
								octree.levelSumsT[level+2].getAndIncrement(); }
							// if octant is offshooting (in detached construction by the thread), it's safe to link it's child to it
							if (octantN.offshoot()) {
									octant1.status |= OFFSHOOT;			// child is also an offshoot
									octantN.octant[oN] = octant1;
									octantN.iteratorAdd(oN);
							// if this is an octant spawned from an existing tree octant, add it as offshoot to the concurrency-safe queue
							} else {
								octant1.status |= OFFSHOOT;
								queue2to1.add(octant1);
								FEM1Octree.offshoots.incrementAndGet();
							}
							octree.splitMade = true;							
						} else {
							neighbourC[a] = octantN = octantN.octant[oN];		// case of child having proper 2:1 neighbour
//							if (((octantN.status>>3) & s2to1flag[17-a]) != 0)	// if it's a true 2:1 relationship (neighbour has no children facing aspect)
//								octantN.int_ext |= 1 << (17 - a + 8);			// flag neighbour's counteraspect as true 2:1 relationship
						}
						
					} else {
						nbr = s2to1aspect[i++];
						if (octant[nbr] == null) {								// case of a 2:1 parent-sibling not existing
							// for in-family work the thread can add directly to local branch
							octant[nbr] = neighbourC[a] = new FEM1Octant(this, nbr, false, true);	
							if (octree.latticeTree) {
								octant[nbr].int_ext = INTERNAL;
								octree.levelSumsT[level+2].getAndIncrement(); }
							octree.splitMade = true;
						}
					}
				} else i += 2;
			}
			
			// adjust child's family & neighbourhood 2:1 recursively
			if (!stopRecursion)	{												// if a 4:1 relationship still can be tested for next level
				if (oEnd - o1 <= 1 || FEM1.freeTasks.get() <= 0) {				// if this is the last call, continue recursion in current thread
					octantO.status |= IN_PROGRESS;
					octantO.split2to1recursor(octree, neighbourC, queue2to1);
				} else {														// otherwise, portion out octants to thread queue
					FEM1.freeTasks.decrementAndGet();
					FEM1.getExecutor(0).execute(new Runnable() {
						@Override public void run() {
							octantO.status |= IN_PROGRESS;
							octantO.split2to1recursor(octree, neighbourC, queue2to1);
							FEM1.freeTasks.incrementAndGet();
					}});
				}
			} else FEM1Octree.processed.addAndGet(octantO.branches);
		}
		FEM1Octree.processed.incrementAndGet();
		//octree.levelSumsT[level].getAndIncrement();
		finish();
	}

	// method mixes/weaves together two branches that share space and have differing offshoots
	static void weaveBranch(FEM1Octant parent, FEM1Octant branch) {
		int o = branch.status & 7;
		if (parent.octant == null) {
			parent.octant = new FEM1Octant[8]; parent.octant[o] = branch; parent.iteratorAdd(o);
			return; }
		if (parent.octant[o] == null) {										// base case: branch directly insertable
			parent.octant[o] = branch; parent.iteratorAdd(o);
			return; }
		if (branch.octant == null) return;									// base case: that octant exists, and inserted  branch is a leaf
		for (FEM1Octant octant : branch.octant)
			if (octant != null) weaveBranch(parent.octant[o], octant);		// this branch existed, try inserting it's subbranches recursively
	}
	
	public void split2to1(FEM1Octree octree, boolean multitask) {
		if (level >= octree.topLevel - 1) return;							// base case: octree had only 2 levels and is already a 2:1 leaf
		FEM1Octant[] neighbour = new FEM1Octant[18];
		
		if (multitask) {
			Queue<FEM1Octant> queue2to1 = new ConcurrentLinkedQueue<FEM1Octant>();
			do {															// propagate 2:1 splits until no split can be made
				octree.splitMade = false;
				FEM1Octree.processed.set(0);
				status |= IN_PROGRESS;
				split2to1recursor(octree, neighbour, queue2to1);
				while (octree.root.branches > FEM1Octree.processed.get());	// wait until every branch is accounted for
				if (queue2to1.size() == 0) {								// base case: all the splits found were only family-local
					process(octree, OCTREE_SUM_BRANCHES|OCTREE_FIX_ITERATORS);
					break;
				}
				FEM1Octant branch2to1;
				while ((branch2to1 = queue2to1.poll()) != null) {			// pop offshoots from head of queue
					int o = branch2to1.status & 7;
					if (branch2to1.parent.octant == null)
						branch2to1.parent.octant = new FEM1Octant[8];
					if (branch2to1.parent.octant[o] != null) {
						if (branch2to1.octant != null)						// if a branch offshoot was spawned, weave it into parent branch
							FEM1Octant.weaveBranch(branch2to1.parent, branch2to1);
					} else { branch2to1.parent.octant[o] = branch2to1; branch2to1.parent.iteratorAdd(o); }				
				}
				process(octree, OCTREE_SUM_BRANCHES|OCTREE_FIX_ITERATORS);
			} while (octree.splitMade);										// propagate 2:1 splits until no split can be made (so-called "ripples")
		} else {
			while (split2to1recursor(octree, neighbour))					// propagate 2:1 splits until no split can be made
				process(octree, OCTREE_SUM_BRANCHES|OCTREE_FIX_ITERATORS);
		}
	}
	
	
	// method keeps only children that are needed to satisfy 2:1 Weak Condition towards face&edge neighbourhood
	// note: method used to cheaply generate internal gradation levels without full-fledged 2:1 balancing
	int cull2to1IST(FEM1Octree octree, boolean recurse) {
		int octantIt = iterator, oEnd = iterator >> 24, oldLeaves = leaves;
		leaves = iterator = 0; status &= 0xF807;								// reset leaves count & iterator & status, will be recounted	
		FEM1Octant[] newFamily = new FEM1Octant[8]; int redeemed = 0;

		for (int o1 = 0; o1 < oEnd; o1++, octantIt >>= 3) {
			int o = octantIt & 7;
			FEM1Octant octantO = octant[o];
			//if (octantO == null) continue;									// no child in that position, nothing to check against
			if (recurse && !octantO.culled() && level < octree.maxLevel - 3) {	// on request, recursively cull children before culling parent
				octantO.cull2to1IST(octree, recurse);
				//leaves += octantO.leaves;
			}
			FEM1Octant[] neighbourC = octantO.getFaceEdgeNeighbours(octree, null, 0, 0x3FFFF);
			boolean culled = true;
			
			for (int i=o*36, a=0, iEnd = i+36; i < iEnd; i+=2, a++) {			// collect or create child level's neighbourhood for this child
				short nbr = s2to1aspect[i];
				if (nbr >= 0) {													// we're only interested in neighbour families
					if (neighbourC[a] == null)
						continue; 												// there was no neighbour to test
					FEM1Octant octantN = neighbourC[a];
					if (octantO.level == octantN.level &&						// if neighbour is of same level (check only against equal children)
						(octantN.status>>3&s2to1flag[17-a])!=0) {				// and if neighbour child has a grandchild towards this octant's child
						newFamily[o] = octantO; redeemed++;						// then keep this child
						iteratorAdd(o);	culled = false; break; }				// update all octant data	
				}
			}
			if (culled) octantO.nodes = -1;
		}
		status |= CULLED;
		octree.levelSums[level+1] -= (oldLeaves - redeemed);
		octant = redeemed == 0 ? null : newFamily;
		return redeemed;
	}
	

	final static int FLG_BFC_MMM = 1|1<<8|1<<9|1<<10|1<<13|1<<14|1<<16|1<<17,		FLG_BFC_PMM = 1<<8|1<<1|1<<10|1<<11|1<<14|1<<15|1<<17|1<<18;
	final static int FLG_BFC_MPM = 1<<9|1<<10|1<<2|1<<12|1<<16|1<<17|1<<19|1<<20,	FLG_BFC_PPM = 1<<10|1<<11|1<<12|1<<3|1<<17|1<<18|1<<20|1<<21;
	final static int FLG_BFC_MMP = 1<<13|1<<14|1<<16|1<<17|1<<4|1<<22|1<<23|1<<24,	FLG_BFC_PMP = 1<<14|1<<15|1<<17|1<<18|1<<22|1<<5|1<<24|1<<25;
	final static int FLG_BFC_MPP = 1<<16|1<<17|1<<19|1<<20|1<<23|1<<24|1<<6|1<<26,	FLG_BFC_PPP = 1<<17|1<<18|1<<20|1<<21|1<<24|1<<25|1<<26|1<<7;
	// array maps child's int_ext flags of upper or lower 4 bits to parent's 27-point int_ext flags
	final static int fMap_bFC[] = { 0, 1, 1<<8, 1|1<<8, 1<<9, 1|1<<9, 1<<8|1<<9, 1|1<<8|1<<9,
		1<<10, 1|1<<10, 1<<8|1<<10, 1|1<<8|1<<10, 1<<9|1<<10, 1|1<<9|1<<10, 1<<8|1<<9|1<<10, 1|1<<8|1<<9|1<<10,
		0, 1<<8, 1<<1, 1<<8|1<<1, 1<<10, 1<<8|1<<10, 1<<1|1<<10, 1<<8|1<<1|1<<10,															// +16
		1<<11, 1<<8|1<<11, 1<<1|1<<11, 1<<8|1<<1|1<<11, 1<<10|1<<11, 1<<8|1<<10|1<<11, 1<<1|1<<10|1<<11, 1<<8|1<<1|1<<10|1<<11,	
		0, 1<<9, 1<<10, 1<<9|1<<10, 1<<2, 1<<9|1<<2, 1<<10|1<<2, 1<<9|1<<10|1<<2,															// +32
		1<<12, 1<<9|1<<12, 1<<10|1<<12, 1<<9|1<<10|1<<12, 1<<2|1<<12, 1<<9|1<<2|1<<12, 1<<10|1<<2|1<<12, 1<<9|1<<10|1<<2|1<<12,
		0, 1<<10, 1<<11, 1<<10|1<<11, 1<<12, 1<<10|1<<12, 1<<11|1<<12, 1<<10|1<<11|1<<12,													// +48
		1<<3, 1<<10|1<<3, 1<<11|1<<3, 1<<10|1<<11|1<<3, 1<<12|1<<3, 1<<10|1<<12|1<<3, 1<<11|1<<12|1<<3, 1<<10|1<<11|1<<12|1<<3,
		
		0, 1<<13, 1<<14, 1<<13|1<<14, 1<<16, 1<<13|1<<16, 1<<14|1<<16, 1<<13|1<<14|1<<16,													// +64
		1<<17, 1<<13|1<<17, 1<<14|1<<17, 1<<13|1<<14|1<<17, 1<<16|1<<17, 1<<13|1<<16|1<<17, 1<<14|1<<16|1<<17, 1<<13|1<<14|1<<16|1<<17,
		0, 1<<14, 1<<15, 1<<14|1<<15, 1<<17, 1<<14|1<<17, 1<<15|1<<17, 1<<14|1<<15|1<<17,													// +80
		1<<18, 1<<14|1<<18, 1<<15|1<<18, 1<<14|1<<15|1<<18, 1<<17|1<<18, 1<<14|1<<17|1<<18, 1<<15|1<<17|1<<18, 1<<14|1<<15|1<<17|1<<18,
		0, 1<<16, 1<<17, 1<<16|1<<17, 1<<19, 1<<16|1<<19, 1<<17|1<<19, 1<<16|1<<17|1<<19,													// +96
		1<<20, 1<<16|1<<20, 1<<17|1<<20, 1<<16|1<<17|1<<20, 1<<19|1<<20, 1<<16|1<<19|1<<20, 1<<17|1<<19|1<<20, 1<<16|1<<17|1<<19|1<<20,
		0, 1<<17, 1<<18, 1<<17|1<<18, 1<<20, 1<<17|1<<20, 1<<18|1<<20, 1<<17|1<<18|1<<20,													// +112
		1<<21, 1<<17|1<<21, 1<<18|1<<21, 1<<17|1<<18|1<<21, 1<<20|1<<21, 1<<17|1<<20|1<<21, 1<<18|1<<20|1<<21, 1<<17|1<<18|1<<20|1<<21,

		0, 1<<4, 1<<22, 1<<4|1<<22, 1<<23, 1<<4|1<<23, 1<<22|1<<23, 1<<4|1<<22|1<<23,														// +128
		1<<24, 1<<4|1<<24, 1<<22|1<<24, 1<<4|1<<22|1<<24, 1<<23|1<<24, 1<<4|1<<23|1<<24, 1<<22|1<<23|1<<24, 1<<4|1<<22|1<<23|1<<24,
		0, 1<<22, 1<<5, 1<<22|1<<5, 1<<24, 1<<22|1<<24, 1<<5|1<<24, 1<<22|1<<5|1<<24,														// +144
		1<<25, 1<<22|1<<25, 1<<5|1<<25, 1<<22|1<<5|1<<25, 1<<24|1<<25, 1<<22|1<<24|1<<25, 1<<5|1<<24|1<<25, 1<<22|1<<5|1<<24|1<<25,
		0, 1<<23, 1<<24, 1<<23|1<<24, 1<<6, 1<<23|1<<6, 1<<24|1<<6, 1<<23|1<<24|1<<6,														// +160
		1<<26, 1<<23|1<<26, 1<<24|1<<26, 1<<23|1<<24|1<<26, 1<<6|1<<26, 1<<23|1<<6|1<<26, 1<<24|1<<6|1<<26, 1<<23|1<<24|1<<6|1<<26,
		0, 1<<24, 1<<25, 1<<24|1<<25, 1<<26, 1<<24|1<<26, 1<<25|1<<26, 1<<24|1<<25|1<<26,													// +176
		1<<7, 1<<24|1<<7, 1<<25|1<<7, 1<<24|1<<25|1<<7, 1<<26|1<<7, 1<<24|1<<26|1<<7, 1<<25|1<<26|1<<7, 1<<24|1<<25|1<<26|1<<7 };

	// method resolves internal/external enumerated flags of "boundary" parameter for parent octant from children's boundaries
	void int_extFromChildren(boolean inclusive) {
		if (octant == null) return;		// cannot resolve anything if there are no children
		if (!inclusive) int_ext = 0;
		if (octant[0]!=null) {	if (octant[0].internal()) int_ext |= FLG_BFC_MMM;
								else { int inEx=octant[0].int_ext&0xFF; int_ext |= fMap_bFC[inEx&15]|fMap_bFC[64+(inEx>>4)]; }}
		if (octant[1]!=null) {	if (octant[1].internal()) int_ext |= FLG_BFC_PMM;
								else { int inEx=octant[1].int_ext&0xFF; int_ext |= fMap_bFC[16+(inEx&15)]|fMap_bFC[80+(inEx>>4)]; }}
		if (octant[2]!=null) {	if (octant[2].internal()) int_ext |= FLG_BFC_MPM;
								else { int inEx=octant[2].int_ext&0xFF; int_ext |= fMap_bFC[32+(inEx&15)]|fMap_bFC[96+(inEx>>4)]; }}
		if (octant[3]!=null) {	if (octant[3].internal()) int_ext |= FLG_BFC_PPM;
								else { int inEx=octant[3].int_ext&0xFF; int_ext |= fMap_bFC[48+(inEx&15)]|fMap_bFC[112+(inEx>>4)]; }}
		if (octant[4]!=null) {	if (octant[4].internal()) int_ext |= FLG_BFC_MMP;
								else { int inEx=octant[4].int_ext&0xFF; int_ext |= fMap_bFC[64+(inEx&15)]|fMap_bFC[128+(inEx>>4)]; }}
		if (octant[5]!=null) {	if (octant[5].internal()) int_ext |= FLG_BFC_PMP;
								else { int inEx=octant[5].int_ext&0xFF; int_ext |= fMap_bFC[80+(inEx&15)]|fMap_bFC[144+(inEx>>4)]; }}
		if (octant[6]!=null) {	if (octant[6].internal()) int_ext |= FLG_BFC_MPP;
								else { int inEx=octant[6].int_ext&0xFF; int_ext |= fMap_bFC[96+(inEx&15)]|fMap_bFC[160+(inEx>>4)]; }}
		if (octant[7]!=null) {	if (octant[7].internal()) int_ext |= FLG_BFC_PPP;
								else { int inEx=octant[7].int_ext&0xFF; int_ext |= fMap_bFC[112+(inEx&15)]|fMap_bFC[176+(inEx>>4)]; }}
	}
	
	// method resolves child's internal/external boundary corner flags from parent, if inclusive=true, it will OR into int_ext rather than replace it
	// note: since only corner points can be inherited, the rest of the existent boundary bits are kept untouched
	void int_extFromParent(boolean inclusive) {
		if (parent == null) return;
		if (!inclusive) int_ext &= 0xFFFFFF00;
		int inExP = parent.int_ext;
		switch (status & 7) {
		case OCT_MMM: int_ext |= inExP&1|inExP>>7&2|inExP>>7&4|inExP>>7&8|inExP>>9&16|inExP>>9&32|inExP>>10&64|inExP>>10&128; break;
		case OCT_PMM: int_ext |= inExP>>8&1|inExP&2|inExP>>8&4|inExP>>8&8|inExP>>10&16|inExP>>10&32|inExP>>11&64|inExP>>11&128; break;
		case OCT_MPM: int_ext |= inExP>>9&1|inExP>>9&2|inExP&4|inExP>>9&8|inExP>>12&16|inExP>>12&32|inExP>>13&64|inExP>>13&128; break;
		case OCT_PPM: int_ext |= inExP>>10&1|inExP>>10&2|inExP>>10&4|inExP&8|inExP>>13&16|inExP>>13&32|inExP>>14&64|inExP>>14&128; break;
		case OCT_MMP: int_ext |= inExP>>13&1|inExP>>13&2|inExP>>14&4|inExP>>14&8|inExP&16|inExP>>17&32|inExP>>17&64|inExP>>17&128; break;
		case OCT_PMP: int_ext |= inExP>>14&1|inExP>>14&2|inExP>>15&4|inExP>>15&8|inExP>>18&16|inExP&32|inExP>>18&64|inExP>>18&128; break;
		case OCT_MPP: int_ext |= inExP>>16&1|inExP>>16&2|inExP>>17&4|inExP>>17&8|inExP>>19&16|inExP>>19&32|inExP&64|inExP>>19&128; break;
		case OCT_PPP: int_ext |= inExP>>17&1|inExP>>17&2|inExP>>18&4|inExP>>18&8|inExP>>20&16|inExP>>20&32|inExP>>20&64|inExP&128;
		}
		// DEBUG: for IST method, creation happens for internal octants of (maxLevel-1) gradation level, so centroid will clearly be internal
		if ((int_ext&0xFF)==0xFF) int_ext |= 1<<17;
	}
	
	
	// method will use nodesX as gradation enumerator & flagger for IST-internal suboctants, if all those flags become set
	// then this octant will become an IST-internal octant itself (which simplifies filtering out IST octants)
	// this process must be top-down from the leaves, and is an EXCLUSIVE process as one progresses down the gradations
	boolean internOfISTFromChildren() {
		nodesX = 0;
		if (octant == null) return false;		// cannot resolve anything if there are no children
		if (octant[0]!=null && (octant[0].nodesX&0xFFFF)>0) nodesX |= 1<<16;
		if (octant[1]!=null && (octant[1].nodesX&0xFFFF)>0) nodesX |= 2<<16;
		if (octant[2]!=null && (octant[2].nodesX&0xFFFF)>0) nodesX |= 4<<16;
		if (octant[3]!=null && (octant[3].nodesX&0xFFFF)>0) nodesX |= 8<<16;
		if (octant[4]!=null && (octant[4].nodesX&0xFFFF)>0) nodesX |= 16<<16;
		if (octant[5]!=null && (octant[5].nodesX&0xFFFF)>0) nodesX |= 32<<16;
		if (octant[6]!=null && (octant[6].nodesX&0xFFFF)>0) nodesX |= 64<<16;
		if (octant[7]!=null && (octant[7].nodesX&0xFFFF)>0) nodesX |= 128<<16;
		// if all children were IST-internal, make octant IST-internal by giving it a gradation number
		if ((nodesX & 0xFF0000) == 0xFF0000) { nodesX |= level; return true; }
		return false;
	}

	// method needs volume internal/external data flags from boundary lattice
	boolean internOfISTFromBLatticeData(short fMMM, short fPMM, short fMPM, short fPPM, short fMMP, short fPMP, short fMPP, short fPPP) {
		nodesX = 0;
		if (fMMM==0x1FF) nodesX |= 1<<16;	if (fPMM==0x1FF) nodesX |= 2<<16;	if (fMPM==0x1FF) nodesX |= 4<<16;	if (fPPM==0x1FF) nodesX |= 8<<16;
		if (fMMP==0x1FF) nodesX |= 16<<16;	if (fPMP==0x1FF) nodesX |= 32<<16;	if (fMPP==0x1FF) nodesX |= 64<<16;	if (fPPP==0x1FF) nodesX |= 128<<16;
		// if all children were IST-internal, make octant IST-internal by giving it a gradation number
		if ((nodesX & 0xFF0000) == 0xFF0000) { nodesX |= level; return true; }
		return false;
	}
	

	
	// method locates the best fitting octant container of the supplied coordinate, returning EITHER a branch OR a leaf octant
	// note: method will search THIS octant object as a tree
	// method will stop if level is reached, setting level to maxLevel is a way to "ignore" the level limit
	public FEM1Octant locateCoordinate(double x, double y, double z, short level) {
		FEM1Octant oct = this, oct2 = this;
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
	public FEM1Octant locateCoordinate(double x, double y, double z) {
		FEM1Octant oct = this, oct2 = this;
		while (oct2.octant != null) {
			if (x <= oct2.xC) {	if (y <= oct2.yC) {	if (z <= oct2.zC)	{ oct2 = oct2.octant[OCT_MMM]; } else { oct2 = oct2.octant[OCT_MMP]; }}
								else {				if (z <= oct2.zC)	{ oct2 = oct2.octant[OCT_MPM]; } else { oct2 = oct2.octant[OCT_MPP]; }}
			} else {			if (y <= oct2.yC) {	if (z <= oct2.zC)	{ oct2 = oct2.octant[OCT_PMM]; } else { oct2 = oct2.octant[OCT_PMP]; }}
								else {				if (z <= oct2.zC)	{ oct2 = oct2.octant[OCT_PPM]; } else { oct2 = oct2.octant[OCT_PPP]; }}}
			if (oct2 == null) break;
			oct = oct2;
		}
		return oct;
	}
	
	
	
	// finds symmetric distant neighbour of this octant along the three x/y/z planes defined by the three bits in diff, the specific plane masked by dirs
	// we must trace down and up again, saving the interim x,y or x,z or y,z statuses while symmetrically walking around x/y/z plane
	// diff: x/bit 0 -> find Pxx neighbour if =1, find Mxx neighbour if =0, and same for y/bit 1 and z/bit 2
	// basically, as long as bits differ (meaning the type of octant, masked by dirs, differs), the XOR will flag for continuing towards root
	// if create=true, the missing neighbor and it's ancestry will be instantiated, they will be marked IN_PROGRESS
	private FEM1Octant symmetricDistantNeighbour(int diff, int dirs, byte[] retrace, boolean create) {
		FEM1Octant octant1 = this, octant1b;
		int r = 0, o;
		while (((octant1.status ^ diff) & dirs) != 0) {
			retrace[r++] = (byte)(octant1.status & 7); octant1 = octant1.parent;
			if (octant1.parent == null) return null;					// this case attainable only if seeking neighbour outside root octant
		}
		octant1b = octant1.parent;
		octant1 = octant1b.octant[o = (octant1.status&7)^dirs];			// mirror across the found branch divider octant (the lowest will always be root)
		
		if (create) {													// construction of symmetric neighbour & ancestry requested?			
			if (octant1 == null) {										// the case of mirroring octant not having the mirrored branch at all
				octant1 = octant1b.octant[o] = new FEM1Octant(octant1b, (short)o, true, true);
			} while (r > 0) {
				if (octant1.octant == null) octant1.octant = new FEM1Octant[8];
				if (octant1.octant[o = retrace[--r] ^ dirs] == null) {
					//if (forIST && octant1.internal()) return octant1;	// for Isosurface Stuffing, do not create internal leaves
					octant1.octant[o] = new FEM1Octant(octant1, (short)o, true, true);
				}
				// TODO: can add a check if the created octant is a lattice leaf created internally, then aborting the creation
				octant1 = octant1.octant[o];
			}
			return octant1;
		} else {
			if (octant1 == null) return null;							// the case of the plane-mirrored octant not existing
			while (r > 0 && octant1.octant != null) {
				octant1b = octant1;
				octant1 = octant1.octant[retrace[--r] ^ dirs];			// toggle to enforce plane-mirrored octant during retracing
				if (octant1 == null) return octant1b;
			}
			return octant1;
		}
	}

	
	
	// method used by Isosurface Stuffing Tetrahedronisation algorithm
	// method utilises the boundary encounter datagrid of FEM1 colution to tell if octant's centerpoint is in/out of the mesh boundary
	public boolean octantCentroidInternal(FEM1Octree octree) {
		// some basic sanity checks that are only needed in a more flexible environment
		//if (!octree.latticeTree || octree.fem.bLatticeNode == null || octree.fem.bLatticeStatus == null) return false;
		int i = (int)((yC - octree.root.yM) / octree.latticeSubdivs), j = (int)((xC - octree.root.xM) / octree.latticeSubdivs);

		if (octree.enforcedLevel == level) i += octree.latticeSubdivs + 1;		// if this is a leaf, use the leaf centroid calculations datafield instead
		double[] fEncounter = octree.fem.fEncounterGrid[i][j];
		int hits = octree.fem.fStatusGrid[i][j] * 3, h = 2;
		int internal = 0;
		while (h < hits && fEncounter[h] < zC) { h += 3; internal = ~internal; }	// flip internal with every encounter less than centroid
		return internal != 0;				// thus, if we hit 1 wall, we're inside, if we hit 2 (enter & exit), we're outside, and so on
	}
	
	
	// method locates the best fitting octant container of the supplied coordinate, then adds the coordinate as the supplied node index
	// if the found topmost container was a branch, method creates a suboctant, inserting node as the only member
	// method returns the octant that became final container of the inserted node
	// the nodeX flag indicates whether to insert a real or extraneous node
	public FEM1Octant addNode(FEM1Octree octree, double x, double y, double z, int n, boolean nodeX) {
		int topmost = -1;
		FEM1Octant oct = this, oct2 = this;
		while (oct2.octant != null && oct2.level < level) {
			if (x <= oct2.xC) {
				if (y <= oct2.yC) {
						if (z<=oct2.zC)	{ oct2.nodes++; oct2 = oct2.octant[OCT_MMM]; } else { oct2.nodes++; oct2 = oct2.octant[OCT_MMP]; }}
				else {	if (z<=oct2.zC)	{ oct2.nodes++; oct2 = oct2.octant[OCT_MPM]; } else { oct2.nodes++; oct2 = oct2.octant[OCT_MPP]; }}
			} else {
				if (y <= oct2.yC) {
						if (z<=oct2.zC)	{ oct2.nodes++; oct2 = oct2.octant[OCT_PMM]; } else { oct2.nodes++; oct2 = oct2.octant[OCT_PMP]; }}
				else {	if (z<=oct2.zC)	{ oct2.nodes++; oct2 = oct2.octant[OCT_PPM]; } else { oct2.nodes++; oct2 = oct2.octant[OCT_PPP]; }}}
			if (oct2 == null) break;
			oct = oct2;
		}
		if (oct2 == null) {															// if containing child octant within branch octant isn't initialised
			topmost = status & 7;
			oct2 = oct.octant[topmost] = new FEM1Octant(oct, (short)topmost, 0, new int[extra_space(1)], 0, null, 0, null, true, true);
			if (oct2.level > octree.topLevel) octree.topLevel = oct2.level;
			if (nodeX) oct2.nodesX++; else oct2.nodes++;	
			oct2.nodeI[0] = n;
		} else {																	// if target leaf container octant was found (or this is octree was a leaf)
			oct2.octantAddNode(n, nodeX);
			if (oct2.nodes > octree.maxNodes) oct2.split(octree, true, (short)0);	// if node count exceeds maxItems, split up the octant
		}
		oct2.finish();
		return oct2;
	}
	
	// method directly inserts node into a known octant, whether as extraneous or as a real node
	public void octantAddNode(int n, boolean nodeX) {
		int nodesTot = nodes + nodesX;
		if (nodeI == null) {
			nodeI = new int[8]; nodeI[nodes++] = n;		
		} else if (nodesTot + 1 >= nodeI.length) {													// need to expand octant's node array?
			int[] nodeNew = new int[FEM1Octant.extra_space(nodesTot + 1)];
			
			if (nodeX) {for (int n1 = 0; n1 < nodesTot; n1++) nodeNew[n1] = nodeI[n1];			// if adding an extraneous node, copy everything straight-on
						nodeNew[nodesTot] = n; 													// add extraneous node at end
						nodesX++; }
			else {		for (int n1 = 0; n1 < nodes; n1++) nodeNew[n1] = nodeI[n1];				// if adding real node, leave gap between nodes & nodesX
						for (int n1 = nodes, n2 = n1 + 1; n1 < nodesTot; n1++, n2++) nodeNew[n2] = nodeI[n1];
						nodeNew[nodes++] = n; }													// add real node in gap
			nodeI = nodeNew;
		} else {
			if (nodeX)	{	nodeI[nodesTot] = n; nodesX++; } 									// regular append at end of array, in nodesX interval
			else {			nodeI[nodesX] = nodeI[nodes]; nodeI[nodes++] = n; }					// move out one node from nodes interval, insert new node
		}
	}

	
	
	// method divides down an edge into the 8 octants, bitflagging the intersected octants, using recursive calls
	// method records splitted edge segments and shifts them toward their proper octant offsets within split[] array
	// the segments are then accessible as 2 coord triplets by steps of (octantNo*6) and their presense is indicated by the bitflags
	// octant bitflags (high to low): PPP,MPP,PMP,MMP,PPM,MPM,PMM,MMM
	int edgeMembership(double xa, double ya, double za, double xb, double yb, double zb, double[] split, int offs, int dim, int bits) {	
		switch (dim) {
		case 0:
			if (xa <= xC) {
				if (xC < xb) {																		// point a in lower x-octants, b in higher
						double xbaD = 1. / (xb - xa), f1 = (xC - xa) * xbaD, f2 = (xb - xC) * xbaD;
						double ys = f1 > f2 ? ya+(yb-ya)*f1 : yb-(yb-ya)*f2, zs = f1 > f2 ? za+(zb-za)*f1 : zb-(zb-za)*f2;
						split[0]=xa; split[1]=ya; split[2]=za; split[3]=xC; split[4]=ys; split[5]=zs;
						split[6]=xC; split[7]=ys; split[8]=zs; split[9]=xb; split[10]=yb; split[11]=zb;
						return 	edgeMembership(xa, ya, za, xC, ys, zs, split, 0, dim+1, 0x55) | edgeMembership(xC, ys, zs, xb, yb, zb, split, 6, dim+1, 0xAA);
				} else if (xa < xM && xb < xM) return 0;											// edge outside octant
				else {	split[0]=xa; split[1]=ya; split[2]=za; split[3]=xb; split[4]=yb; split[5]=zb;
						return 	edgeMembership(xa, ya, za, xb, yb, zb, split, 0, dim+1, 0x55);		// edge completely in OCT_M** domain
				}
			} else if (xb <= xC) {
				if (xC < xa) {																		// reverse: point b in lower x-octants, a in higher
						double xbaD = 1. / (xa - xb), f1 = (xC - xb) * xbaD, f2 = (xa - xC) * xbaD;
						double ys = f1 > f2 ? yb+(ya-yb)*f1 : ya-(ya-yb)*f2, zs = f1 > f2 ? zb+(za-zb)*f1 : za-(za-zb)*f2;
						split[0]=xb; split[1]=yb; split[2]=zb; split[3]=xC; split[4]=ys; split[5]=zs;
						split[6]=xC; split[7]=ys; split[8]=zs; split[9]=xa; split[10]=ya; split[11]=za;
						return 	edgeMembership(xb,yb,zb,xC,ys,zs, split, 0, dim+1, 0x55) | edgeMembership(xC,ys,zs,xa,ya,za, split, 6, dim+1, 0xAA);
				} else if (xa < xM && xb < xM) return 0;											// edge outside octant
				else {	split[0]=xb; split[1]=yb; split[2]=zb; split[3]=xa; split[4]=ya; split[5]=za;
						return 	edgeMembership(xb, yb, zb, xa, ya, za, split, 0, dim+1, 0x55);		// edge completely in OCT_M** domain
				}
			} else if (xa > xP && xb > xP) return 0;												// edge outside octant		
			else {		split[6]=xa; split[7]=ya; split[8]=za; split[9]=xb; split[10]=yb; split[11]=zb;
						return	edgeMembership(xa, ya, za, xb, yb, zb, split, 6, dim+1, 0xAA);		// edge completely in OCT_P** domain	
			}
		case 1:
			if (ya <= yC) {
				if (yC < yb) {																		// point a in lower y-octants, b in higher
						double ybaD = 1. / (yb - ya), f1 = (yC - ya) * ybaD, f2 = (yb - yC) * ybaD;
						double xs = f1 > f2 ? xa+(xb-xa)*f1 : xb-(xb-xa)*f2, zs = f1 > f2 ? za+(zb-za)*f1 : zb-(zb-za)*f2;
						int offY=offs+12; split[offs++]=xa; split[offs++]=ya; split[offs++]=za; split[offs++]=xs; split[offs++]=yC; split[offs++]=zs;
						split[offY++]=xs; split[offY++]=yC; split[offY++]=zs; split[offY++]=xb; split[offY++]=yb; split[offY++]=zb;
						return 	bits & (edgeMembership(xa,ya,za,xs,yC,zs, split, offs-6, dim+1, 0x33) | edgeMembership(xs,yC,zs,xb,yb,zb, split, offY-6, dim+1, 0xCC));
				} else if (ya < yM && yb < yM) return 0;											// edge outside octant
				else	return 	bits & edgeMembership(xa,ya,za,xb,yb,zb, split, offs, dim+1, 0x33);	// edge completely in OCT_*M* domain
			} else if (yb <= yC) {
				if (yC < ya) {																		// reverse: point b in lower y-octants, a in higher
						double ybaD = 1. / (ya - yb), f1 = (yC - yb) * ybaD, f2 = (ya - yC) * ybaD;
						double xs = f1 > f2 ? xb+(xa-xb)*f1 : xa-(xa-xb)*f2, zs = f1 > f2 ? zb+(za-zb)*f1 : za-(za-zb)*f2;
						int offY=offs+12; split[offs++]=xb; split[offs++]=yb; split[offs++]=zb; split[offs++]=xs; split[offs++]=yC; split[offs++]=zs;
						split[offY++]=xs; split[offY++]=yC; split[offY++]=zs; split[offY++]=xa; split[offY++]=ya; split[offY++]=za;
						return 	bits & (edgeMembership(xb,yb,zb,xs,yC,zs, split, offs-6, dim+1, 0x33) | edgeMembership(xs,yC,zs,xa,ya,za, split, offY-6, dim+1, 0xCC));
				} else if (ya < yM && yb < yM) return 0;											// edge outside octant
				else	return 	bits & edgeMembership(xb,yb,zb,xa,ya,za, split, offs, dim+1, 0x33);	// edge completely in OCT_*M* domain
			} else if (ya > yP && yb > yP) return 0;												// edge outside octant		
			else {		int offY=offs+12; split[offY++]=split[offs++]; split[offY++]=split[offs++]; split[offY++]=split[offs++];
						split[offY++]=split[offs++]; split[offY++]=split[offs++]; split[offY++]=split[offs++];
						return	bits & edgeMembership(xa,ya,za,xb,yb,zb, split, offY-6, dim+1,0xCC);// edge completely in OCT_*P* domain	
			}
		case 2:
			if (za <= zC) {
				if (zC < zb) {
						double zbaD = 1. / (zb - za), f1 = (zC - za) * zbaD, f2 = (zb - zC) * zbaD;
						double xs = f1 > f2 ? xa+(xb-xa)*f1 : xb-(xb-xa)*f2, ys = f1 > f2 ? ya+(yb-ya)*f1 : yb-(yb-ya)*f2;
						int offZ=offs+24; split[offs++]=xa; split[offs++]=ya; split[offs++]=za; split[offs++]=xs; split[offs++]=ys; split[offs++]=zC;
						split[offZ++]=xs; split[offZ++]=ys; split[offZ++]=zC; split[offZ++]=xb; split[offZ++]=yb; split[offZ++]=zb;
						return 	bits;																// point a in lower z-octants, b in higher	
				} else if (za < zM && zb < zM) return 0;											// edge outside octant
				else	return 	bits & 0x0F;														// edge completely in OCT_**M domain
			} else if (zb <= zC) {
				if (zC < za) {
						double zbaD = 1. / (za - zb), f1 = (zC - zb) * zbaD, f2 = (za - zC) * zbaD;
						double xs = f1 > f2 ? xb+(xa-xb)*f1 : xa-(xa-xb)*f2, ys = f1 > f2 ? yb+(ya-yb)*f1 : ya-(ya-yb)*f2;
						int offZ=offs+24; split[offs++]=xb; split[offs++]=yb; split[offs++]=zb; split[offs++]=xs; split[offs++]=ys; split[offs++]=zC;
						split[offZ++]=xs; split[offZ++]=ys; split[offZ++]=zC; split[offZ++]=xa; split[offZ++]=ya; split[offZ++]=za;
						return 	bits;																// reverse: point b in lower x-octants, a in higher	
				} else if (za < zM && zb < zM) return 0;											// edge outside octant
				else	return 	bits & 0x0F;														// edge completely in OCT_**M domain
			} else if (za > zP && zb > zP) return 0;												// edge outside octant		
			else {		int offZ=offs+24; split[offZ++]=split[offs++]; split[offZ++]=split[offs++]; split[offZ++]=split[offs++];
						split[offZ++]=split[offs++]; split[offZ++]=split[offs++]; split[offZ++]=split[offs++];
						return	bits & 0xF0;														// edge completely in OCT_**P domain	
			}
		}
		return 0;																					// never ends up here
	}
	
	// method inserts an edge segment into octree by recursive splitting & membership testing
	public void addEdge(FEM1Octree octree, double xa, double ya, double za, double xb, double yb, double zb, int e, int maxEdges, double disbalanceGrowth) {
		
		if (octant == null) {															// if a leaf was reached
			octantAddEdge(e);
			// TODO: splitting demands systemic overhaul for dynamic edge operations, a step only necessary for a properly interactive FEA system
			// TODO: for true dynamic action, split() would need an new inserted edge and either existing edge segments or resplit edges from FEM1 model
			// TODO: currently, the supplied edge index is simply inserted into the correct leaf containers
//			if (level < maxLevel &&														// split this octant if the 2 criterions are fulfilled
//					(nodes > fem.octreeMaxItemsN || edges > fem.octreeMaxItemsE || facets > fem.octreeMaxItemsF)) {
//				int newLeaves = split(fem, true) - 1;
//				if (octree.topLevel < level + 1) octree.topLevel = (short)(level + 1);
//				FEM1Octree octant1 = this;
//				while (octant1 != null) {												// backpropagate new leaf & branch counts
//					octant1.leaves += newLeaves; octant1.branches++; octant1 = octant1.parent;
//				}
//			}
			return;
		}
		double[] split = new double[6 * 8];
		int bits = edgeMembership(xa, ya, za, xb, yb, zb, split, 0, 0, 0xFF);			// what suboctants does the edge belong to, set the proper bits
		
		// DEBUG: wasteful to do rebalancing along with every insertion, better to do infrequent calling of a balanceOctree() method
//		int[] eDisb = { 0,0,0,0,0,0,0,0 };
//		for (int o = 0, bits2 = bits; o < 8; o++, bits2 >>= 1)							// precalculate expected branch disbalance with the added edges
//			if ((bits2 & 1) != 0 && octant[o] != null)
//				eDisb[o] = octant[o].nodes + octant[o].edges + 1;
//		disbalance = octantNodeDisbalance(edges, eDisb[0], eDisb[1], eDisb[2], eDisb[3], eDisb[4], eDisb[5], eDisb[6], eDisb[7]);

		int o = bits < 0x0F ? 0 : ((bits & 0x0F)==0 ? 4 : 0); if (o == 4) bits >>= 4;
		while (bits != 0) {																// start adding edges
			if ((bits & 1) == 1) {
				if (octant[o] == null) {												// if containing octant doesn't exist
					if (disbalanceGrowth < octree.disbalanceCap) {						// if disbalanceCap allows dividing edge containment further
						FEM1Octant octant1 = octant[o] = new FEM1Octant(this, (short)o, true, true);
						if (octant1.level > octree.topLevel) octree.topLevel = octant1.level;
						octant1.octant = new FEM1Octant[8];								// create transitive branch
						octant1.disbalance = 1;											// enforce transitive branch creation until disbalance limit
						int oC = o * 6;
						octant1.addEdge(octree,split[oC++],split[oC++],split[oC++],split[oC++],split[oC++],split[oC++],e,maxEdges, (1+disbalanceGrowth)*.5);
						octant1.finish();
					} else {															// case of disbalanceCap reached, doing one leaf only and returning
						octant[o] = new FEM1Octant(this, (short)o, true, true);			// make leaf octant with one edge
						if (octant[o].level > octree.topLevel) octree.topLevel = octant[o].level;
						octant[o].octantAddEdge(e);
						FEM1Octant octant1 = this;
						synchronized (this) {											// note: concurrent edge adding WILL mess up statistics if unsynchronised
							while (octant1 != null) {									// backpropagate new leaf & branch counts
								octant1.leaves++; octant1.branches++; octant1 = octant1.parent; }
						}
						octant[o].finish();
					}
				} else {																// we're not at leaf level yet, continue
					int oC = o * 6;
					octant[o].addEdge(octree, split[oC++], split[oC++], split[oC++], split[oC++], split[oC++], split[oC++], e, maxEdges, disbalance);
				}
			}
			bits >>= 1; o++;
		}
		edges++;		// note: even though one edge can land in several leaves, it is still counted as a single edge at all branching levels	
	}
	
	
	// method directly inserts edge into a known octant
	public void octantAddEdge(int e) {
		if (edgeI == null) edgeI = new int[8];											// need to initialise coctant's edges?
		else if (edges + 1 >= edgeI.length) {											// need to expand octant's edge array?
			int[] edgeNew = new int[FEM1Octant.extra_space(edges + 1)];		
			for (int n1 = 0; n1 < edges; n1++) edgeNew[n1] = edgeI[n1];					// copy everything straight-on
			edgeNew[edges++] = e; edgeI = edgeNew;										// add extraneous edge at end
			return; }
		edgeI[edges++] = e;
	}

	
	
	// method returns array of all leaves within supplied subtree, or supplied subtree as single member if it's a leaf itself
	// method utilises nonrecursive tree walk with expansion into a stack
	// if forIST=true, picks only octants within IST gradation range an enumerate them by "nodes" parameter
	public FEM1Octant[] octantArray(FEM1Octree octree, int levelSought, int[] count, boolean forIST) {
		
		if (octant == null) {															// octant was a leaf itself
			FEM1Octant[] leafArray = { this };
			if (forIST) this.nodes = 0;													// need to enumerate?
			return leafArray;
		}
		int gradLwr = octree.maxLevel - octree.gradations + 1, indexIST = 0;
		if (gradLwr <= 0) gradLwr = 1;
		FEM1Octant[] leafArray; // = new FEM1Octant[this.leaves];						// worst-case allocation
		if (forIST) {																	// if generating for Isosurface Stuffing, sum gradation levels counts
				int l = gradLwr, totLAsize = 0;											// this gives a worst-case estimate
				while(l <= octree.maxLevel) totLAsize += octree.levelSums[l++];
				leafArray = new FEM1Octant[totLAsize + octree.levelSums[octree.maxLevel]];
		} else	leafArray = new FEM1Octant[levelSought < 0 ? this.leaves : octree.levelSums[levelSought]];
		
		FEM1Octant[] stackO = new FEM1Octant[octree.topLevel - level + 1];				// allocate for octant stack
		stackO[0] = this;
		int[] stack = new int[octree.topLevel - level + 1];								// stack allocation
		stack[0] = iterator;
		int s = 0, leafCnt = 0;
		while (s >= 0) {
			int c = stack[s] & 0xFF000000;												// get remaining children of current octant on stack
			if (c > 0) {						
				FEM1Octant octantChild = stackO[s].octant[stack[s] & 7];				// get outstanding child
				stack[s] = ((stack[s]&0xFFFFFF)>>3) | (c-0x1000000);					// recompose octant stack entry: pop child, make next outstanding
				stack[++s] = octantChild.iterator;
				stackO[s] = octantChild;
				// on IST request, gather up any octant within the gradations range that is either a leaf or whose centroid is internal
				if (forIST && octantChild.level >= gradLwr && (octantChild.level==octree.maxLevel || octantChild.centroid_internal())) {
					leafArray[leafCnt++] = octantChild;
					octantChild.nodes = indexIST++; octantChild.enumerator = 0;
				}
			} else {
				// if it was a leaf of correct level on top of octant stack (ignore this block if IST collection was requested)
				if (!forIST && stackO[s].octant == null && (levelSought < 0 || stackO[s].level == levelSought)) {
					leafArray[leafCnt++] = stackO[s--];									// confirm latest leaf octant, backtrack to parent
				} else s--;
			}
		}
		if (count != null) count[0] = leafCnt;
		return leafArray;
	}

	// recursive version of leafOctantArray() method, microbenchmarking indicates this method is about 30 percent slower, despite lower code complexity
	public FEM1Octant[] octantArrayR(FEM1Octree octree, int levelSought) {
		if (octant == null) {
			if (levelSought < 0 || level == levelSought) { FEM1Octant[] leafArray = {this}; return leafArray; }
			else return new FEM1Octant[1];
		}
		FEM1Octant[] leafArray = new FEM1Octant[this.leaves];							// worst-case allocation
		octantArrayRcall(leafArray, new int[1], levelSought);
		return leafArray;
	}
	public void octantArrayRcall(FEM1Octant[] leafArray, int[] leafCnt, int levelSought) {
		for (int o = 0; o < 8; o++)
			if (octant[o] != null) {
				if (octant[o].octant == null) {
					if (levelSought < 0 || octant[o].level == levelSought) leafArray[leafCnt[0]++] = octant[o];
				} else octant[o].octantArrayRcall(leafArray, leafCnt, levelSought);
			}
	}

	
	// method returns array of arrays of octants of same level, sorted from low to high, from any branch of an octree
	// method utilises nonrecursive tree walk with expansion into a stack
	public FEM1Octant[][] layerOctantArray(FEM1Octree octree) {
		
		if (octant == null) { FEM1Octant[][] layerArray = { {this} }; return layerArray; }	// octant was a leaf itself
		
		int levels = octree.topLevel - level + 1;
		octree.root.process(octree, OCTREE_SUM_LEVELS);
		FEM1Octant[][] layerArray = new FEM1Octant[levels][];							// create the array of octant arrays
		for (int l = level; l <= octree.topLevel; l++) layerArray[l] = new FEM1Octant[octree.levelSums[l]];
		int[] layerCount = new int[levels];												// sum up octant counts per level here
			
		FEM1Octant[] stackO = new FEM1Octant[levels];									// allocate for octant stack
		stackO[0] = this;
		int[] stack = new int[levels];													// stack allocation
		stack[0] = iterator;
		int s = 0;
		while (s >= 0) {
			int c = stack[s] & 0xFF000000;												// get remaining children of current octant on stack
			if (c > 0) {						
				FEM1Octant octantChild = stackO[s].octant[stack[s] & 7];				// get outstanding child
				stack[s] = ((stack[s]&0xFFFFFF)>>3) | (c-0x1000000);					// recompose octant stack entry: pop child, make next outstanding
				stack[++s] = octantChild.iterator;
				stackO[s] = octantChild;
			} else {																	// this octant has no more children (or was leaf)
				FEM1Octant octant1 = stackO[s];
				layerArray[octant1.level][layerCount[octant1.level]++] = octant1;		// store at correct array layer, backtrack to parent
				s--;
			}
		}
		return layerArray;
	}



	// method returns array of all leaves close enough to a coordinate within supplied octant, or supplied octant as single member if it's a leaf itself
	// any octant must be within "dist" distance, method utilises nonrecursive tree walk with DFS expansion into a stack
	public FEM1Octant[] octantArrayByDistance(FEM1Octree octree, double x, double y, double z, double dist) {
		
		if (octant == null) { FEM1Octant[] leafArray = { this }; return leafArray; }	// octant was a leaf itself
		
		FEM1Octant[] leafArray = new FEM1Octant[leaves];								// worst-case-allocate for octant leaves
		FEM1Octant[] stackO = new FEM1Octant[octree.topLevel - level + 1];				// worst-case-allocate for octant stack
		stackO[0] = this;
		int[] stack = new int[octree.topLevel - level + 1];								// worst-case stack allocation
		stack[0] = iterator;
		int s = 0, leafCnt = 0;
		while (s >= 0) {
			int c = stack[s] & 0xFF000000;												// get remaining children of current octant on stack
			if (c > 0) {						
				FEM1Octant octantChild = stackO[s].octant[stack[s] & 7];				// get outstanding child
				stack[s] = ((stack[s]&0xFFFFFF)>>3) | (c-0x1000000);					// recompose octant stack entry: pop child, make next outstanding
				if (octantChild.octantDistance(x, y, z, false) < dist) {				// stack child octant only if it fulfills closeness criterion
					stack[++s] = octantChild.iterator;
					stackO[s] = octantChild;
				}
			} else {																	// if no children remain/exist
				if (stackO[s].octant == null)											// if it was a leaf on top of octant stack
					leafArray[leafCnt++] = stackO[s]; 									// add to leaf collection										
				s--;																	// backtrack
			}
		}
		return leafArray;
	}

	
	
	// method finds out if a bounding box overlaps octant, returning bitfield of the suboctants that are overlapped (=0xFF if overlapping entire octant)
	// method will indicate a corner overlap condition if the three bits are set in upper half of bitfield
	public int octantBBoxOverlap(double xMbb, double yMbb, double zMbb, double xPbb, double yPbb, double zPbb) {
		int status = 0;
		if (xPbb <= xC) {		if (xPbb >= xM) status = 0x10055; else return 0; }
		else if (xMbb >= xC) {	if (xMbb <= xP) status = 0x100AA; else return 0; }
		else status = 0xFF;
		if (yPbb <= yC) {		if (yPbb >= yM) status = (status & 0x33) | 0x20000; else return 0; }
		else if (yMbb >= yC) {	if (yMbb <= yP) status = (status & 0xCC) | 0x20000; else return 0; }
		if (zPbb <= zC) {		if (zPbb >= zM) status = (status & 0x0F) | 0x30000; else return 0; }
		else if (zMbb >= zC) {	if (zMbb <= zP) status = (status & 0xF0) | 0x30000; else return 0; }
		return status;
	}
	
	
	// method returns array of all leaves/octants overlapped by supplied bounding box, or supplied octree as single member if it's a leaf
	// method utilises nonrecursive tree walk with DFS expansion into a stack and dynamic adjustment of overlapArray[] allocator
	// method can be used for fast insertion of SMALLER objects: faces, tetrahedra, without caring about their INTRINSIC geometry within octree
	// TODO: in mixed data octrees, good idea to supply a parameter on what type of data to look for, early-on capping branches that lack this data
	public FEM1Octant[] octantArrayByBBox(FEM1Octree octree,
							double xMbb,double yMbb,double zMbb,double xPbb,double yPbb,double zPbb, boolean leavesOnly, boolean skipCorners) {
		
		if (octant == null) { FEM1Octant[] overlapArray = { this }; return overlapArray; }	// octant was a leaf
		
		FEM1Octant[] overlapArray = new FEM1Octant[octree.oOBB_allocSize];	// assume max 8 overlap results from beginning
		FEM1Octant[] stackO = new FEM1Octant[octree.topLevel - level + 1];	// allocate for octant stack
		stackO[0] = this;
		int[] stack = new int[octree.topLevel - level + 1];					// stack allocation
		stack[0] = octantBBoxOverlap(xMbb, yMbb, zMbb, xPbb, yPbb, zPbb);
		int s = 0, overlaps = 0;
		while (s >= 0) {
			int v = stack[s] & 0xFF;										// get remaining children of current octant on stack
			if (v != 0) {													// if any overlapped children remain	
				int c = stack[s] & 0xFF00;									// extract the counter kept in 0xFF00
				while ((v&1) == 0) { v >>= 1; c += 0x0100; }				// skip nonexistent octants
				FEM1Octant octantChild = stackO[s].octant[c >> 8];			// get overlapped child (counter c indexes it)
				int cornerStatus = stack[s] & 0xFFFF0000;
				stack[s] = cornerStatus | (c+0x0100) | (v>>1);				// recompose octant stack entry: iterator + remaining bitflags
				if (octantChild == null) continue;							// note: octantBBoxOverlap() does not check existence of a child
				if ((!leavesOnly || octantChild.octant==null) && 			// if we accept any octants or we found a a leaf
						(!skipCorners || cornerStatus != 0x30000))	{		// if we're skipping corners and corner indicator wasn't set
					
					if (overlaps >= overlapArray.length) {					// does overlapArray[] need to be extended?
						FEM1Octant[] leafArrayNew = new FEM1Octant[overlapArray.length * 2];	// double the size
						for (int i = 0; i < overlaps; i++) leafArrayNew[i] = overlapArray[i];
						overlapArray = leafArrayNew;
						// do heuristic on running ratio of allocation-free to reallocated calls to find out if a new alloc.size needs to be adopted
						octree.oOBB_lAllocs_successRate = (float)((octree.oOBB_lAllocs_successRate +
																	octree.oOBB_lAllocs - ++octree.oOBB_lAllocs_fail) / 2.);
						if (octree.oOBB_lAllocs_successRate < 100)	octree.oOBB_allocSize <<= 1;
						else if (octree.oOBB_lAllocs_successRate > 1000 && octree.oOBB_allocSize > 4)	octree.oOBB_allocSize >>= 1; 
						if (DEBUG_LEVEL > 2) System.out.println("FEM1Octant.octantsOverlappingBBox() readjusted alloc.size: "
																	+ octree.oOBB_allocSize + " " + octree.oOBB_lAllocs_successRate);
						octree.oOBB_lAllocs_fail = octree.oOBB_lAllocs = 0;
					} else ++octree.oOBB_lAllocs;							// this was a "success" (no reallocation)

					overlapArray[overlaps++] = octantChild; 				// add to overlaps collection
				} else {													// if it was a branch, find out it's children are overlapped
					stack[++s] = octantChild.octantBBoxOverlap(xMbb, yMbb, zMbb, xPbb, yPbb, zPbb);
					if (stack[s] != 0) stackO[s] = octantChild; else s--;	// store only if at least one child branch is overlapped
				}
			} else s--;														// if no children remain/exist, backtrack
		}
		return overlapArray;
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
	final static double OCT_MARGIN = 1e-10;
	public int closestNode(FEM1Octree octree, int n, double[] distance, FEM1Octant[] octants) {
		
		double[] distC = {0}, nodeG = octree.fem.node;
		FEM1Octant octant1;
		if (octants == null || octants[0] == null) {
			int n3 = n * 3;
			if (n3 > nodeG.length) {
					n3 -= nodeG.length;
					octant1 = locateCoordinate(octree.fem.nodeWork[n3++], octree.fem.nodeWork[n3++], octree.fem.nodeWork[n3], octree.maxLevel);
			} else	octant1 = locateCoordinate(nodeG[n3++], nodeG[n3++], nodeG[n3], octree.maxLevel);
		} else
			octant1 = octants[0];
		
		int n3 = n * 3, nClosest;
		double x, y, z;
		if (n3 > nodeG.length) {	n3 -= nodeG.length; x = octree.fem.nodeWork[n3++]; y = octree.fem.nodeWork[n3++]; z = octree.fem.nodeWork[n3]; }
		else				 	  {	x = nodeG[n3++]; y = nodeG[n3++]; z = nodeG[n3]; }		
		FEM1Octant octantC = null;
		
		// this is the case of only one node in container, start by finding the first ancestral relation branch besides direct parentage
		if (octant1.nodes == 1) {
			FEM1Octant octantRB = octant1;
			while ((octantRB.iterator >> 24) <= 1 || octantRB.nodes <= 1)		// backtrack until some other relation branch (with nodes) found
				if (octantRB == this) {											// special case: entire octree has only one node, so return itself
					if (distance != null) distance[0] = 0;
					if (octants != null) octants[1] = null;
					return n;
				}
				else octantRB = octantRB.parent;
			FEM1Octant[] octantsN = octantRB.octantArray(octree, -1, null, false);	// gather all leaf octants of ancestor with relation branches

			double dist = Double.MAX_VALUE;
			for (FEM1Octant octantN : octantsN) {								// find the closest leaf octant in array
				if (octantN == null) break;
				if (octantN == octant1 || octantN.nodes == 0) continue;			// skip node's octant, skip nodeless octants
				double distN = octantN.octantDistance(x, y, z, false);
				if (distN < dist) { dist = distN; octantC = octantN; }
			}
			nClosest = octree.fem.closestNode(x, y, z, octantC.nodeI, distC);	// find the closest node in closest leaf octant
		} else {
			// this is the base case of there being more than one node within container
			nClosest = octree.fem.closestNode(x, y, z, octant1.nodeI, distC);	// find out what local node that is closest
			if (octant1 == this) {												// if this is "archparent" octant (has no neighbours), we're done 
				if (distance != null) distance[0] = Math.sqrt(distC[0]);
				if (octants != null) octants[1] = null;
				return nClosest;
			}
		}
			
		// in both cases we have found a closest node candidate and the distance that will define a testing bounding box
		double dist = distC[0];
		// expand to the smallest ancestor octant that encompasses the bounding box
		double distRt = Math.sqrt(dist);
		double xMn = x - distRt, xPn = x + distRt, yMn = y - distRt, yPn = y + distRt, zMn = z - distRt, zPn = z + distRt;
		FEM1Octant octantE = octant1;
		// expand scope in the direction that isn't overshooting global bounding box and that still is farther than matching octree wall
		if (xMn > octree.fem.bBox[0]) while (octantE.xM >= xMn && octantE != this) octantE = octantE.parent;
		if (xPn < octree.fem.bBox[3]) while (octantE.xP <= xPn && octantE != this) octantE = octantE.parent;
		if (yMn > octree.fem.bBox[1]) while (octantE.yM >= yMn && octantE != this) octantE = octantE.parent;
		if (yPn < octree.fem.bBox[4]) while (octantE.yP <= yPn && octantE != this) octantE = octantE.parent;
		if (zMn > octree.fem.bBox[2]) while (octantE.zM >= zMn && octantE != this) octantE = octantE.parent;
		if (zPn < octree.fem.bBox[5]) while (octantE.zP <= zPn && octantE != this) octantE = octantE.parent;
		if (octantE == octant1)	{													// if smallest encompasser is node's octant itself, then we're done
			if (distance != null) distance[0] = distRt;
			if (octants != null) octants[1] = null;
			return nClosest;
		}
		// call method that collects only the octant's leaves that lie close enough to distance sphere of node's coordinate (criterion neighbours)
		FEM1Octant[] octantsN = octantE.octantArrayByDistance(octree, x, y, z, dist);
		FEM1Octant octantC2 = null;
		
		for (FEM1Octant octantN : octantsN) {										// time to check closeness of nodes of each close-enough leaf octant
			if (octantN == null) break;
			if (octantN == octant1 || octantN == octantC || octantN.nodes == 0)
				continue;															// skip the already tested octant of the node, skip nodeless octants
			int nClosestN = octree.fem.closestNode(x, y, z, octantN.nodeI, distC);	// find out closest node in criterion-neighbour octant
			if (distC[0] < dist) {													// if it was closer than closest node in this octant, assign to it
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
	static float octantNodeDisbalance(int total, int nMMM, int nPMM, int nMPM, int nPPM, int nMMP, int nPMP, int nMPP, int nPPP) {
		double totalD = 0.5773 / total, vDisbX = 0, vDisbY = 0, vDisbZ = 0;
		if (nMMM > 0) { double nMMMf = (double)nMMM * totalD; vDisbX -= nMMMf; vDisbY -= nMMMf; vDisbZ -= nMMMf; }
		if (nPMM > 0) { double nPMMf = (double)nPMM * totalD; vDisbX += nPMMf; vDisbY -= nPMMf; vDisbZ -= nPMMf; }
		if (nMPM > 0) { double nMPMf = (double)nMPM * totalD; vDisbX -= nMPMf; vDisbY += nMPMf; vDisbZ -= nMPMf; }
		if (nPPM > 0) { double nPPMf = (double)nPPM * totalD; vDisbX += nPPMf; vDisbY += nPPMf; vDisbZ -= nPPMf; }
		if (nMMP > 0) { double nMMPf = (double)nMMP * totalD; vDisbX -= nMMPf; vDisbY -= nMMPf; vDisbZ += nMMPf; }
		if (nPMP > 0) { double nPMPf = (double)nPMP * totalD; vDisbX += nPMPf; vDisbY -= nPMPf; vDisbZ += nPMPf; }
		if (nMPP > 0) { double nMPPf = (double)nMPP * totalD; vDisbX -= nMPPf; vDisbY += nMPPf; vDisbZ += nMPPf; }
		if (nPPP > 0) { double nPPPf = (double)nPPP * totalD; vDisbX += nPPPf; vDisbY += nPPPf; vDisbZ += nPPPf; }
		return (float)Math.sqrt(vDisbX * vDisbX + vDisbY * vDisbY + vDisbZ * vDisbZ);
	}
	
	
	// method checks if supplied coordinate lies inside or outside a closed mesh boundary (method assumes nodeWork[] array was integrated)
	// method spans octree either with a user segment (vDir=3) or along an ordinate in +x (vDir=0), +y (vDir=1) or +z (vDir=2) direction
	// method returns null if no mesh penetrations happened and an array of facet encounters in closest-farthest order otherwise
	// method will return total polarity sum and total encounter count in supplied status[] array
	// method can be used for systematic gridded coverage of a mesh topology along a certain ordinate (x/y/z)
	// for domain-spanning ray polarity consistency, the segment should begin outside the mesh, mesh volume must be watertight for expected behaviour
	// status bits (on vDir): 16 = keep polygon check-off table, 32 = find only one intersection
	public double[] facetEncounters(FEM1Octree geocTree, double x1, double y1, double z1, double x2, double y2, double z2, int vDir, int[] status) {
		if (geocTree.fem.nodeWork != null) throw new RuntimeException("FEM1Octree.facetEncounters(): nodeWork[] must be integrated.");
		// DEBUG: a simple ordinate cuboid object will perfectly span the root octant, so checking the coordinates to be inside the global boundbox
		// will make impossible to start checking from outside the mesh, since one then has to step outside the global boundbox
		// might as well rely on the fact that leavesOverlappingBBox() will return null on being outside the domain anyway
		//if (x1 < xM || x1 > xP || y1 < yM || y1 > yP || z1 < zM || z1 > zP) return null;	// base case, coordinate outside octree's boundbox

		boolean clearPolygonCheck = (vDir&16)!=16;											// flag in vDir specifies if visits flagarray is to be kept
		boolean singleEncounter = (vDir&32)==32;											// flag specifies whether only one encounter is necessary
		vDir &= 0x7;
		double x1b = x1, y1b = y1, z1b = z1, x2b = x2, y2b = y2, z2b = z2;
		if (vDir == 3) {
			if (x1 > x2) { x2b = x1; x1b = x2; } else { x1b = x1; x2b = x2; } 				// corrected inputs for boundbox test
			if (y1 > y2) { y2b = y1; y1b = y2; } else { y1b = y1; y2b = y2; } 
			if (z1 > z2) { z2b = z1; z1b = z2; } else { z1b = z1; z2b = z2; } 
			vDir = (x2b-x1b > y2b-y1b ? (x2b-x1b>z2b-z1b?0:2) : (y2b-y1b>z2b-z1b?1:2));		// find out which ordinate to use for sorting order
		} else {
			double vLen = xP - xM + 1;
			switch (vDir) {																	// vDir = 0/1/2 will throw a domain-covering x/y/z ray
				case 0: x2 = x2b = x1 + vLen; y2 = y2b = y1; z2 = z2b = z1; break;			// this case ignores coordinate 2
				case 1: y2 = y2b = y1 + vLen; x2 = x2b = x1; z2 = z2b = z1; break;
				case 2: z2 = z2b = z1 + vLen; x2 = x2b = x1; y2 = y2b = y1; }
		}
		
		FEM1Octant[] octant1 = octantArrayByBBox(geocTree, x1b, y1b, z1b, x2b, y2b, z2b, true, false);	// find leaves only, include corners
		status[0] = status[1] = 0;
		if (octant1[0] == null) return null;
		if (clearPolygonCheck) geocTree.fem.polygonCheck.clearVisits();			// we need a clean checkoff list for facets (unless "do not clear" flagged)
		
		double[] isect = {0,0,0};
		if (singleEncounter) {													// if any encounter requested, find one and return
			for (FEM1Octant o : octant1) {
				if (o == null) break;
				for (int f = 0; f < o.facets; f++) {							// check facet intersections, sum inside/outside status as polarity sum
					int fI = o.facetI[f];
					if (!geocTree.fem.polygonCheck.visited(fI)) {				// do not revisit checked facets
						int polarity = geocTree.fem.facetSegmentIntersection(fI, x1, y1, z1, x2, y2, z2, isect);
						if (polarity != 0) {									// found intersection ?
							status[2] = fI;										// store the intersected facet to give caller another try on it
							status[0] = polarity; status[1]++; return isect;	// return the encounter
						}
						geocTree.fem.polygonCheck.visit(fI);					// check off facet, because some facets are duplicated across octant boundaries
					}
				}
			}
			return null;
		}
		
		double[] encounter = new double[3 * geocTree.maxFacetEncounters];
		int encounters = 0;
		for (FEM1Octant o : octant1) {
			if (o == null) break;
			for (int f = 0; f < o.facets; f++) {								// check facet intersections, sum inside/outside status as polarity sum
				int fI = o.facetI[f];
				if (!geocTree.fem.polygonCheck.visited(fI)) {					// do not revisit checked facets
					int polarity = geocTree.fem.facetSegmentIntersection(fI, x1, y1, z1, x2, y2, z2, isect);
					if (polarity != 0) {										// found intersection ?
						status[0] += polarity;
						//if (encounters < 2) status[2 + encounters] = fI;		// store down facet index for future revisits (limited to 2 facets)
						if (encounter.length < ++encounters * 3) {
							double[] encounterNew = new double[(geocTree.maxFacetEncounters *= 2) * 3];	// double max.encounters for this FEM system
							for (int i = 0; i < encounter.length; i++) encounterNew[i] = encounter[i];
							encounter = encounterNew;
						}
						int i3 = 0;												// add coordinate in order of encounter distance
						for (int i = 1; i < encounters; i++, i3 += 3) {
							if (encounter[i3+vDir] > isect[vDir]) {				// insert in ordinate-decided order, pushing coordinates ahead
								for (int m3t = encounters * 3 - 1, m3f = m3t - 3; m3f > i3;) {
									encounter[m3t--] = encounter[m3f--]; encounter[m3t--] = encounter[m3f--]; encounter[m3t--] = encounter[m3f--];
								}
								break;
							}
						}
						encounter[i3++] = isect[0]; encounter[i3++] = isect[1]; encounter[i3] = isect[2];
					}
					geocTree.fem.polygonCheck.visit(fI);						// check off facet, because some facets are duplicated across octant boundaries
				}
			}
		}
		status[1] = encounters;
		return encounters > 0 ? encounter : null;
	}
	
	
	// method tells how many extraneous nodes to allocate additional space for
	public static int extra_space(int n) { return (n >> 1) == 0 ? n + 1 : n + (n >> 1); }
	// registers octant o in iterator and status bits
	void iteratorAdd(int o) {
		iterator = ((iterator & 0x00FFFFFF) << 3) | o | ((iterator & 0xFF000000) + 0x01000000);
		status |= 1 << (o + 3);
	}
	// following bitflags belong in status shortword, OFFSHOOT is used to flag detached branches (currently only by split2to1() concurrent method)
	// INSIDE flags leaf octant as being inside a mesh boundary, CENTROID_INSIDE flags octant's centroid as being inside a mesh boundary
	final static short IN_PROGRESS = (short)0x8000, OFFSHOOT = (short)0x4000, CULLED = (short)0x2000;
	final static int CENTROID_INTERNAL = 1<<17, INTERNAL = 0xFF|CENTROID_INTERNAL, FULL_INTERNAL = 0x7FFFFFF;
	void finish() { status &= (0xFFFF - IN_PROGRESS); }
	boolean finished() { return (status & IN_PROGRESS) == 0; }
	boolean in_progress() { return (status & IN_PROGRESS) != 0; }
	void finish_offshoot() { status &= (0xFFFF - OFFSHOOT); }
	boolean offshoot() { return (status & OFFSHOOT) != 0; }
	boolean culled() { return (status & CULLED) != 0; }
	boolean internal() { return (int_ext & INTERNAL) == INTERNAL; }		// note: only valid for leaves, ancestor can have a "finer" structure
	boolean centroid_internal() { return (int_ext & CENTROID_INTERNAL) != 0; }

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	

	void indentToLevel(StringBuilder sb) { for (int t = 0; t < level; t++) sb.append("   "); } 
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("subtree, total branches: " + branches + ", total leaves: " + leaves + "\n");
		sb.append(appendToString());
		sb.append("\n");
		return sb.toString();
	}
	
	private static int maxItemPrint = 10;
	static int match_int_ext = 0;				// gets ANDed with int_ext, if given flags are all set -> output octant (set =0 to turn off)
	static int max_level = 128;
	
	StringBuilder appendToString() {
		StringBuilder sb2 = new StringBuilder();
		
		if ((match_int_ext==0||(match_int_ext&int_ext)==match_int_ext) && level <= max_level) {
			indentToLevel(sb2);
			switch(status & 7) {
			case OCT_MMM: sb2.append("MMM: "); break; case OCT_MPM: sb2.append("MPM: "); break;
			case OCT_PMM: sb2.append("PMM: "); break; case OCT_PPM: sb2.append("PPM: "); break;
			case OCT_MMP: sb2.append("MMP: "); break; case OCT_MPP: sb2.append("MPP: "); break;
			case OCT_PMP: sb2.append("PMP: "); break; case OCT_PPP: sb2.append("PPP: "); break;
			}
			sb2.append(String.format(" x(%.3f,%.3f)y(%.3f,%.3f)z(%.3f,%.3f) ", xM,xP,yM,yP,zM,zP));
			if (in_progress() || offshoot() || (nodesX&0xFF0000)==0xFF0000)
				sb2.append((in_progress()?"IN PROGRESS ":"")+(offshoot()?"OFFSHOOT ":"")+((nodesX&0xFF0000)==0xFF0000?"iIST ":"")+"\n");
			else sb2.append("\n");
			indentToLevel(sb2);
			sb2.append((internal()?"I ":"") + (centroid_internal()?"CI ":"") + "Lvl: " + level + ", n: " + nodes + ", e: " + edges + ", f: " + facets);
			if (octant != null) {
				sb2.append(", br: " + branches + ", lv: " + leaves + ", disb: " + String.format("%.3f", disbalance));
				sb2.append(" i: " + (iterator>>24) + "|" +
						((iterator>>21)&7)+((iterator>>18)&7)+((iterator>>15)&7)+((iterator>>12)&7)+
						((iterator>>9)&7)+((iterator>>6)&7)+((iterator>>3)&7)+(iterator&7));
				sb2.append(" e: " + (enumerator>>24) + "|" +
						((enumerator>>21)&7)+((enumerator>>18)&7)+((enumerator>>15)&7)+((enumerator>>12)&7)+
						((enumerator>>9)&7)+((enumerator>>6)&7)+((enumerator>>3)&7)+(enumerator&7)+"\n");
			} else 	sb2.append("\n");
			if (int_ext != 0) {
				indentToLevel(sb2);
				sb2.append(((int_ext&1<<6)!=0? "x":".") + ((int_ext&1<<23)!=0?"x":".") + ((int_ext&1<<4)!=0? "x  ":".  "));
				sb2.append(((int_ext&1<<26)!=0?"x":".") + ((int_ext&1<<24)!=0?"x":".") + ((int_ext&1<<22)!=0?"x  ":".  "));
				sb2.append(((int_ext&1<<7)!=0? "x":".") + ((int_ext&1<<25)!=0?"x":".") + ((int_ext&1<<5)!=0? "x      Z x\n":".      Z x\n"));
				indentToLevel(sb2);
				sb2.append(((int_ext&1<<19)!=0?"x":".") + ((int_ext&1<<16)!=0?"x":".") + ((int_ext&1<<13)!=0?"x  ":".  "));
				sb2.append(((int_ext&1<<20)!=0?"x":".") + ((int_ext&1<<17)!=0?"x":".") + ((int_ext&1<<14)!=0?"x  ":".  "));
				sb2.append(((int_ext&1<<21)!=0?"x":".") + ((int_ext&1<<18)!=0?"x":".") + ((int_ext&1<<15)!=0?"x      ^/\n":".      ^/\n"));
				indentToLevel(sb2);
				sb2.append(((int_ext&1<<2)!=0? "x":".") + ((int_ext&1<<9)!=0? "x":".") + ((int_ext&1)!=0?    "x  ":".  "));
				sb2.append(((int_ext&1<<12)!=0?"x":".") + ((int_ext&1<<10)!=0?"x":".") + ((int_ext&1<<8)!=0? "x  ":".  "));
				sb2.append(((int_ext&1<<3)!=0? "x":".") + ((int_ext&1<<11)!=0?"x":".") + ((int_ext&1<<1)!=0? "x   Y<-o\n":".   Y<-o\n"));
			}
			
			if (nodeI != null && nodes > 0) {								// note: nodes can be > 0 even if node == null, for branches
				indentToLevel(sb2);
				sb2.append("n:[");
				if (maxItemPrint < nodes)
						for (int n = 0; n < maxItemPrint; n++) sb2.append(nodeI[n] + (n == maxItemPrint - 1 ? "...]\n" : ","));
				else	for (int n = 0; n < nodes; n++) sb2.append(nodeI[n] + (n == nodes - 1 ? "]\n" : ","));
			}		
			if (edgeI != null && edges > 0) {								// note: edges can be > 0 even if edge == null, for branches
				indentToLevel(sb2);
				sb2.append("e:[");
				if (maxItemPrint < edges)
						for (int e = 0; e < maxItemPrint; e++) sb2.append(edgeI[e] + (e == maxItemPrint - 1 ? "...]\n" : ","));
				else	for (int e = 0; e < edges; e++) sb2.append(edgeI[e] + (e == edges - 1 ? "]\n" : ","));
			}
			if (facetI != null && facets > 0) {								// note: facets can be > 0 even if facet == null, for branches
				indentToLevel(sb2);
				sb2.append("f:[");
				if (maxItemPrint < facets)
						for (int f = 0; f < maxItemPrint; f++) sb2.append(facetI[f] + (f == maxItemPrint - 1 ? "...]\n" : ","));
				else	for (int f = 0; f < facets; f++) sb2.append(facetI[f] + (f == facets - 1 ? "]\n" : ","));
			}
			sb2.append("\n");
		}
		if (octant != null) {
			if (octant[0] != null) sb2.append(octant[0].appendToString()); if (octant[1] != null) sb2.append(octant[1].appendToString());
			if (octant[2] != null) sb2.append(octant[2].appendToString()); if (octant[3] != null) sb2.append(octant[3].appendToString());
			if (octant[4] != null) sb2.append(octant[4].appendToString()); if (octant[5] != null) sb2.append(octant[5].appendToString());
			if (octant[6] != null) sb2.append(octant[6].appendToString()); if (octant[7] != null) sb2.append(octant[7].appendToString());
		}
		return sb2;
	}
	
	
	// method writes octree structure to OBJ file, output can specify to output only octants holding specific elements: DO_NODES/DO_EDGES/DO_FACETS
	// lvlS = from what level to start output, lvlE = at what level to end output, if both =0 the entire tree written out
	StringBuilder appendVerticesToOBJ(FEM1Octree octree, int lvlS, int lvlE, String precFormat, boolean leavesOnly, boolean internalOnly, int output) {
		StringBuilder sb2 = new StringBuilder();
		if (octree.latticeTree && level > octree.maxLevel - octree.gradations) {								// writeout of lattice tree
			double[] pCoord = {xM,yM,zM, xP,yM,zM, xM,yP,zM, xP,yP,zM, xM,yM,zP, xP,yM,zP, xM,yP,zP, xP,yP,zP};
			for (int p = 0, c = 0; p < 8; p++)
				sb2.append(	"v  "+String.format(precFormat,pCoord[c++]) + 
							" " + String.format(precFormat,pCoord[c++]) +
							" " + String.format(precFormat,pCoord[c++]) + "\n");
			octree.vOBJcnt += 8;
		} else if (	(!leavesOnly || (!octree.latticeTree&&octant==null||level==octree.maxLevel)) &&				// writeout of ordinary octree
					(!internalOnly || internal()) && (lvlS+lvlE==0||(level>=lvlS&&level<=lvlE))) {
		// DEBUG: loop exports only octants on boundary with an external centroid
			//if ((!leavesOnly || octant == null) && (!internalOnly || (boundary != 0 && !internal() && !centroid_internal()))) {
			if (output==0 || (output&DO_NODES)!=0 && nodes!=0 || (output&DO_EDGES)!=0 && edges!=0 || (output&DO_FACETS)!=0 && facets!=0) {
				double[] pCoord = {xM,yM,zM, xP,yM,zM, xM,yP,zM, xP,yP,zM, xM,yM,zP, xP,yM,zP, xM,yP,zP, xP,yP,zP};
				// DEBUG: line exports dot-like cubes to simulate centroids
//				double[] pCoord = {	xC-.003,yC-.003,zC-.003, xC+.003,yC-.003,zC-.003, xC-.003,yC+.003,zC-.003, xC+.003,yC+.003,zC-.003,
//									xC-.003,yC-.003,zC+.003, xC+.003,yC-.003,zC+.003, xC-.003,yC+.003,zC+.003, xC+.003,yC+.003,zC+.003};
				for (int p = 0, c = 0; p < 8; p++)
					sb2.append(	"v  "+String.format(precFormat,pCoord[c++]) +
								" " + String.format(precFormat,pCoord[c++]) +
								" "  +String.format(precFormat,pCoord[c++]) + "\n");
				octree.vOBJcnt += 8;
			}
		}
		if (octant != null) {
			for (int o = 0; o < 8; o++)
				if (octant[o] != null) sb2.append(octant[o].appendVerticesToOBJ(octree, lvlS, lvlE, precFormat, leavesOnly, internalOnly, output)); }
		return sb2;
	}

	StringBuilder appendPolygonsToOBJ(FEM1Octree octree, int lvlS, int lvlE, boolean leavesOnly, boolean internalOnly, int output) {
		StringBuilder sb2 = new StringBuilder();
		if (octree.latticeTree && level > octree.maxLevel - octree.gradations) {								// writeout of lattice tree
			sb2.append("f " + octree.pOBJidx + " " + (octree.pOBJidx+1) + " " + (octree.pOBJidx+5) + " " + (octree.pOBJidx+4) + "\n");
			sb2.append("f " + (octree.pOBJidx+1) + " " + (octree.pOBJidx+3) + " " + (octree.pOBJidx+7) + " " + (octree.pOBJidx+5) + "\n");
			sb2.append("f " + (octree.pOBJidx+3) + " " + (octree.pOBJidx+2) + " " + (octree.pOBJidx+6) + " " + (octree.pOBJidx+7) + "\n");
			sb2.append("f " + octree.pOBJidx + " " + (octree.pOBJidx+4) + " " + (octree.pOBJidx+6) + " " + (octree.pOBJidx+2) + "\n");
			sb2.append("f " + octree.pOBJidx + " " + (octree.pOBJidx+2) + " " + (octree.pOBJidx+3) + " " + (octree.pOBJidx+1) + "\n");
			sb2.append("f " + (octree.pOBJidx+4) + " " + (octree.pOBJidx+5) + " " + (octree.pOBJidx+7) + " " + (octree.pOBJidx+6) + "\n");
			octree.pOBJidx += 8;
		} else if (	(!leavesOnly || (!octree.latticeTree&&octant==null||level==octree.maxLevel)) &&
					(!internalOnly || internal()) && (lvlS+lvlE==0||(level>=lvlS&&level<=lvlE))) {
		// DEBUG: loop exports only octants on boundary with an external centroid
		//if ((!leavesOnly || octant == null) && (!internalOnly || (boundary != 0 &&  !internal() && !centroid_internal()))) {
			if (output==0 || (output&DO_NODES)!=0 && nodes!=0 || (output&DO_EDGES)!=0 && edges!=0 || (output&DO_FACETS)!=0 && facets!=0) {
				sb2.append("f " + octree.pOBJidx + " " + (octree.pOBJidx+1) + " " + (octree.pOBJidx+5) + " " + (octree.pOBJidx+4) + "\n");
				sb2.append("f " + (octree.pOBJidx+1) + " " + (octree.pOBJidx+3) + " " + (octree.pOBJidx+7) + " " + (octree.pOBJidx+5) + "\n");
				sb2.append("f " + (octree.pOBJidx+3) + " " + (octree.pOBJidx+2) + " " + (octree.pOBJidx+6) + " " + (octree.pOBJidx+7) + "\n");
				sb2.append("f " + octree.pOBJidx + " " + (octree.pOBJidx+4) + " " + (octree.pOBJidx+6) + " " + (octree.pOBJidx+2) + "\n");
				sb2.append("f " + octree.pOBJidx + " " + (octree.pOBJidx+2) + " " + (octree.pOBJidx+3) + " " + (octree.pOBJidx+1) + "\n");
				sb2.append("f " + (octree.pOBJidx+4) + " " + (octree.pOBJidx+5) + " " + (octree.pOBJidx+7) + " " + (octree.pOBJidx+6) + "\n");
				octree.pOBJidx += 8;
			}
		}
		if (octant != null) {
			for (int o = 0; o < 8; o++) if (octant[o] != null) sb2.append(octant[o].appendPolygonsToOBJ(octree, lvlS, lvlE, leavesOnly, internalOnly, output)); }
		return sb2;
	}


}
