package lkr74.fem1;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.atomic.AtomicInteger;

public class FEM1Octree {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			TOP LEVEL DEFINITIONS FOR A FEM1 OCTREE														//
	//			Leonard Krylov 2017																			//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////	

	FEM1 fem = null;						// this is the FEM system that spawned this octree
	public FEM1Octant root = null;			// this is the first octant of this octree
	
	final short maxLevel;					// limiter of subdivision level
	public short topLevel;					// maximal attained level
	public final static short FACTOR_NODESX = 128;
	public static final int DO_NODES = 1, DO_EDGES = 2, DO_FACETS = 4;	// flagging for processing of particular container data
	
	double[] node = null;					// TODO: inelegant hack to reuse buildOctree() with the boundary lattice cubes

	short enforcedLevel = 0;				// set > 0 to keep splitting until enforced level is attained
	float disbalanceCap = 0.90f;			// the max. level of disbalance allowed before stopping further branching
	int threadSpawnItemsCap = 70;			// number of items a branch should contain at least, to motivate spawning a thread to process it
	public int maxNodes = 60;				// the max. number of nodes per octant leaf we allow
	public int maxEdges = 30;				// the max. number of edges per octant leaf we allow
	public int maxFacets = 10;				// the max. number of facets per octant leaf we allow
	int maxFacetEncounters = 8;				// how big the coordinate array for facet encounters should be (dynamically readjusted upwards)
	static final AtomicInteger processed = new AtomicInteger(0);	// keeps count of how many octants have been processed by threads
	static final AtomicInteger offshoots = new AtomicInteger(0);	// keeps count of how many offshoots were generated
	int enumerator = 0;						// keeps track of value during enumeration of octants
	int[] levelSums = null;					// octant count sums per level, incremented by instantiators or by process() method
	AtomicInteger[] levelSumsT = null;		// incremented by threaded methods that must track how many octants were processed in a certain level
	boolean splitMade = false;				// for the threadbased 2:1 splitter to flag that more method calls are needed
	boolean latticeTree = false;			// is this a boundary lattice tree?
	int latticeSubdivs = 1;					// how many final cubes per ordinate the mesh volume was subdivided into (defaults to minLatticeSubdivs)
	int gradations = 4;						// specifies how many gradations of level of detail that is requested for the inner volume
	double subdivDtreeWidth = 0;

	// dynamic allocation readjuster for leavesOverlappingBBox() method
	int oOBB_lAllocs = 0, oOBB_lAllocs_fail = 0, oOBB_allocSize = 4;
	float oOBB_lAllocs_successRate = 0;
	// gets relative timings of the three mesh operations (divided by the slowest one), for thread workload balancing
	// note: this is probably irrelevant unless a mechanism for enforcing start of the three node/edge/facet ops all at once is created
	float splitNtiming = 1, splitEtiming = 1, splitFtiming = 1;
	int threadSplitF = 28;											// how many items needed at least to multithread splitting ops
	float splitNf = 29.2f, splitEf = 1, splitFf = 1.59f;			// the factors that modulate threadSplitF
	
	int vOBJcnt = 0, pOBJidx = 1;									// running vertex count & poly.index during OBJ export of octree visualisation

	static int DEBUG_LEVEL;											// debug level at first gotten from the FEM1 object

	// instantiates a skeleton tree with default settings
	public FEM1Octree(FEM1 fem, int doWhat, int maxLevel) {
		this.fem = fem; 
		node = fem.node;
		DEBUG_LEVEL = FEM1.DEBUG_LEVEL;
		levelSums = new int[maxLevel + 1];
		levelSums[0] = 1;
		root = new FEM1Octant(this, doWhat);
		if ((doWhat&DO_FACETS) != 0) {
			int maxBin = fem.edgeSpread.maxBin();					// find out what two groups of facet sizes are the most numerous
			// find out what octant size that groups size will motivate (average size of the two groups)
			maxLevel = maxLevelFromGranularity((fem.edgeSpread.div[maxBin+1] + fem.edgeSpread.div[maxBin]) / 2);
			int divisions = 1; for (int m = maxLevel; m-- > 0; divisions *= 2);
			while (divisions < fem.minLatticeSubdivs) { divisions *=2; maxLevel++; }
		}
		this.maxLevel = (short)maxLevel;
	}
	
	// instantiates first octant of subdivision of a lattice tree (finds bounding box and sets up node reference array)
	public FEM1Octree(FEM1Octree sourceTree, int gradations) {
		fem = sourceTree.fem;
		DEBUG_LEVEL = FEM1.DEBUG_LEVEL;
		latticeTree = true;
		this.gradations = gradations;
		node = sourceTree.fem.bLatticeNode;
		disbalanceCap = 1;					// do not stop subdivision until attaining the lattice cubes as leaves
		maxNodes = 1;						// we want subdivision to one lattice cube per leaf
		// since arbitrary non-power-of-2 subdivs can be specified, we need to expand lattice tree to a power-of-2 subdivision & scale
		// from the subdivision & scale we're given, fitting that subdivision (+1 extra lattice cube layer per side) inside a power-of-2 scale
		while (latticeSubdivs < fem.latticeSubdivs + 2) latticeSubdivs *= 2;
		FEM1Octant stRoot = sourceTree.root;
		double lStep = fem.latticeStep, extent = fem.latticeStep * latticeSubdivs;
		// extend the lattice tree from source tree's zero point (+ space for 1 lattice layer per ordinate) to power-of-2 coverage of
		// the requested lattice subdivision (+ space for 1 lattice layer per ordinate)
		root = new FEM1Octant(sourceTree,	stRoot.xM-lStep, stRoot.yM-lStep, stRoot.zM-lStep,
											stRoot.xM-lStep+extent, stRoot.yM-lStep+extent, stRoot.zM-lStep+extent);
		subdivDtreeWidth = (double)latticeSubdivs / (root.xP - root.xM);	// useful for calculating integer coordinates from an octant
		maxLevel = enforcedLevel = maxLevelFromGranularity(fem.latticeStep + FEM1Octant.OCT_MARGIN);
		levelSums = new int[maxLevel + 1];
		levelSums[0] = 1;
		levelSumsT = new AtomicInteger[maxLevel + 1];
		for (int c = 0; c <= maxLevel; c++) levelSumsT[c] = new AtomicInteger(0);
	}
	
	// user wants a maximal subdivision granularity, method calculates what branch level to stop subdividing at to reach that granularity
	// assumingly the return value will be assigned to maxLevel parameter
	public short maxLevelFromGranularity(double granularity) {
		double subDim = root.xP - root.xM; int level = 0;
		while (subDim > granularity) { level++; subDim *= .5; }
		// find out if final or previous subdivision level lies closer to granularity
		return (granularity - subDim) < (subDim * 2 - granularity) ? (short)level : (short)(level - 1);
	}

	// method calculates what granularity is attainable from current maxLevel subdivision limit, starting at world bounding box
	public double granularityFromMaxLevel() {
		double granularity = root.xP - root.xM; int level = maxLevel;
		while (level-- > 0) { granularity *= .5; }
		return granularity;
	}
	
	
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(latticeTree ? "Lattice tree" : "Octree");
		sb.append(", total branches: " + root.branches + ", total leaves: " + root.leaves + ", top level: " + topLevel + "\n");
		sb.append(root.appendToString(topLevel));
		sb.append("\n");
		return sb.toString();
	}

	
	public void toOBJ(int lvlS, int lvlE, boolean leavesOnly, boolean internal, int output) {
		File file = new File(fem.name + "_octree" + ((output&DO_NODES)!=0?"N":"") + ((output&DO_EDGES)!=0?"E":"") + ((output&DO_FACETS)!=0?"F":"") + ".obj");
		if (!file.exists()) { try {	file.createNewFile(); } catch (IOException e) { e.printStackTrace(); }}
		BufferedWriter bw = null;
		try {		bw = new BufferedWriter(new FileWriter(file));
		} catch (IOException e) { e.printStackTrace(); }
		StringBuilder sb = new StringBuilder();
		String precFormat = "%." + fem.precision_OBJ + "f";
		vOBJcnt = 0; pOBJidx = 1;
		if (lvlE > topLevel) lvlE = topLevel;
		if (lvlS > lvlE) lvlS = lvlE;
		
		sb.append("# Visualisation of octree\n#\n");
		sb.append(root.appendVerticesToOBJ(this, lvlS, lvlE, precFormat, leavesOnly, internal, output));
		sb.append("# " + vOBJcnt + " vertices\n\ng " + fem.name + "_octree\n");
		sb.append(root.appendPolygonsToOBJ(this, lvlS, lvlE, leavesOnly, internal, output));
		sb.append("# " + (pOBJidx - 1) / 6 + " faces\n\ng\n");
		
		try {	bw.write(sb.toString());									// write out & close file
				bw.flush(); bw.close();
		} catch (IOException e) { e.printStackTrace(); }
	}
	
	public void toOBJ(boolean leavesOnly, boolean internal, int output) { toOBJ(0, 0, leavesOnly, internal, output); }

}
