package lkr74.fem1;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import lkr74.matrixlib.BinBitImage;
import lkr74.matrixlib.Matrix;
import lkr74.matrixlib.NSPMatrix;
import lkr74.matrixlib.NspArray;
import lkr74.matrixlib.NspNode;

public class FEM1 {

	// note: for more programmer-friendly access, methods for accessing the following 1D-arrays
	// can be implemented, to reference specific coordinates and elements
	String name;
	// node, element, freed element, polygon & smooth patch counts for this FEM1 object
	public int nodes=0, nodesWork=0, deletedNodes=0, elements=0, elements2=0, elements2Work=0, elements2Free=0, polygons=0, patches=0;
	int precision_OBJ = 5;							// decimal places of the inout OBJ file
	
	public double[] node = null;					// node coordinates come in (x,y,z)-triples within this 1D-array
	public double[] nodeWork = null;				// the work array for newly added nodes is also processed, and is reintegrated after reaching some length
	byte nodeFlag[]=null, nodeWorkFlag[]=null;		// up to 8 multipurpose flags for each node
	int[] deletedNode=null;							// stack of previously deleted nodes
	double[][] nodeCoordSort = null;				// arrays holding coordinates sorted by (x, y, z)
	int[][] nodeCoordSortI = null;					// reindexing tables of the sorted coordinates
	
	int[] element = null;							// each element's node indexes come here (ex. tetrahedral indexes queued in a,b,c,d quartets)
	byte[] nodeCount = null;						// each element's node count is also specifying the type of element

	int[] polygon = null;							// polygons come in N-tuplets, N = no. of nodes (ex. faces come in a,b,c triplets)
	int[] polygonOffset = null;						// polygon node offsets work both as specifiers of poly node count and give position of polygon in array
	int[] polygonEdgeIndex = null;					// indexes matching each shared edge to particular polygons
	double[] polygonEdge = null;					// array specifying edge length of every polygonal edge
	int[] polygonEdgeN = null;						// convenient back-reference from edge index to the underlying nodes
	int polygonEdges;								// total nonredundant count of polygon edges
	//int[] polygonPatch = null;					// every polygon's patch membership is indexed here

	int[] patch = null;								// the patch/smoothing group array holds offsets into face array for the start of every patch's faces
	
	FEM1Element[] element2 = null;					// object-oriented element types go into this structure
	FEM1Element[] element2Work = null;				// the work array for newly added elements is also processed, and is merged after reaching some length
	int[] element2Free = null;
	
	int[][] elementSupports = null;					// indexes every node's supported elements
	int elementSupportMaxArrayL = 0;				// maximum acquired length of nodal element support array for temp.array worst-case allocation
	int[][] elementNeighbours = null;				// indexes neighbours of every element (this array is constructed from elementSupports[])
	int elementNeighboursMaxArray = 0;				// holds largest neighbourhood for temp.array worst-case allocation

	int[][] polygonSupports = null;					// indexes every node's supported polygons
	int polygonSupportMaxArrayL = 0;				// maximum acquired length of nodal polygon support array
	int[][] polygonNeighbours = null;				// indexes neighbours of every polygon (this array is constructed from polygonSupports[])
	int[] borderPolygons = null;					// array optimises iteration over polygons bordering a patch by keeping an indexing over them
	int[][] patchBorders = null;					// every patch can have more than one edge border cycle, cycles lie sequentially in each patch's subarray

	int nodeBitSets = 0;
	long[] nodeChecks = null;						// bitflagging for nodes (used in conjunction with node processing)
	int[] encapCycle = null;
	
	public double[] bBox = null;					// bounding box of dataset kept here

	// parameters pertaining to tetraherdonisation lie here
	double maxTsize = 1;							// maximal bound on size of any constructed tetrahedron
	public int octreeMaxItems = 60;					// the max. number of items per octant leaf we allow
	public int nodeworkFactor = 4;					// let nodeWork[] array be maxium 1/4th of node[] array
	
	// the datatype requested from the object instancer
	public static final int MESH_HANDMADE = 0, MESH_HANDMADE_OBJECTIFIED = 1, MESH_PSC = 2;
	// the states for constructing tetra-elements during file reading
	private static final int TETRA_NEW = 0, TETRA_SEEK = 1, TETRA_CHECK_UNMATCHED = 2, TETRA_SEEK_UNMATCHED = 3;
	// these are the offsets into pertinent fieldvalues of tetrahedral precalculated data
	static final int EDGE01 = 0, EDGE02 = 1, EDGE03 = 2, EDGE12 = 3, EDGE13 = 4, EDGE23 = 5;
	static final int FACE012 = 6, FACE023 = 7, FACE031 = 8, FACE132 = 9;
	
	
	
	static ExecutorService executor;				// thread pool for FEM ops
	static int taskNum, procNum = 1;
	static List<Future<Double>> futureDoubleList;	// list of return doubles from Callable threads
	static Thread taskList[];
	
	// initialise Executor (left to the user, since one can't get Runtime's processor count during class loading)
	public static void initExecutor(int processors) {
        // define the Callable executor pool of reusable threads
		if (processors == 0) 	executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		else					executor = Executors.newFixedThreadPool(processors);
		procNum = (processors == 0 ? Runtime.getRuntime().availableProcessors() : processors);
  		// define a Future list of Callable return values constituting doubles (it will expand dynamically)
		futureDoubleList = new ArrayList<Future<Double>>();	
	}

	// initialise Executor (left to the user, since one can't get Runtime's processor count during class loading)
	public static void initTaskList(int processors) {
        // define the Callable executor pool of reusable threads
		if (processors == 0) 	taskList = new Thread[Runtime.getRuntime().availableProcessors()];
		else					taskList = new Thread[processors];
		procNum = (processors == 0 ? Runtime.getRuntime().availableProcessors() : processors);
		taskNum = 0;
	}

	public static void finishTaskList() {
		for (Thread task: taskList) {
			if (task == null) continue;
			try { task.join();
			} catch (InterruptedException e) { e.printStackTrace(); }
		}
		taskNum = 0;
	}

	// skeleton FEM solution instantiator
	public FEM1(String name) { this.name = name; }
	
	// method instantiates a single FEM1 object from a Lightwave OBJ file, from current buffer position
	// note: this is code accepting a ready-made object subdivided into tetrahedra purely for small handmade test meshes
	// the constraint on the OBJ file: the OBJ face indexes must come in groups of fours, specifying the 4 faces of a tetrahedron
	// if OBJ file is exported from a 3D program, this means that the object in the program must consist of individual tetrahedrons,
	// each tetrahedral having it's own four unique vertexes (meaning, each tetrahedron is a distinct, detached element within the object)
	// since this means most vertexes will be duplicates, a post-optimisation has to be performed to make the tetrahedra share nodes
	public FEM1(BufferedReader br, int dataType) {
				
		int state = 1, tetraCnt = 0, polyNodes = 0;
		boolean keepReading = true, gotNormals = false;
		String s = "";
		String[] dataRow = null;
		
		try {
			while (keepReading)
				switch (state) {
				// state 1, scanning from current position to count the incoming data, mark it to return for reading the data
				case 1:
					br.mark(10000000);
					state = 2; break;

				// state 2 skips the header and anything after the last face index row (for consecutive objects)
				case 2:
					if ((s = br.readLine()) == null) { state = 0; break; }
					if (s.startsWith("v ")) {
						nodes++;
						String vS = s.split(" ")[2];
						double vD = Double.valueOf(vS);									// find out the precision of data
						int integrals = 1;
						while (vD >= 10.0 || vD <= -10) { vD /= 10.0; integrals++; }
						precision_OBJ = vS.length() - (vD < 0 ? 1 : 0) - integrals - 1;
						state = 3; break; }	// if first vertex found, progress to vertex counter
					break;
					
				// state 3 counts the incoming vertexes
				case 3:
					if (br.readLine().startsWith("v ")) nodes++;
					else { state = 4; break; }	// if non-vertex row found, progress
					break;
					
				// state 4 skips the gap between vertexes & polygons
				case 4:
					//if ((s = br.readLine()) == null) { state = 0; break; }
					s = br.readLine();
					if (s.startsWith("s ")) {
						int patchNr = Integer.valueOf(s.split(" ")[1]);						// find out the highest patch/smooth group index
						if (patchNr > patches) patches = patchNr;
					}
					else if (s.startsWith("vn ")) gotNormals = true;
					else if (s.startsWith("f ")) {
						tetraCnt++; polygons++;
						int[] poly = decodePolygon(s, gotNormals);
						polyNodes += poly.length;
						state = 5; break; }													// if first poly found, progress to vertex counter
					break;
					
				// state 5 will increment element count for every 4 polygons read in (which puts a constraint on how the indexes are ordered)
				case 5:
					// note: if element count was zero, do not progress to proper reading in of the data
					//if ((s = br.readLine()) == null) { state = elements > 0 ? 10 : 0; break; }
					s = br.readLine();
					if (s.startsWith("s ")) {
						int patchNr = Integer.valueOf(s.split(" ")[1]);						// find out the highest patch/smooth group index
						if (patchNr > patches) patches = patchNr;
					}
					else if (s.startsWith("f ")) {
						tetraCnt++; polygons++;
						int[] poly = decodePolygon(s, gotNormals);
						polyNodes += poly.length; }
					else { state = 10; break; }
					if ((dataType == MESH_HANDMADE || dataType == MESH_HANDMADE_OBJECTIFIED) && tetraCnt >= 4) { elements++; tetraCnt = 0; }
					break;
					
				// state 10 restarts file reading, initialises data arrays, reads in all data in one state
				case 10:
					br.reset();
					
					node = new double[nodes * 3];
					byte[] nodeCheck = new byte[nodes];
					if (dataType != MESH_PSC) {									// for PSC meshes, elements are created in separate method
						element = new int[elements * 4];
						nodeCount = new byte[elements];
					}
					polygon = new int[polyNodes];
					polygonOffset = new int[polygons + 1];
					//polygonPatch = new int[polygons];
					patch = new int[patches == 0 ? 2 : patches + 1];
					int[] polyUnmatched = new int[polygons];
					
					do { s = br.readLine();	} while (!s.startsWith("v "));
					dataRow = s.split(" ");
					int n3 = 0;
					node[n3++] = Double.valueOf(dataRow[2]);
					node[n3++] = Double.valueOf(dataRow[3]);
					node[n3++] = Double.valueOf(dataRow[4]);
					
					for (int n = 1; n < nodes; n++) {
						dataRow = br.readLine().split(" ");
						node[n3++] = Double.valueOf(dataRow[2]); node[n3++] = Double.valueOf(dataRow[3]); node[n3++] = Double.valueOf(dataRow[4]);
					}
					// find out the name of the 3D-object and acquire it
					do { s = br.readLine();	} while (!s.startsWith("g"));
					name = s.split(" ")[1];
					
					if (dataType == MESH_PSC) {
						int pch = 0, pC = 0;
						for (int f = 0, sGroup = 0; f < polygons;) {
							
							s = br.readLine();
							if (s.startsWith("usemtl")) continue;
							int new_sGroup = catchSmoothingGroup(s); 		// if caught a smoothgroup separator, register facecount for that patch
							if (new_sGroup >= 0) {							// TODO: DEBUG: the patch index from OBJ file is not utilised
								sGroup = new_sGroup;						// update smooth group
								patch[pch + 1] = patch[pch++];				// move on polygon offset
								continue;
							}
							int[] poly = decodePolygon(s, gotNormals);
							
							polygonOffset[f++] = pC;						// note: polygonOffset[offset + 1] - polygonOffset[offset] = polyNodeCount
							for (int ni1 = 0; ni1 < poly.length; ni1++) polygon[pC++] = poly[ni1];

							//polygonPatch[f] = sGroup;						// record face's patch/smooth group
							patch[pch]++;									// increment patch's polycount by incrementing offset
						}
						polygonOffset[polygons] = pC;						// not to forget the final polygon's node count
						
						keepReading = false; break;							// current PSC object is read, we're done with the file
					}
					
					int buildStatus = TETRA_NEW;
					int f = 0, fU = 0, fUseek = 0, a=0, b=0, c=0, d=0, a2 = -1, b2 = -1, c2 = -1;
					
					for (int e = 0, f3 = 0, e4 = 0; e < elements || f < polygons;) {
						
						int buildStatusPrev = buildStatus, faceUIdx;
						int[] poly;
						
						switch (buildStatus) {
						// case loads in the base face of a new tetrahedron
						case TETRA_NEW:
							// take care of situation when a matching face for a tetrahedron can't be found
							if (f >= polygons) throw new InvalidParameterException("FEM1(): Element " + e + " lacks faces.");
							s = br.readLine();
							if (s.startsWith("usemtl") || s.startsWith("s ")) s = br.readLine();		// we're not interested in smoothing groups for handmade meshes
							poly = decodePolygon(s, gotNormals);	
							polygonOffset[f++] = f3;
							polygon[f3++] = a = poly[0];					// record the incoming face					
							polygon[f3++] = b = poly[1];					// only faces accepted for handmade objects, not full polygons
							polygon[f3++] = c = poly[2]; d = -1;
							if (nodeCheck[a] > 0) continue;					// if node has been registered, we've already constructed the element
							if (fU > 0) { buildStatus = TETRA_CHECK_UNMATCHED; fUseek = 0; continue; }
							buildStatus = TETRA_SEEK;

						// case seeks for matching face and 4th node within datafile, storing incoming faces
						case TETRA_SEEK:
							if (f >= polygons) throw new InvalidParameterException("FEM1(): Element " + e + " lacks faces.");
							// we have first 3 nodes of a tetrahedron, consecutive faces must be checked for the fourth node
							s = br.readLine();
							if (s.startsWith("usemtl") || s.startsWith("s ")) s = br.readLine();		// we're not interested in smoothing groups for handmade meshes
							poly = decodePolygon(s, gotNormals);
							polygonOffset[f++] = f3;
							polygon[f3++] = a2 = poly[0];					// only faces accepted for handmade objects, not full polygons
							polygon[f3++] = b2 = poly[1];
							polygon[f3++] = c2 = poly[2];
							if (nodeCheck[a2] > 0) continue;				// if node has been registered, we've already constructed the element
							break;

						// case looks through list of yet unmatched faces to find 4th node for current base face
						case TETRA_CHECK_UNMATCHED:
							faceUIdx = polyUnmatched[fUseek++] * 3;
							a2 = polygon[faceUIdx++]; b2 = polygon[faceUIdx++]; c2 = polygon[faceUIdx]; d = -1;
							if (fUseek >= fU) buildStatus = TETRA_SEEK;
							break;

						// case takes a face from unmatched list for a base face of next tetrahedron
						case TETRA_SEEK_UNMATCHED:
							faceUIdx = polyUnmatched[fUseek = 0] * 3;
							polyUnmatched[fUseek++] = polyUnmatched[--fU];
							a = polygon[faceUIdx++]; b = polygon[faceUIdx++]; c = polygon[faceUIdx]; d = -1;
							buildStatus = fUseek >= fU ? TETRA_SEEK : TETRA_CHECK_UNMATCHED;
							continue;					
						}

						if (a==a2) {		if (b==c2) d = b2; else if (c==b2) d = c2; }
						else if (a==b2) {	if (b==a2) d = c2; else if (c==c2) d = a2; }							
						else if (a==c2) {	if (b==b2) d = a2; else if (c==a2) d = b2; }
						else if (b==b2) {	if (c==a2) d = c2; else if (a==c2) d = a2; }
						else if (b==c2 && c==b2) d = a2;
						else if (c==c2 && b==a2) d = b2;
						
						if (d != -1) {								// if fourth node found
							element[e4++] = a; element[e4++] = b; element[e4++] = c; element[e4++] = d;
							nodeCheck[a]++; nodeCheck[b]++; nodeCheck[c]++; nodeCheck[d]++;
							// if unmatched faces exist, start new tetra. construction from one of them
							buildStatus = fU > 0 ? TETRA_SEEK_UNMATCHED : TETRA_NEW;
							nodeCount[e++] = 4; continue;			// a tetrahedron's node count = 4
						}
						
						if (buildStatusPrev == TETRA_NEW || buildStatusPrev == TETRA_SEEK)
							polyUnmatched[fU++] = f - 1;			// a file-loaded facet was unmatched, put it in list of unmatched faces
					}
					polygonOffset[polygons] = polygonOffset[polygons - 1] + 3;	// not to forget the final polygon's node count
				
				// state 0 reached when EOF was reached, quit while loop
				case 0:
					keepReading = false; break;
				}
		} catch (IOException e) { e.printStackTrace(); }
		
		
		// if object-oriented element allocation was requested, reconfigure dataset
		if (dataType == MESH_HANDMADE_OBJECTIFIED) {
			// collapse nodes that are identical or close enough (since every loaded tetrahedron touches others perfectly)
			uniqueNodes(0.001);

			element2 = new FEM1Element[elements2 = elements];
			for (int e = 0, e4 = 0; e < elements; e++) {
				element2[e] = new FEM1Element(this, element[e4++], element[e4++], element[e4++], element[e4++], null);
			}
			//element = null;
			
			// find datasharing neighbour sets of every tetrahedral element, face-to-face and edge-to-edge
			encapCycle = new int[64];
			elementsNodeSupport();
			tetraNeighbours(null, true);
		}
		else if (dataType == MESH_PSC) {
			polygonsNodeSupport();
			facetNeighbours(false, false);			// find each face's neighbours, accept only facets, skip corner neighbours
			facetEdgeLengths();
			//creaseBoundariesOfFacets();
		}
	}
	
	
	// method allows a seamless indexing of both operative (element2[]) and constructive (element2Work[]) arrays
	public FEM1Element getElement2(int e) { return (e >= element2.length ? element2Work[e - element2.length] : element2[e]); }
	// method occupies a known location with supplied element object
	
	public void setElement2(int e, FEM1Element elem) {
		if (e > element2.length) element2Work[e - element2.length] = elem;
		else  element2[e] = elem;
	}
	public void deleteElement(int e) {
		setElement2(e, null);
		if (element2Free == null) element2Free = new int[elements2];		// if array of free element slots is uninitialised
		element2Free[elements2Free++] = e;
	}
	
	// method will find a free element slot itself, or reallocate work array, first doing a merge if necessary
	public void putElement2(FEM1Element elem) {
		if (elements2Free > 0)
			setElement2(element2Free[--elements2Free], elem);
		else {
			if (element2Work.length >= (element2.length >> 1)) mergeElement2Work();
			if (element2Work == null) {
				element2Work = new FEM1Element[16];
			} else {
				FEM1Element[] element2WorkNew = new FEM1Element[element2Work.length + 16];
				int e = 0, e2WorkLEven = element2Work.length & 0xFFFFFFFE;	// loop unroll test (even-ising count and doing two ops per iteration)
				for (; e < e2WorkLEven;) { element2WorkNew[e] = element2Work[e++]; element2WorkNew[e] = element2Work[e++]; }
				if (element2Work.length > e2WorkLEven) element2WorkNew[e] = element2Work[e];
				element2WorkNew = element2Work;
			}
			elements2 += 16;
			element2Free = new int[elements2];			// buffer increase means buffers were full means no free elements, so we can simply reallocate
			for (int f = 0, f2 = elements2 - 1; f < 16; f++, f2--) { element2Free[f] = f2; element2Free[f] = f2; }
		}
	}
			
	void mergeElement2Work() {
		if (element2Work == null) return;
		FEM1Element[] element2New = new FEM1Element[elements2];
		int e = 0, eEnd = element2.length, eEndEven = eEnd & 0xFFFFFFFE;
		while (e < eEndEven) { element2New[e] = element2[e++]; element2New[e] = element2[e++]; }
		if (eEnd > eEndEven) element2New[e] = element2[e++];
		eEndEven = elements2Work & 0xFFFFFFFE;
		for (int eW = 0; e < eEndEven;) { element2New[e++] = element2Work[eW++]; element2New[e++] = element2Work[eW++]; }
		if (elements2Work > eEndEven) element2New[e] = element2Work[elements2Work - 1];
		element2 = element2New; element2Work = null;
	}
	
	
	// method adds new node, returns it's index
	public int addNode(double x, double y, double z, byte flag) {
		int n;
		// the case of starting a new node collection
		if (node == null) {
			node = new double[(nodes = 64) * NCOORD]; nodeWork  = new double[(n = nodes / nodeworkFactor) * NCOORD];
			node[0] = x; node[1] = y; node[2] = z;
			deletedNode = new int[nodes + n]; deletedNodes = 63;
			for (int nD = 1; nD < nodes; nD++) deletedNode[nD - 1] = nodes - nD;
			return 0;
		}
		// the case of having freed nodes
		if (deletedNodes > 0) {
			int n3 = deletedNode[--deletedNodes] * NCOORD;
			if (n3 >= node.length) {
				n3 -= node.length; n = n3 / 3; nodeWork[n3++] = x; nodeWork[n3++] = y; nodeWork[n3] = z; nodesWork++;
				if (nodeFlag != null) nodeWorkFlag[n] = flag;
			} else {
				n = n3 / NCOORD; node[n3++] = x; node[n3++] = y; node[n3] = z;
				if (nodeFlag != null) nodeFlag[n] = flag; }
			return n;
		}
		// the case of adding to nodeWork[] or increasing it (nodesWork counter can be less than nodeWork[] size)
		n = nodes;
		int n3 = nodesWork * NCOORD;
		if (nodeWork != null && n3 < nodeWork.length) {
			nodeWork[n3++] = x; nodeWork[n3++] = y; nodeWork[n3] = z; nodes++; nodesWork++;
			if (nodeFlag != null) nodeWorkFlag[nodesWork] = flag;
		} else {														// if we filled nodeWork[] array
			mergeNodeWork();											// merge it into node[] array
			int nodesWork2 = nodes / nodeworkFactor;
			if (nodesWork2 < 8) nodesWork2 = 8;
			nodeWork = new double[nodesWork2 * NCOORD];					// create new nodeWork[] array
			nodeWork[0] = x; nodeWork[1] = y; nodeWork[2] = z;
			if (nodeFlag != null) nodeWorkFlag[0] = flag;
			nodesWork = 1;
			deletedNode = new int[nodes++ + nodesWork2];				// assign new array for deleted nodes
		}
		return n;
	}
	
	// since elements refer to nodes by index and reindexing elements is an extra expense, this class will maintain a node's index
	// up until it's complete deletion (only allowed to happen if no element refers to it anymore) and put new nodes in place of deleted ones
	public void deleteNode(int n) {
		if (node == null) return;
		int n3 = n * NCOORD;
		if (n3 < node.length) {
					if (node[n3] == Double.MAX_VALUE) return;						// if node already deleted, do nothing more
					node[n3++] = Double.MAX_VALUE; node[n3++] = Double.MAX_VALUE; node[n3] = Double.MAX_VALUE;
		} else {	n3 -= node.length;
					if (nodeWork[n3] == Double.MAX_VALUE) return;					// if node already deleted, do nothing more
					nodeWork[n3++] = Double.MAX_VALUE; nodeWork[n3++] = Double.MAX_VALUE; nodeWork[n3] = Double.MAX_VALUE; }
		if (deletedNodes < deletedNode.length) deletedNode[deletedNodes++] = n;
	}
	
	void mergeNodeWork() {
		if (nodeWork == null) return;
		int nodesTot3 = nodes * 3;
		double[] nodeNew = new double[nodesTot3];
		byte[] nodeFlagNew = nodeFlag != null ? new byte[nodes] : null;
		int n = 0, nW = 0, n3 = 0, n3W = 0, nEnd = nodes - nodesWork;
		if (nodeFlag != null)
		while (n < nEnd) {		nodeNew[n3] = node[n3++]; nodeNew[n3] = node[n3++]; nodeNew[n3] = node[n3++];
								nodeFlagNew[n] = nodeFlag[n++]; }
		else while (n < nEnd) {	nodeNew[n3] = node[n3++]; nodeNew[n3] = node[n3++]; nodeNew[n3] = node[n3++]; n++; }
		nEnd += nodesWork;
		if (nodeFlag != null)
			while (n < nEnd) {	nodeNew[n3++] = nodeWork[n3W++]; nodeNew[n3++] = nodeWork[n3W++]; nodeNew[n3++] = nodeWork[n3W++];
								nodeFlagNew[n++] = nodeWorkFlag[nW++]; }
		else while (n < nEnd) {	nodeNew[n3++] = nodeWork[n3W++]; nodeNew[n3++] = nodeWork[n3W++]; nodeNew[n3++] = nodeWork[n3W++]; n++; }
		node = nodeNew; nodeWork = null; nodesWork = 0;
		nodeFlag = nodeFlagNew; nodeWorkFlag = null;
	}
	
	
	
	// method eliminates duplicate vertexes lying in same coordinates, found by x/y/z-delta distance error given by "precision"
	// note: method is for INITIALISATION, solely used to eliminate vertex duplicates in handcrafted tetrahedral OBJ files
	private void uniqueNodes(double precision) {
		
		double[] node2 = new double[nodes * 3];				// node2[] contains node coordinates moved according to the remapping
		int nodes2 = 0, nAcc2 = 0;
		int[] nodeRef = new int[nodes];						// nodeRef[] refers every old node index to it's new node index
		for (int n = 0; n < nodes; n++) nodeRef[n] = -1;

		// if elements are unitialised, we only have a shell mesh to operate on
		if (element == null) {
			int[] polygonNew = new int[polygonOffset[polygons]];
			
			for (int p = 0; p < polygons; p++) {			// iterate over current N-tuple of polygon nodes within global polygon array
				for (int pN = polygonOffset[p], pNo = polygonOffset[p + 1]; pN < pNo; pN++) {

					int pIdx = polygon[pN];
					if (nodeRef[pIdx] != -1) {				// if tested node ref. is already remapped, refer to it again
						polygonNew[pN] = nodeRef[pIdx];
						continue;
					}

					int n = pIdx * 3;
					double x1 = node[n++], y1 = node[n++], z1 = node[n++];
					boolean foundNode = false;
					
					int n2 = 0;
					for (int n3 = 0; n2 < nodes2; n2++) {
						double dlt1 = x1 - node2[n3]; if ((dlt1<0?-dlt1:dlt1) > precision) { n3 += 3; continue; } else n3++;
						double dlt2 = y1 - node2[n3]; if ((dlt2<0?-dlt2:dlt2) > precision) { n3 += 2; continue; } else n3++;
						double dlt3 = z1 - node2[n3++]; if ((dlt3<0?-dlt3:dlt3) > precision) continue;
						foundNode = true; break;
					}
					if (!foundNode) {						// if no matching node existed in accumulation array node2[], add this node to array
						polygonNew[pN] = nodes2;			// the current new node ref. in polygonNew[] replacing the old ref. in polygon[]
						nodeRef[pIdx] = nodes2;				// and nodeRef[] holds the mapping from old node ref. to new node ref.
						node2[nAcc2++] = x1; node2[nAcc2++] = y1; node2[nAcc2++] = z1;
						nodes2++;
					} else {
						polygonNew[pN] = n2;				// copy found new node ref. into elementNew[]
						nodeRef[pIdx] = n2;	
					}
				}
			}
		} else {

			// this part operates directly on tetrahedral elements
			// TODO: two close nodes set equal could collapse an element to a zero volume, construct a checker that
			// either keeps such nodes intact, or eliminates zero-volumed elements altogether
			int elements4 = elements * 4;
			int[] elementNew = new int[elements4];
					
			for (int e4 = 0; e4 < elements4; e4++) {
				
				int eIdx = element[e4];
				if (nodeRef[eIdx] != -1) {						// if tested node ref. is already remapped, refer to it again
					elementNew[e4] = nodeRef[eIdx];
					continue;
				}
				
				int n = eIdx * 3;
				double x1 = node[n++], y1 = node[n++], z1 = node[n++];
				boolean foundNode = false;
				
				int n2 = 0;
				for (int n3 = 0; n2 < nodes2; n2++) {
					double dlt1 = x1 - node2[n3]; if ((dlt1<0?-dlt1:dlt1) > precision) { n3 += 3; continue; } else n3++;
					double dlt2 = y1 - node2[n3]; if ((dlt2<0?-dlt2:dlt2) > precision) { n3 += 2; continue; } else n3++;
					double dlt3 = z1 - node2[n3++]; if ((dlt3<0?-dlt3:dlt3) > precision) continue;
					foundNode = true; break;
				}
				if (!foundNode) {						// if no matching node existed in accumulation array node2[], add this node to array
					elementNew[e4] = nodes2;			// the current new node ref. in elementNew[] replacing the old ref. in element[]
					nodeRef[eIdx] = nodes2;				// and nodeRef[] holds the mapping from old node ref. to new node ref.
					node2[nAcc2++] = x1; node2[nAcc2++] = y1; node2[nAcc2++] = z1;
					nodes2++;
				} else {
					elementNew[e4] = n2;				// copy found new node ref. into elementNew[]
					nodeRef[eIdx] = n2;	
				}
			}
					
			// elementNew contains reassigned references
			element = elementNew;
		}

		// reassign all facet node references according to old->new node ref. mapping
		for (int p = 0, p3 = 0; p < polygons; p++) {
			int pN = polygonOffset[p], pNo = polygonOffset[p + 1];
			if (pNo - pN == 3) {												// facet case
				polygon[p3] = nodeRef[polygon[p3++]];
				polygon[p3] = nodeRef[polygon[p3++]];
				polygon[p3] = nodeRef[polygon[p3++]];
			} else
				for (; pN < pNo; pN++) polygon[pN] = nodeRef[polygon[pN]];		// generic polygon case
		}

		if (nodes2 < nodes) {
			node = new double[nodes2 * 3];				// rescale node array to the reduced array size
			for (int n = 0, nEnd = nodes2 * 3; n < nEnd;) { node[n] = node2[n++]; node[n] = node2[n++]; node[n] = node2[n++]; }
			nodes = nodes2;
		} else
			node = node2;
	}
	
	
	
	// method generates for every node the array of elements that are supported by that node
	// this is an INITIALISATION method, while during updates localised methods are better used
	void elementsNodeSupport() {
		
		int elements4 = elements2 * 4;
		elementSupports = new int[nodes][];
		for (int e4 = 0; e4 < elements4; e4++) {
			int e = e4 >> 2, n = getElement2(e).nodeRef[e4 & 3];
			int[] support = elementSupports[n];
			
			if (support == null) {									// optimally, every node is shared by 20 tetrahedra (+ 1 array length integer)
				support = elementSupports[n] = new int[31];			// but as an optimal mesh cannot be guaranteed, allocate defaultly for 30 (+1)
				support[0] = 1; support[1] = e;	continue;			// new support array started, insert it's first element
			}
			
			int nodesS = support[0];
			if (nodesS + 1 > support.length) {						// check if support array needs resizing
				int nodesS2 = nodesS + (nodesS >> 1);
				int[] support2 = new int[nodesS2];
				for (int n2 = 0; n2 < nodesS; n2++) support2[n2] = support[n2];
				support = elementSupports[n] = support2;
			}
			support[0] = ++nodesS;
			support[nodesS] = e;
		}
		
		// scale down arrays to the final support counts
		for (int n = 0; n < nodes; n++) {
			int[] support = elementSupports[n];
			int nodesS = support[0] + 1;
			if (elementSupportMaxArrayL < nodesS) elementSupportMaxArrayL = nodesS;
			if (nodesS < support.length) {
				int[] support2 = new int[nodesS + 1];
				for (int n2 = 0; n2 < nodesS; n2++) support2[n2] = support[n2];
				//for (int n2 = 0; n2 < nodesS; n2++) System.out.print(support2[n2] + (n2 == nodesS - 1 ? "\n" : (n2 == 0 ? ": " : ",")));
				elementSupports[n] = support2;
			} else
				elementSupports[n] = support;
		}
	}
	
	
	// method finds the nodes fully encapsulated into the volume (those that can be moved with full degree of freedom inside the domain)
	// method is only useful for a parsed input volume with already generated tetrahedrons whose internal nodes must be found
	public int[] internalNodes() {
		
		if (element2 == null) return null;							// object-oriented elements must be initialised
		boolean[] externalFlag = new boolean[nodes];
		if (elementSupports == null) {
			elementsNodeSupport();									// we need array of elements supported by every node
			tetraNeighbours(null, true);							// and we need neighbourhood arrays of every element, WITH neighbour interfaces
		}
		
		for (int e = 0; e < elements2; e++) {						// search every element
			FEM1Element elem = getElement2(e);
			if (elem == null || !elem.hasInterfaces()) continue;
			if (elem.neighboursF < 4) {								// if element has facets on the boundary (one or more facets non-interfacing)
				int[] neighbour = elem.neighbour, nodeRef = elem.nodeRef;
				int nff1 = 0;										// OR together all interfacing facets
				for (int f = 1, fEnd = elem.neighboursF * 2; f < fEnd; f += 2) nff1 |= neighbour[f];
				// first, any node indexes associated to boundary faces are set to -1 in internal[] array
				if ((nff1 & NFF012) == 0) { externalFlag[nodeRef[0]]=true; externalFlag[nodeRef[1]]=true; externalFlag[nodeRef[2]]=true; }
				if ((nff1 & NFF023) == 0) { externalFlag[nodeRef[0]]=true; externalFlag[nodeRef[2]]=true; externalFlag[nodeRef[3]]=true; }
				if ((nff1 & NFF031) == 0) { externalFlag[nodeRef[0]]=true; externalFlag[nodeRef[1]]=true; externalFlag[nodeRef[3]]=true; }
				if ((nff1 & NFF132) == 0) { externalFlag[nodeRef[1]]=true; externalFlag[nodeRef[2]]=true; externalFlag[nodeRef[3]]=true; }
			}
		}
		int nIcount = nodes;
		if (nodeFlag == null) {													// if node flags don't exist, initialise them
			nodeFlag = new byte[node.length / 3];
			if (nodesWork > 0) nodeWorkFlag = new byte[nodesWork];
		}
		for (int nF = 0, nodes1 = node.length / 3; nF < nodes1; nF++)			// iterate over node[] array
			if (externalFlag[nF]) nIcount--; else nodeFlag[nF] |= 1;			// count final number of internal nodes
		for (int nF = node.length / 3, n2 = 0; nF < nodes; nF++, n2++)			// iterate over nodeWork[] array
			if (externalFlag[nF]) nIcount--; else nodeWorkFlag[n2] |= 1;		// count final number of internal nodes
		
		int[] internal = new int[nIcount];										// produce indexes of internal nodes into internal[] array
		for (int nI = 0, nF = 0; nF < nodes; nF++) if (!externalFlag[nF]) internal[nI++] = nF;
		return internal;
	}
	
	
	
	// method generates for every node the array of polygons that are supported by that node
	// if registerEdgeNeighbours = false, only edge-edge neighbourhoods will be registered
	void polygonsNodeSupport() {
		
		int polygonsTot = polygonOffset[polygons], nCount = polygonOffset[1] - polygonOffset[0];	// polygonsTot = total node count of all polygons
		polygonSupports = new int[nodes][];
		for (int p = 0, p3 = 0, p3cnt = 0; p3 < polygonsTot; p3++) {

			if (p3cnt >= nCount) { p++; p3cnt = 0; nCount = -polygonOffset[p] + polygonOffset[p + 1]; }
			int n = polygon[p3];
			int[] support = polygonSupports[n];
			
			if (support == null) {										// optimally, every node is shared by 6 faces or by 4 quartics (+ 1 array length integer)
				support = polygonSupports[n] = new int[31];				// but as an optimal mesh cannot be guaranteed, allocate defaultly for 30 (+1)
				support[0] = 1;
				support[1] = p;											// new support array started, insert it's first element
			} else {			
				int nodesS = support[0];
				if (nodesS + 1 > support.length) {						// check if support array needs resizing
					int nodesS2 = nodesS + (nodesS >> 1);				// increase array 1.5 times
					int[] support2 = new int[nodesS2];
					for (int n2 = 0; n2 < nodesS; n2++) support2[n2] = support[n2];
					support = polygonSupports[n] = support2;
				}
				
				// for facetNeighbours() method to work properly, we must make sure the indexes are sorted
				if (support[nodesS] <= p) support[++nodesS] = p;		// if this index is higher than last one, just append
				//else if (support[nodesS] == p) throw new InvalidParameterException("polygonsNodeSupport(): node supports polygon " + p + "twice.");
				else {
					for (int s = 1; s <= nodesS; s++) {					// otherwise insert polygon index in sorted order
						if (p < support[s]) {
							for (int s1 = nodesS + 1; s1 > s; s1++) support[s1] = support[s];
							support[s] = p; break;
						}
					}
				}
				support[0]++;
			}
			p3cnt++;
		}
		
		// scale down arrays to the final support counts
		for (int n = 0; n < nodes; n++) {
			int[] support = polygonSupports[n];
			int nodesS = support[0] + 1;
			if (polygonSupportMaxArrayL < nodesS) polygonSupportMaxArrayL = nodesS;
			if (nodesS < support.length) {
				int[] support2 = new int[nodesS + 1];	// note: an extra slot on end added for optimisation of facetNeighbours() method
				for (int n2 = 0; n2 < nodesS; n2++) support2[n2] = support[n2];
				//for (int n2 = 0; n2 < nodesS; n2++) System.out.print(support2[n2] + (n2 == nodesS - 1 ? "\n" : (n2 == 0 ? ": " : ",")));
				polygonSupports[n] = support2;
			}
		}
	}

	
	
	final static int FL_UNIQUE = 0x80000000, FL_UN_CLR = 0x7FFFFFFF;
	
	// method is optimised for finding neighbourhoods of FACETS, method accepts only tri-faces as polygons
	// the edge & corner neighbours are put in sorted order: edges, border edges, corners, border corners (leading 4 ints are counts for those)
	// if skipNodeNeighbours = true, the corner neighbours will be ignored
	// if onlyTrifaces = true, method will pass over polygons that are not facets, otherwise everything MUST be facets or method fails
	
	// TODO: take care of the unallowed  case of non-closed volumes, a facet can border empty space, that case can be registered
	// in arrays called n123ee[] && n123ce[] for example
	void facetNeighbours(boolean onlyTrifaces, boolean skipCornerNeighbours) {
		
		polygonNeighbours = new int[polygons][];							// write to global face neighbour array
		borderPolygons = new int[polygons + patches];						// store down facets found to be border polygons
		
		// we'll first compare the nodes pairwise: min(n1, n2), then pick the lowest index: min(min(n1, n2), n3)
		int[] n12 = new int[polygonSupportMaxArrayL * 4 + 1];				// worst-case allocation: all uniques in both arrays -> x2 integers
		int[] n123e = new int[3];											// n123e[] carries only edge neighbours
		int[] n123eb = new int[3];											// n123ec[] carries only patch border edge neighbours
		int[] n123c = new int[polygonSupportMaxArrayL * 2];					// n123c[] carries only corner neighbours
		int[] n123cb = new int[polygonSupportMaxArrayL * 2];				// n123cc[] carries only patch border corner neighbours
		int pch = 0, pchStart = patch[pch], pchEnd = patch[pch];			// polygon indexes of current patch maintained by pchStart & pchEnd
			
		for (int p = 0, p3 = 0, pB = 0, pBn = 0; p < polygons; p++) {		// pB indexes current border face, pBn indexes current patch count of border faces
						
			if (polygonOffset[p + 1] - polygonOffset[p] > 3) {
				if (onlyTrifaces) continue;									// if flagged to skip non-facets
				else throw new InvalidParameterException("FEM1.faceNeighbours(): Non-facet encountered.");
			}				
			
			if (p >= pchEnd) {												// if we passed on into next patch
				pchStart = pchEnd;											// update patch delimitations
				pchEnd = patch[++pch];
				pBn = pB++;													// set aside a cell for border polygon count
			}

			// get the facet arrays that form support of all three nodes of facet
			// note: topologically, edge 01 derives from (a1[], a2[]), edge 12 from (a2[], a3[]), edge 20/02 from (a3[], a1[])
			// to maintain facet's edge identity, any unique stored in n12[] will be followed by a facet's node code: 0, 1, 2
			int[] a1 = polygonSupports[polygon[p3++]], a2 = polygonSupports[polygon[p3++]], a3 = polygonSupports[polygon[p3++]];

			int i12 = 0, i1 = 1, i2 = 1, c1 = a1[0], c2 = a2[0];
			int po1 = a1[i1], po2 = a2[i2];
			// TODO: small optimisation can be made, which stops checking self-referencing in consecutive iteration blocks
			while (i1 <= c1 && i2 <= c2) {
				if (po1 == p)		{ po1 = a1[++i1]; continue; }			// skip reference to itself
				if (po2 == p)		{ po2 = a2[++i2]; continue; }
				if (po1 < po2)		{ n12[i12++] = po1 | FL_UNIQUE; n12[i12++] = 0; po1 = a1[++i1]; }
				else if (po1 > po2)	{ n12[i12++] = po2 | FL_UNIQUE; n12[i12++] = 1; po2 = a2[++i2]; }
				else				{ n12[i12++] = po1; po1 = a1[++i1]; po2 = a2[++i2]; }
			}
			// the remains are all uniques, but skip self-refs
			while (i1 <= c1) { if (a1[i1] == p) i1++; else { n12[i12++] = a1[i1++] | FL_UNIQUE; n12[i12++] = 0; }}
			while (i2 <= c2) { if (a2[i2] == p) i2++; else { n12[i12++] = a2[i2++] | FL_UNIQUE; n12[i12++] = 1; }}

			int i123e = 0, i123eb = 0, i123c = 0, i123cb = 0;				// initiate the four counters
			int c12 = i12, c3 = a3[0], i3 = 1;
			i12 = 0;
			int po12 = n12[i12], po3 = a3[i3];
			n123e[0]=-1; n123e[1]=-1; n123e[2]=-1; n123eb[0]=-1; n123eb[1]=-1; n123eb[2]=-1;			// -1 identifies a nonexistent edge of either type
			
			if (skipCornerNeighbours) {																	// if only accepting non-uniques (= edges)
				while (i12 < c12 && i3 <= c3) {
					if (po3 == p)		{ po3 = a3[++i3]; continue; }									// skip references to current facet itself
					if (po12 < 0) {																		// if 1st comparison produced a unique
						int po12c = po12 & FL_UN_CLR;
						if (po12c < po3)	{	i12++; po12 = n12[++i12]; }								// if it's lower -> skip it
						else if (po12c > po3)	po3 = a3[++i3];											// if higher -> a3[] has an unique, skip second
						else { 																			// if they equal -> non-unique -> store whichever
							if (po12c < pchStart || po12c >= pchEnd) {									// if belonging to another patch, put in n123eb[]
									{ n123eb[n12[++i12] + 1] = po12c; i123eb++; }						// (following edge code decides array position)
							} else	{ n123e[n12[++i12] + 1] = po12c; i123e++; }							// else store as edge in n123e[]
							po12 = n12[++i12]; po3 = a3[++i3];
						}
						continue;
					}
					if (po12 < po3) {																	// 1st comparison non-unique lower -> store it
						if (po12 < pchStart || po12 >= pchEnd) {	 									// if belonging to another patch, put in n123eb[]
								{ n123eb[0] = po12; i123eb++; }											// (this being edge01, store at position 0)
						} else	{ n123e[0] = po12; i123e++; }
						po12 = n12[++i12]; }
					else if (po12 > po3)	po3 = a3[++i3];												// a3[] lower -> it's unique -> skip it
					// following case is topologically illegal, since two facets can only have two nodes equal, unless they mirror each other
					else throw new InvalidParameterException("FEM1.facetNeighbours(): 2 facets neighbours on all three nodes.");
//					else {																				// if they equal -> non-unique -> store whichever
//						if (po12 < pchStart || po12 >= pchEnd)	n123ec[i123ec++] = po12;				// if belonging to another patch, put in n123ec[]
//						else									n123e[i123e++] = po12;
//					}
				}
				// the remains of a3[] are by definition all uniques -> skipped
				while (i12 < c12)
					if (n12[i12] >= 0) {																// skip uniques from remains of second comparison
						if (n12[i12] < pchStart || n12[i12] >= pchEnd)
								{ n123eb[0] = n12[i12++]; i123eb++; }									// if belonging to another patch, put in n123eb[]
						else	{ n123e[0] = n12[i12++]; i123e++; }
					} else i12 += 2;																	// jump over node index & node code
				
			// the case of accepting both edge & corner neighbours
			} else {
				while (i12 < c12 && i3 <= c3) {
					if (po3 == p)		{ po3 = a3[++i3]; continue; }									// skip references to current facet itself
					if (po12 < 0) {																		// if 1st comparison produced a unique
						int po12c = po12 & FL_UN_CLR;
						if (po12c < po3)		{														// if it's lower -> unique -> put as corners
							if (po12c < pchStart || po12 >= pchEnd) n123cb[i123cb++] = po12c;			// if belonging to another patch, put in n123eb[]
							else									n123c[i123c++] = po12c;
							i12++; po12 = n12[++i12];
						} else if (po12c > po3){														// if higher, a3[] has an unique -> put as corners
							if (po3 < pchStart || po3 >= pchEnd)	n123cb[i123cb++] = po3;				// if belonging to another patch, put in n123cb[]
							else									n123c[i123c++] = po3;
							po3 = a3[++i3];
						} else {																		// if they equal -> non-unique -> put as edges
							if (po12c < pchStart || po12c >= pchEnd)
									{ n123eb[n12[++i12] + 1] = po12c; i123eb++; }						// if belonging to another patch, put in n123eb[]
							else	{ n123e[n12[++i12] + 1] = po12c; i123e++; }
							po12 = n12[++i12]; po3 = a3[++i3];
						}
						continue;
					}
					if (po12 < po3)			{
						if (po12 < pchStart || po12 >= pchEnd)
								{ n123eb[0] = po12; i123eb++; }											// if belonging to another patch, put in n123eb[]
						else	{ n123e[0] = po12; i123e++; }											// (this is edge 01 of facet)
						po12 = n12[++i12];
					} else if (po12 > po3)	{
						if (po3 < pchStart || po3 >= pchEnd)	n123cb[i123cb++] = po3;					// if belonging to another patch, put in n123cb[]
						else									n123c[i123c++] = po3;
						po3 = a3[++i3];
					// following case is topologically illegal, since two facets can only have two nodes equal, unless they mirror each other
					} else throw new InvalidParameterException("FEM1.facetNeighbours(): 2 facets neighbour on all three nodes.");
//					} else					{
//						if (po12 < pchStart || po12 >= pchEnd)	n123ec[i123eb++] = po12;				// if belonging to another patch, put in n123eb[]
//						else									n123e[i123e++] = po12;
//						po12 = n12[++i12]; po3 = a3[++i3];
//					}
				}
				while (i12 < c12)
					if (n12[i12] >= 0) {																// if a non-unique
						if (n12[i12] < pchStart || n12[i12] >= pchEnd)
								{ n123eb[0] = n12[i12++]; i123eb++; }									// if belonging to another patch, put in n123cc[]
						else	{ n123e[0] = n12[i12++]; i123e++; }
					} else {
						int po12c = n12[i12] & FL_UN_CLR;
						if (po12c < pchStart || po12c >= pchEnd) n123cb[i123cb++] = po12c;				// if belonging to another patch, put in n123cc[]
						else n123c[i123c++] = po12c;
						i12 += 2;
					}
				while (i3 <= c3) {																		// a3[] remains are all uniques
					if (a3[i3] == p) { i3++; continue;	}												// but self-references are skipped
					if (a3[i3] < pchStart || a3[i3] >= pchEnd)	n123cb[i123cb++] = a3[i3++];			// if belonging to another patch, put in n123cc[]
					else										n123c[i123c++] = a3[i3++];
				}	
			}
			
			// if at least one patch-bordering edge found, store in border polygon indexer
			if (i123eb > 0) { borderPolygons[pB++] = p; borderPolygons[pBn]++; }

			// allocate to correct array size and move all data sequentially
			int nTotal = 5 + 3 + 3 + (skipCornerNeighbours ? 0 : i123c + i123cb);
			int[] neighbours = polygonNeighbours[p] = new int[nTotal];
			neighbours[0] = nTotal;
			int p1 = 5;
			neighbours[p1++] = n123e[0]; neighbours[p1++] = n123e[1]; neighbours[p1++] = n123e[2];
			neighbours[p1++] = n123eb[0]; neighbours[p1++] = n123eb[1]; neighbours[p1++] = n123eb[2];
			if (!skipCornerNeighbours) {
				for (int p2 = 0; p2 < i123c; p1++, p2++) neighbours[p1] = n123c[p2];
				for (int p2 = 0; p2 < i123cb; p1++, p2++) neighbours[p1] = n123cb[p2];
			}
			neighbours[1] = i123e; neighbours[2] = i123eb; neighbours[3] = i123c; neighbours[4] = i123cb;
			// DEBUG: array printout
//			for (int pN = 0; pN < nTotal; pN++)
//				if (neighbours[pN] >= 0) System.out.print(neighbours[pN] + (pN == nTotal - 1 ? "\n" : (pN == 4 ? ": " : ",")));
//				else if (pN == nTotal - 1) System.out.print("\n");
		}
	}
	
	
	// method generates nonredundant mapping of every facet edge to it's length, where facet edges can share a length
	// the mapping is from a circular definition of an edge as point(n) to point(n+1), it index-mirrors the polygon nodes array polygon[] exactly,
	// and points to an array of edge lengths polygonEdge[] and an array of node index duplets polygonEdgeN[] (for back-reference to underlying edge data)
	void facetEdgeLengths() {
		int polyNodes = polygonOffset[polygons];
		polygonEdgeIndex = new int[polyNodes];									// polygonEdgeIndex[] is the mapper to the edge lengths
		for (int i = 0; i < polyNodes; i++) polygonEdgeIndex[i] = -1;			// set all entries to unmapped, -1
		double[] polygonEdgeT = new double[polyNodes];							// worst-case allocation
		int[] polygonEdgeNT = new int[polyNodes * 2];							// worst-case allocation
		int eC = 0, eNP = 0;													// the edge counter & the node pair counter
		
		// every consecutive facet is checked if it's edges are already remapped by a previously checked neighbour
		// a new edge length is generated case an edge hasn't been mapped yet
		for (int f = 0; f < polygons; f++) {
			int n = polygonOffset[f];
			if (polygonEdgeIndex[n] == -1) {									// check if current facet edge n0-n1 hasn't been mapped yet
				polygonEdgeT[eC] = distance(polygonEdgeNT[eNP++] = polygon[n++], polygon[polygonEdgeNT[eNP++] = n--], true);
				polygonEdgeIndex[n] = eC++; } n++;								// map to unique position in polygonEdge[]		
			if (polygonEdgeIndex[n] == -1) {									// check if current facet edge n1-n2 hasn't been mapped yet
				polygonEdgeT[eC] = distance(polygonEdgeNT[eNP++] = polygon[n++], polygonEdgeNT[eNP++] = polygon[n--], true);
				polygonEdgeIndex[n] = eC++; } n++;								// map to unique position in polygonEdge[]				
			if (polygonEdgeIndex[n] == -1) {									// check if current facet edge n2-n0 hasn't been mapped yet
				polygonEdgeT[eC] = distance(polygonEdgeNT[eNP++] = polygon[n], polygonEdgeNT[eNP++] = polygon[n - 2], true);
				polygonEdgeIndex[n] = eC++; }									// map to unique position in polygonEdge[]	

			// check neighbours and remap them to the same indexes as current facet in polygonEdgeT[]
			int[] facetNgb = polygonNeighbours[f];
			int fN01 = facetNgb[5] == -1 ? facetNgb[8] : facetNgb[5], fN12 = facetNgb[6] == -1 ? facetNgb[9] : facetNgb[6];
			int fN20 = facetNgb[7] == -1 ? facetNgb[10] : facetNgb[7];	// collect indexes of neighbour facets, from both normal & border edges
			n = polygonOffset[f];
			
			if (fN01 > f) {		// neighbours of lower index than current facet have already been fully mapped
				int[] ngb01 = polygonNeighbours[fN01];
				if		(ngb01[5] == f || ngb01[8] == f)	polygonEdgeIndex[polygonOffset[fN01]] = polygonEdgeIndex[n];
				else if (ngb01[6] == f || ngb01[9] == f)	polygonEdgeIndex[polygonOffset[fN01] + 1] = polygonEdgeIndex[n];
				else if (ngb01[7] == f || ngb01[10] == f)	polygonEdgeIndex[polygonOffset[fN01] + 2] = polygonEdgeIndex[n];
			}	n++;
			if (fN12 > f) {
				int[] ngb12 = polygonNeighbours[fN12];
				if		(ngb12[5] == f || ngb12[8] == f)	polygonEdgeIndex[polygonOffset[fN12]] = polygonEdgeIndex[n];
				else if (ngb12[6] == f || ngb12[9] == f)	polygonEdgeIndex[polygonOffset[fN12] + 1] = polygonEdgeIndex[n];
				else if (ngb12[7] == f || ngb12[10] == f)	polygonEdgeIndex[polygonOffset[fN12] + 2] = polygonEdgeIndex[n];
			}	n++;
			if (fN20 > f) {
				int[] ngb20 = polygonNeighbours[fN20];
				if		(ngb20[5] == f || ngb20[8] == f)	polygonEdgeIndex[polygonOffset[fN20]] = polygonEdgeIndex[n];
				else if (ngb20[6] == f || ngb20[9] == f)	polygonEdgeIndex[polygonOffset[fN20] + 1] = polygonEdgeIndex[n];
				else if (ngb20[7] == f || ngb20[10] == f)	polygonEdgeIndex[polygonOffset[fN20] + 2] = polygonEdgeIndex[n];
			}
		}
		
		// time to rescale polygonEdgeT[] & polygonEdgeNT[] to their final size
		polygonEdge = new double[polygonEdges = eC];
		polygonEdgeN = new int[polygonEdges * 2];
		for (int e = 0, e2 = 0; e < eC; e++) {
			polygonEdge[e] = polygonEdgeT[e];
			polygonEdgeN[e2] = polygonEdgeNT[e2++]; polygonEdgeN[e2] = polygonEdgeNT[e2++];
		}
	}
	
	
	
	// method locates the crease edges of the patches/smoothing groups, laying them out in cyclic order
	// method is specifically operating on facets, NOT generic polygons
	void creaseBoundariesOfFacets() {
		
		for (int pch = 0, pB = 0; pch < patches; pch++) {							// iterate over patches
			// worst case allocation of temp. array for border cycles
			int pchSize = (pch + 1 < patches ? polygons : patch[pch + 1]) - patch[pch], pchSize4 = pchSize * 4;
			// worst case is if every facet in patch is bordering other patches with every edge (repeated cycles of 3 nodes -> topologically discontinuous)
			// since every cycle is preceded by a cycle length count, allocate one extra slot -> 1 + 3 = 4
			int[] borderCycles = new int[pchSize * 4];
			int cyCnt = 0, cyOffs = 1;												// note: at any given time, nodes are sought between cyCnt+1 & cyOffs
			
			for (int pBEnd = pB + borderPolygons[pB]; pB < pBEnd; pB++) {			// iterate over this patch's polygons (border poly.count lying in same array)
				int p1 = borderPolygons[pB], n3 = p1 * 3;							// get current polygon's index (p1) & it's node offset (n3)
				int[] neighbours = polygonNeighbours[p1];							// get current polygon's neighbour list
				int n10 = polygon[n3++], n11 = polygon[n3++], n12 = polygon[n3++];	// get current polygon's nodes
				
				if (borderCycles[cyCnt] == 0) {										// if cycle construction has just begun, store current polygon's nodes
					borderCycles[cyOffs++] = n10; borderCycles[cyOffs++] = n11; borderCycles[cyOffs++] = n12;
					borderCycles[cyCnt] = 3;
					if (neighbours[1] == polygonOffset[p1] - polygonOffset[p1 + 1]) {	// cover the case of all polygon edges bordering other patches
						cyCnt = cyOffs++;
						continue;
					}
				}
				
				if (neighbours[4] >= 0) {													// if edge01-neighbour exists
					int[] neighbours2 = polygonNeighbours[neighbours[4]];
					int n3n = neighbours[4] * 3;											// get neighbour's offset into node array
					int n20 = polygon[n3n++], n21 = polygon[n3n++], n22 = polygon[n3n];		// get indexes of nodes of neighbour polygon
					if (neighbours2[1] > 0) {												// if neighbour polygon also has a boundary neighbour edge
						if (neighbours2[7] >= 0) {											// if that is edge01
							if (n10 == n20)	{												// if current poly's node 0 is chained to neighbour's node 0
								borderCycles[cyOffs - 2] = n21;								// chain on neighbour's node 1
							} else if (n10 == n22) {										// if current poly's node 0 is chained to neighbour's node 2
								borderCycles[cyOffs - 2] = n22;								// chain on neighbour's node 2
							}
						}
					}
						
				}
				// iterate over current facet's neighbour edges that are in same patch
				for (int nbr = 4/*, nbrEnd = neighbours[0] + 4*/; nbr < /*nbrEnd*/7; nbr++) {
					
					if (neighbours[nbr] < 0) continue;
					int[] neighbours2 = polygonNeighbours[neighbours[nbr]];					// get the neighbour polygon's neighbour list
					if (neighbours2[1] > 0) {												// if neighbour polygon also has a boundary neighbour edge					
						int n3n = neighbours[nbr] * 3;										// get neighbour's offset into node array
						int n20 = polygon[n3n++], n22 = polygon[n3n++], n23 = polygon[n3n];	// get indexes of nodes of neighbour polygon
						
						for (int cy = cyCnt + 1; cy < cyOffs; cy++) {						// find neighbour's node cycle chain to node of current polygon's
							if (borderCycles[cy] == n20) {
								borderCycles[++cy] = n22;
							}
						}
					}
					
				}
				
			}
		}
	}
	
	
	
	final static int NFF012=1, NFF023=2, NFF031=4, NFF132=8, NEE01=16, NEE02=32, NEE03=64, NEE12=128, NEE13=256, NEE23=512;
	static final int FULL_NEIGHBOURHOOD = NFF012+NFF023+NFF031+NFF132+NEE01+NEE02+NEE03+NEE12+NEE13+NEE23;
	static final int NFF_BITS = NFF012+NFF023+NFF031+NFF132;
	static final int FULL_0_NEIGHBOURHOOD = NFF012+NFF023+NFF031+NEE01+NEE02+NEE03, FULL_1_NEIGHBOURHOOD = NFF012+NFF031+NFF132+NEE01+NEE12+NEE13;
	static final int FULL_2_NEIGHBOURHOOD = NFF012+NFF023+NFF132+NEE02+NEE12+NEE23, FULL_3_NEIGHBOURHOOD = NFF023+NFF031+NFF132+NEE03+NEE13+NEE23;

	// method is optimised for finding neighbourhood of tetrahedrons
	// if an element list is supplied, the method assumes that new elements were inserted, and will call found neighbours to alter their neighbourhoods
	// findInterfaces = true means the array will consist of element index & interface bitcode pairs, method finds the edge or facet interfaces
	// towards the neighbours and employs a flagging array for found neighbour interfaces with following bits (0,1,2,3 being tetrahedral vertexes):
	//		facet-facet: 012, 023, 031, 132
	//		edge-edge: 01, 02, 03, 12, 13, 23			(all in all 10 bits for a tetrahedron, 10 + 10(neighbour) bits packed in an integer)
	void tetraNeighbours(int[] elemList, boolean findInterfaces) {

		if (elementSupports == null) elementsNodeSupport();				// we need array of elements supported by each node (unless previously calculated)
		
		// we'll first compare the nodes pairwise: min(n1, n2) & min (n3, n4), then pick the lowest index: min(min(n1, n2), min (n3, n4))
		// pairwise comparisons go into these temporary arrays, then again compared into neighboursTmp[]
		int[] n12 = new int[elementSupportMaxArrayL * 4 + 1];			// comparison of nodes 0 and 1
		int[] n34 = new int[elementSupportMaxArrayL * 4 + 1];			// comparison of nodes 2 and 3
		int[] n1234F = new int[elementSupportMaxArrayL * 8];			// comparisons 0&1 and 2&3, triple equalities go here as faces
		int[] n1234e = new int[elementSupportMaxArrayL * 8];			// comparisons 0&1 and 2&3, single equalities go here as edges

		int elemCount = elemList == null ? elements2 : elemList.length;
		for (int e = 0; e < elemCount; e++) {
			
			int[] a1, a2, a3, a4;
			FEM1Element elem1 = elemList == null ? getElement2(e) : getElement2(elemList[e]);
			if (elem1 == null) continue;								// skip recently deleted elements
			int[] nodeRef1 = elem1.nodeRef;
			a1 = elementSupports[nodeRef1[0]]; a2 = elementSupports[nodeRef1[1]];
			a3 = elementSupports[nodeRef1[2]]; a4 = elementSupports[nodeRef1[3]];

			// we need to compare four indexes from four parallel sorted arrays, picking the lowest one every time
			// duplicates are stored normally, unique indexes stored intermediately with 31st bit set
			// a final unique index is a tetrahedron touching another tetrahedron with just one node, we're only interested in edges & faces
			int i12 = 0, i1 = 1, i2 = 1, c1 = a1[0], c2 = a2[0];
			int el1 = a1[i1], el2 = a2[i2];
			while (i1 <= c1 && i2 <= c2) {
				if (el1 == e)		{ el1 = a1[++i1]; continue; }											// skip reference to current element itself
				if (el2 == e)		{ el2 = a2[++i2]; continue; }											// (it will always be it's own neighbour)
				if (el1 < el2)		{ n12[i12++] = el1 | FL_UNIQUE; n12[i12++] = 0; el1 = a1[++i1]; }		// for uniques, additionally store the node index
				else if (el1 > el2)	{ n12[i12++] = el2 | FL_UNIQUE; n12[i12++] = 1; el2 = a2[++i2]; }
				else				{ n12[i12++] = el1;	el1 = a1[++i1]; el2 = a2[++i2]; }					// a neighbour on edge01 found										
			}
			// add remains as uniques, but skip self-references
			while (i1 <= c1) { if (a1[i1] == e) i1++; else { n12[i12++] = a1[i1++] | FL_UNIQUE; n12[i12++] = 0; }}
			while (i2 <= c2) { if (a2[i2] == e) i2++; else { n12[i12++] = a2[i2++] | FL_UNIQUE; n12[i12++] = 1; }}

			int i34 = 0, i3 = 1, i4 = 1, c3 = a3[0], c4 = a4[0];
			int el3 = a3[i3], el4 = a4[i4];
			while (i3 <= c3 && i4 <= c4) {
				if (el3 == e)		{ el3 = a3[++i3]; continue; }											// skip reference to current element itself
				if (el4 == e)		{ el4 = a4[++i4]; continue; }											// (it will always be it's own neighbour)
				if (el3 < el4)		{ n34[i34++] = el3 | FL_UNIQUE; n34[i34++] = 2; el3 = a3[++i3]; }		// for uniques, additionally store the node index
				else if (el3 > el4)	{ n34[i34++] = el4 | FL_UNIQUE; n34[i34++] = 3; el4 = a4[++i4]; }
				else				{ n34[i34++] = el3; el3 = a3[++i3]; el4 = a4[++i4]; }					// a neighbour on edge23 found
			}
			// add remains as uniques, but skip self-references
			while (i3 <= c3) { if (a3[i3] == e) i3++; else { n34[i34++] = a3[i3++] | FL_UNIQUE; n34[i34++] = 2; }}
			while (i4 <= c4) { if (a4[i4] == e) i4++; else { n34[i34++] = a4[i4++] | FL_UNIQUE; n34[i34++] = 3; }}

			int i1234F = 0, i1234e = 0, c12 = i12, c34 = i34;
			i12 = 0; i34 = 0;
			int el12 = n12[i12], el34 = n34[i34];
			while (i12 < c12 && i34 < c34) {
				
				if (el12 < 0) {																				// if first comparison produced a unique
					int el12c = el12 & FL_UN_CLR;
					if (el34 < 0) {																			// case of first & second being unique
						int el34c = el34 & FL_UN_CLR;
						if (el12c < el34c)		{ i12++; el12=n12[++i12]; }									// first is unique & lower -> skip first
						else if (el12c > el34c)	{ i34++; el34=n34[++i34]; }									// second is unique & lower -> skip second
						else {
							n1234e[i1234e++] = el12c;														// first = second -> not unique, store as edge
							if (findInterfaces) n1234e[i1234e++] = elem1.tetraEdgeInterface(n12[++i12], n34[++i34], getElement2(el12c).nodeRef);
							else { i12++; i34++; }
						 	el12=n12[++i12]; el34=n34[++i34]; }
					} else {																				// case of first unique, second not unique
						if (el12c < el34)		{ i12++; el12=n12[++i12]; }									// first is unique & lower -> skip it
						else if (el12c > el34) {
							n1234e[i1234e++] = el34;														// first unique & higher -> store 2nd as edge23
							if (findInterfaces) n1234e[i1234e++] = elem1.tetraEdgeInterface(2, 3, getElement2(el34).nodeRef);
							el34=n34[++i34]; }
						else {																				// first u. = second = triple eq. -> store as facet
							n1234F[i1234F++] = el34;
							if (findInterfaces) n1234F[i1234F++] = elem1.tetraNeighbourInterface(getElement2(el34).nodeRef);
							i12++; el12=n12[++i12]; el34=n34[++i34]; }
					}
					continue;
				}
				if (el34 < 0) {																				// if second comparison produced a unique
					int el34c = el34 & FL_UN_CLR;
					if (el12 < 0) {																			// case of first & second being unique
						int el12c = el12 & FL_UN_CLR;
						if (el12c < el34c)		{ i12++; el12=n12[++i12]; }									// first is unique & lower -> skip first
						else if (el12c > el34c)	{ i34++; el34=n34[++i34]; }									// second is unique & lower -> skip second
						else {
							n1234e[i1234e++] = el12c;														// first = second -> not unique, store as edge
							if (findInterfaces) n1234e[i1234e++] = elem1.tetraEdgeInterface(n12[++i12], n34[++i34], getElement2(el12c).nodeRef);
							else { i12++; i34++; }
							el12=n12[++i12]; el34=n34[++i34]; }
					} else {																				// case of second unique, first not unique
						if (el12 < el34c)		{
							n1234e[i1234e++] = el12;														// second is unique & higher -> store 1st as edge01
							if (findInterfaces) n1234e[i1234e++] = elem1.tetraEdgeInterface(0, 1, getElement2(el12).nodeRef);
							el12=n12[++i12]; }
						else if (el12 > el34c)	{ i34++; el34=n34[++i34]; }									// second is unique & lower -> skip it
						else {
							n1234F[i1234F++] = el12;														// first = second u. = triple eq. -> store as face
							if (findInterfaces) n1234F[i1234F++] = elem1.tetraNeighbourInterface(getElement2(el12).nodeRef);
							el12=n12[++i12]; i34++;  el34=n34[++i34]; }
					}
					continue;
				}
				if (el12 < el34)		{																	// nonunique 1st less than nonunique 2nd -> store as edge01
					n1234e[i1234e++] = el12;
					if (findInterfaces) n1234e[i1234e++] = elem1.tetraEdgeInterface(0, 1, getElement2(el12).nodeRef);
					el12=n12[++i12]; }
				else if (el12 > el34)	{																	// nonunique 2nd less than nonunique 1st -> store as edge23
					n1234e[i1234e++] = el34;
					if (findInterfaces) n1234e[i1234e++] = elem1.tetraEdgeInterface(2, 3, getElement2(el34).nodeRef);
					el34=n34[++i34]; }
				// essentially a topologically invalid state, since it would mean that an element neighbours another element on all four nodes
				//else					{ n1234F[i1234F++] = el12; el12=n12[++i12]; el34=n34[++i34]; }
				else throw new InvalidParameterException("FEM1.tetraNeighbours(): 2 tetrahedral elements neighbours on all four nodes.");
			}
			// the remains are all either edges or uniques/corner neighbours
			while (i12 < c12)
				if (n12[i12] >= 0) {
					n1234e[i1234e++] = n12[i12];
					if (findInterfaces) n1234e[i1234e++] = elem1.tetraEdgeInterface(0, 1, getElement2(n12[i12++]).nodeRef);
					else i12++; }
				else i12 += 2;
			while (i34 < c34)
				if (n34[i34] >= 0) {
					n1234e[i1234e++] = n34[i34];
					if (findInterfaces) n1234e[i1234e++] = elem1.tetraEdgeInterface(2, 3, getElement2(n34[i34++]).nodeRef);
					else i34++; }
				else i34 += 2;
			
			// resize to correct array size and store the rehashed indexes as neighbours
			int nTotal = i1234F + i1234e;
			if (elementNeighboursMaxArray < nTotal) elementNeighboursMaxArray = nTotal;
			int[] neighbours = elem1.neighbour = new int[nTotal];
			if (findInterfaces) {
					elem1.neighbours = nTotal / 2; elem1.neighboursF = i1234F / 2; }
			else {	elem1.neighbours = nTotal; elem1.neighboursF = i1234F; }
						
			int eN = 0;
			if (elemList == null) {						// this block does regular initialisation where every element's neighbourhood is made from scratch
				for (int eF = 0; eF < i1234F; eF++) neighbours[eN++] = n1234F[eF];
				for (int eE = 0; eE < i1234e; eE++) neighbours[eN++] = n1234e[eE];
			} else {									// this block assumes this being a new element, and updates it's found neighbourhood
				for (int eF = 0; eF < i1234F; eF++) {
					if ((eF & 1) == 0) {
						neighbours[eN++] = n1234F[eF];
						FEM1Element elem2 = getElement2(n1234F[eF]);
						if (elem2 == null) continue;
						if (elem2.hasInterfaces()) {
								int flags = n1234F[eF + 1];
								elem2.insertNeighbour(elemList[e], (flags >> 16) | (flags << 16));
						} else	elem2.insertNeighbour(elemList[e], 0);								
					}
				}
				for (int eE = 0; eE < i1234e; eE++) {
					neighbours[eN++] = n1234e[eE];
					if ((eE & 1) == 0) {
						FEM1Element elem2 = getElement2(n1234e[eE]);
						if (elem2 == null) continue;
						if (elem2.hasInterfaces()) {
								int flags = n1234F[eE + 1];
								elem2.insertNeighbour(elemList[e], (flags >> 16) | (flags << 16));
						} else	elem2.insertNeighbour(elemList[e], 0);								
					}
				}
			}
			if (findInterfaces) {
				elem1.flagHasInterfaces();
				elem1.sortInterfaces();
			}

			//for (eN = 0; eN < nTotal; eN+=2) System.out.print(neighbours[eN] + (eN == nTotal - 2 ? "\n" : (eN == 0 ? ": " : ",")));
		}
	}
	
	
	
	public FEM1Element[] deleteElements(int[] deletionRef) {
		
		FEM1Element[] deletion = new FEM1Element[deletionRef.length];
		for(int e = 0; e < deletionRef.length; e++) {							// move deleted element pointers to temp.array while nulling them
			int dRef = deletionRef[e];
			deletion[e] = element2[dRef]; element2[dRef] = null;
		}
		for(int e = 0; e < deletionRef.length; e++) {
			FEM1Element elem1 = deletion[e];
			int[] neighbour = elem1.neighbour;
			boolean hasInterfaces = elem1.hasInterfaces();
			for (int n = 0, ngbrs = elem1.neighbours; n < ngbrs; n++) {			// remove referrals to element from it's neighbourhood
				FEM1Element elem2 = getElement2(neighbour[hasInterfaces ? n * 2 : n]);
				if (elem2 != null) elem2.deleteNeighbour(deletionRef[e]);
			}
		}
		return deletion;
	}
	
	
	
	// method arranges elements in worst-to-best quality sorted element index array, whose length will be wMax (if wMax = 0, ALL elements are sorted)
	// goodCriterion & badCriterion are top & bottom limits for an element, all above goodCriterion are good, all below badCriterion are bad
	public int[] tetraWorstSort(int wMax, double goodCriterion, double badCriterion, boolean localise) {

		if (element2 == null) return null;
		if (badCriterion >= goodCriterion) throw new InvalidParameterException("FEM1.tetraWorstSort(): Bad threshold must be less than good threshold"); 
		int wN = 0, cWorst = 0;
		if (wMax == 0) wMax = elements;
		int[] worst = new int[wMax];
		double[] quality = new double[wMax];
		
		for (int e = 0; e < elements2; e++) {
			
			FEM1Element elem = getElement2(e);
			if (elem == null) continue;								// skip recently deleted elements
			double q = elem.tetraQuality();							// get quality of element, propagate data to neighbours
			elem.propagateInterfaces(this, FEM1Element.PROPAGATE_AREAS | FEM1Element.PROPAGATE_EDGES, false);
			if (q >= goodCriterion || (wN == wMax && q > quality[wN - 1])) continue;
			
			boolean sortedInsert = false;
			for (int w = 0; w < wN; w++) {
				if (q < quality[w]) {
					// while moving elements ahead, less worse elements will be pushed out of queue on far end
					for (int w1 = wN < wMax ? wN : wN - 1; w1 > w; w1--) {
						quality[w1] = quality[w1 - 1]; 
						worst[w1] = worst[w1 - 1]; }
					if (wN < wMax) wN++;							// the queue grows being bounded by wMax
					quality[w] = q;									// insert quality into queue
					worst[w] = e;									// insert element into queue
					sortedInsert = true; break;
				}
			}
			if (!sortedInsert && wN < wMax) {
				quality[wN] = q;
				worst[wN++] = e;
			}
			if (q <= badCriterion && ++cWorst >= wMax)				// if wMax guaranteedly bad tetrahedrons were found, we're done
				return worst;
		}
		return worst;
	}
	
	
	
	// given three existing nodes, spawns a tetrahedron with optimal 4th node position, nodes could be from spawning tetrahedron or boundary facet
	// method demands locally-proper-ordered input nodes: anticlockwise from point of view of the new tetrahedron
	// the 4th node (n3) must be allocated and supplied, and it's coordinates will be replaced
	private static double SQRT6D3 = 0.8164965809277259;
	public FEM1Element optimalTetrahedron(int n0, int n1, int n2, int n3) {
		int n0_3 = n0 * FEM1.NCOORD, n1_3 = n1 * FEM1.NCOORD, n2_3 = n2 * FEM1.NCOORD, n3_3 = n3 * FEM1.NCOORD;
		double[] node0, node1, node2, node3;
		if (n0_3 >= node.length) { node0 = nodeWork; n0_3 -= node.length; } else node0 = node;
		if (n1_3 >= node.length) { node1 = nodeWork; n1_3 -= node.length; } else node1 = node;
		if (n2_3 >= node.length) { node2 = nodeWork; n2_3 -= node.length; } else node2 = node;
		if (n3_3 >= node.length) { node3 = nodeWork; n3_3 -= node.length; } else node3 = node;
		// get normal from vectors n1-n0, n2-n0
		double x10 = node1[n1_3++] - node0[n0_3++], y10 = node1[n1_3++] - node0[n0_3++], z10 = node1[n1_3--] - node0[n0_3--]; n1_3--; n0_3--;
		double x20 = node2[n2_3++] - node0[n0_3++], y20 = node2[n2_3++] - node0[n0_3++], z20 = node2[n2_3--] - node0[n0_3--]; n2_3--; n0_3--;
		double xn = y20 * z10 - z20 * y10, yn = z20 * x10 - x20 * z10, zn = x20 * y10 - y20 * x10;
		double nl = 1. / Math.sqrt(xn*xn + yn*yn + zn*zn); xn *= nl; yn *= nl; zn *= nl;
		// get midpoint of triangle formed by n0, n1, n2
		double xMf = (node0[n0_3++] + node1[n1_3++] + node2[n2_3++])*0.3333333;
		double yMf = (node0[n0_3++] + node1[n1_3++] + node2[n2_3++])*0.3333333, zMf = (node0[n0_3] + node1[n1_3] + node2[n2_3])*0.3333333;
		// insert calculated coordinates into node3
		node3[n3_3++] = xMf + xn * SQRT6D3; node3[n3_3++] = yMf + yn * SQRT6D3; node3[n3_3++] = zMf + zn * SQRT6D3;
		// return instantiated tetrahedron, give it same material property as spawner
		FEM1Element elem = new FEM1Element(this, n0, n1, n2, n3, null);
		return elem;
	}
	

	
	double[] facetNormal(int n1, int n2, int n3) {
		int n1_3 = n1 * 3, n2_3 = n2 * 3, n3_3 = n3 * 3;
		double x21, y21, z21, x31, y31, z31;
		double[] node2, node3, normal = new double[3];
		if (n2_3 >= node.length) { n2_3 -= node.length; node2 = nodeWork; } else node2 = node;
		if (n3_3 >= node.length) { n3_3 -= node.length; node3 = nodeWork; } else node3 = node;
		if (n1_3 < node.length) {
			x21 = node2[n2_3++]-node[n1_3++]; y21 = node2[n2_3++]-node[n1_3++]; z21 = node2[n2_3]-node[n1_3];
			x31 = node3[n3_3++]-node[n1_3++]; y31 = node3[n3_3++]-node[n1_3++]; z31 = node3[n3_3]-node[n1_3];
		} else {
			n1_3 -= node.length; x21 = node2[n2_3++]-nodeWork[n1_3++]; y21 = node2[n2_3++]-nodeWork[n1_3++]; z21 = node2[n2_3]-nodeWork[n1_3];
			n3_3 -= node.length; x31 = node3[n3_3++]-nodeWork[n1_3++]; y31 = node2[n3_3++]-nodeWork[n1_3++]; z31 = node2[n3_3]-nodeWork[n1_3];
		}
		double xn = y21*z31 - z21*y31, yn = x21*z31 - z21*x31, zn = x21*y31 - y21*x31, nLinv = 1. / Math.sqrt(xn*xn + yn*yn + zn*zn);
		normal[0] = xn * nLinv; normal[1] = yn * nLinv; normal[2] = zn * nLinv;
		return normal;
	}
	
	
	
	// method calculates the optimal position of an internal node within it's 1-ring neighbourhood, using the formula of
	// the weighted barycenter of circumcenters divided by total 1-ring volume
	// TODO: debug this method
	double[] nodeOptimal1RingPosition(int n) {
		int[] elSup = elementSupports[n];
		if (elSup.length < 2) return null;					// at least 2 elements needed
		double[] optPos = new double[3], fNormalSum = new double[3];
		double volume1ring = 0, volume1ringI;
		boolean isExternal = externalNode(n);
		
		double sumB = 0;									// boundary facet sum
		
		for (int e = 0; e < elSup.length; e++) {
			FEM1Element elem = getElement2(elSup[e]);
			if (elem == null) continue;
			volume1ring += elem.tetraVolume(true);
			double[] cCenter = elem.tetraCircumcenter();
			optPos[0] += cCenter[0] * elem.volume; optPos[1] += cCenter[1] * elem.volume; optPos[2] += cCenter[2] * elem.volume;
			
			if (isExternal) {								// if node is bordering some boundary facets (perhaps it does on current element)			
				// need to find all the external facets attached to this node by exclusion: flag the neighbouring facets and find the unflagged one/ones
				int internalAreas = 0;						
				for (int i = 1, neighboursF2 = elem.neighboursF * 2; i < neighboursF2; i += 2) internalAreas |= elem.neighbour[i];

				if (internalAreas < 15) {					// precalculate data on tetrahedron only if it has facet/facets on the boundary
					elem.evaluateAreas(false);				// we will need the facet areas
					elem.evaluateEdges();					// we will need the edge lengths
					elem.propagateInterfaces(this, FEM1Element.PROPAGATE_AREAS | FEM1Element.PROPAGATE_EDGES, false);

					int[] nodeRef = elem.nodeRef;
					double[] fNormal1 = {0,0,0}, fNormal2 = {0,0,0}, fNormal3 = {0,0,0};
					// get sums from the found boundary facets, depending on what node n corresponds to on tetrahedron
					// fNormal1/2/3 is the inward-pointing normal of the boundary facet, scaled by area of that facet, N(i,p,q)
					// this sums up boundary factors N(i,p,q) * (||x(p)-x(i)||^2 - ||x(q)-x(i)||^2) from every neighbour's boundary facet
					if (nodeRef[0] == n) {
						int ref0 = nodeRef[0];
						if ((internalAreas&1) == 0) { double sum = (elem.data[EDGE01] + elem.data[EDGE02]) * elem.data[FACE012];	// facet 012
							fNormal1 = facetNormal(ref0, nodeRef[2], nodeRef[1]); fNormal1[0] *= sum; fNormal1[1] *= sum; fNormal1[2] *= sum; }
						if ((internalAreas&2) == 0) { double sum = (elem.data[EDGE02] + elem.data[EDGE03]) * elem.data[FACE023];	// facet 023
							fNormal2 = facetNormal(ref0, nodeRef[3], nodeRef[2]); fNormal2[0] *= sum; fNormal2[1] *= sum; fNormal2[2] *= sum; }
						if ((internalAreas&4) == 0) { double sum = (elem.data[EDGE01] + elem.data[EDGE03]) * elem.data[FACE031];	// facet 031
							fNormal3 = facetNormal(ref0, nodeRef[1], nodeRef[3]); fNormal3[0] *= sum; fNormal3[1] *= sum; fNormal3[2] *= sum; }
					} else if (nodeRef[1] == n) {
						int ref1 = nodeRef[1];
						if ((internalAreas&1) == 0) { double sum = (elem.data[EDGE01] + elem.data[EDGE12]) * elem.data[FACE012];	// facet 012
							fNormal1 = facetNormal(nodeRef[0], nodeRef[2], ref1); fNormal1[0] *= sum; fNormal1[1] *= sum; fNormal1[2] *= sum; }
						if ((internalAreas&4) == 0) { double sum = (elem.data[EDGE01] + elem.data[EDGE13]) * elem.data[FACE031];	// facet 031
							fNormal2 = facetNormal(nodeRef[0], ref1, nodeRef[3]); fNormal2[0] *= sum; fNormal2[1] *= sum; fNormal2[2] *= sum; }
						if ((internalAreas&8) == 0) { double sum = (elem.data[EDGE13] + elem.data[EDGE12]) * elem.data[FACE132];	// facet 132
							fNormal3 = facetNormal(ref1, nodeRef[2], nodeRef[3]); fNormal3[0] *= sum; fNormal3[1] *= sum; fNormal3[2] *= sum; }
					} else if (nodeRef[2] == n) {
						int ref2 = nodeRef[2];
						if ((internalAreas&1) == 0) { double sum = (elem.data[EDGE02] + elem.data[EDGE12]) * elem.data[FACE012];	// facet 012
							fNormal1 = facetNormal(nodeRef[0], ref2, nodeRef[1]); fNormal1[0] *= sum; fNormal1[1] *= sum; fNormal1[2] *= sum; }
						if ((internalAreas&2) == 0) { double sum = (elem.data[EDGE02] + elem.data[EDGE23]) * elem.data[FACE023];	// facet 023
							fNormal2 = facetNormal(nodeRef[2], nodeRef[3], ref2); fNormal2[0] *= sum; fNormal2[1] *= sum; fNormal2[2] *= sum; }
						if ((internalAreas&8) == 0) { double sum = (elem.data[EDGE23] + elem.data[EDGE12]) * elem.data[FACE132];	// facet 132
							fNormal3 = facetNormal(nodeRef[1], ref2, nodeRef[3]); fNormal3[0] *= sum; fNormal3[1] *= sum; fNormal3[2] *= sum; }
					} else /*if (nodeRef[3] == n)*/ {
						int ref3 = nodeRef[3];
						if ((internalAreas&2) == 0) { double sum = (elem.data[EDGE23] + elem.data[EDGE03]) * elem.data[FACE023];	// facet 023
							fNormal1 = facetNormal(nodeRef[0], ref3, nodeRef[2]); fNormal1[0] *= sum; fNormal1[1] *= sum; fNormal1[2] *= sum; }
						if ((internalAreas&4) == 0) { double sum = (elem.data[EDGE03] + elem.data[EDGE13]) * elem.data[FACE031];	// facet 031
							fNormal2 = facetNormal(nodeRef[0], nodeRef[1], ref3); fNormal2[0] *= sum; fNormal2[2] *= sum; fNormal2[2] *= sum; }
						if ((internalAreas&8) == 0) { double sum = (elem.data[EDGE13] + elem.data[EDGE23]) * elem.data[FACE132];	// facet 132
							fNormal3 = facetNormal(nodeRef[1], nodeRef[2], ref3); fNormal3[0] *= sum; fNormal3[3] *= sum; fNormal3[2] *= sum; }
					}				
					fNormalSum[0] += fNormal1[0] + fNormal2[0] + fNormal3[0];
					fNormalSum[1] += fNormal1[1] + fNormal2[1] + fNormal3[1];
					fNormalSum[2] += fNormal1[2] + fNormal2[2] + fNormal3[2];
				}
				
			}
		}
		volume1ringI = 6.0 / volume1ring;
		sumB *= (1./6.);										// divide by 6
		optPos[0] += sumB * fNormalSum[0]; optPos[1] += sumB * fNormalSum[1]; optPos[2] += sumB * fNormalSum[2];
		optPos[0] *= volume1ringI; optPos[1] *= volume1ringI; optPos[2] *= volume1ringI;

		return optPos;
	}
	
	
	
	// method generates first tetrahedron of FEM1 instance as construction seed, FEM1 expected to contain a PSC (Piecewise Smooth Complex) of facets
	FEM1Element firstTetrahedron() {
		if (borderPolygons == null) facetNeighbours(true, true);
		
		maxTsize = 0;														// do simple average edge length calculation
		for (int e = 0; e < polygonEdge.length; e++) maxTsize += polygonEdge[e];
		maxTsize /= (double)polygonEdge.length;
		
		element2 = new FEM1Element[64];										// start at 64 tetrahedrons as default
		int pBo = polygonOffset[borderPolygons[0]];							// just pick first available border facet as basis
		int i0 = polygonEdgeIndex[pBo], i1 = polygonEdgeIndex[pBo+1], i2 = polygonEdgeIndex[pBo+2];
		// find shortest edge of this polygon, using the nonredundant reference array polygonEdgeIndex[]
		double shortestE01 = polygonEdge[i0] < polygonEdge[i1] ? polygonEdge[i0] : polygonEdge[i1];
		double shortestE = shortestE01 < polygonEdge[i2] ? shortestE01 : polygonEdge[i2];
		
		if (shortestE > maxTsize) {											// is seed facet still overall larger than maximum allowed tetrahedron size?
			
		} else {
			
		}
		return element2[0];
	}
	
	
	
	
	// method quicksorts supplied indexed doubles according to increasing value, also creating a reindexing table called sortI[]
	public void sortIndexedDoubles(double[] sortD, int[] sortI, int[] idxI) {
				
		double tmp, d3 = 1. / 3.;
		int rangeC = 2, rangeCN = 0, tmpI, values = sortD.length;
		int[] stack = new int[2];		// stack will hold the ranges of the quicksort subdomains
		int[] stackN = new int[4];		// stackN is the next BFS stack (since the approach is equal to a Breadth First Search expansion)
		stack[1] = values - 1;
		
		while(true) {
			for (int r = 0, rN = 0; r < rangeC;) {
				int rS = stack[r++], rE = stack[r++];								// rS = start idx, rE = end idx
				double pivot = (rE - rS > 2) ? (sortD[rS] + sortD[rS + (rE - rS)/2] + sortD[rE]) * d3 : sortD[rS + 1];
				int rML = rS, rMR = rE;												// rML/rMR = left & right indexes of elements equal to pivot

				while (true) {
					while (sortD[rML] < pivot) rML++;
					while (sortD[rMR] > pivot) rMR--;
					if (rML >= rMR) break;
					tmp=sortD[rML]; sortD[rML]=sortD[rMR]; sortD[rMR]=tmp;			// swap found values (even if equal)
					tmpI=sortI[rML]; sortI[rML++]=sortI[rMR]; sortI[rMR--]=tmpI;	// and progress to next element pair
				}
				if (rML > rMR) { rML = rMR; rMR++; } else rMR++;					// correct left & right pivot delimiters
				// these two lines cut out any equal elements around the pivot
				while (sortD[rML] == pivot && rS < rML) rML--;						// eliminate equal pivot elements in left part of subdomain
//				int rMR2 = rMR;
//				while (sortD[rMR] == pivot && rMR < rE) rMR++;						// leave only one element = pivot in right part of subdomain
//				if (rMR2 < rMR) rMR--;
				
				if (rS < rML) {
					if (rS + 1 == rML) {											// only 2 elements left base-case, just swap
						if (sortD[rS] > sortD[rML]) {
							tmp = sortD[rS]; sortD[rS] = sortD[rML]; sortD[rML] = tmp;
							tmpI =sortI[rS]; sortI[rS] = sortI[rML]; sortI[rML] = tmpI; } 
					} else {	stackN[rN++] = rS; stackN[rN++] = rML;				// store left subpartition for next stack
								rangeCN += 2; }
				}
				if (rMR < rE) {
					if (rMR + 1 == rE) {											// only 2 elements left base-case, just swap
						if (sortD[rMR] > sortD[rE]) {
							tmp = sortD[rE]; sortD[rE] = sortD[rMR]; sortD[rMR] = tmp;
							tmpI =sortI[rE]; sortI[rE] = sortI[rMR]; sortI[rMR] = tmpI; } 
					} else {	stackN[rN++] = rMR; stackN[rN++] = rE; 				// store right subpartition for next stack
								rangeCN += 2; }
				}
			}
			if (rangeCN == 0) break;												// if all ranges have been subdivided & sorted, stop
			stack = stackN;
			stackN = new int[rangeCN * 2];											// worst-case next-stack will be twice the size
			rangeC = rangeCN;
			rangeCN = 0;
		}
		
		for (int v = 0; v < values; v++) idxI[sortI[v]] = v;						// make proper original -> mutated index reindexer		
		// DEBUG check propriety of sorting
		//for (int n = 0; n < nodes - 1; n++) if (sortD[n] > sortD[n+1]) System.out.println("Faulty sort on nodes " + n + " & " + (n+1));
	}
	

	// method utilises an octree and closest node octree search to fill three supplied arrays:
	// closest[] (nodeIdx1,nodeIdx2)-tuple array, distance[] array and closestOctants[] (octantPtr1,octantPtr2)-tuple array
	// creating an indexation of all closest-node pairs, their distances and their container octants
	public int closestNodePairs(FEM1Octree octree, int[] closest, double[] distance, FEM1Octree[] closestOctants) {
		double[] distC = {0};
		FEM1Octree[] octantsC = {null, null};
		
		// flag every node tuple as "unassigned", for skipping symmetric closeness cases (n1 closest to n2 & n2 closest to n1 -> do not store (n2,n1))
		for (int n2 = 0, nodes2 = nodes * 2; n2 < nodes2; n2++) closest[n2++] = -1;
		
		int nPairs = 0;
		FEM1Octree[] leafOctant = octree.leafOctantArray();							// collect leaf octants

		for (FEM1Octree octant : leafOctant) {
			if (octant == null) break;
			octantsC[0] = octant;
			for (int nI = 0; nI < octant.nodes; nI++) {								// iterate over total number of nodes in octant
				octantsC[1] = null;
				int n = octant.node[nI], n2 = n * 2;
				int nClosest = octree.closestNode(this, n, distC, octantsC);
				if (closest[nClosest * 2] == -1 || closest[nClosest * 2 + 1] != n) {// if closest node's closest node isn't this one (symmetric closeness)
					closest[n2] = n; closest[n2 + 1] = nClosest;					// then add it to closest pairs register
					if (distance != null) distance[n] = distC[0];
					if (closestOctants != null) { closestOctants[n2] = octantsC[0]; closestOctants[n2 + 1] = octantsC[1]; }
					nPairs++;
				}
			}
		}

		for (int n2 = 0, d = 0, n2b = 0, nodes2 = nodes * 2; n2 < nodes2;) {		// remove gaps of symmetric pairs in arrays
			if (closest[n2] != -1) {
				if (distance != null) distance[d++] = distance[n2 >> 1];
				closest[n2b++] = closest[n2++]; closest[n2b--] = closest[n2--];
				if (closestOctants != null) { closestOctants[n2b++] = closestOctants[n2++]; closestOctants[n2b++] = closestOctants[n2++]; }
				else { n2 += 2; n2b +=2; }
			} else n2 += 2;
		}
		return nPairs;																// return resultant number of node pairs
	}
	
	
	// method calculates median points of every closest-node pair, adding them into the octree on top of the true nodes
	// method will return a coordinate array of the calculated medians and the containing octants (if an array is supplied)
	// method utilises the extraneous nodes space in octree for quick addition without necessity of tree reconstruction
	// reallocation will happen if extraneous space isn't large enough
	double[] addMediansOfClosestNodes(FEM1Octree octree, int nPairs, int[] closest, FEM1Octree[] closestOctants, FEM1Octree[] medianOctants) {
		
		double[] median = new double[nPairs * 3];
		
		// nP2 = node's neighbour index, np1_3 & npo2_3 = node's & neighbour node's indexes into coord.array, nM3 = the running coord.index of added medians
		for (int nP = 0, nP2 = 1, nP1_3 = 0, nM3 = 0; nP < nPairs; nP++, nP2 += 2, nP1_3 += 3) {
			int nP2_3 = closest[nP2] * 3;
			double xm, ym, zm;												// the median coordinate
			if (nP2_3 >= node.length) { nP2_3 -= node.length; xm = nodeWork[nP2_3++]; ym = nodeWork[nP2_3++]; zm = nodeWork[nP2_3]; }
			else { xm = node[nP2_3++]; ym = node[nP2_3++]; zm = node[nP2_3]; }
			if (nP1_3 >= node.length) {
				nP1_3 -= node.length;
				xm = median[nM3++] = nodeWork[nP1_3] + (xm - nodeWork[nP1_3++]) * .5;
				ym = median[nM3++] = nodeWork[nP1_3] + (ym - nodeWork[nP1_3++]) * .5;
				zm = median[nM3++] = nodeWork[nP1_3] + (zm - nodeWork[nP1_3++]) * .5;
			} else {
				xm = median[nM3++] = node[nP1_3] + (xm - node[nP1_3++]) * .5;
				ym = median[nM3++] = node[nP1_3] + (ym - node[nP1_3++]) * .5;
				zm = median[nM3++] = node[nP1_3] + (zm - node[nP1_3++]) * .5;
			}
			
			// find out if median belongs to the two (or one) octants of the closest nodes
			FEM1Octree octant = null;
			if (closestOctants[nP2].octantDistance(xm, ym, zm, false) == 0)	{
				octant = closestOctants[nP2];
				octant.octantAddNode(xm, ym, zm, nM3 / 3, true);
			} else if (closestOctants[nP2 + 1].octantDistance(xm, ym, zm, false) == 0) {
				octant = closestOctants[nP2 + 1];
				octant.octantAddNode(xm, ym, zm, nM3 / 3, true);
			}
			// the case of the median crossing over to a corner octant to the two neighbour octants:
			// find the topmost container and do an extraneous node insertion
			else	octant = octree.addNode(this, xm, ym, zm, nM3 / 3, octreeMaxItems, true);	
			
			if (medianOctants != null) medianOctants[nP] = octant;
		}
		return median;
	}

		
	// method returns the node in supplied index array that is closest to supplied node index, distance returned in supplied distance[] array
	// utilising the possible existence of additional nodes in nodeWork[]
	// note: comparison happens with squared distances, and the distance is returned squared!
	int closestNode(int n1, int[] nL, double[] distance) {	
		int n1_3 = n1 * 3;
		double[] node1;
		if (n1_3 >= node.length) { node1 = nodeWork; n1_3 -= node.length; } else node1 = node;	
		return closestNode(node1[n1_3++], node1[n1_3++], node1[n1_3], nL, distance);
	}

	// method returns the node in supplied array that is closest to supplied position, returns node index and distance in returned distance[] array
	// note: comparison happens with squared distances, and the distance is returned squared!
	int closestNode(double n1x, double n1y, double n1z, int[] nL, double[] distance) {
		int nClosest = -1;
		double[] node2;
		double dClosest = Double.MAX_VALUE;
		for (int n2 = 0, n2_3 = 0; n2 < nL.length; n2++) {
			n2_3 = nL[n2] * 3;
			if (n2_3 >= node.length) { node2 = nodeWork; n2_3 -= node.length; } else node2 = node;
			double d12x = n1x - node2[n2_3++], d12y = n1y - node2[n2_3++], d12z = n1z - node2[n2_3++];
			double d12 = d12x*d12x + d12y*d12y + d12z*d12z;
			// if position is not equal to a node (meaning, a node position compared to itself) and closer than previous closest
			if (d12 != 0 && d12 < dClosest) { dClosest = d12; nClosest = nL[n2]; }
		}
		if (distance != null) distance[0] = dClosest;
		return nClosest;
	}

	
	// methof calculates distance between nodes, utilising the possible existence of additional nodes in nodeWork[]
	double distance(int n1, int n2, boolean root) {
		double[] node1, node2;
		if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
		if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
		int n3_1 = n1 * 3, n3_2 = n2 * 3;
		double d12x, d12y, d12z;
		if (n3_1 < n3_2) {	d12x = node1[n3_1++]; d12y = node1[n3_1++]; d12z = node1[n3_1];
							d12x -= node2[n3_2++]; d12y -= node2[n3_2++]; d12z -= node2[n3_2]; }
		else {				d12x = node2[n3_2++]; d12y = node2[n3_2++]; d12z = node2[n3_2];
							d12x -= node1[n3_1++]; d12y -= node1[n3_1++]; d12z -= node1[n3_1]; }	
		if (root)	return Math.sqrt(d12x * d12x + d12y * d12y + d12z * d12z);
		else		return d12x * d12x + d12y * d12y + d12z * d12z;
	}

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			STIFFNESS MATRIX ASSEMBLER 1
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	

	static final double DIV6 = 1./6.;
	
	// global stiffness matrix assembler according to Long Chen, "Programming of Finite Elements Method in Matlab", adapted for Java
	// method assembles the tetrahedral simplexes into an NSP sparse matrix
	public int[] assemble(int[] ii, int[] jj, double[] sA, boolean writeRxC) {
		
		int elements4 = elements * 4;
		
		//face = null;
		if (polygon == null) {											// if the input file didn't have tetrahedrons already defined into faces
			polygon = new int[elements4 * 3];
			polygons = elements4;
			for (int f = 0, e4 = 0; e4 < elements4;) {					// gather faces into the 4*elements*3 face array
				int t1 = element[e4++], t2 = element[e4++], t3 = element[e4++], t4 = element[e4++];
				int[] fc;
				polygon[f++] = t2; polygon[f++] = t4; polygon[f++] = t3;
				polygon[f++] = t1; polygon[f++] = t3; polygon[f++] = t4;
				polygon[f++] = t1; polygon[f++] = t4; polygon[f++] = t2;
				polygon[f++] = t1; polygon[f++] = t2; polygon[f++] = t3;
			}
		}
		
		double[] v12 = new double[elements4 * 3], v13 = new double[elements4 * 3];
		double aN[] = new double[elements4 * 3];
		// v12 = node(face(:, 2), :) - node(face(:, 1), :);
		// v13 = node(face(:, 3), :) - node(face(:, 1), :);
		int polygons3 = polygons * 3;
		for (int f3 = 0, v3a = 0, v3b = 0, v3aN = 0; f3 < polygons3; f3++) {			// from the faces, form the edge vectors v12 & v13
			int fv1 = polygon[f3++] * 3, fv2 = polygon[f3++] * 3, fv3 = polygon[f3++] * 3;
			double nfv11 = node[fv1++], nfv12 = node[fv1++], nfv13 = node[fv1++];
			double a1 = v12[v3a++] = node[fv2++] - nfv11;
			double a2 = v12[v3a++] = node[fv2++] - nfv12;
			double a3 = v12[v3a++] = node[fv2++] - nfv13;
			double b1 = v13[v3b++] = node[fv3++] - nfv11;
			double b2 = v13[v3b++] = node[fv3++] - nfv12;
			double b3 = v13[v3b++] = node[fv3++] - nfv13;
			// aN = crossproduct(v12, v13, 2)
			aN[v3aN++] = a2 * b3 - a3 * b2; aN[v3aN++] = a3 * b1 - a1 * b3; aN[v3aN++] = a1 * b2 - a2 * b1;
		}
		
		double[][] normal = new double[elements][12];
		// normal(1:NT, :, 1) = aN(1:NT, :);
		// normal(1:NT, :, 2) = aN(NT+1:2*NT, :);
		// normal(1:NT, :, 3) = aN(2*NT+1:3*NT, :);
		// normal(1:NT, :, 4) = aN(3*NT+1:4*NT, :);
		// DEBUG: seems like one needs to pick the interspersed v12xv13 crossproducts, into normal[]
		for (int e = 0, e1 = 0, e2 = 3, e3 = 6, e4 = 9; e < elements; e++) {
			double[] normal_e = normal[e];
			normal_e[0] = aN[e1++]; normal_e[1] = aN[e2++]; normal_e[2] = aN[e3++]; normal_e[3] = aN[e4++]; 
			normal_e[4] = aN[e1++]; normal_e[5] = aN[e2++]; normal_e[6] = aN[e3++]; normal_e[7] = aN[e4++]; 
			normal_e[8] = aN[e1]; normal_e[9] = aN[e2]; normal_e[10] = aN[e3]; normal_e[11] = aN[e4];
			e1 += 10; e2 += 10; e3 += 10; e4 += 10;
		}
//		int elements3 = elements * 3;
//		for (int e = 0, e1 = 0, e2 = elements3, e3 = elements3 * 2, e4 = elements3 * 3; e < elements; e++) {
//			// note: test if changing increment order speeds up loop
//			double[] normal_e = normal[e];
//			normal_e[0] = aN[e1++]; normal_e[4] = aN[e1++]; normal_e[8] = aN[e1++];
//			normal_e[1] = aN[e2++]; normal_e[5] = aN[e2++]; normal_e[9] = aN[e2++];
//			normal_e[2] = aN[e3++]; normal_e[6] = aN[e3++]; normal_e[10] = aN[e3++];
//			normal_e[3] = aN[e4++]; normal_e[7] = aN[e4++]; normal_e[11] = aN[e4++];
//		}
		
		// v12 = v12(3*NT+1:4*NT, :); v13 = v13(3*NT+1:4*NT, :);
		// cross(v12, v13, 2) (into v12)
		// DEBUG: seems like one needs to pick the interspersed, first face crossproducts, of every tetrahedron
		for (int ea = 0, eb = ea, e2 = ea, eaEnd = elements4 * 3; ea < eaEnd;) {
			double a1 = v12[ea++], a2 = v12[ea++], a3 = v12[ea++], b1 = v13[eb++], b2 = v13[eb++], b3 = v13[eb++];
			v12[e2++] = a2 * b3 - a3 * b2; v12[e2++] = a3 * b1 - a1 * b3; v12[e2++] = a1 * b2 - a2 * b1;			
			ea += 9; eb += 9;
		}
//		for (int ea = elements3 * 3, eb = ea, e2 = ea, eEnd = ea + elements3; ea < eEnd;) {
//			double a1 = v12[ea++], a2 = v12[ea++], a3 = v12[ea++], b1 = v13[eb++], b2 = v13[eb++], b3 = v13[eb++];
//			v12[e2++] = a2 * b3 - a3 * b2; v12[e2++] = a3 * b1 - a1 * b3; v12[e2++] = a1 * b2 - a2 * b1;			
//		}
		
		// v14 = node(elem(:, 4), :) - node (elem(:, 1), :);
		// volume = abs(dot(cross(v12, v13, 2), v14, 2)) / 6;
		double[] volume = new double[elements];
		for (int e = 0, v3 = 0, e1 = 0, e4 = 3; e < elements; e++, e1 += 4, e4 += 4) {
			int n4_3 = element[e4] * 3, n1_3 = element[e1] * 3;							// this is elem(:, 4) and elem(:, 1)
			volume[e] =  v12[v3++] * (node[n4_3++] - node[n1_3++]);
			volume[e] += v12[v3++] * (node[n4_3++] - node[n1_3++]);
			volume[e] += v12[v3++] * (node[n4_3] - node[n1_3]);
			volume[e] *= volume[e] < 0 ? -DIV6 : DIV6;
		}
		
		int[] rCounts = null;	
		
		if (writeRxC) {
			for (int i = 0, idx = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					for (int e = 0, e4 = 0; e4 < elements4; e++, e4 += 4) {
						idx = element[e4 + i] * nodes + element[e4 + j];
						double volD = 1.0 / (volume[e] * 36.0);
						int ni = i, nj = j;
						// sA(idx + 1 : idx + elements) = dot(normal(:, :, i), normal(:, :, j), 2) ./ (36 * volume)
						double[] normal_e = normal[e];
						sA[idx] += normal_e[ni] * normal_e[nj] * volD;
						sA[idx] += normal_e[ni += 4] * normal_e[nj += 4] * volD;
						sA[idx] += normal_e[ni += 4] * normal_e[nj += 4] * volD;
					}
				}
			}			
		} else {
			BinBitImage bitImage = new BinBitImage(new NSPMatrix("", nodes, nodes));	// bitImage will hold flags for nonzeroes
			rCounts = new int[nodes];													// rCounts will hold nonzeroes counts of every matrix row
			
			for (int i = 0, idx = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					for (int e = 0, e4 = 0; e4 < elements4; e++, e4 += 4) {
						int i2 = ii[idx] = element[e4 + i];
						int j2 = jj[idx] = element[e4 + j];
						if (!bitImage.bitIsSet(i2, j2)) {			// if this is an original nonzero entry, increment proper row of rCounts
							rCounts[i2]++;
							bitImage.setBit(i2, j2);
						}
						double volD = 1.0 / (volume[e] * 36.0);
						int ni = i, nj = j;
						// sA(idx + 1 : idx + elements) = dot(normal(:, :, i), normal(:, :, j), 2) ./ (36 * volume)
						double[] normal_e = normal[e];
						sA[idx] = normal_e[ni] * normal_e[nj] * volD;
						sA[idx] += normal_e[ni += 4] * normal_e[nj += 4] * volD;
						sA[idx++] += normal_e[ni += 4] * normal_e[nj += 4] * volD;
					}
				}
			}
		}
		return rCounts;
	}
	
	
	
	public NSPMatrix assembleNSPMatrix() {
		
		// ii[] holds row indexes, jj[] holds column indexes, sA[] holds data, time to assemble the sparse matrix
		// do search runs through column indexes in jj[] for every column we want to add, since appending from right
		// is the quickest way to assemble the sparse matrix		
		int[] ii = new int[elements * 16], jj = new int[elements * 16];
		double[] sA = new double[16 * elements];
		int[] rCounts = this.assemble(ii, jj, sA, false);

		NSPMatrix M = new NSPMatrix(this.name + "_NSP", nodes, nodes);
		int elements16 = elements * 16;
		int[] cOffsets = new int[nodes];								// hold column offset counts here
		
		for (int j = 0; j < nodes; j++)	{								// iterate over sought columns
			int purged = 0;
			boolean purge = false;
			for (int e = 0, eW = 0; e < elements16; e++) {
				if (jj[e] == j) {										// found the column we're currently assembling?
					int i = ii[e], cOffP = cOffsets[i] - 1;
					NspArray aHsp = M.Hsp[i];
					if (aHsp.array == null) {
						aHsp.array = new NspNode[rCounts[i]];
						aHsp.nodes = aHsp.size = rCounts[i];
					}
					if (cOffP >= 0 && aHsp.array[cOffP].c() == j) {
						aHsp.array[cOffP].add(sA[e]);
					} else {
						NspNode node = aHsp.array[cOffsets[i]] = new NspNode(i, j, sA[e], 0, cOffsets[i]++, 0);
						M.nNZ++;
						if (i == j) M.pivotNsp[i] = node;
					}
					
					purged++;
					if (!purge) {
						eW = e;												// since we've read a value out of ii[], jj[] & sA[],
						purge = true;										// set the purging flag & position
					}
				} else
					if (purge) { ii[eW] = ii[e];  jj[eW] = jj[e];  sA[eW++] = sA[e]; } 
			}
			elements16 -= purged;
		}
		
		M.crosslinkVsp();	
		return M;
	}
	
	// the assembleMatrix() method requests the FEM1 assembler to write straight into the allocated matrix array, supplied through getDataRef()
	public Matrix assembleMatrix() {		
		int[] ii = new int[elements * 16], jj = new int[elements * 16];
		Matrix M = new Matrix(this.name, nodes, nodes, Matrix.Type.Null, 1);
		this.assemble(ii, jj, M.getDataRef()[0], true);
		return M;
	}
		
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			FINITE ELEMENT ASSEMBLER 2 (code adapted from Allan F. Bower, solidmechanics.org)
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	

	// NDOF = degrees of freedom, NCOORDS = 3 coordinates for 3D, OFFSPRECOMP = an offset for precomputed arrays, MAXNODES = max. integration nodes used
	static final int NDOF = 3, NCOORD = 3, OFFSPRECOMP = 3, MAXNODES = 20;
	static final double DIV4 = 1./8., DIV8 = 1./8.;
	// arrays return number of integration points for specific element's node count, by direct indexation
	static final int[] numIntegPoints = {0,0,0,0,1,0,0,0,8,0,4,0,0,0,0,0,0,0,0,0,27};
	static final int[] computedArrayFromNodes = {0,0,0,0,0,0,0,0,1,0,2,0,0,0,0,0,0,0,0,0,3};
	static final int[] numIntegPoints2D = {0,0,0,1,4,0,3,0,9};
	static final int[] numIntegPointsMass = {0,0,0,0,4,0,0,0,27,0,5,0,0,0,0,0,0,0,0,0,27};
	
	static final int[] nFaceNodes = {0, 0, 0, 0, 3, 0, 0, 0, 4, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8};
	// array returns indexation of the nodes pertinent to different faces, depending on node count
	static final int[][][] faceNodes = {{},{},{},{},
			{{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {2, 3, 0}}, {}, {}, {},
			{{0, 1, 2, 3}, {4, 7, 6, 5}, {0, 4, 5, 1}, {1, 2, 6, 5}, {2, 6, 7, 3}, {3, 7, 4, 0}}, {},
			{{0, 1, 2, 4, 5, 6}, {0, 3, 1, 7, 8, 4}, {1, 3, 2, 8, 9, 5}, {2, 3, 0, 9, 7, 6}}, {}, {}, {}, {}, {}, {}, {}, {}, {},
			{{0, 1, 2, 3, 8, 9, 10, 11}, {4, 7, 6, 5, 15, 14, 13, 12}, {0, 4, 5, 1, 16, 12, 17, 8},
				{1, 5, 6, 2, 17, 13, 18, 9}, {2, 6, 7, 3, 18, 14, 19, 10}, {3, 7, 4, 0, 19, 15, 16, 11}}};

	// array integPoints[][][] returns by direct indexation the integration coordinates 1D-array for specific point & node counts
	double[][] integPoints_4_10 = {
		{},	{.25,.25,.25}, {}, {},
		{0.58541,0.1382,0.1382, 0.1382,0.58541,0.1382, 0.1382,0.1382,0.58541, 0.1382,0.1382,0.1382},
		{.25,.25,.25, .5,1./6.,1./6., 1./6.,.5,1./6., 1./6.,1./6.,.5, 1./6.,1./6.,1./6.}};
	double[][] integPoints_8_20 = {{}, {0,0,0}, {}, {}, {}, {}, {}, {},
		{-.57735,-.57735,-.57735,.57735,-.57735,-.57735,-.57735,.57735,-.57735,.57735,.57735,-.57735,
			-.57735,-.57735,.57735,.57735,-.57735,.57735,-.57735,.57735,.57735,.57735,.57735,.57735}, {}, {}, {}, {}, {},
		{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {},
		{-.7746,-.7746,-.7746,0,-.7746,-.7746,.7746,-.7746,-.7746,-.7746,0,-.7746,0,0,-.7746,.7746,0,-.7746,-.7746,.7746,-.7746,0,
			.7746,-.7746,.7746,.7746,-.7746,-.7746,-.7746,0,0,-.7746,0,.7746,-.7746,0,-.7746,0,0,0,0,0,.7746,0,0,-.7746,.7746,0,0,
			.7746,0,.7746,.7746,0,-.7746,-.7746,.7746,0,-.7746,.7746,.7746,-.7746,.7746,-.7746,0,.7746,0,0,.7746,.7746,0,.7746,-.7746,
			.7746,.7746,0,.7746,.7746,.7746,.7746,.7746}};
	// data layout: integPoints[nNodes][nPoints][(the points)]
	double[][][] integPoints = {{}, {}, {}, {}, integPoints_4_10, {}, {}, {}, integPoints_8_20, {},
								integPoints_4_10, {}, {}, {}, {}, {}, {}, {}, {}, {}, integPoints_8_20};
	
	// integration points over faces, data layout: integPoints_2D[nFaceN][nPoints][(the points)]
	double[][][] integPoints_2D = {{}, {}, {},
			{{}, {1./3.,1./3.}, {}, {.6,.2, .2,.6, .2,.2}, {1./3.,1./3., .6,.2, .2,.6, .2,.2}, {}, {},
				{1./3.,1./3., .059716,.47014, .47014,.059716, .47014,.47014, .79743,.10129, .10129,.79743, .10129,.10129}},
			{{}, {0,0}, {}, {}, {-.57735,-.57735, .57735,-.57735, -.57735,.57735, .57735,.57735}, {}, {}, {}, {},
				{-.7746,-.7746, 0,-.7746, .7746,-.7746, -.7746,0, 0,0, .7746,0, -.7746,.7746, 0,.7746, .7746,.7746}}, {},
			{{}, {1./3.,1./3.}, {}, {.6,.2, .2,.6, .2,.2}, {1./3.,1./3., .6,.2, .2,.6, .2,.2}, {}, {},
				{1./3.,1./3., .059716,.47014, .47014,.059716, .47014,.47014, .79743,.10129, .10129,.79743, .10129,.10129}}, {},
			{{}, {0,0}, {}, {}, {-.57735,-.57735, .57735,-.57735, -.57735,.57735, .57735,.57735}, {}, {}, {}, {},
					{-.7746,-.7746, 0,-.7746, .7746,-.7746, -.7746,0, 0,0, .7746,0, -.7746,.7746, 0,.7746, .7746,.7746}}};
	
	// array integWeights[][][] returns by direct indexation the mass integration coordinates 1D-array for specific point & node counts
	double[][] integWeights_4_10 = {{}, {1./6.}, {}, {}, {1./24.,1./24.,1./24.,1./24.}};
	double[][] integWeights_8_20 = {{}, {8}, {}, {}, {}, {}, {}, {}, {1,1,1,1,1,1,1,1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {},
			{}, {}, {}, {}, {}, {}, {}, {}, {.17147,.27435,.17147,.27435,.43896,.27435,.17147,.27435,.17147,.27435,.43896,.27435,
				.43896,.70233,.43896,.27435,.43896,.27435,.17147,.27435,.17147,.27435,.43896,.27435,.17147,.27435,.17147}};
	double[][][] integWeights = {{}, {}, {}, {}, integWeights_4_10, {}, {}, {}, integWeights_8_20,
								{}, integWeights_4_10, {}, {}, {}, {}, {}, {}, {}, {}, {}, integWeights_8_20};
	
	double[][][] integWeights_2D = {{}, {}, {},
			{{}, {.5}, {}, {1./6., 1./6., 1./6.}, {-27./96., 25./96., 25./96., 25./96.}},
			{{}, {4}, {}, {}, {1,1,1,1}, {}, {}, {}, {},
				{.30864, .49383, .30864, .49383, .79012, .49383, .30864, .49383, .30864}}, {},
			{{}, {.5}, {}, {1./6., 1./6., 1./6.}, {-27./96., 25./96., 25./96., 25./96.}}, {}, 
			{{}, {4}, {}, {}, {1,1,1,1}, {}, {}, {}, {},
				{.30864, .49383, .30864, .49383, .79012, .49383, .30864, .49383, .30864}}};
	
	
	
	// method constructs a mass matrix for a single finite element and expects an extra preallocated computed[][] array
	// computed[][] stores computations of xiList[], w[], computed point count, dNdxi[] (shape f.derivs), Jacobian inverses and Jacobians for reusage
	// shape function derivatives, determinants and Jacobians come sequentially stored in triples at end of array)
	// whichever calling method can check if they've been calculated already and thus can be reused
	// note: the precomputed data only useful if it is reused for same amount of integration points as the precomputing method's
	// therefore, the sameIntegCount = false will enforce full recomputation
	double[] elementMassMatrix(int nNodes, double coords[], double[] matProps, int elemFlags, boolean lumpedMass, double[][] computed) {
		
		// this is the (partial deriv. of shape function) / (partial deriv. of reference coordinate) matrix
		double[] dxdxi = new double[NCOORD * NCOORD];
		int nDOFn = NDOF * nNodes;
		double[] Mme = new double[nDOFn * nDOFn];						// this is the element mass matrix
		
		// setting up integration points & weights
		int nPoints = numIntegPointsMass[nNodes];
		boolean sameIntegCount = nPoints == computed[3][0] ? true : false;
		double[] xiList = computed[0] != null ? computed[0] : (computed[0] = integrationPoints(nNodes, nPoints));
		double[] w = computed[1] != null ? computed[1] : (computed[1] = integrationWeights(nNodes, nPoints));
		
		for (int intPt = 0, pF = 0, pD = 0, offsPrecomp = OFFSPRECOMP; intPt < nPoints; intPt++, offsPrecomp += 3) {
			
			double[] N = shapeFunctions(nNodes, xiList[pF++], xiList[pF++], xiList[pF++]);
			double det = 0;
		
			// get determinant of Jacobian, |J|, if it's been calculated, from Jacobian's inverse (Jacobian assumed to be already calculated)
			if (sameIntegCount && computed[offsPrecomp + 1] != null)
				det = computed[offsPrecomp + 1][9];							// determinant stored in 9th position
			
			// determinant not calculated (or integ.points mismatch), we must compute Jacobian for this point and calculate it's determinant
			else {
				if (sameIntegCount && computed[offsPrecomp + 2] != null) {
					dxdxi = computed[offsPrecomp + 2];						// retrieve precalculated Jacobian
				} else {
					double[] dNdxi = sameIntegCount && computed[offsPrecomp] != null ? 
							computed[offsPrecomp] : (computed[offsPrecomp] = shapeFunctionDerivs(nNodes, xiList[pD++], xiList[pD++], xiList[pD++]));

					// compute Jacobian & it's determinant
					for (int i = 0, aEnd = NCOORD * nNodes; i < NCOORD; i++) {
						int iC = i * NCOORD;
						for (int j = 0; j < NCOORD; j++) {
							int iCj = iC + j;
							dxdxi[iCj] = 0;
							for (int a = 0; a < aEnd; a += NCOORD)
								dxdxi[iCj] += coords[a + i] * dNdxi[a + j];
						}
					}
					computed[offsPrecomp + 2] = dxdxi;						// store Jacobian for this point sequentially in compute[][]
				}
				// calculate determinant (1st step of invert3x3dualCall())
				computed[offsPrecomp + 1] = Matrix.invert3x3dualCall(dxdxi, null);
				det = computed[offsPrecomp + 1][9];
			}
			
			double rho = matProps[3];										// get density constant of material
			double rho_w_det = rho * w[intPt] * det;
			
			for (int a = 0; a < nNodes; a++) {
				int rowOffs = nDOFn * a;
				double rho_Na_w_det = N[a] * rho_w_det;
				for (int i = 0, iRowOffs = rowOffs; i < NDOF; i++, iRowOffs += nDOFn) {
					for (int b = 0, offsetMme = iRowOffs + i; b < nNodes; b++, offsetMme += NDOF) {
							Mme[offsetMme] += rho_Na_w_det * N[b];
					}
				}
			}
		}
		
		// evaluate lumped mass using row summing method (rows summed into pivots, everything else zeroed)
		if (lumpedMass) {
			for (int a = 0; a < nDOFn; a++) {
				int rowOffs = nDOFn * a, pivotOffs = rowOffs + a;
				for (int b = rowOffs, bEnd = b + nDOFn; b < bEnd; b++)
					if (a != b) {
						Mme[pivotOffs] += Mme[b]; Mme[b] = 0;
					}
			}
		}		
		return Mme;
	}
	
	
	// method constructs a stiffness matrix for a single finite element and expects an extra preallocated computed[][] array
	// the usage of precomputations is equal to elementMassMatrix() method
	double[] elementStiffnessMatrix(int nNodes, double coords[], double[] dsde, int elemFlags, double[] displacement, boolean computeStrain, 
			double[][] computed) {
		
		int nPoints = numIntegPoints[nNodes], offsPrecomp = OFFSPRECOMP;
		// sameIntegCount = true if method is supplied with precalculated arrays for same no. of integ.points as we're calculating for this element
		// recalcJ flag stops recalculation of dNdxi[] and Jacobian if element is linear (tetrahedron/parallellogram), since Jacobian is a constant
		boolean sameIntegCount = nPoints == computed[2][0] ? true : false, recalcJ = true;
		double[] dNdx = new double[nNodes * NCOORD];
		double[] dxdxi = new double[NCOORD * NCOORD];
		double[] strain = new double[NDOF * NCOORD];
		int nDOFn = NDOF * nNodes;
		double[] K = new double[nDOFn * nDOFn], dxidx = null, dNdxi = null;
		
		double[] xiList = computed[0] != null ? computed[0] : (computed[0] = integrationPoints(nNodes, nPoints));
		double[] w = computed[1] != null ? computed[1] : (computed[1] = integrationWeights(nNodes, nPoints));
		
		for (int intPt = 0, pD = 0; intPt < nPoints; intPt++, offsPrecomp += 3) {
			
			//double[] N = shapeFunctions(nNodes, xiList[pF++], xiList[pF++], xiList[pF++]);
			double det = 0;
			
			// for a linear element, this whole Jacobian-processing block needs calculating only once
			if (recalcJ == true) {
				// check if shape function derivatives have been previously calculated for the same integer point count
				dNdxi = sameIntegCount && computed[offsPrecomp] != null ? 
						computed[offsPrecomp] :
						(computed[offsPrecomp] = shapeFunctionDerivs(nNodes, xiList[pD++], xiList[pD++], xiList[pD++]));
		
				// get determinant of Jacobian, |J|, if it's been calculated, from Jacobian's inverse (Jacobian assumed to be already calculated)
				if (sameIntegCount && computed[offsPrecomp + 1] != null) {
					dxidx = computed[offsPrecomp + 1];
					det = dxidx[9];												// determinant stored in 9th position
				
				// determinant not calculated (or integ.points mismatch), we must compute Jacobian for this point and calculate it's determinant
				} else {
					if (sameIntegCount && computed[offsPrecomp + 2] != null) {
						dxdxi = computed[offsPrecomp + 2];						// retrieve precalculated Jacobian
					} else {
								
						// compute Jacobian & it's determinant
						for (int i = 0, aEnd = NCOORD * nNodes; i < NCOORD; i++) {
							int iC = i * NCOORD;
							for (int j = 0; j < NCOORD; j++) {
								int iCj = iC + j;
								dxdxi[iCj] = 0;									// clear previous point's calculations
								for (int a = 0; a < aEnd; a += NCOORD)
									dxdxi[iCj] += coords[a + i] * dNdxi[a + j];
							}
						}
						computed[offsPrecomp + 2] = dxdxi;						// store Jacobian for this point sequentially in compute[][]
					}
					
					// calculate full inverse if necessary
					int offsPrecomp1 = offsPrecomp + 1;
					if (computed[offsPrecomp1] == null)
						computed[offsPrecomp1] = dxidx = Matrix.invert3x3(dxdxi);
					// or calculate full inverse (2nd step of invert3x3dualCall()) if necessary
					else if (computed[offsPrecomp1][10] == 0)
						computed[offsPrecomp1] = dxidx = Matrix.invert3x3dualCall(dxdxi, computed[offsPrecomp1]);
					det = computed[offsPrecomp1][9];
				}
				if (nNodes == 4) recalcJ = false;
			}
			
			// convert shape function derivatives to derivatives with respect to global coordinates
//			for (int a = 0; a < nNodes; a++) {
//				int aN = a * NCOORD;
//				for (int i = 0; i < NCOORD; i++) {
//					int aNi = aN + i;
//					dNdx[aNi] = 0;
//					for (int j = 0, aNj = aN, jN = i; j < NCOORD; j++, aNj++, jN += NCOORD)
//						dNdx[aNi] += dNdxi[aNj] * dxidx[jN];
//				}
//			}
			for (int a = 0; a < nNodes; a++) {								// optimised, unrolled loop for 3D
				int aN = a * NCOORD;
				double v1 = dNdxi[aN], v2 = dNdxi[aN + 1], v3 = dNdxi[aN + 2];
				dNdx[aN++] = v1 * dxidx[0] + v2 * dxidx[NCOORD] + v3 * dxidx[NCOORD + NCOORD];
				dNdx[aN++] = v1 * dxidx[1] + v2 * dxidx[1+NCOORD] + v3 * dxidx[1+NCOORD + NCOORD];
				dNdx[aN++] = v1 * dxidx[2] + v2 * dxidx[2+NCOORD] + v3 * dxidx[2+NCOORD + NCOORD];
			}
			
			// the computation of infinitisemal strain from differentiated displacements is only necessary
			// for nonlinear materials, as linear elasticity's stiffness is independent of strain
			if (computeStrain) {
				for (int i = 0, s = 0; i < NCOORD; i++) {
					int iN = i * nNodes;
					for (int j = 0; j < NCOORD; j++) {
						strain[s] = 0;
						for (int a = 0, iNa = iN, aNj = j, jNa = j * nNodes, aNi = i; a < nNodes; a++, iNa++, aNj += nNodes, jNa++, aNi += nNodes)
							strain[s] += .5 * (displacement[iNa] * dNdx[aNj] + displacement[jNa] * dNdx[aNi]);
						s++;
					}
				}
			}
			
			// compute material tangent stiffness dSigma/dEpsilon tensor, which is equal to linear elasticity's C[i,j,k,l]
			//double[] dsde = matStiffness1D(matProps);			
			double w_intPt_det = w[intPt] * det;

			// compute element's stiffness matrix
			for (int a = 0; a < nNodes; a++) {
				int aNC = a * NCOORD;
				for (int i = 0; i < NDOF; i++) {
					int rowCol = (a * NDOF + i) * nDOFn, iNNN = i*NCOORD*NDOF*NCOORD;
					for (int b = 0; i < nNodes; b++) {
						int bNC = b * NCOORD;
						for (int k = 0; k < NDOF; k++) {
							int iNNNkN = iNNN + k*NCOORD;
							for (int j = 0, bNC2 = bNC, aNC2 = aNC, jNN = 0; j < NCOORD; j++, aNC2++, jNN += NDOF*NCOORD) {
								// ijkl composed to access a 1D-array emulating a 4-D array C[i,j,k,l]: ((i * NCOORD + j) * NDOF + k) * NCOORD + l
								int ijkl = iNNNkN + jNN;
								// last loop, l over 0:NCOORD is unrolled
								K[rowCol++] += dsde[ijkl++] * dNdx[bNC2++] * dNdx[aNC2] * w_intPt_det;				
								K[rowCol++] += dsde[ijkl++] * dNdx[bNC2++] * dNdx[aNC2] * w_intPt_det;				
								K[rowCol++] += dsde[ijkl] * dNdx[bNC2++] * dNdx[aNC2] * w_intPt_det;				
							}							
						}
					}
				}
			}
		}
		return K;
	}
	
	
	double[] elementDistLoadV(int nFaceN, int elemFlags, double[] coords, double t1, double t2, double t3) {
		int nPoints = numIntegPoints2D[nFaceN], nCoordM1 = NCOORD - 1;
		double[] dxdxi = new double[NCOORD * (NCOORD-1)], r = new double[NDOF * nFaceN];
		
		double[] xiList = integPoints_2D[nFaceN][nPoints];
		double[] w = integWeights_2D[nFaceN][nPoints];
		
		for (int intPt = 0, iP = 0; intPt < nPoints; intPt++) {
			double xI0 = xiList[iP++], xI1 = xiList[iP++];
			double[] N = shapeFunctions_2D(nFaceN, xI0, xI1);
			double[] dNdxi = shapeFunctionDerivs_2D(nFaceN, xI0, xI1);

			// Jacobian & determinant calculation
			for (int i = 0, aEnd = NCOORD * nFaceN; i < NCOORD; i++) {
				int iC = i * NCOORD;
				for (int j = 0; j < nCoordM1; j++) {
					int iCj = iC + j;
					dxdxi[iCj] = 0;									// clear previous point's calculations
					for (int a = 0; a < aEnd; a += NCOORD)
						dxdxi[iCj] += coords[a + i] * dNdxi[a + j];
				}
			}
			double dlt3 = dxdxi[0]*dxdxi[3]-dxdxi[1]*dxdxi[2], dlt2 = dxdxi[0]*dxdxi[5]-dxdxi[1]*dxdxi[4], dlt1 = -dxdxi[3]*dxdxi[4]+dxdxi[2]*dxdxi[5];
			double det = Math.sqrt(dlt1*dlt1 + dlt2*dlt1 + dlt3*dlt3), w_det = w[intPt] * det;
			
			for (int a = 0, row = 0; a < nFaceN; a++) {
				double vN_w_det = N[a] * w_det;
				r[row++] += vN_w_det * t1; r[row++] += vN_w_det * t2; r[row++] += vN_w_det * t3; 
			}
		}
		return r;
	}
	
	
	double[] globalMass(int totNodes, double[] coords, int nElements, int[] elemFlags, boolean lumpedMass,
						int[] nElemNodes, int[][] connect, double[] matProps) {
		int totDOFn = NDOF * totNodes;
		double[] M = new double[totDOFn * totDOFn];
		double[] elemCoord = new double[NCOORD * MAXNODES];
		// prepare empty array space for precomputed arrays of differing integration node counts
		double[][][] computedN = new double[4][][];
		computedN[0] = new double[OFFSPRECOMP + 4 * 3][]; computedN[1] = new double[OFFSPRECOMP + 8 * 3][];
		computedN[2] = new double[OFFSPRECOMP + 10 * 3][]; computedN[3] = new double[OFFSPRECOMP + 20 * 3][];
		
		for (int e = 0; e < nElements; e++) {
			
			int nNodes = nElemNodes[e];
			for (int a = 0, c = 0; a < nNodes; a++) {
				int conn = connect[e][a] * NCOORD;
				elemCoord[c++] = coords[conn++]; elemCoord[c++] = coords[conn++]; elemCoord[c++] = coords[conn++];
			}
			
			double[][] computed = null;
			switch(nNodes) {								// we might have precomputed data depending on node count of current element
			case 4: computed = computedN[0]; break;
			case 8: computed = computedN[1]; break;
			case 10: computed = computedN[2]; break;
			case 20: computed = computedN[3]; break;
			}
			double[] Mme = elementMassMatrix(nNodes, elemCoord, matProps, elemFlags[e], lumpedMass, computed);
			
			// add element's stiffness into global stiffness matrix (both simulated in 1D-arrays)
			for (int a = 0, nDOFn = NDOF * nNodes; a < nNodes; a++) {
				int row = connect[e][a] * totDOFn, rowMme = a * nDOFn;
				for (int i = 0; i < NDOF; i++) {
					int row2 = row + i * totDOFn, rowMme2 = rowMme + i * nDOFn;
					for (int b = 0; b < nNodes; b++) {
						int offs = row2 + connect[e][b] * NDOF, offsMme = rowMme2 + b * NDOF;
						M[offs++] += Mme[offsMme++];
						M[offs++] += Mme[offsMme++];
						M[offs++] += Mme[offsMme++];
					}
				}
			}
		}
		return M;
	}
	
	
	double[] globalStiffness(	int totNodes, double[] coords, int nElements, int[] elemFlags,
								int[] nElemNodes, int[][] connect, double[] matProps, double[] nDOFs) {
		double[] S = new double[NDOF * totNodes * NDOF * totNodes];
		double[] elemCoord = new double[NCOORD * MAXNODES];
		double[] elemDOF = new double[NDOF * MAXNODES];
		double[][][] computedA = new double[4][OFFSPRECOMP + MAXNODES * 3][];		// computedA will be filled up for different node counts
		double[] dsde = null;
		
		for (int e = 0; e < nElements; e++) {
			// find out node coordinates and DOF (only 3D supported) of current element
			int nNodes = nElemNodes[e];
			for (int n = 0, c = 0, d = 0; n < nNodes; n++) {
				int conn = connect[e][n] * NCOORD, connDOF = connect[e][n] * NDOF;
				elemCoord[c++] = coords[conn++]; elemCoord[c++] = coords[conn++]; elemCoord[c++] = coords[conn++];			
				elemDOF[d++] = nDOFs[connDOF + 1]; elemDOF[d++] = nDOFs[connDOF + 2]; elemDOF[d++] = nDOFs[connDOF + 3];
			}
			
			if (dsde == null) dsde = matStiffness_1D(matProps);
			double[] K = elementStiffnessMatrix(nNodes, coords, dsde, elemFlags[e], elemDOF, true, computedA[computedArrayFromNodes[nNodes]]);
			int totDOFn = NDOF * totNodes, totDOF2n = totDOFn * NDOF;
			int nDOFn = NDOF * nNodes, nDOF2n = nDOFn * NDOF;
			
			// add element stiffness matrix to global stiffness matrix
			for (int na = 0; na < nNodes; na++) {
				int offsRow = totDOF2n * connect[e][na], offsRowK = nDOF2n * na;
				for (int i = 0; i < NDOF; i++, offsRow += totDOFn, offsRowK += nDOFn) {
					for (int nb = 0; nb < nNodes; nb++) {
						int offsK = offsRow + NDOF * connect[e][nb];
						// add three degrees of freedom in unrolled loop
						S[offsK++] += K[offsRowK++]; S[offsK++] += K[offsRowK++]; S[offsK++] += K[offsRowK++];
					}
				}
			}
		}
		return S;
	}
	
	
	double[] globalTraction(int nNodes, int nDloads, double[] coords, int[] nElemNodes, int[] elemFlags, int[][] connect, int[][] dLoads) {

		double[] r = new double[NDOF * nNodes];			// r is a vector
		
		for (int load = 0; load < nDloads; load++) {
			// extract coordinates of appropriate element face
			int[] dLoadF = dLoads[load];
			int elem = dLoadF[0], face = dLoadF[1];
			int nN = nElemNodes[elem];
			int nFNodes = nFaceNodes[nN];
			int[] nodeList = faceNodes[nN][face];
			double[] elemCoords = new double[NCOORD * nFNodes];
			int[] connectElem = connect[dLoadF[0]];
			
			for (int nF = 0, c = 0; nF < nFNodes; nF++) {
				int offsC = NCOORD * connectElem[nodeList[nF]];
				// transfer a 3D-coordinate to a localised array from indexes inside dLoads array
				elemCoords[c++] = coords[offsC++]; elemCoords[c++] = coords[offsC++]; elemCoords[c++] = coords[offsC++];
			}
			
			double[] rElem = elementDistLoadV(nFNodes, elemFlags[elem], elemCoords, dLoadF[2], dLoadF[3], dLoadF[4]);
			
			// assemble element load vector into global load vector
			for (int nF = 0, row_rEL = 0; nF < nFNodes; nF++) {
				int row = connectElem[nodeList[nF]] * NDOF;
				// unrolled loop over degrees of freedom (NDOF)
				r[row++] += rElem[row_rEL++]; r[row++] += rElem[row_rEL++]; r[row++] += rElem[row_rEL++];
			}
		}
		return r;
		
	}
	
	
	
	// method returns the coordinates of integration points in 1D-array, given the node count & requested integration point count of element
	// note: this method is redundant through usage of integPoints[][][] array
	public static double[] integrationPoints(int nNodes, int nPoints) {
		switch (nNodes) {
		case 4: case 10:
			switch (nPoints) {
			case 1: double[] iPts1 = {.25,.25,.25}; return iPts1;
			case 4: double[] iPts4 = {0.58541,0.1382,0.1382, 0.1382,0.58541,0.1382, 0.1382,0.1382,0.58541, 0.1382,0.1382,0.1382}; return iPts4; 
			case 5: double[] iPts5 = {.25,.25,.25, .5,1./6.,1./6., 1./6.,.5,1./6., 1./6.,1./6.,.5, 1./6.,1./6.,1./6.}; return iPts5;
			default: throw new InvalidParameterException("FEM1.integrationPoints(): Invalid point count.");
			}
		case 8: case 20:
			switch (nPoints) {
			case 1: double[] iPts1 = {0,0,0}; return iPts1;
			case 8: double[] iPts8 = new double[24], x1D = {-0.57735,0.57735};
				for (int n=0,k=0;k<2;k++) for (int j=0;j<2;j++) for (int i=0;i<2;i++) {
					iPts8[n++] = x1D[i]; iPts8[n++] = x1D[j]; iPts8[n++] = x1D[k];
				}
				return iPts8;
			case 27: double[] iPts27 = new double[81], x1D2 = {-0.7746,0,0.7746};
				for (int n=0,k=0;k<3;k++) for (int j=0;j<3;j++) for (int i=0;i<3;i++) {
					iPts27[n++] = x1D2[i]; iPts27[n++] = x1D2[j]; iPts27[n++] = x1D2[k];
				}
				return iPts27;
			default: throw new InvalidParameterException("FEM1.integrationPoints(): Invalid point count.");
			}
		default: throw new InvalidParameterException("FEM1.integrationPoints(): Invalid node count.");
		}
	}
	
	// method returns the coordinates of mass integration points, given the node count & requested mass integration point count of element
	// note: this method is redundant through usage of integWeights[][][] array
	public static double[] integrationWeights(int nNodes, int nPoints) {
		switch (nNodes) {
		case 4: case 10:
			switch (nPoints) {
			case 1: double[] w1 = { 1.0/6.0 }; return w1;
			case 4: double[] w4 = { 1./24.,1./24.,1./24.,1./24. }; return w4;
			default: throw new InvalidParameterException("FEM1.integrationWeights(): Invalid point count.");
			}
		case 8: case 20:
			switch (nPoints) {
			case 1: double[] w1 = { 8 }; return w1;
			case 8: double[] w8 = { 1,1,1,1,1,1,1,1 }; return w8;
			case 27: double[] w27 = new double[nPoints], w1D = {.55556,.88889,.55556};
				for (int n=0,k=0;k<3;k++) for (int j=0;j<3;j++) for (int i=0;i<3;i++) {
					w27[n++] = w1D[i] * w1D[j] * w1D[k];
				}
				return w27;
			default: throw new InvalidParameterException("FEM1.integrationWeights(): Invalid point count.");
			}
		default: throw new InvalidParameterException("FEM1.integrationWeights(): Invalid node count.");
		}
	}
	
	
	// generates shear modulus & Poisson ratio into an elasticity tensor
	public static double[][][][] matStiffness(double matProps[]) {
		double mu = matProps[0], nu = matProps[1];
		double[][][][] tensor = new double[NDOF][NCOORD][NDOF][NCOORD];
		double v = 2 * mu * nu / (1 - 2 * nu);
		
		for (int i = 0; i < NDOF; i++)
			for (int j = 0; j < NCOORD; j++) {
				double[][] tensor2 = tensor[i][j];
				for (int k = 0; k < NDOF; k++) {
					double[] tensor3 = tensor2[k];
					for (int l = 0; l < NCOORD; l++) {
						if (i == j && k == l) 								tensor3[l] += v;
						else if ((i == k && j == l) || (i == l && j == k)) 	tensor3[l] += mu;
					}
				}
			}
		return tensor;
	}
	
	// generates shear modulus & Poisson ratio into an elasticity tensor within a 1D-array
	public static double[] matStiffness_1D(double[] matProps) {
		
		int jBlock = NDOF * NCOORD, iBlock = NCOORD * jBlock, tBlock = NDOF * iBlock;
		double[] tensor = new double[tBlock];
		
		double mu = matProps[0], nu = matProps[1];
		double v = 2.0 * mu * nu / (1.0 - 2.0 * nu);
		
		for (int t = 0, i = 0; i < NDOF; i++)
			for (int j = 0; j < NCOORD; j++) {
				for (int k = 0; k < NDOF; k++) {
					for (int l = 0; l < NCOORD; l++) {
						if (i == j && k == l) 								tensor[t++] += v;
						else if ((i == k && j == l) || (i == l && j == k)) 	tensor[t++] += mu;
						else												t++;
					}
				}
			}
		return tensor;
	}

	
	// computes stress sigma matrix from a strain epsilon matrix
	public static double[][] matStress(double strain[][], double matProps[]) {
		double[][][][] tensor = matStiffness(matProps);
		double[][] stress = new double[NDOF][NCOORD];
		
		for (int i = 0; i < NDOF; i++)
			for (int j = 0; j < NCOORD; j++)
				for (int k = 0; k < NDOF; k++)
					for (int l = 0; l < NCOORD; l++)
						stress[i][j] += tensor[i][j][k][l] * strain[k][l];
		return stress;
	}

	// computes a 1D array stress sigma matrix from a strain epsilon 1D-array matrix
	public static double[] matStress_1D(double strain[], double matProps[]) {
		double[] tensor = matStiffness_1D(matProps);
		int ndnc = NDOF * NCOORD;
		double[] stress = new double[ndnc];
		
		for (int t = 0, i = 0; i < ndnc; i += NCOORD)
			for (int j = i; j < NCOORD; j++) {
				double s = 0;
				for (int k = 0; k < ndnc; k += NCOORD)
					for (int l = k; l < NCOORD; l++)		s += tensor[t++] * strain[l];
				stress[j] = s;
			}
		return stress;
	}
	
	
	double[] shapeFunctions(int nNodes, double xI0, double xI1, double xI2) {
		switch (nNodes) {
		case 4: double[] N4 = {xI0, xI1, xI2, 1.-xI0-xI1-xI2}; return N4;
		case 10: double xI4 = 1.-xI0-xI1-xI2, xI04 = xI0*4, xI14 = xI1*4; 
			double[] N10 = {(xI0*2-1)*xI0, (xI1*2-1)*xI1, (xI2*2-1)*xI2, (xI4*2-1)*xI4, xI04*xI1, xI14*xI2, xI2*xI04, xI04*xI4, xI14*xI4, xI2*4*xI4};
			return N10;
		case 8: { double xI0p1 = xI0+1, xI1p1 = xI1+1, xI2p1d8 = (xI2+1)*DIV8, xI0m1 = 1.-xI0, xI1m1 = 1.-xI1, xI2m1d8 = (1.-xI2)*DIV8;
			double[] N8 = {	xI0m1*xI1m1*xI2m1d8, xI0p1*xI1m1*xI2m1d8, xI0p1*xI1p1*xI2m1d8, xI0m1*xI1p1*xI2m1d8,
							xI0m1*xI1m1*xI2p1d8, xI0p1*xI1m1*xI2p1d8, xI0p1*xI1p1*xI2p1d8, xI0m1*xI1p1*xI2p1d8};
			return N8; }
		case 20: { double xI0p1 = xI0+1, xI1p1 = xI1+1, xI2p1d8 = (xI2+1)*DIV8, xI0m1 = 1.-xI0, xI1m1 = 1.-xI1, xI2m1d8 = (1.-xI2)*DIV8;
			double xI2p1d4 = (xI2+1)*DIV4, xI2m1d4 = (1.-xI2)*DIV4, xI0m1p2 = (1.-xI0*xI0), xI1m1p2 = (1.-xI1*xI1), xI2m1p2d4 = (1.-xI2*xI2)*DIV4;
			double[] N20 = {	xI0m1*xI1m1*xI2m1d8 * (-xI0-xI1-xI2-2), xI0p1*xI1m1*xI2m1d8 * (xI0-xI1-xI2-2),
								xI0p1*xI1p1*xI2m1d8 * (xI0+xI1-xI2-2), xI0m1*xI1p1*xI2m1d8 * (-xI0+xI1-xI2-2),
								xI0m1*xI1m1*xI2p1d8 * (-xI0-xI1+xI2-2), xI0p1*xI1m1*xI2p1d8 * (xI0-xI1+xI2-2),
								xI0p1*xI1p1*xI2p1d8 * (xI0+xI1+xI2-2), xI0m1*xI1p1*xI2p1d8 * (-xI0+xI1+xI2-2),
								xI0m1p2*xI1m1*xI2m1d4, xI0p1*xI1m1p2*xI2m1d4, xI0m1p2*xI1p1*xI2m1d4, xI0m1*xI1m1p2*xI2m1d4,
								xI0m1p2*xI1m1*xI2p1d4, xI0p1*xI1m1p2*xI2p1d4, xI0m1p2*xI1p1*xI2p1d4, xI0m1*xI1m1p2*xI2p1d4,
								xI0m1*xI1m1*xI2m1p2d4, xI0p1*xI1m1*xI2m1p2d4, xI0p1*xI1p1*xI2m1p2d4, xI0m1*xI1p1*xI2m1p2d4};
			return N20; }
		default: throw new InvalidParameterException("FEM1.shapeFunctions(): Invalid node count.");
		}
	}
	
	double[] shapeFunctions_2D(int nFaceN, double xI0, double xI1) {
		switch (nFaceN) {
		case 3: double[] N3 = {xI0, xI1, 1.-xI0-xI1}; return N3;
		case 4: { double xI0m1 = (1.-xI0)*.5, xI0p1 = (1.+xI0)*.5, xI1m1 = (1.-xI1)*.5, xI1p1 = (1.+xI1)*.5;	// .5 * .5 = .25
			double[] N4 = {xI0m1*xI1m1, xI0p1*xI1m1, xI0p1*xI1p1, xI0m1*xI1p1}; return N4; }
		case 6: double xI2 = 1.-xI0-xI1, xI04 = xI0*4;
			double[] N6 = {(2*xI0-1)*xI0, (2*xI1-1)*xI1, (2*xI2-1)*xI2, xI04*xI1, 4*xI1*xI2, xI2*xI04}; return N6;
		case 8: { double xI0m1 = (1.-xI0), xI0p1 = (1.+xI0), xI1m1 = (1.-xI1), xI1p1 = (1.+xI1), xI0sqm1 = 1.-xI0*xI0, xI1sqm1 = 1.-xI1*xI1;
			double[] N8 = {	-.25*xI0m1*xI1m1*(xI0p1+xI1), .25*xI0p1*xI1m1*(xI0m1-xI1), .25*xI0p1*xI1p1*(xI0+xI1m1), .25*xI0m1*xI1p1*(xI1m1-xI0),
							.5*xI0sqm1, .5*xI0p1*xI1sqm1, .5*xI0sqm1*xI1p1, .5*xI0m1*xI1sqm1}; return N8;
		
		}
		default: throw new InvalidParameterException("FEM1.shapeFunctions_2D(): Invalid face count.");
		}
	}

	
	double[] shapeFunctionDerivs(int nNodes, double xI0, double xI1, double xI2) {
		switch (nNodes) {
		case 4: double[] dNdxi4 = {1,0,0,0,1,0,0,0,1,-1,-1,-1}; return dNdxi4;
		case 10: double xI4 = 1. - xI0 - xI1 - xI2, xI0t4 = xI0*4, xI1t4 = xI1*4, xI2t4 = xI2*4, xI4t4 = xI4*4, xI4t4m1 = -(xI4t4 - 1);
			double[] dNdxi10 = {xI0t4-1, 0, 0, 0, xI1t4-1, 0, 0, 0,xI2t4-1, xI4t4m1, xI4t4m1, xI4t4m1, xI1t4, xI0t4, 0, 0, xI2t4, xI1t4, xI2t4, 0, xI0t4,
								xI4t4-xI0t4, -xI0t4, -xI0t4, -xI1t4, xI4t4-xI1t4, -xI1t4, -xI2*xI4t4, -xI2t4, xI4t4-xI2t4};
			return dNdxi10;
		case 8: { double xI0p1 = xI0+1, xI1p1 = xI1+1, xI2p1 = xI2+1, xI0m1 = 1.-xI0, xI1m1 = 1.-xI1, xI2m1 = 1.-xI2;
			double xI1m1d8 = xI1m1*DIV8, xI2m1d8 = xI2m1*DIV8, xI1p1d8 = xI1p1*DIV8, xI2p1d8 = xI2p1*DIV8;
			double[] dNdxi8 = { -xI1m1*xI2m1d8, -xI0m1*xI2m1d8, -xI0m1*xI1m1d8, xI1m1*xI2m1d8, -xI0p1*xI2m1d8, -xI0p1*xI1m1d8,
								xI1p1*xI2m1d8, xI0p1*xI2m1d8, -xI0p1*xI1p1d8, -xI1p1*xI2m1d8, xI0m1*xI2m1d8, -xI0m1*xI1p1d8,
								-xI1m1*xI2p1d8, -xI0m1*xI2p1d8, xI0m1*xI1m1d8, xI1m1*xI2p1d8, -xI0p1*xI2p1d8, xI0p1*xI1m1d8,
								xI1p1*xI2p1d8, xI0p1*xI2p1d8, xI0p1*xI1p1d8, -xI1p1*xI2p1d8, xI0m1*xI2p1d8, xI0m1*xI1p1d8};
			return dNdxi8; }
		case 20: { double xI0p1 = xI0+1, xI1p1 = xI1+1, xI2p1 = xI2+1, xI0m1 = 1.-xI0, xI1m1 = 1.-xI1, xI2m1 = 1.-xI2;
			double[] dNdxi20 = new double[nNodes * 3];
			double xI012mmm2 = (-xI0-xI1-xI2-2)*DIV8, xI2m1d8 = xI2m1*DIV8, xI012mmm1 = xI0m1*xI1m1*xI2m1d8;
			dNdxi20[0] = -xI1m1*xI2m1*xI012mmm2 - xI012mmm1; dNdxi20[1] = -xI0m1*xI2m1*xI012mmm2 - xI012mmm1; dNdxi20[2] = -xI0m1*xI1m1*xI012mmm2 - xI012mmm1;
			double xI012pmm2 = (xI0-xI1-xI2-2)*DIV8, xI012pmm1 = xI0p1*xI1m1*xI2m1d8;
			dNdxi20[3] = xI1m1*xI2m1*xI012pmm2 + xI012pmm1; dNdxi20[4] = -xI0p1*xI2m1*xI012pmm2 - xI012pmm1; dNdxi20[5] = -xI0p1*xI1m1*xI012pmm2 - xI012pmm1;
			double xI012ppm2 = (xI0+xI1-xI2-2)*DIV8, xI012ppm1 = xI0p1*xI1p1*xI2m1d8;
			dNdxi20[6] = xI1p1*xI2m1*xI012ppm2 + xI012ppm1; dNdxi20[7] = xI0p1*xI2m1*xI012ppm2 + xI012ppm1; dNdxi20[8] = -xI0p1*xI1p1*xI012ppm2 - xI012ppm1;
			double xI012mpm2 = (-xI0+xI1-xI2-2)*DIV8, xI012mpm1 = xI0m1*xI1p1*xI2m1d8;
			dNdxi20[9] = -xI1p1*xI2m1*xI012mpm2 - xI012mpm1; dNdxi20[10] = xI0m1*xI2m1*xI012mpm2 + xI012mpm1; dNdxi20[11] = -xI0m1*xI1p1*xI012mpm2 - xI012mpm1;
			double xI012mmp2 = (-xI0-xI1+xI2-2)*DIV8, xI2p1d8 = xI2p1*DIV8, xI012mmp1 = xI0m1*xI1m1*xI2p1d8;
			dNdxi20[12] = -xI1m1*xI2p1*xI012mmp2 - xI012mmp1; dNdxi20[13] = -xI0m1*xI2p1*xI012mmp2 - xI012mmp1; dNdxi20[14] = xI0m1*xI1m1*xI012mmp2 + xI012mmp1;
			double xI012pmp2 = (xI0-xI1+xI2-2)*DIV8, xI012pmp1 = xI0p1*xI1m1*xI2p1d8;
			dNdxi20[15] = xI1m1*xI2p1*xI012pmp2 + xI012pmp1; dNdxi20[16] = -xI0p1*xI2p1*xI012pmp2 - xI012pmp1; dNdxi20[17] = xI0p1*xI1m1*xI012pmp2 + xI012pmp1;
			double xI012ppp2 = (xI0+xI1+xI2-2)*DIV8, xI012ppp1 = xI0p1*xI1p1*xI2p1d8;
			dNdxi20[18] = xI1p1*xI2p1*xI012ppp2 + xI012ppp1; dNdxi20[19] = xI0p1*xI2p1*xI012ppp2 + xI012ppp1; dNdxi20[20] = xI0p1*xI1p1*xI012ppp2 + xI012ppp1;
			double xI012mpp2 = (-xI0+xI1+xI2-2)*DIV8, xI012mpp1 = xI0m1*xI1p1*xI2p1d8;
			dNdxi20[21] = -xI1p1*xI2p1*xI012mpp2 - xI012mpp1; dNdxi20[22] = xI0m1*xI2p1*xI012mpp2 + xI012mpp1; dNdxi20[23] = xI0m1*xI1p1*xI012mpp2 + xI012mpp1;
			double xI2m1d4 = xI2m1*DIV4, xI2p1d4 = xI2p1*DIV4, xI0t2 = -xI0*2, xI1t2 = -xI1*2;
			double xI0p2m1d4 = (1-xI0*xI0)*DIV4; dNdxi20[24] = xI0t2*xI1m1*xI2m1d4; dNdxi20[25] = -xI0p2m1d4*xI2m1; dNdxi20[26] = -xI0p2m1d4*xI1m1;
			double xI1p2m1d4 = (1-xI1*xI1)*DIV4; dNdxi20[27] = xI1p2m1d4*xI2m1; dNdxi20[28] = xI1t2*xI0p1*xI2m1d4; dNdxi20[29] = -xI1p2m1d4*xI0p1;
			dNdxi20[30] = xI0t2*xI1p1*xI2m1d4; dNdxi20[31] = xI0p2m1d4*xI2m1; dNdxi20[32] = -xI0p2m1d4*xI1p1;
			dNdxi20[33] = -xI1p2m1d4*xI2m1; dNdxi20[34] = xI1t2*xI0m1*xI2m1d4; dNdxi20[35] = -xI1p2m1d4*xI0m1;
			dNdxi20[36] = xI0t2*xI1m1*xI2p1d4; dNdxi20[37] = -xI0p2m1d4*xI2p1; dNdxi20[38] = xI0p2m1d4*xI1m1;
			dNdxi20[39] = xI1p2m1d4*xI2p1; dNdxi20[40] = xI1t2*xI0p1*xI2p1d4; dNdxi20[41] = xI1p2m1d4*xI0p1;
			dNdxi20[42] = xI0t2*xI1p1*xI2p1d4; dNdxi20[43] = xI0p2m1d4*xI2p1; dNdxi20[44] = xI0p2m1d4*xI1p1;
			dNdxi20[45] = -xI1p2m1d4*xI2p1; dNdxi20[46] = xI1t2*xI0m1*xI2p1d4; dNdxi20[47] = xI1p2m1d4*xI0m1; 
			double xI2p2m1d4 = (1-xI2*xI2)*DIV4; dNdxi20[48] = -xI2p2m1d4*xI1m1; dNdxi20[49] = -xI2p2m1d4*xI0m1; dNdxi20[50] = -xI2*xI0m1*xI1m1*.5;
			dNdxi20[51] = xI1m1*xI2p2m1d4; dNdxi20[52] = -xI0p1*xI2p2m1d4; dNdxi20[53] = -xI2*xI0p1*xI1m1*.5;
			dNdxi20[54] = xI1p1*xI2p2m1d4; dNdxi20[55] = xI0p1*xI2p2m1d4; dNdxi20[56] = -xI2*xI0p1*xI1p1*.5;
			dNdxi20[57] = -xI1p1*xI2p2m1d4; dNdxi20[58] = xI0m1*xI2p2m1d4; dNdxi20[59] = -xI2*xI0m1*xI1p1*.5;		
			return dNdxi20; }
		default: throw new InvalidParameterException("FEM1.shapeFunctionDerivs(): Invalid node count.");
		}
	}

	double[] shapeFunctionDerivs_2D(int nNodes, double xI0, double xI1) {
		switch(nNodes) {
		case 3: double[] dNdxi3 = {1,0, 0,1, -1,-1}; return dNdxi3;
		case 4: { double xI0p1 = .25*(xI0+1), xI1p1 = .25*(xI1+1), xI0m1 = .25*(1.-xI0), xI1m1 = .25*(1.-xI1);
			double[] dNdxi4 = {-xI1m1,-xI0m1, xI1m1,-xI0p1, xI1p1,xI0p1, -xI1p1,xI0m1};
			return dNdxi4; }
		case 6: double xI2 = 1.-xI0-xI1, xI0t4 = xI0*4, xI1t4 = xI1*4, xI2t4 = xI2*4, xI2t4m1 = xI2t4-1;
			double[] dNdxi6 = {xI0t4-1,0, 0,xI1t4-1, xI2t4m1,xI2t4m1, xI1t4,xI0t4, -xI1t4,-xI0t4, xI2t4-xI0t4,xI2t4-xI1t4};
			return dNdxi6;
		case 8: { double xI0p1 = .25*(xI0+1), xI1p1 = .25*(xI1+1), xI0m1 = .25*(1.-xI0), xI1m1 = .25*(1.-xI1), xI0t2 = xI0*2, xI1t2 = xI1*2;
			double xI0msqm1 = .5*(1.-xI0*xI0), xI1sqm1 = .5*(1.-xI1*xI1);
			double[] dNdxi8 = {	xI1m1*(xI0t2+xI1),xI0m1*(xI1t2+xI0), xI1m1*(xI0t2-xI1),xI0p1*(xI1t2-xI0),
								xI1p1*(xI0t2+xI1),xI0p1*(xI1t2+xI0), xI1p1*(xI0t2-xI1),xI0m1*(xI1t2-xI0),
								-xI0*(1.-xI1),-xI0msqm1, xI1sqm1,-(1.+xI0)*xI1, -xI0*(1.+xI1),xI0msqm1, -xI1sqm1,-(1.-xI0)*xI1};
			return dNdxi8;}
		default: throw new InvalidParameterException("FEM1.shapeFunctionDerivs_2D(): Invalid node count.");
		}
	}
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			HELPER / INLINE METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	public void clearNodeChecks() { nodeChecks = new long[nodeBitSets]; }
	public boolean nodeChecked(int node) { return (nodeChecks[node >> 6] & (0x1L << (node & 63))) != 0; }
	public void checkNode(int node) { nodeChecks[node >> 6] |= (0x1L << (node & 63)); }
	public boolean externalNode(int n) { return n * 3 > node.length ? (nodeWorkFlag[n - nodes] & 1) == 1 : (nodeFlag[n] & 1) == 1; }
	
	int catchSmoothingGroup(String s) {
		if (s.startsWith("s")) {							// if caught a smoothing group separator, register start of patch
			int sGroup = Integer.valueOf(s.split(" ")[1]) - 1;
			return sGroup;
		}
		return -1;
	}

	static int[] decodePolygon(String s, boolean gotNormals) {
		int[] polygon;
		String[] dataRow;
		dataRow = s.split(" ");							// OBJ file with normals has a different format for polygon indexes
		int len = dataRow.length - 1;
		polygon = new int[len];
		for (int i = 1; i <= len; i++) {				// decode each polygon into it's node integer indexes
			if (gotNormals) dataRow[i] = dataRow[i].split("//")[0];
			polygon[i - 1] = Integer.valueOf(dataRow[i]) - 1;
		}
		return polygon;
	}

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int e = 0; e < elements2; e++) {
			FEM1Element elem = getElement2(e);
			if (elem == null) { sb.append("[" + e + " deleted]\n\n"); continue; }
			if (e == elements2) sb.append("<<<<<<<< element2Work[] >>>>>>>>\n");
			switch (elem.flags & 15) {
			case FEM1Element.TETRAHEDRON:
				sb.append("Tetrahedron T" + e);
				sb.append(", refs: n" + elem.nodeRef[0] + ",n" + elem.nodeRef[1] + ",n" + elem.nodeRef[2] + ",n" + elem.nodeRef[3]);
				if (elem.volume > 0) sb.append(" V");
				if (elem.data != null) {
					sb.append(" A:");
					sb.append((elem.data[6]>0?"A":".")+(elem.data[7]>0?"A":".")+(elem.data[8]>0?"A":".")+(elem.data[9]>0?"A":"."));
					sb.append(" e:");
					sb.append(	(elem.data[0]>0?"e":".") + (elem.data[1]>0?"e":".") + (elem.data[2]>0?"e":".")+
								(elem.data[3]>0?"e":".") + (elem.data[4]>0?"e":".") + (elem.data[5]>0?"e":"."));
				}
				sb.append("\n");

				if (elem.hasInterfaces()) {
					sb.append("Interfaces: [");
					for (int n = 0, neighbours2 = elem.neighbours * 2; n < neighbours2; n += 2)
						sb.append(((elem.neighbour[n + 1] & 15) != 0  ? "F" : "e") + elem.neighbour[n] + (n == neighbours2 - 2 ? "]\n" : ", "));
				}
				break;
			case FEM1Element.HEXAHEDRON: sb.append("Hexahedron(E" + e + ")"); break;
			case FEM1Element.ISOTETRAHEDRON: sb.append("Isotetrahedron(E" + e + ")"); break;
			case FEM1Element.ISOHEXAHEDRON: sb.append("Isohexahedron(E" + e + ")"); break;
			}
			
			sb.append("Neighbours: [");
			if (elem.neighbours == 0 || elem.neighbour == null) sb.append("n/a]\n");
			else {
				if (elem.neighboursF > 0) sb.append("F: ");
				for (int n = 0, nF = elem.neighboursF; n < elem.neighbours; n++, nF--) {
					if (elem.hasInterfaces()) {
						if (nF == 0 && elem.neighboursF < elem.neighbours) sb.append("e: ");
						sb.append("T" + elem.neighbour[n * 2] + (n == elem.neighbours - 1 ? "]\n" : ", "));
					} else {
						if (nF == 0 && elem.neighboursF < elem.neighbours) sb.append("e: ");
						sb.append("T" + elem.neighbour[n] + (n == elem.neighbours - 1 ? "]\n" : ", "));
					}
				}
			}
			sb.append("\n");
		}
		return sb.toString();
	}
	
	
	
	public void closestNodesToOBJ(FEM1Octree octree) {
		File file = new File(name + "_processed.obj");
		if (!file.exists()) {
			try {	file.createNewFile();
			} catch (IOException e) { e.printStackTrace(); }
		}
		BufferedWriter bw = null;
		try {		bw = new BufferedWriter(new FileWriter(file));
		} catch (IOException e) { e.printStackTrace(); }
		StringBuilder sb = new StringBuilder();
		sb.append("# Visualisation of result of closest nodes octree search algorithm\n#\n");
		String precFormat = "%." + precision_OBJ + "f"; //(precision < 1 ? 1 : precision > 10 ? 10 : precision) + "f";
		
		int[] closest = new int[nodes * 2];					// worst case allocation
		double[] distance = new double[nodes];				// worst case allocation
		int nPairs = closestNodePairs(octree, closest, distance, null);
		double smallestD = distance[0];						// find smallest distance of two closest nodes to adjust the width of output strips
		for (int nP = 1; nP < nPairs; nP++) if (distance[nP] < smallestD) smallestD = distance[nP];
		mergeNodeWork();									// avoid dealing with extraneous nodes in nodeWork[]
		
		for (int nP = 0; nP < nPairs; nP++) {				// output vertex definitions
			int n1_3 = closest[nP * 2] * 3, n2_3 = closest[nP * 2 + 1] * 3;
			double x1 = node[n1_3++], y1 = node[n1_3++], z1 = node[n1_3++]; 
			double x2 = node[n2_3++], y2 = node[n2_3++], z2 = node[n2_3++];
			// create a polygon strip between 2 closest nodes, whose width if 1/5th of the smallest distance found
			double scaleFactor = smallestD * (1. / 5.) / distance[nP];
			double ny = (x2 - x1) * scaleFactor, nx = -(y2 - y1) * scaleFactor, nz = (z2 - z1) * scaleFactor;
			sb.append("v  " + String.format(precFormat, x1-nx) + " " + String.format(precFormat, y1-nx) + " " + String.format(precFormat, z1-nz) + "\n");
			sb.append("v  " + String.format(precFormat, x1+nx) + " " + String.format(precFormat, y1+ny) + " " + String.format(precFormat, z1+nz) + "\n");
			sb.append("v  " + String.format(precFormat, x2+nx) + " " + String.format(precFormat, y2+ny) + " " + String.format(precFormat, z2+nz) + "\n");
			sb.append("v  " + String.format(precFormat, x2-nx) + " " + String.format(precFormat, y2-ny) + " " + String.format(precFormat, z2-nz) + "\n");
		}
		sb.append("# " + nPairs * 4 + " vertices\n\ng " + name + "_processed\n");

		for (int nP = 0, f = 1; nP < nPairs; nP++) {		// output polygon definitions (vertex count starts with 1!)
			sb.append("f " + f++ + " " + f++ + " " + f++ + " " + f++ + "\n");
		}
		sb.append("# " + nPairs + " faces\n\ng\n");
		
		try {												// write out & close file
			bw.write(sb.toString());
			bw.flush();
			bw.close();
		} catch (IOException e) { e.printStackTrace(); }
	}

}


class Processor implements Runnable {
	
	private int id;
	public Processor(int id) {
		this.id = id;
	}
	
	public void run() {
		System.out.println("Starting: " + id);
		
		try {
			Thread.sleep(2000);
		} catch (InterruptedException e) {}
		
		System.out.println("Completed: " + id);
	}
	
	
}
