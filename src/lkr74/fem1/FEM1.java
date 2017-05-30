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
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import javax.management.RuntimeErrorException;

import lkr74.mathgenerics.SpreadBin;
import lkr74.mathgenerics.VisitBitArray;
import lkr74.matrixlib.BinBitImage;
import lkr74.matrixlib.Matrix;
import lkr74.matrixlib.NSPMatrix;
import lkr74.matrixlib.NSPArray;
import lkr74.matrixlib.NSPNode;

public class FEM1 {
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			FINITE ELEMENT CONSTRUCTOR CLASS 1															//
	//			Leonard Krylov 2017																			//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////	

	// note: for more programmer-friendly access, methods for accessing the following 1D-arrays
	// can be implemented, to reference specific coordinates and elements
	String name;
	// node, element, freed element, polygon & smooth patch counts for this FEM1 object
	public int nodes=0, nodesWork=0, deletedNodes=0, elements=0, elementsWork=0;
	public int elements2=0, elements2Work=0, elements2Free=0, polygons=0, patches=0, maxNodesInPoly = 3;
	int precision_OBJ = 5;							// decimal places of the input OBJ file
	
	public double[] node = null;					// node coordinates come in (x,y,z)-triples within this 1D-array
	public double[] nodeWork = null;				// the work array for newly added nodes is also processed, and is reintegrated after reaching some length
	byte nodeFlag[]=null, nodeWorkFlag[]=null;		// up to 8 multipurpose flags for each node
	int[] deletedNode=null;							// stack of previously deleted nodes
	double[][] nodeCoordSort = null;				// arrays holding coordinates sorted by (x, y, z)
	int[][] nodeCoordSortI = null;					// reindexing tables of the sorted coordinates
	
	int[] element = null, elementWork = null;		// each element's node indexes come here (ex. tetrahedral indexes queued in a,b,c,d quartets)
	byte[] nodeCount = null;						// each element's node count is also specifying the type of element

	int[] polygon = null;							// polygons come in N-tuplets, N = no. of nodes (ex. faces come in a,b,c triplets)
	float[] polygonNormal = null;					// polygon normals come in triplets (for now assumed flat so only 3 nodes used for calculation)
	int[] polygonOffset = null;						// polygon node offsets work both as specifiers of poly node count and give position of polygon in array
	VisitBitArray polygonCheck = null;
	
	int[] edgeIndex = null;							// indexes matching each shared edge to particular polygons
	double[] edge = null;							// array specifying edge length of every polygonal edge
	int[] edgeNode = null;							// convenient back-reference from (edge index * 2) to the underlying nodes
	int edges;										// total nonredundant count of polygon edges

	int[] patch = null;								// the patch/smoothing group array holds offsets into face array for the start of every patch's faces
	
	FEM1Element[] element2 = null;					// object-oriented element types go into this structure
	FEM1Element[] element2Work = null;				// the work array for newly added elements is also processed, and is merged after reaching some length
	int[] element2Free = null;
	
	int[][] elementSupports = null;					// indexes every node's supported elements
	int elementSupportMaxArrayL = 9;				// maximum acquired length of nodal element support array for temp.array worst-case allocation
	int[][] elementNeighbours = null;				// indexes neighbours of every element (this array is constructed from elementSupports[])
	int elementNeighboursMaxArray = 0;				// holds largest neighbourhood for temp.array worst-case allocation

	int[][] polygonSupports = null;					// indexes every node's supported polygons
	int polygonSupportMaxArrayL = 0;				// maximum acquired length of nodal polygon support array
	int[][] polygonNeighbours = null;				// indexes neighbours of every polygon (this array is constructed from polygonSupports[])
	int[] borderPolygons = null;					// array optimises iteration over polygons bordering a patch by keeping an indexing over them
	int[][] patchBorders = null;					// every patch can have more than one edge border cycle, cycles lie sequentially in each patch's subarray

	int nodeBitSets = 0;
	int[] encapCycle = null;
	
	public double[] bBox = null;					// bounding box of dataset kept here

	// global parameters of the FEM solution
	double maxTsize = 1;							// maximal bound on size of any constructed tetrahedron
	public int nodeworkFactor = 8;					// let nodeWork[] array be maxium 1/8th of node[] array
	double avgEdgeLength = 0, minEdgeLength = 0, maxEdgeLength = 0;	// the average and min/max edge lengths found are stored here
	SpreadBin edgeSpread = null;
	static int DEBUG_LEVEL = 1;
	
	// the datatype requested from the object instancer
	public static final int MESH_HANDMADE = 0, MESH_HANDMADE_OBJECTIFIED = 1, MESH_PSC = 2;
	// the states for constructing tetra-elements during file reading
	private static final int TETRA_NEW = 0, TETRA_SEEK = 1, TETRA_CHECK_UNMATCHED = 2, TETRA_SEEK_UNMATCHED = 3;
	// these are the offsets into pertinent fieldvalues of tetrahedral precalculated data
	static final int EDGE01 = 0, EDGE02 = 1, EDGE03 = 2, EDGE12 = 3, EDGE13 = 4, EDGE23 = 5;
	static final int FACE012 = 6, FACE023 = 7, FACE031 = 8, FACE132 = 9;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			CONCURRENCY DEFINITIONS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	

	private static ExecutorService executor = null;			// thread pool for FEM ops
	static int taskNum, procNum = 1;
	static final AtomicInteger freeTasks = new AtomicInteger(procNum);
	public static List<Future<Void>> futureVoidList;		// list of return doubles from Callable threads
	public static List<Future<Double>> futureDoubleList;	// list of return doubles from Callable threads
	public static List<Future<Integer>> futureIntegerList;	// list of return integers from Callable threads
	static Thread taskList[];
	
	// initialise Executor (left to the user, since one can't get Runtime's processor count during class loading)
	public static ExecutorService getExecutor(int processors) {
		if (executor != null) return executor;
        // define the Callable executor pool of reusable threads
		if (processors == 0) 	executor = Executors.newFixedThreadPool(procNum = Runtime.getRuntime().availableProcessors() - 1);
		else					executor = Executors.newFixedThreadPool(procNum = processors);
		freeTasks.set(procNum);
  		// define a Future list of Callable return values constituting doubles (it will expand dynamically)
		futureVoidList = new ArrayList<Future<Void>>();	
		futureDoubleList = new ArrayList<Future<Double>>();	
		futureIntegerList = new ArrayList<Future<Integer>>();
		
		Runtime.getRuntime().addShutdownHook(new Thread() {			// schedule executor shutdown on exiting application
			public void run() {
				FEM1.getExecutor(0).shutdown();
				try { while (!FEM1.getExecutor(0).awaitTermination(10, TimeUnit.SECONDS));
				} catch (InterruptedException e) { e.printStackTrace(); }				
			}
		});		
		return executor;
	}

	// initialise Executor (left to the user, since one can't get Runtime's processor count during class loading)
	public static void initTaskList(int processors) {
        // define the Callable executor pool of reusable threads
		if (processors == 0) 	taskList = new Thread[procNum = Runtime.getRuntime().availableProcessors()];
		else					taskList = new Thread[procNum = processors];
		taskNum = 0;
	}

	public static void finishTaskList() {
		for (Thread task: taskList) {
			if (task == null) continue;
			try { task.join(); } catch (InterruptedException e) { e.printStackTrace(); }
		}
		taskNum = 0;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			CONSTRUCTORS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	

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
						if (maxNodesInPoly < poly.length) maxNodesInPoly = poly.length;
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
							boolean found1stIdentical = false;				// Search polygon for identical node indexes (signifies topology problems)
							for (int n = 0, nEnd = poly.length-1; n < nEnd; n++) {
								for (int n2 = n+1; n2 < poly.length; n2++)
									if (poly[n]==poly[n2]) {
										if (found1stIdentical)
												System.out.print("FEM1.FEM1(): polygon " + f + "contains duplicate node indexes " + n2);
										else	System.out.print(", " + n2);
										found1stIdentical = true;
									}
								if (found1stIdentical) System.out.println();
							}
							
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
			tetraNeighbours(null, true, false);
		}
		else if (dataType == MESH_PSC) {
			polygonsNodeSupport();					// group the polygons supported by every node
			facetNeighbours(false, false);			// find each face's neighbours, accept only facets, skip corner neighbours
			facetEdgeLengths();						// calculate edge lengths, generate nonredundant edge identities
			polygonNormals();						// calculate polygonal normals
		}
		
		getExecutor(0);								// start executor service
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
				n3 -= node.length; n = n3 / 3; nodeWork[n3++] = x; nodeWork[n3++] = y; nodeWork[n3] = z;
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
			if (nodeWorkFlag != null) nodeWorkFlag[nodesWork] = flag;
		} else {														// if we filled nodeWork[] array
			mergeNodeWork();											// merge it into node[] array
			int nodesWork2 = nodes / nodeworkFactor;
			if (nodesWork2 < 8) nodesWork2 = 8;
			nodeWork = new double[nodesWork2 * NCOORD];					// create new nodeWork[] array
			nodeWork[0] = x; nodeWork[1] = y; nodeWork[2] = z;
			if (nodeWorkFlag != null) nodeWorkFlag[0] = flag;
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
			
			int elemsS = support[0];
			if (elemsS + 1 >= support.length) {						// check if support array needs resizing
				int elemsS2 = 1 + elemsS + (elemsS >> 1);			// resize 1.5 times
				int[] support2 = new int[elemsS2];
				for (int n2 = 0; n2 <= elemsS; n2++) support2[n2] = support[n2];
				support = elementSupports[n] = support2;
			}
			support[0] = elemsS;
			support[elemsS] = e;
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
			tetraNeighbours(null, true, false);						// and we need neighbourhood arrays of every element, WITH neighbour interfaces
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
	void polygonsNodeSupport() {		
		int polygonsTot = polygonOffset[polygons], nCount = polygonOffset[1] - polygonOffset[0];	// polygonsTot = total node count of all polygons
		polygonSupports = new int[nodes][];
		for (int p = 0, p3 = 0, p3cnt = 0; p3 < polygonsTot; p3++) {

			if (p3cnt >= nCount) { p++; p3cnt = 0; nCount = -polygonOffset[p] + polygonOffset[p + 1]; }	// progress to next polygon
			int n = polygon[p3];
			int[] support = polygonSupports[n];
			
			if (support == null) {										// optimally, every node is shared by 6 faces or by 4 quartics (+ 1 array length integer)
				support = polygonSupports[n] = new int[9];				// allocate by defaul for 8 (+1) polygons adjoining this node
				support[0] = 1; support[1] = p;							// new support array started, insert it's first element
			} else {			
				int polys = support[0];
				if (polys + 2 >= support.length) {						// check if support array needs resizing
					int polys2 = polys + (polys >> 1);				// increase array 1.5 times
					int[] support2 = new int[polys2];
					for (int n2 = 0; n2 < polys; n2++) support2[n2] = support[n2];
					support = polygonSupports[n] = support2;
				}
				
				// for facetNeighbours() method to work properly, we must make sure the indexes are sorted
				if (support[polys] <= p) support[++polys] = p;		// if this index is higher than last one, just append
				//else if (support[nodesS] == p) throw new InvalidParameterException("FEM1.polygonsNodeSupport(): node supports polygon " + p + "twice.");
				else {
					for (int s = 1; s <= polys; s++) {					// otherwise insert polygon index in sorted order
						if (p < support[s]) {
							for (int s1 = polys + 1; s1 > s; s1++) support[s1] = support[s];
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
			if (support == null) { polygonSupports[n] = new int[2]; continue; }
			int polys = support[0] + 1;
			if (polygonSupportMaxArrayL < polys) polygonSupportMaxArrayL = polys;
			if (polys + 1 < support.length) {
				int[] support2 = new int[polys + 1];	// note: an extra slot on end added for optimisation of facetNeighbours() method
				for (int n2 = 0; n2 < polys; n2++) support2[n2] = support[n2];
				//for (int n2 = 0; n2 < polys; n2++) System.out.print(support2[n2] + (n2 == polys - 1 ? "\n" : (n2 == 0 ? ": " : ",")));
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
		int[] twinFacet = new int[8];										// found facet duplicates stored here (expandable array)
		int supportMult = onlyTrifaces ? 3 : maxNodesInPoly;
		
		// we'll first compare the nodes pairwise: min(n1, n2), then pick the lowest index: min(min(n1, n2), n3)
		int[] n12 = new int[polygonSupportMaxArrayL * 4 + 1];				// worst-case allocation: all uniques in both arrays -> x2 integers
		int[] n123e = new int[3];											// n123e[] carries only edge neighbours
		int[] n123eb = new int[3];											// n123ec[] carries only patch border edge neighbours
		int[] n123c = new int[polygonSupportMaxArrayL*supportMult];			// n123c[] carries only corner neighbours
		int[] n123cb = new int[polygonSupportMaxArrayL*supportMult];		// n123cc[] carries only patch border corner neighbours
		int pch = 0, pchStart = patch[pch], pchEnd = patch[pch];			// polygon indexes of current patch maintained by pchStart & pchEnd
		int twfN = 0, twfE = 0;												// counts number of duplicate facets, number of them erased
			
		// pB indexes current border facet, pBn indexes current patch count of border facets
		for (int p=0, p3=0, pB=0, pBn=0; p < polygons; p++) {	
			
			int ni = polygonOffset[p + 1] - polygonOffset[p];
			if (ni > 3) {
				if (onlyTrifaces) continue;											// if flagged to skip non-facets, jump over polygon
				else throw new InvalidParameterException("FEM1.faceNeighbours(): Non-facet encountered.");
			}				
			
			if (p >= pchEnd) {	pchStart = pchEnd;							// if we passed on into next patch, update patch delimitations
								pchEnd = patch[++pch];
								pBn = pB++;	}								// set aside a cell for border polygon count

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
					else { 
						boolean twinFacetStored = false;
						for (int tw=0; tw<twfN; tw++)													// see if twin facet already been encountered
							if (twinFacet[tw] == p || twinFacet[tw] == po3) { twinFacetStored=true;	break; }	
						if (twinFacetStored) {
						} else {
							polygon[--p3]=-1; polygon[--p3]=-1; polygon[--p3]=-1; p3+=3; twfE++;		// if yes, erase all it's twins, keeping only itself
							twinFacet[twfN++] = po3;													// 1st find of duplicate facet, store in twinFacet[]
							if (DEBUG_LEVEL > 1) System.out.println("FEM1.facetNeighbours(): duplicate facets, " + p + " and " + po12);
							if (twfN >= twinFacet.length) {	int[] twinFacetNew = new int [twfN*2];
															for (int tw2=0; tw2<twfN; tw2++) twinFacetNew[tw2] = twinFacet[tw2]; 
															twinFacet = twinFacetNew; }}
						po12 = n12[++i12]; po3 = a3[++i3];												// skip both element indexes
					}
					//else throw new InvalidParameterException("FEM1.facetNeighbours(): 2 facets neighbours on all three nodes.");
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
						// following case is topologically illegal, since two facets can only have two nodes equal, unless they mirror each other
					} else { 
						boolean twinFacetStored = false;
						for (int tw=0; tw<twfN; tw++)													// see if twin facet already been encountered
							if (twinFacet[tw] == p || twinFacet[tw] == po3) { twinFacetStored=true;	break; }	
						if (twinFacetStored) {
						} else {
							polygon[--p3]=-1; polygon[--p3]=-1; polygon[--p3]=-1; p3+=3; twfE++;		// if yes, erase all it's twins, keeping only itself
							twinFacet[twfN++] = po3;													// 1st find of duplicate facet, store in twinFacet[]
							if (DEBUG_LEVEL > 1) System.out.println("FEM1.facetNeighbours(): duplicate facets, " + p + " and " + po12);
							if (twfN >= twinFacet.length) {	int[] twinFacetNew = new int [twfN*2];
															for (int tw2=0; tw2<twfN; tw2++) twinFacetNew[tw2] = twinFacet[tw2]; 
															twinFacet = twinFacetNew; }}
						po12 = n12[++i12]; po3 = a3[++i3];												// skip both element indexes
					}
					//} else throw new InvalidParameterException("FEM1.facetNeighbours(): 2 facets neighbour on all three nodes.");
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
						if (po12c < pchStart || po12c >= pchEnd) n123cb[i123cb++] = po12c;				// if belonging to another patch, put in n123cb[]
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
		if (DEBUG_LEVEL > 1 && twfN > 0) System.out.println("Duplicate facets found: " + twfN + ", erased: " + twfE);
	}
	
	
	// method generates nonredundant mapping of every facet edge to it's length, where facet edges can share a length
	// the mapping is from a circular definition of an edge as point(n) to point(n+1), it index-mirrors the polygon nodes array polygon[] exactly,
	// and points to an array of edge lengths edge[] and an array of node index duplets edgeNode[] (for back-reference to underlying edge data)
	void facetEdgeLengths() {
		int polyNodes = polygonOffset[polygons];
		edgeIndex = new int[polyNodes];									// edgeIndex[] is the mapper to the edge lengths
		for (int i = 0; i < polyNodes; i++) edgeIndex[i] = -1;			// set all entries to unmapped, -1
		double[] edgeT = new double[polyNodes];							// worst-case allocation
		int[] edgeNT = new int[polyNodes * 2];							// worst-case allocation
		int eC = 0, eNP = 0;											// the edge counter & the node pair counter
		
		// every consecutive facet is checked if it's edges are already remapped by a previously checked neighbour
		// a new edge length is generated in case an edge hasn't been mapped yet
		for (int f = 0; f < polygons; f++) {
			int n = polygonOffset[f];
			if (polygon[n] == -1) continue;								// skip over erased polygons
			if (edgeIndex[n] == -1) {									// check if current facet edge n0-n1 hasn't been mapped yet
				edgeT[eC] = distance(edgeNT[eNP++] = polygon[n++], edgeNT[eNP++] = polygon[n--], true);
				edgeIndex[n] = eC++; } n++;								// map to unique position in polygonEdge[]		
			if (edgeIndex[n] == -1) {									// check if current facet edge n1-n2 hasn't been mapped yet
				edgeT[eC] = distance(edgeNT[eNP++] = polygon[n++], edgeNT[eNP++] = polygon[n--], true);
				edgeIndex[n] = eC++; } n++;								// map to unique position in polygonEdge[]				
			if (edgeIndex[n] == -1) {									// check if current facet edge n2-n0 hasn't been mapped yet
				edgeT[eC] = distance(edgeNT[eNP++] = polygon[n], edgeNT[eNP++] = polygon[n-2], true);
				edgeIndex[n] = eC++; }									// map to unique position in polygonEdge[]
			

			// check neighbours and remap them to the same indexes as current facet in polygonEdgeT[]
			int[] facetNgb = polygonNeighbours[f];
			// collect indexes of neighbour facets, from both normal & border edges
			int fN01 = facetNgb[5] == -1 ? facetNgb[8] : facetNgb[5], fN12 = facetNgb[6] == -1 ? facetNgb[9] : facetNgb[6];
			int fN20 = facetNgb[7] == -1 ? facetNgb[10] : facetNgb[7];
			n = polygonOffset[f];
			
			if (fN01 > f) {		// neighbours of lower index than current facet have already been fully mapped
				int[] ngb01 = polygonNeighbours[fN01];
				if		(ngb01[5] == f || ngb01[8] == f)	edgeIndex[polygonOffset[fN01]] = edgeIndex[n];
				else if (ngb01[6] == f || ngb01[9] == f)	edgeIndex[polygonOffset[fN01] + 1] = edgeIndex[n];
				else if (ngb01[7] == f || ngb01[10] == f)	edgeIndex[polygonOffset[fN01] + 2] = edgeIndex[n];
			}	n++;
			if (fN12 > f) {
				int[] ngb12 = polygonNeighbours[fN12];
				if		(ngb12[5] == f || ngb12[8] == f)	edgeIndex[polygonOffset[fN12]] = edgeIndex[n];
				else if (ngb12[6] == f || ngb12[9] == f)	edgeIndex[polygonOffset[fN12] + 1] = edgeIndex[n];
				else if (ngb12[7] == f || ngb12[10] == f)	edgeIndex[polygonOffset[fN12] + 2] = edgeIndex[n];
			}	n++;
			if (fN20 > f) {
				int[] ngb20 = polygonNeighbours[fN20];
				if		(ngb20[5] == f || ngb20[8] == f)	edgeIndex[polygonOffset[fN20]] = edgeIndex[n];
				else if (ngb20[6] == f || ngb20[9] == f)	edgeIndex[polygonOffset[fN20] + 1] = edgeIndex[n];
				else if (ngb20[7] == f || ngb20[10] == f)	edgeIndex[polygonOffset[fN20] + 2] = edgeIndex[n];
			}
		}
		
		// time to rescale edge[] & edgeNode[] to their final size
		avgEdgeLength = 0;
		minEdgeLength = Double.MAX_VALUE; maxEdgeLength = 0;
		edge = new double[edges = eC];
		edgeNode = new int[edges * 2];
		for (int e = 0, e2 = 0; e < edges; e++) {
			edge[e] = edgeT[e];
			if (edgeT[e] < minEdgeLength) minEdgeLength = edgeT[e];
			if (edgeT[e] > maxEdgeLength) maxEdgeLength = edgeT[e];
			avgEdgeLength += edgeT[e];
			edgeNode[e2] = edgeNT[e2++]; edgeNode[e2] = edgeNT[e2++];
		}
		avgEdgeLength /= (double)edges;
		// separate edges into 8 classes for simple heuristic choice of eventual tetrahedral sizing field
		edgeSpread = new SpreadBin(8, minEdgeLength, maxEdgeLength, 90);
		for (int e = 0; e < edges; e++) edgeSpread.add(edge[e]);
		if (DEBUG_LEVEL > 1) System.out.println("FEM1.facetEdgeLengths() edge lengths spread:\n" + edgeSpread.toString());
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
	int[][] tetraNeighbours(int[] elemList, boolean findInterfaces, boolean forIST) {

		if (!forIST && elementSupports == null) elementsNodeSupport();		// we need array of elements supported by each node (unless previously calculated)
	
		// we'll first compare the nodes pairwise: min(n1, n2) & min (n3, n4), then pick the lowest index: min(min(n1, n2), min (n3, n4))
		// pairwise comparisons go into these temporary arrays, then again compared into neighboursTmp[]
		int[] n12, n34, n1234F, n1234e;
		if (forIST) { n12 = nbrComp12; n34 = nbrComp34; n1234F = nbrComp1234F; n1234e = nbrComp1234e; }
		else {	n12 = new int[elementSupportMaxArrayL * 4 + 1];			// comparison of nodes 0 and 1
				n34 = new int[elementSupportMaxArrayL * 4 + 1];			// comparison of nodes 2 and 3
				n1234F = new int[elementSupportMaxArrayL * 8];			// comparisons 0&1 and 2&3, triple equalities go here as faces
				n1234e = new int[elementSupportMaxArrayL * 8]; }		// comparisons 0&1 and 2&3, single equalities go here as edges

		int elemCount = elemList == null ? elements2 : elemList.length;
		for (int e = 0, eI = 0, e5 = 0; e < elemCount; e++) {
			
			int[] a1, a2, a3, a4;
			FEM1Element elem1 = null;
			if (forIST) {													// for IST mode: get element lists supported by each node
				eI = elemList[e5++];
				a1 = elementNhoodIST[elemList[e5++]]; a2 = elementNhoodIST[elemList[e5++]];
				a3 = elementNhoodIST[elemList[e5++]]; a4 = elementNhoodIST[elemList[e5++]];
			} else {													// for object-oriented elements mode: get elements & related supports from their refs
				if (elemList == null) { elem1 = getElement2(e); eI = e; } else { elem1 = getElement2(elemList[e]); eI = elemList[e]; }
				if (elem1 == null) continue;							// skip recently deleted elements
				int[] nodeRef1 = elem1.nodeRef;
				a1 = elementSupports[nodeRef1[0]]; a2 = elementSupports[nodeRef1[1]];
				a3 = elementSupports[nodeRef1[2]]; a4 = elementSupports[nodeRef1[3]];
			}

			// we need to compare four indexes from four parallel sorted arrays, picking the lowest one every time
			// duplicates are stored normally, unique indexes stored intermediately with 31st bit set
			// a final unique index is a tetrahedron touching another tetrahedron with just one node, we're only interested in edges & faces
			int i12 = 0, i1 = 1, i2 = 1, c1 = a1[0], c2 = a2[0];
			int el1 = a1[i1], el2 = a2[i2];
			while (i1 <= c1 && i2 <= c2) {
				if (el1 == eI)		{ el1 = a1[++i1]; continue; }											// skip reference to current element itself
				if (el2 == eI)		{ el2 = a2[++i2]; continue; }											// (it will always be it's own neighbour)
				if (el1 < el2)		{ n12[i12++] = el1 | FL_UNIQUE; n12[i12++] = 0; el1 = a1[++i1]; }		// for uniques, additionally store the node index
				else if (el1 > el2)	{ n12[i12++] = el2 | FL_UNIQUE; n12[i12++] = 1; el2 = a2[++i2]; }
				else				{ n12[i12++] = el1;	el1 = a1[++i1]; el2 = a2[++i2]; }					// a neighbour on edge01 found										
			}
			// add remains as uniques, but skip self-references
			while (i1 <= c1) { if (a1[i1] == eI) i1++; else { n12[i12++] = a1[i1++] | FL_UNIQUE; n12[i12++] = 0; }}
			while (i2 <= c2) { if (a2[i2] == eI) i2++; else { n12[i12++] = a2[i2++] | FL_UNIQUE; n12[i12++] = 1; }}

			int i34 = 0, i3 = 1, i4 = 1, c3 = a3[0], c4 = a4[0];
			int el3 = a3[i3], el4 = a4[i4];
			while (i3 <= c3 && i4 <= c4) {
				if (el3 == eI)		{ el3 = a3[++i3]; continue; }											// skip reference to current element itself
				if (el4 == eI)		{ el4 = a4[++i4]; continue; }											// (it will always be it's own neighbour)
				if (el3 < el4)		{ n34[i34++] = el3 | FL_UNIQUE; n34[i34++] = 2; el3 = a3[++i3]; }		// for uniques, additionally store the node index
				else if (el3 > el4)	{ n34[i34++] = el4 | FL_UNIQUE; n34[i34++] = 3; el4 = a4[++i4]; }
				else				{ n34[i34++] = el3; el3 = a3[++i3]; el4 = a4[++i4]; }					// a neighbour on edge23 found
			}
			// add remains as uniques, but skip self-references
			while (i3 <= c3) { if (a3[i3] == eI) i3++; else { n34[i34++] = a3[i3++] | FL_UNIQUE; n34[i34++] = 2; }}
			while (i4 <= c4) { if (a4[i4] == eI) i4++; else { n34[i34++] = a4[i4++] | FL_UNIQUE; n34[i34++] = 3; }}

			int i1234F, i1234e, c12 = i12, c34 = i34;
			i12 = 0; i34 = 0;
			if (forIST) { i1234F = 1; i1234e = 1; } else { i1234F = 0; i1234e = 0; }
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
			
			if (forIST) {	n1234F[0] = --i1234F; n1234e[0] = --i1234e;
							int[][] neighbourhood = { n1234F, n1234e }; return neighbourhood; }	// for IST mode: we're done
			
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
				elem1.sortInterfaces(); }
			//for (eN = 0; eN < nTotal; eN+=2) System.out.print(neighbours[eN] + (eN == nTotal - 2 ? "\n" : (eN == 0 ? ": " : ",")));
		}
		return null;
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
	
	
	// returns bounding box of edge e
	public double[] edgeBBox(int e) {
		double[] bbox = new double[6];
		e *= 2; int n0_3 = edgeNode[e++] * NCOORD, n1_3 = edgeNode[e++] * NCOORD;
		if (nodeWork == null) {
			if (node[n0_3] < node[n1_3]) { bbox[0]=node[n0_3++]; bbox[3]=node[n1_3++]; } else { bbox[0]=node[n1_3++]; bbox[3]=node[n0_3++]; }
			if (node[n0_3] < node[n1_3]) { bbox[1]=node[n0_3++]; bbox[4]=node[n1_3++]; } else { bbox[1]=node[n1_3++]; bbox[4]=node[n0_3++]; }
			if (node[n0_3] < node[n1_3]) { bbox[2]=node[n0_3]; bbox[5]=node[n1_3]; } else { bbox[2]=node[n1_3]; bbox[5]=node[n0_3]; }
		} else {
			double[] node0, node1;
			if (n0_3 >= node.length) { node0 = nodeWork; n0_3 -= node.length; } else node0 = node;
			if (n1_3 >= node.length) { node1 = nodeWork; n1_3 -= node.length; } else node1 = node;
			if (node0[n0_3] < node1[n1_3]) { bbox[0]=node0[n0_3++]; bbox[3]=node1[n1_3++]; } else { bbox[0]=node1[n1_3++]; bbox[3]=node0[n0_3++]; }
			if (node0[n0_3] < node1[n1_3]) { bbox[1]=node0[n0_3++]; bbox[4]=node1[n1_3++]; } else { bbox[1]=node1[n1_3++]; bbox[4]=node0[n0_3++]; }
			if (node0[n0_3] < node1[n1_3]) { bbox[2]=node0[n0_3]; bbox[5]=node1[n1_3]; } else { bbox[2]=node1[n1_3]; bbox[5]=node0[n0_3]; }
		}
		return bbox;
	}

	// returns bounding box of facet f
	public double[] facetBBox(int f) {
		double[] bbox = new double[6];
		int fO = polygonOffset[f];
		if (polygon[fO] == -1 || polygonOffset[f + 1] - fO > 3) return null;			// accept only unerased 3-node polygons, fail otherwise
		int n0_3 = polygon[fO++] * NCOORD, n1_3 = polygon[fO++] * NCOORD, n2_3 = polygon[fO] * NCOORD;
		if (nodeWork == null) {
			if (node[n0_3] < node[n1_3]) { bbox[0]=node[n0_3++]; bbox[3]=node[n1_3++]; } else { bbox[0]=node[n1_3++]; bbox[3]=node[n0_3++]; }
			if (node[n0_3] < node[n1_3]) { bbox[1]=node[n0_3++]; bbox[4]=node[n1_3++]; } else { bbox[1]=node[n1_3++]; bbox[4]=node[n0_3++]; }
			if (node[n0_3] < node[n1_3]) { bbox[2]=node[n0_3];   bbox[5]=node[n1_3]; } else {   bbox[2]=node[n1_3];   bbox[5]=node[n0_3]; }
			if (node[n2_3] < bbox[0]) bbox[0]=node[n2_3++]; else if (node[n2_3] > bbox[3]) bbox[3]=node[n2_3++]; else n2_3++;
			if (node[n2_3] < bbox[1]) bbox[1]=node[n2_3++]; else if (node[n2_3] > bbox[4]) bbox[4]=node[n2_3++]; else n2_3++;
			if (node[n2_3] < bbox[2]) bbox[2]=node[n2_3]; else   if (node[n2_3] > bbox[5]) bbox[5]=node[n2_3];
		} else {
			double[] node0, node1, node2;
			if (n0_3 >= node.length) { node0 = nodeWork; n0_3 -= node.length; } else node0 = node;
			if (n1_3 >= node.length) { node1 = nodeWork; n1_3 -= node.length; } else node1 = node;
			if (n2_3 >= node.length) { node2 = nodeWork; n2_3 -= node.length; } else node2 = node;
			if (node0[n0_3] < node1[n1_3]) { bbox[0]=node0[n0_3++]; bbox[3]=node1[n1_3++]; } else { bbox[0]=node1[n1_3++]; bbox[3]=node0[n0_3++]; }
			if (node0[n0_3] < node1[n1_3]) { bbox[1]=node0[n0_3++]; bbox[4]=node1[n1_3++]; } else { bbox[1]=node1[n1_3++]; bbox[4]=node0[n0_3++]; }
			if (node0[n0_3] < node1[n1_3]) { bbox[2]=node0[n0_3];   bbox[5]=node1[n1_3]; } else {   bbox[2]=node1[n1_3];   bbox[5]=node0[n0_3]; }
			if (node2[n2_3] < bbox[0]) bbox[0]=node2[n2_3++]; else if (node2[n2_3] > bbox[3]) bbox[3]=node2[n2_3++]; else n2_3++;
			if (node2[n2_3] < bbox[1]) bbox[1]=node2[n2_3++]; else if (node2[n2_3] > bbox[4]) bbox[4]=node2[n2_3++]; else n2_3++;
			if (node2[n2_3] < bbox[2]) bbox[2]=node2[n2_3]; else   if (node2[n2_3] > bbox[5]) bbox[5]=node2[n2_3];
		}
		return bbox;
	}

	
	
	// given three existing nodes (n0,n1,n2), spawns a tetrahedron with optimal 4th node position
	// method demands locally-proper-ordered input nodes: anticlockwise from point of view of the new tetrahedron
	// nodes could be from spawning tetrahedron or boundary facet, 4th node (n3) must exist, and it's coordinates will be replaced
	private static double SQRT6D3 = 0.8164965809277259;
	public FEM1Element optimalTetrahedron(int n0, int n1, int n2, int n3) {
		int n0_3 = n0 * NCOORD, n1_3 = n1 * NCOORD, n2_3 = n2 * FEM1.NCOORD, n3_3 = n3 * FEM1.NCOORD;
		double[] node0, node1, node2, node3;
		if (nodeWork == null) {
			node0 = node1 = node2 = node3 = node;
		} else {
			if (n0_3 >= node.length) { node0 = nodeWork; n0_3 -= node.length; } else node0 = node;
			if (n1_3 >= node.length) { node1 = nodeWork; n1_3 -= node.length; } else node1 = node;
			if (n2_3 >= node.length) { node2 = nodeWork; n2_3 -= node.length; } else node2 = node;
			if (n3_3 >= node.length) { node3 = nodeWork; n3_3 -= node.length; } else node3 = node;
		}
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
	

	
	// method calculates a triangles (planes) normal, supplied nodes must exist, normal destination array must be supplied and offset into it
	// returns false if triangle area was zero
	// float precision is more than enough for defining a normal, and takes less space
	boolean triangleNormal(int n1, int n2, int n3, float[] normal, int offs) {
		int n1_3 = n1 * 3, n2_3 = n2 * 3, n3_3 = n3 * 3;
		double[] node1, node2, node3;
		if (nodeWork == null) {
					node1 = node2 = node3 = node;
		} else {	if (n1_3 >= node.length) { n1_3 -= node.length; node1 = nodeWork; } else node1 = node;
					if (n2_3 >= node.length) { n2_3 -= node.length; node2 = nodeWork; } else node2 = node;
					if (n3_3 >= node.length) { n3_3 -= node.length; node3 = nodeWork; } else node3 = node; }
		double x1 = node1[n1_3++], y1 = node1[n1_3++], z1 = node1[n1_3];
		double x21 = node2[n2_3++] - x1, y21 = node2[n2_3++] - y1, z21 = node2[n2_3] - z1;
		double x31 = node3[n3_3++] - x1, y31 = node3[n3_3++] - y1, z31 = node3[n3_3] - z1;
		double xn = y21*z31 - z21*y31, yn = z21*x31 - x21*z31, zn = x21*y31 - y21*x31, nLinv = Math.sqrt(xn*xn + yn*yn + zn*zn);
		if (nLinv == 0) return false; else nLinv = 1. / nLinv;					// if face area = zero, return null
		normal[offs++] = (float)(xn * nLinv); normal[offs++] = (float)(yn * nLinv); normal[offs] = (float)(zn * nLinv);
		return true;
	}
	
	boolean facetNormal(int f, float[] normal, int offs) { return triangleNormal(polygonOffset[f++], polygonOffset[f++], polygonOffset[f], normal, offs); }
	
	
	final static double FACET_MARGIN = 1e-8;
	// method calculates intersection point of line segment passing through a facet, method demands that polygonNormal[] array is constructed
	// if pa behind facet, returns 2, if pa outside, returns -2, if no intersection exists, returns 0 and zeros in isect[]
	// if isect[] = null, does not calculate intersection point
	int facetSegmentIntersection(int f, double xa, double ya, double za, double xb, double yb, double zb, double[] isect) {
		
		int fN = f * 3, status = 0;
		f = polygonOffset[f];
		int n0_3 = polygon[f++] * 3, n1_3 = polygon[f++] * 3, n2_3 = polygon[f] * 3;
		double x0 = node[n0_3++], y0 = node[n0_3++], z0 = node[n0_3];
		double x0a = x0 - xa, y0a = y0 - ya, z0a = z0 - za, xba = xb-xa, yba = yb-ya, zba = zb-za, rI012;
		double xN = polygonNormal[fN++], yN = polygonNormal[fN++], zN = polygonNormal[fN];
		if ((rI012 = xN * x0a + yN * y0a + zN * z0a) >= 0) status = 1; else status = -1;
		double xb0 = xb - x0, yb0 = yb - y0, zb0 = zb - z0;
		if (xN * xb0 + yN * yb0 + zN * zb0 <= 0) status--; else status++;
		// status = 0 -> pa & pb outside facet plane, a non-intersecting case
		if (status == 0) return 0;

		double rI012D = xN*xba + yN*yba + zN*zba;
		if (rI012D != 0) {					// if segment & facet are not coplanar/parallel
			rI012 = rI012 / rI012D;
			double xv012 = xa + rI012*xba - x0, yv012 = ya + rI012*yba - y0, zv012 = za + rI012*zba - z0;	// plane isect point - node0 vector
			double x10 = node[n1_3++]-x0, y10 = node[n1_3++]-y0, z10 = node[n1_3]-z0, x20 = node[n2_3++]-x0, y20 = node[n2_3++]-y0, z20 = node[n2_3]-z0;
			double v10v20 =	x10*x20 + y10*y20 + z10*z20;
			double v012v20 = xv012*x20 + yv012*y20 + zv012*z20, v20v20 = x20*x20 + y20*y20 + z20*z20, v012v10 = xv012*x10 + yv012*y10 + zv012*z10;
			double v10v10 = x10*x10 + y10*y10 + z10*z10, st012D = 1. / (v10v20*v10v20 - v10v10*v20v20);
			double sI012 = (v10v20 * v012v20 - v20v20 * v012v10) * st012D, tI012 = (v10v20 * v012v10 - v10v10 * v012v20) * st012D;
			//if (sI012 >= 0 && tI012 >= 0 && sI012 + tI012 <= 1) {	// if intersection point is inside parameterised triangle
			// DEBUG: parameterised triangle area somewhat expanded
			if (sI012 >= -FACET_MARGIN && tI012 >= -FACET_MARGIN && sI012 + tI012 <= 1+FACET_MARGIN) {
				if (isect != null) { isect[0] = x0+sI012*x10+tI012*x20; isect[1] = y0+sI012*y10+tI012*y20; isect[2] = z0+sI012*z10+tI012*z20; }
				return status;
			} 
		}
		return 0;
	}

	
	// method generates normals for all the FEM1-constituent polygons (failed/zero-area polygons return zeroed normals)
	public float[] polygonNormals() {
		polygonNormal = new float[polygons * 3];
		for (int p = 0, pN3 = 0; p < polygons; p++, pN3 += 3) {
			// if the three chosen nodes form a zero area, try next node as third coordinate
			int pS = polygonOffset[p];
			if (polygon[pS] == -1) continue;
			for (int pE = pS + 2, p2End = polygonOffset[p + 1]; pE < p2End; pE++)
				if (triangleNormal(polygon[pS], polygon[pS + 1], polygon[pE], polygonNormal, pN3)) break;
		}
		// since normals are maent for fast checking of interferences with polygons,
		// this is a fitting place to instantiate a checklist for visited polygons
		if (polygonCheck == null) polygonCheck = new VisitBitArray(polygons);
		return polygonNormal;
	}
	
	
	// method calculates the optimal position of an internal node within its 1-ring neighbourhood, using the formula of
	// the weighted barycenter of circumcenters divided by total 1-ring volume (an energy minimisation formulation)
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
					float[] fNormal1 = {0,0,0}, fNormal2 = {0,0,0}, fNormal3 = {0,0,0};
					// get sums from the found boundary facets, depending on what node n corresponds to on tetrahedron
					// fNormal1/2/3 is the inward-pointing normal of the boundary facet, scaled by area of that facet, N(i,p,q)
					// this sums up boundary factors N(i,p,q) * (||x(p)-x(i)||^2 - ||x(q)-x(i)||^2) from every neighbour's boundary facet
					if (nodeRef[0] == n) {
						int ref0 = nodeRef[0];
						if ((internalAreas&1) == 0) { double sum = (elem.data[EDGE01] + elem.data[EDGE02]) * elem.data[FACE012];	// facet 012
							triangleNormal(ref0, nodeRef[2], nodeRef[1], fNormal1, 0); fNormal1[0] *= sum; fNormal1[1] *= sum; fNormal1[2] *= sum; }
						if ((internalAreas&2) == 0) { double sum = (elem.data[EDGE02] + elem.data[EDGE03]) * elem.data[FACE023];	// facet 023
							triangleNormal(ref0, nodeRef[3], nodeRef[2], fNormal2, 0); fNormal2[0] *= sum; fNormal2[1] *= sum; fNormal2[2] *= sum; }
						if ((internalAreas&4) == 0) { double sum = (elem.data[EDGE01] + elem.data[EDGE03]) * elem.data[FACE031];	// facet 031
							triangleNormal(ref0, nodeRef[1], nodeRef[3], fNormal3, 0); fNormal3[0] *= sum; fNormal3[1] *= sum; fNormal3[2] *= sum; }
					} else if (nodeRef[1] == n) {
						int ref1 = nodeRef[1];
						if ((internalAreas&1) == 0) { double sum = (elem.data[EDGE01] + elem.data[EDGE12]) * elem.data[FACE012];	// facet 012
							triangleNormal(nodeRef[0], nodeRef[2], ref1, fNormal1, 0); fNormal1[0] *= sum; fNormal1[1] *= sum; fNormal1[2] *= sum; }
						if ((internalAreas&4) == 0) { double sum = (elem.data[EDGE01] + elem.data[EDGE13]) * elem.data[FACE031];	// facet 031
							triangleNormal(nodeRef[0], ref1, nodeRef[3], fNormal2, 0); fNormal2[0] *= sum; fNormal2[1] *= sum; fNormal2[2] *= sum; }
						if ((internalAreas&8) == 0) { double sum = (elem.data[EDGE13] + elem.data[EDGE12]) * elem.data[FACE132];	// facet 132
							triangleNormal(ref1, nodeRef[2], nodeRef[3], fNormal3, 0); fNormal3[0] *= sum; fNormal3[1] *= sum; fNormal3[2] *= sum; }
					} else if (nodeRef[2] == n) {
						int ref2 = nodeRef[2];
						if ((internalAreas&1) == 0) { double sum = (elem.data[EDGE02] + elem.data[EDGE12]) * elem.data[FACE012];	// facet 012
							triangleNormal(nodeRef[0], ref2, nodeRef[1], fNormal1, 0); fNormal1[0] *= sum; fNormal1[1] *= sum; fNormal1[2] *= sum; }
						if ((internalAreas&2) == 0) { double sum = (elem.data[EDGE02] + elem.data[EDGE23]) * elem.data[FACE023];	// facet 023
							triangleNormal(nodeRef[2], nodeRef[3], ref2, fNormal2, 0); fNormal2[0] *= sum; fNormal2[1] *= sum; fNormal2[2] *= sum; }
						if ((internalAreas&8) == 0) { double sum = (elem.data[EDGE23] + elem.data[EDGE12]) * elem.data[FACE132];	// facet 132
							triangleNormal(nodeRef[1], ref2, nodeRef[3], fNormal3, 0); fNormal3[0] *= sum; fNormal3[1] *= sum; fNormal3[2] *= sum; }
					} else /*if (nodeRef[3] == n)*/ {
						int ref3 = nodeRef[3];
						if ((internalAreas&2) == 0) { double sum = (elem.data[EDGE23] + elem.data[EDGE03]) * elem.data[FACE023];	// facet 023
							triangleNormal(nodeRef[0], ref3, nodeRef[2], fNormal1, 0); fNormal1[0] *= sum; fNormal1[1] *= sum; fNormal1[2] *= sum; }
						if ((internalAreas&4) == 0) { double sum = (elem.data[EDGE03] + elem.data[EDGE13]) * elem.data[FACE031];	// facet 031
							triangleNormal(nodeRef[0], nodeRef[1], ref3, fNormal2, 0); fNormal2[0] *= sum; fNormal2[2] *= sum; fNormal2[2] *= sum; }
						if ((internalAreas&8) == 0) { double sum = (elem.data[EDGE13] + elem.data[EDGE23]) * elem.data[FACE132];	// facet 132
							triangleNormal(nodeRef[1], nodeRef[2], ref3, fNormal3, 0); fNormal3[0] *= sum; fNormal3[3] *= sum; fNormal3[2] *= sum; }
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
	public int closestNodePairs(FEM1Octree octree, int[] closest, double[] distance, FEM1Octant[] closestOctants) {
		double[] distC = {0};
		FEM1Octant[] octantsC = {null, null};
		
		// flag every node tuple as "unassigned", for skipping symmetric closeness cases (n1 closest to n2 & n2 closest to n1 -> do not store (n2,n1))
		for (int n2 = 0, nodes2 = nodes * 2; n2 < nodes2; n2++) closest[n2++] = -1;
		
		int nPairs = 0;
		FEM1Octant[] leafOctant = octree.root.octantArray(octree, -1, null, false);	// collect leaf octants

		for (FEM1Octant octant : leafOctant) {
			if (octant == null) break;
			octantsC[0] = octant;
			for (int nI = 0; nI < octant.nodes; nI++) {								// iterate over total number of nodes in octant
				octantsC[1] = null;
				int n = octant.nodeI[nI], n2 = n * 2;
				int nClosest = octree.root.closestNode(octree, n, distC, octantsC);
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
	double[] addMediansOfClosestNodes(FEM1Octree octree, int nPairs, int[] closest, FEM1Octant[] closestOctants, FEM1Octant[] medianOctants) {
		
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
			FEM1Octant octant = null;
			if (closestOctants[nP2].octantDistance(xm, ym, zm, false) == 0)	{
				octant = closestOctants[nP2];
				octant.octantAddNode(nM3 / 3, true);
			} else if (closestOctants[nP2 + 1].octantDistance(xm, ym, zm, false) == 0) {
				octant = closestOctants[nP2 + 1];
				octant.octantAddNode(nM3 / 3, true);
			}
			// the case of the median crossing over to a corner octant to the two neighbour octants:
			// find the topmost container and do an extraneous node insertion
			else	octant = octree.root.addNode(octree, xm, ym, zm, nM3 / 3, true);	
			
			if (medianOctants != null) medianOctants[nP] = octant;
		}
		return median;
	}

		
	// method returns the node in supplied index array that is closest to supplied node index, distance returned in supplied distance[] array
	// utilising the possible existence of additional nodes in nodeWork[]
	// note: comparison happens with squared distances, and the distance is returned squared!
	int closestNode(int n1, int[] nL, double[] distance) {	
		if ((n1 *= 3) >= node.length) {
			n1 -= node.length;
				return closestNode(nodeWork[n1++], nodeWork[n1++], nodeWork[n1], nL, distance); }
		else 	return closestNode(node[n1++], node[n1++], node[n1], nL, distance);
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
	//			ISOSURFACE STUFFING TETRAHEDRAL SUBDIVIDER
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	

	// global parameters
	boolean newSolution = true;
	double latticeStep = 0;							// the cube width of boundary lattice
	double minElemQuality = 0.2;					// minimum defectivity allowed for an element to pass final quality check
	double minUnsnapQuality = 0.4;					// minimum quality to attain when unsnapping a bad quality element
	double minSliverQuality = 0.2;					// minimum defectivity allowed for a sliver (element with all nodes snapped to surface)
	double alphaS = 0.25, alphaL = /*0.412*/.35;	// these factors tell whether a node is "violated" by a cut node and a node snap is merited
	double maxZordinDev = 1e-8;						// max deviation of an edge from Z-ordinate to consider using the Z-intersection data for the edge

	boolean discardAllSliversIST = false;			// discards any sliver found (can help remove SOME straddling elements)
	// note: following switches are a test-only attempt to improve volume fit for low-res, and will never give a 100% satisfactory solution,
	// complex heuristics would have to be applied to nearly all elements, the computations would be costlier than generating a higher resolution mesh,
	// therefore: the only thing that will give a guaranteedly satisfactory solution is to increase resolution!
	boolean discardStraddlers = false;				// discards elements straddling gaps in the mesh (note: will leave holes & increase computation)
	boolean discardBadSliversIST = true;			// flags for discarding elements fully snapped onto surface that are bad
	// debug switches
	boolean skipVolumeInternalsIST = false;			// DEBUG: skips elements inside of the hull (good for surface quality previews)
	int maxOutputElementsIST = 0;					// DEBUG: write out up to this number number of elements (set = 0 to ignore)
	
	// lattice structures
	int latticeSubdivs=0, minLatticeSubdivs=16;		// the specified and the minimum allowed number of lattice cubes along all three ordinates
	double subsPerSize;								// subdivs divided by lattice size and "zero" point of lattice (=xM/yM/zM of underlying octree)
	double[][][] fEncounterGrid = null;				// holds a latticed grid of facet encounters, sorted coordinate-wise by z-distance
	int[][] fStatusGrid = null;						// holds encounter counts and polarities of the abovementioned
	double[] bLatticeNode = null;					// holds coordinates of the boundary lattice, to be formed into an octree
	short[] bLatticeStatus = null;					// holds inside/outside flags for every boundary lattice cube, the array index-matches bLatticeNode[]
	
	// datastructures
	int tetIST0, tetIST1, tetIST2, tetIST3;			// the 4 current/running node indexes for a tetrahedron element	
	int tetRef0, tetRef1, tetRef2, tetRef3;			// the 4 current/running node's element origins for a tetrahedron element
	int slotcodeIST=0, edgeCodeIST=0, mirrorFlagIST=0;	// the data held for an element writedown in processISTnode()

	int nodesIST = 0, streamNodesIST = 0;			// how many unique nodes were accumulated, how many nodes expected in a node stream
	double[] nodeIST = null;						// final array for unique IST solution nodes
	NodeHashTable nodeHTableIST = null;
	int[] cutNodes = null;							// indexes to nodes created during edge cut tests are put here
	byte[] snapNodeStatus = null;					// a status set when a node is snapped to a cut point
	
	int elementsIST = 0, elementsIntIST = 0;		// counter for total number of IST internal and boundary elements & number of internal elements
	int cutEdgesIST = 0;							// counter for total number of unique cutedges
	final int elementPgSize=1000*6;					// page sizes of elementPg[] (must be evenly dividable by 5 and 6)
	int ePg = 0;									// page accumulation counters
	int elementPgI = 0, elementPgN=400;				// page counter for elementPg[] & total pages count
	int[][] elementPg = null;						// paging structure for element array
	int[] elementIST = null;						// paged boundary elements are regathered into this array
	short[] edgeOctFlags = null;					// holds child flags for every edge neighbour of an IST octant (bit set if children exist on that edge)
	int[][] elementNhoodIST = null;
	
	double[] elementStream = null;					// final element stream
	long[] elementStreamCode = null;				// final element bitpacked compaction code stream
	
	// arrays used for neighbourhood comparison during nonthreaded stage 4 boundary slivers analysis
	int[] nbrComp12=null, nbrComp34=null, nbrComp1234F=null, nbrComp1234e=null;
	
	// IST debugging counters
	private int debugIST_nonsharedEdge = 0, debugIST_nonsharedEdge2 = 0, debugIST_locallyUnique = 0;
	private int debugIST_pISTe_visits = 0, debugIST_pISTe_vainSnaps = 0, debugIST_pISTe_fastVisits = 0, debugIST_pISTe_fastFacetVisits = 0;
	private int debugIST_pISTe_skipVisits = 0;
	private long debugIST_phase1time=0, debugIST_phase2time=0, debugIST_phase3time=0, debugIST_phase4time=0, debugIST_phase5time=0;
	private long debugIST_pISTe_time1 = 0, debugIST_pISTe_time2 = 0, debugIST_pISTn_time = 0;
	private long debugIST_bLatt_time = 0, debugIST_aBOIST_time = 0, debugIST_analysisTime = 0;
	
	public FEM1Octree latticeTree(FEM1Octree geocTree, int minSubdivs, boolean multitask) {
		boundaryLattice(geocTree, minSubdivs, true, true);					// lattice data with centroid encounters and inner volume cubes generated
		FEM1Octree latticeTree = new FEM1Octree(geocTree, 4);
		latticeTree.root.build(latticeTree, 0, multitask);
		bLatticeNode = null; bLatticeStatus = null;							// these arrays are not needed anymore
		extendBCCdoNhoodsIST(latticeTree);									// extra leaf octants made, ensuring only BCC tetrahedra intersect isosurface
		latticeTree.root.process(latticeTree, FEM1Octant.OCTREE_SUM_BRANCHES);
		return latticeTree;
	}
	
	// method generates an isosurface stuffing tetrahedral volume mesh
	// subdivScale gives the size of the lattice that will subdivide the mesh in metric scale, and need to be at most half of the scale
	// of the smallest features that need to be captured
	public FEM1Octree volumeMeshIST(FEM1Octree geocTree, double subdivScale, int gradations, boolean multitask) {
		
		if (!newSolution) clearISTsolution();
		if (minUnsnapQuality > 0.64) { System.out.println("FEM1.volumeMeshIST(): minUnsnapQuality too high, setting to " + (minUnsnapQuality=0.64)); }

		FEM1Octree latticeTree = null;
		int subdivisions = 1, divFactor;
		if (subdivScale > 0) {
			subdivisions = (int)((geocTree.root.xP - geocTree.root.xM) / subdivScale);
		} else {
			//int maxBin = edgeSpread.maxBin();											// find out what group of facet sizes is the most numerous
			long bin0 = edgeSpread.bin[0], bin1 = edgeSpread.bin[1];
			divFactor = geocTree.maxLevelFromGranularity((edgeSpread.div[0]*bin0 + edgeSpread.div[1]*bin1) / ((bin0+bin1)*2));
			while (divFactor-- > 0) subdivisions *= 2;
			if (DEBUG_LEVEL > 1) System.out.println("FEM1.volumeMeshIST(): calculated lattice subdivision from largest edge group: " + subdivisions);			
		}
		if (subdivisions < minLatticeSubdivs) {
			if (DEBUG_LEVEL > 1)
				System.out.println("FEM1.volumeMeshIST(): subdivision of "+subdivisions+" too low, setting to minimum value: "+minLatticeSubdivs);
			subdivisions = minLatticeSubdivs; }
		latticeSubdivs = subdivisions;
		boundaryLattice(geocTree, subdivisions, true, true);					// create the marching cubes Z-ordinate intersection data
		latticeTree = new FEM1Octree(geocTree, 4);
		bBox[0] = latticeTree.root.xM; bBox[1] = latticeTree.root.yM; bBox[2] = latticeTree.root.zM;
		bBox[3] = latticeTree.root.xP; bBox[4] = latticeTree.root.yP; bBox[5] = latticeTree.root.zP;
		double granularity = (bBox[3] - bBox[0]) / latticeSubdivs;
		divFactor = geocTree.maxLevelFromGranularity(granularity);
		latticeTree.gradations = gradations<=0 || gradations >= divFactor ? 4 : gradations;
		maxZordinDev = granularity / 20.0;										// how far an edge can deviate from Z-ordinate to be considered to be on it
		latticeTree.root.build(latticeTree, 0, multitask);						// build the lattice tree
		bLatticeNode = null; bLatticeStatus = null; 							// these arrays are not needed anymore
		// this call makes extra octant leaves according to the corners vs. centroid int/ext criterion
		// and also forms the face neighbourhoods and bitfields for fast checking of children at edges
		FEM1Octant[] nhoodIST = extendBCCdoNhoodsIST(latticeTree);				// extra leaf octants made, ensuring only BCC tets intersect isosurface
		sortStencilPermCode();
		//FEM1Octant testOctant = latticeTree.root.locateCoordinate(0.398, 0.407, -0.757);			// DEBUG: finds an erroneous BCC octant
		elementGeneratorIST(latticeTree, geocTree, nhoodIST, 0, 0);	
		newSolution = false;
		//int errorElement = findElementIndex(element, nodeIST, 2.0373, -0.6068, -1.6933, 0, 3);	// DEBUG: finds an erroneous final element
		return latticeTree;
	}

	// TODO: the version that generates a mesh volume towards a node or element count target
	public FEM1Octree volumeMeshIST(FEM1Octree geocTree, double subdivScale, boolean multitask, int elementTarget, int nodeTarget) {
		
		if (!newSolution) clearISTsolution();
		if (minUnsnapQuality > 0.64) { System.out.println("FEM1.volumeMeshIST(): minUnsnapQuality too high, setting to " + (minUnsnapQuality=0.64)); }

		FEM1Octree latticeTree = null;
		int subdivisions = 1;
		if (subdivScale > 0) {
			subdivisions = (int)((geocTree.root.xP - geocTree.root.xM) / subdivScale);
		} else {
			int maxBin = edgeSpread.maxBin();											// find out what group of facet sizes is the most numerous
			long bin0 = edgeSpread.bin[0], bin1 = edgeSpread.bin[1];
			int divFactor = geocTree.maxLevelFromGranularity((edgeSpread.div[0]*bin0 + edgeSpread.div[1]*bin1) / ((bin0+bin1)*2));
			while (divFactor-- > 0) subdivisions *= 2;
			if (DEBUG_LEVEL > 1) System.out.println("FEM1.volumeMeshIST(): calculated lattice subdivision from largest edge group: " + subdivisions);			
		}
		if (subdivisions < minLatticeSubdivs) {
			if (DEBUG_LEVEL > 1)
				System.out.println("FEM1.volumeMeshIST(): subdivision of "+subdivisions+" too low, setting to minimum value: "+minLatticeSubdivs);
			subdivisions = minLatticeSubdivs; }
		latticeSubdivs = subdivisions;
		// define upper and lower bounds for search of an optimal subdivision to match element & node targets
		int lwrSubdiv = (subdivisions>>1)<minLatticeSubdivs ? minLatticeSubdivs : subdivisions>>1;
		int uprSubdiv = latticeSubdivs==minLatticeSubdivs ? minLatticeSubdivs*2 : latticeSubdivs;
		if (elementTarget < 0) elementTarget = 0; if (nodeTarget < 0) nodeTarget = 0;
		boolean approximateTargets = (elementTarget > 0 || nodeTarget > 0) ? true : false;		
		boolean elementTargetOK = true, nodeTargetOK = true;
		int[] enReturn = null;
		
		do {
			if (!newSolution) clearISTsolution();
			boundaryLattice(geocTree, subdivisions, true, true);					// create the marching cubes Z-ordinate intersection data
			latticeTree = new FEM1Octree(geocTree, 4);
			bBox[0] = latticeTree.root.xM; bBox[1] = latticeTree.root.yM; bBox[2] = latticeTree.root.zM;
			bBox[3] = latticeTree.root.xP; bBox[4] = latticeTree.root.yP; bBox[5] = latticeTree.root.zP;
			maxZordinDev = (bBox[3] - bBox[0]) / (latticeSubdivs * 20);				// how far an edge can deviate from Z-ordinate to be considered to be on it
			latticeTree.root.build(latticeTree, 0, multitask);
			bLatticeNode = null; bLatticeStatus = null; 							// these arrays are not needed anymore
			// this call makes extra octant leaves according to the corners vs. centroid int/ext criterion
			// and also forms the face neighbourhoods and bitfields for fast checking of children at edges
			FEM1Octant[] arrayIST = extendBCCdoNhoodsIST(latticeTree);				// extra leaf octants made, ensuring only BCC tetrahedra intersect isosurface
			sortStencilPermCode();
			enReturn = elementGeneratorIST(latticeTree, geocTree, arrayIST, elementTarget, nodeTarget);
			newSolution = false;
			if (enReturn[0] == 0) {
				if (!approximateTargets) break;
				subdivisions=(uprSubdiv-=(uprSubdiv-lwrSubdiv)>>1); 
			} else if (enReturn[0] < 0) {
				if (subdivisions<=minLatticeSubdivs) { if (DEBUG_LEVEL>1) System.out.println("FEM1.volumeMeshIST(): element target unattainable."); break; }
				subdivisions=(uprSubdiv-=(uprSubdiv-lwrSubdiv)>>1); elementTargetOK = false; }
			else if ((lwrSubdiv+(uprSubdiv-lwrSubdiv)>>1) >= uprSubdiv - 1) {
				 if (DEBUG_LEVEL>1) System.out.println("FEM1.volumeMeshIST(): element target attained."); break;
			} else elementTargetOK = true;
			if (enReturn[1] < 0) {
				if (subdivisions<=minLatticeSubdivs) { if (DEBUG_LEVEL>1) System.out.println("FEM1.volumeMeshIST(): element target unattainable."); break; }
				subdivisions=(uprSubdiv>>=1); elementTargetOK = false; }
			else elementTargetOK = true;
		} while (approximateTargets && (!elementTargetOK || !nodeTargetOK) && lwrSubdiv < uprSubdiv);
		return latticeTree;
	}

	
	
	public void clearISTsolution() {
		latticeStep = subsPerSize = 0;
		latticeSubdivs = 0; minLatticeSubdivs=16;
		fEncounterGrid = null; fStatusGrid = null;
		bLatticeNode = null; bLatticeStatus = null; elementNhoodIST = null;
		
		tetIST0 = tetIST1 = tetIST2 = tetIST3 = tetRef0 = tetRef1 = tetRef2 = tetRef3 = slotcodeIST = nodesIST = streamNodesIST = 0;
		nodeHTableIST = null;
		snapNodeStatus = null;		
		elementsIST = elementsIntIST = cutEdgesIST = ePg = elementPgI = 0;
		elementPg = null; cutNodes = elementIST = null; elementPgN=400;
		edgeOctFlags = null;	
		nodeIST = elementStream = null; elementStreamCode = null;		
		debugIST_nonsharedEdge = debugIST_nonsharedEdge2 = debugIST_locallyUnique = 0;
		debugIST_pISTe_visits = debugIST_pISTe_vainSnaps = 0;
		debugIST_phase1time = debugIST_phase2time = debugIST_phase3time = debugIST_phase4time = debugIST_phase5time = 0;
		debugIST_pISTe_time1=debugIST_pISTe_time2=debugIST_pISTn_time=debugIST_bLatt_time=debugIST_pISTe_fastVisits=debugIST_pISTe_fastFacetVisits=0;
		debugIST_aBOIST_time = debugIST_analysisTime = 0; debugIST_pISTe_skipVisits = 0;
		readElementPg = true; finalArrayPg = false; printSlotcodes1 = false; printSlotcodes2 = false; printSnapNodes = false;
		newSolution = true;
	}
	
	
	// method finds mesh boundary-intersecting lattice cubes (specified by their centroids) with a facet size distribution based resolution
	// method casts rays from xy-plane across domain's entire z-span, storing intersections into a gridded datastructure
	// method then traverses the intersections by marching cubes, storing cube centroid for every cube found to be intersecting the boundary
	// and if addInternals = true, also the cubes fully inside are added
	// the in/out status byteflags of cube corners are stored, if centroids=true, the centroid in/out status is checked and stored
	public void boundaryLattice(FEM1Octree geocTree, int division, boolean centroidStatus, boolean addInternals) {
		
		debugIST_bLatt_time = System.nanoTime();
		FEM1Octant geocRoot = geocTree.root;
		
		double geocSize = geocRoot.xP - geocRoot.xM, bLw = latticeStep = geocSize/(double)division, bLwD2 = bLw*.5, bLwI = 1/bLw;
		// on encountering mesh holes/defects, the salvaging step factor will retest intersections at intervals of 1/100000th of max volume dimension
		double salvF = geocSize / 100000.0;
		subsPerSize = (double)division / geocSize;
		// division signifying subdivided cubes count, their defining edges are one more than no. of lattice cubes
		int dataEndXY = ++division, dataWidth = centroidStatus ? division * 2 - 1 : division;
		// holds grid of facet encounters along z,y,x-ordinates + extra sidegrids for centroids, if requested
		fEncounterGrid = new double[dataWidth][division][];						// intersection points for Z-rays		
		fStatusGrid = new int[dataWidth][division];								// intersections counts	
		int[] status = {0, 0, 0, 0}, fStatusC = null;
		boolean iLast = false;
		double[][] fEncounterC = null;
		
		// generate rays along Z-ordinate
		double baseZ = geocRoot.zM - FEM1Octant.OCT_MARGIN;
		for (int i = 0; i < dataEndXY; i++) {									// i is the Y-ordinate lattice step
			if  (!iLast && i == dataEndXY - 1) iLast = true;					// flag if this is last row of rays
			else {	fEncounterC = fEncounterGrid[dataEndXY + i];				// assign for centroidal encounter if not on last row
					fStatusC = fStatusGrid[dataEndXY + i]; }
			double yStep = geocRoot.yM + (double)i * bLw, xStep = geocRoot.xM;
			
			for (int j = 0; j < division; j++) {								// j is the X-ordinate lattice step
				fEncounterGrid[i][j] = geocRoot.facetEncounters(geocTree, xStep, yStep, baseZ, 0, 0, 0, 2, status);
				if ((status[1] & 1) != 0) {										// if there was a hole (odd no. of intersections), try sampling around area
					if (DEBUG_LEVEL > 1) System.out.println("FEM1.boundaryLattice(): mesh hole at ray: x: "+xStep+"y: "+yStep+", salvaging.");
					fEncounterGrid[i][j] = geocRoot.facetEncounters(geocTree,xStep+salvF,yStep,baseZ,0,0,0,2,status);
					if ((status[1] & 1) != 0) fEncounterGrid[i][j] = geocRoot.facetEncounters(geocTree,xStep-salvF,yStep,baseZ,0,0,0,2,status);
					if ((status[1] & 1) != 0) fEncounterGrid[i][j] = geocRoot.facetEncounters(geocTree,xStep,yStep+salvF,baseZ,0,0,0,2,status);
					if ((status[1] & 1) != 0) fEncounterGrid[i][j] = geocRoot.facetEncounters(geocTree,xStep,yStep-salvF,baseZ,0,0,0,2,status);
					if ((status[1] & 1) != 0) if (DEBUG_LEVEL > 1) System.out.println("FEM1.boundaryLattice(): salvage failed."); }
				fStatusGrid[i][j] = status[1];
				
				if (centroidStatus && !iLast && j < division - 1) {
					fEncounterC[j] = geocRoot.facetEncounters(geocTree, xStep + bLwD2, yStep + bLwD2, baseZ, 0, 0, 0, 2, status);
					if ((status[1] & 1) != 0) {									// if there was a hole in the mesh, try sampling around the area
						if (DEBUG_LEVEL > 1)
							System.out.println("FEM1.boundaryLattice(): mesh hole at ray: x: "+(xStep+bLwD2)+" y: "+(yStep+bLwD2)+", salvaging.");
						fEncounterC[j] = geocRoot.facetEncounters(geocTree,xStep+bLwD2+salvF,yStep+bLwD2,baseZ,0,0,0,2,status);
						if ((status[1]&1)!=0) fEncounterC[j]=geocRoot.facetEncounters(geocTree,xStep+bLwD2-salvF,yStep+bLwD2,baseZ,0,0,0,2,status);
						if ((status[1]&1)!=0) fEncounterC[j]=geocRoot.facetEncounters(geocTree,xStep+bLwD2,yStep+bLwD2+salvF,baseZ,0,0,0,2,status);
						if ((status[1]&1)!=0) fEncounterC[j]=geocRoot.facetEncounters(geocTree,xStep+bLwD2,yStep+bLwD2-salvF,baseZ,0,0,0,2,status);
						if ((status[1]&1)!=0) if (DEBUG_LEVEL > 1) System.out.println("FEM1.boundaryLattice(): salvage failed."); }
					fStatusC[j] = status[1];
				}
				xStep += bLw;
			}
		}
		
		// generate envelopes, allocate envelope array pagewise
		division--;
		double[][] bLatticePg = new double[division][];
		short[][] bLatticeStPg = new short[division][];
		int pageSize = division * division * FEM1.NCOORD, pageSizeSt = division * division;
		bLatticePg[0] = new double[pageSize];									// allocate first page
		bLatticeStPg[0] = new short[pageSizeSt];

		int e = 0, e3 = 0, eP = 0, eSt = 0, eStP = 0;							// envelopes & page counters
		for (int i = 0; i < division; i++) {
			int[] fStatus_i = fStatusGrid[i], fStatus_i1 = fStatusGrid[i + 1];
			double[][] fEncounter_i = fEncounterGrid[i], fEncounter_i1 = fEncounterGrid[i + 1];
			double bLw_iD2 = bLw * (double)i + bLwD2;
			for (int j = 0; j < division; j++) {
	
				int hits_ij = fStatus_i[j], hits_ij1 = fStatus_i[j+1], hits_i1j = fStatus_i1[j], hits_i1j1 = fStatus_i1[j+1], hitsSum;
				int hits_C = centroidStatus ? fStatusGrid[i + division + 1][j] : 0;
				if ((hitsSum = hits_ij + hits_ij1 + hits_i1j + hits_i1j1 + hits_C) == 0) continue;	// do not harvest empty space
				hits_ij *= 3; hits_ij1 *= 3; hits_i1j *= 3; hits_i1j1 *= 3;  hits_C *= 3;
				int stepMM = 2, stepPM = 2, stepMP = 2, stepPP = 2, stepC = 2;
				int envelope = 0, centroid = 0;
				double bLw_jD2 = bLw * (double)j + bLwD2;
				
				double[] fEnc_ij = fEncounter_i[j], fEnc_ij1 = fEncounter_i[j+1], fEnc_i1j = fEncounter_i1[j], fEnc_i1j1 = fEncounter_i1[j+1];
				double[] fEnc_C = centroidStatus ? fEncounterGrid[i + division + 1][j] : null;
				// calculate first intersection encounter to start iteration at
				double en_ij = hits_ij == 0 ? Double.MAX_VALUE : fEnc_ij[stepMM], en_ij1 = hits_ij1 == 0 ? Double.MAX_VALUE : fEnc_ij1[stepPM];
				double en_i1j = hits_i1j == 0 ? Double.MAX_VALUE : fEnc_i1j[stepMP], en_i1j1 = hits_i1j1 == 0 ? Double.MAX_VALUE : fEnc_i1j1[stepPP];
				double en_T1 = en_ij < en_ij1 ? en_ij : en_ij1, en_T2 = en_i1j < en_i1j1 ? en_i1j : en_i1j1, en = en_T1 < en_T2 ? en_T1 : en_T2;
				int kStart = 1 + (int)((en - geocRoot.zM) * bLwI);
				double testPoint = baseZ + bLw*(double)kStart, testPointM = testPoint - bLwD2;
				
				for (int k = kStart; k <= division; k++) {
					// k signifies the iteration through the facet encounters coord arrays
					// note: the while() guarantees to pass cubewalker upto testPoint through microstructures that otherwise can produce "cube blinking"
					// TODO: is it motivated to do something about the features smaller than the lattice step?
					// TODO: maybe statistically accumulate number of microfeature "misses" and decide to redo lattice at increased subdiv. level?
					while (stepMM < hits_ij &&		fEnc_ij[stepMM] < testPoint) {	stepMM+=3; envelope^=16; hitsSum--; }
					while (stepPM < hits_ij1 &&		fEnc_ij1[stepPM] < testPoint) {	stepPM+=3; envelope^=32; hitsSum--; }
					while (stepMP < hits_i1j &&		fEnc_i1j[stepMP] < testPoint) {	stepMP+=3; envelope^=64; hitsSum--; }
					while (stepPP < hits_i1j1 &&	fEnc_i1j1[stepPP] < testPoint) {stepPP+=3; envelope^=128; hitsSum--; }
					// en will hold the next closest boundary encounter point for this cube
					if (!addInternals) {
						if (en > fEnc_ij[stepMM]) en = fEnc_ij[stepMM];		if (en > fEnc_ij1[stepPM]) en = fEnc_ij1[stepPM];
						if (en > fEnc_i1j[stepMP]) en = fEnc_i1j[stepMP];	if (en > fEnc_i1j1[stepPP]) en = fEnc_i1j1[stepPP]; }
					if (centroidStatus)
						while (stepC < hits_C &&	fEnc_C[stepC] < testPointM) { stepC += 3; centroid ^= 256; }

					// every step forward in a facet encounter array means entry or exit into the mesh volume, therefore:
					// unless all bits of envelope are set (fully inside) or cleared (fully outside), the cube is intersecting the boundary
					// OR, if internals are requested, all bits must be set for inclusion (note that centroid can potentially by outside despite that)
					if (addInternals && envelope == 0xFF || envelope != 0 && envelope != 0xFF) {
						e++;
						double[] ePage = bLatticePg[eP];
						if (ePage.length <= e3) { ePage = bLatticePg[++eP] = new double[pageSize]; e3 = 0; }
						if (bLatticeStPg[eStP].length <= eSt) { bLatticeStPg[++eStP] = new short[pageSizeSt]; eSt = 0; }
						ePage[e3++] = geocRoot.xM + bLw_jD2; ePage[e3++] = geocRoot.yM + bLw_iD2; ePage[e3++] = baseZ + bLw*(double)k - bLwD2;
						// store down envelope's inside/outside bitflags + optional centroid in/out flag
						// note: naturally, for an insider, all in/out bits will be set
						bLatticeStPg[eStP][eSt++] = (short)(envelope | centroid);	
						if (hitsSum < 1) break;									// was this the last status change (last predicted cube in row)?
						testPoint += bLw; testPointM += bLw;
					} else {
						if (!addInternals) {									// if internal cubes weren't requested
							k = (int)((en - geocRoot.zM) * bLwI);				// jump past inner volume to next encounter (attention of the k++ above)
							testPoint = baseZ + bLw*k;
							testPointM = testPoint - bLwD2;
						} else { testPoint += bLw; testPointM += bLw; }			// if internal cubes requested, increment one cube length
					}		
					envelope = (envelope >> 4) | (envelope & 0xF0);				// shift intersection status rightwards, which become cube's 4 lower flags
					//testPoint += bLw; 
				}				
			}
		}
		
		// collect pages into a single bLattice[] array
		bLatticeNode = new double[e * NCOORD];
		for (int eP1 = 0, e1 = 0, e3b = 0; eP1 <= eP; eP1++) {
			double[] envelopeP = bLatticePg[eP1];
			for (int e2 = 0; e2 < pageSize;) {
				bLatticeNode[e3b++] = envelopeP[e2++]; bLatticeNode[e3b++] = envelopeP[e2++]; bLatticeNode[e3b++] = envelopeP[e2++];
				if (++e1 >= e) break;
			}
		}
		bLatticeStatus = new short[e];
		for (int eStP1 = 0, e1 = 0, eStb = 0; eStP1 <= eStP; eStP1++) {
			short[] envelopeStP = bLatticeStPg[eStP1];
			if (eStP1 < eStP)
					for (int e2 = 0; e2 < pageSizeSt; e1++) bLatticeStatus[eStb++] = envelopeStP[e2++];		// no need to test end of data until last page
			else	for (int e2 = 0; e2 < pageSizeSt;) { bLatticeStatus[eStb++] = envelopeStP[e2++]; if (++e1 >= e) break; }
		}
		debugIST_bLatt_time = System.nanoTime() - debugIST_bLatt_time;
	}
	
	
	// method completes the octants that have the centroid differing in inside/outside status from the edge vertexes
	// this is made to ensure that only BCC tetrahedra intersect the surface, a leaf centroid boundary int/ext status differing from any given
	// octant vertex int/ext status will enforce 3 octants to be created at the 3 faces that meet on that vertex
	public FEM1Octant[] extendBCCdoNhoodsIST(FEM1Octree latticeTree) {

		int[] culledCount = new int[latticeTree.gradations];
		debugIST_aBOIST_time = System.nanoTime();
		// get ALL octants pertaining to Isosurface Stuffing solution, but we'll only iterate over highest level
		int[] countO = {0};											// gets true count of array, since leafOctantArray[] does worst-case allocation
		FEM1Octant[] octantIST = latticeTree.root.octantArray(latticeTree, latticeTree.maxLevel, countO, true);	
		FEM1Octant[] neighbourhoodIST = new FEM1Octant[octantIST.length * 7];
		int oLast=countO[0], oNonLeafC=0, maxLevel=latticeTree.maxLevel, firstLevel=maxLevel-latticeTree.gradations+1, leafCnt=0;
		int[] gradCount = new int[latticeTree.gradations];

		// iteration creates leaf octants according to Continuation Criterion and generates neighbourhoods both for existing and new leaves
		// neighbourhood backreferences are written and the internal/external status is resolved from the parents
		for (int o = 0, n7 = 0; o < oLast; o++) {
			FEM1Octant oct = octantIST[o];
			if (oct.level < maxLevel) { n7+=7; continue; }
			neighbourhoodIST[n7] = oct;								// first octant of a neighbourhood is neighbourhood's spawning octant
			
			int nFlags =(neighbourhoodIST[++n7]==null ? 1 : 0)  | (neighbourhoodIST[++n7]==null ? 4 : 0)   | (neighbourhoodIST[++n7]==null ? 16 : 0) |
						(neighbourhoodIST[++n7]==null ? 64 : 0) | (neighbourhoodIST[++n7]==null ? 256 : 0) | (neighbourhoodIST[++n7]==null ? 1024 : 0);
			if (nFlags==0) { n7++; continue; }						// a fully assigned octant (mostly interior ones), nothing more to do
			boolean leafProcessed = false; n7 -= 5;
			
			if (o >= oLast) 									// get neighbourhood for the new leaves (do not generate neighbours)
				oct.generateGetFaceNeighbours(latticeTree, neighbourhoodIST, nFlags, n7, true);
			else {												// get neighbourhood for leaves (creation on centroid/corner internal/external criterion)
				boolean centroidIn = oct.centroid_internal();	// flag creation of neighbours that fulfill int/ext criterion, centroid versus corners
				nFlags |= (centroidIn != ((oct.int_ext&0x55)!=0) ? 32 : 0) | (centroidIn != ((oct.int_ext&0xAA)!=0) ? 128 : 0);		// x-, x+ (bits 6, 8)
				nFlags |= (centroidIn != ((oct.int_ext&0x33)!=0) ? 8 : 0)  | (centroidIn != ((oct.int_ext&0xCC)!=0) ? 512 : 0);		// y-, y+ (bits 4, 10)
				nFlags |= (centroidIn != ((oct.int_ext&0x0F)!=0) ? 2 : 0)  | (centroidIn != ((oct.int_ext&0xF0)!=0) ? 2048 : 0);	// z-, z+ (bits 2, 12)
				// the criterions for octant centroid versus each face's vertexes are sent to generateGetFaceNeighbours()
				oct.generateGetFaceNeighbours(latticeTree, neighbourhoodIST, nFlags, n7, true);
				leafProcessed = true; }
			leafCnt++;
			int n7end = n7 + 6, level = oct.level, oEnum = oct.nodes, nbrRevIdx = 6;
			
			while (n7 < n7end) {									// assign found neighbourhoods to pertinent neighbours of octant
				FEM1Octant octN = neighbourhoodIST[n7];
				if (octN == null) { nbrRevIdx--; n7++; continue; }
				
				if (leafProcessed) {								// if that was a leaf processed, process eventual additional created leaves
					if (oLast >= octantIST.length) {				// created leaves will be added, check neighbourhoods&octants arrays overflow
						int newLength = octantIST.length + octantIST.length/4;
						FEM1Octant[] octantISTnew = new FEM1Octant[newLength];
						FEM1Octant[] neighbourhoodISTnew = new FEM1Octant[newLength * 7];
						for (int o1 = 0; o1 < oLast; o1++) octantISTnew[o1] = octantIST[o1];
						for (int o7 = 0, o7end = oLast*7; o7 < o7end; o7++) neighbourhoodISTnew[o7] = neighbourhoodIST[o7];
						octantIST = octantISTnew;
						neighbourhoodIST = neighbourhoodISTnew;
					}
					if (octN.in_progress()) {							// if neighbour had to be created, append to leafArray
						octantIST[oLast] = octN;
						octN.nodes = oLast++;							// continue the octant enumeration						
						octN.int_extFromParent(false);					// resolve int/ext status from parent's int. status (some int. leaves removed earlier)
						if (octN.internal()) octN.int_ext |= FEM1Octant.CENTROID_INTERNAL;		// an internal leaf has internal centroid by default
						// process() will halt on unfinished nodes, so clear in_progress flag for the entire new octant chain
						FEM1Octant octN1 = octN;
						while (octN1.in_progress()) { octN1.finish(); octN1 = octN1.parent; }
					}
				}
				
				// note: IST octant can have a level N-1 neighbour, but NOT a level N+1 neighbour, therefore if levels aren't equal,
				// a back-reference is only done from neighbours that are same level or level N+1
				if (octN.nodes > oEnum)	{									// only interested in assigning future neighbourhoods (past ones already done)
					if (octN.level == level) neighbourhoodIST[octN.nodes * 7 + nbrRevIdx] = oct; 
					// also if same-level back-referencing already was found, ignore N to N-1 level backreferences
					else { /* flag N-1 to N neighbourhood */ }}
				nbrRevIdx--; n7++;
			}
		}
		gradCount[latticeTree.gradations-1] = leafCnt;

		// the iteration creates neighbourhoods on one hand, and shaves off children on each consecutive level on the other
		// IST-constituent octants of lower levels are handled, it checks face&edge neighbourhoods and removes
		// unnecessary children if full gradation requested, face-neighbourhoods are created with backreferences
		int lastLevel = latticeTree.maxLevel-1;
		
		for (int l=lastLevel; l >= firstLevel; l--)
			for (int o = 0, n7 = 0; o < oLast; o++) {
				FEM1Octant oct = octantIST[o];
				if (oct.level != l)	{ n7+=7; continue; }						// skip other levels than l
				int grad = oct.level - firstLevel;
				gradCount[grad]++;
				
				if (l==lastLevel && oct.iterator != 0) oct.nodesX = 0;				// if leaves were generated in next highest level, clear IST-internal status
					
				// if octant is IST-internal with children, try culling children that don't upset a 2:1 balance with neighbours
				else if (oct.level<=latticeTree.maxLevel-2 && oct.octant!= null && (oct.nodesX & 0xFF0000) == 0xFF0000)	
					culledCount[grad+1] += (oct.iterator>>24) - oct.cull2to1IST(latticeTree, false);
				
				neighbourhoodIST[n7] = oct;										// first octant of a neighbourhood is neighbourhood's spawning octant				
				int nFlags =(neighbourhoodIST[++n7]==null ? 1 : 0)  | (neighbourhoodIST[++n7]==null ? 4 : 0)   | (neighbourhoodIST[++n7]==null ? 16 : 0) |
							(neighbourhoodIST[++n7]==null ? 64 : 0) | (neighbourhoodIST[++n7]==null ? 256 : 0) | (neighbourhoodIST[++n7]==null ? 1024 : 0);
				if (nFlags==0) { n7++; continue; }								// a fully assigned octant (mostly interior ones), nothing more to do
				n7 -= 5;

				oct.generateGetFaceNeighbours(latticeTree, neighbourhoodIST, nFlags, n7, true);
				//octantIST[oNonLeafC++] = oct;									// reform octantIST[] to only hold non-"leaf" (=non-maxLevel) octants	
				int n7end = n7 + 6, level = oct.level, oEnum = oct.nodes, nbrRevIdx = 6;
				
				while (n7 < n7end) {											// assign backreferences to octant's neighbours
					FEM1Octant octN = neighbourhoodIST[n7];
					if (octN == null || octN.nodes == -1) { nbrRevIdx--; n7++; continue; }
					
					// note: IST octant can have a level N-1 neighbour, but NOT a level N+1 neighbour, therefore if levels aren't equal,
					// a back-reference is only done for (N -> N) or (N -> N-1)
					if (octN.nodes > oEnum)	{								// only interested in assigning future neighbourhoods (past ones already done)
						if (octN.level == level) neighbourhoodIST[octN.nodes * 7 + nbrRevIdx] = oct; 
						// also if same-level back-referencing already was found, ignore N to N-1 level backreferences
						else { /* flag N-1 to N neighbourhood */ }}
					nbrRevIdx--; n7++;
				}
			}
		
		if (DEBUG_LEVEL > 1) {
			System.out.println("******************* FEM1.extendBCCdoNhoodsIST() *******************");
			System.out.println("Gradation levels:");
			for (int g=0; g < latticeTree.gradations; g++) System.out.println("Gradation " + g + ": " + gradCount[g] + " octants.");
			System.out.println("Gradation culls:");
			for (int g=0; g < latticeTree.gradations; g++) System.out.println("Gradation " + g + ": " + culledCount[g] + " culled octants.");
			System.out.print("\n\n"); }
		gradCount = new int[latticeTree.gradations];
		
		// since a neighbour octant might have been created AFTER a neighbourhood was found, need to recheck the cases
		// of level N octants neighbouring level N-1 octants, to see if a generated child can be made neighbour
		// also any children culled must be purged from neighbourhoods
		for (int o = 0, n7 = 0; o < oLast; o++) {
			FEM1Octant oct = neighbourhoodIST[n7];
			if (oct.nodes == -1) { neighbourhoodIST[n7] = null; n7+=7; continue; }		// skip neighbourhoods of culled octants
			gradCount[oct.level - firstLevel]++;
			if (oct.level == latticeTree.maxLevel) { n7+=7; continue; }					// maxLevel octants can't find any children
			octantIST[oNonLeafC++] = oct;												// recollect non-leaves into octantIST[]
			int n7end = ++n7 + 6, n7b = n7;
			short level = oct.level;
			while (n7 < n7end) {
				FEM1Octant octN = neighbourhoodIST[n7];
				if (octN != null) {
					if (octN.nodes == -1) neighbourhoodIST[n7] = null;
					else if (octN.level < level && octN.octant != null) {
						switch (n7 - n7b) {
						case 0: neighbourhoodIST[n7] = octN.locateCoordinate(oct.xC, oct.yC, oct.zM - FEM1Octant.OCT_MARGIN, level); break;	// z-
						case 1: neighbourhoodIST[n7] = octN.locateCoordinate(oct.xC, oct.yM - FEM1Octant.OCT_MARGIN, oct.zC, level); break;	// y-
						case 2: neighbourhoodIST[n7] = octN.locateCoordinate(oct.xM - FEM1Octant.OCT_MARGIN, oct.yC, oct.zC, level); break; // x-
						case 3: neighbourhoodIST[n7] = octN.locateCoordinate(oct.xP + FEM1Octant.OCT_MARGIN, oct.yC, oct.zC, level); break;	// x+
						case 4: neighbourhoodIST[n7] = octN.locateCoordinate(oct.xC, oct.yP + FEM1Octant.OCT_MARGIN, oct.zC, level); break;	// y+
						case 5: neighbourhoodIST[n7] = octN.locateCoordinate(oct.xC, oct.yC, oct.zP + FEM1Octant.OCT_MARGIN, level); }		// z+
					}}
				n7++;
			}
		}

		if (DEBUG_LEVEL > 2) {
			System.out.println("Gradation levels, post-2:1-culling:");
			for (int g=0; g < latticeTree.gradations; g++) System.out.println("Gradation " + g + ": " + gradCount[g] + " octants.");
			System.out.print("\n\n"); }
		
		int[] edgeOctCheck = new int[oLast];
		// unflag face-neighbour search for the neighbourhood checks
		for (int o = 0; o < oLast; o++) edgeOctCheck[o] = 0x3FFFF - (1<<2|1<<6|1<<8|1<<9|1<<11|1<<15);
		edgeOctFlags = new short[oLast];
		FEM1Octant[] neighbourhood18IST = new FEM1Octant[18];
		
		// get the 18-neighbourhood for every IST octant that isn't a leaf, but only those 12 neighbours that haven's been found yet,
		// and flag any edge that carries children towards current octant
		for (int o = 0; o < oNonLeafC; o++) {
			FEM1Octant oct = octantIST[o];
			if (oct == null) continue;
			int level = oct.level, oEnum = oct.nodes;
			int nCheck = edgeOctCheck[oEnum], n6 = oEnum*7+1;
			oct.getFaceEdgeNeighbours(latticeTree, neighbourhood18IST, 0, nCheck);
			// need to integrate the face neighbours to get their childflags
			neighbourhood18IST[2] = neighbourhoodIST[n6++]; neighbourhood18IST[6] = neighbourhoodIST[n6++];
			neighbourhood18IST[8] = neighbourhoodIST[n6++]; neighbourhood18IST[9] = neighbourhoodIST[n6++];
			neighbourhood18IST[11] = neighbourhoodIST[n6++]; neighbourhood18IST[15] = neighbourhoodIST[n6++];
			short flagsEN = 0;
			
			for (int n18 = 0; n18 < 18 && nCheck != 0; n18++, nCheck>>=1) {
				FEM1Octant octN = neighbourhood18IST[n18];
				if (octN == null || octN.level != level) continue;		// child-bitflagging ONLY relevant for same-level neighbourhoods
				int statusN = octN.status>>3&0xFF;
				switch (n18) {											// flags are stored in bitfields within edgeOctFlag[] for every edge-neighbourhood
				case 0:		if ((statusN&FLG_XPP) != 0) flagsEN |= 1; break;
				case 1:		if ((statusN&FLG_PXP) != 0) flagsEN |= 2; break;
				case 2:		if ((statusN&FLG_XMP) != 0)	flagsEN |= 1; if ((statusN&FLG_MXP) != 0) flagsEN |= 2;
							if ((statusN&FLG_PXP) != 0)	flagsEN |= 4; if ((statusN&FLG_XPP) != 0) flagsEN |= 8; break;
				case 3:		if ((statusN&FLG_MXP) != 0) flagsEN |= 4; break;
				case 4:		if ((statusN&FLG_XMP) != 0) flagsEN |= 8; break;
				case 5:		if ((statusN&FLG_PPX) != 0) flagsEN |= 16; break;
				case 6:		if ((statusN&FLG_XPM) != 0)	flagsEN |= 1; if ((statusN&FLG_MPX) != 0) flagsEN |= 16;
							if ((statusN&FLG_PPX) != 0)	flagsEN |= 32; if ((statusN&FLG_XPP) != 0) flagsEN |= 256; break;
				case 7:		if ((statusN&FLG_MPX) != 0) flagsEN |= 32; break;
				case 8:		if ((statusN&FLG_PXM) != 0)	flagsEN |= 2; if ((statusN&FLG_PMX) != 0) flagsEN |= 16;
							if ((statusN&FLG_PPX) != 0)	flagsEN |= 64; if ((statusN&FLG_PXP) != 0) flagsEN |= 512; break;
				case 9:		if ((statusN&FLG_MXM) != 0)	flagsEN |= 4; if ((statusN&FLG_MMX) != 0) flagsEN |= 32;
							if ((statusN&FLG_MPX) != 0)	flagsEN |= 128; if ((statusN&FLG_MXP) != 0) flagsEN |= 1024; break;
				case 10:	if ((statusN&FLG_PMX) != 0) flagsEN |= 64; break;
				case 11:	if ((statusN&FLG_XMM) != 0)	flagsEN |= 8; if ((statusN&FLG_MMX) != 0) flagsEN |= 64;
							if ((statusN&FLG_PMX) != 0)	flagsEN |= 128; if ((statusN&FLG_XMP) != 0) flagsEN |= 2048; break;
				case 12:	if ((statusN&FLG_MMX) != 0) flagsEN |= 128; break;
				case 13:	if ((statusN&FLG_XPM) != 0) flagsEN |= 256; break;
				case 14:	if ((statusN&FLG_PXM) != 0) flagsEN |= 512; break;
				case 15:	if ((statusN&FLG_XMM) != 0)	flagsEN |= 256; if ((statusN&FLG_MXM) != 0) flagsEN |= 512;
							if ((statusN&FLG_PXM) != 0)	flagsEN |= 1024; if ((statusN&FLG_XPM) != 0) flagsEN |= 2048; break;
				case 16: 	if ((statusN&FLG_MXM) != 0) flagsEN |= 1024; break;
				case 17: 	if ((statusN&FLG_XMM) != 0) flagsEN |= 2048;
				}
			}
			edgeOctFlags[oEnum] = flagsEN;
		}
		octantIST = null;
		debugIST_aBOIST_time = System.nanoTime() - debugIST_aBOIST_time;
		return neighbourhoodIST;										// return neighbourhoods collection
	}
	
	
	// accumulates IST elements supported by every node in individual arrays
	void elementsNodeSupportIST() {
		elementNhoodIST = new int[nodesIST][];

		for (int e = 0, e6 = 0; e < elementsIST; e++, e6+=2) {
			for (int n = 0; n < 4; n++) {
				int node = elementIST[e6++] & 0x3FFFFFFF;	
				int[] support = elementNhoodIST[node];
				if (support == null) {	support = elementNhoodIST[node] = new int[9];			// start allocating for 8 elements
										support[0] = 1; support[1] = e;	continue; }					// new support array started, insert it's first element
				
				int elemsS = support[0];
				if (++elemsS >= support.length-1) {													// check if support array needs resizing
					int elemsS2 = 1 + elemsS + (elemsS >> 1);										// resize 1.5 times
					int[] support2 = new int[elemsS2];
					if (elementSupportMaxArrayL < elemsS2) elementSupportMaxArrayL = elemsS2;
					for (int e2 = 0; e2 < elemsS; e2++) support2[e2] = support[e2];
					support = elementNhoodIST[node] = support2; }
				
				// lists must be sorted for tetraNeighbours() method to work correctly, therefore, insert in ascending order
				if (support[--elemsS] < e) { support[0] = ++elemsS; support[elemsS] = e; }			// if this index is higher than LAST one, just append
				else {	for (int s = 1; s <= elemsS; s++) {											// otherwise insert element index in sorted order
							if (e == support[s]) break;												// if element exists, stop
							if (e < support[s]) {													// is THIS index larger than inserted one?
								for (int s1 = elemsS + 1; s1 > s; s1++) support[s1] = support[s];	// yes, push it & followers one slot forward
								support[s] = e; support[0] = ++elemsS; break; }						// insert index
					}
				}
			}
		}
	}
	
	void nodeSupportAddElementIST(int node, int e) {
		int[] support = elementNhoodIST[node];
		if (support == null) {
			support = elementNhoodIST[node] = new int[elementSupportMaxArrayL];				// start allocating for 8 elements
			support[0] = 1; support[1] = e;	return; }										// new support array started, insert it's first element
		
		int elemsS = support[0];
		if (++elemsS >= support.length-1) {													// check if support array needs resizing
			int elemsS2 = 1 + elemsS + (elemsS >> 1);										// resize 1.5 times
			int[] support2 = new int[elemsS2];
			if (elementSupportMaxArrayL < elemsS2) elementSupportMaxArrayL = elemsS2;
			for (int e2 = 0; e2 < elemsS; e2++) support2[e2] = support[e2];
			support = elementNhoodIST[node] = support2; }
		
		// list must be sorted for tetraNeighbours() method to work correctly, therefore, insert in ascending order
		if (support[--elemsS] < e) { support[0] = ++elemsS; support[elemsS] = e; }			// if this index is higher than LAST one, just append
		else {	for (int s = 1; s <= elemsS; s++) {											// otherwise insert element index in sorted order
					if (e == support[s]) break;												// if element exists, stop
					if (e < support[s]) {													// is THIS index larger than inserted one?
						for (int s1 = elemsS + 1; s1 > s; s1++) support[s1] = support[s];	// yes, push it & followers one slot forward
						support[s] = e; support[0] = ++elemsS; break; }						// insert index
			}
		}		
	}

	

	// the edge status bitsets are merged into a bitpattern that will be matched against all possible permutations (x12) of that pattern
	final static int EC_P_P=9, EC_P_M=1, EC_M_P=8, EC_P_0=17, EC_0_P=10, EC_M_M=0, EC_M_0=16, EC_0_M=2, EC_0_0=18, EC_PxM=5, EC_MxP=12, EC_Mx0=20, EC_0xM=6;
	final static int EC_C01=1<<2, EC_C02=1<<7, EC_C03=1<<12, EC_C12=1<<17, EC_C13=1<<22, EC_C23=1<<27;

	// the stencil permutation codes are pregenerated to test for a specific stencil from an arbitrary starting node and aspect,
	// since a tetrahedron can potentially be generated from any starting node and in any allowed positive-volume index permutation
	// note: this code was autogenerated by stencilPermutations() & stencilPermToJavaCode() methods
	final static int[] stencilPermCode = {
			// the 3 stencils of arbitrary orientation & reflection
			EC_P_0|EC_PxM<<5|EC_P_0<<10|EC_0_M<<15|EC_0_0<<20|EC_M_0<<25, 0|1<<4|2<<8|3<<12 | 1<<16|8<<20|8<<24|8<<28,	// 0, stencil 1
			EC_PxM|EC_P_0<<5|EC_P_0<<10|EC_M_0<<15|EC_M_0<<20|EC_0_0<<25, 0|2<<4|3<<8|1<<12 | 0<<16|8<<20|8<<24|8<<28,
			EC_P_0|EC_P_0<<5|EC_PxM<<10|EC_0_0<<15|EC_0_M<<20|EC_0_M<<25, 0|3<<4|1<<8|2<<12 | 2<<16|8<<20|8<<24|8<<28,
			EC_0_0|EC_0_M<<5|EC_0_P<<10|EC_0_M<<15|EC_0_P<<20|EC_MxP<<25, 1|3<<4|2<<8|0<<12 | 5<<16|8<<20|8<<24|8<<28,
			EC_0_M|EC_0_P<<5|EC_0_0<<10|EC_MxP<<15|EC_M_0<<20|EC_P_0<<25, 1|2<<4|0<<8|3<<12 | 3<<16|8<<20|8<<24|8<<28,
			EC_0_P|EC_0_0<<5|EC_0_M<<10|EC_P_0<<15|EC_PxM<<20|EC_0_M<<25, 1|0<<4|3<<8|2<<12 | 4<<16|8<<20|8<<24|8<<28,
			EC_M_0|EC_M_0<<5|EC_MxP<<10|EC_0_0<<15|EC_0_P<<20|EC_0_P<<25, 2|1<<4|3<<8|0<<12 | 2<<16|8<<20|8<<24|8<<28,
			EC_MxP|EC_M_0<<5|EC_M_0<<10|EC_P_0<<15|EC_P_0<<20|EC_0_0<<25, 2|0<<4|1<<8|3<<12 | 0<<16|8<<20|8<<24|8<<28,
			EC_M_0|EC_MxP<<5|EC_M_0<<10|EC_0_P<<15|EC_0_0<<20|EC_P_0<<25, 2|3<<4|0<<8|1<<12 | 1<<16|8<<20|8<<24|8<<28,
			EC_0_M|EC_0_0<<5|EC_0_P<<10|EC_M_0<<15|EC_MxP<<20|EC_0_P<<25, 3|2<<4|1<<8|0<<12 | 4<<16|8<<20|8<<24|8<<28,
			EC_0_0|EC_0_P<<5|EC_0_M<<10|EC_0_P<<15|EC_0_M<<20|EC_PxM<<25, 3|1<<4|0<<8|2<<12 | 5<<16|8<<20|8<<24|8<<28,
			EC_0_P|EC_0_M<<5|EC_0_0<<10|EC_PxM<<15|EC_P_0<<20|EC_M_0<<25, 3|0<<4|2<<8|1<<12 | 3<<16|8<<20|8<<24|8<<28,
			EC_PxM|EC_PxM<<5|EC_P_0<<10|EC_M_M<<15|EC_M_0<<20|EC_M_0<<25, 0|1<<4|2<<8|3<<12 | 0<<16|1<<20|8<<24|8<<28,	// 12, stencil 2
			EC_PxM|EC_P_0<<5|EC_PxM<<10|EC_M_0<<15|EC_M_M<<20|EC_0_M<<25, 0|2<<4|3<<8|1<<12 | 2<<16|0<<20|8<<24|8<<28,
			EC_P_0|EC_PxM<<5|EC_PxM<<10|EC_0_M<<15|EC_0_M<<20|EC_M_M<<25, 0|3<<4|1<<8|2<<12 | 1<<16|2<<20|8<<24|8<<28,
			EC_M_0|EC_M_M<<5|EC_MxP<<10|EC_0_M<<15|EC_0_P<<20|EC_MxP<<25, 1|3<<4|2<<8|0<<12 | 2<<16|5<<20|8<<24|8<<28,
			EC_M_M|EC_MxP<<5|EC_M_0<<10|EC_MxP<<15|EC_M_0<<20|EC_P_0<<25, 1|2<<4|0<<8|3<<12 | 1<<16|3<<20|8<<24|8<<28,
			EC_MxP|EC_M_0<<5|EC_M_M<<10|EC_P_0<<15|EC_PxM<<20|EC_0_M<<25, 1|0<<4|3<<8|2<<12 | 0<<16|4<<20|8<<24|8<<28,
			EC_M_M|EC_M_0<<5|EC_MxP<<10|EC_M_0<<15|EC_MxP<<20|EC_0_P<<25, 2|1<<4|3<<8|0<<12 | 4<<16|2<<20|8<<24|8<<28,
			EC_MxP|EC_M_M<<5|EC_M_0<<10|EC_PxM<<15|EC_P_0<<20|EC_M_0<<25, 2|0<<4|1<<8|3<<12 | 3<<16|0<<20|8<<24|8<<28,
			EC_M_0|EC_MxP<<5|EC_M_M<<10|EC_0_P<<15|EC_0_M<<20|EC_PxM<<25, 2|3<<4|0<<8|1<<12 | 5<<16|1<<20|8<<24|8<<28,
			EC_0_M|EC_0_M<<5|EC_0_P<<10|EC_M_M<<15|EC_MxP<<20|EC_MxP<<25, 3|2<<4|1<<8|0<<12 | 5<<16|4<<20|8<<24|8<<28,
			EC_0_M|EC_0_P<<5|EC_0_M<<10|EC_MxP<<15|EC_M_M<<20|EC_PxM<<25, 3|1<<4|0<<8|2<<12 | 3<<16|5<<20|8<<24|8<<28,
			EC_0_P|EC_0_M<<5|EC_0_M<<10|EC_PxM<<15|EC_PxM<<20|EC_M_M<<25, 3|0<<4|2<<8|1<<12 | 4<<16|3<<20|8<<24|8<<28,
			EC_PxM|EC_PxM<<5|EC_PxM<<10|EC_M_M<<15|EC_M_M<<20|EC_M_M<<25, 0|1<<4|2<<8|3<<12 | 0<<16|1<<20|2<<24|8<<28,	// 24, stencil 3
			EC_PxM|EC_PxM<<5|EC_PxM<<10|EC_M_M<<15|EC_M_M<<20|EC_M_M<<25, 0|2<<4|3<<8|1<<12 | 2<<16|0<<20|1<<24|8<<28,
			EC_PxM|EC_PxM<<5|EC_PxM<<10|EC_M_M<<15|EC_M_M<<20|EC_M_M<<25, 0|3<<4|1<<8|2<<12 | 1<<16|2<<20|0<<24|8<<28,
			EC_M_M|EC_M_M<<5|EC_MxP<<10|EC_M_M<<15|EC_MxP<<20|EC_MxP<<25, 1|3<<4|2<<8|0<<12 | 2<<16|5<<20|4<<24|8<<28,
			EC_M_M|EC_MxP<<5|EC_M_M<<10|EC_MxP<<15|EC_M_M<<20|EC_PxM<<25, 1|2<<4|0<<8|3<<12 | 1<<16|3<<20|5<<24|8<<28,
			EC_MxP|EC_M_M<<5|EC_M_M<<10|EC_PxM<<15|EC_PxM<<20|EC_M_M<<25, 1|0<<4|3<<8|2<<12 | 0<<16|4<<20|3<<24|8<<28,
			EC_M_M|EC_M_M<<5|EC_MxP<<10|EC_M_M<<15|EC_MxP<<20|EC_MxP<<25, 2|1<<4|3<<8|0<<12 | 4<<16|2<<20|5<<24|8<<28,
			EC_MxP|EC_M_M<<5|EC_M_M<<10|EC_PxM<<15|EC_PxM<<20|EC_M_M<<25, 2|0<<4|1<<8|3<<12 | 3<<16|0<<20|4<<24|8<<28,
			EC_M_M|EC_MxP<<5|EC_M_M<<10|EC_MxP<<15|EC_M_M<<20|EC_PxM<<25, 2|3<<4|0<<8|1<<12 | 5<<16|1<<20|3<<24|8<<28,
			EC_M_M|EC_M_M<<5|EC_MxP<<10|EC_M_M<<15|EC_MxP<<20|EC_MxP<<25, 3|2<<4|1<<8|0<<12 | 5<<16|4<<20|2<<24|8<<28,
			EC_M_M|EC_MxP<<5|EC_M_M<<10|EC_MxP<<15|EC_M_M<<20|EC_PxM<<25, 3|1<<4|0<<8|2<<12 | 3<<16|5<<20|1<<24|8<<28,
			EC_MxP|EC_M_M<<5|EC_M_M<<10|EC_PxM<<15|EC_PxM<<20|EC_M_M<<25, 3|0<<4|2<<8|1<<12 | 4<<16|3<<20|0<<24|8<<28,
			// 2 stencils limited to LONG/SHORT edge orientation & mirror permutation along assymetric axis
			EC_M_0|EC_MxP<<5|EC_MxP<<10|EC_0_P<<15|EC_0_P<<20|EC_P_P<<25, 0|1<<4|2<<8|3<<12 | 1<<16|2<<20|8<<24|8<<28,	// 36, stencil 4
			EC_P_P|EC_P_0<<5|EC_PxM<<10|EC_P_0<<15|EC_PxM<<20|EC_0_M<<25, 3|2<<4|1<<8|0<<12 | 4<<16|2<<20|8<<24|8<<28,
			EC_P_P|EC_PxM<<5|EC_P_0<<10|EC_PxM<<15|EC_P_0<<20|EC_M_0<<25, 2|3<<4|0<<8|1<<12 | 1<<16|3<<20|8<<24|8<<28,
			EC_0_M|EC_0_P<<5|EC_0_P<<10|EC_MxP<<15|EC_MxP<<20|EC_P_P<<25, 1|0<<4|3<<8|2<<12 | 4<<16|3<<20|8<<24|8<<28,
			EC_MxP|EC_M_0<<5|EC_MxP<<10|EC_P_0<<15|EC_P_P<<20|EC_0_P<<25, 0|1<<4|2<<8|3<<12 | 0<<16|2<<20|8<<24|8<<28,	// 40, stencil 4 mirrored
			EC_P_0|EC_P_P<<5|EC_PxM<<10|EC_0_P<<15|EC_0_M<<20|EC_PxM<<25, 3|2<<4|1<<8|0<<12 | 5<<16|2<<20|8<<24|8<<28,
			EC_0_P|EC_0_M<<5|EC_0_P<<10|EC_PxM<<15|EC_P_P<<20|EC_MxP<<25, 2|3<<4|0<<8|1<<12 | 5<<16|3<<20|8<<24|8<<28,
			EC_PxM|EC_P_P<<5|EC_P_0<<10|EC_MxP<<15|EC_M_0<<20|EC_P_0<<25, 1|0<<4|3<<8|2<<12 | 0<<16|3<<20|8<<24|8<<28,
			EC_PxM|EC_P_P<<5|EC_PxM<<10|EC_MxP<<15|EC_M_M<<20|EC_PxM<<25, 0|1<<4|2<<8|3<<12 | 0<<16|2<<20|3<<24|5<<28,	// 44, stencil 5
			EC_MxP|EC_M_M<<5|EC_MxP<<10|EC_PxM<<15|EC_P_P<<20|EC_MxP<<25, 3|2<<4|1<<8|0<<12 | 5<<16|2<<20|3<<24|0<<28,
			EC_PxM|EC_P_P<<5|EC_PxM<<10|EC_MxP<<15|EC_M_M<<20|EC_PxM<<25, 2|3<<4|0<<8|1<<12 | 5<<16|3<<20|2<<24|0<<28,
			EC_MxP|EC_M_M<<5|EC_MxP<<10|EC_PxM<<15|EC_P_P<<20|EC_MxP<<25, 1|0<<4|3<<8|2<<12 | 0<<16|3<<20|2<<24|5<<28,
			EC_P_P|EC_PxM<<5|EC_PxM<<10|EC_PxM<<15|EC_PxM<<20|EC_M_M<<25, 0|1<<4|2<<8|3<<12 | 1<<16|2<<20|3<<24|4<<28,	// 48, stencil 5 mirrored
			EC_M_M|EC_MxP<<5|EC_MxP<<10|EC_MxP<<15|EC_MxP<<20|EC_P_P<<25, 3|2<<4|1<<8|0<<12 | 4<<16|2<<20|3<<24|1<<28,
			EC_M_M|EC_MxP<<5|EC_MxP<<10|EC_MxP<<15|EC_MxP<<20|EC_P_P<<25, 2|3<<4|0<<8|1<<12 | 1<<16|3<<20|2<<24|4<<28,
			EC_P_P|EC_PxM<<5|EC_PxM<<10|EC_PxM<<15|EC_PxM<<20|EC_M_M<<25, 1|0<<4|3<<8|2<<12 | 4<<16|3<<20|2<<24|1<<28,
			// 3 stencils limited to LONG/SHORT edge orientation & to Parity Rule (matching COMPOUND tet. edge orientation neighbour elements)
			EC_MxP|EC_MxP<<5|EC_M_0<<10|EC_P_P<<15|EC_P_0<<20|EC_P_0<<25, 0|1<<4|2<<8|3<<12 | 0<<16|1<<20|8<<24|8<<28,	// 52, stencil 6
			EC_0_P|EC_0_P<<5|EC_0_M<<10|EC_P_P<<15|EC_PxM<<20|EC_PxM<<25, 3|2<<4|1<<8|0<<12 | 5<<16|4<<20|8<<24|8<<28,
			EC_P_0|EC_PxM<<5|EC_P_P<<10|EC_0_M<<15|EC_0_P<<20|EC_MxP<<25, 2|3<<4|0<<8|1<<12 | 5<<16|1<<20|8<<24|8<<28,
			EC_PxM|EC_P_0<<5|EC_P_P<<10|EC_M_0<<15|EC_MxP<<20|EC_0_P<<25, 1|0<<4|3<<8|2<<12 | 0<<16|4<<20|8<<24|8<<28,
			EC_MxP|EC_MxP<<5|EC_MxP<<10|EC_P_P<<15|EC_P_P<<20|EC_P_P<<25, 0|1<<4|2<<8|3<<12 | 0<<16|1<<20|2<<24|8<<28,	// 56, stencil 7
			EC_P_P|EC_P_P<<5|EC_PxM<<10|EC_P_P<<15|EC_PxM<<20|EC_PxM<<25, 3|2<<4|1<<8|0<<12 | 5<<16|4<<20|2<<24|8<<28,
			EC_P_P|EC_PxM<<5|EC_P_P<<10|EC_PxM<<15|EC_P_P<<20|EC_MxP<<25, 2|3<<4|0<<8|1<<12 | 5<<16|1<<20|3<<24|8<<28,
			EC_PxM|EC_P_P<<5|EC_P_P<<10|EC_MxP<<15|EC_MxP<<20|EC_P_P<<25, 1|0<<4|3<<8|2<<12 | 0<<16|4<<20|3<<24|8<<28,
			EC_MxP|EC_MxP<<5|EC_M_M<<10|EC_P_P<<15|EC_PxM<<20|EC_PxM<<25, 0|1<<4|2<<8|3<<12 | 0<<16|1<<20|4<<24|5<<28,	// 60, stencil 8
			EC_MxP|EC_MxP<<5|EC_M_M<<10|EC_P_P<<15|EC_PxM<<20|EC_PxM<<25, 3|2<<4|1<<8|0<<12 | 5<<16|4<<20|1<<24|0<<28,
			EC_PxM|EC_PxM<<5|EC_P_P<<10|EC_M_M<<15|EC_MxP<<20|EC_MxP<<25, 2|3<<4|0<<8|1<<12 | 5<<16|1<<20|4<<24|0<<28,
			EC_PxM|EC_PxM<<5|EC_P_P<<10|EC_M_M<<15|EC_MxP<<20|EC_MxP<<25, 1|0<<4|3<<8|2<<12 | 0<<16|4<<20|1<<24|5<<28 };
	
	final byte[] stencilPermSwitch = {	1,1,1,1,1,1,1,1,1,1,1,1,	2,2,2,2,2,2,2,2,2,2,2,2,	3,3,3,3,3,3,3,3,3,3,3,3,
										4,4,4,4,4,4,4,4,	5,5,5,5,50,50,50,50,	6,6,6,6,	7,7,7,7,	8,8,8,8};
	
	// octant codes & flag masks for checking children within parents
	final static short OCT_MMM=0, OCT_PMM=1, OCT_MPM=2, OCT_PPM=3, OCT_MMP=4, OCT_PMP=5, OCT_MPP=6, OCT_PPP=7;
	static final byte FLG_MXX = (byte)0x55, FLG_PXX = (byte)0xAA, FLG_XMX = (byte)0x33, FLG_XPX = (byte)0xCC, FLG_XXM = (byte)0x0F, FLG_XXP = (byte)0xF0;
	final static byte FLG_XMM = FLG_XMX&FLG_XXM, FLG_MXM = FLG_MXX&FLG_XXM, FLG_PXM = FLG_PXX&FLG_XXM, FLG_XPM = FLG_XPX&FLG_XXM;
	final static byte FLG_MMX = FLG_MXX&FLG_XMX, FLG_PMX = FLG_PXX&FLG_XMX, FLG_MPX = FLG_MXX&FLG_XPX, FLG_PPX = FLG_PXX&FLG_XPX;
	final static byte FLG_XMP = FLG_XMX&FLG_XXP, FLG_MXP = FLG_MXX&FLG_XXP, FLG_PXP = FLG_PXX&FLG_XXP, FLG_XPP = FLG_XPX&FLG_XXP;	
	
	// program codes specifying which corner/center/halfway coordinate within octant to make a node out of
	final static byte Cd_MMM=0, Cd_PMM=1, Cd_MPM=2, Cd_PPM=3, Cd_MMP=4, Cd_PMP=5, Cd_MPP=6, Cd_PPP=7, Cd_CCC=8;
	final static byte Cd_CMM=9, Cd_CPM=10, Cd_CMP=11, Cd_CPP=12, Cd_MCM=13, Cd_PCM=14, Cd_MCP=15, Cd_PCP=16, Cd_MMC=17, Cd_PMC=18, Cd_MPC=19, Cd_PPC=20;
	final static byte Cd_MCC=21, Cd_PCC=22, Cd_CMC=23, Cd_CPC=24, Cd_CCM=25, Cd_CCP=26, Cd_CCCn=27;
	
	// defines the edge vector direction in space (limited to the possible movements withinin the BCC grid), forms a bitfield with only 1 bit allowed 
	final static int	EV_X=2, EV_Y=2, EV_Z=4, EV_XYZ=8, EV_xYZ=16, EV_xyZ=32, EV_XyZ=64,
						EV_x=128+2, EV_y=128+2, EV_z=128+4, EV_XYz=128+8, EV_xYz=128+16, EV_xyz=128+32, EV_Xyz=128+64;
	// array specifies the edge vectors for every possible element shape (enumerated along the shapes in faceSeq_IST[])
//	final static int[] elemEdgeVec = {
//			EV_xYz|EV_xYZ<<8|EV_x<<16|EV_Z<<24|EV_xyZ<<32|EV_xyz<<40, EV_xyz|EV_xyZ<<8|EV_y<<16|EV_Z<<24|EV_XyZ<<16|EV_Xyz<<24};
	
	// these codes specify the type of edge a tetrahedron will have between each node, ex: Long(0-1), Short(1-2), Short(2-3), Long(3-0)
	// "long" edges are those running along the BCC ordinates, "short" edges are the diagonal ones
	// "Te_M" flags that a capped tetrahedron of that position should have its consituent stencil tetrahedra mirrored
	// TODO: unclear which status to apply to which edge for all other types of tetrahedra besides the obvious BCC cases
	final static byte Te_SSSSSS=0, Te_LSSSSS=2, Te_SLSSSS=2, Te_SSLLSS = 12, Te_SSSLSS = 8, Te_SSSSLS=16, Te_SSSSSL=32, Te_M = (byte)0x80;
	
	// tetra slots manipulation codes: Sl_ChgA: change all, Sl_FtoL: move 0th to 3rd change 0th, Sl_Hld2: hold 2nd change rest, Sl_Ch01: change 0th&1st,
	// Sl_Mv01: move 0th to 1st change 0th,  Sl_Chg0/1/2/3: put coordinate in slot 0/1/2/3, Sl_Ch03: change nodes 0/3
	final static int CLR28=0xEFFFFFFF, SET28=0x10000000, CLR29=0xDFFFFFFF, SET29=0x20000000;
	final static int CLR30=0xBFFFFFFF, SET30=0x40000000, CLR31=0x7FFFFFFF, SET31=0x80000000;
	final static byte Sl_ChgA=0, Sl_FtoL=1, Sl_Hld2=2, Sl_Ch01=3, Sl_Mv01=4, Sl_Ch03=5, Sl_ShMB=6, Sl_Chg0=8, Sl_Chg1=9, Sl_Chg2=10, Sl_Chg3=11, Sl_Ch02=12;
	final static byte slotcodeSumIST[] = { 4, 1, 3, 2, 1, 2, 2, 0, 1, 1, 1, 1, 2, 0, 0, 0 };
	final static int COMPACT_NEXT_ChgA = SET31, COMPACT_internal = SET30;
	
	// following is a metaprogram for investigating and taking decisions on tetrahedra creation from an octant's face/edge neighbourhood perspective
	// the metaprogram will mark tetrahedral coordinates for a compact bitcoded (slot codes) stream where duplicate coordinates are rehashed
	// note: aspects are enumerated ANTICLOCKWISE as viewed from INSIDE from central point (Cd_CCC) of octant (Cd_CCCn is neighbour's center)
	final static byte[] faceSeq_IST = {
		0, 2, 1, 3, 5, 4,
		// [Program 1] code for NEIGHBOUR 0 enumerated face of same size, the layout description (example):
		//   edge nbr, slotcode, Long/Short edges,     (picked 4x octant coords),          (alternative program B if checked corner nbr is of higher level),
		//          0,  Sl_ChgA,          Te_SLSL,  Cd_PMM,Cd_CCCn,Cd_CCC,Cd_MMM, Sl_ChgA,Cd_PCM,Cd_CCCn,Cd_CCC,Cd_MMM,Sl_FtoL,Cd_PMM,Cd_CCCn,Cd_CCC,Cd_PCM,
		FLG_XXP,	// the mask for child-bits of neighbour
/*7*/	0, Sl_ChgA,     Te_SSLLSS, Cd_PMM,Cd_CCCn,Cd_CCC,Cd_MMM, Sl_ChgA,Te_SSSLSS,Cd_CMM,Cd_CCCn,Cd_CCC,Cd_MMM,Sl_FtoL,Te_SSSLSS,Cd_PMM,Cd_CCCn,Cd_CCC,Cd_CMM,
/*26*/	2, Sl_FtoL|Te_M,Te_SSLLSS, Cd_PPM,Cd_CCCn,Cd_CCC,Cd_PMM, Sl_FtoL,Te_SSSLSS,Cd_PCM,Cd_CCCn,Cd_CCC,Cd_PMM,Sl_FtoL,Te_SSSLSS,Cd_PPM,Cd_CCCn,Cd_CCC,Cd_PCM,
/*45*/	3, Sl_FtoL,     Te_SSLLSS, Cd_MPM,Cd_CCCn,Cd_CCC,Cd_PPM, Sl_FtoL,Te_SSSLSS,Cd_CPM,Cd_CCCn,Cd_CCC,Cd_PPM,Sl_FtoL,Te_SSSLSS,Cd_MPM,Cd_CCCn,Cd_CCC,Cd_CPM,
/*64*/	1, Sl_FtoL|Te_M,Te_SSLLSS, Cd_MMM,Cd_CCCn,Cd_CCC,Cd_MPM, Sl_FtoL,Te_SSSLSS,Cd_MCM,Cd_CCCn,Cd_CCC,Cd_MPM,Sl_FtoL,Te_SSSLSS,Cd_MMM,Cd_CCCn,Cd_CCC,Cd_MCM,
		// [Program 2] code for NEIGHBOUR 0 enumerated face with facing children
/*83*/	0, Sl_ChgA,		Te_SSSSSS, Cd_PMM,Cd_CCM,Cd_CCC,Cd_MMM, Sl_ChgA,Te_SSSSSL,Cd_CMM,Cd_CCM,Cd_CCC,Cd_MMM,Sl_FtoL,Te_SLSSSS,Cd_PMM,Cd_CCM,Cd_CCC,Cd_CMM,
/*102*/	2, Sl_FtoL|Te_M,Te_SSSSSS, Cd_PPM,Cd_CCM,Cd_CCC,Cd_PMM, Sl_FtoL,Te_SSSSSL,Cd_PCM,Cd_CCM,Cd_CCC,Cd_PMM,Sl_FtoL,Te_SLSSSS,Cd_PPM,Cd_CCM,Cd_CCC,Cd_PCM,
/*121*/	3, Sl_FtoL,		Te_SSSSSS, Cd_MPM,Cd_CCM,Cd_CCC,Cd_PPM, Sl_FtoL,Te_SSSSSL,Cd_CPM,Cd_CCM,Cd_CCC,Cd_PPM,Sl_FtoL,Te_SLSSSS,Cd_MPM,Cd_CCM,Cd_CCC,Cd_CPM,
/*140*/	1, Sl_FtoL|Te_M,Te_SSSSSS, Cd_MMM,Cd_CCM,Cd_CCC,Cd_MPM, Sl_FtoL,Te_SSSSSL,Cd_MCM,Cd_CCM,Cd_CCC,Cd_MPM,Sl_FtoL,Te_SLSSSS,Cd_MMM,Cd_CCM,Cd_CCC,Cd_MCM,
		// [Program 3] code for NEIGHBOUR 0 enumerated face of double size, the layout:
		// (first program with diagonal slicing for octants like MMM,MMP,PPM,PPP), (second program with diagonal for octants like PMM,PMP,MPM,MPP)
/*159*/	Sl_Hld2,Te_LSSSSS,Cd_PPM,Cd_MMM,Cd_CCC,Cd_PMM,Sl_FtoL,Te_SSSSLS,Cd_MPM,Cd_MMM,Cd_CCC,Cd_PPM,
		Sl_Hld2,Te_LSSSSS,Cd_PMM,Cd_MPM,Cd_CCC,Cd_MMM,Sl_FtoL,Te_SSSSLS,Cd_PPM,Cd_MPM,Cd_CCC,Cd_PMM,

		// [Program 1] code for NEIGHBOUR 2 enumerated face of same size
		FLG_PXX,
/*184*/	6, Sl_Ch01|Te_M,Te_SSLLSS, Cd_MPP,Cd_CCCn,Cd_CCC,Cd_MPM, Sl_Hld2,Te_SSSLSS,Cd_MPC,Cd_CCCn,Cd_CCC,Cd_MPM,Sl_FtoL,Te_SSSLSS,Cd_MPP,Cd_CCCn,Cd_CCC,Cd_MPC,
/*203*/	9, Sl_FtoL,		Te_SSLLSS, Cd_MMP,Cd_CCCn,Cd_CCC,Cd_MPP, Sl_FtoL,Te_SSSLSS,Cd_MCP,Cd_CCCn,Cd_CCC,Cd_MPP,Sl_FtoL,Te_SSSLSS,Cd_MMP,Cd_CCCn,Cd_CCC,Cd_MCP,
/*222*/	4, Sl_FtoL|Te_M,Te_SSLLSS, Cd_MMM,Cd_CCCn,Cd_CCC,Cd_MMP, Sl_FtoL,Te_SSSLSS,Cd_MMC,Cd_CCCn,Cd_CCC,Cd_MMP,Sl_FtoL,Te_SSSLSS,Cd_MMM,Cd_CCCn,Cd_CCC,Cd_MMC,
/*241*/	1, Sl_FtoL,		Te_SSLLSS, Cd_MPM,Cd_CCCn,Cd_CCC,Cd_MMM, Sl_FtoL,Te_SSSLSS,Cd_MCM,Cd_CCCn,Cd_CCC,Cd_MMM,Sl_FtoL,Te_SSSLSS,Cd_MPM,Cd_CCCn,Cd_CCC,Cd_MCM,
		// [Program 2] code for NEIGHBOUR 2 enumerated face with facing children
/*260*/	6, Sl_Ch01|Te_M,Te_SSSSSS, Cd_MPP,Cd_MCC,Cd_CCC,Cd_MPM, Sl_Hld2,Te_SSSSSL,Cd_MPC,Cd_MCC,Cd_CCC,Cd_MPM,Sl_FtoL,Te_SLSSSS,Cd_MPP,Cd_MCC,Cd_CCC,Cd_MPC,
/*279*/	9, Sl_FtoL,		Te_SSSSSS, Cd_MMP,Cd_MCC,Cd_CCC,Cd_MPP, Sl_FtoL,Te_SSSSSL,Cd_MCP,Cd_MCC,Cd_CCC,Cd_MPP,Sl_FtoL,Te_SLSSSS,Cd_MMP,Cd_MCC,Cd_CCC,Cd_MCP,
/*298*/	4, Sl_FtoL|Te_M,Te_SSSSSS, Cd_MMM,Cd_MCC,Cd_CCC,Cd_MMP, Sl_FtoL,Te_SSSSSL,Cd_MMC,Cd_MCC,Cd_CCC,Cd_MMP,Sl_FtoL,Te_SLSSSS,Cd_MMM,Cd_MCC,Cd_CCC,Cd_MMC,
/*317*/	1, Sl_FtoL,		Te_SSSSSS, Cd_MPM,Cd_MCC,Cd_CCC,Cd_MMM, Sl_FtoL,Te_SSSSSL,Cd_MCM,Cd_MCC,Cd_CCC,Cd_MMM,Sl_FtoL,Te_SLSSSS,Cd_MPM,Cd_MCC,Cd_CCC,Cd_MCM,
		// [Program 3] code for NEIGHBOUR 2 enumerated face of double size
/*336*/	Sl_Hld2,Te_LSSSSS,Cd_MPP,Cd_MMM,Cd_CCC,Cd_MPM,Sl_FtoL,Te_SSSSLS,Cd_MMP,Cd_MMM,Cd_CCC,Cd_MPP,
		Sl_Hld2,Te_LSSSSS,Cd_MPM,Cd_MMP,Cd_CCC,Cd_MMM,Sl_FtoL,Te_SSSSLS,Cd_MPP,Cd_MMP,Cd_CCC,Cd_MPM,

		// [Program 1] code for NEIGHBOUR 1 enumerated face of same size
		FLG_XPX,
/*361*/	4, Sl_Ch01,		Te_SSLLSS, Cd_MMP,Cd_CCCn,Cd_CCC,Cd_MMM, Sl_Hld2,Te_SSSLSS,Cd_MMC,Cd_CCCn,Cd_CCC,Cd_MMM,Sl_FtoL,Te_SSSLSS,Cd_MMP,Cd_CCCn,Cd_CCC,Cd_MMC,
/*380*/	8, Sl_FtoL|Te_M,Te_SSLLSS, Cd_PMP,Cd_CCCn,Cd_CCC,Cd_MMP, Sl_FtoL,Te_SSSLSS,Cd_CMP,Cd_CCCn,Cd_CCC,Cd_MMP,Sl_FtoL,Te_SSSLSS,Cd_PMP,Cd_CCCn,Cd_CCC,Cd_CMP,
/*399*/	5, Sl_FtoL,		Te_SSLLSS, Cd_PMM,Cd_CCCn,Cd_CCC,Cd_PMP, Sl_FtoL,Te_SSSLSS,Cd_PMC,Cd_CCCn,Cd_CCC,Cd_PMP,Sl_FtoL,Te_SSSLSS,Cd_PMM,Cd_CCCn,Cd_CCC,Cd_PMC,
/*418*/	0, Sl_FtoL|Te_M,Te_SSLLSS, Cd_MMM,Cd_CCCn,Cd_CCC,Cd_PMM, Sl_FtoL,Te_SSSLSS,Cd_CMM,Cd_CCCn,Cd_CCC,Cd_PMM,Sl_FtoL,Te_SSSLSS,Cd_MMM,Cd_CCCn,Cd_CCC,Cd_CMM,
		// [Program 2] code for NEIGHBOUR 1 enumerated face with facing children
/*437*/	4, Sl_Ch01,		Te_SSSSSS, Cd_MMP,Cd_CMC,Cd_CCC,Cd_MMM, Sl_Hld2,Te_SSSSSL,Cd_MMC,Cd_CMC,Cd_CCC,Cd_MMM,Sl_FtoL,Te_SLSSSS,Cd_MMP,Cd_CMC,Cd_CCC,Cd_MMC,
/*456*/	8, Sl_FtoL|Te_M,Te_SSSSSS, Cd_PMP,Cd_CMC,Cd_CCC,Cd_MMP, Sl_FtoL,Te_SSSSSL,Cd_CMP,Cd_CMC,Cd_CCC,Cd_MMP,Sl_FtoL,Te_SLSSSS,Cd_PMP,Cd_CMC,Cd_CCC,Cd_CMP,
/*475*/	5, Sl_FtoL,		Te_SSSSSS, Cd_PMM,Cd_CMC,Cd_CCC,Cd_PMP, Sl_FtoL,Te_SSSSSL,Cd_PMC,Cd_CMC,Cd_CCC,Cd_PMP,Sl_FtoL,Te_SLSSSS,Cd_PMM,Cd_CMC,Cd_CCC,Cd_PMC,
/*494*/	0, Sl_FtoL|Te_M,Te_SSSSSS, Cd_MMM,Cd_CMC,Cd_CCC,Cd_PMM, Sl_FtoL,Te_SSSSSL,Cd_CMM,Cd_CMC,Cd_CCC,Cd_PMM,Sl_FtoL,Te_SLSSSS,Cd_MMM,Cd_CMC,Cd_CCC,Cd_CMM,
		// [Program 3] code for NEIGHBOUR 1 enumerated face of double size
/*513*/	Sl_Hld2,Te_LSSSSS,Cd_PMP,Cd_MMM,Cd_CCC,Cd_MMP,Sl_FtoL,Te_SSSSLS,Cd_PMM,Cd_MMM,Cd_CCC,Cd_PMP,
		Sl_Hld2,Te_LSSSSS,Cd_MMP,Cd_PMM,Cd_CCC,Cd_MMM,Sl_FtoL,Te_SSSSLS,Cd_PMP,Cd_PMM,Cd_CCC,Cd_MMP,

		// [Program 1] code for NEIGHBOUR 3 enumerated face of same size (this one cannot compact up with the previous on topological grounds)
		FLG_MXX,
/*538*/	7, Sl_Hld2|Te_M,Te_SSLLSS, Cd_PPM,Cd_CCCn,Cd_CCC,Cd_PPP, Sl_Hld2,Te_SSSLSS,Cd_PPC,Cd_CCCn,Cd_CCC,Cd_PPP,Sl_FtoL,Te_SSSLSS,Cd_PPM,Cd_CCCn,Cd_CCC,Cd_PPC,
/*576*/	2, Sl_FtoL,		Te_SSLLSS, Cd_PMM,Cd_CCCn,Cd_CCC,Cd_PPM, Sl_FtoL,Te_SSSLSS,Cd_PCM,Cd_CCCn,Cd_CCC,Cd_PPM,Sl_FtoL,Te_SSSLSS,Cd_PMM,Cd_CCCn,Cd_CCC,Cd_PCM,
/*557*/	5, Sl_FtoL|Te_M,Te_SSLLSS, Cd_PMP,Cd_CCCn,Cd_CCC,Cd_PMM, Sl_FtoL,Te_SSSLSS,Cd_PMC,Cd_CCCn,Cd_CCC,Cd_PMM,Sl_FtoL,Te_SSSLSS,Cd_PMP,Cd_CCCn,Cd_CCC,Cd_PMC,
/*595*/	10,Sl_FtoL,		Te_SSLLSS, Cd_PPP,Cd_CCCn,Cd_CCC,Cd_PMP, Sl_FtoL,Te_SSSLSS,Cd_PCP,Cd_CCCn,Cd_CCC,Cd_PMP,Sl_FtoL,Te_SSSLSS,Cd_PPP,Cd_CCCn,Cd_CCC,Cd_PCP,
		// [Program 2] code for NEIGHBOUR 3 enumerated face with facing children
/*614*/	7, Sl_Hld2|Te_M,Te_SSSSSS, Cd_PPM,Cd_PCC,Cd_CCC,Cd_PPP, Sl_Hld2,Te_SSSSSL,Cd_PPC,Cd_PCC,Cd_CCC,Cd_PPP,Sl_FtoL,Te_SLSSSS,Cd_PPM,Cd_PCC,Cd_CCC,Cd_PPC,
/*633*/	2, Sl_FtoL,		Te_SSSSSS, Cd_PMM,Cd_PCC,Cd_CCC,Cd_PPM, Sl_FtoL,Te_SSSSSL,Cd_PCM,Cd_PCC,Cd_CCC,Cd_PPM,Sl_FtoL,Te_SLSSSS,Cd_PMM,Cd_PCC,Cd_CCC,Cd_PCM,
/*652*/	5, Sl_FtoL|Te_M,Te_SSSSSS, Cd_PMP,Cd_PCC,Cd_CCC,Cd_PMM, Sl_FtoL,Te_SSSSSL,Cd_PMC,Cd_PCC,Cd_CCC,Cd_PMM,Sl_FtoL,Te_SLSSSS,Cd_PMP,Cd_PCC,Cd_CCC,Cd_PMC,
/*671*/	10,Sl_FtoL,		Te_SSSSSS, Cd_PPP,Cd_PCC,Cd_CCC,Cd_PMP, Sl_FtoL,Te_SSSSSL,Cd_PCP,Cd_PCC,Cd_CCC,Cd_PMP,Sl_FtoL,Te_SLSSSS,Cd_PPP,Cd_PCC,Cd_CCC,Cd_PCP,
		// [Program 3] code for NEIGHBOUR 3 enumerated face of double size
/*690*/	Sl_Hld2,Te_LSSSSS,Cd_PPP,Cd_PMM,Cd_CCC,Cd_PMP,Sl_FtoL,Te_SSSSLS,Cd_PPM,Cd_PMM,Cd_CCC,Cd_PPP,
		Sl_Hld2,Te_LSSSSS,Cd_PMP,Cd_PPM,Cd_CCC,Cd_PMM,Sl_FtoL,Te_SSSSLS,Cd_PPP,Cd_PPM,Cd_CCC,Cd_PMP,

		// [Program 1] code for NEIGHBOUR 5 enumerated face of same size
		FLG_XXM,
/*715*/	8, Sl_Ch01,		Te_SSLLSS, Cd_MMP,Cd_CCCn,Cd_CCC,Cd_PMP, Sl_Hld2,Te_SSSLSS,Cd_CMP,Cd_CCCn,Cd_CCC,Cd_PMP,Sl_FtoL,Te_SSSLSS,Cd_MMP,Cd_CCCn,Cd_CCC,Cd_CMP,
/*734*/	9, Sl_FtoL|Te_M,Te_SSLLSS, Cd_MPP,Cd_CCCn,Cd_CCC,Cd_MMP, Sl_FtoL,Te_SSSLSS,Cd_MCP,Cd_CCCn,Cd_CCC,Cd_MMP,Sl_FtoL,Te_SSSLSS,Cd_MPP,Cd_CCCn,Cd_CCC,Cd_MCP,
/*753*/	11,Sl_FtoL,		Te_SSLLSS, Cd_PPP,Cd_CCCn,Cd_CCC,Cd_MPP, Sl_FtoL,Te_SSSLSS,Cd_CPP,Cd_CCCn,Cd_CCC,Cd_MPP,Sl_FtoL,Te_SSSLSS,Cd_PPP,Cd_CCCn,Cd_CCC,Cd_CPP,
/*772*/	10,Sl_FtoL|Te_M,Te_SSLLSS, Cd_PMP,Cd_CCCn,Cd_CCC,Cd_PPP, Sl_FtoL,Te_SSSLSS,Cd_PCP,Cd_CCCn,Cd_CCC,Cd_PPP,Sl_FtoL,Te_SSSLSS,Cd_PMP,Cd_CCCn,Cd_CCC,Cd_PCP,
		// [Program 2] code for NEIGHBOUR 5 enumerated face with facing children
/*791*/	8, Sl_Ch01,		Te_SSSSSS, Cd_MMP,Cd_CCP,Cd_CCC,Cd_PMP, Sl_Hld2,Te_SSSSSL,Cd_CMP,Cd_CCP,Cd_CCC,Cd_PMP,Sl_FtoL,Te_SLSSSS,Cd_MMP,Cd_CCP,Cd_CCC,Cd_CMP,
/*810*/	9, Sl_FtoL|Te_M,Te_SSSSSS, Cd_MPP,Cd_CCP,Cd_CCC,Cd_MMP, Sl_FtoL,Te_SSSSSL,Cd_MCP,Cd_CCP,Cd_CCC,Cd_MMP,Sl_FtoL,Te_SLSSSS,Cd_MPP,Cd_CCP,Cd_CCC,Cd_MCP,
/*829*/	11,Sl_FtoL,		Te_SSSSSS, Cd_PPP,Cd_CCP,Cd_CCC,Cd_MPP, Sl_FtoL,Te_SSSSSL,Cd_CPP,Cd_CCP,Cd_CCC,Cd_MPP,Sl_FtoL,Te_SLSSSS,Cd_PPP,Cd_CCP,Cd_CCC,Cd_CPP,
/*848*/	10,Sl_FtoL|Te_M,Te_SSSSSS, Cd_PMP,Cd_CCP,Cd_CCC,Cd_PPP, Sl_FtoL,Te_SSSSSL,Cd_PCP,Cd_CCP,Cd_CCC,Cd_PPP,Sl_FtoL,Te_SLSSSS,Cd_PMP,Cd_CCP,Cd_CCC,Cd_PCP,
		// [Program 3] code for NEIGHBOUR 5 enumerated face of double size
/*867*/	Sl_Hld2,Te_LSSSSS,Cd_PPP,Cd_MMP,Cd_CCC,Cd_MPP,Sl_FtoL,Te_SSSSLS,Cd_PMP,Cd_MMP,Cd_CCC,Cd_PPP,
		Sl_Hld2,Te_LSSSSS,Cd_MPP,Cd_PMP,Cd_CCC,Cd_MMP,Sl_FtoL,Te_SSSSLS,Cd_PPP,Cd_PMP,Cd_CCC,Cd_MPP,

		// [Program 1] code for NEIGHBOUR 4 enumerated face of same size, note that the coordinate compaction (code Sl_FtoL) comes in seamlessly
		// from neighbour 2 if the two programs are executed in order
		FLG_XMX,
/*892*/	11,Sl_Ch01|Te_M,Te_SSLLSS, Cd_MPP,Cd_CCCn,Cd_CCC,Cd_PPP, Sl_Hld2,Te_SSSLSS,Cd_CPP,Cd_CCCn,Cd_CCC,Cd_PPP,Sl_FtoL,Te_SSSLSS,Cd_MPP,Cd_CCCn,Cd_CCC,Cd_CPP,
/*911*/	6, Sl_FtoL,		Te_SSLLSS, Cd_MPM,Cd_CCCn,Cd_CCC,Cd_MPP, Sl_FtoL,Te_SSSLSS,Cd_MPC,Cd_CCCn,Cd_CCC,Cd_MPP,Sl_FtoL,Te_SSSLSS,Cd_MPM,Cd_CCCn,Cd_CCC,Cd_MPC,
/*930*/	3, Sl_FtoL|Te_M,Te_SSLLSS, Cd_PPM,Cd_CCCn,Cd_CCC,Cd_MPM, Sl_FtoL,Te_SSSLSS,Cd_CPM,Cd_CCCn,Cd_CCC,Cd_MPM,Sl_FtoL,Te_SSSLSS,Cd_PPM,Cd_CCCn,Cd_CCC,Cd_CPM,
/*949*/	7, Sl_FtoL,		Te_SSLLSS, Cd_PPP,Cd_CCCn,Cd_CCC,Cd_PPM, Sl_FtoL,Te_SSSLSS,Cd_PPC,Cd_CCCn,Cd_CCC,Cd_PPM,Sl_FtoL,Te_SSSLSS,Cd_PPP,Cd_CCCn,Cd_CCC,Cd_PPC,
		// [Program 2] code for NEIGHBOUR 4 enumerated face with facing children
/*968*/	11,Sl_Ch01|Te_M,Te_SSSSSS, Cd_MPP,Cd_CPC,Cd_CCC,Cd_PPP, Sl_Hld2,Te_SSSSSL,Cd_CPP,Cd_CPC,Cd_CCC,Cd_PPP,Sl_FtoL,Te_SLSSSS,Cd_MPP,Cd_CPC,Cd_CCC,Cd_CPP,
/*987*/	6, Sl_FtoL,		Te_SSSSSS, Cd_MPM,Cd_CPC,Cd_CCC,Cd_MPP, Sl_FtoL,Te_SSSSSL,Cd_MPC,Cd_CPC,Cd_CCC,Cd_MPP,Sl_FtoL,Te_SLSSSS,Cd_MPM,Cd_CPC,Cd_CCC,Cd_MPC,
/*1006*/3, Sl_FtoL|Te_M,Te_SSSSSS, Cd_PPM,Cd_CPC,Cd_CCC,Cd_MPM, Sl_FtoL,Te_SSSSSL,Cd_CPM,Cd_CPC,Cd_CCC,Cd_MPM,Sl_FtoL,Te_SLSSSS,Cd_PPM,Cd_CPC,Cd_CCC,Cd_CPM,
/*1025*/7, Sl_FtoL,		Te_SSSSSS, Cd_PPP,Cd_CPC,Cd_CCC,Cd_PPM, Sl_FtoL,Te_SSSSSL,Cd_PPC,Cd_CPC,Cd_CCC,Cd_PPM,Sl_FtoL,Te_SLSSSS,Cd_PPP,Cd_CPC,Cd_CCC,Cd_PPC,
		// [Program 3] code for NEIGHBOUR 4 enumerated face of double size
/*1044*/Sl_Hld2,Te_LSSSSS,Cd_PPP,Cd_MPM,Cd_CCC,Cd_PPM,Sl_FtoL,Te_SSSSLS,Cd_MPP,Cd_MPM,Cd_CCC,Cd_PPP,
		Sl_Hld2,Te_LSSSSS,Cd_PPM,Cd_MPP,Cd_CCC,Cd_MPM,Sl_FtoL,Te_SSSSLS,Cd_PPP,Cd_MPP,Cd_CCC,Cd_PPM
	};
	
	// these masks clear and set pertinent node types in the stencil-matching pattern of an element (they only need to set a node to be a 0-node)
	final static int NINT = 0xFFFFFFFF;
	final static int CLR_n0=NINT-(3|3<<5|3<<10), CLR_n1=NINT-(3<<3|3<<15|3<<20), CLR_n2=NINT-(3<<8|3<<18|3<<25), CLR_n3=NINT-(3<<13|3<<23|3<<28);
	final static int FLG_n0_0=(2|2<<5|2<<10), FLG_n1_0=(2<<3|2<<15|2<<20), FLG_n2_0=(2<<8|2<<18|2<<25), FLG_n3_0=(2<<13|2<<23|2<<28);
	// these masks clear and set pertinent edge status bits in the status word, FLG_EBND = mask to test if all edges became internal/external
	final static int FLG_EST0=42, FLG_EST1=642, FLG_EST2=2184, FLG_EST3=2592, FLG_EBND = 2+8+32+128+512+2048;
	final static int CLR_EST0=NINT-FLG_EST0, CLR_EST1=NINT-FLG_EST1, CLR_EST2=NINT-FLG_EST2, CLR_EST3=NINT-FLG_EST3;
	final static int eStatClr[] = { CLR_EST0, CLR_EST1, CLR_EST2, CLR_EST3 };
	// these are the bits within pattern word that equal to the cut flags
	final static int FLG_CUT01=1<<2, FLG_CUT02=1<<7, FLG_CUT03=1<<12, FLG_CUT12=1<<17, FLG_CUT13=1<<22, FLG_CUT23=1<<27;
	final static int edgeCut[] = { FLG_CUT01, FLG_CUT02, FLG_CUT03, FLG_CUT12, FLG_CUT13, FLG_CUT23 };
	// the TETREFx bits lie in the status word and flag whether a node is a node index or a reference to another element's node index 
	final static int TETREF0 = 1<<22, TETREF1 = 1<<23, TETREF2 = 1<<24, TETREF3 = 1<<25, TETREF_CLR = NINT - (1<<22|1<<23|1<<24|1<<25);
	
	// method generates the internal & external tetrahedra for adjustment towards the boundary from a supplied latticeTree
	// tetrahedra are generated and separated into an internal paged array elementIntPg[] and a boundary page array elementBndPg[]
	// unique nodes and edges are generated from octant coordinates, the boundary element array is further used to do the edge intersection tests
	// to find cut and snap points: nodes close to boundary are snapped, the rest are stored as cut nodes for generating the boundary tetra stencils
	// method finally assembles the nodes into a compact tetrahedron stream, minimising redundant coordinate loads
	private int[] elementGeneratorIST(FEM1Octree latticeTree, FEM1Octree geocTree, FEM1Octant[] neighbourhoodIST, int elementTarget, int nodeTarget) {
		
		if (DEBUG_LEVEL > 1) System.out.println("**************** FEM1.elementGeneratorIST() stage 1 ****************");
		debugIST_phase1time = System.nanoTime();
		int lSubD1 = latticeTree.latticeSubdivs*2 + 1;				// note: since centroids also form nodes, the grid has to be  doubled for hashing
		int[] enReturn = {0,0};
		nodeHTableIST = new NodeHashTable(0, bBox[3] - bBox[0], lSubD1, bBox[0], bBox[1], bBox[2]);
		// allocate first pages of paging structures for internal & boundary elements
		elements = -1;
		elementPg = new int[elementPgN][];
		elementPg[0] = new int[elementPgSize];
		int count_ChgA=0, count_FtoL=0, count_Ch03=0, count_Hld2=0, count_Ch01=0;	// DEBUG: slotcode usage counters
				
		// this block investigates every octant's facing neighbours according to the byte-program, generating tetrahedra depending on neighbourhood
		// ----------------------------------------------------------------------------------------------------------------------------------------
		for (int n7 = 0; n7 < neighbourhoodIST.length; n7+=6) {
			FEM1Octant oct = neighbourhoodIST[n7++];
			if (oct == null) continue;
			int octStatus = oct.status>>3&0xFF;
			// don't need to process a lower level gradation if it's fully subdivided by it's children (the children will be doing the processing)
			if (octStatus == 0xFF) continue;
			oct.nodeI = new int[66];
			int nbrEflgs = edgeOctFlags[oct.nodes];

			int o = oct.status & 7;
			byte slotcode = Sl_ChgA;
			boolean forceChgA = true, forceHld2 = false, facedLarger = false, facedEqual = false, didProgramB = false;

			for (int seq=0, i1=1; seq < 6; seq++) {									// for each face-neighbour of this octant

				startPrintElement = elementsIST > 40 ? elementsIST - 40 : 0;		// tell toString() we want to print 40 elements in hindsight
				FEM1Octant nbrF = neighbourhoodIST[n7 + faceSeq_IST[seq]];
				int seqP = 6 + seq * 177;											// choose particular code group depending on neighbour-enumerated face
				byte facedFlg = faceSeq_IST[seqP++], facingFlg = (byte)(~facedFlg);

				if ((nbrF==null||nbrF.level<oct.level)&&(facingFlg&octStatus)==0) {	// if neighbour is boundary/lower level, DEBUG: supposedly 2:1 balanced tree, unnecessary test?
					i1 = 2;															// progress to 2:1 tetra generation program 3 for that face
					switch (faceSeq_IST[seq]) {										// choose DIAGONAL SPLIT DIRECTION according to position of this octant
					case 0: seqP += (o==OCT_MMM || o==OCT_PPM || o==OCT_MMP || o==OCT_PPP) ? 152 : 164; break;
					case 5: seqP += (o==OCT_MMP || o==OCT_PPP || o==OCT_MMM || o==OCT_PPM) ? 152 : 164; break;
					case 4: seqP += (o==OCT_PPM || o==OCT_MPP || o==OCT_PMM || o==OCT_MMP) ? 164 : 152; break;
					case 1: seqP += (o==OCT_PMM || o==OCT_MMP || o==OCT_PPM || o==OCT_MPP) ? 164 : 152; break;
					case 2: seqP += (o==OCT_MMP || o==OCT_MPM || o==OCT_PMP || o==OCT_PPM) ? 164 : 152; break;
					case 3: seqP += (o==OCT_PMP || o==OCT_PPM || o==OCT_MMP || o==OCT_MPM) ? 164 : 152;
					}
					while (i1-- > 0) {												// the program writes two tetrahedra

						facedLarger = true;											// signal to equal-octants generator that we've been here
						if (forceChgA) { slotcode=Sl_ChgA; seqP++; }				// always change all nodes on new octant
						// if change of 3 nodes was flagged or previous program was 1(B)/2(B), hold centroid node only
						else if (forceHld2 || facedEqual) { slotcode = Sl_Hld2; seqP++; facedEqual = false; }
						else slotcode = faceSeq_IST[seqP++];
						int edgeCode = faceSeq_IST[seqP++] << 16;
						// check internal status of the expected nodes, skip fully external tetrahedra (Sl_ChgA breaks compaction sequence)
						if (!oct.centroid_internal() && externalTetrahedron(oct, nbrF, seqP)) { seqP += 4; forceHld2 = true; continue; }
						int seq24 = seq<<24;

						switch (slotcode&0x7F) {
						case Sl_ChgA:	processISTnode(oct, nbrF, slotcode&0xFF|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
										processISTnode(oct, nbrF, Sl_Chg1|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
										processISTnode(oct, nbrF, Sl_Chg2|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
										processISTnode(oct, nbrF, Sl_Chg3|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
										forceChgA=false; count_ChgA++; break;
						case Sl_FtoL:	processISTnode(oct, nbrF, slotcode&0xFF|faceSeq_IST[seqP]<<8|edgeCode|seq24);
										seqP += 4; count_FtoL++; break;
						case Sl_Hld2:	processISTnode(oct, nbrF, slotcode&0xFF|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
										processISTnode(oct, nbrF, Sl_Chg1|faceSeq_IST[seqP++]<<8|edgeCode|seq24); seqP++;
										processISTnode(oct, nbrF, Sl_Chg3|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
										forceHld2 = false; count_Hld2++; break;
						}
					}
					
				} else if (nbrF!= null && nbrF.level == oct.level) {				// if it's an existing neighbour of same size
					
					facedEqual = true;												// signal to hi-vs-low octant generator that we've been here
					boolean program2 = false;
					// if this octant's aspect-facing children are forming a complete wall versus that neighbour, do premature abortion of edge tests
					if ((facingFlg & octStatus) == facingFlg) { forceHld2 = true; continue; }

					// if neighbour's children face octant OR octant's children face neighbour (=there's a node at face center)
					if ((facedFlg & nbrF.status>>3&0xFF) != 0 || (facingFlg & octStatus) != 0)
						{ seqP += 76; program2 = true; }							// skip to program 2, tetra-generation versus neighbour child/children

					// if there's no node at face center, then INTRUSION into neighbour must be done -> only allowed towards x/y/z-POSITIVE neighbours
					else if (!(oct.xC<nbrF.xC || oct.yC<nbrF.yC || oct.zC<nbrF.zC))	// process only neighbours above oneself in ordinate-increasing order
						{ forceHld2 = true; continue; }								// which eliminates duplicate work (only for equal-level construction!)
//					else if ((nbrF.facets&(1<<(5-faceSeq_IST[seq]))) != 0)			// if neighbour had processed this face aspect, avoid double job
//						{ forceHld2 = true; continue; }								// note: nodepacking sequence is interrupted -> Sl_Hld2

					int e = 4;
					boolean axisLoaded = false, programB = false;
					
					while (e-- > 0) {												// test the 4 edges of this face
						int n = faceSeq_IST[seqP++];								// get the neighbour of current edge
						i1 = 1;

						// test whether an equal-size corner neighbour has children toward this octant (thus a 2:1 relation)
						if ((nbrEflgs&1<<n) != 0) {									// and if it does...
							seqP += 6; i1++; programB = true; }						// ...skip to program 1B or 2B (which outputs two tetrahedra)
						else programB = false;

						while (i1-- > 0) {											// write one or two tetrahedra, depending on program
							// do not write a quadrisected BCC (program 2B) for current corner if octant has a child there
							if (program2 && programB && oct.octant != null) {
								if      (i1==1 && oct.octant[faceSeq_IST[seqP+5]] != null) { seqP += 6; forceHld2 = true; continue; }
								else if (i1==0 && oct.octant[faceSeq_IST[seqP+2]] != null) { seqP += 6; forceHld2 = true; break; }
							}
							// on starting a new octant, forceChgA will always be true
							if (forceChgA) { slotcode = (byte)(Sl_ChgA|(faceSeq_IST[seqP++]&0x80)); 
							// if some element of octant was done by SOME program, the worst case node replacement is to keep only centroid
							} else if (forceHld2 || facedLarger) {
								slotcode = (byte)((axisLoaded ? Sl_Ch03 : Sl_Hld2)|(faceSeq_IST[seqP++]&0x80));
							// on changing back to main program from program B, if axis was loaded, need to replace only 2 nodes, 0/3
							} else if (!programB && didProgramB) {
								slotcode = (byte)((axisLoaded ? Sl_Ch03 : Sl_Hld2)|(faceSeq_IST[seqP++]&0x80));
							// if the program sequence was unbroken, follow the sequence's compaction codes
							} else slotcode = faceSeq_IST[seqP++];
							int edgeCode = (int)faceSeq_IST[seqP++] << 16;
							
							// check internal status of the expected nodes, skip fully external tetrahedra (keep centroid node)
							if (!oct.centroid_internal() && externalTetrahedron(oct, nbrF, seqP))
								{ seqP += programB ? 4 : 16; forceHld2 = true; continue; }
							int seq24 = seq>>24;
							
							switch (slotcode & 0x7F) {
							case Sl_ChgA:	processISTnode(oct, nbrF, slotcode&0xFF|(faceSeq_IST[seqP++]<<8)|edgeCode|seq24);
											processISTnode(oct, nbrF, Sl_Chg1|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
											processISTnode(oct, nbrF, Sl_Chg2|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
											processISTnode(oct, nbrF, Sl_Chg3|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
											forceChgA = false; if (!programB) seqP += 12;
											count_ChgA++; break;
							case Sl_FtoL:	processISTnode(oct, nbrF, slotcode&0xFF|faceSeq_IST[seqP]<<8|edgeCode|seq24);
											seqP += programB ? 4 : 16;
											count_FtoL++; break;
							case Sl_Hld2:	processISTnode(oct, nbrF, slotcode&0xFF|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
											processISTnode(oct, nbrF, Sl_Chg1|faceSeq_IST[seqP++]<<8|edgeCode|seq24); seqP++;
											processISTnode(oct, nbrF, Sl_Chg3|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
											forceHld2 = false; if (!programB) seqP += 12;
											count_Hld2++; break;
							case Sl_Ch01:	processISTnode(oct, nbrF, slotcode&0xFF|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
											processISTnode(oct, nbrF, Sl_Chg1|faceSeq_IST[seqP--]<<8|edgeCode|seq24);
											seqP += programB ? 4 : 16;
											count_Ch01++; break;
							case Sl_Ch03:	processISTnode(oct, nbrF, slotcode&0xFF|faceSeq_IST[seqP]<<8|edgeCode|seq24); seqP+=3;
											processISTnode(oct, nbrF, Sl_Chg3|faceSeq_IST[seqP++]<<8|edgeCode|seq24);
											if (!programB) seqP += 12;
											count_Ch03++;
							}
							axisLoaded = true;							// flag that centroidal nodes (1&2) were stored
							facedLarger = false;
							didProgramB = programB;
						}
					}
					//oct.facets |= 1 << faceSeq_IST[seq];				// flag this face aspect as processed
				}
			}
			oct.nodeI = null;											// deallocate local node retriever
		}
		processISTnode(null, null, SET31);								// write out the last element

		// gather paged boundary elements into elementIST[] array
		elementIST = new int[elementsIST * 6];
		for (int elemNo = 0, p = 0, elems = elementsIST * 6; p < elementPgN; p++, elems -= elementPgSize) {
			int[] elemPage = elementPg[p];
			if (elemPage == null) break;
			for (int e6 = 0, eEnd = elems < elementPgSize ? elems: elementPgSize; e6 < eEnd; e6++, elemNo++)
				elementIST[elemNo] = elemPage[e6];
		}
		readElementPg = false;											// DEBUG: stop toString() of elementPg structure
		elementPg = null;												// kill the boundary elements paging array

		nodeIST = nodeHTableIST.arrayInterleaved();						// collect hashed nodes into single coordinate array
		nodesIST = nodeHTableIST.count();
		debugIST_phase1time = System.nanoTime() - debugIST_phase1time;	// time phase 1
		
		if (DEBUG_LEVEL > 1) {
			System.out.print("Primary tetrahedra generated: "+ elementsIST);
			System.out.print("\nUnique nodes accumulated: " + nodeHTableIST.count());
			System.out.print("\nTime spent in FEM1.processISTnode(): " + debugIST_pISTn_time + " ns\n");
			System.out.print("\nLocally (prematurely) unique node cases: " + debugIST_locallyUnique);
			System.out.print("\nSlotCode counts:\nChgA: " + count_ChgA + "\nFtoL: " + count_FtoL);
			System.out.print("\nCh03: " + count_Ch03 + "\nHld2: " + count_Hld2 + "\nCh01: " + count_Ch01+"\n\n");
			System.out.println("**************** FEM1.elementGeneratorIST() stage 2 ****************");
		}
		
		// time to do the close-cuts node snapping procedure (doing intersection tests and moving nodes violated by cuts to the cuts)
		// --------------------------------------------------------------------------------------------------------------------------
		debugIST_phase2time = System.nanoTime();
		snapNodeStatus = new byte[nodesIST];							// storage of node snapping statuses
		printSnapNodes = true;
		int[] testElemsIST = new int[(elementsIST-elementsIntIST)*5];	// fast-access of boundary element nodes for post-stage 2 degenerate element test
		int elemsBnd = 0;
		
		lStepI_pISTe = 2/latticeStep; lStepI2_pISTe = 1/latticeStep;	// gridsteps for processISTedge() lattice precalculated data access
		divCofs_pISTe = latticeSubdivs + 1;								// offset to centroid lattice ray data for processISTedge()
		intOfsX_pISTe = -geocTree.root.xM + maxZordinDev; intOfsY_pISTe = -geocTree.root.yM + maxZordinDev;
		for (int s = 0; s < nodesIST; s++) snapNodeStatus[s] = 127;		// set them all to "max distance" (an unattainable value)
		int[] eStatus = {0, -1};										// indirect values (snap type, fast facet) returned from processISTedge()

		int cyclicBufSizeM1 = 15, seekAvg = 0, seekAvgIters = 0; 		// size of cyclic buffer - 1, average elements searched over N runs, runs counter
		long[] nonsharedCutNode = new long[16];							// fast storage area for lookup of non-shared edges
		int ceIdx = cyclicBufSizeM1;									// index pointer into nonsharedCutNode[] & counter of stored cutnodes				
		int count_internal=0, count_boundary=0, count_cyclic_finds=0, count_cyclic_misses=0;	// DEBUG counters
		int[] count_edgetype = {0,0,0,0};								// DEBUG edge type counters
		
		for (int e6 = 0, e = 0; e < elementsIST; e++) {
			startPrintElement = e > 6 ? e - 6 : 0;											// tell toString() we want to print from 9 elements back
			int nX0 = elementIST[e6++], status = elementIST[e6+4];
			if ((status & FLG_EBND) == 0) { e6+=5; count_internal++; continue; }			// skip internal elements

			int nX1 = elementIST[e6++], nX2 = elementIST[e6++], nX3 = elementIST[e6++], edgeGenFlags = status & 0xFFF;
			int edge = 0;
			// nonunique nodes can/will be referred to be fetched from their origin element, the TETREFx flags tell which need to be retrieved
			int n0, n1, n2, n3, eX0=-1, eX1=-1, eX2=-1, eX3=-1;
			if ((status&TETREF0)!=0) { n0=elementIST[(eX0 = nX0&0x3FFFFFFF)*6 + (nX0=nX0>>30&3)]; } else n0=nX0;
			if ((status&TETREF1)!=0) { n1=elementIST[(eX1 = nX1&0x3FFFFFFF)*6 + (nX1=nX1>>30&3)]; } else n1=nX1;
			if ((status&TETREF2)!=0) { n2=elementIST[(eX2 = nX2&0x3FFFFFFF)*6 + (nX2=nX2>>30&3)]; } else n2=nX2;
			if ((status&TETREF3)!=0) { n3=elementIST[(eX3 = nX3&0x3FFFFFFF)*6 + (nX3=nX3>>30&3)]; } else n3=nX3;
			int n0F = n0&0x3FFFFFFF, n1F = n1&0x3FFFFFFF, n2F = n2&0x3FFFFFFF, n3F = n3&0x3FFFFFFF;

			eStatus[1] = -1;																// specify fast facet nonexistent at start of edge processing
			while (edge < 6) {
				int eGF = edgeGenFlags & 3;
				count_edgetype[eGF]++;
				boolean nonshared = false, doCutTest = false;

				if (eGF == 3) doCutTest = true;
				else if (eGF == 2) {
					switch (edge) {
					// check for previously unknown unique edge (one that uses 2 nodes from 2 elements, but forms first occurrence in space of that edge)
					case 0:	nonshared=eX0!=eX1&&(eX0<eX1?(eX0<0?-2:transitiveRef(eX0,eX1,nX0,nX1)):(eX1<0?-2:transitiveRef(eX1,eX0,nX1,nX0)))==-1; break;
					case 1:	nonshared=eX0!=eX2&&(eX0<eX2?(eX0<0?-2:transitiveRef(eX0,eX2,nX0,nX2)):(eX2<0?-2:transitiveRef(eX2,eX0,nX2,nX0)))==-1; break;
					case 2:	nonshared=eX0!=eX3&&(eX0<eX3?(eX0<0?-2:transitiveRef(eX0,eX3,nX0,nX3)):(eX3<0?-2:transitiveRef(eX3,eX0,nX3,nX0)))==-1; break;
					case 3:	nonshared=eX1!=eX2&&(eX1<eX2?(eX1<0?-2:transitiveRef(eX1,eX2,nX1,nX2)):(eX2<0?-2:transitiveRef(eX2,eX1,nX2,nX1)))==-1; break;
					case 4: nonshared=eX1!=eX3&&(eX1<eX3?(eX1<0?-2:transitiveRef(eX1,eX3,nX1,nX3)):(eX3<0?-2:transitiveRef(eX3,eX1,nX3,nX1)))==-1; break;
					case 5:	nonshared=eX2!=eX3&&(eX2<eX3?(eX2<0?-2:transitiveRef(eX2,eX3,nX2,nX3)):(eX3<0?-2:transitiveRef(eX3,eX2,nX3,nX2)))==-1;
					}
					if (nonshared) debugIST_nonsharedEdge2++;
				} else { edge++; edgeGenFlags >>= 2; continue; }							// skip nonboundary edge

				if (nonshared) {														// if edge unique but formed by two old nodes
					long edgeKey = 0; boolean foundEdge = false;
					switch (edge) {
					case 0: edgeKey = n0F < n1F ? (long)n1F<<32|n0F : (long)n0F<<32|n1F; break;
					case 1: edgeKey = n0F < n2F ? (long)n2F<<32|n0F : (long)n0F<<32|n2F; break;
					case 2: edgeKey = n0F < n3F ? (long)n3F<<32|n0F : (long)n0F<<32|n3F; break;
					case 3: edgeKey = n1F < n2F ? (long)n2F<<32|n1F : (long)n1F<<32|n2F; break;
					case 4: edgeKey = n1F < n3F ? (long)n3F<<32|n1F : (long)n1F<<32|n3F; break;
					case 5: edgeKey = n2F < n3F ? (long)n3F<<32|n2F : (long)n2F<<32|n3F;
					}
					// see if the supposed nonshared edge is stored in fast cyclic array, if so, the edge was already processed
					//for (int iStart = (ceIdx + 1) & cyclicBufSizeM1, i = iStart; i != ceIdx; i = ++i & cyclicBufSizeM1)
					for (int i = (ceIdx + 1) & cyclicBufSizeM1; i != ceIdx; i = ++i & cyclicBufSizeM1)
						if (nonsharedCutNode[i] == edgeKey) {
							foundEdge = true;
							//seekAvg += (i - iStart)&cyclicBufSizeM1;				// DEBUG: finds out a fitting size (=power of 2!) for cyclic buffer
							//if (++seekAvgIters >= 8) { System.out.println(seekAvg/8); seekAvgIters = seekAvg = 0; }
							count_cyclic_finds++; break; }

					if (!foundEdge) {
						nonsharedCutNode[ceIdx] = edgeKey; ceIdx = --ceIdx&cyclicBufSizeM1;		// edge not found, store it's key
						doCutTest = true; count_cyclic_misses++; }								// edge will be processed
				}

				if (doCutTest) {		// cut test also receives the long/short edge status bits deciding the proximity alpha
					eStatus[0] = 0;
					switch (edge) {
					case 0:	processISTedge(geocTree, n0, n1, eStatus, (status&1<<16)!=0); break;		// edge 01
					case 1:	processISTedge(geocTree, n0, n2, eStatus, (status&1<<17)!=0); break;		// edge 02
					case 2:	processISTedge(geocTree, n0, n3, eStatus, (status&1<<18)!=0); break;		// edge 03
					case 3:	processISTedge(geocTree, n1, n2, eStatus, (status&1<<19)!=0); break;		// edge 12
					case 4:	processISTedge(geocTree, n1, n3, eStatus, (status&1<<20)!=0); break;		// edge 13
					case 5:	processISTedge(geocTree, n2, n3, eStatus, (status&1<<21)!=0); }				// edge 23
				}
				edge++; edgeGenFlags >>= 2;
			}
			// store down boundary elements in direct nodeindexed array for degeneracy check
			testElemsIST[elemsBnd++]=n0; testElemsIST[elemsBnd++]=n1; testElemsIST[elemsBnd++]=n2; testElemsIST[elemsBnd++]=n3; testElemsIST[elemsBnd++]=e;
			e6+=2;
		}
			
		// DEBUG: finds an erroneous element whose centroid is closest to the supplied coordinate
		//startPrintElement = testElemsIST[findElementIndex(testElemsIST, nodeIST, 2.402, 0.275, 1.813, 1, 6)*5+4];
		
		// do post-stage 2 degenerate element test
		debugIST_analysisTime = System.nanoTime();
		if (!discardAllSliversIST) elementNhoodIST = new int[nodesIST][];
		int sliverCount = 0;
		int count_unsnaps=0, count_degeneracy_checks = 0, count_degenerates = 0, count_defectives = 0;
		
		for (int e5 = 0, e5S = 0; e5 < elemsBnd;) {
			double Q = 0;
			int n0 = testElemsIST[e5++], n1 = testElemsIST[e5++], n2 = testElemsIST[e5++], n3 = testElemsIST[e5++], eIdx = testElemsIST[e5++];
			startPrintElement = eIdx;		
			int n0F = n0 & 0x3FFFFFFF, n1F = n1 & 0x3FFFFFFF, n2F = n2 & 0x3FFFFFFF, n3F = n3 & 0x3FFFFFFF;
			byte snap0=snapNodeStatus[n0F], snap1=snapNodeStatus[n1F], snap2=snapNodeStatus[n2F], snap3=snapNodeStatus[n3F];
			// flag boundary slivers, then flag them as degenerate or defective for discarding
			if (snap0<0 && snap1<0 && snap2<0 && snap3<0) {
				
				if (discardAllSliversIST) {
					// all slivers were requested for discarding
					elementIST[eIdx*6+5] |= SET28|SET29;
					
				} else if (tetraVolumePositivityIST(n0F, n1F, n2F, n3F, snap0, snap1, snap2, snap3) < 0.000001) {
					// find and discard degenerate/negative volume slivers (a normal case for externalised boundary elements at lor resolution)
					if (DEBUG_LEVEL > 2) System.out.println("Boundary sliver element " + eIdx + " negative/zero-volume, deleted.");
					elementIST[eIdx*6+5] |= SET28|SET29; count_degenerates++;
				} else {
					// quality sift according to 3*inner/circumradius criterion, stacked slivers cases will now note missing neighbours
					if ((Q = tetraQualityIST(n0F, n1F, n2F, n3F, snap0, snap1, snap2, snap3)) < minSliverQuality) {
						
						if (DEBUG_LEVEL > 1) System.out.println("Boundary sliver element "+eIdx+" of too low quality: "+
																	String.format("%.3f", Q)+(discardBadSliversIST ? ", deleted." : "."));
						if (discardBadSliversIST) elementIST[eIdx*6+5] |= SET28;
						count_defectives++;
						
					} else if (!sliverCentroidInternal(geocTree, n0F, n1F, n2F, n3F)) {
						if (DEBUG_LEVEL > 2) System.out.println("Boundary sliver element "+eIdx+" a straddler, deleted.");
						elementIST[eIdx*6+5] |= SET28; count_defectives++;
					} else {
						// store found nondefective slivers in same array, we'll have to test them for stacking with defective ones
						testElemsIST[e5S++]=n0F; testElemsIST[e5S++]=n1F; testElemsIST[e5S++]=n2F; testElemsIST[e5S++]=n3F; testElemsIST[e5S++]=eIdx;
						sliverCount++; }
					// note down sliver-sliver neighbourhoods to eliminate the stacked sliver cases 
					nodeSupportAddElementIST(n0F, eIdx); nodeSupportAddElementIST(n1F, eIdx);
					nodeSupportAddElementIST(n2F, eIdx); nodeSupportAddElementIST(n3F, eIdx);
					elementIST[eIdx*6+5] |= SET29;
					count_boundary++;
				}
				continue;
			}
			int sC0=snap0<0?1:0, sC1=snap1<0?1:0, sC2=snap2<0?1:0, sC3=snap3<0?1:0, snapSum = sC0 + sC1 + sC2 + sC3;
			
			// we're interested in elements with at least two snapped nodes and at least one internal node
			if (snapSum > 1 && ((n0&SET30)!=0&&sC0==0||(n1&SET30)!=0&&sC1==0||(n2&SET30)!=0&&sC2==0||(n3&SET30)!=0&&sC3==0)) {
				if ((Q = tetraQualityIST(n0F, n1F, n2F, n3F, snap0, snap1, snap2, snap3)) < minUnsnapQuality) {
					if (DEBUG_LEVEL > 1) {
						System.out.println("Boundary element "+eIdx+" defective: "+String.format("%.3f", Q)+". Element is unsnapped");
						if (DEBUG_LEVEL > 2) {
							int n03 = n0F * 6, n13 = n1F * 6, n23 = n2F * 6, n33 = n3F * 6;
							if (snap0 < 0) n03+=3; System.out.print(String.format("      n0: %.3f,%.3f,%.3f, ",nodeIST[n03++],nodeIST[n03++],nodeIST[n03]));
							if (snap1 < 0) n13+=3; System.out.print(String.format("      n1: %.3f,%.3f,%.3f, ",nodeIST[n13++],nodeIST[n13++],nodeIST[n13]));
							if (snap2 < 0) n23+=3; System.out.print(String.format("      n2: %.3f,%.3f,%.3f, ",nodeIST[n23++],nodeIST[n23++],nodeIST[n23]));
							if (snap3 < 0) n33+=3; System.out.print(String.format("      n3: %.3f,%.3f,%.3f  ",nodeIST[n33++],nodeIST[n33++],nodeIST[n33]));
							System.out.println();
						}
						count_defectives++;
					}
					
					boolean success = false;
					// prioritise snapnodes according to snap proximity, unsnap the node that has to move the LEAST distance for quality approval
					int[] prio = {3, 2, 1, 0};									// DEBUG: reverse count to start unsnapping the MOST distant snaps instead
					int val0 = snap0&127, val1 = snap1&127, val2 = snap2&127, val3 = snap3&127;
					if (val0 > val3) { int tmp = prio[0]; prio[0]=prio[3]; prio[3]=tmp; tmp = val0; val0=val3; val3=tmp; }	
					if (val1 > val2) { int tmp = prio[1]; prio[1]=prio[2]; prio[2]=tmp; tmp = val1; val1=val2; val2=tmp; }
					if (val0 > val1) { int tmp = prio[0]; prio[0]=prio[1]; prio[1]=tmp; }
					if (val2 > val3) { int tmp = prio[2]; prio[2]=prio[3]; prio[3]=tmp; }
				
					for (int p = 0; !success && p < 4; p++)						// unsnap in order until the criterion is fulfilled
						switch (prio[p]) {
						case 0: if (snap0 < 0) { snap0 = snapNodeStatus[n0F]&=127; count_unsnaps++;
								if (tetraQualityIST(n0F,n1F,n2F,n3F, snap0, snap1, snap2, snap3) >= minUnsnapQuality) success=true;} break;
						case 1: if (snap1 < 0) { snap1 = snapNodeStatus[n1F]&=127; count_unsnaps++;
								if (tetraQualityIST(n0F,n1F,n2F,n3F, snap0, snap1, snap2, snap3) >= minUnsnapQuality) success=true;} break;
						case 2: if (snap2 < 0) { snap2 = snapNodeStatus[n2F]&=127; count_unsnaps++;
								if (tetraQualityIST(n0F,n1F,n2F,n3F, snap0, snap1, snap2, snap3) >= minUnsnapQuality) success=true;} break;
						case 3: if (snap3 < 0) { snap3 = snapNodeStatus[n3F]&=127; count_unsnaps++;
								if (tetraQualityIST(n0F,n1F,n2F,n3F, snap0, snap1, snap2, snap3) >= minUnsnapQuality) success=true;} }
					count_boundary++;
				// DEBUG: straddling discarding needs much better (and intensive) heuristics, avoid using
				} else if (discardStraddlers && snapSum > 2 && !elementCentroidInternal3S(geocTree, n0F, n1F, n2F, n3F, snap0, snap1, snap2, snap3)) {
					if (DEBUG_LEVEL > 1) System.out.println("Boundary element "+eIdx+" a 3-snap straddler, deleted.");
					elementIST[eIdx*6+5] |= SET28; count_defectives++;				
				} else if (discardStraddlers && snapSum > 1 && !elementCentroidInternal2S(geocTree, n0F, n1F, n2F, n3F, snap0, snap1, snap2, snap3)) {
					if (DEBUG_LEVEL > 1) System.out.println("Boundary element "+eIdx+" a 2-snap straddler, deleted.");
					elementIST[eIdx*6+5] |= SET28; count_defectives++; }
				else count_boundary++;
				count_degeneracy_checks++;
			}
		}
				
		// this iteration discards those stacked slivers whose neighbours have been discarded in previous iteration
		// note: the sliver neighbourhoods had been accumulated during the sliver testing iteration
		if (!discardAllSliversIST) {
			nbrComp12 = new int[elementSupportMaxArrayL * 4 + 1];							// prepare fixed arrays (nonthreaded) for neighbourhood analysis
			nbrComp34 = new int[elementSupportMaxArrayL * 4 + 1];
			nbrComp1234F = new int[elementSupportMaxArrayL * 8];
			nbrComp1234e = new int[elementSupportMaxArrayL * 8]; }

		for (int e5 = 0, e5end = sliverCount * 5; e5 < e5end;) {
			int n0=testElemsIST[e5++], n1=testElemsIST[e5++], n2=testElemsIST[e5++], n3=testElemsIST[e5++], eIdx=testElemsIST[e5++], eSoffs=eIdx*6+5;
			startPrintElement = eIdx;		
			if (snapNodeStatus[n0]>=0||snapNodeStatus[n1]>=0||snapNodeStatus[n2]>=0||snapNodeStatus[n3]>=0) {
				elementIST[eSoffs] &= CLR29; continue; }							// set unsnapped sliver as boundary element
			
			if (discardBadSliversIST) {
				if ((elementIST[eSoffs]&SET28) != 0) continue;						// check if some sliver already eliminated the group
				int[] nList = {eIdx, n0, n1, n2, n3};
				int[][] ngbrs = tetraNeighbours(nList, false, true);				// test if sliver is a "hanger", with <=1 facet neighbours
				int[] ngbrsF = ngbrs[0];	
				for (int i = 1, iEnd = ngbrsF[0]; i <= iEnd; i++)
					if ((elementIST[ngbrsF[i]*6 + 5] & SET28) != 0) {				// was this sliver in this neighbourhood discarded? 
						if (DEBUG_LEVEL > 1) {
							System.out.print("Sliver stack discarded: "+eIdx+", ");
							for (int i1=1; i1 <= iEnd; i1++) System.out.print(ngbrsF[i1]+(i1==iEnd?".\n":", "));
						}
						elementIST[eSoffs] |= SET28; count_boundary--;
						while (--i >= 1) {											// make sure to discard entire group
							elementIST[ngbrsF[i]*6 + 5] |= SET28;
						}
						break;
					}
			}
		}
		nbrComp12 = nbrComp34 = nbrComp1234F = nbrComp1234e = null;
		testElemsIST = null; elementNhoodIST = null;

		// if element target count exceeded, return it's negative value
		if (elementTarget>0 && count_internal+count_boundary > elementTarget) {
			if (DEBUG_LEVEL > 1) System.out.println("FEM1.elementGeneratorIST(): Element target exceeded, aborting.");
			enReturn[0] = -(count_internal+count_boundary); return enReturn; }

		debugIST_analysisTime = System.nanoTime() - debugIST_analysisTime;				// time sliver analysis time
		debugIST_phase2time = System.nanoTime() - debugIST_phase2time;					// time phase 2

		if (DEBUG_LEVEL > 1) {
			System.out.print("\nInternal tetrahedra skipped: " + count_internal);
			System.out.print("\nEdge type counts:\n    int/ext uniques: " + count_edgetype[3] + "\n    boundary uniques: " + count_edgetype[1]);
			System.out.print("\n    int/ext nonuniques: " + count_edgetype[2] + "\n    boundary nonuniques: " + count_edgetype[0]);
			System.out.print("\nTotal edges: " + (count_edgetype[0]+count_edgetype[1]+ count_edgetype[2]+ count_edgetype[3]));
			System.out.print(", of those, nonshared edges: " + debugIST_nonsharedEdge2);
			System.out.print("\nVisits to FEM1.processISTedge(): " + debugIST_pISTe_visits + ", fast Z-ordinate visits: " + debugIST_pISTe_fastVisits);
			System.out.print(", fast facets visits: " + debugIST_pISTe_fastFacetVisits + ", skip visits: " + debugIST_pISTe_skipVisits);
			System.out.print("\nTime spent in FEM1.processISTedge(): " + debugIST_pISTe_time2 + " ns\n");
			System.out.print("Wasted cut functions in FEM1.processISTedge(): " + debugIST_pISTe_vainSnaps);
			System.out.print("\nNonshared edges cyclic buffer cases: " + (count_cyclic_misses + count_cyclic_finds));
			System.out.print(" of which directly resolved: " + count_cyclic_finds);
			System.out.print("\nSlivers found: " + sliverCount);
			System.out.print("\nUnsnapped elements: " + count_unsnaps);
			System.out.print("\nDegenerated checks: "+count_degeneracy_checks+", degenerated: "+count_degenerates+", defectives: "+count_defectives+"\n\n");
			System.out.println("**************** FEM1.elementGeneratorIST() stage 3 ****************");
		}
		
		// after snapping nodes, it remains to find the cutnodes and register them in node hash, patterns and cutNodes[] array
		// -------------------------------------------------------------------------------------------------------------------
		debugIST_phase3time = System.nanoTime();
		nodeHTableIST = new NodeHashTable(nodesIST, bBox[3] - bBox[0], lSubD1, bBox[0], bBox[1], bBox[2]);	// clear node hash table for the cut nodes
		cutNodes = new int[elementsIST * 6];							// storage of the cut nodes in edge-indexed clusters
		nonsharedCutNode = new long[(cyclicBufSizeM1 = 15)+1];
		ceIdx = cyclicBufSizeM1 - 1;
		count_cyclic_finds=0; count_cyclic_misses=0;
		int count_becameIntExt=0, count_cutnodes=0;
		count_edgetype = new int[4];
		debugIST_pISTe_visits = debugIST_pISTe_fastVisits = debugIST_pISTe_fastFacetVisits = 0;
		debugIST_pISTe_time1 = debugIST_pISTe_time2; debugIST_pISTe_time2 = 0;		
		
		for (int e6 = 0, e = 0; e < elementsIST; e++) {
			
			startPrintElement = e > 7 ? e - 7 : 0;			
			boolean nonSharedEdge = false;
			int n0 = elementIST[e6++], status = elementIST[e6+4];
			if ((status&SET28) != 0) { e6+=5; continue; }										// skip defective elements
			int n1 = elementIST[e6++], n2 = elementIST[e6++], n3 = elementIST[e6++], pattern = elementIST[e6++]; e6++;
			
			// nonunique nodes can/will be dereferenced from their origin element, the TETREFx flags tell which need to be retrieved
			// the previously internal nodes that became snapped (bit 31 set) are deflagged from being internal (bit 30 cleared)
			int eX0=-1, eX1=-1, eX2=-1, eX3=-1, nX0=-1, nX1=-1, nX2=-1, nX3=-1;
			if ((status&TETREF0)!=0) n0=elementIST[(eX0=n0&0x3FFFFFFF)*6 + (nX0=n0>>30&3)]; 
			if ((status&TETREF1)!=0) n1=elementIST[(eX1=n1&0x3FFFFFFF)*6 + (nX1=n1>>30&3)];
			if ((status&TETREF2)!=0) n2=elementIST[(eX2=n2&0x3FFFFFFF)*6 + (nX2=n2>>30&3)];
			if ((status&TETREF3)!=0) n3=elementIST[(eX3=n3&0x3FFFFFFF)*6 + (nX3=n3>>30&3)];
			int n0F = n0 & 0x3FFFFFFF, n1F = n1 & 0x3FFFFFFF, n2F = n2 & 0x3FFFFFFF, n3F = n3 & 0x3FFFFFFF;
			
			// reconfigure pattern & status according to snapnode results, clear internal status of snapped nodes
			if ((status&SET29)!=0) {
				pattern = pattern & (CLR_n0&CLR_n1&CLR_n2&CLR_n3) | (FLG_n0_0|FLG_n1_0|FLG_n2_0|FLG_n3_0);
				status &= CLR_EST0&CLR_EST1&CLR_EST2&CLR_EST3; n0=(n0|SET31)&CLR30; n1=(n1|SET31)&CLR30; n2=(n2|SET31)&CLR30; n3=(n3|SET31)&CLR30;
			} else {
				if (snapNodeStatus[n0F] < 0) { pattern = pattern&CLR_n0|FLG_n0_0; status &= CLR_EST0; n0 = (n0|SET31) & CLR30; }
				if (snapNodeStatus[n1F] < 0) { pattern = pattern&CLR_n1|FLG_n1_0; status &= CLR_EST1; n1 = (n1|SET31) & CLR30; }
				if (snapNodeStatus[n2F] < 0) { pattern = pattern&CLR_n2|FLG_n2_0; status &= CLR_EST2; n2 = (n2|SET31) & CLR30; }
				if (snapNodeStatus[n3F] < 0) { pattern = pattern&CLR_n3|FLG_n3_0; status &= CLR_EST3; n3 = (n3|SET31) & CLR30; }
			}
			e6 -= 6;
			elementIST[e6++] = n0; elementIST[e6++] = n1; elementIST[e6++] = n2; elementIST[e6++] = n3;	
			elementIST[e6++] = pattern; elementIST[e6] = status & TETREF_CLR;

			// if snapped element became internal/external (checked by edges status), no cuts can be made, do nothing more
			if ((status&FLG_EBND) == 0) { e6++; count_becameIntExt++; continue; }
			int eM6 = e * 6, edgeGenFlags = status & 0xFFF, edge = 0;
			
			while (edge < 6) {
				int eGF = edgeGenFlags & 3;
				count_edgetype[eGF]++;
				switch (eGF) {
				case 3: // if unique flag is set & boundary-edge flag is set
					eStatus[0] = 1+2;												// flag cut point calculation & usage of existing snapnodes
					double[] cutNodeO = null;
					switch (edge) {
					case 0:	cutNodeO = processISTedge(geocTree, n0, n1, eStatus, false); pattern|=FLG_CUT01; break;	// edge 01
					case 1:	cutNodeO = processISTedge(geocTree, n0, n2, eStatus, false); pattern|=FLG_CUT02; break;	// edge 02
					case 2:	cutNodeO = processISTedge(geocTree, n0, n3, eStatus, false); pattern|=FLG_CUT03; break;	// edge 03
					case 3:	cutNodeO = processISTedge(geocTree, n1, n2, eStatus, false); pattern|=FLG_CUT12; break;	// edge 12
					case 4:	cutNodeO = processISTedge(geocTree, n1, n3, eStatus, false); pattern|=FLG_CUT13; break;	// edge 13							
					case 5:	cutNodeO = processISTedge(geocTree, n2, n3, eStatus, false); pattern|=FLG_CUT23;			// edge 23
					}
					// note: since an exhaustive edge uniqueness check is not done, "indicate" flag is set so node hash WILL check for node non-uniqueness
					int nIdx = (int)(nodeHTableIST.add(cutNodeO[0], cutNodeO[1], cutNodeO[2], e, true)&0x7FFFFFFF);
					count_cutnodes++;
					// signal was set for adding this cut of an assumed nonshared edge to fast cutnode index search array
					if (nonSharedEdge) {
						nonsharedCutNode[ceIdx] = nIdx;
						ceIdx = (ceIdx-3)&cyclicBufSizeM1; nonSharedEdge = false;
						count_cyclic_misses++; }
					cutNodes[eM6 + edge] = nIdx;
					break;

				// the case of a non-unique boundary edge, meaning that we need to look for nodes referenced to other elements
				// to first find the matching (or assumingly unmatched) edge, and to recover it's attached cutpoint
				case 2:
					int cutIdx = -1;
					long edgeKey = 0;
					// abstract any edge to it's two endpoints: nXa & nXb, and it's two referenced elements: eXa & eXb
					switch (edge) {
					case 0: cutIdx = eX0<eX1 ? cutNodeRef(eX0,eX1,nX0,nX1,n0) : cutNodeRef(eX1,eX0,nX1,nX0,n1); pattern|=FLG_CUT01; break;
					case 1: cutIdx = eX0<eX2 ? cutNodeRef(eX0,eX2,nX0,nX2,n0) : cutNodeRef(eX2,eX0,nX2,nX0,n2); pattern|=FLG_CUT02; break;
					case 2: cutIdx = eX0<eX3 ? cutNodeRef(eX0,eX3,nX0,nX3,n0) : cutNodeRef(eX3,eX0,nX3,nX0,n3); pattern|=FLG_CUT03; break;
					case 3: cutIdx = eX1<eX2 ? cutNodeRef(eX1,eX2,nX1,nX2,n1) : cutNodeRef(eX2,eX1,nX2,nX1,n2); pattern|=FLG_CUT12; break;
					case 4: cutIdx = eX1<eX3 ? cutNodeRef(eX1,eX3,nX1,nX3,n1) : cutNodeRef(eX3,eX1,nX3,nX1,n3); pattern|=FLG_CUT13; break;
					case 5: cutIdx = eX2<eX3 ? cutNodeRef(eX2,eX3,nX2,nX3,n2) : cutNodeRef(eX3,eX2,nX3,nX2,n3); pattern|=FLG_CUT23;
					}
					if (cutIdx < 0) {													// edge references to two elements that don't share an edge
						boolean foundEdge = false;
						switch (edge) {
						case 0: edgeKey = n0F < n1F ? (long)n1F<<32|n0F : (long)n0F<<32|n1F; break;
						case 1: edgeKey = n0F < n2F ? (long)n2F<<32|n0F : (long)n0F<<32|n2F; break;
						case 2: edgeKey = n0F < n3F ? (long)n3F<<32|n0F : (long)n0F<<32|n3F; break;
						case 3: edgeKey = n1F < n2F ? (long)n2F<<32|n1F : (long)n1F<<32|n2F; break;
						case 4: edgeKey = n1F < n3F ? (long)n3F<<32|n1F : (long)n1F<<32|n3F; break;
						case 5: edgeKey = n2F < n3F ? (long)n3F<<32|n2F : (long)n2F<<32|n3F;
						}
						// see if the supposed nonshared edge's cutnode is stored in fast cyclic array, if so, assign to it immediately
						//for (int iStart = (ceIdx+2)&cyclicBufSizeM1, i = iStart; i != ceIdx; i = ++i&cyclicBufSizeM1)
						for (int i = (ceIdx+2)&cyclicBufSizeM1; i != ceIdx; i = ++i&cyclicBufSizeM1)
							if (nonsharedCutNode[i++] == edgeKey) {
								cutNodes[eM6 + edge] = (int)nonsharedCutNode[i]; foundEdge=true;
//								seekAvg += (i - iStart)&cyclicBufSizeM1;				// DEBUG: finds out a fitting size (=power of 2!) for cyclic buffer
//								if (++seekAvgIters >= 16) { System.out.println(seekAvg/16); seekAvgIters = seekAvg = 0; }
								count_cyclic_finds++; break; }

						if (!foundEdge) {												// if cut index not found in fast cyclic array
							nonsharedCutNode[ceIdx++] = edgeKey; nonSharedEdge = true;
							edgeGenFlags|=3;  continue;									// enforce a cut, signal for it's storage in fast array
						} else {
							edge++; edgeGenFlags >>= 2; continue; }						// if found, continue to next edge
					}
					cutNodes[eM6 + edge] = cutNodes[cutIdx];							// edge reference found, store it's precomputed cutnode
					break;
				}
				edge++; edgeGenFlags >>= 2;
			}
			elementIST[--e6] = pattern;	e6 +=2;											// matching pattern has been updated	
		}
		
		// the snap nodes were separately updated to closest possible snap points and need now to replace their counterparts
		printSnapNodes = false;
		int count_snapped = 0;
		for (int n = 0, n3t = 0; n < nodesIST; n++) {
			int n3f = n * 6; if (snapNodeStatus[n] < 0) { n3f += 3; count_snapped++; }
			nodeIST[n3t++] = nodeIST[n3f++]; nodeIST[n3t++] = nodeIST[n3f++]; nodeIST[n3t++] = nodeIST[n3f];
		}
		snapNodeStatus = null;

		debugIST_phase3time = System.nanoTime() - debugIST_phase3time;					// time phase 3

		if (DEBUG_LEVEL > 1) {
			System.out.print("Internalised/externalised elements: " + count_becameIntExt);
			System.out.print("\nTotal nodes: " + (nodesIST + count_cutnodes) + ", of which cutnodes: " + count_cutnodes);
			System.out.print("\nTotal unique snaps: " + count_snapped);
			System.out.print("\nEdge type counts:\n    boundary uniques: " + count_edgetype[3] + "\n    int/ext uniques: " + count_edgetype[1]);
			System.out.print("\n    boundary nonuniques: " + count_edgetype[2] + "\n    int/ext nonuniques: " + count_edgetype[0]);
			System.out.print("\nof those, nonshared edges: " + debugIST_nonsharedEdge);
			System.out.print("\nVisits to FEM1.processISTedge(): " + debugIST_pISTe_visits + ", fast Z-ordinate visits: " + debugIST_pISTe_fastVisits);
			System.out.print(", fast facet visits: " + debugIST_pISTe_fastFacetVisits);
			System.out.print("\nTime spent in FEM1.processISTedge(): " + debugIST_pISTe_time2 + " ns\n");
			System.out.print("Nonshared edges cyclic buffer cases: " + (count_cyclic_misses + count_cyclic_finds));
			System.out.print(" of which directly resolved: " + count_cyclic_finds +"\n\n");
		}
	
		// if zero elements were produced, return zero
		if (count_becameIntExt == elementsIST) {
			if (DEBUG_LEVEL > 1) System.out.println("FEM1.elementGeneratorIST(): No stuffable elements generated, aborting.");
			return enReturn; }

		// time to match the patterns to the permutated stencils and generate the final stencil tetrahedral combos
		// the stencil elements will be added to the internal elements paged array
		// -------------------------------------------------------------------------------------------------------
		if (DEBUG_LEVEL > 1) System.out.println("**************** FEM1.elementGeneratorIST() stage 4 ****************");
		debugIST_phase4time = System.nanoTime();		
		elements = 0;
		elementPg = new int[elementPgN=400][]; elementPg[0] = new int[elementPgSize];	// we will regather elements & stencil elements in elementPg
		elementPgI = 0; ePg = 0;
		
		int[] nodeRemap = new int[nodesIST+nodeHTableIST.count()];						// node indexes remapped to a condensed indexing
		int newNodeIdx = 1;																// the new running node index
		boolean force_ChgA = false;
		int count_ChgAforces=0, writeElements = maxOutputElementsIST > 0 ? maxOutputElementsIST : Integer.MAX_VALUE;
		int[] count_stencil = {0,0,0,0,0,0,0,0};
		SpreadBin angleSpread = null;
		if (DEBUG_LEVEL > 2) angleSpread = new SpreadBin(30, 0, 180, 80);
		
		for (int e = 0, e6 = 0; e < elementsIST && elements < writeElements; e++) {
			
			// stop if element or node target count was exceeded
			if (elementTarget > 0 && elements > elementTarget) {
				if (DEBUG_LEVEL > 1) System.out.println("FEM1.elementGeneratorIST(): Element target exceeded, aborting.");
				enReturn[0] = -elements; return enReturn; }
			if (nodeTarget > 0 && newNodeIdx > nodeTarget) {
				if (DEBUG_LEVEL > 1) System.out.println("FEM1.elementGeneratorIST(): Node target exceeded, aborting.");
				enReturn[1] = -newNodeIdx; return enReturn; }
			
			startPrintElement = e;
			int s0 = elementIST[e6++], s1 = elementIST[e6++], s2 = elementIST[e6++], s3 = elementIST[e6++];
			int n0 = s0&0x3FFFFFFF, n1 = s1&0x3FFFFFFF, n2 = s2&0x3FFFFFFF,n3 = s3&0x3FFFFFFF;
			int pattern = elementIST[e6++], status = elementIST[e6++];
			
			// DEBUG: skip volume-internal gradation elements
			if (skipVolumeInternalsIST && (status&SET30)!=0 && ((s0|s1|s2|s3)&SET31)==0) { force_ChgA = true; continue; }
			if ((status&SET28)!=0) { force_ChgA = true; continue; }								// skip previously discarded
			if ((status&SET29)==0 && ((s0|s1|s2|s3)&SET30)==0) { force_ChgA = true; continue; }	// skip externals, except approved slivers
				
			double Q=0;
			if (DEBUG_LEVEL > 2) {														// do angle spread analysis before stencilling							
				double[] angles = tetraDihedralAngles(nodeIST, n0, n1, n2, n3);
				if ((Q = tetraSmallestDihedral(angles)) < 10) System.out.println("Boundary element " + e + " with too small dihedral angle: " + Q);
				if ((Q = tetraLargestDihedral(angles)) > 170) System.out.println("Boundary element " + e + " with too large dihedral angle: " + Q);
				angleSpread.add(angles[0]); angleSpread.add(angles[1]); angleSpread.add(angles[2]); 
				angleSpread.add(angles[3]); angleSpread.add(angles[4]); angleSpread.add(angles[5]);
			}
			
			int[] elemPage = elementPg[elementPgI];
			if ((status & FLG_EBND) == 0) {								// directly store internal elements (those lacking boundary edges)
				elemPage = verifySpaceElementPg();						// verify that there is space in the paging array
				elemPage[ePg++] = nodeRemap[n0] > 0 ? nodeRemap[n0] : (nodeRemap[n0] = newNodeIdx++);
				elemPage[ePg++] = nodeRemap[n1] > 0 ? nodeRemap[n1] : (nodeRemap[n1] = newNodeIdx++);
				elemPage[ePg++] = nodeRemap[n2] > 0 ? nodeRemap[n2] : (nodeRemap[n2] = newNodeIdx++);
				elemPage[ePg++] = nodeRemap[n3] > 0 ? nodeRemap[n3] : (nodeRemap[n3] = newNodeIdx++);
				// compaction code is still valid (unless previous element has force_ChgA set)
				if (force_ChgA) {
					elemPage[ePg++] = Sl_ChgA; streamNodesIST += 4; force_ChgA = false;
					count_ChgAforces++;
				} else {
					elemPage[ePg] = status>>12 & 0xF;
					switch (elemPage[ePg++] & 0x3FFFFFFF) {
					case Sl_ChgA:								streamNodesIST+=4; break;
					case Sl_Ch01: case Sl_Ch03: case Sl_ShMB:	streamNodesIST+=2; break;
					case Sl_Hld2:								streamNodesIST+=3; break;
					default: streamNodesIST++;
					}
				}
				elements++; continue;
			}
			
			int slotcode = status>>12 & 0xF;
			boolean mirror = (status&SET31) != 0;
			int cn = e*6, seekS=0, seekE=126, seekM=62, countDown = 7;

			while (countDown-- > 0) {
				if (pattern < stencilPermCode[seekM]) seekE = seekM; else if (pattern > stencilPermCode[seekM]) seekS = seekM; else break;
				seekM = seekS + (0xFE&(seekE-seekS)>>1);
			}
			if (countDown==0) throw new RuntimeException("FEM1.elementGeneratorIST(): Couldn't match stencil for element " + e);
			int s = stencilPermSwitch[seekM>>1];
			int idxMix = stencilPermCode[++seekM];

			int cutMix = 0, nL0=0, nL1=0, nL2=0, nP0=0, nP1=0, nP2=0, nP3=0, cn1=0, cn2=0, cn3=0, cn4=0;
			if (s >= 4) {
				// the node-compacted order is lost with the compounded stencils 4/5/6/7/8
				cutMix = idxMix >> 16;								// get the cutnode indexes part
				nL0 = idxMix&15; nL1 = idxMix>>4&15; nL2 = idxMix>>8&15;
				// do the permutation of nodes from aspect of the incoming node order to aspect of basis stencil
				nP0 = nL0==0?n0:(nL1==0?n1:(nL2==0?n2:n3)); nP1 = nL0==1?n0:(nL1==1?n1:(nL2==1?n2:n3));
				nP2 = nL0==2?n0:(nL1==2?n1:(nL2==2?n2:n3)); nP3 = nL0==3?n0:(nL1==3?n1:(nL2==3?n2:n3));
				cn1 = cn+(15&cutMix); cn2 = cn+(15&cutMix>>4); cn3 = cn+(15&cutMix>>8); cn4 = cn+(15&cutMix>>12);
				force_ChgA = true;
			} else {
				// for stencils 1/2/3, snap down the proper nodes to the cutnodes, internalising & finalising the element, order is preserved
				if ((pattern&(1<<2|1<<7|1<<12)) != 0) {
					if ((pattern&1<<2)!=0)	{ if ((pattern&3)==1) n1 = cutNodes[cn]; else n0 = cutNodes[cn]; }
					if ((pattern&1<<7)!=0)	{ if (((pattern>>5)&3)==1) n2 = cutNodes[cn+1]; else n0 = cutNodes[cn+1]; }
					if ((pattern&1<<12)!=0) { if (((pattern>>10)&3)==1) n3 = cutNodes[cn+2]; else n0 = cutNodes[cn+2]; }}
				if ((pattern&(1<<17|1<<22|1<<27)) != 0) {
					if ((pattern&1<<17)!=0)	{ if (((pattern>>15)&3)==1) n2 = cutNodes[cn+3]; else n1 = cutNodes[cn+3]; }
					if ((pattern&1<<22)!=0)	{ if (((pattern>>20)&3)==1) n3 = cutNodes[cn+4]; else n1 = cutNodes[cn+4]; }
					if ((pattern&1<<27)!=0)	{ if (((pattern>>25)&3)==1) n3 = cutNodes[cn+5]; else n2 = cutNodes[cn+5]; }}
			}

			switch (s) {
			case 1: case 2: case 3:	{								// found stencil 1/2/3?
				streamNodesIST += slotcodeSumIST[force_ChgA?Sl_ChgA:slotcode]; elements++;
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				elemPage[ePg++] = nodeRemap[n0] > 0 ? nodeRemap[n0] : (nodeRemap[n0] = newNodeIdx++);
				elemPage[ePg++] = nodeRemap[n1] > 0 ? nodeRemap[n1] : (nodeRemap[n1] = newNodeIdx++);
				elemPage[ePg++] = nodeRemap[n2] > 0 ? nodeRemap[n2] : (nodeRemap[n2] = newNodeIdx++);
				elemPage[ePg++] = nodeRemap[n3] > 0 ? nodeRemap[n3] : (nodeRemap[n3] = newNodeIdx++);
				// DEBUG: edgesharing tetra will also be cut on same edge (unless externalised/internalised), nodestream CAN continue
				elemPage[ePg++] = force_ChgA?Sl_ChgA:slotcode;
				force_ChgA = false;									// only stencils 1/2/3 can continue compaction sequence
				count_stencil[s-1]++; }
				break;
					
			case 4:	{												// found stencil 4?
				streamNodesIST += 4 + 1; elements+=2;
				int c1 = cutNodes[cn1], c2 = cutNodes[cn2];
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				elemPage[ePg++] = nP1 = nodeRemap[nP1] > 0 ? nodeRemap[nP1] : (nodeRemap[nP1] = newNodeIdx++);
				elemPage[ePg++] = c2 = nodeRemap[c2] > 0 ? nodeRemap[c2] : (nodeRemap[c2] = newNodeIdx++);
				elemPage[ePg++] = nP2 = nodeRemap[nP2] > 0 ? nodeRemap[nP2] : (nodeRemap[nP2] = newNodeIdx++);
				elemPage[ePg++] = nodeRemap[c1] > 0 ? nodeRemap[c1] : (nodeRemap[c1] = newNodeIdx++);
				elemPage[ePg++] = Sl_ChgA;
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				elemPage[ePg++] = nodeRemap[nP3] > 0 ? nodeRemap[nP3] : (nodeRemap[nP3] = newNodeIdx++);
				elemPage[ePg++] = c2;
				elemPage[ePg++] = nP2;
				elemPage[ePg++] = nP1;
				elemPage[ePg++] = Sl_FtoL|COMPACT_NEXT_ChgA;		// can be compacted by First-To-Last-mode with previous 4 nodes
				count_stencil[3]++; } break;

			case 5:	{												// found stencil 5?
				streamNodesIST += 4 + 1 + 1; elements+=3;
				int c1 = cutNodes[cn1], c2 = cutNodes[cn2], c3 = cutNodes[cn3], c4 = cutNodes[cn4];
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				elemPage[ePg++] = c3 = nodeRemap[c3] > 0 ? nodeRemap[c3] : (nodeRemap[c3] = newNodeIdx++);
				elemPage[ePg++] = nP2 = nodeRemap[nP2] > 0 ? nodeRemap[nP2] : (nodeRemap[nP2] = newNodeIdx++);
				elemPage[ePg++] = c2 = nodeRemap[c2] > 0 ? nodeRemap[c2] : (nodeRemap[c2] = newNodeIdx++);
				elemPage[ePg++] = nodeRemap[c4] > 0 ? nodeRemap[c4] : (nodeRemap[c4] = newNodeIdx++);
				elemPage[ePg++] = Sl_ChgA;
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				elemPage[ePg++] = nP0 = nodeRemap[nP0] > 0 ? nodeRemap[nP0] : (nodeRemap[nP0] = newNodeIdx++);
				elemPage[ePg++] = nP2;
				elemPage[ePg++] = c2;
				elemPage[ePg++] = c3;
				elemPage[ePg++] = Sl_FtoL;							// can be compacted by First-To-Last-mode with previous 4 nodes
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				elemPage[ePg++] = nodeRemap[c1] > 0 ? nodeRemap[c1] : (nodeRemap[c1] = newNodeIdx++);
				elemPage[ePg++] = nP0;
				elemPage[ePg++] = c2;
				elemPage[ePg++] = c3;
				elemPage[ePg++] = Sl_Mv01|COMPACT_NEXT_ChgA;		// can be compacted by Move-1st-To-2nd with previous 4 nodes
				count_stencil[4]++; } break;
				
			case 50: {												// found stencil 5 mirrored?
				streamNodesIST += 4 + 1 + 1; elements+=3;
				int c1 = cutNodes[cn1], c2 = cutNodes[cn2], c3 = cutNodes[cn3], c4 = cutNodes[cn4];
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				elemPage[ePg++] = nP1 = nodeRemap[nP1] > 0 ? nodeRemap[nP1] : (nodeRemap[nP1] = newNodeIdx++);
				elemPage[ePg++] = nodeRemap[c4] > 0 ? nodeRemap[c4] : (nodeRemap[c4] = newNodeIdx++);
				elemPage[ePg++] = c3 = nodeRemap[c3] > 0 ? nodeRemap[c3] : (nodeRemap[c3] = newNodeIdx++);
				elemPage[ePg++] = c2 = nodeRemap[c2] > 0 ? nodeRemap[c2] : (nodeRemap[c2] = newNodeIdx++);
				elemPage[ePg++] = Sl_ChgA;
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				elemPage[ePg++] = nP0 = nodeRemap[nP0] > 0 ? nodeRemap[nP0] : (nodeRemap[nP0] = newNodeIdx++);
				elemPage[ePg++] = nP1;
				elemPage[ePg++] = c3;
				elemPage[ePg++] = c2;
				elemPage[ePg++] = Sl_Mv01;							// can be compacted by First-To-Last-mode with previous 4 nodes
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				elemPage[ePg++] = nodeRemap[c1] > 0 ? nodeRemap[c1] : (nodeRemap[c1] = newNodeIdx++);
				elemPage[ePg++] = nP0;
				elemPage[ePg++] = c3;
				elemPage[ePg++] = c2;
				elemPage[ePg++] = Sl_Mv01|COMPACT_NEXT_ChgA;		// can be compacted by Move-1st-To-2nd with previous 4 nodes
				count_stencil[4]++; } break;

			case 6: {												// found stencil 6?
				streamNodesIST += 4 + 1; elements+=2;
				int c1 = cutNodes[cn1], c2 = cutNodes[cn2];
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				if (mirror) {										// was parity rule flag set, enforcing a mirror?
					elemPage[ePg++] = nP1 = nodeRemap[nP1] > 0 ? nodeRemap[nP1] : (nodeRemap[nP1] = newNodeIdx++);
					elemPage[ePg++] = nP3 = nodeRemap[nP3] > 0 ? nodeRemap[nP3] : (nodeRemap[nP3] = newNodeIdx++);
					elemPage[ePg++] = c2 = nodeRemap[c2] > 0 ? nodeRemap[c2] : (nodeRemap[c2] = newNodeIdx++);
					elemPage[ePg++] = nodeRemap[c1] > 0 ? nodeRemap[c1] : (nodeRemap[c1] = newNodeIdx++);
					elemPage[ePg++] = Sl_ChgA;
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nodeRemap[nP2] > 0 ? nodeRemap[nP2] : (nodeRemap[nP2] = newNodeIdx++);
					elemPage[ePg++] = nP3;
					elemPage[ePg++] = c2;
					elemPage[ePg++] = nP1;
				} else {
					elemPage[ePg++] = c1 = nodeRemap[c1] > 0 ? nodeRemap[c1] : (nodeRemap[c1] = newNodeIdx++);
					elemPage[ePg++] = nP3 = nodeRemap[nP3] > 0 ? nodeRemap[nP3] : (nodeRemap[nP3] = newNodeIdx++);
					elemPage[ePg++] = nP2 = nodeRemap[nP2] > 0 ? nodeRemap[nP2] : (nodeRemap[nP2] = newNodeIdx++);
					elemPage[ePg++] = nodeRemap[c2] > 0 ? nodeRemap[c2] : (nodeRemap[c2] = newNodeIdx++);
					elemPage[ePg++] = Sl_ChgA;
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nodeRemap[nP1] > 0 ? nodeRemap[nP1] : (nodeRemap[nP1] = newNodeIdx++);
					elemPage[ePg++] = nP3;
					elemPage[ePg++] = nP2;
					elemPage[ePg++] = c1;
				}
				elemPage[ePg++] = Sl_FtoL|COMPACT_NEXT_ChgA;		// can be compacted by First-To-Last-mode with previous 4 nodes
				count_stencil[5]++; } break;
				
			case 7: {												// found stencil 7?
				streamNodesIST += 4 + 1 + 1; elements+=3;
				int c1 = cutNodes[cn1], c2 = cutNodes[cn2], c3 = cutNodes[cn3];
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				if (mirror) {
					elemPage[ePg++] = c2 = nodeRemap[c2] > 0 ? nodeRemap[c2] : (nodeRemap[c2] = newNodeIdx++);
					elemPage[ePg++] = nP1 = nodeRemap[nP1] > 0 ? nodeRemap[nP1] : (nodeRemap[nP1] = newNodeIdx++);
					elemPage[ePg++] = c3 = nodeRemap[c3] > 0 ? nodeRemap[c3] : (nodeRemap[c3] = newNodeIdx++);
					elemPage[ePg++] = nodeRemap[c1] > 0 ? nodeRemap[c1] : (nodeRemap[c1] = newNodeIdx++);
					elemPage[ePg++] = Sl_ChgA;
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nP2 = nodeRemap[nP2] > 0 ? nodeRemap[nP2] : (nodeRemap[nP2] = newNodeIdx++);
					elemPage[ePg++] = nP1;
					elemPage[ePg++] = c3;
					elemPage[ePg++] = c2;
					elemPage[ePg++] = Sl_FtoL;						// can be compacted by First-To-Last-mode with previous 4 nodes
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nodeRemap[nP3] > 0 ? nodeRemap[nP3] : (nodeRemap[nP3] = newNodeIdx++);
					elemPage[ePg++] = nP1;
					elemPage[ePg++] = c3;
					elemPage[ePg++] = nP2;
					elemPage[ePg++] = Sl_FtoL|COMPACT_NEXT_ChgA;	// can be compacted by First-To-Last-mode with previous 4 nodes
				} else {
					elemPage[ePg++] = nP2 = nodeRemap[nP2] > 0 ? nodeRemap[nP2] : (nodeRemap[nP2] = newNodeIdx++);
					elemPage[ePg++] = c1 = nodeRemap[c1] > 0 ? nodeRemap[c1] : (nodeRemap[c1] = newNodeIdx++);
					elemPage[ePg++] = c3 = nodeRemap[c3] > 0 ? nodeRemap[c3] : (nodeRemap[c3] = newNodeIdx++);
					elemPage[ePg++] = nodeRemap[c2] > 0 ? nodeRemap[c2] : (nodeRemap[c2] = newNodeIdx++);
					elemPage[ePg++] = Sl_ChgA;
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nP1 = nodeRemap[nP1] > 0 ? nodeRemap[nP1] : (nodeRemap[nP1] = newNodeIdx++);
					elemPage[ePg++] = c1;
					elemPage[ePg++] = c3;
					elemPage[ePg++] = nP2;
					elemPage[ePg++] = Sl_FtoL;						// can be compacted by First-To-Last-mode with previous 4 nodes
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nodeRemap[nP3] > 0 ? nodeRemap[nP3] : (nodeRemap[nP3] = newNodeIdx++);
					elemPage[ePg++] = nP1;
					elemPage[ePg++] = c3;
					elemPage[ePg++] = nP2;
					elemPage[ePg++] = Sl_Mv01|COMPACT_NEXT_ChgA;	// can be compacted by Move-1st-To-2nd with previous 4 nodes
				}
				count_stencil[6]++; } break;

			case 8: {												// found stencil 8?
				streamNodesIST += 4 + 1 + 1; elements+=3;
				int c1 = cutNodes[cn1], c2 = cutNodes[cn2], c3 = cutNodes[cn3], c4 = cutNodes[cn4];
				elemPage = verifySpaceElementPg();					// verify that there is space in the paging array
				if (mirror) {
					elemPage[ePg++] = c2 = nodeRemap[c2] > 0 ? nodeRemap[c2] : (nodeRemap[c2] = newNodeIdx++);
					elemPage[ePg++] = nP1 = nodeRemap[nP1] > 0 ? nodeRemap[nP1] : (nodeRemap[nP1] = newNodeIdx++);
					elemPage[ePg++] = c3 = nodeRemap[c3] > 0 ? nodeRemap[c3] : (nodeRemap[c3] = newNodeIdx++);
					elemPage[ePg++] = nodeRemap[c1] > 0 ? nodeRemap[c1] : (nodeRemap[c1] = newNodeIdx++);
					elemPage[ePg++] = Sl_ChgA;
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nP2 = nodeRemap[nP2] > 0 ? nodeRemap[nP2] : (nodeRemap[nP2] = newNodeIdx++);
					elemPage[ePg++] = nP1;
					elemPage[ePg++] = c3;
					elemPage[ePg++] = c2;
					elemPage[ePg++] = Sl_FtoL;					// can be compacted by First-To-Last-mode with previous 4 nodes
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nodeRemap[c4] > 0 ? nodeRemap[c4] : (nodeRemap[c4] = newNodeIdx++);
					elemPage[ePg++] = nP2;
					elemPage[ePg++] = c3;
					elemPage[ePg++] = c2;
					elemPage[ePg++] = Sl_Mv01|COMPACT_NEXT_ChgA;	// can be compacted by First-To-Last-mode with previous 4 nodes
				} else {
					elemPage[ePg++] = nP2 = nodeRemap[nP2] > 0 ? nodeRemap[nP2] : (nodeRemap[nP2] = newNodeIdx++);
					elemPage[ePg++] = c1 = nodeRemap[c1] > 0 ? nodeRemap[c1] : (nodeRemap[c1] = newNodeIdx++);
					elemPage[ePg++] = c4 = nodeRemap[c4] > 0 ? nodeRemap[c4] : (nodeRemap[c4] = newNodeIdx++);
					elemPage[ePg++] = nodeRemap[c2] > 0 ? nodeRemap[c2] : (nodeRemap[c2] = newNodeIdx++);
					elemPage[ePg++] = Sl_ChgA;
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nP1 = nodeRemap[nP1] > 0 ? nodeRemap[nP1] : (nodeRemap[nP1] = newNodeIdx++);
					elemPage[ePg++] = c1;
					elemPage[ePg++] = c4;
					elemPage[ePg++] = nP2;
					elemPage[ePg++] = Sl_FtoL;						// can be compacted by First-To-Last-mode with previous 4 nodes
					elemPage = verifySpaceElementPg();				// verify that there is space in the paging array
					elemPage[ePg++] = nodeRemap[c3] > 0 ? nodeRemap[c3] : (nodeRemap[c3] = newNodeIdx++);
					elemPage[ePg++] = c1;
					elemPage[ePg++] = c4;
					elemPage[ePg++] = nP1;
					elemPage[ePg++] = Sl_FtoL|COMPACT_NEXT_ChgA;	// can be compacted by First-To-Last-mode with previous 4 nodes
				}
				count_stencil[7]++; } break;
			}
			if (force_ChgA) count_ChgAforces++;
		}
		cutNodes = null;											// cut nodes only needed for stencil generation, kill array
		elementIST = null;											// kill the boundary elements gathering array
		
		// recombine snapped nodes & cutnodes into a single condensed array, unused nodes eliminated			
		double[] nodeIST1 = nodeHTableIST.array(0, nodeRemap, newNodeIdx);
		for (int i=0, i3f = 0; i < nodesIST; i++) {
			int i3t = nodeRemap[i];
			if (i3t == 0) { i3f += 3; continue; } else i3t *= 3;
			nodeIST1[i3t++] = nodeIST[i3f++]; nodeIST1[i3t++] = nodeIST[i3f++]; nodeIST1[i3t] = nodeIST[i3f++];
		}
		nodeRemap = null;
		nodesIST = newNodeIdx;
		nodeIST = nodeIST1;
		nodeHTableIST = null;
		debugIST_phase4time = System.nanoTime() - debugIST_phase4time;		// time phase 4

		if (DEBUG_LEVEL > 1) {
			System.out.print(elements + " final tetrahedra.\n");
			System.out.print("Condensed nodes final count: " + newNodeIdx + "\nStencils encountered:\n");
			for (int s = 0; s < 8; s++) System.out.println("Stencil " + s + ": " + count_stencil[s]);
			System.out.println("\nCompaction, forced all-node changes: " + count_ChgAforces);
			if (DEBUG_LEVEL > 2) System.out.print("Angle distribution prior to stencilling:\n" + angleSpread.toString() + "\n\n");
			else System.out.print("\n\n");
			System.out.println("**************** FEM1.elementGeneratorIST() stage 5 ****************");
		}
		
		debugIST_phase5time = System.nanoTime();							// start time phase 5
		readElementPg = finalArrayPg = true;
		startPrintElement=0; //printSlotcodes1 = true;
		int[] count_foundCompactions = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		// time to do the final collection of nodes into a single compact stream
		// ---------------------------------------------------------------------
		elementStream = new double[streamNodesIST * 3];
		elementStreamCode = new long[elements / 16 + 1];
		int[] elemPage = elementPg[0];
		int n0=0, n1=0, n2=0, n3=0, eS = 0, n03, n13, n23, n33, eSC = 0;
		
		// preload first element for one-ahead comparison
		n0=elemPage[0]; n1=elemPage[1]; n2=elemPage[2]; n3=elemPage[3]; 
		elementStream[eS++] = nodeIST[n03=n0*3]; elementStream[eS++] = nodeIST[++n03]; elementStream[eS++] = nodeIST[++n03];	// store first element
		elementStream[eS++] = nodeIST[n13=n1*3]; elementStream[eS++] = nodeIST[++n13]; elementStream[eS++] = nodeIST[++n13];
		elementStream[eS++] = nodeIST[n23=n2*3]; elementStream[eS++] = nodeIST[++n23]; elementStream[eS++] = nodeIST[++n23];
		elementStream[eS++] = nodeIST[n33=n3*3]; elementStream[eS++] = nodeIST[++n33]; elementStream[eS++] = nodeIST[++n33];
		long slotcodes = 0, slotCnt = 1;
		
		for (int e = 1, ePage = 0, e5 = 5; e < elements; e++, slotCnt++) {
			startPrintElement = e;
			if (e5 >= elementPgSize) { e5 = 0; elemPage = elementPg[++ePage]; }			// go to next page of elements if last one was read
			if (slotCnt >= 16) { elementStreamCode[eSC++] = slotcodes; slotcodes = slotCnt = 0; }
			
			int n0N = elemPage[e5++], n1N = elemPage[e5++], n2N = elemPage[e5++], n3N = elemPage[e5++], slotcode = elemPage[e5++];
			
			// on changing all nodes, do a check if new nodes match up according to SOME slot coding scheme, employ scheme if they do
			slotcode &= 0xF;
			//if (slotcode == Sl_ChgA) {
			int oldSlotcode = slotcode;
				if (n1N==n1 && n2N==n2)		{
									if (n3N==n0)		slotcode = Sl_FtoL;				// found a First-to-Last transition?
									else 				slotcode = Sl_Ch03;	}			// found Change-3rd transition?
				else if (n2N==n2 && n3N==n3) {
									if (n1N==n0)		slotcode = Sl_Mv01;				// found Move 0 to 1?
									else if (n0N==n0)	slotcode = Sl_Chg1;				// found Change-1st transition?
									else 				slotcode = Sl_Ch01; }			// found Change 0 & 1?
				else if (n0N==n0 && n1N==n1 && n2N==n2)	slotcode = Sl_Chg3;				// found Change-3rd transition?
				else if (n0N==n0 && n1N==n1 && n3N==n3)	slotcode = Sl_Chg2;				// found Change-2nd transition?
				else if (n2N==n1 && n3N==n2)			slotcode = Sl_ShMB;				// found Shift mid-two to back transition?
				else if (n1N==n1 && n3N==n3)			slotcode = Sl_Ch02;				// found Change 0 & 2?
				else if (n2N==n2)						slotcode = Sl_Hld2;				// found Hold-2nd transition?
				else 									slotcode = Sl_ChgA;
				if (oldSlotcode != slotcode) count_foundCompactions[slotcode]++;
			//}

			switch (slotcode) {
			case Sl_FtoL:
			case Sl_Mv01:	elementStream[eS++] = nodeIST[n03=n0N*3]; elementStream[eS++] = nodeIST[++n03]; elementStream[eS++] = nodeIST[++n03]; break;
			case Sl_ShMB:
			case Sl_Ch01:	elementStream[eS++] = nodeIST[n03=n0N*3]; elementStream[eS++] = nodeIST[++n03]; elementStream[eS++] = nodeIST[++n03];
			case Sl_Chg1: 	elementStream[eS++] = nodeIST[n13=n1N*3]; elementStream[eS++] = nodeIST[++n13]; elementStream[eS++] = nodeIST[++n13]; break;
			case Sl_Ch02:	elementStream[eS++] = nodeIST[n03=n0N*3]; elementStream[eS++] = nodeIST[++n03]; elementStream[eS++] = nodeIST[++n03];
			case Sl_Chg2:	elementStream[eS++] = nodeIST[n23=n2N*3]; elementStream[eS++] = nodeIST[++n23]; elementStream[eS++] = nodeIST[++n23]; break;
			case Sl_ChgA: case Sl_Ch03:
			case Sl_Hld2:	elementStream[eS++] = nodeIST[n03=n0N*3]; elementStream[eS++] = nodeIST[++n03]; elementStream[eS++] = nodeIST[++n03];
						if (slotcode==Sl_Hld2 || slotcode==Sl_ChgA) {
							elementStream[eS++] = nodeIST[n13=n1N*3]; elementStream[eS++] = nodeIST[++n13]; elementStream[eS++] = nodeIST[++n13]; }
						if (slotcode==Sl_ChgA) {
							elementStream[eS++] = nodeIST[n23=n2N*3]; elementStream[eS++] = nodeIST[++n23]; elementStream[eS++] = nodeIST[++n23]; }
			case Sl_Chg3:	elementStream[eS++] = nodeIST[n33=n3N*3]; elementStream[eS++] = nodeIST[++n33]; elementStream[eS++] = nodeIST[++n33]; break;
			}
			n0 = n0N; n1 = n1N; n2 = n2N; n3 = n3N;
			slotcodes |= (long)slotcode << (slotCnt*4);
		}
		if (slotcodes != 0) elementStreamCode[eSC++] = slotcodes;
		
		startPrintElement = 0; printSlotcodes1 = false; printSlotcodes2 = true;
		int oldStreamNodesIST = streamNodesIST;
		streamNodesIST = eS / 3;													// the true final count of streamed nodes	
		element = new int[elements * 4];											// gather elements into a final index-quartet array
		elemPage = elementPg[0];
				
		for (int e = 0, ePage = 0, e4 = 0, e5 = 0; e < elements; e++) {
			startPrintElement=e;
			if (e5 >= elementPgSize) { e5 = 0; elemPage = elementPg[++ePage]; }		// go to next page if last element was read
			element[e4++] = elemPage[e5++]; element[e4++] = elemPage[e5++];
			element[e4++] = elemPage[e5++]; element[e4++] = elemPage[e5++];
			e5++;
		}
		debugIST_phase5time = System.nanoTime() - debugIST_phase5time;				// time phase 5

		if (DEBUG_LEVEL > 1) {
			System.out.print(oldStreamNodesIST + " premature streamnodes, " + streamNodesIST + " final streamnodes\nStream compactions found:\n");
			for (int s = 1; s < 16; s++)
				if (count_foundCompactions[s] != 0) System.out.print(slotcodeName[s] + ": " + count_foundCompactions[s] + (s==15 ? "\n\n" : "\n"));
			System.out.print("Total stage timings:\nBndOcts: " + debugIST_aBOIST_time + " ns\n");
			System.out.print("Lattice: " + debugIST_bLatt_time + " ns\n");
			System.out.print("Stage 1: " + debugIST_phase1time + " ns\n");
			System.out.print("Stage 2: "+debugIST_phase2time+" ns (cut function): "+debugIST_pISTe_time1+" ns, (quality analysis): "+debugIST_analysisTime+" ns\n");
			System.out.print("Stage 3: " + debugIST_phase3time + " ns, (cut function): " + debugIST_pISTe_time2 + " ns\nStage 4: " + 
								debugIST_phase4time + " ns\nStage 5: " + debugIST_phase5time + " ns\n");
			double timeSum = debugIST_aBOIST_time + debugIST_bLatt_time + 
					debugIST_phase1time + debugIST_phase2time + debugIST_phase3time + debugIST_phase4time + debugIST_phase5time;
			System.out.print("Percentage time:\nBndOcts: " + String.format("%.1f", 100*(double)debugIST_aBOIST_time/timeSum));
			System.out.print("\nLattice: " + String.format("%.1f", 100*(double)debugIST_bLatt_time/timeSum)+"%");
			System.out.print("\nStage 1: " + String.format("%.1f", 100*(double)debugIST_phase1time/timeSum)+"%");
			System.out.print("\nStage 2: " + String.format("%.1f", 100*(double)debugIST_phase2time/timeSum)+"%");
			System.out.print(" (cut function): " + String.format("%.1f", 100*(double)debugIST_pISTe_time1/timeSum)+"%");
			System.out.print(" (quality analysis): " + String.format("%.1f", 100*(double)debugIST_analysisTime/timeSum)+"%");
			System.out.print("\nStage 3: " + String.format("%.1f", 100*(double)debugIST_phase3time/timeSum)+"%");
			System.out.print(" (cut function): " + String.format("%.1f", 100*(double)debugIST_pISTe_time2/timeSum)+"%");
			System.out.print("\nStage 4: " + String.format("%.1f", 100*(double)debugIST_phase4time/timeSum)+"%");
			System.out.print("\nStage 5: " + String.format("%.1f", 100*(double)debugIST_phase5time/timeSum) + "%\n\n");
		}
		
		if (DEBUG_LEVEL > 1) {
			count_defectives = 0;
			angleSpread = new SpreadBin(30, 0, 180, 80);
			for (int e = 0, e4 = 0; e < elements; e++) {
				n0 = element[e4++]; n1 = element[e4++]; n2 = element[e4++]; n3 = element[e4++];
				double Q = 0;
				if ((Q = tetraQuality(nodeIST, n0, n1, n2, n3)) < minElemQuality) {
					System.out.println("Element " + e + " has bad quality: " + Q);
					if (DEBUG_LEVEL > 2) {
						System.out.print(String.format  ("      n0: %.3f,%.3f,%.3f, ", nodeIST[n0*3], nodeIST[n0*3+1], nodeIST[n0*3+2]));
						System.out.print(String.format  ("      n1: %.3f,%.3f,%.3f, ", nodeIST[n1*3], nodeIST[n1*3+1], nodeIST[n1*3+2]));
						System.out.print(String.format  ("      n2: %.3f,%.3f,%.3f, ", nodeIST[n2*3], nodeIST[n2*3+1], nodeIST[n2*3+2]));
						System.out.println(String.format("      n3: %.3f,%.3f,%.3f, ", nodeIST[n3*3], nodeIST[n3*3+1], nodeIST[n3*3+2]));
					}
					count_defectives++;
				}				
				double[] angles = tetraDihedralAngles(nodeIST, n0, n1, n2, n3);
				if ((Q = tetraSmallestDihedral(angles)) < 10) System.out.println("Boundary element " + e + " with too small dihedral angle: " + Q);
				if ((Q = tetraLargestDihedral(angles)) > 170) System.out.println("Boundary element " + e + " with too large dihedral angle: " + Q);
				angleSpread.add(angles[0]); angleSpread.add(angles[1]); angleSpread.add(angles[2]); 
				angleSpread.add(angles[3]); angleSpread.add(angles[4]); angleSpread.add(angles[5]);
			}
			if (count_defectives > 0) System.out.print("Bad quality tetrahedra: " + count_defectives);
			System.out.println("\nAngle distribution, final:\n" + angleSpread.toString() + "\n\n");
		}
		printSlotcodes2 = false;
		enReturn[0] = elements; enReturn[1] = nodesIST; return enReturn;
	}

	
	// methods abstract any edge to it's two endpoints: nXa & nXb, and it's two referenced elements: eXa (farthest) & eXb (nearest), and finds out
	// which node within the nearest element of endpoint nXb is equal to the node na in the tested element,
	// returning the cut node position pertinent to that edge within the nearest element
	// note: it can in fact be that two node refs to two elements exist but no shared edge exist, method returns -1 if so
	private int cutNodeRef(int eXa, int eXb, int nXa, int nXb, int na) {

		if (eXa==eXb) {												// the simple case of both endpoints referring to one single neighbour tetra-element
			return eXa*6 + (nXa==0||nXb==0?nXa+nXb-1:nXa+nXb);
		} else {
			// we have situation where this element refers to a neighbour who in turn has a node referral to same neighbour as this element
			int eXb6 = eXb * 6;
			switch (nXb) {
			case 0:	// check that on neighbour's edge 01/02/03, if it's node ref. = this node ref.
				if (elementIST[eXb6+1]==na) return eXb*6;			// nbrs e01 cut node
				if (elementIST[eXb6+2]==na) return eXb*6+1;			// nbrs e02 cut node
				if (elementIST[eXb6+3]==na) return eXb*6+2;	break;	// nbrs e02 cut node
			case 1:	// check that on neighbour's edge 01/12/13, if it's  the node ref. = this node ref.
				if (elementIST[eXb6]==na)   return eXb*6;			// nbrs e01 cut node
				if (elementIST[eXb6+2]==na) return eXb*6+3;			// nbrs e12 cut node
				if (elementIST[eXb6+3]==na) return eXb*6+4;	break;	// nbrs e13 cut node
			case 2:	// check that on neighbour's edge 02/12/23, if it's  the node ref. = this node ref.
				if (elementIST[eXb6]==na)   return eXb*6+1;			// nbrs e02 cut node
				if (elementIST[eXb6+1]==na) return eXb*6+3;			// nbrs e12 cut node
				if (elementIST[eXb6+3]==na) return eXb*6+5;	break;	// nbrs e23 cut node
			case 3:	// check that on neighbour's edge 03/13/23, if it's  the node ref. = this node ref.
				if (elementIST[eXb6]==na)   return eXb*6+2;			// nbrs e03 cut node
				if (elementIST[eXb6+1]==na) return eXb*6+4;			// nbrs e13 cut node
				if (elementIST[eXb6+2]==na) return eXb*6+5;			// nbrs e23 cut node
			}
		}
		debugIST_nonsharedEdge++;
		return -1;
	}
	
	// method returns 1 if transient edge reference exists, operates with the references stored in the elements, rather than node indexes
	private int transitiveRef(int eXa, int eXb, int nXa, int nXb) {
		// we have situation where this element refers to a neighbour who in turn has a node referral to same neighbour as this element
		int eXb6 = eXb * 6, ngbS = elementIST[eXb6 + 5], refa = nXa<<30|eXa;
		switch (nXb) {
		case 0:	// check that on neighbour's edge 01/02/03, if it's node ref. = this node ref.
			if ((ngbS&TETREF1)!=0 && elementIST[eXb6+1]==refa) return 1;
			if ((ngbS&TETREF2)!=0 && elementIST[eXb6+2]==refa) return 1;
			if ((ngbS&TETREF3)!=0 && elementIST[eXb6+3]==refa) return 1;
		case 1:	// check that on neighbour's edge 01/12/13, if it's  the node ref. = this node ref.
			if ((ngbS&TETREF0)!=0 && elementIST[eXb6]==refa)   return 1;
			if ((ngbS&TETREF2)!=0 && elementIST[eXb6+2]==refa) return 1;
			if ((ngbS&TETREF3)!=0 && elementIST[eXb6+3]==refa) return 1;
		case 2:	// check that on neighbour's edge 02/12/23, if it's  the node ref. = this node ref.
			if ((ngbS&TETREF0)!=0 && elementIST[eXb6]==refa)   return 1;
			if ((ngbS&TETREF1)!=0 && elementIST[eXb6+1]==refa) return 1;
			if ((ngbS&TETREF3)!=0 && elementIST[eXb6+3]==refa) return 1;
		case 3:	// check that on neighbour's edge 03/13/23, if it's  the node ref. = this node ref.
			if ((ngbS&TETREF0)!=0 && elementIST[eXb6]==refa)   return 1;
			if ((ngbS&TETREF1)!=0 && elementIST[eXb6+1]==refa) return 1;
			if ((ngbS&TETREF2)!=0 && elementIST[eXb6+2]==refa) return 1;
		}
		debugIST_nonsharedEdge2++;
		return -1;
	}

	
	// checks int_ext bits pertinent to the formed tetrahedron, if any one is internal, the method fails
	private boolean externalTetrahedron(FEM1Octant o, FEM1Octant  n, int seqP) {
		int inExT = 0, c = 4, inEx = o.int_ext;
		while (c-- > 0) {
			int code = faceSeq_IST[seqP++];
			switch (code) {
			case Cd_MMM: inExT|=inEx; break;		case Cd_PMM: inExT|=inEx>>1; break;
			case Cd_MPM: inExT|=inEx>>2; break;		case Cd_PPM: inExT|=inEx>>3; break;
			case Cd_MMP: inExT|=inEx>>4; break;		case Cd_PMP: inExT|=inEx>>5; break;
			case Cd_MPP: inExT|=inEx>>6; break;		case Cd_PPP: inExT|=inEx>>7; break;
			case Cd_CMM: inExT|=inEx>>8; break;		case Cd_CPM: inExT|=inEx>>12; break;
			case Cd_CMP: inExT|=inEx>>22; break;	case Cd_CPP: inExT|=inEx>>26; break;
			case Cd_MCM: inExT|=inEx>>9; break;		case Cd_PCM: inExT|=inEx>>11; break;
			case Cd_MCP: inExT|=inEx>>23; break;	case Cd_PCP: inExT|=inEx>>25; break;
			case Cd_MMC: inExT|=inEx>>13; break;	case Cd_PMC: inExT|=inEx>>15; break;
			case Cd_MPC: inExT|=inEx>>19; break;	case Cd_PPC: inExT|=inEx>>21; break;
			case Cd_CCC: inExT|=inEx>>17; break;	case Cd_CCCn:inExT|=n.int_ext>>17; break;
			case Cd_MCC: inExT|=inEx>>16; break;	case Cd_PCC: inExT|=inEx>>18; break;
			case Cd_CMC: inExT|=inEx>>14; break;	case Cd_CPC: inExT|=inEx>>20; break;
			case Cd_CCM: inExT|=inEx>>10; break;	case Cd_CCP: inExT|=inEx>>24; break;
			default:	throw new RuntimeException("FEM1.processISTnode(): Invalid aspect code received."); }
			if ((inExT&1) != 0) return false;
		}
		return true;
	}
	

	private double lStepI_pISTe = 0, lStepI2_pISTe = 0, intOfsX_pISTe = 0, intOfsY_pISTe = 0;
	private int divCofs_pISTe = 0;

	// method does the boundary intersection test on a "short" or "long" IST edge, if any end of edge is within criterion closeness to the
	// boundary, that node-end of edge is snapped to voundary and turned into a 0-node, and if no end fulfills the snapping criterion
	// the resulting cut node is hash-stored together with the other unique IST nodes, it's approx. distance & status stored in snapNodeStatus[] array
	// note: method works on the supplied geometry octree that carries the boundary facets
	// TODO: method could benefit from a simple salvaging strategy for non-watertight meshes, since the boundaryLattice() method has one
	private double[] processISTedge(FEM1Octree geocTree, int na, int nb, int[] statusE, boolean longEdge) {

		debugIST_pISTe_visits++;
		long timerS = DEBUG_LEVEL > 1 ? System.nanoTime() : 0;
		int naF = na & 0x3FFFFFFF, nbF = nb & 0x3FFFFFFF, n3a = naF * 6, n3b = nbF * 6, snapa = snapNodeStatus[naF], snapb = snapNodeStatus[nbF];
		int status = statusE[0];
		// test case of both edges being snapped within an "good enough" margin 
		if (status == 0 && (snapa & snapb & 0x80) !=0) {
			if ((snapa&127)<=15 && (snapb&127)<=15) { statusE[0] = 3; debugIST_pISTe_skipVisits++; return null; }
		}
		
		/***** Cut function block starts ******************************************/
		// cut function must return the intersection point (or null) in the
		// cutNodeO[] array, naF & nbF give indexes into the nArray[] node array
		
		double[] cutNodeO = null;
		int[] statusG = {-1, 0, -1};
		// do intersection test within the geometry octree
		double nax=0, nay=0, naz=0, nbx=0, nby=0, nbz=0;
		if ((status&2)==2) {							// flag in status integer tells to use EXISTING snap nodes instead of original nodes
			if (snapa < 0) n3a += 3; nax = nodeIST[n3a++]; nay = nodeIST[n3a++]; naz = nodeIST[n3a];
			if (snapb < 0) n3b += 3; nbx = nodeIST[n3b++]; nby = nodeIST[n3b++]; nbz = nodeIST[n3b];
		} else {
			nax = nodeIST[n3a++]; nay = nodeIST[n3a++]; naz = nodeIST[n3a]; nbx = nodeIST[n3b++]; nby = nodeIST[n3b++]; nbz = nodeIST[n3b]; }
		// this edge-shrinking code expects to take care of certain topological defects when a "straddling" element gets snapped across empty space
		// if that happens, a defective cutnode might return, right at the cusp where the edge got snapped onto the surface, this can be solved
		// by shrinking the edge somewhat inward, and is ONLY done for stage 3, the cutnode search iterator
		if ((status&1)==1) {
			double shrinkx = (nax-nbx)*maxZordinDev, shrinky = (nay-nby)*maxZordinDev, shrinkz = (naz-nbz)*maxZordinDev;
			nax -= shrinkx; nay -= shrinky; naz -= shrinkz; nbx += shrinkx; nby += shrinky; nbz += shrinkz;
		}
		double xdba = nbx<nax ? nax-nbx : nbx-nax, ydba = nby<nay ? nay-nby : nby-nay, zdba = nbz<naz ? naz-nbz : nbz-naz;
		
		// since we've already calculated volume intersections along Z ordinate, for Z ordinate aligned edges that data will be used,
		// IF the edge is (approximately) aligned with Z-ordinate
		if (xdba<maxZordinDev && ydba<maxZordinDev) {
			int inx = (int)((nax + intOfsX_pISTe) * lStepI2_pISTe), iny = (int)((nay + intOfsY_pISTe) * lStepI_pISTe);	// convert X,Y to array offsets
			iny = (iny&1)==1 ? (iny>>1) + divCofs_pISTe : iny>>1;				// if edge goes through centroid ray, offset to centroid rays data part
			int enc3 = fStatusGrid[iny][inx] * NCOORD;							// how many values (of coord triplets) we need to iterate over
			cutNodeO = new double[3];
			double[] fEncounterXY = fEncounterGrid[iny][inx];
			double z1, z2; if (naz < nbz) { z1 = naz; z2 = nbz; } else { z2 = naz; z1 = nbz; }
			for (int e3 = 2; e3 < enc3; e3+=3) {
				double zEnc = fEncounterXY[e3];
				if (z1 < zEnc && zEnc < z2) { cutNodeO[0] = nax; cutNodeO[1] = nay; cutNodeO[2] = zEnc; break; }}
			debugIST_pISTe_fastVisits++;
		} else if (statusE[1] >= 0) {											// do we have a fast facet stored away?
			cutNodeO = new double[3];
			if (geocTree.fem.facetSegmentIntersection(statusE[1], nax, nay, naz, nbx, nby, nbz, cutNodeO) != 0) {
				debugIST_pISTe_fastFacetVisits++;
			} else {
				geocTree.fem.polygonCheck.clearVisits();						// fast facet failed to intersect, fallback to facetEncounters()
				geocTree.fem.polygonCheck.visit(statusE[1]);
				cutNodeO = geocTree.root.facetEncounters(geocTree, nax, nay, naz, nbx, nby, nbz, 3+16+32, statusG);
				statusE[1] = statusG[2]; }
		} else {
			cutNodeO = geocTree.root.facetEncounters(geocTree, nax, nay, naz, nbx, nby, nbz, 3+32, statusG);
			statusE[1] = statusG[2];
		}
		// handle numeric border cases (node exactly at a facet plane, etc) by miniscule expansion of edge length
//		if (statusG[0] == 0) {
//			double growx = (nax-nbx)*1e-9, growy = (nay-nby)*1e-9, growz = (naz-nbz)*1e-9;
//			nax += growx; nay += growy; naz += growz; nbx -= growx; nby -= growy; nbz -= growz;
//			cutNodeO = geocTree.root.facetEncounters(geocTree, nax, nay, naz, nbx, nby, nbz, 3+32, statusG);
//			statusE[1] = statusG[2];
//		}
		
		/***** Cut function block ends ********************************************/
		
		if ((status&1) == 1) {
			if (DEBUG_LEVEL > 1) debugIST_pISTe_time2 += System.nanoTime() - timerS;
			return cutNodeO; }											// return cut point if request was only for it's calculation
		
		double xdca = cutNodeO[0] - nax, ydca = cutNodeO[1] - nay, zdca = cutNodeO[2] - naz, eL=0, cL=0;
		// test violation criterion from endpoint "a" of the edge, do direct length comparison if edge lies along x/y/z-ordinate
		if (xdba < maxZordinDev) {
			if (ydba<maxZordinDev) eL = zdba; else if (zdba<maxZordinDev) eL = ydba;
			else eL = Math.sqrt(xdba*xdba + ydba*ydba + zdba*zdba);
		} else if (ydba < maxZordinDev) {
			if (zdba<maxZordinDev) eL = xdba;
			else eL = Math.sqrt(xdba*xdba + ydba*ydba + zdba*zdba); }
		else eL = Math.sqrt(xdba*xdba + ydba*ydba + zdba*zdba);			// this is the edge length
		
		if (xdca < maxZordinDev && xdca > -maxZordinDev) {
			if (ydca<maxZordinDev && ydca>-maxZordinDev) cL = zdca; else if (zdca<maxZordinDev && zdca>-maxZordinDev) cL = ydca;
			else cL = Math.sqrt(xdca*xdca + ydca*ydca + zdca*zdca);
		} else if (ydca < maxZordinDev && ydca > -maxZordinDev) {
			if (zdca<maxZordinDev && zdca>-maxZordinDev) cL = xdca;
			else cL = Math.sqrt(xdca*xdca + ydca*ydca + zdca*zdca); }
		else cL = Math.sqrt(xdca*xdca + ydca*ydca + zdca*zdca);			// this is the length from end a to the snap/cutpoint c
		if (eL < 0) eL = -eL; if (cL < 0) cL = -cL; 
		double eLf = eL * (longEdge ? alphaL : alphaS);
		
		if (cL < eLf) {													// cut node violated node "a" of edge?
			byte factor = (byte)((cL * 127) / eL);						// compress violation factor to 7 bits of precision (quite enough for our needs)
			if (factor < (snapa & 127)) {								// is it closer than previous snap distance?
				nodeIST[n3a=naF*6+3] = cutNodeO[0]; nodeIST[++n3a] = cutNodeO[1]; nodeIST[++n3a] = cutNodeO[2];
				snapNodeStatus[naF] = (byte)(0x80|factor);				// will flag change to snapnode
			} else debugIST_pISTe_vainSnaps++;
			statusE[0] = 1;												// return that we snapped to node a
			//if ((na&SET30)==0) snapNodeStatus[nbF]&=127;				// a snap of an external node obliterates a snap of an internal node on other end
			if (DEBUG_LEVEL > 1) debugIST_pISTe_time2 += System.nanoTime() - timerS;
			return cutNodeO;
			
		} else if ((eL - cL) < eLf) {									// cut node violated node "b" of edge?
			byte factor = (byte)(((eL-cL) * 127) / eL);
			if (factor < (snapb & 127)) {								// is it closer than previous snap distance?
				nodeIST[n3b=nbF*6+3] = cutNodeO[0]; nodeIST[++n3b] = cutNodeO[1]; nodeIST[++n3b] = cutNodeO[2];
				snapNodeStatus[nbF] = (byte)(0x80|factor);				// will flag change to snapnode
			} else debugIST_pISTe_vainSnaps++;
			statusE[0] = 2;												// return that we snapped to node b
			//if ((nb&SET30)==0) snapNodeStatus[naF]&=127;				// a snap of an external node obliterates a snap of an internal node on other end
			if (DEBUG_LEVEL > 1) debugIST_pISTe_time2 += System.nanoTime() - timerS;
			return cutNodeO;
		}
		
		statusE[0] = 3;													// return that neither end was violated by cutnode
		if (DEBUG_LEVEL > 1) debugIST_pISTe_time2 += System.nanoTime() - timerS;
		return cutNodeO;
	}
	
	
	
	final static int[] nodeVisitEnumF = { 1,1<<2,1<<6,1<<8,1<<18,1<<20,1<<24,1<<26,
			1<<13,1<<1,1<<7,1<<19,1<<25,1<<3,1<<5,1<<21,1<<23,1<<9,1<<11,1<<15,1<<17,1<<12,1<<14,1<<10,1<<16,1<<4,1<<22,0/*CCCn*/};

	// helping method for Isosurface Stuffing method backgroundTetraGridBCC() extracts nodes from octant according to a supplied code scheme
	// the scheme aids in compacting the nodes later (avoiding repititions) but also avoids reading them more than once
	// method registers the running internal/external coordinate status in the intExt bitfield
	private void processISTnode(FEM1Octant o, FEM1Octant n, int codes) {
		
		long timerS = DEBUG_LEVEL > 1 ? System.nanoTime() : 0;
		byte slotcode = (byte)(codes & 0x7F), aspectCode = (byte)((codes>>8) & 0xFF);
		int edgeCode = codes & 0x3F0000, seq2 = 2 * (codes>>24 & 7);
		int mirrorFlag = (codes & 0x80) != 0 ? SET31 : 0;
		// flush previous tetrahedron either to the internal or the boundary element arrays
		if (slotcode < 8) {
			switch (slotcodeIST) {													// need to mark the  vertexes that remained in place as nonunique
			case Sl_FtoL: tetIST1 |= SET31; tetIST2 |= SET31; tetIST3 |= SET31; break;
			case Sl_Hld2: tetIST2 |= SET31; break;
			case Sl_Ch01: tetIST2 |= SET31; tetIST3 |= SET31; break;
			case Sl_Mv01: tetIST1 |= SET31; tetIST2 |= SET31; tetIST3 |= SET31; break;
			case Sl_Ch03: tetIST1 |= SET31; tetIST2 |= SET31; break;
			}
			int tF0=tetIST0>>30&3, tF1=tetIST1>>30&3, tF2=tetIST2>>30&3, tF3=tetIST3>>30&3;

			if (elements < 0) elements++;											// avoid writing an empty element by a 1st increment of "elements"
			// an internal element has, per entry: 1 global order index, 4 node indexes, 1 int for the slot code
			
			else if (((tF0|tF1|tF2|tF3)&1) != 0) {									// if at least one node is internal
				// the boundary element array is special, having per entry:
				// 1 global order index, 4 node indexes/element refs, 1 pattern for stencil matching, 1 edge status & slotcode & noderef. status int
				int[] elemPage = verifySpaceElementPg();

				int edgeFlags=0;
				// if a node is nonunique, flag node index as being node's original element index (bits 22-25)
				switch (tF0) {
				case 0:	edgeFlags=21;										elemPage[ePg++]=tetIST0&CLR31; break;	// if node is external & unique
				case 1:	edgeFlags=21+FLG_EST0;								elemPage[ePg++]=tetIST0&CLR31; break;	// if node is internal & unique
				case 2: edgeFlags|=TETREF0; tF0&=1;							elemPage[ePg++]=tetRef0; break;		// external & nonunique
				case 3:	edgeFlags|=TETREF0; edgeFlags^=FLG_EST0; tF0&=1;	elemPage[ePg++]=tetRef0; }			// internal & nonunique, saving element ref.
				switch (tF1) {
				case 0:	edgeFlags|=321;										elemPage[ePg++]=tetIST1&CLR31; break;
				case 1:	edgeFlags|=321; edgeFlags^=FLG_EST1;				elemPage[ePg++]=tetIST1&CLR31; break;
				case 2: edgeFlags|=TETREF1; tF1&=1;							elemPage[ePg++]=tetRef1; break;
				case 3:	edgeFlags|=TETREF1; edgeFlags^=FLG_EST1; tF1&=1;	elemPage[ePg++]=tetRef1; }
				switch (tF2) {
				case 0:	edgeFlags|=1092;									elemPage[ePg++]=tetIST2&CLR31; break;
				case 1:	edgeFlags|=1092; edgeFlags^=FLG_EST2;				elemPage[ePg++]=tetIST2&CLR31; break;
				case 2: edgeFlags|=TETREF2; tF2&=1;							elemPage[ePg++]=tetRef2; break;
				case 3:	edgeFlags|=TETREF2; edgeFlags^=FLG_EST2; tF2&=1;	elemPage[ePg++]=tetRef2; }
				switch (tF3) {
				case 0:	edgeFlags|=1296;									elemPage[ePg++]=tetIST3&CLR31; break;
				case 1:	edgeFlags|=1296; edgeFlags^=FLG_EST3;				elemPage[ePg++]=tetIST3&CLR31; break;
				case 2: edgeFlags|=TETREF3; tF3&=1;							elemPage[ePg++]=tetRef3; break;
				case 3:	edgeFlags|=TETREF3; edgeFlags^=FLG_EST3; tF3&=1;	elemPage[ePg++]=tetRef3; }

				if ((edgeFlags & FLG_EBND) == 0) {										// if element is internal, pattern unnecessary
					ePg++;
					elemPage[ePg++] = SET30|edgeFlags|slotcodeIST<<12|edgeCodeIST;			// bit 31 -> volume-internal non-boundary element
					elementsIntIST++;
				} else {
					// generate & store the pattern
					elemPage[ePg++] = tF0|tF1<<3|tF0<<5|tF2<<8|tF0<<10|tF3<<13|tF1<<15|tF2<<18|tF1<<20|tF3<<23|tF2<<25|tF3<<28;
					// generate & store the edge testing flags (two bits must be set for a legitimate testing edge)
					elemPage[ePg++] = mirrorFlagIST|edgeFlags|slotcodeIST<<12|edgeCodeIST;	// tag on compaction code & edges long/short code
				}
				slotcodeIST=slotcode; edgeCodeIST=edgeCode; mirrorFlagIST=mirrorFlag;		// store data for the next element writedown
				elementsIST++;																// we've made another entry in the boundary element array
			} else {
				// DEBUG: should not end up here since tetrahedron external-ness was already tested!
				throw new RuntimeException("FEM1.processISTnode(): Tetrahedron was external despite being checked as boundary.");
			}
			if ((codes&SET31) != 0) return;						// if last element writedown was flagged, do nothing more
		}
		
		// before it gets flagged locally, check if node (as specified by it's in-octant position) is locally unique for octant (or its neighbour centroid)
		boolean locallyUnique = aspectCode==Cd_CCCn && o.nodeI[60+seq2] == 0 || (o.enumerator & nodeVisitEnumF[aspectCode]) == 0;
		
		double x = 0, y = 0, z = 0;
		int intExt = 0;
		switch (aspectCode) {
		case Cd_MMM: x=o.xM; y=o.yM; z=o.zM; o.enumerator|=1;     intExt = o.int_ext; break;
		case Cd_PMM: x=o.xP; y=o.yM; z=o.zM; o.enumerator|=1<<2;  intExt = o.int_ext>>1; break;
		case Cd_MPM: x=o.xM; y=o.yP; z=o.zM; o.enumerator|=1<<6;  intExt = o.int_ext>>2; break;
		case Cd_PPM: x=o.xP; y=o.yP; z=o.zM; o.enumerator|=1<<8;  intExt = o.int_ext>>3; break;
		case Cd_MMP: x=o.xM; y=o.yM; z=o.zP; o.enumerator|=1<<18; intExt = o.int_ext>>4; break;
		case Cd_PMP: x=o.xP; y=o.yM; z=o.zP; o.enumerator|=1<<20; intExt = o.int_ext>>5; break;
		case Cd_MPP: x=o.xM; y=o.yP; z=o.zP; o.enumerator|=1<<24; intExt = o.int_ext>>6; break;
		case Cd_PPP: x=o.xP; y=o.yP; z=o.zP; o.enumerator|=1<<26; intExt = o.int_ext>>7; break;
		case Cd_CMM: x=o.xC; y=o.yM; z=o.zM; o.enumerator|=1<<1;  intExt = o.int_ext>>8; break;
		case Cd_CPM: x=o.xC; y=o.yP; z=o.zM; o.enumerator|=1<<7;  intExt = o.int_ext>>12; break;
		case Cd_CMP: x=o.xC; y=o.yM; z=o.zP; o.enumerator|=1<<19; intExt = o.int_ext>>22; break;
		case Cd_CPP: x=o.xC; y=o.yP; z=o.zP; o.enumerator|=1<<25; intExt = o.int_ext>>26; break;
		case Cd_MCM: x=o.xM; y=o.yC; z=o.zM; o.enumerator|=1<<3;  intExt = o.int_ext>>9; break;
		case Cd_PCM: x=o.xP; y=o.yC; z=o.zM; o.enumerator|=1<<5;  intExt = o.int_ext>>11; break;
		case Cd_MCP: x=o.xM; y=o.yC; z=o.zP; o.enumerator|=1<<21; intExt = o.int_ext>>23; break;
		case Cd_PCP: x=o.xP; y=o.yC; z=o.zP; o.enumerator|=1<<23; intExt = o.int_ext>>25; break;
		case Cd_MMC: x=o.xM; y=o.yM; z=o.zC; o.enumerator|=1<<9;  intExt = o.int_ext>>13; break;
		case Cd_PMC: x=o.xP; y=o.yM; z=o.zC; o.enumerator|=1<<11; intExt = o.int_ext>>15; break;
		case Cd_MPC: x=o.xM; y=o.yP; z=o.zC; o.enumerator|=1<<15; intExt = o.int_ext>>19; break;
		case Cd_PPC: x=o.xP; y=o.yP; z=o.zC; o.enumerator|=1<<17; intExt = o.int_ext>>21; break;
		case Cd_CCC: x=o.xC; y=o.yC; z=o.zC; o.enumerator|=1<<13; intExt = o.int_ext>>17; break;
		case Cd_CCCn:x=n.xC; y=n.yC; z=n.zC; 					  intExt = n.int_ext>>17; break;
		case Cd_MCC: x=o.xM; y=o.yC; z=o.zC; o.enumerator|=1<<12; intExt = o.int_ext>>16; break;
		case Cd_PCC: x=o.xP; y=o.yC; z=o.zC; o.enumerator|=1<<14; intExt = o.int_ext>>18; break;
		case Cd_CMC: x=o.xC; y=o.yM; z=o.zC; o.enumerator|=1<<10; intExt = o.int_ext>>14; break;
		case Cd_CPC: x=o.xC; y=o.yP; z=o.zC; o.enumerator|=1<<16; intExt = o.int_ext>>20; break;
		case Cd_CCM: x=o.xC; y=o.yC; z=o.zM; o.enumerator|=1<<4;  intExt = o.int_ext>>10; break;
		case Cd_CCP: x=o.xC; y=o.yC; z=o.zP; o.enumerator|=1<<22; intExt = o.int_ext>>24; break;
		default: throw new RuntimeException("FEM1.processISTnode(): Invalid aspect code received."); }
				
		int nodeN = 0, nodeRefN = elementsIST, local = aspectCode == Cd_CCCn ? 60 + seq2 : aspectCode<<1;
		// this block utilises the quick nodecode-based octant-local node index storage, a node exists locally if locallyUnique=false
		if (locallyUnique) {															// if node appears unique for octant, insert into hash table
			debugIST_locallyUnique++;
			switch (slotcode) {															// prepare node specification in reference to this node
			case Sl_Chg1: nodeRefN |= 1<<30; break; case Sl_Chg2: nodeRefN |= 2<<30; break; case Sl_Chg3: nodeRefN |= 3<<30; }
			
			long nodeLN = nodeHTableIST.add(x, y, z, nodeRefN, true);					// note: "indicate" flag set, duplicate node will be checked for & it's index returned
			nodeN=(int)(nodeLN&0xFFFFFFFFL)|(intExt=(intExt&1)<<30);
			nodeRefN=(int)(nodeLN>>32&NINT);											// if node wasn't unique globally, existing indexes will be returned					
			o.nodeI[local++] = nodeN | intExt; o.nodeI[local] = nodeRefN;
		} else {																		// direct-retrieve index locally
			nodeN = o.nodeI[local++] | SET31; nodeRefN = o.nodeI[local];
		}

		// the slotcode decides how to shuffle the four tetrahedral node indexes
		switch (slotcode) {
		case Sl_Mv01: tetIST1 = tetIST0; tetRef1 = tetRef0; tetIST0 = nodeN; tetRef0 = nodeRefN; break;
		case Sl_FtoL: tetIST3 = tetIST0; tetRef3 = tetRef0;
		case Sl_ChgA: case Sl_Hld2: case Sl_Chg0: case Sl_Ch01: case Sl_Ch03: tetIST0 = nodeN; tetRef0 = nodeRefN; break;
		case Sl_Chg1: tetIST1 = nodeN; tetRef1 = nodeRefN; break;
		case Sl_Chg2: tetIST2 = nodeN; tetRef2 = nodeRefN; break;
		case Sl_Chg3: tetIST3 = nodeN; tetRef3 = nodeRefN; break;
		}
		if (DEBUG_LEVEL > 1) debugIST_pISTn_time += System.nanoTime() - timerS;
	}
	
	
	// method unpacks an element node stream into coordinate quadruplets signifying tetrahedral elements
	public double[] elementStreamToArray(double[] strN, long[] strC, int elements, int[] arrN, boolean unpackIndexes) {
		long slotcodes = 0;
		double x0=0, y0=0, z0=0, x1=0, y1=0, z1=0, x2=0, y2=0, z2=0, x3=0, y3=0, z3=0;
		double[] nA = new double[elements * 4 * 3];

		int iStrm = 1, i4 = 0, idx0=-1, idx1=-1, idx2=-1, idx3=-1;
		for (int s3 = 0, sC = 0, s16 = 0, n3 = 0, e = 0; e < elements; e++) {
			
			if (--s16 <= 0) { slotcodes = strC[sC++]; s16 = 16; } else slotcodes >>= 4;
			int slotcode = (int)(slotcodes & 0xF);
			
			switch (slotcode) {												// find out what nodes to get to complete an element
			case Sl_FtoL:	x3 = x0; y3 = y0; z3 = z0;
							x0 = strN[s3++]; y0 = strN[s3++]; z0 = strN[s3++]; break;
			case Sl_Mv01:	x1 = x0; y1 = y0; z1 = z0;
							x0 = strN[s3++]; y0 = strN[s3++]; z0 = strN[s3++]; break;
			case Sl_ShMB:	x3 = x2; y3 = y2; z3 = z2; x2 = x1; y2 = y1; z2 = z1;
			case Sl_Ch01:	x0 = strN[s3++]; y0 = strN[s3++]; z0 = strN[s3++];
			case Sl_Chg1: 	x1 = strN[s3++]; y1 = strN[s3++]; z1 = strN[s3++]; break;
			case Sl_Ch02:	x0 = strN[s3++]; y0 = strN[s3++]; z0 = strN[s3++];
			case Sl_Chg2:	x2 = strN[s3++]; y2 = strN[s3++]; z2 = strN[s3++]; break;
			case Sl_ChgA: case Sl_Ch03:
			case Sl_Hld2:	x0 = strN[s3++]; y0 = strN[s3++]; z0 = strN[s3++];
						if (slotcode==Sl_Hld2 || slotcode==Sl_ChgA) {
							x1 = strN[s3++]; y1 = strN[s3++]; z1 = strN[s3++]; }
						if (slotcode==Sl_ChgA) {
							x2 = strN[s3++]; y2 = strN[s3++]; z2 = strN[s3++]; }
			case Sl_Chg3:	x3 = strN[s3++]; y3 = strN[s3++]; z3 = strN[s3++]; break;
			}
			if (unpackIndexes) {
				switch (slotcode) {
				case Sl_FtoL: idx3 = idx0; idx0 = iStrm++; break; case Sl_Mv01: idx1 = idx0; idx0 = iStrm++; break;
				case Sl_ShMB: idx3 = idx2; idx2 = idx1; case Sl_Ch01: idx0 = iStrm++; idx1 = iStrm++; break;
				case Sl_Chg1: idx1 = iStrm++; break; case Sl_Chg2: idx2 = iStrm++; break;
				case Sl_Ch02: idx0 = iStrm++; idx2 = iStrm++; break; case Sl_Ch03: idx0 = iStrm++; case Sl_Chg3: idx3 = iStrm++; break;
				case Sl_Hld2: idx0 = iStrm++; idx1 = iStrm++; idx3 = iStrm++; break;
				case Sl_ChgA: idx0 = iStrm++; idx1 = iStrm++; idx2 = iStrm++; idx3 = iStrm++; }
				arrN[i4++] = idx0; arrN[i4++] = idx1; arrN[i4++] = idx2; arrN[i4++] = idx3;
			}		
			nA[n3++] = x0; nA[n3++] = y0; nA[n3++] = z0; nA[n3++] = x1; nA[n3++] = y1; nA[n3++] = z1; 	// write down element nodes
			nA[n3++] = x2; nA[n3++] = y2; nA[n3++] = z2; nA[n3++] = x3; nA[n3++] = y3; nA[n3++] = z3; 
		}
		return nA;
	}
	
	// internal element paging array space checker
	private int[] verifySpaceElementPg() {
		if (elementPgSize <= ePg) {											// if we've run out of space in current element array page
			if (elementPgI + 1 >= elementPgN) {								// if we've run out of element array pages, add 400 more
				int[][] elementPgNew = new int[elementPgN += 400][];
				for (int i = 0; i <= elementPgI; i++) elementPgNew[i] = elementPg[i];
				elementPg = elementPgNew;
			}
			ePg = 0;
			return elementPg[++elementPgI] = new int[elementPgSize]; 			// allocate new page
		}
		return elementPg[elementPgI];
	}
	
	
	// method takes any one of the element arrays generated by FEM1 class, a "skip" provided will jump past any status data,
	// the node array must be provided in nArray[] and the multiplier for node offsets must be given in "multN" (usually = 3)
	static int findElementIndex(int[] eArray, double[] nArray, double x, double y, double z, int skip, int multN) {
		double shortestDist = Double.MAX_VALUE;
		int closestE = -1;
		for (int e = 0, eA = 0; eA < eArray.length; e++) {
			int n0 = (eArray[eA++]&0x3FFFFFFF)*multN, n1 = (eArray[eA++]&0x3FFFFFFF)*multN;
			int n2 = (eArray[eA++]&0x3FFFFFFF)*multN, n3 = (eArray[eA++]&0x3FFFFFFF)*multN;
			eA += skip;
			double x0=nArray[n0++], y0=nArray[n0++], z0=nArray[n0++], x1=nArray[n1++], y1=nArray[n1++], z1=nArray[n1++];
			double x2=nArray[n2++], y2=nArray[n2++], z2=nArray[n2++], x3=nArray[n3++], y3=nArray[n3++], z3=nArray[n3++];
			double xc = (x0+x1+x2+x3)*0.25, yc = (y0+y1+y2+y3)*0.25, zc = (z0+z1+z2+z3)*0.25;
			double dx = xc - x, dy = yc - y, dz = zc - z, dist = Math.sqrt(dx*dx+dy*dy+dz*dz);
			if (dist < shortestDist) { shortestDist = dist; closestE = e; }
		}
		return closestE;
	}
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			STENCIL PERMUTATIONS GENERATION CODE FOR ISOSURFACE STUFFING
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	

	// these are the stencil seeds
	private final static int[] stencilSeed = {	1,2,0,2,		// the +0-0 stencil 1 (all aspects)
												1,0,0,2,		// the +--0 stencil 2 (all aspects)
												1,0,0,0,		// the +--- stencil 3 (all aspects)
												0,2,1,1,		// the -0++ stencil 4 (long/short edge aligned)
												0,1,2,1,		// the -0++ stencil 4 mirrored -> -+0+
												1,0,1,0,		// the +-+- stencil 5 (long/short edge aligned)
												1,1,0,0,		// the +-+- stencil 5 mirrored -> ++--
												0,1,1,2,		// the -++0 stencil 6 (long/short edge aligned)
												0,1,1,1,		// the -+++ stencil 7 (long/short edge aligned, parity rule applies)
												0,1,1,0};		// the -++- stencil 8 (long/short edge aligned, parity rule applies)
	private final static int[] stencilEdgeSeed = {
			0,2,-1,-1,-1,-1,-1,-1,	0,1,0,2,-1,-1,-1,-1,	0,1,0,2,0,3,-1,-1,
			0,2,0,3,-1,-1,-1,-1,	0,1,0,3,-1,-1,-1,-1,
			0,1,0,3,1,2,2,3,		0,2,0,3,1,2,1,3,
			0,1,0,2,-1,-1,-1,-1,	0,1,0,2,0,3,-1,-1,		0,1,0,2,1,3,2,3 };
	private final static int[] stencilPerm = {
			0,1,2,3, 0,2,3,1, 0,3,1,2, 1,3,2,0, 1,2,0,3, 1,0,3,2, 2,1,3,0, 2,0,1,3, 2,3,0,1, 3,2,1,0, 3,1,0,2, 3,0,2,1,	// the all-aspect permutations
			0,1,2,3, 3,2,1,0, 2,3,0,1, 1,0,3,2 };																		// long/short edge mirror permutations

	public static int[] stencilPermutations() {
		int[] stencil = new int[(12*3 + 4*4 + 4*3)*2];			// 3 stencils permutate along any positive node combo, 2 mirror&reflect, 3 only reflecting
		for (int pt = 0, s = 0, es = 0; pt < 10*4; pt += 4, es += 8) {
			int pmEnd = (pt >= 3*4 ? 16*4 : 12*4);				// choose the 3 first or 7 last permutation combos
			for (int pm = pt >= 3*4 ? 48 : 0; pm < pmEnd;) {
				int[] eSeed = stencilEdgeSeed.clone();
				for (int e = es; e < es+8;) {
					if (eSeed[e] < 0) { e += 2; continue; }
					for (int c = 0; c < 4; c++) if (stencilPerm[pm + c] == eSeed[e]) { eSeed[e] = c; break; } e++;
					for (int c = 0; c < 4; c++) if (stencilPerm[pm + c] == eSeed[e]) { eSeed[e] = c; break; }
					if (eSeed[e] < eSeed[e-1]) { int tmp = eSeed[e]; eSeed[e] = eSeed[e-1]; eSeed[e-1] = tmp; } e++;
				}
				int pm0 = stencilPerm[pm++], pm1 = stencilPerm[pm++], pm2 = stencilPerm[pm++], pm3 = stencilPerm[pm++];
				int n0 = stencilSeed[pt+pm0], n1 = stencilSeed[pt+pm1], n2 = stencilSeed[pt+pm2], n3 = stencilSeed[pt+pm3];
				int pattern = n0|n1<<3|n0<<5|n2<<8|n0<<10|n3<<13|n1<<15|n2<<18|n1<<20|n3<<23|n2<<25|n3<<28;
				// indexMix bitpacks the intended order of node picking for an element permutation (4 bits per element-local node index)
				// and intended order of picking the edge cutnodes (4 bits per element-local edge index)
				int indexMix = pm0 | pm1<<4 | pm2<<8 | pm3<<12;
				for (int e = es, bOfs=16; e < es + 8; e += 2, bOfs+=4) {
					if (eSeed[e]<0) { indexMix |= 8 << bOfs; continue; }
					int eSeedA = eSeed[e], eSeedB = eSeed[e+1];
					if (eSeedA==0 && eSeedB==1) { pattern|=EC_C01; } else if (eSeedA==0 && eSeedB==2) { pattern|=EC_C02; } else
					if (eSeedA==0 && eSeedB==3) { pattern|=EC_C03; } else if (eSeedA==1 && eSeedB==2) { pattern|=EC_C12; } else
					if (eSeedA==1 && eSeedB==3) { pattern|=EC_C13; } else if (eSeedA==2 && eSeedB==3) { pattern|=EC_C23; }
					indexMix |= (eSeedA==0||eSeedB==0?eSeedA+eSeedB-1:eSeedA+eSeedB) << bOfs;
				}
				stencil[s++] = pattern; stencil[s++] = indexMix;
			}
		}
		return stencil;
	}
	
	boolean stencilPermCodeSortedIST = false;		// flags whether stencil permutations codes has been sorted (so one can use binary search)
	void sortStencilPermCode() {
		if (stencilPermCodeSortedIST) return; boolean done = false;
		while (!done) {
			done = true;
			for (int i = 0, i2 = 0; i < stencilPermCode.length-2; i+=2, i2++) {
				if (stencilPermCode[i] > stencilPermCode[i+2]) {
					int tmp = stencilPermCode[i]; stencilPermCode[i] = stencilPermCode[i+2]; stencilPermCode[i+2] = tmp;
					tmp = stencilPermCode[i+1]; stencilPermCode[i+1] = stencilPermCode[i+3]; stencilPermCode[i+3] = tmp;
					byte tmpB = stencilPermSwitch[i2]; stencilPermSwitch[i2] = stencilPermSwitch[i2+1]; stencilPermSwitch[i2+1] = tmpB; done = false; }
			}
		} 
		stencilPermCodeSortedIST = true;
	}
	
	public static String stencilPermToJavaCode(int[] stencil) {
		System.out.println();
		String stencilCode = "final static int[] stencilPermCode = {\n\t// the 3 stencils of arbitrary orientation & reflection\n";
		for (int c = 0; c < stencil.length;) {
			if (c==12*3*2)		 stencilCode += "\t// 2 stencils limited to LONG/SHORT edge orientation & mirror permutation along assymetric axis\n";
			if (c==12*3*2+4*4*2) stencilCode += "\t// 3 stencils limited to LONG/SHORT edge orientation & to Parity Rule (matching COMPOUND tet. edge orientation neighbour elements)\n";
			stencilCode += "\t" + FEM1.stencilToString(stencil[c++], true)+", "+FEM1.indexMixToString(stencil[c++]) + (c<stencil.length-1?",\n":" };\n");
		}
		return stencilCode;
	}
	
	public static String stencilToString(int s, boolean outputCode) {
		String s0 = (s&3)==0?"-":((s&3)==1?"+":"0"), s1 = (s>>3&3)==0?"-":((s>>3&3)==1?"+":"0");
		String s2 = (s>>8&3)==0?"-":((s>>8&3)==1?"+":"0"), s3 = (s>>13&3)==0?"-":((s>>13&3)==1?"+":"0");
		String c01 = (s&4)==4?"o":" ", c02 = (s&128)==128?"o":" ", c03 = (s&4096)==4096?"o":" ", c12 = (s&131072)==131072?"o":" ";
		String c13 = (s&4194304)==4194304?"o":" ", c23 = (s&134217728)==134217728?"o":" ";
		String test = "["+s0+c01+s1+"]["+s0+c02+s2+"]["+s0+c03+s3+"]["+s1+c12+s2+"]["+s1+c13+s3+"]["+s2+c23+s3+"]\n";
		String code = "";
		if (outputCode) for (int i = 0, shft = 0; i < 6; i++, s>>=5, shft += 5)
				code+="EC_"+((s&3)==0?"M":((s&3)==1?"P":"0"))+((s&4)==4?"x":"_")+((s&24)==0?"M":((s&24)==8?"P":"0"))+(shft>0?"<<"+shft:"")+(i<5?"|":"");		
		return outputCode ? code : test;
	}
	
	public static String indexMixToString(int im) {
		return (im&15)+"|"+(im>>4&15)+"<<4|"+(im>>8&15)+"<<8|"+(im>>12&15)+"<<12 | "+(im>>16&15)+"<<16|"+(im>>20&15)+"<<20|"+(im>>24&15)+"<<24|"+(im>>28&15)+"<<28";
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
					NSPArray aHsp = M.Hsp[i];
					if (aHsp.array == null) {
						aHsp.array = new NSPNode[rCounts[i]];
						aHsp.nodes = aHsp.size = rCounts[i];
					}
					if (cOffP >= 0 && aHsp.array[cOffP].c() == j) {
						aHsp.array[cOffP].add(sA[e]);
					} else {
						NSPNode node = aHsp.array[cOffsets[i]] = new NSPNode(i, j, sA[e], 0, cOffsets[i]++, 0);
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
	
	public boolean externalNode(int n) { return n * NCOORD > node.length ? (nodeWorkFlag[n - nodes] & 1) == 1 : (nodeFlag[n] & 1) == 1; }
	
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

	
	static int[] findNode(double[] node, double x, double y, double z, double error) {
		int[] finds = new int[16]; int f=0;
		for (int i3 = 0; i3 < node.length;) {
			double dx = node[i3++] - x, dy = node[i3++] - y, dz = node[i3++] - z;
			if (dx<error && dx>-error && dy<error && dy>-error && dz<error && dz>-error) { finds[f++] = (i3-3)/3; if (f>=16) return finds; }}
		return f==0 ? null : finds;
	}
	
	
	public static double tetraQuality(double[] node, int n0, int n1, int n2, int n3) {
		n0 *= NCOORD; n1 *= NCOORD; n2 *= NCOORD; n3 *= NCOORD;
		double x0=node[n0++], y0=node[n0++], z0=node[n0], x1=node[n1++], y1=node[n1++], z1=node[n1];
		double x2=node[n2++], y2=node[n2++], z2=node[n2], x3=node[n3++], y3=node[n3++], z3=node[n3];
		double x10 = x1-x0, y10 = y1-y0, z10 = z1-z0, e01 = Math.sqrt(x10*x10+y10*y10+z10*z10);
		double x20 = x2 - x0, y20 = y2 - y0, z20 = z2 - z0, e02 = Math.sqrt(x20*x20+y20*y20+z20*z20);
		double x30 = x3 - x0, y30 = y3 - y0, z30 = z3 - z0, e03 = Math.sqrt(x30*x30+y30*y30+z30*z30);
		double x21 = x2 - x1, y21 = y2 - y1, z21 = z2 - z1, e12 = Math.sqrt(x21*x21+y21*y21+z21*z21);
		double x31 = x3 - x1, y31 = y3 - y1, z31 = z3 - z1, e13 = Math.sqrt(x31*x31+y31*y31+z31*z31);
		double x32 = x3 - x2, y32 = y3 - y2, z32 = z3 - z2, e23 = Math.sqrt(x32*x32+y32*y32+z32*z32);
		double sq1 = y10 * z20 - z10 * y20, sq2 = z10 * x20 - x10 * z20, sq3 = x10 * y20 - y10 * x20;
		double area012 = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);						// facet 012 area	
		double Vt6 = (sq1 * x30 + sq2 * y30 + sq3 * z30);									// Vt6 is tetrahedron volume * 6
		if (Vt6 < 0) Vt6 = -Vt6;
		sq1 = y10 * z30 - z10 * y30; sq2 = z10 * x30 - x10 * z30; sq3 = x10 * y30 - y10 * x30;
		double area023 = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);						// facet 023 area
		sq1 = y20 * z30 - z20 * y30; sq2 = z20 * x30 - x20 * z30; sq3 = x20 * y30 - y20 * x30;
		double area031 = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);						// facet 031 area
		sq1 = y21 * z31 - z21 * y31; sq2 = z21 * x31 - x21 * z31; sq3 = x21 * y31 - y21 * x31;
		double area132 = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);						// facet 132 area
		double e01e23 = e01 * e23, e02e13 = e02 * e13, e03e12 = e03 * e12;
		double M5 = (e01e23 + e02e13 + e03e12) * (-e01e23 + e02e13 + e03e12) * (e01e23 - e02e13 + e03e12) * (e01e23 + e02e13 - e03e12) * 0.0625;
		//double cRadius = Math.sqrt(M5 < 0 ? -M5 : M5) / Vt6, iRadius = 0.5 * Vt6 / (area012 + area023 + area031 + area132);
		return (3./2.) * Vt6*Vt6 / (Math.sqrt(M5 < 0 ? -M5 : M5) * (area012 + area023 + area031 + area132));
	}

	public double tetraQualityIST(int n0, int n1, int n2, int n3, byte snap0, byte snap1, byte snap2, byte snap3) {
		double x0=nodeIST[n0=snap0<0?n0*6+3:n0*6], y0=nodeIST[++n0], z0=nodeIST[++n0];
		double x1=nodeIST[n1=snap1<0?n1*6+3:n1*6], y1=nodeIST[++n1], z1=nodeIST[++n1];
		double x2=nodeIST[n2=snap2<0?n2*6+3:n2*6], y2=nodeIST[++n2], z2=nodeIST[++n2];
		double x3=nodeIST[n3=snap3<0?n3*6+3:n3*6], y3=nodeIST[++n3], z3=nodeIST[++n3];
		double x10 = x1-x0, y10 = y1-y0, z10 = z1-z0, e01 = Math.sqrt(x10*x10+y10*y10+z10*z10);
		double x20 = x2 - x0, y20 = y2 - y0, z20 = z2 - z0, e02 = Math.sqrt(x20*x20+y20*y20+z20*z20);
		double x30 = x3 - x0, y30 = y3 - y0, z30 = z3 - z0, e03 = Math.sqrt(x30*x30+y30*y30+z30*z30);
		double x21 = x2 - x1, y21 = y2 - y1, z21 = z2 - z1, e12 = Math.sqrt(x21*x21+y21*y21+z21*z21);
		double x31 = x3 - x1, y31 = y3 - y1, z31 = z3 - z1, e13 = Math.sqrt(x31*x31+y31*y31+z31*z31);
		double x32 = x3 - x2, y32 = y3 - y2, z32 = z3 - z2, e23 = Math.sqrt(x32*x32+y32*y32+z32*z32);
		double sq1 = y10 * z20 - z10 * y20, sq2 = z10 * x20 - x10 * z20, sq3 = x10 * y20 - y10 * x20;
		double area012 = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);						// facet 012 area	
		double Vt6 = (sq1 * x30 + sq2 * y30 + sq3 * z30);									// Vt6 is tetrahedron volume * 6
		if (Vt6 < 0) Vt6 = -Vt6;
		sq1 = y10 * z30 - z10 * y30; sq2 = z10 * x30 - x10 * z30; sq3 = x10 * y30 - y10 * x30;
		double area023 = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);						// facet 023 area
		sq1 = y20 * z30 - z20 * y30; sq2 = z20 * x30 - x20 * z30; sq3 = x20 * y30 - y20 * x30;
		double area031 = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);						// facet 031 area
		sq1 = y21 * z31 - z21 * y31; sq2 = z21 * x31 - x21 * z31; sq3 = x21 * y31 - y21 * x31;
		double area132 = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);						// facet 132 area
		double e01e23 = e01 * e23, e02e13 = e02 * e13, e03e12 = e03 * e12;
		double M5 = (e01e23 + e02e13 + e03e12) * (-e01e23 + e02e13 + e03e12) * (e01e23 - e02e13 + e03e12) * (e01e23 + e02e13 - e03e12) * 0.0625;
		//double cRadius = Math.sqrt(M5 < 0 ? -M5 : M5) / Vt6, iRadius = 0.5 * Vt6 / (area012 + area023 + area031 + area132);
		return (3./2.) * Vt6*Vt6 / (Math.sqrt(M5 < 0 ? -M5 : M5) * (area012 + area023 + area031 + area132));
	}

	public static double tetraVolumePositivity(double[] node, int n0, int n1, int n2, int n3) {
		n0 *= NCOORD; n1 *= NCOORD; n2 *= NCOORD; n3 *= NCOORD;
		double x0=node[n0++], y0=node[n0++], z0=node[n0], x1=node[n1++], y1=node[n1++], z1=node[n1], x2=node[n2++], y2=node[n2++], z2=node[n2];
		double x10 = x1 - x0, y10 = y1 - y0, z10 = z1 - z0, x20 = x2 - x0, y20 = y2 - y0, z20 = z2 - z0;
		double n012x = y10 * z20 - z10 * y20, n012y = z10 * x20 - x10 * z20, n012z = x10 * y20 - y10 * x20;
		double nx = node[n3++] - x0, ny = node[n3++] - y0, nz = node[n3] - z0;
		return -(n012x * nx + n012y * ny + n012z * nz);
	}

	public double tetraVolumePositivityIST(int n0, int n1, int n2, int n3, byte snap0, byte snap1, byte snap2, byte snap3) {
		double x0 = nodeIST[n0=snap0<0?n0*6+3:n0*6], y0 = nodeIST[++n0], z0 = nodeIST[++n0];
		double x1 = nodeIST[n1=snap1<0?n1*6+3:n1*6], y1 = nodeIST[++n1], z1 = nodeIST[++n1];
		double x2 = nodeIST[n2=snap2<0?n2*6+3:n2*6], y2 = nodeIST[++n2], z2 = nodeIST[++n2];
		double x10 = x1 - x0, y10 = y1 - y0, z10 = z1 - z0, x20 = x2 - x0, y20 = y2 - y0, z20 = z2 - z0;
		double n012x = y10 * z20 - z10 * y20, n012y = z10 * x20 - x10 * z20, n012z = x10 * y20 - y10 * x20;
		double nx = nodeIST[n3=snap3<0?n3*6+3:n3*6] - x0, ny = nodeIST[++n3] - y0, nz = nodeIST[++n3] - z0;
		return -(n012x * nx + n012y * ny + n012z * nz);
	}
	
	public boolean sliverCentroidInternal(FEM1Octree geocTree, int n0, int n1, int n2, int n3) {
		double x0=nodeIST[n0=n0*6+3], y0=nodeIST[++n0], z0=nodeIST[++n0], x1=nodeIST[n1=n1*6+3], y1=nodeIST[++n1], z1=nodeIST[++n1];
		double x2=nodeIST[n2=n2*6+3], y2=nodeIST[++n2], z2=nodeIST[++n2], x3=nodeIST[n3=n3*6+3], y3=nodeIST[++n3], z3=nodeIST[++n3];
		double xc = (x0+x1+x2+x3)*0.25, yc = (y0+y1+y2+y3)*0.25, zc = (z0+z1+z2+z3)*0.25;
		int[] status = {0, 0};
		geocTree.root.facetEncounters(geocTree, xc, yc, zc, xc, yc, geocTree.root.zP+0.001, 3, status);
		return status[0] != 0;
	}
	
	final static double bias_eCI3S = .31, intl_eCI3S = 1.0-(3*bias_eCI3S);		// bias_eCI3S should not exceed .3333
	// method checks if a tetrahedral centroid is internal or external, method expects 3 snapped nodes and one can bias versus or away from them
	public boolean elementCentroidInternal3S(FEM1Octree geocTree, int n0, int n1, int n2, int n3, byte snap0, byte snap1, byte snap2, byte snap3) {
		double f0 = bias_eCI3S, f1 = bias_eCI3S, f2 = bias_eCI3S, f3 = bias_eCI3S;
		if (snap0>=0) f0=intl_eCI3S; else if (snap1>=0) f1=intl_eCI3S; else if (snap2>=0) f2=intl_eCI3S; else f3=intl_eCI3S;
		double x0=nodeIST[n0=snap0<0?n0*6+3:n0*6], y0=nodeIST[++n0], z0=nodeIST[++n0];
		double x1=nodeIST[n1=snap1<0?n1*6+3:n1*6], y1=nodeIST[++n1], z1=nodeIST[++n1];
		double x2=nodeIST[n2=snap2<0?n2*6+3:n2*6], y2=nodeIST[++n2], z2=nodeIST[++n2];
		double x3=nodeIST[n3=snap3<0?n3*6+3:n3*6], y3=nodeIST[++n3], z3=nodeIST[++n3];
		double xc = x0*f0+x1*f1+x2*f2+x3*f3, yc = y0*f0+y1*f1+y2*f2+y3*f3, zc = z0*f0+z1*f1+z2*f2+z3*f3; int[] status = {0, 0};
		geocTree.root.facetEncounters(geocTree, xc, yc, zc, xc, yc, geocTree.root.zP+0.001, 3, status);
		return status[0] != 0;
	}

	final static double bias_eCI2S = .40, intl_eCI2S = 1.0-(2*bias_eCI2S);		// bias_eCI2S should not exceed .4999
	// method expects two nodes to be snapped, one can bias versus or away from them
	public boolean elementCentroidInternal2S(FEM1Octree geocTree, int n0, int n1, int n2, int n3, byte snap0, byte snap1, byte snap2, byte snap3) {
		double f0 = bias_eCI2S, f1 = bias_eCI2S, f2 = bias_eCI2S, f3 = bias_eCI2S;
		if (snap0>=0) f0=intl_eCI2S*.5; if (snap1>=0) f1=intl_eCI2S*.5; if (snap2>=0) f2=intl_eCI2S*.5; if (snap3>=0) f3=intl_eCI2S*.5;
		double x0=nodeIST[n0=snap0<0?n0*6+3:n0*6], y0=nodeIST[++n0], z0=nodeIST[++n0];
		double x1=nodeIST[n1=snap1<0?n1*6+3:n1*6], y1=nodeIST[++n1], z1=nodeIST[++n1];
		double x2=nodeIST[n2=snap2<0?n2*6+3:n2*6], y2=nodeIST[++n2], z2=nodeIST[++n2];
		double x3=nodeIST[n3=snap3<0?n3*6+3:n3*6], y3=nodeIST[++n3], z3=nodeIST[++n3];
		double xc = x0*f0+x1*f1+x2*f2+x3*f3, yc = y0*f0+y1*f1+y2*f2+y3*f3, zc = z0*f0+z1*f1+z2*f2+z3*f3; int[] status = {0, 0};
		geocTree.root.facetEncounters(geocTree, xc, yc, zc, xc, yc, geocTree.root.zP+0.001, 3, status);
		return status[0] != 0;
	}


	public static double[] tetraDihedralAngles(double[] node, int n0, int n1, int n2, int n3) {
		n0 *= NCOORD; n1 *= NCOORD; n2 *= NCOORD; n3 *= NCOORD;
		double[] angles = {0,0,0,0,0,0};
		double x0=node[n0++], y0=node[n0++], z0=node[n0], x1=node[n1++], y1=node[n1++], z1=node[n1];
		double x2=node[n2++], y2=node[n2++], z2=node[n2], x3=node[n3++], y3=node[n3++], z3=node[n3];
		double x10 = x1 - x0, y10 = y1 - y0, z10 = z1 - z0, x20 = x2 - x0, y20 = y2 - y0, z20 = z2 - z0;
		double n012x = y10 * z20 - z10 * y20, n012y = z10 * x20 - x10 * z20, n012z = x10 * y20 - y10 * x20;
		double n012D = 1.0/Math.sqrt(n012x*n012x+n012y*n012y+n012z*n012z); n012x *= n012D; n012y *= n012D; n012z *= n012D;
		double x30 = x3 - x0, y30 = y3 - y0, z30 = z3 - z0;
		double n023x = y20 * z30 - z20 * y30, n023y = z20 * x30 - x20 * z30, n023z = x20 * y30 - y20 * x30;
		double n023D = 1.0/Math.sqrt(n023x*n023x+n023y*n023y+n023z*n023z); n023x *= n023D; n023y *= n023D; n023z *= n023D;
		double n031x = y30 * z10 - z30 * y10, n031y = z30 * x10 - x30 * z10, n031z = x30 * y10 - y30 * x10;
		double n031D = 1.0/Math.sqrt(n031x*n031x+n031y*n031y+n031z*n031z); n031x *= n031D; n031y *= n031D; n031z *= n031D;
		double x21 = x2 - x1, y21 = y2 - y1, z21 = z2 - z1, x31 = x3 - x1, y31 = y3 - y1, z31 = z3 - z1;
		double n132x = y31 * z21 - z31 * y21, n132y = z31 * x21 - x31 * z21, n132z = x31 * y21 - y31 * x21;
		double n132D = 1.0/Math.sqrt(n132x*n132x+n132y*n132y+n132z*n132z); n132x *= n132D; n132y *= n132D; n132z *= n132D;
		double toDeg = 180.0 / Math.PI;
		angles[0] = Math.acos(-(n012x*n023x + n012y*n023y + n012z*n023z)) * toDeg;	// theta012v023
		angles[1] = Math.acos(-(n012x*n031x + n012y*n031y + n012z*n031z)) * toDeg;	// theta012v031
		angles[2] = Math.acos(-(n012x*n132x + n012y*n132y + n012z*n132z)) * toDeg;	// theta012v132
		angles[3] = Math.acos(-(n023x*n031x + n023y*n031y + n023z*n031z)) * toDeg;	// theta023v031
		angles[4] = Math.acos(-(n023x*n132x + n023y*n132y + n023z*n132z)) * toDeg;	// theta023v132
		angles[5] = Math.acos(-(n031x*n132x + n031y*n132y + n031z*n132z)) * toDeg;	// theta031v132
		return angles;
	}
	
	public static double tetraSmallestDihedral(double[] angles) {
		double min01 = angles[0] < angles[1] ? angles[0] : angles[1], min23 = angles[2] < angles[3] ? angles[2] : angles[3];
		double min45 = angles[4] < angles[5] ? angles[4] : angles[5], min01_a2 = min01 < min23 ? min01 : min23;
		return min45 < min01_a2 ? min45 : min01_a2; }

	public static double tetraLargestDihedral(double[] angles) {
		double max01 = angles[0] > angles[1] ? angles[0] : angles[1], max23 = angles[2] > angles[3] ? angles[2] : angles[3];
		double max45 = angles[4] > angles[5] ? angles[4] : angles[5], max01_a2 = max01 > max23 ? max01 : max23;
		return max45 > max01_a2 ? max45 : max01_a2; }

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	final private static String[] slotcodeName = {"ChgA","FtoL","Hld2","Ch01","Mv01","Ch03","ShMB","**7*","Chg0","Chg1","Chg2","Chg3","Ch02","*13*","*14*","*15*"};
	final private static String[] edgeCodeName = {"SSSSSS","","LSSSSS/SLSSSS","","","","","","SSSLSS","","","","SSLLSS","","","","SSSSLS",
													"","","","","","","","","","","","","","","","SSSSSL"};
	// note: START_PRINT_ELEMENT can be assigned by any algorithm that processes IST elements to tell toString() what element to start printing from
	private static int maxPrintElements = 100, startPrintElement = 0;
	private boolean readElementPg = true, finalArrayPg = false, printSlotcodes1 = false, printSlotcodes2 = false, printSnapNodes = false;

	public static void setDebugLevel(int l) { DEBUG_LEVEL = l > 0 ? l : 0; }
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		if (printSlotcodes1 ||printSlotcodes2) { String slotcodes = toStringSlotcodes(); return slotcodes; }
		
		// print Isosurface Stuffing data first, as it is first to be generated
		if (elementIST != null || elementPg != null) {
			sb.append("<<<<<<<< Isosurface Stuffing Tetrahedral Data >>>>>>>>\n");
			int printElems = startPrintElement+maxPrintElements < elementsIST ? startPrintElement+maxPrintElements : elementsIST;
			double[] nodeA = nodeIST;
			if (nodeA==null) nodeA = nodeHTableIST.array();
			for (int e = 0, e6 = 0, e6b = 0, eP = 0; e < printElems; e++) {
				int n0=0, n1=0, n2=0, n3=0, pattern=0, status=0;
				int [] elA = null;
				
				if (readElementPg) {
					e6b = e6 - eP * elementPgSize;
					if (e6b >= elementPgSize) { eP++; e6b -= elementPgSize; }
					elA =  elementPg[eP];
					n0 = elA[e6b++]; n1 = elA[e6b++]; n2 = elA[e6b++]; n3 = elA[e6b++];
					if (!finalArrayPg) { pattern = elA[e6b++]; status = elA[e6b++]; e6 += 6; } else { status = elA[e6b++]; e6 += 5; }
				} else {
					e6b = e6;
					n0 = elementIST[e6b++]; n1 = elementIST[e6b++]; n2 = elementIST[e6b++]; n3 = elementIST[e6b++];
					pattern = elementIST[e6b++]; status = elementIST[e6b++];	
					e6 += 6;
				}
				
				if (e < startPrintElement) continue;						// do not start printing until getting to the target element index
				
				boolean elemInternal = (status&FLG_EBND)==0||(status&SET30)==SET30, elemSliver=(status&SET29)==SET29, elemDiscard=(status&SET28)==SET28;
				String elemType = (elemDiscard?"Discarded ":(elemInternal?"Internal ":"Boundary ")) + (elemSliver?"sliver ":"")+"element ";
				int eR0=0, eR1=0, eR2=0, eR3=0, nX0=0, nX1=0, nX2=0, nX3=0;
				boolean ref0=false, ref1=false, ref2=false, ref3=false;
				if (!finalArrayPg) {
					if ((status&TETREF0)!=0) {
						eR0 = n0&0x3FFFFFFF; n0 = n0>>30&3; ref0 = true;
						nX0=readElementPg ? elementPg[(eR0*6)/elementPgSize][eR0*6-elementPgSize*((eR0*6)/elementPgSize)+n0] : elementIST[eR0*6+n0]; }
					if ((status&TETREF1)!=0) {
						eR1 = n1&0x3FFFFFFF; n1 = n1>>30&3; ref1 = true;
						nX1=readElementPg ? elementPg[(eR1*6)/elementPgSize][eR1*6-elementPgSize*((eR1*6)/elementPgSize)+n1] : elementIST[eR1*6+n1]; }
					if ((status&TETREF2)!=0) {
						eR2 = n2&0x3FFFFFFF; n2 = n2>>30&3; ref2 = true;
						nX2=readElementPg ? elementPg[(eR2*6)/elementPgSize][eR2*6-elementPgSize*((eR2*6)/elementPgSize)+n2] : elementIST[eR2*6+n2]; }
					if ((status&TETREF3)!=0) {
						eR3 = n3&0x3FFFFFFF; n3 = n3>>30&3; ref3 = true;
						nX3=readElementPg ? elementPg[(eR3*6)/elementPgSize][eR3*6-elementPgSize*((eR3*6)/elementPgSize)+n3] : elementIST[eR3*6+n3]; }
				}
				sb.append(elemType + e + ((status & SET31)!=0&&!finalArrayPg ? " [M]" : ""));
				sb.append(", nodes: (");
				
				String nX0pfx = ((nX0&SET31)!=0?"s":"")+((nX0&SET30)!=0?"i":""), nX1pfx = ((nX1&SET31)!=0?"s":"")+((nX1&SET30)!=0?"i":"");
				String nX2pfx = ((nX2&SET31)!=0?"s":"")+((nX2&SET30)!=0?"i":""), nX3pfx = ((nX3&SET31)!=0?"s":"")+((nX3&SET30)!=0?"i":"");
				String n0pfx = ((n0&SET31)!=0?"s":"")+((n0&SET30)!=0?"i":""), n1pfx = ((n1&SET31)!=0?"s":"")+((n1&SET30)!=0?"i":"");
				String n2pfx = ((n2&SET31)!=0?"s":"")+((n2&SET30)!=0?"i":""), n3pfx = ((n3&SET31)!=0?"s":"")+((n3&SET30)!=0?"i":"");
				sb.append((ref0?nX0pfx+(nX0&0x3FFFFFFF)+"[E"+eR0+"|n"+n0+"]":n0pfx+(n0&0x3FFFFFFF)) + ", ");
				sb.append((ref1?nX1pfx+(nX1&0x3FFFFFFF)+"[E"+eR1+"|n"+n1+"]":n1pfx+(n1&0x3FFFFFFF)) + ", ");
				sb.append((ref2?nX2pfx+(nX2&0x3FFFFFFF)+"[E"+eR2+"|n"+n2+"]":n2pfx+(n2&0x3FFFFFFF)) + ", ");
				sb.append((ref3?nX3pfx+(nX3&0x3FFFFFFF)+"[E"+eR3+"|n"+n3+"]":n3pfx+(n3&0x3FFFFFFF)) + "), ");
				sb.append("s.code: "+(finalArrayPg?slotcodeName[status&0xF]:slotcodeName[status>>12&0xF])+((status&SET31)!=0?"[CA]":"")+", ");
				if (!finalArrayPg) sb.append("e.code: " + edgeCodeName[status>>16 & 0xF] + "\n"); else sb.append("\n");
				
				if (nodeA != null) {
					int m = printSnapNodes ? 6 : 3;
					int n0_3 = (ref0 ? nX0&0x3FFFFFFF : n0&0x3FFFFFFF) * m, n1_3 = (ref1 ? nX1&0x3FFFFFFF : n1&0x3FFFFFFF) * m;
					int n2_3 = (ref2 ? nX2&0x3FFFFFFF : n2&0x3FFFFFFF) * m, n3_3 = (ref3 ? nX3&0x3FFFFFFF : n3&0x3FFFFFFF) * m;
					if (n0_3<nodeA.length) sb.append(String.format("      n0: %.3f,%.3f,%.3f, ", nodeA[n0_3++], nodeA[n0_3++], nodeA[n0_3++]));
					else sb.append("      n0: n/a, ");
					if (n1_3<nodeA.length) sb.append(String.format(      "n1: %.3f,%.3f,%.3f, ", nodeA[n1_3++], nodeA[n1_3++], nodeA[n1_3++]));
					else sb.append("      n1: n/a, ");
					if (n2_3<nodeA.length) sb.append(String.format(      "n2: %.3f,%.3f,%.3f, ", nodeA[n2_3++], nodeA[n2_3++], nodeA[n2_3++]));
					else sb.append("      n2: n/a, ");
					if (n3_3<nodeA.length) sb.append(String.format(      "n3: %.3f,%.3f,%.3f\n", nodeA[n3_3++], nodeA[n3_3++], nodeA[n3_3++]));
					else sb.append("      n3: n/a, ");
				}
				if (!elemInternal && printSnapNodes) {
					int n0_3 = (ref0 ? nX0&0x3FFFFFFF : n0&0x3FFFFFFF) * 6 + 3, n1_3 = (ref1 ? nX1&0x3FFFFFFF : n1&0x3FFFFFFF) * 6 + 3;
					int n2_3 = (ref2 ? nX2&0x3FFFFFFF : n2&0x3FFFFFFF) * 6 + 3, n3_3 = (ref3 ? nX3&0x3FFFFFFF : n3&0x3FFFFFFF) * 6 + 3;
					if (!(	nodeIST[n0_3]==0&&nodeIST[n0_3+1]==0&&nodeIST[n0_3+2]==0&&nodeIST[n1_3]==0&&nodeIST[n1_3+1]==0&&nodeIST[n1_3+2]==0&&
							nodeIST[n2_3]==0&&nodeIST[n2_3+1]==0&&nodeIST[n2_3+2]==0&&nodeIST[n3_3]==0&&nodeIST[n3_3+1]==0&&nodeIST[n3_3+2]==0)) {
						if (!(nodeIST[n0_3]==0&&nodeIST[n0_3+1]==0&&nodeIST[n0_3+2]==0))
							sb.append(String.format("snap: n0: %.3f,%.3f,%.3f, ",nodeIST[n0_3++],nodeIST[n0_3++],nodeIST[n0_3]));
						else sb.append("snap:                        ");
						if (!(nodeIST[n1_3]==0&&nodeIST[n1_3+1]==0&&nodeIST[n1_3+2]==0))
							sb.append(String.format("n1: %.3f,%.3f,%.3f, ",nodeIST[n1_3++],nodeIST[n1_3++],nodeIST[n1_3]));
						else sb.append("                       ");
						if (!(nodeIST[n2_3]==0&&nodeIST[n2_3+1]==0&&nodeIST[n2_3+2]==0))
							sb.append(String.format("n2: %.3f,%.3f,%.3f, ",nodeIST[n2_3++],nodeIST[n2_3++],nodeIST[n2_3]));
						else sb.append("                       ");
						if (!(nodeIST[n3_3]==0&&nodeIST[n3_3+1]==0&&nodeIST[n3_3+2]==0))
							sb.append(String.format("n3: %.3f,%.3f,%.3f\n",nodeIST[n3_3++],nodeIST[n3_3++],nodeIST[n3_3]));
						else sb.append("\n");
					}
				}
				if (!elemInternal && !finalArrayPg) {
					sb.append("       e01   e02   e03   e12   e13   e23\n      ");
					for (int p = 0; p < 6; p++) {
						sb.append((pattern&3)==0?"-=":((pattern&3)==1?"+=":"0=")); pattern >>= 2;
						sb.append((pattern&1)!=0?"o":"="); pattern >>= 1;
						sb.append((pattern&3)==0?"=- ":((pattern&3)==1?"=+ ":"=0 ")); pattern >>= 2;
					} 
					sb.append("\n       ");
					for (int s = 0; s < 6; s++) {
						sb.append(((status&1)==0?" ":"U") + (cutNodes!=null&&cutNodes[e*6+s]!=0?"C":" ") + ((status&2)==0?"    ":"B   "));
						status >>= 2; }
					sb.append("\n");
				}
				sb.append("\n");
			}
		}
		
		int printElems = maxPrintElements < elements2 ? maxPrintElements : elements2;
		for (int e = 0; e < printElems; e++) {
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
	
	
	
	public void toOBJ(boolean IST, boolean stream) {
		File file = new File(name + "_tet.obj");
		if (!file.exists()) { try {	file.createNewFile(); } catch (IOException e) { e.printStackTrace(); }}
		BufferedWriter bw = null;
		try { bw = new BufferedWriter(new FileWriter(file));
		} catch (IOException e) { e.printStackTrace(); }
		StringBuilder sb = new StringBuilder();
		
		sb.append("# Visualisation of Isosurface Stuffing Tetrahedral mesh volume\n#\n");
		String precFormat = "%." + precision_OBJ + "f";
		int[] nodeIdx = null;
		
		if (IST) {
			double[] nA;
			if (stream)	{
				//nA = elementStreamToArray(elementStream, elementStreamCode, elements, null, FALSE);	// temporarily decompact element stream for writeout
				nodeIdx = new int[elements * 4];
				nA = elementStreamToArray(elementStream, elementStreamCode, elements, nodeIdx, true);
				for (int e = 0, e3 = 0; e < elements; e++) {
					String n0S = String.format(precFormat,nA[e3++])+" "+String.format(precFormat,nA[e3++])+" "+String.format(precFormat,nA[e3++]);
					String n1S = String.format(precFormat,nA[e3++])+" "+String.format(precFormat,nA[e3++])+" "+String.format(precFormat,nA[e3++]);
					String n2S = String.format(precFormat,nA[e3++])+" "+String.format(precFormat,nA[e3++])+" "+String.format(precFormat,nA[e3++]);
					String n3S = String.format(precFormat,nA[e3++])+" "+String.format(precFormat,nA[e3++])+" "+String.format(precFormat,nA[e3++]);
					sb.append("v "+n0S+"\nv "+n1S+"\nv "+n2S+"\nv "+n3S+"\n"); }
			} else {
				nA = nodeIST;
				for (int e = 0, e4 = 0; e < elements; e++) {
					int n03 = element[e4++]*3, n13 = element[e4++]*3, n23 = element[e4++]*3, n33 = element[e4++]*3;
					String n0S = String.format(precFormat,nA[n03++])+" "+String.format(precFormat,nA[n03++])+" "+String.format(precFormat,nA[n03]);
					String n1S = String.format(precFormat,nA[n13++])+" "+String.format(precFormat,nA[n13++])+" "+String.format(precFormat,nA[n13]);
					String n2S = String.format(precFormat,nA[n23++])+" "+String.format(precFormat,nA[n23++])+" "+String.format(precFormat,nA[n23]);
					String n3S = String.format(precFormat,nA[n33++])+" "+String.format(precFormat,nA[n33++])+" "+String.format(precFormat,nA[n33]);
					sb.append("v "+n0S+"\nv "+n1S+"\nv "+n2S+"\nv "+n3S+"\n"); }
			}
		} else {
			if (nodeWork != null) mergeNodeWork();
			for (int e = 0, e4 = 0; e < elements; e++) {
				int n03 = element[e4++]*3, n13 = element[e4++]*3, n23 = element[e4++]*3, n33 = element[e4++]*3;
				String n0S = String.format(precFormat,node[n03++])+" "+String.format(precFormat,node[n03++])+" "+String.format(precFormat,node[n03]);
				String n1S = String.format(precFormat,node[n13++])+" "+String.format(precFormat,node[n13++])+" "+String.format(precFormat,node[n13]);
				String n2S = String.format(precFormat,node[n23++])+" "+String.format(precFormat,node[n23++])+" "+String.format(precFormat,node[n23]);
				String n3S = String.format(precFormat,node[n33++])+" "+String.format(precFormat,node[n33++])+" "+String.format(precFormat,node[n33]);
				sb.append("v "+n0S+"\nv "+n1S+"\nv "+n2S+"\nv "+n3S+"\n");
			}
		}
		
		sb.append("# " + 12*elements + " vertices\n\ng " + name + "_tet\n");
		
		int n0=0, n1=0, n2=0, n3=0;
		if (IST) {
			for (int e = 0, e4 = 1; e < elements; e++, e4+=4) {
				sb.append(	  "f "+(e4)+  " "+(e4+1)+" "+(e4+2)+									// write out tetrahedral facets
							"\nf "+(e4)+  " "+(e4+2)+" "+(e4+3)+
							"\nf "+(e4)+  " "+(e4+3)+" "+(e4+1)+
							"\nf "+(e4+1)+" "+(e4+3)+" "+(e4+2)+"\n"); }
		} else {
			if (element2 != null) {
				if (element2Work != null) mergeElement2Work();
				for (int e = 0, e4 = 0; e < elements2; e++) {
					n0 = element2[e4].nodeRef[0]+1; n1 = element2[e4].nodeRef[1]+1; n2 = element2[e4].nodeRef[2]+1; n3 = element2[e4++].nodeRef[3]+1;
					// write out tetrahedral facets
					sb.append("f "+n0+" "+n1+" "+n2+"\nf "+n0+" "+n3+" "+n1+"\nf "+n0+" "+n2+" "+n3+"\nf "+n1+" "+n3+" "+n2+"\n");
				}
			} else {
				for (int e = 0, e4 = 0; e < elements; e++) {
					n0 = element[e4++]+1; n1 = element[e4++]+1; n2 = element[e4++]+1; n3 = element[e4++]+1;
					// write out tetrahedral facets
					sb.append("f "+n0+" "+n1+" "+n2+"\nf "+n0+" "+n3+" "+n1+"\nf "+n0+" "+n2+" "+n3+"\nf "+n1+" "+n3+" "+n2+"\n");
				}
			}
		}
		sb.append("# " + (IST ? elements*4 : (elements2>0 ? elements2*4 : elements*4)) + " faces\n\ng\n");

		try {	bw.write(sb.toString());									// write out & close file
		bw.flush(); bw.close();
		} catch (IOException e) { e.printStackTrace(); }
	}
	
	
	public String toStringSlotcodes() {
		StringBuilder sb = new StringBuilder();
		if (printSlotcodes2) {
			long slotcodes = elementStreamCode[0]; 
			sb.append(slotcodeName[(int)slotcodes & 0xF] + " "); slotcodes>>=4;
			for (int c = 0, cEnd = elementStreamCode.length-1, c16 = 15; c < cEnd;) {
				sb.append(slotcodeName[(int)slotcodes & 0xF] + " ");
				if (--c16 <= 0) { slotcodes = elementStreamCode[++c]; c16 = 16; sb.append("\n"); } else slotcodes>>=4;
			}
		}
		if (printSlotcodes1) {
			int[] elemPage = elementPg[0];
			for (int e = 1, e5 = 4, ePage = 0; e <= elements; e++, e5+=5) {
				if (e5 >= elementPgSize) { e5 = 4; elemPage = elementPg[++ePage]; }
				sb.append(slotcodeName[elemPage[e5] & 0xF] + ((e&15)==0 ? "\n" : " "));
			}
		}
		return sb.toString();
	}
	
	
	// method provies a visual check on closestNodePairs() method, exporting the result to OBJ 3D-file
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
	
	public static String toStringBCCtetraInfo() {
		double[] tetBCCcoords = {	0,1,0.5, 0.5,0.5,0, -0.5,0.5,0, 0,0,0.5,		0,1,0.5, 0.5,0.5,0, 0,0.5,0, 0,0,0.5,
									0,0.5,0.5, 0.5,0.5,0, -0.5,0.5,0, 0,0,0.5,		0,0.5,0.5, 0.5,0.5,0, 0,0.5,0, 0,0,0.5,
									0,0.5,0.5, 0.5,0,0, -0.5,1,0, -0.5,0,0};
		String tetraInfo =	"BCC tetrahedron type 1 quality " + tetraQuality(tetBCCcoords, 0, 1, 2, 3)
							+ ", smallest dihedral angle " + tetraSmallestDihedral(tetraDihedralAngles(tetBCCcoords, 0, 1, 2, 3)) + "\n";
		tetraInfo += 		"BCC tetrahedron type 2 quality " + tetraQuality(tetBCCcoords, 4, 5, 6, 7)
							+ ", smallest dihedral angle " + tetraSmallestDihedral(tetraDihedralAngles(tetBCCcoords, 4, 5, 6, 7)) + "\n";
		tetraInfo += 		"BCC tetrahedron type 3 quality " + tetraQuality(tetBCCcoords, 8, 9, 10, 11)
							+ ", smallest dihedral angle " + tetraSmallestDihedral(tetraDihedralAngles(tetBCCcoords, 8, 9, 10, 11)) + "\n";
		tetraInfo += 		"BCC tetrahedron type 4 quality " + tetraQuality(tetBCCcoords, 12, 13, 14, 15)
							+ ", smallest dihedral angle " + tetraSmallestDihedral(tetraDihedralAngles(tetBCCcoords, 12, 13, 14, 15)) + "\n";
		return tetraInfo +	"BCC tetrahedron type 5 quality " + tetraQuality(tetBCCcoords, 16, 17, 18, 19)
							+ ", smallest dihedral angle " + tetraSmallestDihedral(tetraDihedralAngles(tetBCCcoords, 16, 17, 18, 19)) + "\n";		
	}

}
