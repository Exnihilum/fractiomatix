package lkr74.fem1;

import lkr74.matrixlib.Matrix;

public class FEM1Element {
	// the flags parameter is set to one of these defintions
	public static final int TETRAHEDRON = 0, HEXAHEDRON = 1, ISOTETRAHEDRON = 2, ISOHEXAHEDRON = 3;
	// the named offsets into a tetrahedron's data[] array
	public static final int TD_EDGE01 = 0, TD_EDGE02 = 1, TD_EDGE03 = 2, TD_EDGE12 = 3, TD_EDGE13 = 4, TD_EDGE23 = 5;
	public static final int TD_AREA012 = 6, TD_AREA023 = 7, TD_AREA031 = 8, TD_AREA132 = 9;
	public static final int TD_VGRAD0 = 10, TD_VGRAD1 = 13, TD_VGRAD2 = 16, TD_VGRAD3 = 19;
	public static final int TD_NORM012 = 22, TD_NORM023 = 25, TD_NORM031 = 28, TD_NORM132 = 31;

	int nodes=0, flags=0;					// how many integration nodes, flags for this element
	int neighboursF = 0, neighbours = 0;	// no. of facet-neighbours facet-facet & total neighbours for this element
	int patch=0;							// index of the smoothing patch that this element belongs to
	int[] nodeRef;							// indexing into global coordinate list of what nodes that belong to this element
	double[] node=null, nodeWork=null;		// global nodes reference
	//double[] node0, node1, node2, node3;	// location assignments (calculation-temporary) of where the four node coordinates are stored
	double[] delta = null;					// the often calculated v1-v0, v2-v0, v3-v0, v3-v1, v2-v1 tetrahedral delta-(x,y,z) components can be stored here 
	public double[] data = null;			// edge lengths [6], areas [4], volume gradients [3x4], crossproduct normals [3x4] go here
	double[] matProperty = null;			// material properties of this element contained here
	double[][] computed = null;				// this is a link to possibly precomputed data for this type of element
	//double[] dsde = null;					// the precomputed stress/strain derivative for this element type
	double volume=0, iRadius=0, cRadius=0;	// volume, inscribed radius & circumscribed radius of element
	int[] neighbour;						// holds neighbour elements
	
	static byte[] iSortSlot = new byte[513];
	{ iSortSlot[2]=1; iSortSlot[4]=2; iSortSlot[8]=3; iSortSlot[16]=4; iSortSlot[32]=5; iSortSlot[64]=6; iSortSlot[128]=7; iSortSlot[256]=8; iSortSlot[512]=9; }
	static byte[] iAreaSlot = new byte[9]; { iAreaSlot[1]=TD_AREA012; iAreaSlot[2]=TD_AREA023; iAreaSlot[4]=TD_AREA031; iAreaSlot[8]=TD_AREA132; }
	static byte[] iAreaSlot2 = new byte[9]; { iAreaSlot2[1]=0; iAreaSlot2[2]=1; iAreaSlot2[4]=2; iAreaSlot2[8]=3; }
	static byte[] iEdgeSlot = new byte[513];
	{ iEdgeSlot[16]=TD_EDGE01; iEdgeSlot[32]=TD_EDGE02; iEdgeSlot[64]=TD_EDGE03; iEdgeSlot[128]=TD_EDGE12; iEdgeSlot[256]=TD_EDGE13; iEdgeSlot[512]=TD_EDGE23; }
	
	// instantiates a tetrahedron
	public FEM1Element(FEM1 fem, int cRef1, int cRef2, int cRef3, int cRef4, double[] matProperty) {
		this.nodes = 4;
		flags = TETRAHEDRON;
		this.nodeRef = new int[4];
		nodeRef[0] = cRef1; nodeRef[1] = cRef2; nodeRef[2] = cRef3; nodeRef[3] = cRef4;
		node = fem.node; nodeWork = fem.nodeWork;
		this.matProperty = matProperty == null ? null : matProperty.clone();
	}

	public int elementType() { return flags & 15; }
	public void flagHasDeltas() { flags |= 16; }
	public void clearHasDeltas() { flags &= 0xFFFFFFFF - 16; }
	public boolean hasDeltas() { return (flags & 16) != 0; }
	public void flagHasAreas() { flags |= 32; }
	public void clearHasAreas() { flags &= 0xFFFFFFFF - 32; }
	public boolean hasAreas() { return (flags & 32) != 0; }
	public void flagHasEdges() { flags |= 64; }
	public void clearHasEdges() { flags &= 0xFFFFFFFF - 64; }
	public boolean hasEdges() { return (flags & 64) != 0; }
	public void flagHasVGradients() { flags |= 128; }
	public void clearHasVGradients() { flags &= 0xFFFFFFFF - 128; }
	public boolean hasVGradients() { return (flags & 128) != 0; }
	public void flagHasInterfaces() { flags |= 256; }
	public void clearHasInterfaces() { flags &= 0xFFFFFFFF - 256; }
	public boolean hasInterfaces() { return (flags & 256) != 0; }
	public void flagHasXNormals() { flags |= 512; }
	public void clearHasXNormals() { flags &= 0xFFFFFFFF - 512; }
	public boolean hasXNormals() { return (flags & 512) != 0; }
	public void flagHasNormals() { flags |= 1024; }
	public void clearHasNormals() { flags &= 0xFFFFFFFF - 1024; }
	public boolean hasNormals() { return (flags & 1024) != 0; }
	// these bits will flag existence of individual normals
	final static int FLAG_NORM012 = 2048, FLAG_NORM023 = 4096, FLAG_NORM031 = 8192, FLAG_NORM132 = 16384;
	final static int FLAG_NORMALS = FLAG_NORM012 + FLAG_NORM023 + FLAG_NORM031 + FLAG_NORM132;
		
	
	// method adds fields to the precalculated data array, depending on what field is requested
	// the algorithm blocks are naturally arranged according to how far offsetted the new data is in array
	final static int EDGE_FIELD = 1, AREA_FIELD = 2, VOLGRADIENT_FIELD = 4, NORMAL_FIELD = 8;
	void setupData(int field) {
		if ((field & NORMAL_FIELD) != 0) {								// case of request for a normals datafield
			if (data == null) { data = new double[TD_NORM132 + FEM1.NCOORD]; return; }
			if (data.length <= 6+4+3*4) {								// is data[] short of a normals datafield?
				double[] dataNew = new double[TD_NORM132 + FEM1.NCOORD];// copy over the eventual edge, area & gradients data
				for (int d = 0; d <= TD_NORM132; d++) dataNew[d]=data[d];
				data = dataNew; return;
			}
			return;
		}

		if ((field & VOLGRADIENT_FIELD) != 0) {							// case of request for volume gradient datafield
			if (data == null) { data = new double[TD_VGRAD3 + FEM1.NCOORD]; return; }
			if (data.length <= 6 + 4) {									// is data[] short of a volume gradient datafield?
				double[] dataNew = new double[TD_VGRAD3 + FEM1.NCOORD];	// copy over the eventual edge & area data
				dataNew[0]=data[0]; dataNew[1]=data[1]; dataNew[2]=data[2]; dataNew[3]=data[3]; dataNew[4]=data[4]; dataNew[5]=data[5];
				dataNew[6]=data[6]; dataNew[7]=data[7]; dataNew[8]=data[8]; dataNew[9]=data[9];
				data = dataNew; return;
			}
			return;
		}
		
		if ((field & (EDGE_FIELD + AREA_FIELD)) != 0) {					// case of request for edge or area datafields
			if (data == null) { data = new double[6 + 4]; return; }
			if (data.length > 6 + 4) {									// do we have volume gradient data?
				double[] dataNew = new double[TD_AREA132 + 1];			// copy over the eventual volume gradient data
				dataNew[10]=data[10]; dataNew[11]=data[11]; dataNew[12]=data[12]; dataNew[13]=data[13]; dataNew[14]=data[14]; dataNew[15]=data[15];
				dataNew[16]=data[16]; dataNew[17]=data[17]; dataNew[18]=data[18]; dataNew[19]=data[19]; dataNew[20]=data[20]; dataNew[21]=data[21];
				data = dataNew; return;
			}
		}
	}
	
	
	// method copies over element's calculated boundary-shared datafields to neighbours of element if absorb = false
	// if absorb = true, element absorbs data from the neighbours instead
	// TODO: debug the absorption additions
	public final static int PROPAGATE_AREAS = 1, PROPAGATE_EDGES = 2, PROPAGATE_NORMALS = 4;
	public void propagateInterfaces(FEM1 fem, int propType, boolean absorb) {
		if (data == null) return;																// sanity check: do we have calculated data at all?
		switch (flags & 15) {
		case TETRAHEDRON:
			if (!hasInterfaces()) return;														// sanity check: does this element have interfaces?
			if ((propType & PROPAGATE_AREAS) != 0 && !(absorb && hasAreas())) {
				
				if (absorb) setupData(AREA_FIELD);
				for (int i = 0, neighboursF2 = neighboursF * 2; i < neighboursF2; i += 2) {				
					FEM1Element elem2 = fem.getElement2(neighbour[i]);
					int intf_flag = neighbour[i + 1];
					if (absorb) {
						double area2 = elem2.data[iAreaSlot[(intf_flag >> 16) & 15]];
						if (area2 > 0) data[iAreaSlot[intf_flag & 15]] = area2;					// does neighbour have a calculated area to share?
							
						if (data[TD_AREA012] > 0 && data[TD_AREA023] > 0 && data[TD_AREA031] > 0 && data[TD_AREA132] > 0) {
							flagHasAreas(); break; } // if all areas set, we're done absorbing
					} else {
						if (elem2.hasAreas()) continue;											// neighbour's areas are up to date
						double area1 = data[iAreaSlot[intf_flag & 15]];
						if (area1 > 0) {														// do we have a calculated area to share?
							elem2.setupData(AREA_FIELD);										// check if area field needs setting up in neighbour
							double[] data2 = elem2.data;
							data2[iAreaSlot[(intf_flag >> 16) & 15]] = area1;
							if (data2[TD_AREA012] > 0 && data2[TD_AREA023] > 0 && data2[TD_AREA031] > 0 && data2[TD_AREA132] > 0) elem2.flagHasAreas();
						}
					}
				}
			}
			if ((propType & PROPAGATE_NORMALS) != 0 && !(absorb && hasNormals())) {
				
				if (absorb) setupData(NORMAL_FIELD);
				for (int i = 0, neighboursF2 = neighboursF * 2; i < neighboursF2; i += 2) {				
					FEM1Element elem2 = fem.getElement2(neighbour[i]);
					if (absorb) {
						if (hasNormals() || hasXNormals() && elem2.hasXNormals()) continue;		// normals are up-to-date to at least same level as neighbour
						int intf_flag = neighbour[i + 1], dOfs2 = TD_NORM012 + iAreaSlot[(intf_flag >> 16) & 15];
						
						double normal2x = elem2.data[dOfs2++], normal2y = elem2.data[dOfs2++], normal2z = elem2.data[dOfs2];
						if (normal2x != 0 || normal2y != 0 || normal2z != 0) {					// does neighbour have a calculated normal to share?
							int dOfs1 = TD_NORM012 + iAreaSlot[intf_flag & 15];
							switch (dOfs1) {case TD_NORM012: flags|=FLAG_NORM012; break; case TD_NORM023: flags|=FLAG_NORM023; 
											case TD_NORM031: flags|=FLAG_NORM031; break; case TD_NORM132: flags|=FLAG_NORM132;}
							data[dOfs1++] = normal2x; data[dOfs1++] = normal2y; data[dOfs1] = normal2z;
						}
							
						if ((flags & FLAG_NORMALS) == FLAG_NORMALS) {
							// since both normals and crossp.normals might be absorbed, we can only be sure that we've got upto crossp.normals
							flagHasXNormals(); break;
						}
					} else {
						if (elem2.hasNormals() || elem2.hasXNormals() && hasXNormals()) continue;	// neighbour's normals are up-to-date to same level
						int intf_flag = neighbour[i + 1], dOfs1 = TD_NORM012 + iAreaSlot[intf_flag & 15];
						
						double normal1x = data[dOfs1++], normal1y = data[dOfs1++], normal1z = data[dOfs1];
						if (normal1x != 0 || normal1y != 0 || normal1z != 0) {						// do we have a calculated normal to share?
							
							elem2.setupData(NORMAL_FIELD);											// check if normal field needs setting up in neighbour
							int dOfs2 = TD_NORM012 + iAreaSlot[(intf_flag >> 16) & 15];
							switch (dOfs2) {case TD_NORM012: elem2.flags|=FLAG_NORM012; break; case TD_NORM023: elem2.flags|=FLAG_NORM023; 
											case TD_NORM031: elem2.flags|=FLAG_NORM031; break; case TD_NORM132: elem2.flags|=FLAG_NORM132;}
							elem2.data[dOfs2++] = normal1x; elem2.data[dOfs2++] = normal1y; elem2.data[dOfs2] = normal1z;
						}
						if ((elem2.flags & FLAG_NORMALS) == FLAG_NORMALS) elem2.flagHasXNormals();
					}
				}
			}
			if ((propType & PROPAGATE_EDGES) != 0 && !(absorb && hasEdges())) {
				
				if (absorb) setupData(EDGE_FIELD);														// make sure this element has an edge datafield
				int neighboursF2 = neighboursF * 2;
				for (int i = 0, neighbours2 = neighbours * 2; i < neighbours2; i += 2) {				
					int intf_flag = neighbour[i + 1];
					FEM1Element elem2 = fem.getElement2(neighbour[i]);
					if (!absorb) {
						if (elem2.hasEdges()) continue;												// if neighbour's edges are up to date
						elem2.setupData(EDGE_FIELD);												// allocate for edges & facet areas if necessary
					}
					
					// process facet-facet first (a facet-facet interface is equivalent to three edge-edge interfaces)
					if (i < neighboursF2) {
						
						int o1 = iAreaSlot2[intf_flag & 15], e2f = iAreaSlot2[(intf_flag>>16) & 15], o2 = e2f * 3;	// get indexes of faces from status bits
						final int[] e1ni = { 0, 0, 0, 1 }, e2ni = { 0,1,2, 0,2,3, 0,3,1, 1,3,2 };	// initialise indirect face->edge addressing arrays
						final int[] edges1 = { 0,3,1,0,3, 1,5,2,1,5, 2,4,0,2,4, 4,5,3,4,5 };
						final int[] edges2 = { 1,3,0,1,3, 2,5,1,2,5, 0,4,2,0,4, 3,5,4,3,5 };

						int[] r2 = elem2.nodeRef;
						if (nodeRef[e1ni[o1]] == r2[e2ni[o2]]) 	o2 = e2f * 5; else
						if (nodeRef[e1ni[o1]] == r2[e2ni[++o2]]) o2 = e2f * 5 + 2; else
						/*if (nodeRef[e1ni[o1]] == r2[e2ni[++o2]])*/ o2 = e2f * 5 + 1;
						o1 *= 5;
						 double eL1 = 0;
						if (absorb) {
									if ((eL1 = elem2.data[edges2[o2++]]) != 0)	data[edges1[o1]] = eL1;	// for every edge of found facet, check if
							o2++;	if ((eL1 = elem2.data[edges2[o2++]]) != 0)	data[edges1[o1]] = eL1;	// it's calculated before propagation
							o2++;	if ((eL1 = elem2.data[edges2[o2]]) != 0)	data[edges1[o1]] = eL1;
						} else {
							double[] edgeLength2 = elem2.data;
									if ((eL1 = data[edges1[o1++]]) != 0)	edgeLength2[edges2[o2]] = eL1;
							o2++;	if ((eL1 = data[edges1[o1++]]) != 0)	edgeLength2[edges2[o2]] = eL1;
							o2++;	if ((eL1 = data[edges1[o1]]) != 0)		edgeLength2[edges2[o2]] = eL1;
						}
					} else {
						if (absorb) {
							double edgeL = elem2.data[iEdgeSlot[((intf_flag >> 16) & 1008)]];
							if (edgeL > 0)															// does neighbour have a calculated edge length to share?
								data[iEdgeSlot[(intf_flag & 1008)]] = edgeL;
						} else {
							double edgeL = data[iEdgeSlot[(intf_flag & 1008)]];
							if (edgeL > 0)															// do we have a calculated edge length to share?
								elem2.data[iEdgeSlot[((intf_flag >> 16) & 1008)]] = edgeL;
						}
					}
//					double[] data2 = elem2.data;
//					if (data2[0] > 0 && data2[1] > 0 && data2[2] > 0 && data2[3] > 0 && data2[4] > 0 && data2[5] > 0) elem2.flagHasEdges();
				}
			}
			break;
		case ISOTETRAHEDRON:
		case HEXAHEDRON:
		case ISOHEXAHEDRON:
		}
	}	

	
	
	// method sorts element's edge-edge interfaces grouping them according to the edge index
	void sortInterfaces() {
		//if (!hasInterfaces()) return;									// DEBUG: caller method assumed to test if element has interfaces
		switch (flags & 15) {
		case TETRAHEDRON:
			int[][] sorter = new int[10][neighbours * 2];
			int[] eSCount = {0,0,0,0,0,0,0,0,0,0};

			for (int ngb = 0, neighbours2 = neighbours * 2; ngb < neighbours2;) {
				// spread out incoming neighbour interfaces into 10 indexed temporary arrays, for later sorted recollection
				int slot = iSortSlot[neighbour[ngb + 1] & 65535];
				sorter[slot][eSCount[slot]++] = neighbour[ngb++];
				sorter[slot][eSCount[slot]++] = neighbour[ngb++];
			}
			for (int i = 0, e = 0; e < 10; e++)							// regather neighbours in sorted order
				for (int eS = 0, eSC = eSCount[e]; eS < eSC;) {
					neighbour[i++] = sorter[e][eS++];
					neighbour[i++] = sorter[e][eS++];
				}
			break;
		case ISOTETRAHEDRON:
		case HEXAHEDRON:
		case ISOHEXAHEDRON:
		}
	}

	
	// method either inserts a neighbour/interface pair in sorted order or, if element lacks interfaces, inserting just the neighbour
	// in the second case, flags = 0 means a facet neighbour inserted, flags = 1 means an edge neighbour
	public void insertNeighbour(int e2, int flags) {
		boolean hasInterfaces = hasInterfaces();
		int[] neighbourNew;
		// only way this case is enacted is for a free-floating neighbourless tetrahedron
		if (neighbour == null) neighbour = neighbourNew = new int[hasInterfaces ? 8 * 2 : 8];
		else if (neighbours + (hasInterfaces ? 2 : 1) > neighbour.length)			// existing array can't fit another neighbour?
			neighbourNew = new int[neighbour.length + (hasInterfaces ? 8 * 2 : 8)];	// optimisation: allocate for 8 additional neighbours
		else neighbourNew = new int[neighbour.length];								// we CAN fit interface into current array length

		if (hasInterfaces) {														// this block will insert a neighbour & interface pair
			if (neighbour == null) neighbour = neighbourNew = new int[8 * 2];
			else if (neighbours + 2 > neighbour.length)								// existing array can't fit another neighbour?
				neighbourNew = new int[neighbour.length + 8 * 2];					// optimisation: allocate for 8 additional neighbours
			else neighbourNew = new int[neighbour.length];							// we CAN fit interface into current array length

			int flags1 = flags & 65535;
			for (int ngb = 0, ngb2 = 0, neighbours2 = neighbours * 2; ngb < neighbours2;) {
				int flags2 = (neighbour[ngb + 1] & 65535);
				// TODO: an element can currently only have one interface towards another element, on linear element assumption
				if (/*flags1 == flags2 && */neighbour[ngb] == e2) return;			// if that neighbour & interface already exist, do nothing more
				else if (flags1 > flags2) {											// if reached slot that has a lower flag value than inserted flag
					neighbourNew[ngb2++] = e2; neighbourNew[ngb2++] = flags;		// do the insertion here
				}
				neighbourNew[ngb2++] = neighbour[ngb++];							// regular pair copying
				neighbourNew[ngb2++] = neighbour[ngb++];
			}
			if ((flags & 15) != 0) neighboursF++;									// if it was a facet interface, increase facet count
			
		} else {																	// this block only inserts neighbour indexes
			if (neighbour == null) neighbour = neighbourNew = new int[8];
			else if (neighbours + 1 > neighbour.length)								// existing array can't fit another neighbour?
				neighbourNew = new int[neighbour.length + 8];						// optimisation: allocate for 8 additional neighbours
			else neighbourNew = new int[neighbour.length];							// we CAN fit interface into current array length

			for (int ngb = 0, ngb2 = 0; ngb < neighbours; ngb++) {
				if (neighbour[ngb] == e2) return;									// if neighbour already exists, do nothing more
				if (flags == 0 && ngb2 == neighboursF)								// if we're inserting a facet neighbour & we're past last facet neighbour
					neighbourNew[ngb2++] = e2;										// insert
				neighbourNew[ngb2++] = neighbour[ngb];								// regular copying
			}
			if (flags == 0) neighboursF++;											// if it was a facet neighbour, increase facet count
		}
		neighbours++;
		neighbour = neighbourNew;
	}
	
	
	
	public void deleteNeighbour(int e2) {
		boolean foundNeighbour = false;
		if (hasInterfaces()) {														// this block will delete a neighbour & interface pair
			int ngb = 0, neighbours2 = neighbours * 2;
			for (;ngb < neighbours2; ngb += 2)
				if (neighbour[ngb] == e2) { foundNeighbour = true; break; }			// found neighbour? Break out to next loop
			
			if (!foundNeighbour) return;
			if ((ngb>>1) < neighboursF) neighboursF--;								// if it was a facet naighbour, decrease facet count
			int ngb2 = ngb + 2;
			while (ngb2 < neighbours2) {
				neighbour[ngb++] = neighbour[ngb2++];								// shift neighbours one step down
				neighbour[ngb++] = neighbour[ngb2++];
			}					
		} else {																	// this block only inserts neighbour indexes
			int ngb = 0;
			for (; ngb < neighbours; ngb++)
				if (neighbour[ngb] == e2) { foundNeighbour = true; break; }			// found neighbour? Break out to next loop

			if (!foundNeighbour) return;
			if ((ngb>>1) < neighboursF) neighboursF--;								// if it was a facet naighbour, decrease facet count
			int ngb2 = ngb + 1;
			while (ngb2 < neighbours) neighbour[ngb++] = neighbour[ngb2++];			// shift neighbours one step down
		}
		neighbours--;
	}

	
	
	// method, provided with two node indexes and a neighbour tetrahedron's node reference array, returns edge-edge bitshifted interface bitflags
	// note: calling method has already found that elements are edge-neighbours, so a zero shouldn't be returnable
	final static int NFF012=1, NFF023=2, NFF031=4, NFF132=8, NEE01=16, NEE02=32, NEE03=64, NEE12=128, NEE13=256, NEE23=512;
	public int tetraEdgeInterface(int rA, int rB, int[] nodeRef2) {
		int nA = nodeRef[rA], nB = nodeRef[rB], fl1 = 0;
		switch (rA) {
		case 0: switch(rB) {case 1: fl1 = NEE01; break; case 2: fl1 = NEE02; break; case 3: fl1 = NEE03;} break;
		case 1: switch(rB) {case 0: fl1 = NEE01; break; case 2: fl1 = NEE12; break; case 3: fl1 = NEE13;} break;
		case 2: switch(rB) {case 0: fl1 = NEE02; break; case 1: fl1 = NEE12; break; case 3: fl1 = NEE23;} break;
		case 3: switch(rB) {case 0: fl1 = NEE03; break; case 1: fl1 = NEE13; break; case 2: fl1 = NEE23;} break;
		}
		if (nA == nodeRef[0]) {
			if (nB==nodeRef[1]) return fl1|(NEE01<<16); else if (nB==nodeRef[2]) return fl1|(NEE02<<16); else /*if (nB==nodeRef[3])*/ return fl1|(NEE03<<16); //*else return -1;
		} else if (nA == nodeRef[1]) {
			if (nB==nodeRef[0]) return fl1|(NEE01<<16); else if (nB==nodeRef[2]) return fl1|(NEE12<<16); else /*if (nB==nodeRef[3])*/ return fl1|(NEE13<<16); //*else return -1;
		} else if (nA == nodeRef[2]) {
			if (nB==nodeRef[0]) return fl1|(NEE02<<16); else if (nB==nodeRef[1]) return fl1|(NEE12<<16); else /*if (nB==nodeRef[3])*/ return fl1|(NEE23<<16); //else return -1;
		} else if (nA == nodeRef[3]) {
			if (nB==nodeRef[0]) return fl1|(NEE03<<16); else if (nB==nodeRef[1]) return fl1|(NEE13<<16); else /*if (nB==nodeRef[2])*/ return fl1|(NEE23<<16); //else return -1;
		}
		return 0;
	}

	
	// method locates the edge-edge or facet-facet interface of this tetrahedron with node references of a known neighbour
	// result returned as bitshifted bitflags, this tetrahedron's in lower 16 bits, the neighbour's in upper 16 bits
	int tetraNeighbourInterface(int[] nodeRef2) {
		int va0 = nodeRef[0], va1 = nodeRef[1], va2 = nodeRef[2], va3 = nodeRef[3];
		int vb0 = nodeRef2[0], vb1 = nodeRef2[1], vb2 = nodeRef2[2], vb3 = nodeRef2[3], fl1 = 0, fl2 = 0;
		
		// the neighbourhood interfacing rules check for face-face interfaces, with fallback to edge-edge if two references match
		for (int i = 2; i > 0; i--) {
			if (va0==vb0) {
				if (va1==vb1) {if (va2==vb3) {fl1|=NFF012;fl2|=NFF031;} else if (va3==vb2) {fl1|=NFF031;fl2|=NFF012;} else {fl1|=NEE01;fl2|=NEE01;}} else
				if (va1==vb2) {if (va2==vb1) {fl1|=NFF012;fl2|=NFF012;} else if (va3==vb3) {fl1|=NFF031;fl2|=NFF023;} else {fl1|=NEE01;fl2|=NEE02;}} else
				if (va1==vb3) {if (va2==vb2) {fl1|=NFF012;fl2|=NFF023;} else if (va3==vb1) {fl1|=NFF031;fl2|=NFF031;} else {fl1|=NEE01;fl2|=NEE03;}} else
				if (va2==vb2) {if (va1==vb3) {fl1|=NFF012;fl2|=NFF023;} else if (va3==vb1) {fl1|=NFF023;fl2|=NFF012;} else {fl1|=NEE02;fl2|=NEE02;}} else
				if (va2==vb3) {if (va1==vb1) {fl1|=NFF012;fl2|=NFF031;} else if (va3==vb2) {fl1|=NFF023;fl2|=NFF023;} else {fl1|=NEE02;fl2|=NEE03;}} else
				if (va3==vb3) {if (va2==vb1) {fl1|=NFF023;fl2|=NFF031;} else if (va1==vb2) {fl1|=NFF031;fl2|=NFF023;} else {fl1|=NEE03;fl2|=NEE03;}}
			} else if (va0==vb1) {
				if (va1==vb0) {if (va2==vb2) {fl1|=NFF012;fl2|=NFF012;} else if (va3==vb3) {fl1|=NFF031;fl2|=NFF031;} else {fl1|=NEE01;fl2|=NEE01;}} else
				if (va1==vb2) {if (va2==vb3) {fl1|=NFF012;fl2|=NFF132;} else if (va3==vb0) {fl1|=NFF031;fl2|=NFF012;} else {fl1|=NEE01;fl2|=NEE12;}} else
				if (va1==vb3) {if (va2==vb0) {fl1|=NFF012;fl2|=NFF031;} else if (va3==vb2) {fl1|=NFF031;fl2|=NFF132;} else {fl1|=NEE01;fl2|=NEE13;}} else
				if (va2==vb2) {if (va1==vb0) {fl1|=NFF012;fl2|=NFF012;} else if (va3==vb3) {fl1|=NFF023;fl2|=NFF132;} else {fl1|=NEE02;fl2|=NEE12;}} else
				if (va2==vb3) {if (va1==vb2) {fl1|=NFF012;fl2|=NFF132;} else if (va3==vb0) {fl1|=NFF023;fl2|=NFF031;} else {fl1|=NEE02;fl2|=NEE13;}} else
				if (va3==vb2) {if (va2==vb0) {fl1|=NFF023;fl2|=NFF012;} else if (va1==vb3) {fl1|=NFF031;fl2|=NFF132;} else {fl1|=NEE03;fl2|=NEE12;}} else
				if (va3==vb3) {if (va2==vb2) {fl1|=NFF023;fl2|=NFF132;} else if (va1==vb0) {fl1|=NFF031;fl2|=NFF031;} else {fl1|=NEE03;fl2|=NEE13;}}
			} else if (va0==vb2) {
				if (va1==vb0) {if (va2==vb3) {fl1|=NFF012;fl2|=NFF023;} else if (va3==vb1) {fl1|=NFF031;fl2|=NFF012;} else {fl1|=NEE01;fl2|=NEE02;}} else
				if (va1==vb1) {if (va2==vb0) {fl1|=NFF012;fl2|=NFF012;} else if (va3==vb3) {fl1|=NFF031;fl2|=NFF132;} else {fl1|=NEE01;fl2|=NEE12;}} else
				if (va1==vb3) {if (va2==vb1) {fl1|=NFF012;fl2|=NFF132;} else if (va3==vb0) {fl1|=NFF031;fl2|=NFF023;} else {fl1|=NEE01;fl2|=NEE23;}} else
				if (va2==vb0) {if (va1==vb1) {fl1|=NFF012;fl2|=NFF012;} else if (va3==vb3) {fl1|=NFF023;fl2|=NFF023;} else {fl1|=NEE02;fl2|=NEE02;}} else
				if (va2==vb1) {if (va1==vb3) {fl1|=NFF012;fl2|=NFF132;} else if (va3==vb0) {fl1|=NFF023;fl2|=NFF012;} else {fl1|=NEE02;fl2|=NEE12;}} else
				if (va2==vb3) {if (va1==vb0) {fl1|=NFF012;fl2|=NFF023;} else if (va3==vb1) {fl1|=NFF023;fl2|=NFF132;} else {fl1|=NEE02;fl2|=NEE23;}} else
				if (va3==vb1) {if (va2==vb3) {fl1|=NFF023;fl2|=NFF132;} else if (va1==vb0) {fl1|=NFF031;fl2|=NFF012;} else {fl1|=NEE03;fl2|=NEE12;}} else
				if (va3==vb3) {if (va2==vb0) {fl1|=NFF023;fl2|=NFF023;} else if (va1==vb1) {fl1|=NFF031;fl2|=NFF132;} else {fl1|=NEE03;fl2|=NEE23;}}
			} else if (va0==vb3) {
				if (va1==vb0) {if (va2==vb1) {fl1|=NFF012;fl2|=NFF031;} else if (va3==vb2) {fl1|=NFF031;fl2|=NFF023;} else {fl1|=NEE01;fl2|=NEE03;}} else
				if (va1==vb1) {if (va2==vb2) {fl1|=NFF012;fl2|=NFF132;} else if (va3==vb0) {fl1|=NFF031;fl2|=NFF031;} else {fl1|=NEE01;fl2|=NEE13;}} else
				if (va1==vb2) {if (va2==vb0) {fl1|=NFF012;fl2|=NFF023;} else if (va3==vb1) {fl1|=NFF031;fl2|=NFF132;} else {fl1|=NEE01;fl2|=NEE23;}} else
				if (va2==vb0) {if (va1==vb2) {fl1|=NFF012;fl2|=NFF023;} else if (va3==vb1) {fl1|=NFF023;fl2|=NFF031;} else {fl1|=NEE02;fl2|=NEE03;}} else
				if (va2==vb1) {if (va1==vb0) {fl1|=NFF012;fl2|=NFF031;} else if (va3==vb2) {fl1|=NFF023;fl2|=NFF132;} else {fl1|=NEE02;fl2|=NEE13;}} else
				if (va2==vb2) {if (va1==vb1) {fl1|=NFF012;fl2|=NFF132;} else if (va3==vb0) {fl1|=NFF023;fl2|=NFF023;} else {fl1|=NEE02;fl2|=NEE23;}} else
				if (va3==vb0) {if (va2==vb2) {fl1|=NFF023;fl2|=NFF023;} else if (va1==vb1) {fl1|=NFF031;fl2|=NFF031;} else {fl1|=NEE03;fl2|=NEE03;}} else
				if (va3==vb1) {if (va2==vb0) {fl1|=NFF023;fl2|=NFF031;} else if (va1==vb2) {fl1|=NFF031;fl2|=NFF132;} else {fl1|=NEE03;fl2|=NEE13;}} else
				if (va3==vb2) {if (va2==vb1) {fl1|=NFF023;fl2|=NFF132;} else if (va1==vb0) {fl1|=NFF031;fl2|=NFF023;} else {fl1|=NEE03;fl2|=NEE23;}}
			} else if (va1==vb1) {
				if (va2==vb2) {if (va3==vb0) {fl1|=NFF132;fl2|=NFF012;} else if (va0==vb3) {fl1|=NFF012;fl2|=NFF132;} else {fl1|=NEE12;fl2|=NEE12;}} else
				if (va2==vb3) {if (va3==vb2) {fl1|=NFF132;fl2|=NFF132;} else if (va0==vb0) {fl1|=NFF012;fl2|=NFF031;} else {fl1|=NEE12;fl2|=NEE13;}} else
				if (va3==vb3) {if (va2==vb0) {fl1|=NFF132;fl2|=NFF031;} else if (va0==vb2) {fl1|=NFF031;fl2|=NFF132;} else {fl1|=NEE13;fl2|=NEE13;}}
			} else if (va1==vb2) {
				if (va2==vb1) {if (va0==vb0) {fl1|=NFF012;fl2|=NFF012;} else if (va3==vb3) {fl1|=NFF132;fl2|=NFF132;} else {fl1|=NEE12;fl2|=NEE12;}} else
				if (va2==vb3) {if (va0==vb1) {fl1|=NFF012;fl2|=NFF132;} else if (va3==vb0) {fl1|=NFF132;fl2|=NFF023;} else {fl1|=NEE12;fl2|=NEE23;}} else
				if (va3==vb3) {if (va2==vb1) {fl1|=NFF132;fl2|=NFF132;} else if (va0==vb0) {fl1|=NFF031;fl2|=NFF023;} else {fl1|=NEE13;fl2|=NEE23;}}			
			} else if (va1==vb3) {
				if (va2==vb1) {if (va3==vb0) {fl1|=NFF132;fl2|=NFF031;} else if (va0==vb2) {fl1|=NFF012;fl2|=NFF132;} else {fl1|=NEE12;fl2|=NEE13;}} else
				if (va2==vb2) {if (va3==vb1) {fl1|=NFF132;fl2|=NFF132;} else if (va0==vb0) {fl1|=NFF012;fl2|=NFF023;} else {fl1|=NEE12;fl2|=NEE23;}} else
				if (va3==vb1) {if (va2==vb2) {fl1|=NFF132;fl2|=NFF132;} else if (va0==vb0) {fl1|=NFF031;fl2|=NFF031;} else {fl1|=NEE13;fl2|=NEE13;}} else
				if (va3==vb2) {if (va2==vb0) {fl1|=NFF132;fl2|=NFF023;} else if (va0==vb1) {fl1|=NFF031;fl2|=NFF132;} else {fl1|=NEE13;fl2|=NEE23;}}
			} else if (va2==vb2) {
				if (va3==vb3) {if (va0==vb1) {fl1|=NFF023;fl2|=NFF132;} else if (va1==vb0) {fl1|=NFF132;fl2|=NFF023;} else {fl1|=NEE23;fl2|=NEE23;}}
			} else if (va2==vb3) {
				if (va3==vb2) {if (va0==vb0) {fl1|=NFF023;fl2|=NFF023;} else if (va1==vb1) {fl1|=NFF132;fl2|=NFF132;} else {fl1|=NEE23;fl2|=NEE23;}}
			}

			// the rules need to be checked again in reverse aspect (because they are compactified for one-way comparison)
			if (i == 2) {	if (fl1 != 0) break;									// no need to do reverse check if neighbour found
							int tmp = va0; va0 = vb0; vb0 = tmp;
								tmp = va1; va1 = vb1; vb1 = tmp;
								tmp = va2; va2 = vb2; vb2 = tmp;
								tmp = va3; va3 = vb3; vb3 = tmp; }
			int tmp = fl1; fl1 = fl2; fl2 = tmp;
		}
		return fl1 | (fl2 << 16);
	}

	
	// method produces edge lengths of tetrahedron
	// calculation optimised for tetrahedral meshes where neighbour tetrahedrons can share their edge & tri.area calculations
	// zero edge length signifies an uncalculated tetrahedral edge (which means, set to zero to enforce recalculation)
	public void evaluateEdges() {
		boolean setHasDeltas = false;
		if (data == null) { data = new double[6 + 4]; setHasDeltas=true; }	// note: we know we have ALL deltas only if we're evaluating ALL edges
		if (delta == null) delta = new double[5 * FEM1.NCOORD];				// we need 3 (1-0,2-0,3-0) + 2 (2-1,3-1) x3 spaces for deltas
		int n2 = nodeRef[2] * FEM1.NCOORD, n3 = nodeRef[3] * FEM1.NCOORD;
		double[] node0, node1, node2, node3;
		if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
		if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;
		if (hasDeltas()) {
			if (data[TD_EDGE01] == 0) { data[TD_EDGE01] = Math.sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]); }
			if (data[TD_EDGE02] == 0) { data[TD_EDGE02] = Math.sqrt(delta[3]*delta[3] + delta[4]*delta[4] + delta[5]*delta[5]); }
			if (data[TD_EDGE03] == 0) { data[TD_EDGE03] = Math.sqrt(delta[6]*delta[6] + delta[7]*delta[7] + delta[8]*delta[8]); }
			if (data[TD_EDGE12] == 0) { data[TD_EDGE12] = Math.sqrt(delta[9]*delta[9] + delta[10]*delta[10] + delta[11]*delta[11]); }
			if (data[TD_EDGE13] == 0) { data[TD_EDGE13] = Math.sqrt(delta[12]*delta[12] + delta[13]*delta[13] + delta[14]*delta[14]); }
			if (data[TD_EDGE23] == 0) { double x32 = -node2[n2++]+node3[n3++], y32 = -node2[n2++]+node3[n3++], z32 = -node2[n2]+node3[n3];
			data[TD_EDGE23] = Math.sqrt(x32*x32 + y32*y32 + z32*z32); }
		} else {
			int n0 = nodeRef[0] * FEM1.NCOORD, n1 = nodeRef[1] * FEM1.NCOORD;
			if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
			if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
			if (data[TD_EDGE01] == 0) { double x10 = delta[0] = -node0[n0++]+node1[n1++], y10 = delta[1] = -node0[n0++]+node1[n1++], z10 = delta[2] = -node0[n0]+node1[n1];
										data[TD_EDGE01] = Math.sqrt(x10*x10 + y10*y10 + z10*z10); n0 -= 2; n1 -= 2; }
			if (data[TD_EDGE02] == 0) { double x20 = delta[3] = -node0[n0++]+node2[n2++], y20 = delta[4] = -node0[n0++]+node2[n2++], z20 = delta[5] = -node0[n0]+node2[n2];
										data[TD_EDGE02] = Math.sqrt(x20*x20 + y20*y20 + z20*z20); n0 -= 2; n2 -= 2; }
			if (data[TD_EDGE03] == 0) { double x30 = delta[6] = -node0[n0++]+node3[n3++], y30 = delta[7] = -node0[n0++]+node3[n3++], z30 = delta[8] = -node0[n0]+node3[n3];
										data[TD_EDGE03] = Math.sqrt(x30*x30 + y30*y30 + z30*z30); n0 -= 2; n3 -= 2; }
			if (data[TD_EDGE12] == 0) { double x21 = delta[9] = -node1[n1++]+node2[n2++], y21 = delta[10] = -node1[n1++]+node2[n2++], z21 = delta[11] = -node1[n1]+node2[n2];
										data[TD_EDGE12] = Math.sqrt(x21*x21 + y21*y21 + z21*z21); n1 -= 2; n2 -= 2; }
			if (data[TD_EDGE13] == 0) { double x31 = delta[12] = -node1[n1++]+node3[n3++], y31 = delta[13] = -node1[n1++]+node3[n3++], z31 = delta[14] = -node1[n1]+node3[n3];
										data[TD_EDGE13] = Math.sqrt(x31*x31 + y31*y31 + z31*z31); n1 -= 2; n3 -= 2; }
			if (data[TD_EDGE23] == 0) { double x32 = -node2[n2++]+node3[n3++], y32 = -node2[n2++]+node3[n3++], z32 = -node2[n2]+node3[n3];
										data[TD_EDGE23] = Math.sqrt(x32*x32 + y32*y32 + z32*z32); }
		}
		if (setHasDeltas) flagHasDeltas();
		flagHasEdges();
	}

	
	// calculate the radius of a sphere that inscribes/fits into a tetrahedron element
	static final double DIV6 = 1./6.;
	public void evaluateAreas(boolean findInscribedRadius) {

		if (data == null) data = new double[6 + 4];
		if (!findInscribedRadius && hasAreas()) return;						// do nothing more if all areas been calculated & that was the only request

		double x10, y10, z10, x20, y20, z20, x30, y30, z30, x21, y21, z21, x31, y31, z31, Vt6 = 0;
		if (hasDeltas()) {
			x10 = delta[0]; y10 = delta[1]; z10 = delta[2]; x20 = delta[3]; y20 = delta[4]; z20 = delta[5]; x30 = delta[6]; y30 = delta[7]; z30 = delta[8];
			x21 = delta[9]; y21 = delta[10]; z21 = delta[11]; x31 = delta[12]; y31 = delta[13]; z31 = delta[14];
		} else {
			double x0, y0, z0, x1, y1, z1;
			if (delta == null) delta = new double[5 * FEM1.NCOORD];
			int n0 = nodeRef[0] * FEM1.NCOORD, n1 = nodeRef[1] * FEM1.NCOORD, n2 = nodeRef[2] * FEM1.NCOORD, n3 = nodeRef[3] * FEM1.NCOORD; 
			double[] node0, node1, node2, node3;
			if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
			if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
			if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
			if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;
			x0 = node0[n0++]; y0 = node0[n0++]; z0 = node0[n0]; x1 = node1[n1++]; y1 = node1[n1++]; z1 = node1[n1];
			x10 = delta[0] = x1-x0; y10 = delta[1] = y1-y0; z10 = delta[2] = z1-z0;
			x20 = delta[3] = node2[n2++] - x0; y20 = delta[4] = node2[n2++] - y0; z20 = delta[5] = node2[n2--] - z0; n2--;
			x30 = delta[6] = node3[n3++] - x0; y30 = delta[7] = node3[n3++] - y0; z30 = delta[8] = node3[n3--] - z0; n3--;
			x21 = delta[9] = node2[n2++] - x1; y21 = delta[10] = node2[n2++] - y1; z21 = delta[11] = node2[n2] - z1;
			x31 = delta[12] = node3[n3++] - x1; y31 = delta[13] = node3[n3++] - y1; z31 = delta[14] = node3[n3] - z1;
		}
		
		// find tetrahedron's facet areas
		double sq1 = y10 * z20 - z10 * y20, sq2 = z10 * x20 - x10 * z20, sq3 = x10 * y20 - y10 * x20;
		if (data[TD_AREA012] == 0)
			data[6] = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);			// facet 012
		
		if (findInscribedRadius) {
			// Vt6 is tetrahedron volume * 6
			Vt6 = (sq1 * x30 + sq2 * y30 + sq3 * z30);
			if (Vt6 < 0) Vt6 = -Vt6;
			if (volume == 0) volume = Vt6 * DIV6;							// store the volume (if it's not calculated)
		}
		
		if (data[TD_AREA023] == 0) {
			sq1 = y10 * z30 - z10 * y30; sq2 = z10 * x30 - x10 * z30; sq3 = x10 * y30 - y10 * x30;
			data[7] = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);			// facet 031
		}
		if (data[TD_AREA031] == 0) {
			sq1 = y20 * z30 - z20 * y30; sq2 = z20 * x30 - x20 * z30; sq3 = x20 * y30 - y20 * x30;
			data[8] = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);			// facet 023
		}
		if (data[TD_AREA132] == 0) {
			sq1 = y21 * z31 - z21 * y31; sq2 = z21 * x31 - x21 * z31; sq3 = x21 * y31 - y21 * x31;
			data[9] = 0.5 * Math.sqrt(sq1*sq1 + sq2*sq2 + sq3*sq3);			// facet 132
		}
		flagHasAreas();
		// evaluate the inscribed sphere radius
		if (findInscribedRadius)
			iRadius = 0.5 * Vt6 / (data[TD_AREA012] + data[TD_AREA023] + data[TD_AREA031] + data[TD_AREA132]);
	}
	
	
	public void evaluateNormals(boolean normalise) {
		if (hasNormals()) return;
		if (!hasXNormals()) {
			setupData(NORMAL_FIELD);							// we want to save the calculated normals
			int n0, n1, n2, n3;
			
			n0 = nodeRef[0] * FEM1.NCOORD; n1 = nodeRef[1] * FEM1.NCOORD; n2 = nodeRef[2] * FEM1.NCOORD; n3 = nodeRef[3] * FEM1.NCOORD;
			double[] node0, node1, node2, node3;
			if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
			if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
			if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
			if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;

			double x0 = node0[n0++], y0 = node0[n0++], z0 = node0[n0], x1 = node1[n1++], y1 = node1[n1++], z1 = node1[n1];	
			double x10 = x1 - x0, y10 = y1 - y0, z10 = z1 - z0;
			double x2 = node2[n2++], y2 = node2[n2++], z2 = node2[n2];
			double x20 = x2 - x0, y20 = y2 - y0, z20 = z2 - z0, x21 = x2 - x1, y21 = y2 - y1, z21 = z2 - z1;
			if ((flags & FLAG_NORM012) == 0) {		// if facet 012 normal uninitialised
				data[TD_NORM012] = y10 * z20 - z10 * y20; data[TD_NORM012+1] = z10 * x20 - x10 * z20; data[TD_NORM012+2] = x10 * y20 - y10 * x20;
			}
			double x3 = node3[n3++], y3 = node3[n3++], z3 = node3[n3];
			double x30 = x3 - x0, y30 = y3 - y0, z30 = z3 - z0;
			if ((flags & FLAG_NORM023) == 0) {		// if facet 023 normal uninitialised
				data[TD_NORM023] = y20 * z30 - z20 * y30; data[TD_NORM023+1] = z20 * x30 - x20 * z30; data[TD_NORM023+2] = x20 * y30 - y20 * x30;
			}
			if ((flags & FLAG_NORM031) == 0) {		// if facet 031 normal uninitialised
				data[TD_NORM031] = y30 * z10 - z30 * y10; data[TD_NORM031+1] = x30 * z10 - z30 * x10; data[TD_NORM031+2] = x30 * y10 - y30 * x10;
			}
			if ((flags & FLAG_NORM132) == 0) {		// if facet 132 normal uninitialised
				double x31 = x3 - x1, y31 = y3 - y1, z31 = z3 - z1;
				data[TD_NORM132] = y31 * z21 - z31 * y21; data[TD_NORM132+1] = z31 * x21 - x31 * z21; data[TD_NORM132+2] = x31 * y21 - y31 * x21;
			}
			flags |= FLAG_NORM012 + FLAG_NORM023 + FLAG_NORM031 + FLAG_NORM132;
			flagHasXNormals();
		}
		if (normalise) {
			double l012D = 1. / Math.sqrt(data[TD_NORM012]*data[TD_NORM012]+data[TD_NORM012+1]*data[TD_NORM012+1]+data[TD_NORM012+2]*data[TD_NORM012+2]);
			data[TD_NORM012]*=l012D; data[TD_NORM012+1]*=l012D; data[TD_NORM012+2]*=l012D;
			double l023D = 1. / Math.sqrt(data[TD_NORM023]*data[TD_NORM023]+data[TD_NORM023+1]*data[TD_NORM023+1]+data[TD_NORM023+2]*data[TD_NORM023+2]);
			data[TD_NORM023]*=l023D; data[TD_NORM023+1]*=l023D; data[TD_NORM023+2]*=l023D;
			double l031D = 1. / Math.sqrt(data[TD_NORM031]*data[TD_NORM031]+data[TD_NORM031+1]*data[TD_NORM031+1]+data[TD_NORM031+2]*data[TD_NORM031+2]);
			data[TD_NORM031]*=l031D; data[TD_NORM031+1]*=l031D; data[TD_NORM031+2]*=l031D;
			double l132D = 1. / Math.sqrt(data[TD_NORM132]*data[TD_NORM132]+data[TD_NORM132+1]*data[TD_NORM132+1]+data[TD_NORM132+2]*data[TD_NORM132+2]);
			data[TD_NORM132]*=l132D; data[TD_NORM132+1]*=l132D; data[TD_NORM132+2]*=l132D;
			clearHasXNormals();
			flagHasNormals();
		}
	}
	
	
	// method evaluares the radius of a sphere that is enclosed by tetrahedron
	public double tetraInscribedRadius() {
		if (iRadius > 0) return iRadius;						// a zero inscribed radius indicates it's not calculated/bound for recalculation
		evaluateAreas(true); return iRadius;
	}
	
	
	// method calculates the radius of a sphere that circumscribes/envelopes a tetrahedron element
	public double tetraCircumRadius() {

		if (cRadius > 0) return cRadius;						// a zero circumradius indicates it's not calculated/bound for recalculation
		if (!hasEdges()) evaluateEdges();
		double e01e23 = data[TD_EDGE01] * data[TD_EDGE23], e02e13 = data[TD_EDGE02] * data[TD_EDGE13], e03e12 = data[TD_EDGE03] * data[TD_EDGE12];
		double M5 = (e01e23 + e02e13 + e03e12) * (-e01e23 + e02e13 + e03e12) * (e01e23 - e02e13 + e03e12) * (e01e23 + e02e13 - e03e12) * 0.0625;
		double volume6;
		if (volume <= 0) volume = (volume6 = tetraVolume(true)) / 6;
		else volume6 = volume * 6;
		return cRadius = Math.sqrt(M5 < 0 ? -M5 : M5) / volume6;
	}
	
	
	// method determines the barycentric or circumcenter coordinates of tetrahedron by setting up an Ax = b type 3x3-matrix solution,
	// getting det(A) and det(Ai) according to Cramer's rule, thus solving the three coordinates, plus resolving the circumcenter
	// coordinates if requested, returns null if det(A) is near zero
	public double[] tetraBaryCircum(boolean circumcenter) {

		double x0, y0, z0, x10, y10, z10, x20, y20, z20, x30, y30, z30;
		int n0 = nodeRef[0] * FEM1.NCOORD;
		double[] node0;
		if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
		if (hasDeltas()) {
			x10 = delta[0]; y10 = delta[1]; z10 = delta[2]; x20 = delta[3]; y20 = delta[4]; z20 = delta[5]; x30 = delta[6]; y30 = delta[7]; z30 = delta[8];
		} else {
			int n1 = nodeRef[1] * FEM1.NCOORD, n2 = nodeRef[2] * FEM1.NCOORD, n3 = nodeRef[3] * FEM1.NCOORD; 
			double[] node1, node2, node3;
			if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
			if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
			if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;
			x0 = node0[n0++]; y0 = node0[n0++]; z0 = node0[n0--]; n0--;
			x10 = node1[n1++]-x0; y10 = node1[n1++]-y0; z10 = node1[n1]-z0;
			x20 = node2[n2++]-x0; y20 = node2[n2++]-y0; z20 = node2[n2]-z0;
			x30 = node3[n3++]-x0; y30 = node3[n3++]-y0; z30 = node3[n3]-z0;
		}

		double[] m3x3 = new double[9];
		m3x3[0] = x10*x10 + y10*y10 + z10*z10;	m3x3[1] = x10*x20 + y10*y20 + z10*z20;	m3x3[2] = x10*x30 + y10*y30 + z10*z30;
		m3x3[3] = m3x3[1];						m3x3[4] = x20*x20 + y20*y20 + z20*z20;	m3x3[5] = x20*x30 + y20*y30 + z20*z30;
		m3x3[6] = m3x3[2];						m3x3[7] = m3x3[5];						m3x3[8] = x30*x30 + y30*y30 + z30*z30;
		double det_m3x3I = Matrix.determinant3x3(m3x3);
		if (Matrix.nearZero(det_m3x3I)) return null; else det_m3x3I = 1. / det_m3x3I;
		double bX = .5 * m3x3[0], bY = .5 * m3x3[4], bZ = .5 * m3x3[8];
		double[] m3x3X = m3x3.clone(), m3x3Y = m3x3.clone(), m3x3Z = m3x3.clone();
		m3x3X[0] = bX; m3x3X[3] = bY; m3x3X[6] = bZ; m3x3Y[1] = bX; m3x3Y[4] = bY; m3x3Y[7] = bZ; m3x3Z[2] = bX; m3x3Z[5] = bY; m3x3Z[8] = bZ;
		
		if (circumcenter) {
			double[] cCenter = new double[3];
			double b0 = Matrix.determinant3x3(m3x3X)*det_m3x3I, b1 = Matrix.determinant3x3(m3x3Y)*det_m3x3I, b2 = Matrix.determinant3x3(m3x3Z)*det_m3x3I;
			cCenter[0] = node0[n0++] + b0 * x10 + b1 * x20 + b2 * x30;
			cCenter[1] = node0[n0++] + b0 * y10 + b1 * y20 + b2 * y30;
			cCenter[2] = node0[n0]   + b0 * z10 + b1 * z20 + b2 * z30;
			return cCenter;
		}
		
		double[] cBarycentric = new double[4];
		cBarycentric[0] = Matrix.determinant3x3(m3x3X) * det_m3x3I;
		cBarycentric[1] = Matrix.determinant3x3(m3x3Y) * det_m3x3I;
		cBarycentric[2] = Matrix.determinant3x3(m3x3Z) * det_m3x3I;
		cBarycentric[3] = 1.0 - cBarycentric[0] - cBarycentric[1] - cBarycentric[2];
		return cBarycentric;
	}
	
	public double[] tetraCircumcenter() { double[] cCenter = tetraBaryCircum(true); return cCenter; }
	public double[] tetraBarycentric() { double[] cBarycenter = tetraBaryCircum(false); return cBarycenter; }

	
	
	// finds volume of tetrahedron, optionally returns it directly as V*6 since some calculations need it
	public double tetraVolume(boolean times6) {
		if (volume > 0) return volume;							// a zero volume indicates indicates it's not calculated/bound for recalculation		
		double x0, y0, z0, x10, y10, z10, x20, y20, z20, x30, y30, z30;
		if (hasDeltas()) {
			x10 = delta[0]; y10 = delta[1]; z10 = delta[2]; x20 = delta[3]; y20 = delta[4]; z20 = delta[5]; x30 = delta[6]; y30 = delta[7]; z30 = delta[8];
		} else {
			int n0 = nodeRef[0] * FEM1.NCOORD, n1 = nodeRef[1] * FEM1.NCOORD, n2 = nodeRef[2] * FEM1.NCOORD, n3 = nodeRef[3] * FEM1.NCOORD; 
			double[] node0, node1, node2, node3;
			if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
			if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
			if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
			if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;
			x0 = node0[n0++]; y0 = node0[n0++]; z0 = node0[n0]; x10=node1[n1++]-x0; y10=node1[n1++]-y0; z10=node1[n1]-z0;
			x20=node2[n2++]-x0; y20=node2[n2++]-y0; z20=node2[n2]-z0; x30=node3[n3++]-x0; y30=node3[n3++]-y0; z30=node3[n3]-z0;
		}
		// find tetrahedron's triangular areas
		double sq1 = y10 * z20 - z10 * y20, sq2 = z10 * x20 - x10 * z20, sq3 = x10 * y20 - y10 * x20;
		double Vt6 = (sq1 * x30 + sq2 * y30 + sq3 * z30);	// Vt6 is tetrahedron volume * 6		
		if (times6) return (Vt6 < 0 ? -Vt6 : Vt6);				// caller wants V*6 returned, do not store
		return volume = (Vt6 < 0 ? -Vt6 : Vt6) * DIV6;			// store the volume
	}
	
	
	// finds center od mass of tetrahedron
	static final double DIV4 = 1./4.;
	public double[] tetraMassCenter() {
		double[] massC = {0, 0, 0};
		int n0 = nodeRef[0] * FEM1.NCOORD, n1 = nodeRef[1] * FEM1.NCOORD, n2 = nodeRef[2] * FEM1.NCOORD, n3 = nodeRef[3] * FEM1.NCOORD;
		double[] node0, node1, node2, node3;
		if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
		if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
		if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
		if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;
		massC[0] = (node0[n0++] + node1[n1++] + node2[n2++] + node3[n3++]) * DIV4;
		massC[1] = (node0[n0++] + node1[n1++] + node2[n2++] + node3[n3++]) * DIV4;
		massC[2] = (node0[n0]   + node1[n1]   + node2[n2]   + node3[n3])   * DIV4;
		return massC;
	}
	
	// returns quality of this element if it's a tetrahedron: the circumradius to longest edge ratio
	// the constant is a correction factor alpha = 2 * sqrt(6)
	public double tetraQuality() {
		if (!hasEdges()) evaluateEdges();
		double elMaxA = data[TD_EDGE01] > data[TD_EDGE02] ? data[TD_EDGE01] : data[TD_EDGE02];		// compare out longest edge out of 6
		double elMaxB = data[TD_EDGE03] > data[TD_EDGE12] ? data[TD_EDGE03] : data[TD_EDGE12];
		double elMaxC = data[TD_EDGE13] > data[TD_EDGE23] ? data[TD_EDGE12] : data[TD_EDGE23];
		double elMaxAB = elMaxA > elMaxB ? elMaxA : elMaxB, elMaxABC = elMaxAB > elMaxC ? elMaxAB : elMaxC;
		return 4.898979485566356 * tetraInscribedRadius() / elMaxABC;
	}
	
	// returns quality based on circumscribed radius to shortest edge ratio
	public double tetraQuality2() {
		if (!hasEdges()) evaluateEdges();
		double elMinA = data[TD_EDGE01] < data[TD_EDGE02] ? data[TD_EDGE01] : data[TD_EDGE02];		// compare out shortest edge out of 6
		double elMinB = data[TD_EDGE03] < data[TD_EDGE12] ? data[TD_EDGE03] : data[TD_EDGE12];
		double elMinC = data[TD_EDGE13] < data[TD_EDGE23] ? data[TD_EDGE13] : data[TD_EDGE23];
		double elMinAB = elMinA < elMinB ? elMinA : elMinB, elMinABC = elMinAB < elMinC ? elMinAB : elMinC;
		return tetraCircumRadius() / elMinABC;
	}
	
		
	// returns the gradient of tetrahedral volume with respect to each vertex or one or more chosen vertices
	public void tetraVolumeGradient(boolean g0, boolean g1, boolean g2, boolean g3) {

		int g = 10;
		setupData(VOLGRADIENT_FIELD);
		flagHasVGradients();
	
		int n0 = nodeRef[0] * FEM1.NCOORD, n1 = nodeRef[1] * FEM1.NCOORD, n2 = nodeRef[2] * FEM1.NCOORD, n3 = nodeRef[3] * FEM1.NCOORD;
		double[] node0, node1, node2, node3;
		if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
		if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
		if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
		if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;

		if (g0 && g1 && g2 && g3) {
			double x1 = node1[n1++], y1 = node1[n1++], z1 = node1[n1];
			double x31 = node3[n3++] - x1, y31 = node3[n3++] - y1, z31 = node3[n3] - z1; n3 -= 2;
			double x2 = node2[n2++], y2 = node2[n2++], z2 = node2[n2];
			double x21 = x2 - x1, y21 = y2 - y1, z21 = z2 - z1;
			data[g++] = y31*z21 - z31*y21; data[g++] = x31*z21 - z31*x21; data[g++] = x31*y21 - y31*x21;
			double x0 = node0[n0++], y0 = node0[n0++], z0 = node0[n0];
			double x20 = x2 - x0, y20 = y2 - y0, z20 = z2 - z0;
			double x30 = node3[n3++] - x0, y30 = node3[n3++] - y0, z30 = node3[n3] - z0;
			data[g++] = y20*z30 - z20*y30; data[g++] = x20*z30 - z20*x30; data[g++] = x20*y30 - y20*x30;
			double x10 = x1 - x0, y10 = y1 - y0, z10 = z1 - z0;
			data[g++] = y30*z10 - z30*y10; data[g++] = x30*z10 - z30*x10; data[g++] = x30*y10 - y30*x10;
			data[g++] = y10*z20 - z10*y20; data[g++] = x10*z20 - z10*x20; data[g++] = x10*y20 - y10*x20;
			return;
		}
		if (g0) {
			double x1 = node1[n1++], y1 = node1[n1++], z1 = node1[n1];
			double x31 = node3[n3++] - x1, y31 = node3[n3++] - y1, z31 = node3[n3] - z1; n3 -= 2;
			double x2 = node2[n2++], y2 = node2[n2++], z2 = node2[n2];
			double x21 = x2 - x1, y21 = y2 - y1, z21 = z2 - z1;
			data[g++] = y31*z21 - z31*y21; data[g++] = x31*z21 - z31*x21; data[g++] = x31*y21 - y31*x21; }
		if (g1) {
			double ax = node0[n0++], ay = node0[n0++], az = node0[n0];
			double cax = node2[n2++] - ax, cay = node2[n2++] - ay, caz = node2[n2] - az;
			double dax = node3[n3++] - ax, day = node3[n3++] - ay, daz = node3[n3] - az;
			data[g++] = cay*daz - caz*day; data[g++] = cax*daz - caz*dax; data[g++] = cax*day - cay*dax; }
		if (g2) {
			double ax = node0[n0++], ay = node0[n0++], az = node0[n0];
			double bax = node1[n1++] - ax, bay = node1[n1++] - ay, baz = node1[n1] - az;
			double dax = node3[n3++] - ax, day = node3[n3++] - ay, daz = node3[n3] - az;
			data[g++] = day*baz - daz*bay; data[g++] = dax*baz - daz*bax; data[g++] = dax*bay - day*bax; }
		if (g3) {
			double ax = node0[n0++], ay = node0[n0++], az = node0[n0];
			double bax = node1[n1++] - ax, bay = node1[n1++] - ay, baz = node1[n1] - az;
			double cax = node2[n2++] - ax, cay = node2[n2++] - ay, caz = node2[n2] - az;
			data[g++] = bay*caz - baz*cay; data[g++] = bax*caz - baz*cax; data[g++] = bax*cay - bay*cax;
		}
	}
	
	
	// method returns 0 if supplied vertex is inside tetrahedron, and >0 if outside
	// if rigorousTest=true, method returns 1 in bit 0 if supplied vertex is inside tetrahedron, and bits NFF012<<1 or NFF023<<1 or NFF031<<1 or NFF132<<1
	// or a mix of those if vertex is outside the plane/planes defined by facets NFF012 or NFF023 or NFF031 or NFF132
	public int tetraVertexEnclosure(double x, double y, double z, boolean rigorousTest) {

		double[] node0, node1, node2, node3;
		int n0, n1, n2, n3, enclosure = 0;
		double n012x, n012y, n012z, n023x, n023y, n023z, n031x, n031y, n031z, n132x, n132y, n132z;	
		n0 = nodeRef[0] * FEM1.NCOORD; n1 = nodeRef[1] * FEM1.NCOORD; n2 = nodeRef[2] * FEM1.NCOORD; n3 = nodeRef[3] * FEM1.NCOORD;
		if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
		if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
		if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
		if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;
		double x0 = node0[n0++], y0 = node0[n0++], z0 = node0[n0], x1 = node1[n1++], y1 = node1[n1++], z1 = node1[n1];
		
		if (hasXNormals() || hasNormals()) {				// if we already have calculated normals or crossproduct normals
			double nx = x - x0, ny = y - y0, nz = z - z0;	// do immediate dot product comparisons from node0 and node1 aspects
			if (data[22] * nx + data[23] * ny + data[24] * nz > 0) {		enclosure = (NFF012<<1)+1; if (!rigorousTest) return enclosure; }
			if (data[25] * nx + data[26] * ny + data[27] * nz > 0) {		enclosure |=(NFF023<<1)+1; if (!rigorousTest) return enclosure; }
			if (data[28] * nx + data[29] * ny + data[30] * nz > 0) {		enclosure |=(NFF031<<1)+1; if (!rigorousTest) return enclosure; }
			if (data[31]*(x-x1) + data[32]*(y-y1) + data[33]*(z-z1) > 0) {	enclosure |=(NFF132<<1)+1; if (!rigorousTest) return enclosure; }
			return enclosure;
		}

		setupData(NORMAL_FIELD);							// we want to save calculated crossproduct normals
		double x10 = x1 - x0, y10 = y1 - y0, z10 = z1 - z0;
		double x2 = node2[n2++], y2 = node2[n2++], z2 = node2[n2];
		double x20 = x2 - x0, y20 = y2 - y0, z20 = z2 - z0, x21 = x2 - x1, y21 = y2 - y1, z21 = z2 - z1;
		n012x = data[TD_NORM012] = y10 * z20 - z10 * y20; n012y = data[TD_NORM012+1] = z10 * x20 - x10 * z20; n012z = data[TD_NORM012+2] = x10 * y20 - y10 * x20;
		flags |= FLAG_NORM012;
		double nx = x - x0, ny = y - y0, nz = z - z0;
		if (n012x * nx + n012y * ny + n012z * nz > 0) {	
			enclosure = (NFF012<<1)+1; if (!rigorousTest) return enclosure; }
		
		double x3 = node3[n3++], y3 = node3[n3++], z3 = node3[n3];
		double x30 = x3 - x0, y30 = y3 - y0, z30 = z3 - z0;
		n023x = data[TD_NORM023] = y20 * z30 - z20 * y30; n023y = data[TD_NORM023+1] = z20 * x30 - x20 * z30; n023z = data[TD_NORM023+2] = x20 * y30 - y20 * x30;
		flags |= FLAG_NORM023;
		if (n023x * nx + n023y * ny + n023z * nz > 0) {	
			enclosure |=(NFF023<<1)+1; if (!rigorousTest) return enclosure; }
		
		n031x = data[TD_NORM031] = y30 * z10 - z30 * y10; n031y = data[TD_NORM031+1] = x30 * z10 - z30 * x10; n031z = data[TD_NORM031+2] = x30 * y10 - y30 * x10;
		flags |= FLAG_NORM031;
		if (n031x * nx + n031y * ny + n031z * nz > 0) {
			enclosure |=(NFF031<<1)+1; if (!rigorousTest) return enclosure; }
		
		double x31 = x3 - x1, y31 = y3 - y1, z31 = z3 - z1;
		n132x = data[TD_NORM132] = y31 * z21 - z31 * y21; n132y = data[TD_NORM132+1] = z31 * x21 - x31 * z21; n132z = data[TD_NORM132+2] = x31 * y21 - y31 * x21;
		flags |= FLAG_NORM132;
		if (n132x * (x - x1) + n132y * (y - y1) + n132z * (z - z1) > 0) {
			enclosure |=(NFF132<<1)+1; if (!rigorousTest) return enclosure; }
		
		flagHasNormals();
		return enclosure;
	}
	
	
	// method checks if a line segment intersects the tetrahedron, if is assumed that the end coordinates already have
	// been checked with tetraVertexEnclosure(), checkEnds=false will skip testing whether segment ends are inside the tetrahedron
	// if isect[] supplied, the intersection coordinates will be returned inside, an (x,y,z)-triplet for every (affected, otherwise (0,0,0)) facet
	public int tetraSegmentIntersection(double xa, double ya, double za, double xb, double yb, double zb, boolean checkEnds, double[] isect) {
		
		int intersector = 0;
		if (checkEnds) {
			intersector =  tetraVertexEnclosure(xb, yb, zb, true)<<5; 	//if (intersector > 0) return a_enclosure;	// DEBUG: do complete test
			intersector |= tetraVertexEnclosure(xa, ya, za, true); 		//if (intersector > 0) return b_enclosure;	// DEBUG: do complete test
		}
		
		double[] node0, node1, node2, node3;
		int n0, n1, n2, n3;
		n0 = nodeRef[0] * FEM1.NCOORD; n1 = nodeRef[1] * FEM1.NCOORD; n2 = nodeRef[2] * FEM1.NCOORD; n3 = nodeRef[3] * FEM1.NCOORD;
		if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
		if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
		if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
		if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;
		
		double x0 = node0[n0++], y0 = node0[n0++], z0 = node0[n0];
		if (!hasNormals()) evaluateNormals(true);
		
		// test the four facet cases: test if segment is parallel to facet plane, then test if it punctures the plane, then test if it punctures the facet
		double x0a = x0 - xa, y0a = y0 - ya, z0a = z0 - za, xba = xb-xa, yba = yb-ya, zba = zb-za;
		double rI012D = data[TD_NORM012]*xba + data[TD_NORM012+1]*yba + data[TD_NORM012+2]*zba;
		if (rI012D != 0) {
			double rI012 =	(data[TD_NORM012]*x0a + data[TD_NORM012+1]*y0a + data[TD_NORM012+2]*z0a) / rI012D;
			if (0 <= rI012 && rI012 <= 1) {				// if line segment intersects plane of facet 012
				double xv012 = xa + rI012*xba - x0, yv012 = ya + rI012*yba - y0, zv012 = za + rI012*zba - z0;	// plane isect point - node0 vector
				double x10 = node1[n1++]-x0, y10 = node1[n1++]-y0, z10 = node1[n1]-z0, x20 = node2[n2++]-x0, y20 = node2[n2++]-y0, z20 = node2[n2]-z0;
				n1 -= 2; n2 -= 2;
				double v10v20 =	x10*x20 + y10*y20 + z10*z20;
				double v012v20 = xv012*x20 + yv012*y20 + zv012*z20, v20v20 = x20*x20 + y20*y20 + z20*z20, v012v10 = xv012*x10 + yv012*y10 + zv012*z10;
				double v10v10 = x10*x10 + y10*y10 + z10*z10, st012D = 1. / (v10v20*v10v20 - v10v10*v20v20);
				double sI012 = (v10v20 * v012v20 - v20v20 * v012v10) * st012D, tI012 = (v10v20 * v012v10 - v10v10 * v012v20) * st012D;
				if (sI012 >= 0 && tI012 >= 0 && sI012 + tI012 <= 1) {	// if intersection point is inside parameterised triangle
					if (isect != null) { isect[0] = x0+sI012*x10+tI012*x20; isect[1] = y0+sI012*y10+tI012*y20; isect[2] = z0+sI012*z10+tI012*z20; }
					intersector |= 1024;								// bit 10 flags intersection of facet012
				}
		}}
		double rI023D = data[TD_NORM023]*xba + data[TD_NORM023+1]*yba + data[TD_NORM023+2]*zba;
		if (rI023D != 0) {
			double rI023 =	(data[TD_NORM023]*x0a + data[TD_NORM023+1]*y0a + data[TD_NORM023+2]*z0a) / rI023D;
			if (0 <= rI023 && rI023 <= 1) {				// if line segment intersects plane of facet 023
				double xv023 = xa + rI023*xba - x0, yv023 = ya + rI023*yba - y0, zv023 = za + rI023*zba - z0;	// plane isect point - node0 vector
				double x20 = node2[n2++]-x0, y20 = node2[n2++]-y0, z20 = node2[n2]-z0, x30 = node3[n3++]-x0, y30 = node3[n3++]-y0, z30 = node3[n3]-z0;
				n2 -= 2; n3 -= 2;
				double v20v30 =	x20*x30 + y20*y30 + z20*z30;
				double v023v30 = xv023*x30 + yv023*y30 + zv023*z30, v30v30 = x30*x30 + y30*y30 + z30*z30, v023v20 = xv023*x20 + yv023*y20 + zv023*z20;
				double v20v20 = x20*x20 + y20*y20 + z20*z20, st023D = 1. / (v20v30*v20v30 - v20v20*v30v30);
				double sI023 = (v20v30 * v023v30 - v30v30 * v023v20) * st023D, tI023 = (v20v30 * v023v20 - v20v20 * v023v30) * st023D;
				if (sI023 >= 0 && tI023 >= 0 && sI023 + tI023 <= 1) {	// if intersection point is inside parameterised triangle
					if (isect != null) { isect[3] = x0+sI023*x20+tI023*x30; isect[4] = y0+sI023*y20+tI023*y30; isect[5] = z0+sI023*z20+tI023*z30; }
					intersector |= 2048;								// bit 11 flags intersection of facet023
				}
		}}
		double rI031D = data[TD_NORM031]*xba + data[TD_NORM031+1]*yba + data[TD_NORM031+2]*zba;
		if (rI031D != 0) {
			double rI031 =	(data[TD_NORM031]*x0a + data[TD_NORM031+1]*y0a + data[TD_NORM031+2]*z0a) / rI031D;
			if (0 <= rI031 && rI031 <= 1) {				// if line segment intersects plane of facet 031
				double xv031 = xa + rI031*xba - x0, yv031 = ya + rI031*yba - y0, zv031 = za + rI031*zba - z0;	// plane isect point - node0 vector
				double x10 = node1[n1++]-x0, y10 = node1[n1++]-y0, z10 = node1[n1]-z0, x30 = node3[n3++]-x0, y30 = node3[n3++]-y0, z30 = node3[n3]-z0;
				n1 -= 2; n3 -= 2;
				double v30v10 =	x30*x10 + y30*y10 + z30*z10;
				double v031v10 = xv031*x10 + yv031*y10 + zv031*z10, v30v30 = x30*x30 + y30*y30 + z30*z30, v031v30 = xv031*x30 + yv031*y30 + zv031*z30;
				double v10v10 = x10*x10 + y10*y10 + z10*z10, st031D = 1. / (v30v10*v30v10 - v30v30*v10v10);
				double sI031 = (v30v10 * v031v10 - v10v10 * v031v30) * st031D, tI031 = (v30v10 * v031v30 - v30v30 * v031v10) * st031D;
				if (sI031 >= 0 && tI031 >= 0 && sI031 + tI031 <= 1) {	// if intersection point is inside parameterised triangle
					if (isect != null) { isect[6] = x0+sI031*x30+tI031*x10; isect[7] = y0+sI031*y30+tI031*y10; isect[8] = z0+sI031*z30+tI031*z10; }
					intersector |= 4096;								// bit 12 flags intersection of facet031
				}
		}}
		double x1 = node1[n1++], y1 = node1[n1++], z1 = node1[n1];
		double rI132D = data[TD_NORM132]*xba + data[TD_NORM132+1]*yba + data[TD_NORM132+2]*zba;
		if (rI132D != 0) {
			double rI132 =	(data[TD_NORM132]*(x1 - xa) + data[TD_NORM132+1]*(y1 - ya) + data[TD_NORM132+2]*(z1 - za)) / rI132D;
			if (0 <= rI132 && rI132 <= 1) {				// if line segment intersects plane of facet 132
				double xv132 = xa + rI132*xba - x1, yv132 = ya + rI132*yba - y1, zv132 = za + rI132*zba - z1;	// plane isect point - node1 vector
				double x31 = node3[n3++]-x1, y31 = node3[n3++]-y1, z31 = node3[n3]-z1, x21 = node2[n2++]-x1, y21 = node2[n2++]-y1, z21 = node2[n2]-z1;
				double v31v21 =	x31*x21 + y31*y21 + z31*z21;
				double v132v21 = xv132*x21 + yv132*y21 + zv132*z21, v21v21 = x21*x21 + y21*y21 + z21*z21, v132v31 = xv132*x31 + yv132*y31 + zv132*z31;
				double v31v31 = x31*x31 + y31*y31 + z31*z31, st132D = 1. / (v31v21*v31v21 - v31v31*v21v21);
				double sI132 = (v31v21 * v132v21 - v21v21 * v132v31) * st132D, tI132 = (v31v21 * v132v31 - v31v31 * v132v21) * st132D;
				if (sI132 >= 0 && tI132 >= 0 && sI132 + tI132 <= 1) {	// if intersection point is inside parameterised triangle
					if (isect != null) { isect[9] = x1+sI132*x31+tI132*x21; isect[10] = y1+sI132*y31+tI132*y21; isect[11] = z1+sI132*z31+tI132*z21; }
					intersector |= 8192;								// bit 13 flags intersection of facet132
				}
		}}
		return intersector;
	}
	
	
	
	// method calculates circumcenter of one of 4 chosen facets of a tetrahedron: 012, 023, 031, 132
	// TODO: debug routine
	public double[] tetraFacetCircumCenter(int facet) {
		double xa=0, ya=0, za=0, xb=0, yb=0, zb=0, xc=0, yc=0, zc=0;
		int n0 = nodeRef[0] * FEM1.NCOORD, n1 = nodeRef[1] * FEM1.NCOORD, n2 = nodeRef[2] * FEM1.NCOORD, n3 = nodeRef[3] * FEM1.NCOORD;
		double[] node0, node1, node2, node3;
		if (n0 >= node.length) { node0 = nodeWork; n0 -= node.length; } else node0 = node;
		if (n1 >= node.length) { node1 = nodeWork; n1 -= node.length; } else node1 = node;
		if (n2 >= node.length) { node2 = nodeWork; n2 -= node.length; } else node2 = node;
		if (n3 >= node.length) { node3 = nodeWork; n3 -= node.length; } else node3 = node;
		switch (facet) {
		case NFF012: xa=node0[n0++]; ya=node0[n0++]; za=node0[n0]; xb=node1[n1++]; yb=node1[n1++]; zb=node1[n1]; xc=node2[n2++]; yc=node2[n2++]; zc=node2[n2]; break;
		case NFF023: xa=node0[n0++]; ya=node0[n0++]; za=node0[n0]; xb=node2[n2++]; yb=node2[n2++]; zb=node2[n2]; xc=node3[n3++]; yc=node3[n3++]; zc=node3[n3]; break;
		case NFF031: xa=node0[n0++]; ya=node0[n0++]; za=node0[n0]; xb=node3[n3++]; yb=node3[n3++]; zb=node3[n3]; xc=node1[n1++]; yc=node1[n1++]; zc=node1[n1]; break;
		case NFF132: xa=node1[n1++]; ya=node1[n1++]; za=node1[n1]; xb=node3[n3++]; yb=node3[n3++]; zb=node3[n3]; xc=node2[n2++]; yc=node2[n2++]; zc=node2[n2]; break;
		}
		double xba = xb - xa, yba = yb - ya, zba = zb - za, xca = xc - xa, yca = yc - ya, zca = zc - za;
		double xCr1 = yba * zca - zba * yca, yCr1 = xba * zca - zba * xca, zCr1 = xba * yca - yba * xca;
		double baXcaL2D = 1 / (Math.sqrt(xCr1*xCr1 + yCr1*yCr1 + zCr1*zCr1) * 2);
		double xCr3 = yca * zCr1 - zca * yCr1, yCr3 = xca * zCr1 - zca * xCr1, zCr3 = xca * yCr1 - yca * xCr1;
		double xCr2 = yCr1 * zba - zCr1 * yba, yCr2 = xCr1 * zba - zCr1 * xba, zCr2 = xCr1 * yba - yCr1 * xba;
		double caL = Math.sqrt(xca*xca + yca*yca + zca*zca), baL = Math.sqrt(xba*xba + yba*yba + zba*zba);
		double xSumD = caL * xCr2 + baL * xCr3, ySumD = caL * yCr2 + baL * yCr3, zSumD = caL * zCr2 + baL * zCr3;
		double[] cCenter = new double[3];
		cCenter[0] = xa + xSumD * baXcaL2D; cCenter[1] = ya + ySumD * baXcaL2D; cCenter[2] = za + zSumD * baXcaL2D;		
		return cCenter;
	}
		
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			OUTPUT METHODS
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Element: ");
		switch (flags & 15) {
			case TETRAHEDRON: sb.append("tetrahedron\n"); break;
			case HEXAHEDRON: sb.append("hexahedron\n"); break;
			case ISOTETRAHEDRON: sb.append("isotetrahedron\n"); break;
			case ISOHEXAHEDRON: sb.append("isohexahedron\n"); break;
		}
		if (volume > 0) sb.append(String.format("Volume: %.5f m^3\n", volume));
		switch (flags & 15) {
			case TETRAHEDRON:
				for (int n = 0; n < 4; n++) {
					sb.append("Node ref " + n + ": " + "n" + nodeRef[n]);
					int n1 = nodeRef[n]*3;
					double[] node0;
					if (n1 >= node.length) { n1 -= node.length; node0 = nodeWork; } else node0 = node;
					sb.append(", [" + node0[n1++] + "," + node0[n1++] + "," + node0[n1] + "]\n");
				}
				if (data != null) {
					if (data[TD_AREA012] > 0) sb.append("AREA 012: " + String.format("%.5f", data[TD_AREA012]) + " m^2\n");
					if (data[TD_AREA023] > 0) sb.append("AREA 023: " + String.format("%.5f", data[TD_AREA023]) + " m^2\n");
					if (data[TD_AREA031] > 0) sb.append("AREA 031: " + String.format("%.5f", data[TD_AREA031]) + " m^2\n");
					if (data[TD_AREA132] > 0) sb.append("AREA 132: " + String.format("%.5f", data[TD_AREA132]) + " m^2\n");
					if (data[TD_EDGE01] > 0) sb.append("edge 01: " + String.format("%.5f", data[TD_EDGE01]) + " m\n");
					if (data[TD_EDGE02] > 0) sb.append("edge 02: " + String.format("%.5f", data[TD_EDGE02]) + " m\n");
					if (data[TD_EDGE03] > 0) sb.append("edge 03: " + String.format("%.5f", data[TD_EDGE03]) + " m\n");
					if (data[TD_EDGE12] > 0) sb.append("edge 12: " + String.format("%.5f", data[TD_EDGE12]) + " m\n");
					if (data[TD_EDGE13] > 0) sb.append("edge 13: " + String.format("%.5f", data[TD_EDGE13]) + " m\n");
					if (data[TD_EDGE23] > 0) sb.append("edge 23: " + String.format("%.5f", data[TD_EDGE23]) + " m\n");
					if (data.length > 10) {
						if (data[TD_VGRAD0] != 0) sb.append(String.format("Vol.gradient n0: %.5f, %.5f, %.5f\n", data[TD_VGRAD0], data[TD_VGRAD0+1], data[TD_VGRAD0+2]));
						if (data[TD_VGRAD1] != 0) sb.append(String.format("Vol.gradient n1: %.5f, %.5f, %.5f\n", data[TD_VGRAD1], data[TD_VGRAD1+1], data[TD_VGRAD1+2]));
						if (data[TD_VGRAD2] != 0) sb.append(String.format("Vol.gradient n2: %.5f, %.5f, %.5f\n", data[TD_VGRAD2], data[TD_VGRAD2+1], data[TD_VGRAD2+2]));
						if (data[TD_VGRAD3] != 0) sb.append(String.format("Vol.gradient n3: %.5f, %.5f, %.5f\n", data[TD_VGRAD3], data[TD_VGRAD3+1], data[TD_VGRAD3+2]));
						if ((flags&FLAG_NORM012) != 0) sb.append(String.format("Normal f012: %.5f, %.5f, %.5f\n", data[TD_NORM012], data[TD_NORM012+1], data[TD_NORM012+2]));
						if ((flags&FLAG_NORM023) != 0) sb.append(String.format("Normal f023: %.5f, %.5f, %.5f\n", data[TD_NORM023], data[TD_NORM023+1], data[TD_NORM023+2]));
						if ((flags&FLAG_NORM031) != 0) sb.append(String.format("Normal f032: %.5f, %.5f, %.5f\n", data[TD_NORM031], data[TD_NORM031+1], data[TD_NORM031+2]));
						if ((flags&FLAG_NORM132) != 0) sb.append(String.format("Normal f132: %.5f, %.5f, %.5f\n", data[TD_NORM132], data[TD_NORM132+1], data[TD_NORM132+2]));
					}
				}
				if (hasInterfaces())
					for (int nI = 0, neighbours2 = neighbours * 2; nI < neighbours2; nI += 2) {
						int intf_flag = neighbour[nI + 1];
						if ((intf_flag & 15) != 0) {
							sb.append("SHARES FACE ");
							switch (intf_flag & 15) {
							case 1: sb.append("012 [n" + nodeRef[0] + ",n" + nodeRef[1] + ",n" + nodeRef[2] + "]"); break;
							case 2:	sb.append("023 [n" + nodeRef[0] + ",n" + nodeRef[2] + ",n" + nodeRef[3] + "]"); break;
							case 4: sb.append("031 [n" + nodeRef[0] + ",n" + nodeRef[3] + ",n" + nodeRef[1] + "]"); break;
							case 8: sb.append("132 [n" + nodeRef[1] + ",n" + nodeRef[3] + ",n" + nodeRef[2] + "]"); break;
							}
							sb.append(" WITH FACE ");
							switch ((intf_flag >> 16) & 15) {
							case 1: sb.append("012"); break;
							case 2:	sb.append("023"); break;
							case 4: sb.append("031"); break;
							case 8: sb.append("132"); break;
							}
							sb.append(" OF ELEMENT " + neighbour[nI] + "\n");
						} else {
							sb.append("Shares edge ");
							switch (intf_flag & 1008) {
							case 16: 	sb.append("01 [n" + nodeRef[0] + ",n" + nodeRef[1] + "]"); break;
							case 32:	sb.append("02 [n" + nodeRef[0] + ",n" + nodeRef[2] + "]"); break;
							case 64: 	sb.append("03 [n" + nodeRef[0] + ",n" + nodeRef[3] + "]"); break;
							case 128: 	sb.append("12 [n" + nodeRef[1] + ",n" + nodeRef[2] + "]"); break;
							case 256: 	sb.append("13 [n" + nodeRef[1] + ",n" + nodeRef[3] + "]"); break;
							case 512: 	sb.append("23 [n" + nodeRef[2] + ",n" + nodeRef[3] + "]"); break;
							}
							sb.append(" with edge ");
							switch ((intf_flag >> 16) & 1008) {
							case 16: 	sb.append("01"); break;
							case 32:	sb.append("02"); break;
							case 64: 	sb.append("03"); break;
							case 128: 	sb.append("12"); break;
							case 256: 	sb.append("13"); break;
							case 512: 	sb.append("23"); break;
							}
							sb.append(" of element " + neighbour[nI] + "\n");
						}
					}
				break;
			case HEXAHEDRON:
			case ISOTETRAHEDRON:
			case ISOHEXAHEDRON:
		}
		
		sb.append("Neighbours: [");
		if (neighbours == 0 || neighbour == null) sb.append("n/a]\n");
		else {
			if (neighboursF > 0) sb.append("F: ");
			for (int n = 0, nF = neighboursF; n < neighbours; n++, nF--) {
				if (hasInterfaces()) {
					if (nF == 0 && neighboursF < neighbours) sb.append("e: ");
					sb.append("T" + neighbour[n * 2] + (n == neighbours - 1 ? "]\n" : ", "));
				} else {
					if (nF == 0 && neighboursF < neighbours) sb.append("e: ");
					sb.append("T" + neighbour[n] + (n == neighbours - 1 ? "]\n" : ", "));
				}
			}
		}
		sb.append("\n");
		return sb.toString();
	}
	
	public static String toBitsInteger(int i, String bitfield) {
		StringBuilder sb = new StringBuilder();
		sb.append(bitfield);		
		for (int b = 31; b >= 0; b--) if ((i & (1 << b)) != 0) sb.append("1"); else sb.append(" ");
		sb.append("\n"); return sb.toString();
	}
}
