package lkr74.fem1;

import java.security.InvalidParameterException;
import lkr74.mathgenerics.VisitBitArray;

public class NodeHashTable {
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//			HASHTABLE FOR COORDINATES WITH UNIQUE INDEXES AND ATTACHED REFERENCE VALUES 				//
	//			Leonard Krylov 2017																			//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////	

	private final double[][] table;				// coordinate hash table
	private final int[][] tableI;				// matching index & references & counts table (corresponds exactly to table[][] +2 offset of 2 elements)
	private final double[] fastTable;			// fast circular table of "sizeFT" nodes
	private final int[] fastTableI;				// fast circular table of "sizeFT" indexes + references
	private final int size, sizeB, sizeFT3, mask, nodesPerOrdinate, nodesPerOrdinateSq;
	private final double nodesPerVolScale;		// helps quick calculation of a node's integer coordinates
	private final VisitBitArray nodeFlagger;	// quick-access check of presense of a node
	private final int startIndex;				// the starting node index
	private int nodeIndex=0;					// the running node index
	private int fastTableIndex=0;				// the fast table circular index
	private final double xZero, yZero, zZero;	// no coordinate should be lower than this specified zero-point of the collection
	private boolean fromBack = false;			// tells whether scanning should happen from end of arrays or from beginning
	private boolean useFastTable = false;		// tells whether fast table should be used (activating this on nodecount/bucket heuristic?)
	
	final static int NCOORD = FEM1.NCOORD, NCOORD2 = NCOORD*2;
	
	// for instantiation, NodeHashTable needs the volume scale (from 0 to a positive max.number, assumed to define a cube, equal on x/y/z)
	// and how many nodes per each ordinate to divide the volume scale into
	public NodeHashTable(int startIndex, double volScale, int nodesPerOrdinate, double xZero, double yZero, double zZero) {
		if (nodesPerOrdinate < 1 || volScale < 0)
			throw new InvalidParameterException("NodeHashTable.NodeHashTable(): supplied spatial dimensions invalid.");
		this.startIndex = this.nodeIndex = startIndex;
		long nPerO3D = nodesPerOrdinate * nodesPerOrdinate * nodesPerOrdinate;
		// DEBUG: want to make sure there isn't an out of memory issue since the allocation is the cubic of the provided nodes/ordinate value
		// let a test gird of 256x256x256 units be the maximal direct checklist
		while (nPerO3D > 256*256*256) { nodesPerOrdinate>>=1; nPerO3D = nodesPerOrdinate * nodesPerOrdinate * nodesPerOrdinate; }
		this.nodesPerOrdinate = nodesPerOrdinate;
		nodesPerOrdinateSq = nodesPerOrdinate * nodesPerOrdinate;
		nodesPerVolScale = (double)nodesPerOrdinate / volScale;
		int bucketCount = sizeB = (int)Math.sqrt(nPerO3D) / 32;
		int pow2size =0;
		bucketCount *= 64;						// (bucket size)/32 and (bucket count)*2
		// trim bucket count to nearest upper power-of-2
		while (bucketCount != 0) { pow2size = (pow2size<<1) + 1; bucketCount >>= 1; }
		mask = pow2size++;
		int sizeFT = (pow2size>>1) > 8 ? 8 : pow2size>>1;
		fastTable = new double[sizeFT3 = sizeFT * 3];
		fastTableI = new int[sizeFT3];
		table = new double[size = pow2size][];
		tableI = new int[size][];
		nodeFlagger = new VisitBitArray((int)nPerO3D);
		this.xZero = xZero; this.yZero = yZero; this.zZero = zZero; 
	}
	
	public int count() { return nodeIndex - startIndex; }
	
	// a hashing value for 3D coordinates suggested by David Mauro
	private int hashValue(int x, int y, int z) {
		long max = (x > y ? (x > z ? x : z) : (y > z ? y : z)), hash = max * max * max + 2 * max * z + z;
		if (max == z) { int max_xy = x > y ? x : y; hash += max_xy * max_xy; }
		if (y >= x) hash += x + y; else hash += y;
		return (int)(hash & mask);
	}
	
	public boolean contains(double x, double y, double z) {
		int xInt = (int)((x - xZero) * nodesPerVolScale), yInt = (int)((y - yZero) * nodesPerVolScale), zInt = (int)((z - zZero) * nodesPerVolScale);
		int composeIdx = zInt * nodesPerOrdinateSq + yInt * nodesPerOrdinate + xInt;
		if (!nodeFlagger.visited(composeIdx)) return false;
		return true;
	}

	// if unique node, inserted in table and running index returned, if indicate=true, an existing index with top bit set & reference returned
	// with fineCheck=true, the caller forces fine search for a node, assuming that the node is unique
	public long add(double x, double y, double z, int reference, boolean indicate) {
		
		int xInt = (int)((x - xZero) * nodesPerVolScale), yInt = (int)((y - yZero) * nodesPerVolScale), zInt = (int)((z - zZero) * nodesPerVolScale);
		int composeIdx = zInt * nodesPerOrdinateSq + yInt * nodesPerOrdinate + xInt, bkSize;
		int hashIdx = (18397*xInt + 20483*yInt + 29303*zInt) & mask;
		double[] bucket = table[hashIdx];
		int[] idxBucket = tableI[hashIdx];

		if (nodeFlagger.visited(composeIdx)) {				// if coarse uniqueness check fails
			if (indicate) {									// if existing node indication requested, find the node
				int i3 = 0, sFT2 = 0;						// TODO: do not use fast table on a misses criterion or being >max bucket length
				if (useFastTable)
					while (i3 < sizeFT3) {					// check in fast table
						if (fastTable[i3++] == x)
							if (fastTable[i3++] == y) { if (fastTable[i3++] == z)
								return ((long)fastTableI[sFT2+1])<<32 | 0x80000000L | fastTableI[sFT2]; }
							else i3++;
						else i3 += 2;
						sFT2 += 3;
					}
				if (!fromBack) {
					i3 = 0;
					for (int i = 1, iEnd = 2 + idxBucket[i++]*2; i < iEnd; i += 2) {
						if (bucket[i3++] == x)
							if (bucket[i3++] == y) { if (bucket[i3++] == z) return ((long)idxBucket[i+1])<<32 | 0x80000000L | idxBucket[i]; }
							else i3++;
						else i3 += 2;
					}
				} else {
				// heuristic: search from back, since for spatially arranged vertices the locally relevant nodes will always be in the end
				// TODO: check if it's possible to toggle between starting from beginning or end of array on some criterion, or just doing it randomly?
					i3 = idxBucket[1] * NCOORD - 1;
					for (int i = idxBucket[1]*2; i > 1; i -= 2) {
						if (bucket[i3--] == z)
							if (bucket[i3--] == y) { if (bucket[i3--] == x) return ((long)idxBucket[i+1])<<32 | 0x80000000L | idxBucket[i]; }
							else i3--;
						else i3 -= 2;
					}
				}
			}
			// important note: without indication, the failure is under a COARSER criterion of a node as an int coordinate (approcimated to nodes/ordinate)
			else return -1;														// on node existence and no indication requested, return failure as -1
			
		} else {																// node is unique
			nodeFlagger.visit(composeIdx);										// coarsely flag it's integer coordinate as non-unique
			if (bucket == null) {												// control bucket existense and proper size
				bucket = table[hashIdx] = new double[sizeB * NCOORD];
				idxBucket = tableI[hashIdx] = new int[sizeB * 2 + 2];
				idxBucket[0] = sizeB * 2;
			} else if (idxBucket[0] <= idxBucket[1]) {
				double[] newBucket = new double[(bkSize = idxBucket[0]) * 6];	// make bucket twice the size, copy over the data
				int[] newIdxBucket = new int[2 + 4 * bkSize];
				for (int i = 2, iEnd = bkSize + 2, i3 = 0; i < iEnd;) {
					newBucket[i3] = bucket[i3++]; newBucket[i3] = bucket[i3++]; newBucket[i3] = bucket[i3++];
					newIdxBucket[i] = idxBucket[i++]; newIdxBucket[i] = idxBucket[i++];
				}
				newIdxBucket[0] = bkSize * 2;
				table[hashIdx] = newBucket;
				tableI[hashIdx] = newIdxBucket;
			}
		}
		
		// at this point, either a unique node was received (coarse search failed), or a node was not found on indication request
		// in both cases, need to insert the node
		int i3 = idxBucket[1], i = i3 * 2 + 2; i3 *= NCOORD;
		if (useFastTable) {
			idxBucket[i++] = fastTableI[fastTableIndex] = nodeIndex;		// store the inserted node's index
			bucket[i3++] = fastTable[fastTableIndex++] = x;					// store the node
			idxBucket[i++] = fastTableI[fastTableIndex] = reference;		// store the supplied reference
			bucket[i3++] = fastTable[fastTableIndex++] = y;					// store the node
			fastTableI[fastTableIndex] = hashIdx;							// store hash index for fast removal
			bucket[i3++] = fastTable[fastTableIndex++] = z;					// store the node
			if (fastTableIndex >= sizeFT3) fastTableIndex = 0; 
		} else {
			idxBucket[i++] = nodeIndex;	idxBucket[i++] = reference;			// store the inserted node's index & the supplied reference
			bucket[i3++] = x; bucket[i3++] = y; bucket[i3++] = z;			// store the node
		}
		idxBucket[1]++;														// increase count of that bucket
		return ((long)reference)<<32 | (long)nodeIndex++;					// return long combo of running index & reference value
	}
	
	
	
	
	// method removes a recent node entry found in fast table specifically for IST algorithm
	// method expects to be called for ALL recent nodes upto LAST node index, to avoid breaking the incremental node index consistency
	// therefore, the node index can be trivially DECREMENTED during node removal
	// method returns false if the node wasn't found
	public boolean removeRecentIST(int idx) {
		if (!useFastTable) return false;
		for (int iFT = 0; iFT < sizeFT3; iFT += NCOORD) {
			if (fastTableI[iFT] == idx) {
				fastTableI[iFT] = -1;
				fastTable[iFT++] = Double.MAX_VALUE;					// mark fast table's coordinate as "unattainable"	
				int hashIdx = fastTableI[++iFT];
				double[] bucket = table[hashIdx];
				int[] idxBucket = tableI[hashIdx];
				for (int i = 1, iEnd = 2 + idxBucket[i++]*2; i < iEnd; i += 2) {
					if (idxBucket[i] == idx) {
						idxBucket[i] = -1;								// true data deletion too costly and unnecessary, mark as "deleted"
						// mark x-coordinate (z-coordinate for reverse scanning) as "unattainable"
						if (fromBack)	bucket[((i-2)>>1)*NCOORD + 2] = Double.MAX_VALUE;
						else 			bucket[((i-2)>>1)*NCOORD] = Double.MAX_VALUE;
						break;
					}
				}
				nodeIndex--;
				return true;
			}
		}
		return false;
	}
	
	// returns nodes as one array correctly indexed according to matching index per node
	// two types of node collection possible:
	// 1) method can create extra space in front (inFront=true) or back of array
	// 2) method can take a node condensing reindexer array nodeRemap[] that reorders indexes, remapSize will tell the size of a condensed array
	// startNode will tell what node index to start collecting from
	public double[] arrayWithSpace(int extraSpace, boolean inFront, int startNode, int[] nodeRemap, int remapSize, boolean interleave) {
		extraSpace *= NCOORD;
		double[] nodeArray = new double[(nodeRemap==null ? ((interleave?2:1)*nodeIndex*NCOORD) : remapSize*NCOORD) + extraSpace];
		for (int b = 0; b < size; b++) {
			double[] bucket = table[b];
			int[] idxBucket = tableI[b];
			if (bucket == null) continue;
			if (nodeRemap != null)
				for (int i = 0, i3 = 0, ip2 = 2, iEnd = idxBucket[i+1]; i < iEnd; i++, ip2+=2) {
					int nIdx = idxBucket[ip2];
					if (nIdx < 0 || nIdx < startNode) { i3 += NCOORD; continue;	}		// skip "deleted" nodes and node indexes < startNode
					int i3a = nodeRemap[nIdx] * NCOORD;
					nodeArray[i3a++] = bucket[i3++]; nodeArray[i3a++] = bucket[i3++]; nodeArray[i3a] = bucket[i3++];
				}
			else {
				if (interleave) {
					for (int i = 0, i3 = 0, ip2 = 2, iEnd = idxBucket[i+1]; i < iEnd; i++, ip2+=2) {
						int nIdx = idxBucket[ip2];
						if (nIdx < 0 || nIdx < startNode) { i3+=NCOORD; continue;	}	// skip "deleted" nodes and node indexes < startNode
						int i6a = inFront ? nIdx * NCOORD2 + extraSpace : nIdx * NCOORD2;
						nodeArray[i6a++] = bucket[i3++]; nodeArray[i6a++] = bucket[i3++]; nodeArray[i6a] = bucket[i3++];
					}
				} else {
					for (int i = 0, i3 = 0, ip2 = 2, iEnd = idxBucket[i+1]; i < iEnd; i++, ip2+=2) {
						int nIdx = idxBucket[ip2];
						if (nIdx < 0 || nIdx < startNode) { i3+=NCOORD; continue;	}	// skip "deleted" nodes and node indexes < startNode
						int i3a = inFront ? nIdx * NCOORD + extraSpace : nIdx * NCOORD;
						nodeArray[i3a++] = bucket[i3++]; nodeArray[i3a++] = bucket[i3++]; nodeArray[i3a] = bucket[i3++];
					}
				}
			}
		}
		return nodeArray;
	}
	
	public double[] array() { double[] nodeArray = arrayWithSpace(0, false, 0, null, 0, false); return nodeArray; }
	public double[] arrayInterleaved() { double[] nodeArray = arrayWithSpace(0, false, 0, null, 0, true); return nodeArray; }
	public double[] array(int startNode) { double[] nodeArray = arrayWithSpace(0, false, startNode, null, 0, false); return nodeArray; }
	public double[] array(int extraSpace, boolean inFront, int startNode) {
		double[] nodeArray = arrayWithSpace(extraSpace, inFront, startNode, null, 0, false); return nodeArray; }
	public double[] array(int startNode, int[] nodeRemap, int remapSize) {
		double[] nodeArray = arrayWithSpace(0, false, startNode, nodeRemap, remapSize, false); return nodeArray; }
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		int maxBucket = 0, barMaxLength = 60;
		for (int i = 0; i < size; i++) if (tableI[i]!= null && tableI[i][1] > maxBucket) maxBucket = tableI[i][1];	// find largest containing bucket count
		if (maxBucket == 0) return sb.toString();
		sb.append("Node Hash Table distribution:\n");
		for (int i = 0; i < size; i++) {
			if (tableI[i] != null) {
				int barLength = maxBucket < barMaxLength ? tableI[i][1] : (tableI[i][1] * barMaxLength) / maxBucket;
				sb.append(String.format("%-12d", tableI[i][1]));
				for (int bar = 0; bar < barLength; bar++) sb.append("#");
			}
			sb.append("\n");
		}
		return sb.toString();
	}

}
