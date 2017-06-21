package lkr74.mathgenerics;

import java.security.InvalidParameterException;

// RandFill class is supplied with an (integer) range and will guarantee to return a random integer value within that range
// but never the same value as any previous one, thus populating the range evenly with random values
// a useful class for iteratively test-filling arrays & matrices with random values in a non-repeating fashion
// "sectors" holds integer pairs, the first is start index of a sector, the second an end index
// there is always at least one occupied slot between each sector, every new occupied slot can create a new sector partition
//
// Leonard Krylov 2016

public class RandFill {
	private int sectorCnt = 1, slotCount, startRange, endRange;
	private int[] sectors;
	
	public RandFill(int startRange, int endRange) {
		if (startRange > endRange || endRange < 0) throw new InvalidParameterException("RandFill(): Invalid range.");
		this.startRange = startRange;
		this.endRange = endRange;
		slotCount = endRange - startRange == 0 ? 1 : endRange - startRange;
		sectors = new int[slotCount + 1];		// holds definitions of sectors between occupied random value slots
		sectors[1] = slotCount - 1;				// define the first default sector as covering entire range
	}
	
	public int remainingSlots() { return slotCount; }
	
	public void reset() {
		for (int i = 1; i < sectors.length; i++) sectors[i] = 0;
		sectors[1] = sectors.length - 2;
	}
	
	public int getRandom() {
		
		// zero sectors left means we've exhausted all random slots, return 0
		if (sectorCnt == 0) return 0;
		
		int rndSec = (int)(Math.random() * sectorCnt) * 2;			// select a random sector
		int secLen = sectors[rndSec + 1] - sectors[rndSec] + 1;		// get it's length
		int rnd;
		
		// did we find a single-element sector?
		if (secLen == 1) {
			sectorCnt--;
			slotCount--;
			rnd = sectors[rndSec];
			sectors[rndSec] = sectors[sectorCnt * 2];				// delete this 1-element sector by copying last-in-array sector over it
			sectors[rndSec + 1] = sectors[sectorCnt * 2 + 1];
			return startRange + rnd;								// return this random slot
		}
		// select a random slot within this sector
		rnd = (int)(Math.random() * secLen) + sectors[rndSec];
		if (rnd == sectors[rndSec])	sectors[rndSec]++;				// slot is at start of sector, increment sector start point
		else if (rnd == sectors[rndSec + 1]) sectors[rndSec + 1]--;	// slot is at end of sector, decrement sector end point
		else {														// slot is in middle of sector, partition into two sectors
			sectors[sectorCnt * 2] = rnd + 1;						// create new sector at end of sector array
			sectors[sectorCnt * 2 + 1] = sectors[rndSec + 1];		// let it start above the slot
			sectors[rndSec + 1] = rnd - 1;							// decrement current sector below the slot
			sectorCnt++;
		}
		slotCount--;
		return startRange + rnd;
	}
	
	@Override
	public RandFill clone() {
		Object o = null;
		try { o = super.clone(); } catch (CloneNotSupportedException e) { e.printStackTrace(); }
		RandFill rf = (RandFill) o;
		return rf;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Range: " + startRange + "-" + endRange + "\n");
		sb.append("Sectors: ");
		for (int s = 0, s2 = 0; s < sectorCnt; s++)
			sb.append("[" + (sectors[s2++]+startRange) + "-" + (sectors[s2++]+startRange) + "]");
		return sb.toString();
	}
}
