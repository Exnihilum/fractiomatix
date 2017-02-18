package lkr74.kineticstorage;

public class FES {		
	
	String name;
	int rims;
	public double height, speed;
	double[] rimRadius;
	FESRimMaterial[] rimMaterial;
	
	public FES(String name, int rims, double height, double speed) {
		if (rims <= 0) throw new RuntimeException("FlyWheel(): invalid rim count.");
		this.name = name;
		this.rims = rims;
		rimRadius = new double[rims + 1];
		// insert some default radius measurements
		for (int i = 0; i < rims; i++) rimRadius[i] = 0.5 + 0.25 * i;
		
		rimMaterial = new FESRimMaterial[rims];
		
		if (height <= 0) throw new RuntimeException("FlyWheel(): invalid rim height.");
		this.height = height;

		if (speed <= 0) throw new RuntimeException("FlyWheel(): invalid speed.");
		this.speed = speed;
	}
	
	public void setRim(int rim, double radiusI, double radiusO, String materialName) {
		
		if (rim < 0 || rim > rims) throw new RuntimeException("FlyWheel.setRim(): invalid rim.");
		if (rim != 0 && radiusI < rimRadius [rim - 1]) throw new RuntimeException("FlyWheel.setRim(): radius surpasses inner rim.");
		if (radiusI >= radiusO) throw new RuntimeException("FlyWheel.setRim(): invalid radii.");

		rimRadius[rim] = radiusI;
		double radiusD = radiusO - rimRadius[rim + 1];
		// readjust all rims past the updated one
		for (int r = rim + 1; r < rimRadius.length; r++) rimRadius[r] += radiusD;
		
		rimMaterial[rim] = new FESRimMaterial(materialName);
	}
	
	public double inertia() {
		double inertia = 0;
		for (int r = 0; r < rims; r++) {
			double rdsI = rimRadius[r], rdsO = rimRadius[r + 1];
			inertia += rimMaterial[r].density * (rdsO*rdsO*rdsO*rdsO - rdsI*rdsI*rdsI*rdsI);
		}
		return (Math.PI / 2.0) * height * inertia;
	}
	
	public double energy() {
		double energy = 0;
		for (int r = 0; r < rims; r++) {
			double rdsI = rimRadius[r], rdsO = rimRadius[r + 1];
			energy += 	0.25 * rimMaterial[r].density * Math.PI * height
						* (rdsO*rdsO*rdsO*rdsO - rdsI*rdsI*rdsI*rdsI)
						* speed * speed;
		}
		return energy;
	}
	
	// returns the max speed attainable for the material combination
	public double speedHoopStrLimit() {
		// we're interested in outermost rim's hoop strength
		int matIdx = rimMaterial.length - 1;
		return Math.sqrt(rimMaterial[matIdx].tensile / rimMaterial[matIdx].density);
	}
	
	
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("Kinetic storage system " + name + String.format("\nInertia: %.4f kg*m^2\n", inertia()));
		sb.append(String.format("Energy at %.4f rad/s (%.0frpm): %.4f Joules\n", speed, speed*60/(Math.PI*2), energy()));
		double spLim = speedHoopStrLimit();
		sb.append(String.format("Hoop strength limited speed: %.4f rad/s (%.0f rpm)\n", spLim, spLim*60/(Math.PI*2)));
		return sb.toString();
	}
	
}
