package kineticstorage;


// Singleton initialisation of a rim material database
class FESRimMaterial {
	
	final String[] materialData = {
			"carbon T300", "carbon T400HB", "carbon T700SC", "carbon T800HB", "carbon T1000GB",
			"carbon M35JB", "carbon M40JB", "carbon M46JB", "carbon M60JB" };
	// tensile hoop strength in Pascals
	final double[] tensileData = {	3530e6, 4410e6, 4900e6, 5490e6, 6370e6,
									4700e6, 4400e6, 4020e6, 3820e6 };
	// tensile modulus in Pascals
	final double[] modulusData = {	230e9, 250e9, 230e9, 294e9, 294e9,
									343e9, 377e9, 436e9, 588e9 };
	// density in kg/m^3
	final double[] densityData = {	1760, 1800.0, 1800.0, 1810.0, 1800.0,
									1750.0, 1750.0, 1840.0, 1930.0 };
	// elongation in %
	final double[] elongationData = {	101.5, 101.8, 102.1, 101.9, 102.2,
										101.4, 101.2, 100.9, 100.7 };
	
	public double tensile, modulus, density, elongation;
	
	FESRimMaterial(int index) {
		checkMaterialValues(index);
		tensile = tensileData[index];
		modulus = modulusData[index];
		density = densityData[index];
		elongation = elongationData[index];
	}
	
	FESRimMaterial(String name) {
		for (int i = 0; i < materialData.length; i++)
			if (materialData[i].equals(name)) {
				checkMaterialValues(i);
				tensile = tensileData[i];
				modulus = modulusData[i];
				density = densityData[i];
				elongation = elongationData[i];
				return;
			}
		throw new RuntimeException("FESRimMaterial(): database lacks material:" + name);		
	}
	
	private void checkMaterialValues(int index) {
		if (tensileData[index] < 300e6 || tensileData[index] > 1000e12)
			throw new RuntimeException("FESRimMaterial(): faulty tensile strength:" + tensileData[index]);
		if (modulusData[index] < 300e6 || modulusData[index] > 1000e14)
			throw new RuntimeException("FESRimMaterial(): faulty elasticity modulus:" + modulusData[index]);
		if (densityData[index] < 100.0 || densityData[index] > 19200.0)
			throw new RuntimeException("FESRimMaterial(): faulty density:" + densityData[index]);
		if (elongationData[index] < 100.0 || elongationData[index] > 140.0)
			throw new RuntimeException("FESRimMaterial(): faulty elongation:" + elongationData[index]);
	}

}
