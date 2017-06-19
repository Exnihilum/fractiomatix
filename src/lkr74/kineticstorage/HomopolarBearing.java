package lkr74.kineticstorage;

public class HomopolarBearing {

	int magnets = 2;							// how many magnets in a bearing row
	double dIr, dOr;							// rotor inner & outer diameter
	double dImO, dOmO;							// outer magnet inner & outer diameter (width increases stiffness quadratically)
	double tm;									// magnet thickness (thicker -> lower speed operation)
	double dIepO, dOepO, tep;					// outer end plate inner & outer diameter and thickness
	double dIipsO, dOipsO, tips;				// outer pole shoe inner & outer diameter and thickness
	double airgap, airgapE;						// airgap & emergency airgap
	double hpolarFluxComp = 0.244;				// homopolar flux component in Teslas
	
	double conductivityCopper = 5.96*10e-7;		// copper conductivity, siemens/m
	double conductivityAluminium = 3.5*10e-7;	// aluminium conductivity
	double conductivityGraphene = 1*10e-8;		// graphene conductivity
	
	// instantiate single row homopolar bearing model
	public HomopolarBearing(int magnets, double dImO, double dOmO, double tm, double airgap, double lowestRotSpeed) {
		this.magnets = magnets;
		this.dImO = dImO; this.dOmO = dOmO; this.tm = tm;
		this.airgap = airgap * .5;
		this.airgapE = airgap * .5;						// set aside half of the airgap for emergency airgap
		this.dOr = dImO - airgap*2;
		this.dIr = dOr - leastRotorThickness(lowestRotSpeed, conductivityCopper) * 2;
		this.dIepO = dImO; this.dOepO = dOmO;			// let end plats be same width as magnets initially
		this.dIipsO = dImO; this.dOipsO = dOmO*0.9; 	// leave gap on pole shoe outer edge to stop flux leaking towards housing)
	}
	
	static double leastRotorThickness(double lowestOpSpeed, double conductivity) {
		return Math.sqrt(2/(lowestOpSpeed*conductivity*4*Math.PI*10e-7));
	}
	double innerRadialFluxGradient(double remFluxD) {			// remanent flux density Br as parameter
		double rImO = dImO*.5, tr = (dOr-dIr)*.5;				// magnet inner radius & rotor thickness
		return (remFluxD*(rImO-airgap-tr-airgapE) - remFluxD*(rImO-airgapE))/
								(tr+airgap*2);					// divide by the operationg interval
	}
	double outerRadialFluxGradient(double remFluxD) {			// remanent flux density Br as parameter
		double rOmO = dOmO*.5, tr = (dOr-dIr)*.5;				// magnet onner radius & rotor thickness
		return (remFluxD*(rOmO+airgapE) - remFluxD*(rOmO+airgap+airgapE))/
								(tr+airgap*2);					// divide by the operationg interval
	}
	double intermediateRadialFluxGradient(double remFluxD) {	// remanent flux density Br as parameter
		double tr = (dOr-dIr)*.5;
		return -2*remFluxD/(tr+2*(airgap+airgapE));
	}
	// displRotor is rotor's axial displacement from zero along y-axis, displAngle is the rotated angle away from x-axis
	double intermediateNormalFluxTimeDerivative(double rotSpeed, double time, double fluxB0, double displRotor, double displAngle) {
		return rotSpeed*displRotor*-intermediateRadialFluxGradient(fluxB0)*Math.cos(rotSpeed*time+displAngle);
	}
	double intermediateNormalFluxD(double rotSpeed, double time, double fluxB0, double displRotor, double displAngle) {
		return rotSpeed*displRotor*-intermediateRadialFluxGradient(fluxB0)*Math.sin(rotSpeed*time+displAngle) + hpolarFluxComp;		
	}
	double intermediateRemanentFluxD(double rotSpeed, double time, double fluxB0, double displRotor, double displAngle) {
		double fluxBr = intermediateNormalFluxD(rotSpeed, time, fluxB0, displRotor, displAngle); return fluxBr;
	}
	//double effectivePoleArea()
}
