package lkr74.sizeof;

import lkr74.matrixlib.CSRMatrix;
import lkr74.matrixlib.FrontalMatrix;
import lkr74.matrixlib.Matrix;
import lkr74.matrixlib.NSPMatrix;

public class SizeOf {
	private long size = 0;
	
	// Adaptation of Vladimir Roubtsov's method for finding true size of an allocated Object
	// the switch clause covers the measurable object types, their instantiation is made by cloning an existing
	// object that needs measuring
	
	public SizeOf (Object gaugeObject, boolean verbose)
	{
		// Warm up all classes/methods we will use
		runGC ();
		usedMemory ();
		// Array to keep strong references to allocated objects
		int count = 10000, i = -1;
		Object [] objects = new Object [count];

		long maxUseMem = s_runtime.freeMemory() / 8;			// max free memory we allow SizeOf to utilise

		long heap1 = 0;
		String fullClassName = gaugeObject.getClass().getName();
		String className = fullClassName.substring(fullClassName.lastIndexOf('.') + 1);
		
		for (; i < count; i++) {
			if (usedMemory() > maxUseMem) break;				// interrupt on allocating more than 1/8th of free memory
			Object object = null;
			// Instantiate your data here and assign it to object
			switch (className) {
			case "Matrix":			object = ((Matrix)gaugeObject).clone(); break;
			case "CSRMatrix":		object = ((CSRMatrix)gaugeObject).clone(); break;
			case "NSPMatrix":		object = ((NSPMatrix)gaugeObject).clone(); break;
			case "FrontalMatrix":	object = ((FrontalMatrix)gaugeObject).clone(); break;
			default: System.out.println("Unknown SizeOf class.\n"); return;
			}
			if (i >= 0) objects[i] = object;
			else {	object = null;								// Discard the warm up object
					runGC ();
					heap1 = usedMemory(); }						// Take a before heap snapshot
		}
		runGC ();
		long heap2 = usedMemory (); // Take an after heap snapshot:

		size = Math.round (((float)(heap2 - heap1)) / (i - 1));
		if (verbose) {
			System.out.println ("'before' heap: " + heap1 + ", 'after' heap: " + heap2);
			System.out.println ("heap delta: " + (heap2 - heap1) + ", {" + className + "} size = " + size + " bytes");
		}
		for (i = 0; i < count; ++ i) objects [i] = null;
		objects = null;
	}
	
	public long size() { return size; }

	// It helps to call Runtime.gc() using several method calls:
	private static void runGC () { for (int r = 0; r < 4; ++ r) _runGC (); }

	private static void _runGC () {
		long usedMem1 = usedMemory (), usedMem2 = Long.MAX_VALUE;
		for (int i = 0; (usedMem1 < usedMem2) && (i < 500); ++ i) {
			s_runtime.runFinalization ();
			s_runtime.gc ();
			Thread.currentThread ();
			Thread.yield ();
			usedMem2 = usedMem1;
			usedMem1 = usedMemory ();
		}
	}


	private static long usedMemory () { return s_runtime.totalMemory () - s_runtime.freeMemory (); }

	private static final Runtime s_runtime = Runtime.getRuntime ();
}