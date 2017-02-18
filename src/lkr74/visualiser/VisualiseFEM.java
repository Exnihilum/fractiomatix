package lkr74.visualiser;

import processing.core.*;
import processing.opengl.*;

public class VisualiseFEM extends PApplet {

	float gridWidth = 100f, gridStep = 10f;
	float fontSize = 30;
	PFont f;
	
	public static void main(String[] args) { PApplet.main("lkr74.visualiser.VisualiseFEM"); }


	public void settings() {
		size(800, 800, P3D);
	}
	
	
	public void setup() {
		  // Step 3: Create Font
		  f = createFont("Arial", fontSize);
		} 

	public	void draw() {
		background(255);
		
		// draw a point grid
//		for (float x = -gridWidth; x < gridWidth; x += gridStep)
//			for (float y = -gridWidth; y < gridWidth; y += gridStep)
//				point(x, y, 0);

		// draw a line grid
		for (float x = -gridWidth; x < gridWidth; x += gridStep)
			line(x, -gridWidth, x, gridWidth);
		for (float y = -gridWidth; y < gridWidth; y += gridStep)
			line(-gridWidth, y, gridWidth, y);

		camera(60f, 0f, -30f, 0f, 0f, 0f, 0f, 0f, 1f);
		//frustum(-40, -40, 40, 40, 1, 200);
		translate(-80, 0, 0);
		//rotateX(-PI/6);
		//rotateY(PI/3);
		fill(220);
		box(45);
		
		// some text
		rotateY(PI/2);
		textFont(f, fontSize); fill(0); text("Some text here", 0, 0);
		}
		
	public void mousePressed() { saveFrame("output.png"); }		// do something on mouse press here

}
