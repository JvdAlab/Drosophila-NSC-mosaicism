package uk.ac.cam.gurdon.deadpanorama;

import ij.ImagePlus;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.GaussianBlur3D;
import ij.plugin.ImageCalculator;

public class DoG {

	//3D Difference of Gaussians
	public static ImagePlus process(ImagePlus imp, int chan, double sigma, double k){

		Duplicator dup = new Duplicator();
		int Z = imp.getNSlices();
		Calibration cal = imp.getCalibration();
		ImagePlus map = dup.run(imp,chan,chan,1,Z,1,1);
		ImagePlus sub = dup.run(map,1,1,1,Z,1,1);
		
		double sigmaXY = sigma/cal.pixelWidth;
		double sigmaZ = sigma/cal.pixelDepth;
		GaussianBlur3D.blur( map, sigmaXY, sigmaXY, sigmaZ );
		GaussianBlur3D.blur( sub, sigmaXY*k, sigmaXY*k, sigmaZ*k );
		
		ImageCalculator ic = new ImageCalculator();
		ImagePlus dog = ic.run("Subtract create stack", map, sub);
		map.close();
		sub.close();
		
		dog.setTitle("DoG_"+imp.getTitle());
		return dog;
	}
	
}
