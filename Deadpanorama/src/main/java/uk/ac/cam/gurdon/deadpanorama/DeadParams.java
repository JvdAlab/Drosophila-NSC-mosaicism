package uk.ac.cam.gurdon.deadpanorama;

import java.text.SimpleDateFormat;
import java.util.Date;

import ij.IJ;
import ij.Prefs;

public class DeadParams {

	public static final double DEFAULT_K = 4.0;
	public static final int DEFAULT_STARTK = 16;
	public static final int DEFAULT_JOINZ = 1;
	public static final double DEFAULT_JOINOVERLAP = 0.2;
	public static final double DEFAULT_MINFLATNESS = 0.2;
	public static final boolean DEFAULT_WATERSHED = false;
	
	public int deadpanC, startZ, endZ, startK, joinZ;
	public double sigma, K, minA, maxA, threshold, minCirc, joinOverlap, minFlatness;
	public boolean advanced, watershed;

	
	public DeadParams(){
		
	}
	
	public static DeadParams getPrefs(){
		DeadParams params = new DeadParams();
		params.deadpanC = Prefs.getInt("DeadParams.channel", 1);
		params.startZ = Prefs.getInt("DeadParams.startZ", 1);
		params.endZ = Prefs.getInt("DeadParams.endZ", 100);
		params.startK = Prefs.getInt("DeadParams.startK", DEFAULT_STARTK);
		params.sigma = Prefs.getDouble("DeadParams.sigma", 1.0);
		params.K = Prefs.getDouble("DeadParams.K", DEFAULT_K);
		params.minA = Prefs.getDouble("DeadParams.minA", 10.0);
		params.maxA = Prefs.getDouble("DeadParams.maxA", 100.0);
		params.threshold = Prefs.getDouble("DeadParams.threshold", 0.0);
		params.minCirc = Prefs.getDouble("DeadParams.minCirc", 0.6);
		params.joinZ = Prefs.getInt("DeadParams.joinZ", DEFAULT_JOINZ);
		params.joinOverlap = Prefs.getDouble("DeadParams.joinOverlap", DEFAULT_JOINOVERLAP);
		params.minFlatness = Prefs.getDouble("DeadParams.minFlatness", DEFAULT_MINFLATNESS);
		params.advanced = Prefs.getBoolean("DeadParams.advanced", false);
		params.watershed = Prefs.getBoolean("DeadParams.watershed", DEFAULT_WATERSHED);
		
		return params;
	}
	
	public void save(){
		
		Prefs.set("DeadParams.channel", deadpanC);
		Prefs.set("DeadParams.startZ", startZ);
		Prefs.set("DeadParams.endZ", endZ);
		Prefs.set("DeadParams.startK", startK);
		Prefs.set("DeadParams.sigma", sigma);
		Prefs.set("DeadParams.K", K);
		Prefs.set("DeadParams.minA", minA);
		Prefs.set("DeadParams.maxA", maxA);
		Prefs.set("DeadParams.threshold", threshold);
		Prefs.set("DeadParams.minCirc", minCirc);
		Prefs.set("DeadParams.joinZ", joinZ);
		Prefs.set("DeadParams.joinOverlap", joinOverlap);
		Prefs.set("DeadParams.minFlatness", minFlatness);
		Prefs.set("DeadParams.advanced", advanced);
		Prefs.set("DeadParams.watershed", watershed);
		
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		String nl = System.getProperty("line.separator");
		SimpleDateFormat formatter= new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");  
		Date date = new Date(System.currentTimeMillis());  
		sb.append("Deadpanorama Parameters "+formatter.format(date)+nl);
		
		if(IJ.getImage()!=null)	sb.append("Image: "+IJ.getImage().getTitle()+nl);
		sb.append("Deadpan Channel: "+ deadpanC+nl);
		sb.append("Slice Range: "+ startZ +"-"+ endZ+nl);
		sb.append("Sigma: "+ sigma+nl);
		sb.append("Sigma Factor: "+ K+nl);
		sb.append("Starting K: "+ startK+nl);
		sb.append("Watershed: "+ watershed+nl);
		sb.append("Area Range: "+ minA +"-"+ maxA+nl);
		sb.append("Threshold: "+ threshold+nl);
		sb.append("Min Circularity: "+ minCirc+nl);
		sb.append("Min Flatness: "+ minFlatness+nl);
		sb.append("Max Slice Distance: "+ joinZ+nl);
		sb.append("Min Join Overlap: "+ joinOverlap+nl);
		sb.append(nl);
		
		return sb.toString();
	}
	
}
