package uk.ac.cam.gurdon.deadpanorama;

import java.awt.Color;
import java.awt.Font;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;

import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.gui.TextRoi;
import ij.plugin.RoiEnlarger;
import ij.process.ImageProcessor;

public class Nucleus3D {

	private static final Font LABELFONT = new Font(Font.SANS_SERIF, Font.BOLD, 12);
	
	private ImagePlus imp;
	
	private ArrayList<Roi> rois, surroundingRois;
	public int index;
	
	private double[] centroid;
	private ArrayList<ArrayList<Float>> nucleusValues;
	private ArrayList<Float> nucleusMeans;
	private ArrayList<ArrayList<Float>> surroundingValues;
	private ArrayList<Float> surroundingMeans;
	private double volume = 0d;
	
	public Nucleus3D(ImagePlus imp, int index){
		this.imp = imp;
		this.index = index;
		rois = new ArrayList<Roi>();
	}
	
	public ArrayList<Roi> getRois(){
		return rois;
	}
	
	public int size(){
		return rois.size();
	}
	
	public void add(Roi roi){
		rois.add(roi);
	}
	
	public boolean contains(Roi roi){
		return rois.contains(roi);
	}
	
	public void setIndex(int i){
		index = i;
	}
	
	public void addToOverlay(Overlay ol){

		float h = (float) Math.random();
		float s = 1.0f;
		float v = 1.0f;
		Color colour = Color.getHSBColor(h,s,v);
		
		for(Roi roi : rois){
			roi.setStrokeColor(colour);
			ol.add(roi);
			Rectangle roiRect = roi.getBounds();
			TextRoi label = new TextRoi( roiRect.x, roiRect.y-12, ""+index, LABELFONT );
			label.setPosition( -1, roi.getZPosition(), 1 );
			label.setStrokeColor(colour);
			ol.add(label);
		}
		
		if(surroundingRois!=null){	//created when surrounding values are measured
			for(Roi ring : surroundingRois){
				ring.setStrokeColor(colour);
				ol.add(ring);
			}
		}
		
	}

	//calculate centroid, volume and voxel intensity value Lists
	private void calculate() {
		centroid = new double[3];	//uncalibrated values, calibration applied when added to ResultsTable
		volume = 0d;
		int C = imp.getNChannels();
		nucleusValues = new ArrayList<ArrayList<Float>>();
		for(int c=0;c<=C;c++){	//add empty list at 0 for 1-based channel indexing
			ArrayList<Float> vc = new ArrayList<Float>();
			nucleusValues.add(vc);
		}
		for(Roi roi : rois){
			double[] centroid2D = roi.getContourCentroid();
			int roiZ = roi.getZPosition();
			double roiA = roi.getStatistics().area;
			centroid[0] += centroid2D[0];
			centroid[1] += centroid2D[1];
			centroid[2] += roiZ;
			volume += roiA;
			Point[] pixels = roi.getContainedPoints();
			for(int c=1;c<=C;c++){	//get all contained pixel values for each channel
				ImageProcessor ip = imp.getStack().getProcessor( imp.getStackIndex(c,roiZ,1) );
				for(Point p:pixels){
					float f = ip.getf(p.x,p.y);
					nucleusValues.get(c).add(f);
				}
			}
		}
		float nf = rois.size();
		centroid[0] /= nf;
		centroid[1] /= nf;
		centroid[2] /= nf;
	}
	
	public double[] getCentroid(){
		if(centroid==null) calculate();
		return centroid;
	}
	
	public int[] getIntCentroid(){
		if(centroid==null) calculate();
		int[] ic = new int[3];
		ic[0] = (int)centroid[0];
		ic[1] = (int)centroid[1];
		ic[2] = (int)centroid[2];
		return ic;
	}
	
	public double getVolume(){
		if(volume==0d) calculate();
		return volume;
	}
	
	public float getInsideMean(int channel){
		if(nucleusMeans==null){
			nucleusMeans = new ArrayList<Float>();
			int C = imp.getNChannels();
			for(int c=0;c<=C;c++){
				Float[] values = getInsideValues(c);
				float channelMean = (float) (Arrays.stream(values).mapToDouble(f->f).sum() / (float)values.length);
				nucleusMeans.add(channelMean);
			}
		}
		return nucleusMeans.get(channel);
	}
	
	public float getSurroundingMean(int channel){
		if(surroundingMeans==null){
			getSurroundingValues(channel);
		}
		return surroundingMeans.get(channel);
	}
	
	public float getMean(int channel, boolean outer){
		if(outer){
			return getSurroundingMean(channel);
		}
		else{
			return getInsideMean(channel);
		}
	}
	
	public Float[] getInsideValues(int channel){
		if(nucleusValues==null) calculate();
		if(nucleusValues.size()<channel-1) return null;
		return nucleusValues.get(channel).toArray( new Float[nucleusValues.get(channel).size()] );
	}
	
	public Float[] getSurroundingValues(int channel){	//separate from calculate() so only run when needed
		if(surroundingValues==null){
			int C = imp.getNChannels();
			surroundingValues = new ArrayList<ArrayList<Float>>();
			surroundingMeans = new ArrayList<Float>();
			for(int c=0;c<=C;c++){	//add empty list at 0 for 1-based channel indexing
				ArrayList<Float> vc = new ArrayList<Float>();
				surroundingValues.add(vc);
				surroundingMeans.add(0f);
			}
			
			surroundingRois = new ArrayList<Roi>();
			for(Roi roi : rois){
				int roiZ = roi.getZPosition();
				double roiA = roi.getStatistics().area;
				double expand = Math.sqrt(roiA/Math.PI);	//expand by the radius of a circle with the same area
				Roi big = RoiEnlarger.enlarge(roi, expand);
				Roi outerRing = new ShapeRoi(big).not(new ShapeRoi(roi));
				Point[] pixels = outerRing.getContainedPoints();
	
				for(int c=1;c<=C;c++){
					
					ImageProcessor ip = imp.getStack().getProcessor( imp.getStackIndex(c,roiZ,1) );
					int W = ip.getWidth();
					int H = ip.getHeight();
					float mean = 0f;
					for(Point p:pixels){
						if(p.x>=0&&p.y>=0&&p.x<W&&p.y<H){
							float f = ip.getf(p.x,p.y);
							surroundingValues.get(c).add(f);
							mean += f;
						}
					}
					mean /= (float)surroundingValues.get(c).size();
					surroundingMeans.set(c, mean);

				}
				
				outerRing.setPosition( -1, roiZ, 1 );
				surroundingRois.add(outerRing);
				
			}
		}
		return surroundingValues.get(channel).toArray(new Float[surroundingValues.get(channel).size()]);
	}
	
	public Float[] getValues(int channel, boolean outer){
		if(outer){
			return getSurroundingValues(channel);
		}
		else{
			return getInsideValues(channel);
		}
	}
	
}
