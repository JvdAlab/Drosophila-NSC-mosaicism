package uk.ac.cam.gurdon.deadpanorama;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.IntStream;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.plugin.filter.EDM;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.filter.ThresholdToSelection;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.StackStatistics;

public class Segmenter {

	private ExecutorService executor;
	private Set<Future<ArrayList<Roi>>> futures;
	private boolean cancel = false;
	private StatusPanel status;
	
	public Segmenter(StatusPanel status){
		this.status = status;
	}
	
	/**	Maps nuclei - applies Processor, 2D hierarchical K-means segmentation and joins Rois in 3D to create a list of <code>Nucleus3D</code>.
	 * 
	 * @param imp	the raw ImagePlus
	 * @param params	parameters to be used
	 * @param userRoi	the selected area to be analysed, if null the whole image will be used
	 * 
	 * @return	an ArrayList of <code>Nucleus3D</code>
	 */
	public ArrayList<Nucleus3D> getNuclei(ImagePlus imp, DeadParams params, Roi userRoi){
		ArrayList<Nucleus3D> nuclei = null;
		try{
			if(cancel) return null;
status.log("Segmenting nuclei...");
			ImagePlus proc = DoG.process(imp, params.deadpanC, params.sigma, params.K);	//returns a DoG processed copy of the deadpan channel

			/*if(true){
				proc.show();
				return null;
			}*/
			
			HistogramCluster hc = new HistogramCluster(proc);
			double[] levels = hc.getLevels(params.startK);
status.log("Agglomerative K-means gave "+levels.length+" levels");
			if(levels.length==0){
				return null;
			}
			ArrayList<Roi> roiList = getRois(proc, levels, params);
status.log("Got "+roiList.size()+" ROIs");
			if(cancel) return null;
			if(userRoi!=null){
				roiList = filterRois(roiList, userRoi);
			}
			nuclei = joinRois(imp, roiList, params);
			
			proc.close();
			proc.flush();
		}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
		return nuclei;
	}
	
	public synchronized void cancel(){
		cancel = true;
		if(executor!=null){
			executor.shutdownNow();
		}
		if(futures!=null){
			for(Future<ArrayList<Roi>> task : futures){
				task.cancel(true);
			}
		}
	}
	
	/**	Get Rois with a centroid inside the user defined Roi
	 */
	private static ArrayList<Roi> filterRois(ArrayList<Roi> rois, Roi userRoi){
		ArrayList<Roi> filtered = new ArrayList<Roi>();
		for(Roi roi : rois){
			double[] centroid = roi.getContourCentroid();
			if(userRoi.contains((int)centroid[0], (int)centroid[1])){
				filtered.add(roi);
			}
		}
		return filtered;
	}
	
	/** Segments a 2D slice using a hierarchy of intensity levels. Implements Callable to be run in parallel for multiple slices and return Rois from a Future.
	 */
	private class SliceSegmenter implements Callable<ArrayList<Roi>>{

		ImagePlus imp;
		int z;
		double[] levels;
		DeadParams params;
		ImageProcessor procip;
		
		/**Segments a 2D slice using a hierarchy of intensity levels.
		 * 
		 * @param imp	a whole processed z-stack with one channel
		 * @param z	the slice index to segment (1-based)
		 * @param levels	hierarchy of intensity levels
		 * @param params	parameters to be used
		 */
		public SliceSegmenter(ImagePlus imp, int z, double[] levels, DeadParams params){
			this.imp = imp;
			this.z = z;
			this.levels = levels;
			this.params = params;
			this.procip = imp.getStack().getProcessor( imp.getStackIndex(params.deadpanC, z, 1) ).duplicate();
		}
		
		@Override
		public ArrayList<Roi> call() {
			ArrayList<Roi> rois = new ArrayList<Roi>();
			try{
				ThresholdToSelection tts = new ThresholdToSelection();
				
				int W = imp.getWidth();
				int H = imp.getHeight();
				double pixW = imp.getCalibration().pixelWidth;
				ImageStatistics stackStats = new StackStatistics(imp);
				
				ByteProcessor outip = new ByteProcessor(W, H);
				outip.setColor(Color.WHITE);	//to add objects
				procip.setColor(Color.BLACK);	//to remove objects
				for(int m=0; m<levels.length; m++){		//means in increasing order
				//for(int m=levels.length-1; m>=0; m--){		//means in decreasing order
	
					ImageProcessor mask = new ByteProcessor(W, H);
					final double thresh = levels[m];
					IntStream.range(0,W*H).forEach( i -> mask.set(i, (procip.getf(i)>=thresh?255:0)	));
					
					if(params.watershed){
						double tolerance = 0.3;
						FloatProcessor edm = new EDM().makeFloatEDM(mask, 0, false);
						ImageProcessor seg = new MaximumFinder().findMaxima(edm, tolerance, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, false, true);
						mask.copyBits(seg, 0, 0, Blitter.AND);
					}
	
					if(mask.getStatistics().mean==0) continue;
					Roi roi = null;
					try{
						mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE);
						roi = tts.convert(mask);
					}catch(ArrayIndexOutOfBoundsException oob){continue;}	//ignore ArrayIndexOutOfBoundsException from ThresholdToSelection when image is empty
					if(roi==null) continue;
					
					if(roi!=null){
						Roi[] split = new ShapeRoi(roi).getRois();
						
						for(int ri=0;ri<split.length;ri++){
							Roi r = split[ri];
							
							procip.setRoi(r);
							ImageStatistics procStats = ImageStatistics.getStatistics(procip, ImageStatistics.AREA+ImageStatistics.MEAN, imp.getCalibration());
							outip.setRoi(r);
							if(outip.getStatistics().mean > 0){ //already added
								continue;
							}
							
							double perim = r.getLength() * pixW;
							double circ = 4*Math.PI*(procStats.area/(perim*perim));
							//double intensityF = procStats.max / stackStats.max;	//TODO: use procStats.mean instead?
							double intensityF = procStats.mean / stackStats.max;
							
							//filter
							if( procStats.area>=params.minA && procStats.area<=params.maxA && intensityF>=params.threshold && circ>=params.minCirc ){
								r.setPosition(-1,z,1);
								rois.add(r);
								outip.fill(r); //add to the binary mask (white)
								procip.fill(r); //remove from the thresholding image (black)
							}
						}
					}
				}
			}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
			return rois;
		}
		
	}
	
	/**Gets all Rois from a stack using a new <code>SliceSegmenter</code> running on a separate Thread for each slice.
	 * 
	 * @param imp	a processed z-stack with one channel
	 * @param levels	hierarchy of intensity levels
	 * @param params	parameters to be used, passed to <code>SliceSegmenter</code>
	 * 
	 * @return	a combined ArrayList containing all Rois from each <code>SliceSegmenter</code>
	 */
	private ArrayList<Roi> getRois(ImagePlus imp, double[] levels, DeadParams params) {
		ArrayList<Roi> rois = new ArrayList<Roi>();
		try{
			status.log("Hierarchical mapping...");
			int nThreads = Runtime.getRuntime().availableProcessors()-1;
			executor = Executors.newFixedThreadPool(nThreads);
			futures = new HashSet<Future<ArrayList<Roi>>>();
			for(int z=params.startZ;z<=params.endZ;z++){
				SliceSegmenter seg = new SliceSegmenter(imp, z, levels, params);
				Future<ArrayList<Roi>> future = executor.submit(seg);
				futures.add(future);
			}
			for(Future<ArrayList<Roi>> future : futures){
				try {
					if(future.isCancelled()) continue;
					ArrayList<Roi> list = future.get();
					rois.addAll( list );
				} catch (CancellationException ce) {
					System.out.print("SliceSegmenter did not cancel normally: "+ce.toString());
				}catch (Exception e) {
					System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));
				}
			}
		}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
		return rois;
	}
	
	/** Joins Rois into nuclei.
	 * 
	 * @param imp	a processed z-stack with one channel
	 * @param roiList	the Rois to be joined
	 * @param params	parameters to be used, passed to <code>SliceSegmenter</code>
	 * 
	 * @return	an ArrayList of <code>Nucleus3D</code>
	 */
	private ArrayList<Nucleus3D> joinRois(ImagePlus imp, ArrayList<Roi> roiList, DeadParams params) {
		if(cancel) return null;
		
		ArrayList<Nucleus3D> nuclei = new ArrayList<Nucleus3D>();
		
		int[] index = new int[roiList.size()];
		IntStream.range(0,roiList.size()).forEach(i -> index[i] = i);
		
		int n = roiList.size();
		int nPairs = (int) ((n*(n+1))/2f) - n;	//the number of pairs to be tested for merging is the Gaussian sum of n minus n
		int count = 0;
		int step = (int) Math.max(nPairs/10f, 10);
		
		//merge indices for Rois to be joined
		for(int i=0;i<roiList.size();i++){
			if(cancel) return null;
			Roi roii = roiList.get(i);
			int zi = roii.getZPosition();
			double iA = roii.getStatistics().area;
			for(int j=i+1;j<roiList.size();j++){
				count++;
				if(count%step==0){	// 1/10 steps to output progress in linear time
					String percent = IJ.d2s(((count/(float)nPairs)*100), 1)+"%";
					status.log("3D Reconstruction "+percent);
				}
				Roi roij = roiList.get(j);
				int zj = roij.getZPosition();
				double jA = roii.getStatistics().area;
				if(Math.abs(zi-zj)<=params.joinZ){
					Roi overlap = new ShapeRoi(roii).and(new ShapeRoi(roij));
					if(overlap.getLength()==0) continue;
					double overlapA = overlap.getStatistics().area;
					double overlapF = (2*overlapA) / (iA+jA);
					if(overlapF >= params.joinOverlap && index[j] != index[i]){
						int oldi = index[j];
						int newi = index[i];
						for(int r=0;r<index.length;r++){
							if(index[r]==oldi){
								index[r] = newi;
							}
						}
					}
				}
			}
		}
		if(cancel) return null;
		
		//create Nucleus3Ds from Rois with the same index
		HashSet<Integer> iset = new HashSet<Integer>();
		for(int i : index){
			iset.add(i);
		}
		for(Integer i : iset){
			Nucleus3D nuc = new Nucleus3D(imp, i);
			for(int a=0;a<roiList.size();a++){
				if(index[a]==i){
					nuc.add(roiList.get(a));
				}
			}
			if( params.endZ-params.startZ>0 && nuc.size()==1 ) continue;	//ignore if not a single slice preview and only one Roi
			nuclei.add(nuc);
		}
		
		// filter out deadpan labelled trachea - remove whole Nucleus3D if any component Roi is a trachea fragment by density variation scoring
		if(params.minFlatness>0){
			ArrayList<Nucleus3D> filtered = new ArrayList<Nucleus3D>();
			status.log("Trachea Filtering...");
			for(Nucleus3D nuc : nuclei){
				if(cancel) return null;
				boolean good = true;
				for(Roi roi : nuc.getRois()){
					ImageProcessor ip = imp.getStack().getProcessor( imp.getStackIndex(params.deadpanC, roi.getZPosition(), 1) );
					ip.setRoi(roi);
					ImageStatistics stats = ImageStatistics.getStatistics(ip, ImageStatistics.MIN_MAX+ImageStatistics.MEAN+ImageStatistics.STD_DEV, imp.getCalibration());
					double dv = (stats.stdDev*stats.stdDev) / (stats.mean*stats.area);	// variance / integrated density ("Density Variation")
					double flatness = 1d / Math.sqrt(dv);	//reciprocal of root density variation
					if( flatness < params.minFlatness ){	//exclude nucleus if any Roi DV is too high
						good = false;
						break;
					}
				}
				if(good){
					filtered.add(nuc);
				}
			}
			nuclei = filtered;
		}
		
		//sort nuclei by coordinates to make indices deterministic - multithreading makes indices non-deterministic due to order of Roi List
		nuclei.sort(new Comparator<Nucleus3D>(){
			public int compare(Nucleus3D nuc0, Nucleus3D nuc1){
				double[] c0 = nuc0.getCentroid();
				double[] c1 = nuc1.getCentroid();
				int d = Double.compare(c0[2], c1[2]);	//compare by z coordinate
				if(d==0) d = Double.compare(c0[1]+c0[0], c1[1]+c0[0]);	//if z is the same, compare by y+x
				return d;
			}
		});
		int ni = 0;
		for(Nucleus3D nuc : nuclei){
			nuc.setIndex(ni++);	//set sequential indices 
		}
		status.log("Got "+nuclei.size()+" nuclei");	
		return nuclei;
	}
	
}
