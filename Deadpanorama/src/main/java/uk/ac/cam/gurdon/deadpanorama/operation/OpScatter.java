package uk.ac.cam.gurdon.deadpanorama.operation;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.Collectors;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.DefaultXYDataset;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.ResultsTable;
import ij.process.StackStatistics;
import uk.ac.cam.gurdon.deadpanorama.DeadParams;
import uk.ac.cam.gurdon.deadpanorama.Nucleus3D;

public class OpScatter extends Operation {
	
	static{
		TYPE = Operation.Type.SCATTER;
		TIP = "Scatter plot of mean intensity values for each nucleus in two channels.<br>Applies OPTICS to identify clusters in the scatter plot.";
		NCHANNELS = 2;
	}
	
	public OpScatter(){
		super();
	}
	
	@Override
	public void run(ImagePlus imp, ArrayList<Nucleus3D> nuclei, DeadParams params, ResultsTable rt) {
		try{
			int channelA = channelACombo.getSelectedIndex();
			int channelB = channelBCombo.getSelectedIndex();
			boolean outerA = outerTickA.isSelected();
			boolean outerB = outerTickB.isSelected();
			if(channelA==0||channelB==0||channelA>imp.getNChannels()||channelB>imp.getNChannels()) return;
			int bd = imp.getBitDepth();
			if(bd==16&&new StackStatistics(imp).max<=4095) bd = 12;
			double max = Math.pow(2, bd);
			
			ArrayList<Detection> points = new ArrayList<Detection>();
			for(int i=0;i<nuclei.size();i++){
				if(cancel) return;
				Nucleus3D nuc = nuclei.get(i);
				float Amean = nuc.getMean(channelA, outerA);
				float Bmean = nuc.getMean(channelB, outerB);
				rt.setValue("C"+channelA+" Norm. "+(outerA?"Surrounding":""), nuc.index, Amean/max);
				rt.setValue("C"+channelB+" Norm. "+(outerB?"Surrounding":""), nuc.index, Bmean/max);
				
				points.add( new Detection(Amean/max, Bmean/max) );
			}
			DefaultXYDataset dataset = new DefaultXYDataset();
	
			double minX = Double.POSITIVE_INFINITY; double minY = Double.POSITIVE_INFINITY;
			double maxX = Double.NEGATIVE_INFINITY; double maxY = Double.NEGATIVE_INFINITY;
			for(Detection d:points){
				minX = Math.min(minX,d.x); minY = Math.min(minY,d.y);
				maxX = Math.max(maxX,d.x); maxY = Math.max(maxY,d.y);
			}
			double maxR = 0.1;	//maximum linking distance for normalised intensity values
			
			ArrayList<ArrayList<Detection>> clusters = OPTICS(points, maxR);
			ArrayList<Detection> remaining = (ArrayList<Detection>) points.stream().collect(Collectors.toList());
			for(ArrayList<Detection> cluster : clusters){
				int n = cluster.size();
				double[][] series = new double[2][n];
				int i = 0;
				double meanX = 0d;
				double meanY = 0d;
				for(Detection d:cluster){
					if(remaining.contains(d)){
						series[0][i] = d.x;
						series[1][i] = d.y;
						meanX += d.x;
						meanY += d.y;
						remaining.remove(d);
						i++;
					}
				}
				meanX /= (float)n;
				meanY /= (float)n;
				dataset.addSeries(IJ.d2s(meanX,2)+","+IJ.d2s(meanY,2)+" ("+n+")", series);
			}
			double[][] remainingSeries = new double[2][remaining.size()];
			for(int i=0;i<remaining.size();i++){	//series for points not assigned to a cluster
				remainingSeries[0][i] = remaining.get(i).x;
				remainingSeries[1][i] = remaining.get(i).y;
			}
			String name = clusters.size()>0?"outliers":"all";	//either there were no separate clusters or the remaining points are outliers
			dataset.addSeries(name+" ("+remaining.size()+")", remainingSeries);
			
			makeChart(dataset, "C"+channelA+(outerA?" Surrounding":" Nuclei")+" vs C"+channelB+(outerB?" Surrounding":" Nuclei"), "C"+channelA, "C"+channelB);
		}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	}

	private static void makeChart( DefaultXYDataset dataset, String title, String X, String Y) throws Exception{
		JFreeChart chart = ChartFactory.createScatterPlot(title, X, Y, dataset, PlotOrientation.VERTICAL, true, true, true);
		XYPlot plot = chart.getXYPlot();
        plot.setDomainCrosshairVisible(false);
        plot.setRangeCrosshairVisible(false);
        plot.setBackgroundPaint(Color.BLACK);
    	plot.setDomainGridlinePaint(Color.GRAY);
    	plot.setRangeGridlinePaint(Color.GRAY);
    	
    	NumberAxis domain = (NumberAxis) plot.getDomainAxis();	//x-axis
        domain.setRange(0.00, 1.00);
        domain.setTickUnit(new NumberTickUnit(0.1));
        domain.setVerticalTickLabels(false);
        NumberAxis range = (NumberAxis) plot.getRangeAxis();	//y-axis
        range.setRange(0.0, 1.0);
        range.setTickUnit(new NumberTickUnit(0.1));
    	
        Shape shape = new Ellipse2D.Float(-1f, -1f, 2f, 2f);
        XYItemRenderer render = plot.getRenderer();
        for(int i=0;i<plot.getSeriesCount();i++){
        	float h = i/(float)plot.getSeriesCount();
    		float s = 1.0f;
    		float v = 1.0f;
    		Color hsv = Color.getHSBColor(h,s,v);
    		float[] rgb = hsv.getRGBColorComponents(new float[3]);
    		Color colour = new Color(rgb[0], rgb[1], rgb[2]);
	    	render.setSeriesPaint(i, colour);
	    	render.setSeriesShape(i, shape );
        }

        ChartFrame frame = new ChartFrame(title,chart);
        frame.pack();
		frame.setSize( new Dimension(800, 800) );
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}
	
	
	
	
	
	
	
	
	
	
	
	private class Detection{
		double x, y;
		
		private boolean processed;
		private Double coreDistance;
		private Double reachabilityDistance;
		
		
		public Detection(double x, double y){
			this.x = x;
			this.y = y;
			this.processed = false;
		}
		
		public double distance2D(Detection other){
			double dist = Math.sqrt( sqdistance2D(other) );
			return dist;
		}
		
		public double sqdistance2D(Detection other){
			return ((x-other.x)*(x-other.x))+((y-other.y)*(y-other.y));
		}
		
		public Double getCd(){
			return coreDistance;
		}
		public void setCd(double value){
			coreDistance = value;
		}
		
		public Double getRd(){
			return reachabilityDistance;
		}
		public void setRd(double value){
			reachabilityDistance = value;
		}
		
		public boolean isProcessed(){
			return processed;
		}
		public void setProcessed(boolean value){
			processed = value;
		}
		
		@Override
		public boolean equals(Object obj){
			return this==obj;
		}
		@Override
		public int hashCode(){
			return (int) (31 + Math.round(x+y));
		}

	}
	
	
	
	
	int minN = 3;
	ArrayList<Detection> unprocessed;
	ArrayList<Detection> ordered;
	
	private ArrayList<ArrayList<Detection>> OPTICS(ArrayList<Detection> points, double maxR){
		
		unprocessed = (ArrayList<Detection>) points.stream().collect(Collectors.toList());
		ordered = new ArrayList<Detection>();
		
		while(unprocessed.size()>0){
			
			Detection p = unprocessed.get(0);
			processed( p );
			ArrayList<Detection> neighbours = neighbours(p, points, maxR);
			if( coreDistance(p, neighbours) != null){
				ArrayList<Detection> seeds = new ArrayList<Detection>();
				update( neighbours, p, seeds );
				while(seeds.size()>0){
					seeds.sort( Comparator.comparing(d->d.getRd()) );
					Detection n = seeds.get(0);	//seed with smallest reachability distance
					seeds.remove(n);
					processed(n);
					ArrayList<Detection> neighn = neighbours(n, points, maxR);
					if (coreDistance(n, neighn)!=null){
						update(neighn, n, seeds);
					}
				}
			}
		}
		
		ArrayList<Integer> separators = clusterSeparators(ordered, maxR);

        ArrayList<ArrayList<Detection>> clusters = new ArrayList<ArrayList<Detection>>();
        for(int i=0;i<separators.size()-2;i++){
            int start = separators.get(i);
            int end = separators.get(i + 1);
            if(end - start >= minN){
            	ArrayList<Detection> clus = new ArrayList<Detection>();
            	for(int a=start;a<end;a++){
            		clus.add(ordered.get(a));
            	}
            	clusters.add(clus);
            }
        }
        
		return clusters;
	}
	
	private void update(ArrayList<Detection> neighbours, Detection p, ArrayList<Detection> seeds){
		for(Detection neigh:neighbours){
			if(neigh.isProcessed()) continue;
			double reachability = Math.max( p.getCd(), p.distance2D(neigh) );
			if(neigh.getRd()==null){
				neigh.setRd(reachability);
				seeds.add( neigh );
			}
			else if(reachability < neigh.getRd()){
				neigh.setRd(reachability);
			}
		}
	}
	
	private Double coreDistance(Detection p, ArrayList<Detection> neighbours) {
		if(p.getCd()!=null) return p.getCd();
		
		if(neighbours.size() >= minN-1){
			ArrayList<Double> distances = (ArrayList<Double>) neighbours.stream().map( d->d.distance2D(p) ).sorted().collect( Collectors.toList() );
			p.setCd( distances.get(minN-2) );
		}
		
		return p.getCd();
	}

	private ArrayList<Detection> neighbours(Detection p, ArrayList<Detection> points, double maxR) {
		ArrayList<Detection> neigh = (ArrayList<Detection>) points.stream().filter( d-> p.distance2D(d)<maxR  ).collect( Collectors.toList() );
		return neigh;
	}

	private void processed(Detection point) {
		point.setProcessed(true);
		unprocessed.remove(point);
		ordered.add(point);
	}

	private ArrayList<Integer> clusterSeparators(ArrayList<Detection> ordered, double maxR){
		ArrayList<Integer> separators = new ArrayList<Integer>();
		for(int i=0;i<ordered.size();i++){
			double rd = Double.POSITIVE_INFINITY;
			if(ordered.get(i).getRd()!=null){
				rd = ordered.get(i).getRd();
			}
			if( rd > maxR){ //separate the clusters based on maximum linking radius
			    separators.add(i);
			}
		}
	    separators.add(ordered.size());
	    return separators;
	}

}
