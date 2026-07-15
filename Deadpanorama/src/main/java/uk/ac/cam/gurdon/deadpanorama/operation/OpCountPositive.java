package uk.ac.cam.gurdon.deadpanorama.operation;

import java.awt.BasicStroke;
import java.awt.Color;
import java.util.ArrayList;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

import ij.ImagePlus;
import ij.measure.ResultsTable;
import ij.process.AutoThresholder;
import ij.process.StackStatistics;
import uk.ac.cam.gurdon.deadpanorama.DeadParams;
import uk.ac.cam.gurdon.deadpanorama.Nucleus3D;

public class OpCountPositive extends Operation {
	
	static{
		TYPE = Operation.Type.COUNTPOSITIVE;
		TIP = "Count positive cells in the selected channel.";
		NCHANNELS = 1;
	}
	
	public OpCountPositive(){
		super();
	}
	
	@Override
	public void run(ImagePlus imp, ArrayList<Nucleus3D> nuclei, DeadParams params, ResultsTable rt) {
		int channel = channelACombo.getSelectedIndex();
		boolean outer = outerTickA.isSelected();
		String col = outer?" Surrounding":" Nucleus";
		String cols = outer?" Surrounding":" Nuclei";
		if(channel==0||channel>imp.getNChannels()) return;
		double[] means = new double[nuclei.size()];
		float min = Float.MAX_VALUE;
		float max = Float.MIN_VALUE;
		int i = 0;
		for(Nucleus3D nuc : nuclei){
			if(cancel) return;
			float mean = nuc.getMean(channel, outer);
			
			means[i] = mean;
			min = Math.min(min, mean);
			max = Math.max(max, mean);
			i++;
		}
		
		HistogramDataset dataset = new HistogramDataset();
		dataset.setType(HistogramType.FREQUENCY);
		int nbins = 256;
		dataset.addSeries("C"+channel+col, means, nbins);
		int[] histogram = new int[nbins];
		for(int b=0;b<nbins;b++){
			int y = dataset.getY(0, b).intValue();
			histogram[b] = y;
		}
		
		int bd = imp.getBitDepth();
		if(bd==16&&new StackStatistics(imp).max<=4095) bd = 12;
		double bdMax = Math.pow(2, bd);
		int positiveThresh = (int) (bdMax);	//no positive objects
		String ptStr = "";
		
		if(max-min>bdMax/4f){	//only set a threshold to detect positive objects if the intensity range is large enough
			int thresh = new AutoThresholder().getThreshold(AutoThresholder.Method.MaxEntropy, histogram);
			positiveThresh = dataset.getX(0, thresh).intValue();
			ptStr = " (>="+String.valueOf(positiveThresh)+")";
		}
		
		
		int count = 0;
		for(Nucleus3D nuc : nuclei){
			if(cancel) return;
			double mean = nuc.getMean(channel, outer);
			rt.setValue("C"+channel+col+" Mean", nuc.index, mean);	//will overwrite same value if OpStatistics has already been run on this channel
			String pos = "-";
			if(mean>=positiveThresh){
				count++;
				pos = "+";
			}
			rt.setValue("C"+channel+col+" +/-ve"+ptStr, nuc.index, pos );
		}
		
		JFreeChart chart = ChartFactory.createHistogram(count+" C"+channel+cols+" positive out of "+nuclei.size(), "C"+channel+col+" Mean Intensity", "n", dataset, PlotOrientation.VERTICAL, false, true, false);
		XYPlot plot = chart.getXYPlot();
		XYBarRenderer renderer = (XYBarRenderer) plot.getRenderer();
		renderer.setSeriesPaint(0, Color.BLUE);
		XYLineAnnotation threshAnn = new XYLineAnnotation(positiveThresh, 0, positiveThresh, nuclei.size(), new BasicStroke(1f), Color.RED);
		renderer.addAnnotation(threshAnn);
		StandardXYBarPainter painter = new StandardXYBarPainter();
		renderer.setBarPainter(painter);
		
		NumberAxis range = (NumberAxis) plot.getRangeAxis();	//y-axis
	    range.setTickUnit(new NumberTickUnit(1));
		
		ChartFrame frame = new ChartFrame("C"+channel+col+" Intensity Frequency", chart);
		frame.setSize(800, 800);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		
	}
	
}

