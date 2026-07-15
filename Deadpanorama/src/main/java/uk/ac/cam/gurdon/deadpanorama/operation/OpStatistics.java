package uk.ac.cam.gurdon.deadpanorama.operation;

import java.util.ArrayList;
import java.util.Arrays;

import ij.ImagePlus;
import ij.measure.ResultsTable;
import uk.ac.cam.gurdon.deadpanorama.DeadParams;
import uk.ac.cam.gurdon.deadpanorama.Nucleus3D;

public class OpStatistics extends Operation {
	
	static{
		TYPE = Operation.Type.STATISTICS;
		TIP = "Measure intensity mean, standard deviation, min and max for each nucleus.";
		NCHANNELS = 1;
	}
	
	public OpStatistics(){
		super();
	}
	
	@Override
	public void run(ImagePlus imp, ArrayList<Nucleus3D> nuclei, DeadParams params, ResultsTable rt) {
		int channel = channelACombo.getSelectedIndex();
		boolean outer = outerTickA.isSelected();
		if(channel==0||channel>imp.getNChannels()) return;
		for(Nucleus3D nuc : nuclei){
			if(cancel) return;
			Float[] values = nuc.getValues(channel, outer);
			String col = outer?" Surrounding":" Nucleus";
			
			double mean = Arrays.stream(values).mapToDouble(f->f).sum() / (float)values.length;
			rt.setValue("C"+channel+col+" Mean", nuc.index, mean);
			double ssd = Arrays.stream(values).mapToDouble(f-> (f-mean)*(f-mean) ).sum();
			double stdDev = Math.sqrt(ssd/(float)values.length);
			rt.setValue("C"+channel+col+" StdDev", nuc.index, stdDev);
			
			double min = Arrays.stream(values).mapToDouble(f->f).min().getAsDouble();
			double max = Arrays.stream(values).mapToDouble(f->f).max().getAsDouble();
			rt.setValue("C"+channel+col+" Min", nuc.index, min);
			rt.setValue("C"+channel+col+" Max", nuc.index, max);
			
		}
	}
	
}
