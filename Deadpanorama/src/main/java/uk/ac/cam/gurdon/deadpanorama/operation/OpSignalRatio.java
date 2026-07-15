package uk.ac.cam.gurdon.deadpanorama.operation;

import java.util.ArrayList;

import ij.ImagePlus;
import ij.measure.ResultsTable;
import uk.ac.cam.gurdon.deadpanorama.DeadParams;
import uk.ac.cam.gurdon.deadpanorama.Nucleus3D;


public class OpSignalRatio extends Operation {

	static{
		TYPE = Operation.Type.SIGNALRATIO;
		NCHANNELS = 2;
		TIP = "Calculate the ratio of mean signal intensities between two channels in each nucleus.";
	}
	
	public OpSignalRatio(){
		super();
	}

	@Override
	public void run(ImagePlus imp, ArrayList<Nucleus3D> nuclei, DeadParams params, ResultsTable rt) {
		int channelA = channelACombo.getSelectedIndex();
		int channelB = channelBCombo.getSelectedIndex();
		boolean outerA = outerTickA.isSelected();
		boolean outerB = outerTickA.isSelected();
		String colA = outerA?" Surrounding":" Nucleus";
		String colB = outerB?" Surrounding":" Nucleus";
		if(channelA==0||channelB==0||channelA>imp.getNChannels()||channelB>imp.getNChannels()) return;
		for(Nucleus3D nuc : nuclei){
			if(cancel) return;
			
			double meanA = nuc.getMean(channelA, outerA);
			double meanB = nuc.getMean(channelB, outerB);
			double ratio = meanA / meanB;
			
			rt.setValue("C"+channelA+colA+" Mean", nuc.index, meanA);
			rt.setValue("C"+channelB+colB+" Mean", nuc.index, meanB);
			rt.setValue("C"+channelA+colA+":C"+channelB+colB+" Mean Ratio", nuc.index, ratio);
			
		}
	}

}
