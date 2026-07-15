package uk.ac.cam.gurdon.deadpanorama.operation;

import java.util.ArrayList;
import java.util.Arrays;

import ij.ImagePlus;
import ij.measure.ResultsTable;
import uk.ac.cam.gurdon.deadpanorama.DeadParams;
import uk.ac.cam.gurdon.deadpanorama.Nucleus3D;

public class OpCorrelation extends Operation {
	
	static{
		TYPE = Operation.Type.CORRELATION;
		TIP = "Calculate Pearson's Correlation Coefficient, Manders' Overlap Coefficient and <br>Li's Intensity Correlation Quotient for two channels in each nucleus.";
		NCHANNELS = 2;
	}
	
	public OpCorrelation(){
		super();
		outerTickA.setVisible(false);	//both channels have to use inner or outer values, use outerTickB for both
	}
	
	@Override
	public void run(ImagePlus imp, ArrayList<Nucleus3D> nuclei, DeadParams params, ResultsTable rt) {
		int channelA = channelACombo.getSelectedIndex();
		int channelB = channelBCombo.getSelectedIndex();
		boolean outer = outerTickB.isSelected();
		if(channelA==0||channelB==0||channelA>imp.getNChannels()||channelB>imp.getNChannels()) return;
		for(Nucleus3D nuc : nuclei){
			if(cancel) return;
			
			Float[] Avalues = nuc.getValues(channelA, outer);
			Float[] Bvalues = nuc.getValues(channelB, outer);
			
			float nf = (float)Avalues.length;
			
			double Amean = Arrays.stream(Avalues).mapToDouble(f->f).sum() / nf;
			double Bmean = Arrays.stream(Bvalues).mapToDouble(f->f).sum() / nf;

			double numPearson = 0.0;
			double denom1Pearson = 0.0;
			double denom2Pearson = 0.0;
			double numManders = 0.0;
			double denom1Manders = 0.0;
			double denom2Manders = 0.0;
			double ICQ = 0.0;
			for(int i=0;i<nf;i++){
				if(cancel) return;
				numPearson += (Avalues[i]-Amean)*(Bvalues[i]-Bmean);
				denom1Pearson += (Avalues[i]-Amean)*(Avalues[i]-Amean);
				denom2Pearson += (Bvalues[i]-Bmean)*(Bvalues[i]-Bmean);

				numManders += Avalues[i] * Bvalues[i];
				denom1Manders += Avalues[i] * Avalues[i];
				denom2Manders += Bvalues[i] * Bvalues[i];

				if( (Avalues[i]-Amean) * (Bvalues[i]-Bmean) > 0){
					ICQ += 1;
				}
			}

			double pearson = numPearson/Math.sqrt(denom1Pearson*denom2Pearson);
			double manders = numManders/Math.sqrt(denom1Manders*denom2Manders);
			ICQ = ICQ/nf - 0.5;
			
			String col = outer?" Surrounding":" Nucleus";
			rt.setValue("Pearson Correlation C"+channelA+" vs C"+channelB+col, nuc.index, pearson);
			rt.setValue("Manders Overlap C"+channelA+" vs C"+channelB+col, nuc.index, manders);
			rt.setValue("ICQ C"+channelA+" vs C"+channelB+col, nuc.index, ICQ);
		}
	}
	
}
