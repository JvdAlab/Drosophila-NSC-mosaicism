package uk.ac.cam.gurdon.deadpanorama.operation;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.JLabel;
import javax.swing.JPanel;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.plugin.RGBStackMerge;
import ij.process.ByteProcessor;
import ij.process.ImageConverter;
import ij.process.LUT;
import uk.ac.cam.gurdon.deadpanorama.DeadParams;
import uk.ac.cam.gurdon.deadpanorama.Nucleus3D;
import uk.ac.cam.gurdon.deadpanorama.Projectile;

public class OpRender3D extends Operation {

	static{
		TYPE = Operation.Type.RENDER3D;
		TIP = 	"Render in 3D.";
		NCHANNELS = 1;
	}
	
	public OpRender3D() {
		panel = new JPanel();
		panel.add( new JLabel("3D Render Deadpan nuclei with ") );
		
		panel.add( new JLabel("C") );
		channelACombo = new CComboBox( );
		channelACombo.setSelectedIndex(2);
		panel.add( channelACombo );

		panel.setToolTipText(TIPSTART+TIP);
	}

	@Override
	public void run(ImagePlus imp, ArrayList<Nucleus3D> nuclei, DeadParams params, ResultsTable rt) {
		try{
			int channel = channelACombo.getSelectedIndex();
			if(channel==0||channel>imp.getNChannels()) return;
			
			ImageStack impStack = imp.getStack();
			int W = imp.getWidth();
			int H = imp.getHeight();
			ImageStack maskStack = new ImageStack(W,H);
			ImageStack signalStack = new ImageStack(W,H);
			for(int z=0;z<imp.getNSlices();z++){
				maskStack.addSlice( new ByteProcessor(W,H) );
				signalStack.addSlice( impStack.getProcessor( imp.getStackIndex(channel, z+1, 1) ).convertToByteProcessor() );
			}
			
			int nColours = Math.min(255, nuclei.size());
			
			byte[] r = new byte[256];
			byte[] g = new byte[256];
			byte[] b = new byte[256];
			for(int i=1;i<=nColours;i++){	//0 is 0,0,0
				
				float h = i/(float)nColours;
				if(i%2==0) h = 1-h;	//alternate colour wheel directions
				float s = 1.0f;
				float v = 1.0f;
				float alpha = 1f;
				Color hsv = Color.getHSBColor(h,s,v);
				float[] rgb = hsv.getRGBColorComponents(new float[3]);
				Color colour = new Color(rgb[0], rgb[1], rgb[2], alpha);
				
				r[i] = (byte) colour.getRed();
				g[i] = (byte) colour.getGreen();
				b[i] = (byte) colour.getBlue();
			}
			LUT colourful = new LUT(r,g,b);
			
			for(Nucleus3D nuc : nuclei){
				if(cancel) return;
				for(Roi roi : nuc.getRois()){
					int z = roi.getZPosition();
					ByteProcessor bp = (ByteProcessor) maskStack.getProcessor(z);
					bp.setColor(nuc.index<255?nuc.index:255%nuc.index);
					bp.fill(roi);
				}
				
			}
			
			ImagePlus signalImp = new ImagePlus("",signalStack);
			
			LUT[] luts = imp.getLuts();
			signalImp.setLut(  luts[channel-1] );
			
			signalImp.setDisplayRange(0,255);
			ImagePlus maskImp = new ImagePlus("", maskStack);
			maskImp.setLut( colourful );
			ImagePlus hyper = RGBStackMerge.mergeChannels( new ImagePlus[]{signalImp, maskImp}, false );
			hyper.setDisplayMode(IJ.COMPOSITE);
			hyper.setCalibration(imp.getCalibration());
			new ImageConverter(hyper).convertToRGB();
			hyper.setTitle(imp.getTitle()+"_Deadpan+C"+channel);

			new Projectile().run(hyper, Projectile.brightestPoint, Projectile.yAxis, 0, 360, 10);
			
		}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	}
	
}
