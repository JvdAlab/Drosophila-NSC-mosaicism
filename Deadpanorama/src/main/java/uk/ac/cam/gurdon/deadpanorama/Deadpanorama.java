package uk.ac.cam.gurdon.deadpanorama;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;

import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.Roi;
import uk.ac.cam.gurdon.deadpanorama.operation.OperationPanel;


@Plugin(type = Command.class, menuPath = "Plugins>Deadpanorama")
public class Deadpanorama implements Command{

	private Segmenter segmenter;
	private OperationPanel opPanel;
	
	private ArrayList<Nucleus3D> nuclei;
	
	public void run() {
		new DPgui(this);
	}

	public ArrayList<Nucleus3D> run(ImagePlus imp, DeadParams params, OperationPanel opp, StatusPanel status){
		this.opPanel = opp;
		nuclei = null;
		try{
			Roi userRoi = imp.getRoi();
			imp.killRoi();

			segmenter = new Segmenter(status);
			nuclei = segmenter.getNuclei(imp, params, userRoi);
			if(nuclei==null) return null;
			
			if(opPanel!=null){	//null passed for previews
				opPanel.execute(imp, nuclei, params, status);
			}
			
			Overlay ol = new Overlay();
			if(userRoi!=null){
				Roi userRoiShow = (Roi) userRoi.clone();
				userRoiShow.setPosition(-1, -1, -1);
				userRoiShow.setStrokeWidth(3);
				userRoiShow.setStrokeColor(Color.MAGENTA);
				ol.add(userRoiShow);
			}
			for(Nucleus3D nuc : nuclei){
				nuc.addToOverlay(ol);
			}
			imp.setOverlay(ol);
			imp.setRoi(userRoi);
			
status.log("Finished "+imp.getTitle()+"\n");
		}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
		return nuclei;
	}
	
	public ArrayList<Nucleus3D> getNuclei(){
		return nuclei;
	}
	
	public void kill(){
		if(segmenter!=null) segmenter.cancel();
		if(opPanel!=null) opPanel.cancel();
	}
	
	public static void main(String[] arg){
		
		ImageJ.main(arg);
		/*//ImagePlus img = new ImagePlus("E:\\Jelle\\Deadpan\\19-01-15 WOR58 WorG4-GFP-G80ts x BlwTRIP NaOH wander at29 Dpn-GFP-Imp-pH3.lif - Mark_and_Find 001 Position001.tif");
		//ImagePlus img = new ImagePlus("E:\\Jelle\\Deadpan\\18-04-09 WOR51 WorG4-GFP-G80ts x Rbf280 at29 wander Dpn-GFP-Imp-PH3.lif - Mark_and_Find 001 Position001.tif");
		//ImagePlus img = new ImagePlus("E:\\Jelle\\Deadpan\\18-04-09 WOR51 WorG4-GFP-G80ts x Rbf280 at29 wander Dpn-GFP-Imp-PH3.lif - Mark_and_Find 001 Position006.tif");
		ImagePlus img = new ImagePlus("E:\\Jelle\\smalltest_20_09_21_hsflp_F20_DpnFRTSTOPLexA_LexAoptdTom_LeoxOmCherry_Wandering_25_1hHS_Dpn_RFP_VNC1.tif");
		
		final ImagePlus image = HyperStackConverter.toHyperStack(img, img.getNChannels(), img.getNSlices(), img.getNFrames());
		image.setDisplayMode(IJ.COLOR);
		image.show();*/
		
		new Deadpanorama().run();
	}
	
	
}
