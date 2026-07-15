package uk.ac.cam.gurdon.deadpanorama.operation;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.JPanel;

import ij.ImagePlus;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import uk.ac.cam.gurdon.deadpanorama.DeadParams;
import uk.ac.cam.gurdon.deadpanorama.Nucleus3D;
import uk.ac.cam.gurdon.deadpanorama.StatusPanel;


public class OperationPanel extends JPanel{
	private static final long serialVersionUID = 3416215893929047100L;
	
	private Path path = Paths.get(System.getProperty("user.home")+System.getProperty("file.separator")+"Deadpanorama.cfg");
	
	private ArrayList<Operation> opList;
	private boolean cancel = false;
	private ResultsTable rt;
	private String resultsTitle;
	
	public OperationPanel(){
		loadOpList();
		setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
	}
	
	/** Create and add a new Operation with default parameters
	 */
	public void add(Operation.Type type){
		Operation op = Operation.getOperation(type);
		add(op.getPanel());
		opList.add(op);
	}
	
	/** Add an Operation with parameters possibly already set
	 */
	public void add(Operation op){
		add(op.getPanel());
		opList.add(op);
	}
	
	public void remove(){
		if(opList.size()>0){
			remove( opList.get(opList.size()-1).getPanel() );
			opList.remove(opList.size()-1);
		}
	}
	
	public void execute(ImagePlus imp, ArrayList<Nucleus3D> nuclei, DeadParams params, StatusPanel status){
		try{
			Calibration cal = imp.getCalibration();
			String unitV = " ("+cal.getUnit()+"\u00b3)";
			rt = new ResultsTable();
			rt.showRowNumbers(false);
			int row = 0;
			status.log("Measuring...");
			for(Nucleus3D nuc : nuclei){	//create table and get spatial statistics output for everything
				if(cancel) return;
				rt.setValue("Deadpan Nucleus", row, nuc.index);
				double[] centroid = nuc.getCentroid();
				rt.setValue("X", row, centroid[0]*cal.pixelWidth);
				rt.setValue("Y", row, centroid[1]*cal.pixelHeight);
				rt.setValue("Z", row, centroid[2]*cal.pixelDepth);
				double volume = nuc.getVolume();
				rt.setValue("Volume"+unitV, row, volume*cal.pixelWidth*cal.pixelHeight*cal.pixelDepth);
				row++;
			}
			
			for(Operation op : opList){
				if(cancel) return;
				status.log("Operation : "+op.toString());
				op.run(imp, nuclei, params, rt);
			}
			
			if(cancel) return;
			resultsTitle = imp.getTitle()+" Deadpanorama";
			rt.show(resultsTitle);
			saveOpList();
		}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	}
	
	public String getResultsTitle(){
		return resultsTitle;
	}
	
	public void cancel(){
		cancel = true;
		if(opList!=null){
			for(Operation op : opList){
				op.cancel();
			}
		}
	}
	
	public ResultsTable getResultsTable(){
		return rt;
	}
	
	private void saveOpList(){
		try{
			final String nl = System.getProperty("line.separator");
			StringBuilder sb = new StringBuilder();
			for(Operation op:opList){
				sb.append( op.getCode()+nl );
			}
			String str = sb.toString();
			Files.write(path, str.getBytes());
		}catch(IOException e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	}
	
	private void loadOpList(){
		try{
			opList = new ArrayList<Operation>();
			if(!Files.exists(path)){
				//Files.createFile(path);
				return;
			}
			List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
			for(String line:lines){
				add( Operation.fromCode(line) );
			}
		}catch(IOException e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	}
	
}
