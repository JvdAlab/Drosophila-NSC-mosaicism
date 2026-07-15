package uk.ac.cam.gurdon.deadpanorama.operation;

import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

import ij.ImagePlus;
import ij.measure.ResultsTable;
import uk.ac.cam.gurdon.deadpanorama.DeadParams;
import uk.ac.cam.gurdon.deadpanorama.Nucleus3D;


public abstract class Operation{

	static Type TYPE;
	static int NCHANNELS;
	static String TIP;
	
	protected static final String TIPSTART = "<html><p>";
	private static final String TIPEND = "<br>Select \"outer\" to measure the surrounding volume instead of inside nuclei.</p></html>";
	private static final String[] CS = new String[]{ "-", "1", "2", "3", "4", "5", "6", "7", "8" };
	
	public static enum Type{
		STATISTICS("Statistics"), CORRELATION("Correlation"), SCATTER("Scatter Plot"),
		COUNTPOSITIVE("Count Positive"), SIGNALRATIO("Signal Ratio"), RENDER3D("3D Render");
		
		String str;
		Type(String str){
			this.str = str;
		}
		
		public String getName(){
			return str;
		}
		
		public static String[] getNames(){
			return Arrays.stream(values()).map(Type::getName).toArray(String[]::new);
		}
		
		public static Type getType(String str){
			for(Type type : values()){
				if(type.str.equals(str)) return type;
			}
			return null;
		}

	}
	
	/** static Operation factory
	 * @return	a new Operation of the specified Type
	 */
	public static Operation getOperation(Type type){
		if(type == Operation.Type.STATISTICS){
			return new OpStatistics();
		}
		else if(type == Operation.Type.CORRELATION){
			return new OpCorrelation();
		}
		else if(type == Operation.Type.SCATTER){
			return new OpScatter();
		}
		else if(type == Operation.Type.COUNTPOSITIVE){
			return new OpCountPositive();
		}
		else if(type == Operation.Type.SIGNALRATIO){
			return new OpSignalRatio();
		}
		else if(type == Operation.Type.RENDER3D){
			return new OpRender3D();
		}
		return null;
	}
	
	public static Operation fromCode(String str){
		try{
			String[] code = str.split(",");
	
			Type type = Type.getType(code[0]);
			if(type!=null){
				Operation op = getOperation(type);
				op.setChannelA( Integer.valueOf(code[1]) );
				if(code.length>=3){
					op.setChannelB( Integer.valueOf(code[2]) );
				}
				return op;
			}
		}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
		return null;
	}
	
	
	class CComboBox extends JComboBox<String>{
		private static final long serialVersionUID = -8225028998672676010L;

		public CComboBox(){
			super(CS);
			addItemListener(	//don't allow selection of index 0, keep consistent 1-based indexing to match ImageJ
					new ItemListener(){
						@Override
						public void itemStateChanged(ItemEvent ie) {
							if(getSelectedIndex()==0){
								setSelectedIndex(1);
							}
						}
					}
			);
		}
		
	}
	
	
	boolean cancel = false;
	JPanel panel;
	JComboBox<String> channelACombo, channelBCombo;
	JCheckBox outerTickA, outerTickB;
	
	public Operation(){
		
		panel = new JPanel();
		panel.add( new JLabel(TYPE.getName()) );
		
		if(NCHANNELS>=1){
			panel.add( new JLabel("C") );
			channelACombo = new CComboBox( );
			channelACombo.setSelectedIndex(1);
			panel.add( channelACombo );
			
			outerTickA = new JCheckBox("outer", false);
			panel.add(outerTickA);
		}
		if(NCHANNELS>=2){
			panel.add( new JLabel("vs C") );
			channelBCombo = new CComboBox( );
			channelBCombo.setSelectedIndex(2);
			panel.add( channelBCombo );
			
			outerTickB = new JCheckBox("outer", false);
			panel.add(outerTickB);
		}
		
		panel.setToolTipText(TIPSTART+TIP+TIPEND);
	}
	
	
	/** Run the operation and add results to the passed ResultsTable
	*/
	public abstract void run(ImagePlus imp, ArrayList<Nucleus3D> nuclei, DeadParams params, ResultsTable rt);
	
	/** Stop the operation before running on the next Nucleus3D in the list
	 */
	public void cancel() {
		cancel = true;
	}
	
	/**	Get the panel containing GUI components for configuring this Operation
	 */
	public JPanel getPanel() {
		return panel;
	}

	public String toString(){
		if(NCHANNELS==2){
			return TYPE.getName()+" C"+channelACombo.getSelectedIndex()+" vs C"+channelBCombo.getSelectedIndex();
		}
		else if(NCHANNELS==1){
			return TYPE.getName()+" C"+channelACombo.getSelectedIndex();
		}
		else{
			return TYPE.getName();
		}
	}
	
	/** Get a String code that can be saved and used to load this Operation
	 */
	public String getCode() {
		if(NCHANNELS==2){
			return TYPE.getName()+","+channelACombo.getSelectedIndex()+","+channelBCombo.getSelectedIndex();
		}
		else if(NCHANNELS==1){
			return TYPE.getName()+","+channelACombo.getSelectedIndex();
		}
		else{
			return TYPE.getName();
		}
	}
	
	public void setChannelA(int chan) {
		channelACombo.setSelectedIndex(chan);
	}

	public void setChannelB(int chan) throws Exception{
		if(NCHANNELS<2) throw new Exception("This operation does not have a second channel: "+toString());
		channelBCombo.setSelectedIndex(chan);
	}
	
}
