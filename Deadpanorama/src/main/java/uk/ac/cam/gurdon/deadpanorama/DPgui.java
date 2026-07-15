package uk.ac.cam.gurdon.deadpanorama;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Timer;
import java.util.TimerTask;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.ToolTipManager;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.process.AutoThresholder;
import ij.process.StackStatistics;
import ij.text.TextPanel;
import ij.text.TextWindow;
import uk.ac.cam.gurdon.deadpanorama.operation.Operation;
import uk.ac.cam.gurdon.deadpanorama.operation.OperationPanel;
import vos.de.kurt.CellCntrMarker;
import vos.de.kurt.CellCntrMarkerVector;
import vos.de.kurt.CellCounter;


public class DPgui extends JFrame implements ActionListener{
	private static final long serialVersionUID = 2254811415033008281L;
	
	private static final String tag0 = "<html><p>";
	private static final String tag1 = "</p></html>";
	
	private static enum Tip{
		NONE(""),
		ADVANCED("Enable additional controls for advanced parameters. Default values optimised for standard use cases<br>will be applied if advanced controls are not shown. Press the log button to output all current parameters<br>to the ImageJ log."),
		OVERLAY("Toggle image overlay visibility in the original image or the Cell Counter copy if it is open."),
		OPERATIONS("Add or remove analysis operations to be run on the segmented nuclei."),
		
		DEADPANC("Channel number to use for nucleus mapping."),
		SLICES("Range of Z-slices to use for nucleus mapping."),
		SIGMA("Gaussian standard deviation - use to set the smallest nucleus radius."),
		SCALEK("Scaling factor for sigma in difference of Gaussians - use to set the largest nucleus radius.  Default = "+DeadParams.DEFAULT_K),
		K("Initial number of intensity levels for Hierarchical K-Means clustering - increase to detect nuclei with a<br>larger range of intensities.  Default = "+DeadParams.DEFAULT_STARTK),
		WATERSHED("Apply watershed segmentation to mapped objects before filtering. Use to detect more objects in<br>areas of continuous intensity, but often over-segments."),
		
		AREA("Only 2D objects within this area range will be used for nucleus reconstruction."),
		THRESHOLD("Minimum intensity threshold as a scale normalised value 0..1. Press the gear button to set the Otsu<br>auto-threshold for the current deadpan channel or selected area of it."),
		CIRCULARITY("Minimum 2D object circularity, 1 is a circle."),
		FLATNESS("Minimum density flatness \u221A(\u03C3\u00B2/\u00B5A)\u207B\u00B9 - use to exclude trachea fragments. Default = "+DeadParams.DEFAULT_MINFLATNESS),
		//flatness = 1/sqrt(variance/integrated density)
		
		SLICEDIST("Maximum Z-slice distance between fragments of nuclei to be joined. The minimum is 1 to join<br>fragments only in adjacent slices.  Default = "+DeadParams.DEFAULT_JOINZ),
		OVERLAP("Minimum projected overlap proportion for fragments of nuclei to be joined, 1 is complete overlap<br>between slices.  Default = "+DeadParams.DEFAULT_JOINOVERLAP);
		
		String txt;
		private Tip(String txt){
			this.txt = txt;
		}
		
		public String getText(){
			return tag0+txt+tag1;
		}
		
	}
	
	private CardLayout card;
	
	private Deadpanorama dp;
	private KSpinner deadpanCSpinner, startZSpinner, endZSpinner, startKSpinner, joinZSpinner;	//int
	private KSpinner sigmaSpinner, dogKSpinner, minASpinner, maxASpinner, minCircSpinner, joinOverlapSpinner, thresholdSpinner, flatnessSpinner;	//double
	private JCheckBox advancedTick, watershedTick, showOverlayTick, ccOverlayTick;
	private StatusPanel status;
	private JPopupMenu typeMenu;
	
	private ArrayList<KPanel> panels;
	
	private OperationPanel opPanel;
	
	private CellCounter cc;
	private OvalRoi marker;
	
	
	private class KSpinner extends JSpinner{
		private static final long serialVersionUID = -1874314266594266222L;
		private final Dimension DIM = new Dimension(50,20);
		
		public KSpinner(int value, int min, int max, int step){
			super( new SpinnerNumberModel(value, min, max, step) );
		}
		public KSpinner(double value, double min, double max, double step){
			super( new SpinnerNumberModel(value, min, max, step) );
		}
		
		@Override
		public Dimension getPreferredSize(){
			return DIM;
		}
		
	}
	
	private class KPanel extends JPanel{
		private static final long serialVersionUID = -4784577659478968533L;
		public final Color ADVANCED_COLOUR = new Color(196, 172, 172);
		public boolean advanced;
		
		public KPanel(boolean advanced, Tip tip, Object... comps){
			this.advanced = advanced;
			if(advanced) setBackground(ADVANCED_COLOUR);
			setLayout(new FlowLayout(FlowLayout.CENTER, 2, 2));
			if(tip!=Tip.NONE)	setToolTipText(tip.getText());

			add(Box.createVerticalStrut(20));	//min height even without components
			for(Object obj : comps){
				if(obj instanceof JComponent){
					if(tip!=Tip.NONE)	((JComponent)obj).setToolTipText(tip.getText());
					if(advanced) ((JComponent)obj).setBackground(ADVANCED_COLOUR);
					add((JComponent)obj);
				}
				else if(obj instanceof String){
					JLabel label = new JLabel((String)obj);
					if(tip!=Tip.NONE)	label.setToolTipText(tip.getText());
					add(label);
				}
				else if(obj instanceof Integer){
					add(Box.createHorizontalStrut((int) obj));
				}
				else if(obj instanceof Color){
					setBackground((Color) obj);
				}
			}

			panels.add(this);
		}
		
	}
	
	/** A circular button with a choice of simple icons
	 */
	private class KButton extends JButton{
		private static final long serialVersionUID = -747890320829992807L;
		private static final int PLUS=-1, MINUS=-2, AUTO=-3, LOG=-4;
		private static final int SIZE = 21;
		private int p0 = 4, p1 = 10, p2 = 16;
		
		private int icon;
		private Insets ninsets = new Insets(0, 0, 0, 0);
		private Dimension dim = new Dimension(SIZE, SIZE);
		private Shape shape = new Ellipse2D.Float(0, 0, SIZE, SIZE);
		private BasicStroke stroke1 = new BasicStroke(1f);
		private BasicStroke stroke2 = new BasicStroke(2f);
		private BasicStroke stroke4 = new BasicStroke(4f);
		
		public KButton(int icon){
			this.icon = icon;
			setMargin(ninsets);
			setContentAreaFilled(false);
			setFocusPainted(false);
		}

		@Override
		public void paintComponent(Graphics g1d){
			Graphics2D g = (Graphics2D) g1d;
			Color fillColour = Color.LIGHT_GRAY;
			Color drawColour = Color.BLACK;
			if (getModel().isArmed()) {
				fillColour = Color.GRAY;
				drawColour = Color.LIGHT_GRAY;
			}
			g.setColor(fillColour);
			g.fillOval(0, 0, SIZE, SIZE );

			g.setColor(drawColour);
			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			
			if(icon==LOG){
				g.setStroke( stroke2 );
				for(int y=p0;y<=p2;y+=4){
					g.drawLine(p0+2,y,p2-2,y);
				}
				return;
			}
			else{
				g.setStroke( stroke4 );
				g.drawLine(p0,p1,p2,p1);	//horizontal for PLUS, MINUS and AUTO
				if(icon==PLUS||icon==AUTO)	g.drawLine(p1,p0,p1,p2);	//vertical
				if(icon==AUTO){
					g.drawLine(p0+2,p0+2,p2-2,p2-2);	//diagonals
					g.drawLine(p0+2,p2-2,p2-2,p0+2);
					g.setColor(fillColour);
					g.fillOval(p1-2,p1-2,5,5);			//hole
				}
			}
		}
		
		@Override
		public void paintBorder(Graphics g1d) {
			Graphics2D g = (Graphics2D) g1d;
		    g.setColor(Color.DARK_GRAY);
		    g.setStroke(stroke1);
		    g.drawOval(0, 0, SIZE-1, SIZE-1 );
		}

		@Override
		public boolean contains(int x, int y) {
			return shape.contains(x, y);
		}
		
		@Override
		public Dimension getPreferredSize(){
			return dim;
		}
		
	}
	
	public DPgui(Deadpanorama dp){
		super("Deadpanorama");
		this.dp = dp;
		this.panels = new ArrayList<KPanel>();

		setIconImage( Toolkit.getDefaultToolkit().getImage(getClass().getResource("logo_icon.gif")) );
		ToolTipManager.sharedInstance().setDismissDelay(10000);
		
		DeadParams params = DeadParams.getPrefs();
		
		minASpinner = new KSpinner(params.minA, 0.0, 10000.0, 1.0);
		maxASpinner = new KSpinner(params.maxA, 0.0, 10000.0, 1.0);
		
		deadpanCSpinner = new KSpinner(params.deadpanC, 1,16, 1);
		startZSpinner = new KSpinner(params.startZ, 1,512, 1);
		
		startKSpinner = new KSpinner(params.startK, 1,512, 1);
		endZSpinner = new KSpinner(params.endZ, 1,512, 1);
		
		minCircSpinner = new KSpinner(params.minCirc, 0.0,1.0, 0.01);
		joinZSpinner = new KSpinner(params.joinZ, 1,20,1);
		joinOverlapSpinner = new KSpinner(params.joinOverlap, 0.0,1.0, 0.01);
		flatnessSpinner = new KSpinner(params.minFlatness, 0.0,1.0, 0.01);
		thresholdSpinner = new KSpinner(params.threshold, 0.0,1.0, 0.01);
		
		sigmaSpinner = new KSpinner(params.sigma, 0.0, 100.0, 0.1);
		dogKSpinner = new KSpinner(params.K, 0.0, 100.0, 0.1);
		
		
		card = new CardLayout();
		setLayout(card);
		JPanel main = new JPanel();
		main.setLayout( new BoxLayout(main, BoxLayout.Y_AXIS) );
		
		JPanel deadPanel = new JPanel();
		deadPanel.setLayout( new BoxLayout(deadPanel, BoxLayout.Y_AXIS) );
		
		advancedTick = new JCheckBox("Advanced", params.advanced);
		advancedTick.addActionListener(this);
		KButton paramsButton = new KButton(KButton.LOG);
		paramsButton.setActionCommand("Log Parameters");
		paramsButton.addActionListener(this);
		deadPanel.add( new KPanel(false, Tip.ADVANCED, advancedTick, paramsButton) );
		
		deadPanel.add( new KPanel(false, Tip.NONE, "Deadpan Mapping", Color.GRAY) );
		deadPanel.add( new KPanel(false, Tip.DEADPANC, "Deadpan Channel: ", deadpanCSpinner) );
		deadPanel.add( new KPanel(false, Tip.SLICES, "Slice Range: ", startZSpinner, "-", endZSpinner) );
		deadPanel.add( new KPanel(false, Tip.SIGMA, "Sigma: ", sigmaSpinner, "\u00b5m") );
		deadPanel.add( new KPanel(true, Tip.SCALEK, "Sigma Factor: ", dogKSpinner) );
		deadPanel.add( new KPanel(true, Tip.K, "Starting K: ", startKSpinner) );
		watershedTick = new JCheckBox("Watershed", params.watershed);
		deadPanel.add( new KPanel(true, Tip.WATERSHED, watershedTick) );
		
		deadPanel.add( new KPanel(false, Tip.NONE) );
		deadPanel.add( new KPanel(false, Tip.NONE, "Filtering", Color.GRAY) );
		deadPanel.add( new KPanel(false, Tip.AREA, "Area Range: ", minASpinner, "-", maxASpinner, "\u00b5m\u00b2") );
		KButton threshButton = new KButton(KButton.AUTO);
		threshButton.setActionCommand("Auto-Threshold");
		threshButton.addActionListener(this);
		deadPanel.add( new KPanel(false, Tip.THRESHOLD, "Threshold: ", thresholdSpinner, threshButton) );
		deadPanel.add( new KPanel(false, Tip.CIRCULARITY, "Min Circularity: ", minCircSpinner) );
		deadPanel.add( new KPanel(true, Tip.FLATNESS, "Min Flatness: ", flatnessSpinner) );

		deadPanel.add( new KPanel(true, Tip.NONE, getBackground()) );
		deadPanel.add( new KPanel(true, Tip.NONE, "3D Reconstruction", Color.GRAY) );
		deadPanel.add( new KPanel(true, Tip.SLICEDIST, "Max Slice Distance: ", joinZSpinner) );
		deadPanel.add( new KPanel(true, Tip.OVERLAP, "Min Join Overlap: ", joinOverlapSpinner) );
		
		deadPanel.add( new KPanel(false, Tip.NONE) );
		deadPanel.add( new KPanel(false, Tip.NONE, "Results", Color.GRAY) );
		JButton toCC = new JButton("To Cell Counter");
		toCC.addActionListener(this);
		
		JButton target = new JButton("Target from Table");
		target.addActionListener(this);
		
		deadPanel.add(new KPanel(false, Tip.NONE, toCC, target));
		showOverlayTick = new JCheckBox("Original Image", true);
		showOverlayTick.setActionCommand("Overlay in Original Image");
		showOverlayTick.addActionListener(this);
		ccOverlayTick = new JCheckBox("Cell Counter", true);
		ccOverlayTick.setActionCommand("Overlay in Cell Counter");
		ccOverlayTick.addActionListener(this);
		deadPanel.add(new KPanel(false, Tip.OVERLAY, "Overlay ", showOverlayTick, ccOverlayTick));
		
		main.add(deadPanel);
		main.add( new KPanel(false, Tip.NONE) );
		
		main.add( new KPanel(false, Tip.NONE, "Analysis", Color.GRAY) );
		KButton addOpButton = new KButton(KButton.PLUS);
		addOpButton.setActionCommand("add operation");
		addOpButton.addActionListener(this);
		KButton remOpButton = new KButton(KButton.MINUS);
		remOpButton.setActionCommand("remove operation");
		remOpButton.addActionListener(this);
		opPanel = new OperationPanel();	
		main.add( new KPanel(false, Tip.OPERATIONS, addOpButton, remOpButton) );
		main.add(opPanel);
		main.add( new KPanel(false, Tip.NONE) );
		
		JButton run = new JButton("Run");
		run.addActionListener(this);
		JButton preview = new JButton("Preview");
		preview.addActionListener(this);
		JButton cancel = new JButton("Cancel");
		cancel.addActionListener(this);
		main.add(new KPanel(false, Tip.NONE,  run, preview, 20, cancel));
		
		card.addLayoutComponent(main, "main");
		add(main);
		
		JPanel working = new JPanel();
		working.setLayout(new BorderLayout());
		ImageIcon img = new ImageIcon( Toolkit.getDefaultToolkit().getImage(getClass().getResource("work.gif")) );
		working.add( new KPanel(false, Tip.NONE, new JLabel(img)), BorderLayout.NORTH );
		JButton stop = new JButton("Stop");
		stop.setBackground(Color.RED);
		stop.addActionListener(this);
		working.add( new KPanel(false, Tip.NONE, stop), BorderLayout.SOUTH );
		status = new StatusPanel();
		working.add(status, BorderLayout.CENTER );
		
		card.addLayoutComponent(working, "working");
		add(working);
		
		setAdvanced(false);
		
		card.show(getContentPane(), "main");
		
		pack();
		setLocationRelativeTo(null);
		setVisible(true);
	}
	
	private DeadParams getParams(){
		
		DeadParams params = new DeadParams();
		
		params.advanced = advancedTick.isSelected();
		params.deadpanC = (int) deadpanCSpinner.getValue();
		params.startZ = (int) startZSpinner.getValue();
		params.endZ = (int) endZSpinner.getValue();
		
		params.sigma = (double) sigmaSpinner.getValue();
		params.minA = (double) minASpinner.getValue();
		params.maxA = (double) maxASpinner.getValue();
		params.threshold = (double) thresholdSpinner.getValue();
		params.minCirc = (double) minCircSpinner.getValue();

		if(params.advanced){
			params.K = (double) dogKSpinner.getValue();
			params.startK = (int) startKSpinner.getValue();
			params.joinZ = (int) joinZSpinner.getValue();
			params.joinOverlap = (double) joinOverlapSpinner.getValue();
			params.minFlatness = (double) flatnessSpinner.getValue();
			params.watershed = watershedTick.isSelected();
		}
		else{	//use defaults when not in advanced mode
			params.K = DeadParams.DEFAULT_K;
			params.startK = DeadParams.DEFAULT_STARTK;
			params.joinZ = DeadParams.DEFAULT_JOINZ;
			params.joinOverlap = DeadParams.DEFAULT_JOINOVERLAP;
			params.minFlatness = DeadParams.DEFAULT_MINFLATNESS;
			params.watershed = DeadParams.DEFAULT_WATERSHED;
		}
		
		return params;
	}
	
	private void setAdvanced(boolean adv){
		for(KPanel panel : panels){
			panel.setVisible( !panel.advanced || adv );
		}
		pack();
		setLocationRelativeTo(null);
	}
	
	private void target(ImagePlus imp) {
		try{
			if (marker != null)	return; // do nothing if there is already a marker - marker is set to null when the display timer finishes
			String title = opPanel.getResultsTitle(); // ResultsTable title set by OperationPanel
			Frame frame = WindowManager.getFrame(title);
			if (frame == null || frame instanceof TextWindow == false) {
				JOptionPane.showMessageDialog(this, "No ResultsTable.", "Deadpanorama", JOptionPane.ERROR_MESSAGE);
				return;
			}
			TextPanel tp = ((TextWindow) frame).getTextPanel();
			int start = tp.getSelectionStart();
			if (start == -1) {
				JOptionPane.showMessageDialog(this, "No row selected.", "Deadpanorama", JOptionPane.ERROR_MESSAGE);
				return;
			}
			String[] head = tp.getColumnHeadings().split("\\t"); // ResultsTable headers
			String[] line = tp.getLine(start).split("\\t"); 	 // ResultsTable selected line
			Double x = null, y = null, z = null;
			for (int i = 0; i < head.length; i++) {
				if (head[i].equals("X")) {
					x = Double.valueOf(line[i]);
				} else if (head[i].equals("Y")) {
					y = Double.valueOf(line[i]);
				} else if (head[i].equals("Z")) {
					z = Double.valueOf(line[i]);
				}
			}
			if (x != null && y != null && z != null) {
				Calibration cal = imp.getCalibration();
				int xi = (int) (x / cal.pixelWidth);
				int yi = (int) (y / cal.pixelHeight);
				int zi = (int) (z / cal.pixelDepth);
	
				Overlay lol = imp.getOverlay();
				if (lol == null) {
					lol = new Overlay();
					imp.setOverlay(lol);
				}
				final Overlay ol = lol;
				marker = new OvalRoi(xi - 3, yi - 3, 6, 6);
				ol.add(marker);
				
				imp.getWindow().toFront();
				
				Timer timer = new Timer();
				TimerTask task = new TimerTask() {
					int n = 0;
					int maxn = 100;
					Color colour = Color.YELLOW;
	
					public void run() {
						try{
							ol.remove(marker);
							n++;
							if (n >= maxn) {
								marker = null;
								imp.updateAndDraw();
								timer.cancel();
							} else {
								int r = maxn - n;
								marker = new OvalRoi(xi - r, yi - r, 2 * r, 2 * r);
								marker.setStrokeWidth(n);
								marker.setStrokeColor(colour);
								imp.setPosition(imp.getChannel(), zi, imp.getFrame());
								ol.add(marker);
							}
							imp.updateAndDraw();
						}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
					}
				};
				timer.scheduleAtFixedRate(task, 0, 5);
	
			}
		}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	}
	
	@Override
	public void actionPerformed(ActionEvent ae) {
		SwingUtilities.invokeLater(
			new Runnable(){
				public void run(){

					String event = ae.getActionCommand();

					if(event.equals("Cancel")){	//main gui cancel button
						dispose();
						return;
					}
					else if(event.equals("Stop")){
						status.log("Stopping...");
						dp.kill();
						card.show(getContentPane(), "main");
						return;
					}
					else if(event.equals("Advanced")){
						setAdvanced(advancedTick.isSelected());
						return;
					}

					ImagePlus imp = WindowManager.getCurrentImage();
					if(imp==null) {
						JOptionPane.showMessageDialog(null, "No image open.", "No Image", JOptionPane.ERROR_MESSAGE);
						return;
					}
					
					DeadParams params = getParams();
					if(params.deadpanC>imp.getNChannels()){
						params.deadpanC = 1;
						deadpanCSpinner.setValue(1);
					}
					if(params.startZ>imp.getNSlices()){
						params.startZ = imp.getNSlices();
						startZSpinner.setValue(params.startZ);
					}
					if(params.endZ>imp.getNSlices()){
						params.endZ = imp.getNSlices();
						endZSpinner.setValue(params.endZ);
					}
					if(params.startZ>params.endZ){
						params.startZ = params.endZ;
						startZSpinner.setValue(params.startZ);
					}

					if(event.equals("To Cell Counter")){
						try{
							cc = new CellCounter();
							cc.setKeepOriginal(true);
							cc.initializeImage();
							CellCounter.setType("1");	//Stringly typed 1-based type index
							CellCntrMarkerVector markerVector = cc.getCurrentMarkerVector();	//Vector for current type
							ArrayList<Nucleus3D> nuclei = dp.getNuclei();
							if(nuclei!=null){
								for(Nucleus3D nuc : nuclei){
									int[] xyz = nuc.getIntCentroid();
									xyz[2] = imp.getStackIndex(params.deadpanC, xyz[2], 1);	//uses stack index for z, not slice index
									CellCntrMarker mark = new CellCntrMarker(xyz[0], xyz[1], xyz[2]);
									markerVector.add(mark);
								}
							}
							cc.populateTxtFields();	//display count
							ImagePlus ccimp = cc.getCounterImage();
							if(ccimp.getOverlay()!=null){
								ccimp.setHideOverlay(!ccOverlayTick.isSelected());
							}
						}catch(Exception e){System.out.print(e.toString()+"\n~~~~~\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
						return;
					}
					else if(event.equals("Auto-Threshold")){
						ImagePlus dpImp = new Duplicator().run(imp, params.deadpanC,params.deadpanC, params.startZ,params.endZ, 1,1);
						int[] hist = new StackStatistics( dpImp ).histogram;
						dpImp.close();

						int threshi = new AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, hist );
						double threshf = Math.round( (threshi/256d) * 100d ) / 100d;	//threshold as factor rounded to 2 decimal places
						thresholdSpinner.setValue( threshf );
						return;
					}
					else if(event.equals("add operation")){
						if(typeMenu==null){
							ActionListener addListener = new ActionListener(){
								public void actionPerformed(ActionEvent ev){
									String str = ev.getActionCommand();
									Operation.Type type = Operation.Type.getType(str);
									opPanel.add(type);
									pack();
								}
							};
							typeMenu = new JPopupMenu();
							for(String type : Operation.Type.getNames()){
								JMenuItem item = new JMenuItem(type);
								item.addActionListener(addListener);
								typeMenu.add(item);
							}
						}
						Component comp = (Component) ae.getSource();
						typeMenu.show(comp, comp.getWidth(), 0);
					}
					else if(event.equals("remove operation")){
						opPanel.remove();
						pack();
					}
					else if(event.equals("Overlay in Original Image")){
						if(imp.getOverlay()!=null){
							imp.setHideOverlay(!showOverlayTick.isSelected());
						}
					}
					else if(event.equals("Overlay in Cell Counter")){
						if(cc==null) return;
						ImagePlus ccimp = cc.getCounterImage();
						if(ccimp.getOverlay()!=null){
							ccimp.setHideOverlay(!ccOverlayTick.isSelected());
						}
					}
					else if(event.equals("Target from Table")){
						target(imp);							
					}
					else if(event.equals("Log Parameters")){
						IJ.log( params.toString() );
					}
					else if(event.equals("Run")||event.equals("Preview")){
						SwingWorker<Object, Void> worker = new SwingWorker<Object, Void>() {
							public Object doInBackground() {
								card.show(getContentPane(), "working");
								if(event.equals("Run")){
									dp.run(imp, params, opPanel, status);
									params.save();
								}
								else if(event.equals("Preview")){
									params.startZ = imp.getSlice();
									params.endZ = imp.getSlice();
									dp.run(imp, params, null, status);
								}
								if(imp.getOverlay()!=null){
									imp.setHideOverlay(!showOverlayTick.isSelected());
								}
								card.show(getContentPane(), "main");
								return null;
							}
						};
						worker.execute();
					}
				}
			}
		);
	}
	
}
