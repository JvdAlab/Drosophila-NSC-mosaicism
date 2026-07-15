package vos.de.kurt;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.RenderingHints;
import java.awt.color.ColorSpace;
import java.awt.color.ICC_ColorSpace;
import java.awt.color.ICC_Profile;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.prefs.Preferences;

import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JPanel;



/** Replacement for the real class to work around java.lang.IllegalArgumentException: Invalid service: net.imagej.legacy.LegacyConsoleService
 *  - RB
 */
public class CellCounterOptions {
	
	private static int NTYPES = 8;
	
	JDialog gui;
	TypeColourButton[] typeColourButtons;
	CellCntrImageCanvas ic;
	Preferences prefs;
	Color[] colours;
	
	public CellCounterOptions(CellCntrImageCanvas ic){
		this.ic = ic;
		prefs = Preferences.userRoot().node("Deadpanorama.CellCounterOptions");
		colours = new Color[NTYPES+1];	//1-based indexing
		for(int i=1;i<=NTYPES;i++){
			colours[i] = getColour(i);
		}
	}
	
	private class TypeColourButton extends JButton{
		private static final long serialVersionUID = 3341540561861279310L;
		private Dimension DIM = new Dimension(40,40);
		private Font FONT = new Font(Font.SANS_SERIF, Font.BOLD, 16);
		
		private int type;
		private Color colour;
		
		public TypeColourButton(int type){
			this.type = type;
			this.colour = colours[type];
			setMargin(new Insets(0,0,0,0));
			
			addActionListener( new ActionListener(){
				public void actionPerformed(ActionEvent ae){
					chooseColour();
				}
			} );
			
		}
		
		public void chooseColour() {
			Color c = JColorChooser.showDialog(gui, "Type "+type, colours[type]);
			if(c!=null)	setColour(c);
		}
		
		public void setColour(Color c){
			colour = c;
			String str = colour.getRed()+","+colour.getGreen()+","+colour.getBlue()+","+colour.getAlpha();
			prefs.put("C"+type, str);
			colours[type] = colour;
			repaint();
			//ic.repaint();
		}
		
		@Override
		public Dimension getPreferredSize(){
			return DIM;
		}
		
		@Override
		public void paintComponent(Graphics g1d){
			Graphics2D g = (Graphics2D) g1d;
			Dimension size = getSize();
			g.setColor(colour);
			g.fillRect(0,0,size.width+1,size.height+1);
			g.setColor(Color.BLACK);
			g.setFont(FONT);
			FontMetrics fm = g.getFontMetrics();
			String str = ""+type;
			g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
			g.drawString(str, (int)((size.width - fm.stringWidth(str)) / 2d), (int)((size.height - fm.getHeight()) / 2d) + fm.getAscent() );
		}
		
	}
	
	// default colour for given index out of NTYPES
	private Color getDefaultColour(int i){
		float nf = (float)(NTYPES);
		//float h = (i%2==0) ? i/nf : 1f-(i/nf);
		//float h = i/nf;
		//float h = (i%2==0) ? i/nf : ((i/nf)+0.33f)%1f;
			 if(i==1) return Color.RED;
		else if(i==2) return Color.GREEN;
		else if(i==3) return Color.BLUE;
		else if(i==4) return Color.YELLOW;
		else if(i==5) return Color.MAGENTA;
		else if(i==6) return Color.CYAN;
		else if(i==7) return Color.ORANGE;
		else if(i==8) return Color.PINK;
		else{
			float f = (i-1)/nf;
			float h = (i%2==0) ? f : 1f-f;
			float s = 1.0f;
			float v = 1.0f;
			Color hsv = Color.getHSBColor(h,s,v);
			return hsv;
		}
	}
	
	public Color getColour(int id) {
		if(colours[id]!=null) return colours[id];
		String key = "C"+id;
		
		Color hsv = getDefaultColour(id);
		int rdef = hsv.getRed();
		int gdef = hsv.getGreen();
		int bdef = hsv.getBlue();
		int adef = 255;
		
		
		String def = rdef+","+gdef+","+bdef+","+adef;
		String str = prefs.get(key, def);
		//System.out.println(id+" : "+def+" ~ "+str);
		if(!str.matches("\\d+,\\d+,\\d+,\\d+")){	//format is R,G,B,A int
			str = def;
		}
		String[] rgb = str.split(",");
		int r = Integer.parseInt(rgb[0]);
		int g = Integer.parseInt(rgb[1]);
		int b = Integer.parseInt(rgb[2]);
		int a = Integer.parseInt(rgb[3]);
		return new Color(r,g,b,a);
	}
	
	public static Color colour(int id) {
		Preferences prefs = Preferences.userRoot().node("Deadpanorama.CellCounterOptions");
		String key = "C"+id;
		String def = "0,255,0,255";
		String str = prefs.get(key, def);
		if(!str.matches("\\d+,\\d+,\\d+,\\d+")){	//format is R,G,B,A int
			str = def;
		}
		String[] rgb = str.split(",");
		int r = Integer.parseInt(rgb[0]);
		int g = Integer.parseInt(rgb[1]);
		int b = Integer.parseInt(rgb[2]);
		int a = Integer.parseInt(rgb[3]);
		return new Color(r,g,b,a);
	}
	
	public void show(){
		if(gui==null){
			gui = new JDialog(new Frame(), "Cell Counter Options");
			gui.setLayout(new BorderLayout());
			gui.setModal(true);
			JPanel main = new JPanel();
			int cols = 4;
			main.setLayout(new GridLayout(0,cols,2,2));
			typeColourButtons = new TypeColourButton[33];
			for(int i=1;i<=NTYPES;i++){			//TODO: max number of types???
				TypeColourButton chooser = new TypeColourButton(i);
				main.add( chooser );
				typeColourButtons[i] = chooser;
			}
			gui.add(main, BorderLayout.CENTER);
			
			ActionListener buttonListener = new ActionListener(){
				public void actionPerformed(ActionEvent ae){
					String event = ae.getActionCommand();
					if(event.equals("OK")){
						gui.setVisible(false);
						if(ic!=null){ 
							ic.repaint();
							ic.getCellCounter().populateTxtFields();	//update colours
						}
					}
					else if(event.equals("Reset")){
						for(int i=1;i<=NTYPES;i++){
							Color hsv = getDefaultColour(i);
							colours[i] = hsv;
							typeColourButtons[i].setColour(hsv);
						}
					}
				}
			};
			
			JPanel buttonPan = new JPanel();
			JButton ok = new JButton("OK");
			ok.addActionListener(buttonListener);
			buttonPan.add(ok);
			JButton reset = new JButton("Reset");
			reset.addActionListener(buttonListener);
			buttonPan.add(reset);
			gui.add(buttonPan, BorderLayout.SOUTH);
			gui.pack();
			gui.setResizable(false);
		}
		gui.setLocationRelativeTo(null);
		gui.setVisible(true);
	}
	
}
