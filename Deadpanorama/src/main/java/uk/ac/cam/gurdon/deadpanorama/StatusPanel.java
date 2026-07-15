package uk.ac.cam.gurdon.deadpanorama;

import java.awt.Dimension;

import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import javax.swing.text.Document;

import ij.IJ;

public class StatusPanel extends JScrollPane{
	private static final long serialVersionUID = -3701927983210537069L;
	private static final Dimension dim = new Dimension(100,20);
	
	private JTextArea ta;
	

	public StatusPanel(){
		super();
		ta = new JTextArea();
		ta.setEditable(false);
		setViewportView( ta );
	}
	
	/**	Always log to the JTextArea. Optional Boolean args set whether to also log to System.out and IJ.log.
	 */
	public void log(String txt, Boolean... out){
		if(out.length>1&&out[0]){
			System.out.println(txt);
		}
		if(out.length>2&&out[1]){
			IJ.log(txt);
		}
		SwingUtilities.invokeLater(
				new Runnable(){
					public void run(){
						ta.append(txt+"\n");
						Document doc = ta.getDocument();
						int length = doc.getLength();
						ta.setCaretPosition(length);
						repaint();
					}
				}
		);
	}
	
	@Override
	public Dimension getPreferredSize(){
		return dim;
	}
	
}