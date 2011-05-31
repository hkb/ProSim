package tool;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.swing.*;

import dataStructure.CTNode;
import dataStructure.ChainTree;

import edu.j3dScene.J3DScene;

public class BinaryTreePainter extends JFrame{
	private static final long serialVersionUID = 1L;
	public PNode root;
	public J3DScene scene;
	private PaintPane pp;
	
	
	public BinaryTreePainter(ChainTree cTree) {
		this(cTree.getRoot());
	}
	
	public BinaryTreePainter(CTNode node) {		
		super("ChainTree painter");
		super.setSize(600,500);
		
		setup();
		setRoot(new Node(node));
		super.setVisible(true);
		super.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	/**
	 * Wrapper class for CTNodes
	 * 
	 * @author hkb
	 */
	private class Node implements BinNode {
		
		CTNode node;
		
		public Node(CTNode node) {
			this.node = node;
		}

		@Override
		public void click() {
			// dummy
		}

		@Override
		public BinNode getLeft() {
			if(this.node.getLeft() != null) {
				return new Node(this.node.getLeft());
			}
			return null;
		}

		@Override
		public Color getLegColor() {
			return this.node.active ? Color.GREEN : this.node.isLocked ? Color.RED : Color.BLACK;
		}

		@Override
		public Color getNodeColor() {
			return this.node.active ? Color.GREEN : this.node.isLocked ? Color.RED : Color.BLACK;
		}

		@Override
		public BinNode getRight() {
			if(this.node.getRight() != null) {
				return new Node(this.node.getRight());
			}
			return null;
		}
		
		@Override
		public String toString() {
			return this.node.toString();
		}
	}
	
	
	/* --------- lol stuff ---------- */
	
	public BinaryTreePainter(BinNode treeRoot){
		super("ChainTree painter");
		super.setSize(600,500);
		
		setup();
		setRoot(treeRoot);
		super.setVisible(true);
	}
	
	void setup(){
		pp = new PaintPane();
		super.getContentPane().add(pp);
	}
	
	public void setRoot(BinNode treeRoot){
		this.root = new PNode(treeRoot,0);
		pp.updateComponents();
	}
	
	
	private class PaintPane extends JLayeredPane{
		private static final long serialVersionUID = 1L;
		Map<PNode,JComponent> componentMap = new HashMap<PNode,JComponent>();
		Map<JComponent,PNode> nodeMap = new HashMap<JComponent,PNode>();
		PaintPane(){	
			super();
		}
		double[] treeBorder = new double[4];
		int gap = 20;
		
		private void updateComponents(){
			List<PNode> nodes = new LinkedList<PNode>();
			root.collect(nodes);
			for(PNode n: nodes) {
				JComponent comp = new JComponent(){
					private static final long serialVersionUID = 1L;};
				super.add(comp,2000);
				Point p = toPoint(n.x,n.y);
				comp.setBounds(p.x-5,p.y-5,10,10);
				componentMap.put(n, comp);
				nodeMap.put(comp,n);
				comp.addMouseListener(new MouseAdapter(){
					public void mousePressed(MouseEvent evt){
						nodeMap.get(evt.getSource()).n.click();
						repaint();
					}
				});
			}
		}

		public void paint(Graphics g){
			Graphics2D g2d = (Graphics2D)g;
			g2d.setColor(Color.white);
			g2d.fillRect(0, 0, getWidth(), getHeight());
			g2d.setColor(Color.BLACK);
			super.paint(g);
			setTreeBorder();
			paintNode(root,g2d);
		}
		
		private Point paintNode(PNode n, Graphics2D g){
			Point p = toPoint(n.x,n.y);
			if(n.left!=null) {
				Point leftP = paintNode(n.left,g);
				g.setColor(n.n.getLegColor());
				g.drawLine(p.x, p.y, leftP.x, leftP.y);
			}
			if(n.right!=null) {
				Point rightP = paintNode(n.right,g);
				g.setColor(n.n.getLegColor());
				g.drawLine(p.x, p.y, rightP.x, rightP.y);
			}
			g.setColor(n.n.getNodeColor());
			g.fillOval(p.x-3, p.y-3, 6,6);
			g.setColor(Color.black);
			g.drawString(n.n.toString(), p.x, p.y);
			componentMap.get(n).setBounds(p.x-5, p.y-5, 10,10);
			
			return p;
		}
		void setTreeBorder(){
			List<PNode> nodes = new LinkedList<PNode>();
			root.collect(nodes);
			for(PNode n: nodes) {
				if(n.x<treeBorder[0]) treeBorder[0]=n.x;
				if(n.x>treeBorder[1]) treeBorder[1]=n.x;
				if(n.y<treeBorder[2]) treeBorder[2]=n.y;
				if(n.y>treeBorder[3]) treeBorder[3]=n.y;
			}
		}
		
		Point toPoint(double x, double y){
			double canonicalX = (x-treeBorder[0])/(treeBorder[1]-treeBorder[0]);
			double canonicalY = (y-treeBorder[2])/(treeBorder[3]-treeBorder[2]);
			return new Point((int)(canonicalX*(getWidth()-gap-gap))+gap,(int)(canonicalY*(getHeight()-gap-gap))+gap);
		}
	}

	public static class PNode{
		BinNode n;
		PNode left, right;
		double x,y;
		int height, depth;
		public PNode(BinNode n, int depth){
			this.n = n;
			this.depth = depth;
			if(n.getLeft()!=null) left = new PNode(n.getLeft(), depth+1);
			if(n.getRight()!=null) right = new PNode(n.getRight(), depth+1);
			
			if(n.getLeft()==null){
				if(n.getRight()==null) height=0;
				else height = right.height+1;
			}else {
				if(n.getRight()==null) height = left.height+1;
				else height = Math.max(left.height, right.height)+1; 
			}
			y = depth;
			x = 0;
			
			if(left!=null && right!=null){
				right.translate(1-horizontalDistance(left, right));
				double centerDX = (left.x+right.x)/2;
				left.translate(-centerDX);
				right.translate(-centerDX);
			}else if(left!=null){
				left.translate(-left.x-1);
			}else if(right!=null){
				right.translate(-right.x+1);
			}
		}

		private void translate(double dX){
			x+=dX;
			if(left!=null) left.translate(dX);
			if(right!=null) right.translate(dX);
		}
		
		public static double horizontalDistance(PNode left, PNode right){
			List<PNode> leftNodes = new LinkedList<PNode>();
			List<PNode> rightNodes = new LinkedList<PNode>();
			left.collect(leftNodes);
			right.collect(rightNodes);
			double minDist = Double.POSITIVE_INFINITY;
			for(int d=Math.min(left.depth,right.depth);;d++){
				double leftExtremeRight = Double.NEGATIVE_INFINITY;
				double rightExtremeLeft = Double.POSITIVE_INFINITY;
				boolean foundLeft= false, foundRight = false;
				for(PNode n: leftNodes){
					if(n.depth==d){
						foundLeft=true;
						if(n.x>leftExtremeRight) leftExtremeRight = n.x;
					}
				}
				for(PNode n: rightNodes){
					if(n.depth==d){
						foundRight = true;
						if(n.x<rightExtremeLeft) rightExtremeLeft = n.x;
					}
				}
				if(!foundLeft||!foundRight) break;
				minDist = Math.min(minDist, rightExtremeLeft-leftExtremeRight);
			}
			
			return minDist;
//			if(left.depth!=right.depth) throw new Error("Ahem!");
//			double dist = right.x-left.x;
//			PNode leftChild = left.right;
//			if(leftChild==null) leftChild = left.left;
//			PNode rightChild = right.left;
//			if(rightChild==null) rightChild = right.right;
//			if(rightChild==null||leftChild==null) return dist;
//			
//			return Math.min(dist, horizontalDistance(leftChild,rightChild));
		}
		private void collect(List<PNode> collection){
			collection.add(this);
			if(left!=null) left.collect(collection);
			if(right!=null) right.collect(collection);
		}
		public String toString(){
			return "PNode["+x+","+y+"]";
		}
	}
	
	public static interface BinNode{
		public BinNode getLeft();
		public BinNode getRight();
		public Color getNodeColor();
		public Color getLegColor();
		public void click();
		public String toString();
	}
}
