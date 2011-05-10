package tool;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import chemestry.AminoAcid;

public class NeighborDependentRamachandranDistribution {
	
	Double[][][][] rightNeightbour = new Double[AminoAcid.count][AminoAcid.count][360/5][360/5];
	Double[][][][] leftNeightbour = new Double[AminoAcid.count][AminoAcid.count][360/5][360/5];
	
	public NeighborDependentRamachandranDistribution() throws Exception {
		//try{
			FileInputStream fstream = new FileInputStream("/home/hkb/data/bachelor/Neighbor-dependent Ramachandran Distributions/NDRD_TCBIG.txt");
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			
			String strLine;
			while ((strLine = br.readLine()) != null)   {
				//System.out.println(strLine);
				if(strLine.length() > 0 && strLine.charAt(0) != '#') {
					String[] tokens = strLine.split("\\s+");
					
					if(tokens[2].equals("ALL") || tokens[0].equals("CPR") || tokens[2].equals("CPR"))
						continue;
					
					AminoAcid.Type aminoAcid = AminoAcid.Type.valueOf(tokens[0]);
					AminoAcid.Type neightbourAcid = AminoAcid.Type.valueOf(tokens[2]);
					int phi = Integer.valueOf(tokens[3]);
					int psi = Integer.valueOf(tokens[4]);
					double probability = Double.valueOf(tokens[5]);
					
					Double[][][][] neightbour = (tokens[1].equals("left")) ? this.leftNeightbour : this.rightNeightbour;

					neightbour[AminoAcid.typeToInt(aminoAcid)][AminoAcid.typeToInt(neightbourAcid)][(phi + 180)/5][(psi + 180)/5] = probability;
				}
			}
			
			in.close();
		/*} catch (Exception e){
			System.err.println("Error: " + e.getMessage());
		}*/
	}
	
	public double probability(AminoAcid.Type aminoAcid, AminoAcid.Type leftNeighbour, AminoAcid.Type rightNeighbour, double phi, double psi) {
		double left = this.leftNeightbour[AminoAcid.typeToInt(aminoAcid)][AminoAcid.typeToInt(leftNeighbour)][this.toBin(this.radianToDegree(phi))][this.toBin(this.radianToDegree(psi))];
		double right = this.leftNeightbour[AminoAcid.typeToInt(aminoAcid)][AminoAcid.typeToInt(rightNeighbour)][this.toBin(this.radianToDegree(phi))][this.toBin(this.radianToDegree(psi))];
		
		return (left + right) / 2;
	}
	
	private double radianToDegree(double radian) {
		return radian * (180 / Math.PI);
	}
	
	private int toBin(double angle) {
		return (int) (angle + 180)/5;
	}
}
