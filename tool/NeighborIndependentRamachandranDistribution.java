package tool;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.List;

import math.Tuple2;

import chemestry.AminoAcid;
import dataStructure.ChainTree;

public class NeighborIndependentRamachandranDistribution {

	private static int BIN_SIZE = 5;
	
	Double[][][] probabilities = new Double[AminoAcid.count][360/BIN_SIZE][360/BIN_SIZE];

	public NeighborIndependentRamachandranDistribution() throws Exception {

		//try{
			FileInputStream fstream = new FileInputStream("/home/hkb/data/bachelor/Neighbor-dependent Ramachandran Distributions/NDRD_TCBIG.txt");
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			
			String strLine;
			while ((strLine = br.readLine()) != null)   {
				//System.out.println(strLine);
				if(strLine.length() > 0 && strLine.charAt(0) != '#') {
					String[] tokens = strLine.split("\\s+");
					
					if(!tokens[2].equals("ALL") || tokens[0].equals("CPR"))
						continue;
					
					AminoAcid.Type aminoAcid = AminoAcid.Type.valueOf(tokens[0]);
					int phi = Integer.valueOf(tokens[3]);
					int psi = Integer.valueOf(tokens[4]);
					double probability = Double.valueOf(tokens[5]);

					this.probabilities[AminoAcid.typeToInt(aminoAcid)][(phi + 180)/5][(psi + 180)/5] = probability;
				}
			}
			
			in.close();
		/*} catch (Exception e){
			System.err.println("Error: " + e.getMessage());
		}*/
	}
	

	public double probability(AminoAcid.Type aminoAcid, double phi, double psi) {		
		return this.probabilities[AminoAcid.typeToInt(aminoAcid)][this.toBin(this.radianToDegree(phi))][this.toBin(this.radianToDegree(psi))];
	}
	
	public static double radianToDegree(double radian) {
		return radian * (180 / Math.PI);
	}
	
	private int toBin(double angle) {
		return (int) (angle + 180)/BIN_SIZE;
	}
}
