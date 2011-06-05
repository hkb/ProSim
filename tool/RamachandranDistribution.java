package tool;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.HashMap;
import java.util.Map;

import math.Tuple2;

import chemestry.AminoAcid;
import chemestry.AminoAcid.Type;

public class RamachandranDistribution {
	
	private static double RADIAN_TO_DEGREE_FACTOR = 180 / Math.PI;
	
	private Double[][][][] rightNeightbour = new Double[AminoAcid.count][AminoAcid.count][360/5][360/5];
	private Double[][][][] leftNeightbour = new Double[AminoAcid.count][AminoAcid.count][360/5][360/5];
	
	private Double[][][] singleProbabilities = new Double[AminoAcid.count][360/5][360/5];
		
	public RamachandranDistribution() {
		this("/home/hkb/data/bachelor/Neighbor-dependent Ramachandran Distributions/NDRD_TCB.txt");
	}
	
	public RamachandranDistribution(String dataFile) {
		try{
			BufferedReader data = new BufferedReader(new FileReader(dataFile));
			
			String strLine;
			while ((strLine = data.readLine()) != null)   {
				//System.out.println(strLine);
				if(strLine.length() > 0 && strLine.charAt(0) != '#') {
					String[] tokens = strLine.split("\\s+");
					
					if(tokens.length == 0 || tokens[0].equals("CPR") || tokens[2].equals("CPR"))
						continue;
					
					AminoAcid.Type aminoAcid = AminoAcid.Type.valueOf(tokens[0]);
					int phi = Integer.valueOf(tokens[3]);
					int psi = Integer.valueOf(tokens[4]);
					double probability = Double.valueOf(tokens[5]);
					
					if (tokens[2].equals("ALL")) {
						if(this.singleProbabilities[AminoAcid.typeToInt(aminoAcid)][toBin(phi)][toBin(psi)] == null) {
							this.singleProbabilities[AminoAcid.typeToInt(aminoAcid)][toBin(phi)][toBin(psi)] = probability;
						} else {
							this.singleProbabilities[AminoAcid.typeToInt(aminoAcid)][toBin(phi)][toBin(psi)] += probability;
							this.singleProbabilities[AminoAcid.typeToInt(aminoAcid)][toBin(phi)][toBin(psi)] /= 2;
						}
					} else {
						AminoAcid.Type neightbourAcid = AminoAcid.Type.valueOf(tokens[2]);
						
						Double[][][][] neightbour = (tokens[1].equals("left")) ? this.leftNeightbour : this.rightNeightbour;
	
						neightbour[AminoAcid.typeToInt(aminoAcid)][AminoAcid.typeToInt(neightbourAcid)][toBin(phi)][toBin(psi)] = probability;
					}
				}
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public double probability(AminoAcid.Type aminoAcid, AminoAcid.Type leftNeighbour, AminoAcid.Type rightNeighbour, double phi, double psi) {
		double left = this.leftNeightbour[AminoAcid.typeToInt(aminoAcid)][AminoAcid.typeToInt(leftNeighbour)][toBin(radianToDegree(phi))][toBin(radianToDegree(psi))];
		double right = this.rightNeightbour[AminoAcid.typeToInt(aminoAcid)][AminoAcid.typeToInt(rightNeighbour)][toBin(radianToDegree(phi))][toBin(radianToDegree(psi))];
		
		return (left * right) / this.probability(aminoAcid, phi, psi);
	}
	
	public double probability(AminoAcid.Type aminoAcid, double phi, double psi) {		
		return this.singleProbabilities[AminoAcid.typeToInt(aminoAcid)][toBin(radianToDegree(phi))][toBin(radianToDegree(psi))];
	}
	
	public Tuple2<Double,Double> purposeAngle(AminoAcid.Type aminoAcid, AminoAcid.Type leftNeighbour, AminoAcid.Type rightNeighbour) {
		// http://en.wikipedia.org/wiki/Rejection_sampling
		double phi = 0;
		double psi = 0;
		double p;
		
		do {
			phi = (Math.random() - 0.5) * Math.PI * 2;
			psi = (Math.random() - 0.5) * Math.PI * 2;
			
			p = this.probability(aminoAcid, leftNeighbour, rightNeighbour, phi, psi);
			
		} while(Math.random() >= p);
		
		return new Tuple2<Double,Double>(phi, psi);
	}
	
	public Tuple2<Double,Double> purposeAngle(AminoAcid.Type aminoAcid) {
		// http://en.wikipedia.org/wiki/Rejection_sampling
		double phi = 0;
		double psi = 0;
		double p;
		
		do {
			phi = (Math.random() - 0.5) * Math.PI * 2;
			psi = (Math.random() - 0.5) * Math.PI * 2;
			
			p = this.probability(aminoAcid, phi, psi);
			
		} while(Math.random() >= p);

		return new Tuple2<Double,Double>(phi, psi);
	}
	
	private static double radianToDegree(double radian) {
		return radian * RADIAN_TO_DEGREE_FACTOR;
	}
	
	private static int toBin(double angle) {
		return (int) (angle + 180)/5;
	}
}
