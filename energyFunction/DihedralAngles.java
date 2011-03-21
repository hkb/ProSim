package energyFunction;

import java.util.List;

import dataStructure.ChainTree;

public class DihedralAngles implements EnergyFunction {
	
	private ChainTree testing, target;
	private List<Double> targetAngles;

	public DihedralAngles(ChainTree testing, ChainTree target) {
		this.testing = testing;
		this.target = target;
		
		this.targetAngles = target.getDihedralAngles();
	}
	
	@Override
	public double compute() {
		List<Double> testingAngles = this.testing.getDihedralAngles();
		
		double sum = 0;
		
		for (int i = 0, j = testingAngles.size(); i < j; i++) {
			double diff = this.getAngleDiff(testingAngles.get(i), this.targetAngles.get(i));
			
			sum += diff * diff;
		}
		
		return Math.sqrt(sum / testingAngles.size());
	}
	
	private double getAngleDiff(double a, double b) {
		double pi = Math.PI;
		double pipi = pi * pi;
		
		while (a > pi) a -= pipi;
		while (a <= -pi) a += pipi;
		while (b > pi) b -= pipi;
		while (b <= -pi) b += pipi;
		if (a*b >= 0.0) return Math.abs(a-b);
		else {
			double sum = Math.abs(a) + Math.abs(b);
			if (sum <= pi) return sum; else return pipi - sum;
		}
	}
}
