package chemestry;

public class AminoAcid {

	public static enum Type {ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, SEC, PYL};
	public static int count = 22;
	
	public static Type stringToType(String str) {
		return Type.valueOf(str);
	}
	
	public static int typeToInt(Type type) {
		switch(type) {
			case ALA: return 0;
			case ARG: return 1;
			case ASN: return 2;
			case ASP: return 3;
			case CYS: return 4;
			case GLU: return 5;
			case GLN: return 6;
			case GLY: return 7;
			case HIS: return 8;
			case ILE: return 9;
			case LEU: return 10;
			case LYS: return 11;
			case MET: return 12;
			case PHE: return 13;
			case PRO: return 14;
			case SER: return 15;
			case THR: return 16;
			case TRP: return 17;
			case TYR: return 18;
			case VAL: return 19;
			case SEC: return 20;
			case PYL: return 21;
		}
		
		throw new IllegalArgumentException("Unknown amoni acid type!");
	}
}
