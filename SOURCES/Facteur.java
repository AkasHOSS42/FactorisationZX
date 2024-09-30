import java.math.BigInteger;

class Facteur{
    static BigInteger ZERO=BigInteger.ZERO;
    static BigInteger ONE=BigInteger.ONE;
	
    int multiplicite;
    BigInteger[] pol; //les coefficients du facteur
    Facteur(BigInteger[] q, int m){
	multiplicite=m;
	pol=q;
    }
	
    public String toString(){
	String poly = "(";
		
	if(pol.length ==0){
	    return "0";
	}
		
	for(int i=0;i<pol.length;i++){
	    if(pol[i].equals(ZERO))
		continue;
	    if(i==0||!pol[i].equals(ONE))
		poly =poly+ pol[i].toString();
	    if(i!=0){
		poly = poly+ "X";
		if(i!=1)
		    poly=poly+"^"+String.valueOf(i);
	    }
	    if(i!=pol.length-1)
		poly=poly+" + ";
	}
	poly=poly+")";
	if(multiplicite>1)
	    poly=poly+"^"+multiplicite;
	return poly;
    }
}
