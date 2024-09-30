import java.util.*;
import java.math.BigInteger;

class Rationnel{
    BigInteger num;
    BigInteger den;

    static BigInteger ZERO=BigInteger.ZERO;
    static BigInteger ONE=BigInteger.ONE;
    
    public Rationnel(BigInteger b1, BigInteger b2){
	num=b1;
	den=b2;
    }
		     
    //Additionne 2 nombres rationnels
    public static Rationnel addition(Rationnel r1, Rationnel r2){
	Rationnel re =  new Rationnel(r1.num.multiply(r2.den).add(r1.den.multiply(r2.num)), r1.den.multiply(r2.den));
	re.reduit();
	return re;
    }

    //Soustrait 2 nombres rationnels  
    public static Rationnel soustraction(Rationnel r1,Rationnel r2){
        Rationnel re= new Rationnel(r1.num.multiply(r2.den).subtract(r1.den.multiply(r2.num)), r1.den.multiply(r2.den));
	re.reduit();
	return re;
    }

    //Multiplie 2 nombres rationnels
    public static Rationnel multiplication(Rationnel r1,Rationnel r2){
	Rationnel r=new Rationnel(r1.num.multiply(r2.num), r1.den.multiply(r2.den));
	r.reduit();
    	return r;
    }

    //Divise 2 nombres rationnels
    public static Rationnel division(Rationnel r1,Rationnel r2){
	Rationnel r=new Rationnel(r1.num.multiply(r2.den), r1.den.multiply(r2.num));
	r.reduit();
	return r;
    }

    //Donne l'inverse d'un nombre rationnel
    public static Rationnel inverse(Rationnel r){
    	return new Rationnel(r.den, r.num);
    }

    //Verifie si 2 nombres rationnels sont egaux
    public boolean equals(Rationnel r){
	return num.multiply(r.den).equals(den.multiply(r.num));
    }

    //Reduit le numerateur et le denominateur
    void reduit(){
    	if(num.equals(ZERO)){
	    den=ONE;
	    return;
    	}
	BigInteger pgcd = Arithmetique.pgcdBIG(num,den);
	num = num.divide(pgcd);
	den = den.divide(pgcd);
    }
    
    public static Rationnel[] somme(Rationnel[] p1, Rationnel[] p2){
	Rationnel[] ans;
	
    	if(p1.length!=p2.length){
	    if(p1.length>p2.length){
		Rationnel[] p=p1;
		p1=p2;
		p2=p;
	    }
	    ans = new Rationnel [p2.length];
	}else{
	    int count = p1.length;
	    for(int i=count-1;i>-1;i--){
		if(addition(p1[i], p2[i]).num.equals(ZERO)){
		    count--;
		}else{
		    break;
		}
	    }
	    ans = new Rationnel[count];
	}
	int i;
	for(i=0;i<ans.length&&i<p1.length;i++){
	    ans[i]= addition(p1[i], p2[i]);
	}
	for(; i<ans.length; i++){
	    ans[i]=p2[i];
	}
	return ans;
    }

    public static Rationnel[] multiplication(Rationnel[] p1, Rationnel[] p2){
	Rationnel[] ans;

	if(p1.length == 0 || p2.length == 0){
	    ans = new Rationnel[0];
	}else{
	    ans = new Rationnel[p1.length+p2.length-1];
	}
	for(int i=0; i<ans.length; i++)
	    ans[i]=new Rationnel(ZERO, ONE);
	for(int i=0;i<p1.length;i++){
	    for(int j=0;j<p2.length;j++){
		ans[i+j]=addition(multiplication(p1[i], p2[j]), ans[i+j]);
	    }
	}
	return ans;
    }
    
    static Rationnel[] soustraction(Rationnel[] p1, Rationnel[] p2){
    	Rationnel[] tmp=new Rationnel[p2.length];
    	for(int i=0; i<tmp.length; i++)
	    tmp[i]=new Rationnel(p2[i].num.negate(), p2[i].den);
    	return somme(p1, tmp);
    }
    
    static Rationnel[][] division(Rationnel[] dividende, Rationnel[] diviseur){
	if(dividende.length<diviseur.length)
	    return new Rationnel[][]{new Rationnel[0], dividende};
	Rationnel[] quotient=new Rationnel[dividende.length-diviseur.length+1];
	for(int i=0; i<quotient.length; i++)
	    quotient[i]=new Rationnel(ZERO, ONE);
	Rationnel[] reste=dividende;
	while(reste.length>=diviseur.length){
	    quotient[reste.length-diviseur.length]=multiplication(reste[reste.length-1],inverse(diviseur[diviseur.length-1]));
	    reste=soustraction(dividende, multiplication(diviseur, quotient));
	}
	return new Rationnel[][]{quotient, reste};
    }
    
    // l'un des polynomes est unitaire
    static BigInteger[] pgcd(BigInteger[] t1, BigInteger[] t2){
	Rationnel[] p1=new Rationnel[t1.length];
	Rationnel[] p2=new Rationnel[t2.length];
    
	for(int i=0; i<t1.length; i++)      
	    p1[i]=new Rationnel(t1[i], ONE);
	for(int i=0; i<t2.length; i++)
	    p2[i]=new Rationnel(t2[i], ONE);
	return pgcdAux(p1, p2);
    }
    
    static BigInteger[] pgcdAux(Rationnel[] p1, Rationnel[] p2){
    	if(p2.length>p1.length){
    	    Rationnel[] p = p2;
    	    p2 = p1;
    	    p1 = p;
	}
	Rationnel[][] divis = division(p1,p2);
    	if(divis[1].length==0){
    	    BigInteger[] ans=new BigInteger[p2.length];
    	    Rationnel elt;
    	    Rationnel inv=inverse(p2[p2.length-1]);
    	    for(int i=0; i<p2.length; i++){
    	    	elt=multiplication(inv, p2[i]);
    	    	ans[i]=elt.num.divide(elt.den);
	    }
	    return ans;
	}
	return pgcdAux(p2,divis[1]);
    }
    
    //Retourne le determinant d'une matrice
    static Rationnel determinant(Rationnel[][] matrix){
    	int len=matrix.length;
    	Rationnel ans=new Rationnel(ONE, ONE);
    	Rationnel coeff;
    	for(int j=0; j<len; j++){
	    for(int i=j; i<len; i++)
		if(!matrix[i][j].num.equals(ZERO)){
		    Rationnel tmp;
		    if(i!=j)
			ans.num=ans.num.negate();
		    for(int k=j; k<len; k++){
		       	tmp=matrix[j][k];
    			matrix[j][k]=matrix[i][k];
    			matrix[i][k]=tmp;
		    }
		    break;
		}
	    if(matrix[j][j].num.equals(ZERO))
		return new Rationnel(ZERO, ONE);
	    ans=multiplication(ans, matrix[j][j]);
	    coeff=inverse(matrix[j][j]);
	    for(int k=j; k<len; k++)
		matrix[j][k]=multiplication(matrix[j][k], coeff);
	    for(int i=0; i<len; i++){
		if(i==j)
		    continue;
	        coeff=matrix[i][j];
		for(int k=j; k<len; k++)
	            matrix[i][k]=soustraction(matrix[i][k], multiplication(coeff, matrix[j][k]));
            }   
        }
        return ans;
    }
}
