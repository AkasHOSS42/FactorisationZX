import java.math.BigInteger;

/* Classe pour manipuler des polynomes a coefficients dans Z. */
class PolyZ{
    static BigInteger ZERO=BigInteger.ZERO;
    static BigInteger ONE=BigInteger.ONE;

    static BigInteger toBI(int n){
	return new BigInteger(String.valueOf(n));
    }
    
    public static BigInteger[] somme(BigInteger[] p1, BigInteger[] p2){
	BigInteger[] ans;
	if(p1.length!=p2.length){
	    if(p1.length>p2.length){
		BigInteger[] p=p1;
		p1=p2;
		p2=p;
	    }
	    ans = new BigInteger [p2.length];
	}else{
	    int count = p1.length;
	    for(int i=count-1;i>-1;i--){
		if(p1[i].add(p2[i]).equals(ZERO)){
		    count--;
		}else{
		    break;
		}
	    }
	    ans = new BigInteger[count];
	}
	int i;
	for(i=0;i<ans.length&&i<p1.length;i++){
	    ans[i]= p1[i].add(p2[i]);
	}
	for(; i<ans.length; i++){
	    ans[i]=p2[i];
	}
	return ans;
    }

    public static BigInteger[] multiplication(BigInteger[] p1, BigInteger[] p2){
    	BigInteger[] ans;
    	if(p1.length == 0 || p2.length == 0){
	    ans = new BigInteger[0];
    	}else{
	    ans = new BigInteger[p1.length+p2.length-1];
    	}
    	for(int i=0; i<ans.length; i++)
	    ans[i]=ZERO;
    	for(int i=0;i<p1.length;i++){
	    for(int j=0;j<p2.length;j++){
		ans[i+j]=ans[i+j].add(p1[i].multiply(p2[j]));
	    }
    	}
    	return ans;
    }
	
    static BigInteger[] exponentiation(BigInteger[] p, int n){
	if(n==1)
	    return p;
	BigInteger[] half=exponentiation(p, n/2);
	if(n%2==0)
	    return multiplication(half, half);
	return multiplication(half, multiplication(half, p));
    }
	
    static BigInteger[] soustraction(BigInteger[] p1, BigInteger[] p2){
	BigInteger[] tmp=new BigInteger[p2.length];
	for(int i=0; i<p2.length; i++)
	    tmp[i]=p2[i].negate();
	return somme(p1, tmp);
    }
	
    /* Une des deux conditions doit etre respectee :
     * 		diviseur est unitaire
     * ou
     * 		diviseur divise dividende (le reste est 0). */
    static BigInteger[][] division(BigInteger[] dividende, BigInteger[] diviseur){
	if(dividende.length<diviseur.length)
	    return new BigInteger[][]{new BigInteger[0], dividende};
	BigInteger[] quotient=new BigInteger[dividende.length-diviseur.length+1];
	for(int i=0; i<quotient.length; i++)
	    quotient[i]=ZERO;
	BigInteger[] reste=dividende;
	while(reste.length>=diviseur.length){
	    quotient[reste.length-diviseur.length]=reste[reste.length-1].divide(diviseur[diviseur.length-1]);
	    reste=soustraction(dividende, multiplication(diviseur, quotient));
	}
	return new BigInteger[][]{quotient, reste};
    }
    
    static BigInteger[] derivation(BigInteger[] p){
    	BigInteger[] ans=new BigInteger[p.length-1];
    	for(int i=1; i<p.length; i++)
	    ans[i-1]=PolyMod.toBI(i).multiply(p[i]);
    	return ans;
    }
	
    //Donne le contenu d'un polynome
    static int contenu(int[] poly){
    	int nb = 0;
    	if(poly.length ==0){
	    return 1;
    	}
    	for(int i=0;i<poly.length;i++){
	    if(poly[i] ==0) continue;
	    if(nb==0)
		nb=poly[i];
	    else
		nb = Arithmetique.pgcdEntier(nb,poly[i]);
    	}
    	return nb;
    }


    //Donne le contenu d'un polynome
    static BigInteger contenu(BigInteger[] poly){
    	BigInteger nb = ZERO;
    	if(poly.length ==0){
	    return ONE;
    	}
    	for(int i=0;i<poly.length;i++){
	    if(poly[i].equals(ZERO)) continue;
	    if(nb.equals(ZERO))
		nb=poly[i];
	    else
		nb = Arithmetique.pgcdBIG(nb,poly[i]);
    	}
    	return nb;
    }

    static BigInteger majore_norme_euclidienne(BigInteger[] tab){
    	BigInteger ans=ZERO;
    	for(int i=0; i<tab.length; i++){
	    ans=ans.add(tab[i].multiply(tab[i]));
	}
	return Arithmetique.sqrt(ans);
    }
    
    static BigInteger majore_coeff_facteurs(BigInteger [] tab, int degMax){
    	return toBI(2).multiply(majore_norme_euclidienne(tab)).multiply(toBI(Arithmetique.binomial(degMax, degMax/2)));
    }
    
    static BigInteger resultant(BigInteger[] p1, BigInteger[] p2){
    	int len=p1.length+p2.length-2;
    	Rationnel[][] matrix=new Rationnel[len][len];
    	for(int i=0; i<len; i++)
	    for(int j=0; j<len; j++)
		matrix[i][j]=new Rationnel(ZERO, ONE);
    	for(int j=0; j<p2.length-1; j++)
	    for(int i=j; i<j+p1.length; i++)
		matrix[i][j].num=p1[p1.length-1-i+j];
    	for(int j=0; j<p1.length-1; j++)
	    for(int i=j; i<j+p2.length; i++)
		matrix[i][j+p2.length-1].num=p2[p2.length-1-i+j];
    	Rationnel ans=Rationnel.determinant(matrix);
    	return ans.num.divide(ans.den);
    }
}
