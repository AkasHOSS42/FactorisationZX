import java.util.*;
import java.math.BigInteger;

/* Un polynome a coefficients dans Z/nZ.
 * Certaines methodes n'ont de sens que dans Fp. */
class PolyMod{
    static BigInteger ZERO=BigInteger.ZERO;
    static BigInteger ONE=BigInteger.ONE;
	
    BigInteger p; // dans quel corps sont les coefficients
    BigInteger[] coeff;
    
    static BigInteger toBI(int n){
	return new BigInteger(String.valueOf(n));
    }
    
    public PolyMod(BigInteger val,int[] tab){
    	p =val;
    	int len;
    	for(len=tab.length; len!=0&&tab[len-1]%p.intValue()==0; len--){}
	coeff =new BigInteger[len];
    	for(int i=0; i<len; i++)
	    coeff[i]=toBI(tab[i]).mod(p).add(p).mod(p);
    }

    public PolyMod(BigInteger val, BigInteger[] tab){
    	p =val;
    	int len;
    	for(len=tab.length; len!=0&&tab[len-1].mod(p).equals(ZERO); len--){}
	coeff =new BigInteger[len];
    	for(int i=0; i<len; i++)
	    coeff[i]=tab[i];
	reduit();
    }

    /* reduit les coefficients modulo p,
     * pour eviter de manipuler des Biginteger trop grands. */
    void reduit(){
    	for(int i=0; i<coeff.length; i++)
	    coeff[i]=coeff[i].mod(p).add(p).mod(p);
    }
    
    //Additionne 2 polynomes
    public static PolyMod somme(PolyMod p1, PolyMod p2){
    	BigInteger[] tab;
    	if(p1.coeff.length!=p2.coeff.length){
	    if(p1.coeff.length>p2.coeff.length){
		PolyMod p=p1;
		p1=p2;
		p2=p;
	    }
	    tab = new BigInteger[Math.max(p1.coeff.length,p2.coeff.length)];
    	}else{
	    int count = p1.coeff.length;
	    for(int i=count-1;i>-1;i--){
		if((p1.coeff[i].add(p2.coeff[i])).mod(p1.p).equals(ZERO)){
		    count--;
		}else{
		    break;
		}
	    }
	    tab = new BigInteger[count];
    	}
    	int i;
    	for(i=0;i<tab.length&&i<p1.coeff.length;i++){
	    tab[i]= p1.coeff[i].add(p2.coeff[i]);
    	}
    	for(; i<tab.length; i++){
	    tab[i]=p2.coeff[i];
    	}
    	PolyMod ans= new PolyMod(p1.p, tab);
    	ans.reduit();
    	return ans;
    }
    
    //Multiplie 2 polynomes
    public static PolyMod multiplication(PolyMod p1, PolyMod p2){
    	BigInteger[] tab;
    	if(p1.coeff.length == 0 || p2.coeff.length == 0){
	    tab = new BigInteger[0];
    	}else{
	    tab = new BigInteger[p1.coeff.length+p2.coeff.length-1];
    	}
    	for(int i=0; i<tab.length; i++)
	    tab[i]=ZERO;
    	for(int i=0;i<p1.coeff.length;i++){
	    for(int j=0;j<p2.coeff.length;j++){
		tab[i+j]=tab[i+j].add(p1.coeff[i].multiply(p2.coeff[j]));		
	    }
    	}
    	PolyMod ans= new PolyMod(p1.p, tab);
	ans.reduit();
	return ans;
    }
    
    //Multiplie un polynome avec un entier
    void multiplie(int n){
    	for(int i=0; i<coeff.length; i++)
	    coeff[i]=coeff[i].multiply(toBI(n));
    	reduit();
    }
    
    //Soustrait 2 polynomes
    static PolyMod soustraction(PolyMod p1, PolyMod p2){
    	p2.multiplie(-1);
    	PolyMod ans=somme(p1, p2);
    	p2.multiplie(-1);
    	return ans;
    }
    
    //Divise 2 polynomes. renvoie le couple {quotient, reste}
    static PolyMod[] division(PolyMod dividende, PolyMod diviseur){
    	if(dividende.coeff.length<diviseur.coeff.length)
	    return new PolyMod[]{new PolyMod(dividende.p, new int[0]), dividende};
    	BigInteger[] tab=new BigInteger[dividende.coeff.length-diviseur.coeff.length+1];
	for(int i=0; i<tab.length; i++)
	    tab[i]=ZERO;
    	PolyMod quotient=new PolyMod(dividende.p, tab);
    	quotient.coeff=tab;
    	PolyMod reste=dividende;
    	while(reste.coeff.length>=diviseur.coeff.length){
	    quotient.coeff[reste.coeff.length-diviseur.coeff.length]=reste.coeff[reste.coeff.length-1].multiply(Arithmetique.inverse(diviseur.p, diviseur.coeff[diviseur.coeff.length-1]));
	    reste=soustraction(dividende, multiplication(diviseur, quotient));
    	}
    	return new PolyMod[]{quotient, reste};
    }
    
    //Trouve le pgcd de 2 polynomes avec Euclide
    static PolyMod pgcd(PolyMod p1, PolyMod p2){
    	if(p2.coeff.length>p1.coeff.length){
	    PolyMod p = p2;
	    p2 = p1;
	    p1 = p;
    	}
    	PolyMod [] divis = division(p1,p2);
    	if(divis[1].coeff.length==0){
	    p2.multiplie(Arithmetique.inverse(p2.p, p2.coeff[p2.coeff.length-1]).intValue());
	    return p2;
    	}
    	return pgcd(p2,divis[1]);
    }

    /* Troucve u et v tels que :
     * up1+vp2=1
     * degu<degp2, degv<degp1 */
    static PolyMod[] trouveBezout(PolyMod p1, PolyMod p2){
    	PolyMod[] div=division(p1, p2);
    	if(div[1].coeff.length==1){
	    PolyMod ans0=new PolyMod(p1.p, new BigInteger[]{Arithmetique.inverse(p1.p, div[1].coeff[0])});
	    div[0].multiplie(-1);
	    return new PolyMod[]{ans0, multiplication(ans0, div[0])};
    	}
    	div=division(p2, p1);
    	if(div[1].coeff.length==1){
	    PolyMod ans0=new PolyMod(p1.p, new BigInteger[]{Arithmetique.inverse(p2.p, div[1].coeff[0])});
	    div[0].multiplie(-1);
	    return new PolyMod[]{multiplication(ans0, div[0]), ans0};
    	}
	int len1=p1.coeff.length;
	int len2=p2.coeff.length;
	int len=len1+len2-2;
	/* On cherche les solutions d'un systeme
	 * lineaire dont les inconnues sont les coefficients
	 * de u et v. */
	int[][] mat=new int[len][len+1];
	for(int j=0; j<len2-1; j++)
	    for(int i=j; i<j+len1; i++)
		mat[i][j]=p1.coeff[p1.coeff.length-1+j-i].intValue();
	for(int j=0; j<len1-1; j++)
	    for(int i=j; i<j+len2; i++)
		mat[i][j+len2-1]=p2.coeff[p2.coeff.length-1+j-i].intValue();
	mat[len-1][len]=1;
	Arithmetique.pivot(p1.p.intValue(), mat);
	int nbParams=Factorisation.nbFact(mat);
	int[] pivots=new int[len-nbParams];
	for(int i=0; i<pivots.length; i++)
	    for(int j=0; j<len; j++)
		if(mat[i][j]!=0){
		    pivots[i]=j;
		    break;
		}
	int[] ans1=new int[len2-1];
	int[] ans2=new int[len1-1];
	int rank=0;
	for(int i=0; i<pivots.length; i++){
	    if(pivots[i]<len2-1)
		ans1[ans1.length-1-pivots[i]]=mat[i][len];
	    else
		ans2[ans2.length-1-pivots[i]+len2-1]=mat[i][len];
	}
	return new PolyMod[]{new PolyMod(p1.p, ans1), new PolyMod(p2.p, ans2)};
    }
}
