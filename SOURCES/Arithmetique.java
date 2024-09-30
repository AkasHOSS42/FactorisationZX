import java.math.BigInteger;

class Arithmetique{
    static BigInteger ZERO=BigInteger.ZERO;
    static BigInteger ONE=BigInteger.ONE;

    static BigInteger toBI(int n){
	return new BigInteger(String.valueOf(n));
    }
    
    //Trouve une ralation de Bezout pour deux entiers
    private static BigInteger[] trouveBezout(BigInteger n1, BigInteger n2){
    	if(n2.equals(ONE))
	    return new BigInteger[]{ONE, ONE.subtract(n1)};
    	BigInteger r=n1.mod(n2);
    	if(r.equals(ONE))
	    return new BigInteger[]{ONE, n1.divide(n2).negate()};
    	BigInteger q=n1.divide(n2);
    	BigInteger[] tab=trouveBezout(n2, r);
    	return new BigInteger[]{tab[1], tab[0].subtract(q.multiply(tab[1]))};
    }
    
    //Retourne l'inverse de n dans Fp
    static BigInteger inverse(BigInteger p, BigInteger n){
    	return trouveBezout(p, n.mod(p).add(p).mod(p))[1].mod(p);
    }
	
    //Echange la place des lignes i1 et i2 dans tab
    private static void swap(int i1, int i2, int[][]tab){
    	int tmp;
    	for(int j=0; j<tab[0].length; j++){
	    tmp=tab[i1][j];
	    tab[i1][j]=tab[i2][j];
	    tab[i2][j]=tmp;
    	}
    }
    
    //Fait un pivot dans une matrice donnee sur Fp
    private static void pivote(int iPivot, int iTarget, int j, int[][]matrix, int p){
    	int mult=-matrix[iTarget][j];
    	for(int j2=0; j2<matrix[0].length; j2++)
	    matrix[iTarget][j2]=(matrix[iTarget][j2]+mult*matrix[iPivot][j2])%p;
    }
    
    //Fait un pivot de Gauss dans une matrice sur Fp
    static void pivot(int p, int[][]matrix){
    	int len=matrix.length;
    	int j=0;
    	int i;
    	int iMin=0;
    	int inverse;
    	while(j<len&&iMin<len){
	    for(i=iMin; i<len; i++)
		if(matrix[i][j]%p!=0){
		    swap(i, iMin, matrix);
		    inverse=inverse(toBI(p), toBI(matrix[iMin][j]%p)).intValue();
		    for(int j2=j; j2<matrix[0].length; j2++)
			matrix[iMin][j2]*=inverse;
		    for(i=0; i<len; i++)
			if(i!=iMin)
			    pivote(iMin, i, j, matrix, p);
		    iMin++;
		    break;
		}
	    j++;
    	}
    }
    
    /* Approxime la racine carree de n. */
    static BigInteger sqrt(BigInteger n){
	BigInteger b1=ZERO.setBit(n.bitLength()/2);
	BigInteger b2=b1;
	BigInteger b3;
	for(;true;){
	    b3=b1.add(n.divide(b1)).shiftRight(1);
	    if(!((!b3.equals(b1))&&!b3.equals(b2)))
		return b3;
	    b2=b1;
	    b1=b3;
	}
    }
    
    static int binomial(int n, int k){
    	int ans1=1;
    	int ans2=1;
    	for(int i=0; i<n-k; i++){
	    ans1*=k+i+1;
	    ans2*=i+1;
    	}
	int ans=ans1/ans2;
    	return ans;
    }
    
    static boolean isPrime(int p){
    	int max=(int)(Math.ceil(Math.sqrt(p)));
    	for(int i=2; i<=max; i++)
	    if(p%i==0)
    		return false;
    	return true;
    }
    
    /* Retourne un nombre premier plus grand que p. */
    static int nextPrime(int p){
    	int ans=p+2;
    	while(!isPrime(ans))
	    ans+=2;
    	return ans;
    }
    
    public static BigInteger pgcdBIG(BigInteger nb1, BigInteger nb2){
    	if(nb1.compareTo(ZERO)<0)
	    nb1=nb1.negate();
    	if(nb2.compareTo(ZERO)<0)
	    nb2=nb2.negate();
	if(nb1.compareTo(nb2)<0){
	    BigInteger nb = nb1;
	    nb1 =nb2;
	    nb2 = nb;
	}
	BigInteger reste = nb1.mod(nb2);
	if(reste.equals(ZERO))
	    return nb2;
	return pgcdBIG(nb2,reste);
    }

    public static int pgcdEntier(int nb1,int nb2){
    	if(nb1<0)
	    nb1=-nb1;
    	if(nb2<0)
	    nb2=-nb2;
	if(nb1<nb2){
	    int nb = nb1;
	    nb1 =nb2;
	    nb2 = nb;
	}
	int reste = nb1%nb2;
	if(reste==0)
	    return nb2;
	return pgcdEntier(nb2,reste);
    }
}
