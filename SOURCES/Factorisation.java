import java.math.BigInteger;
import java.util.*;

class Factorisation{
    static BigInteger ZERO=BigInteger.ZERO;
    static BigInteger ONE=BigInteger.ONE;

    static BigInteger toBI(int n){
	return new BigInteger(String.valueOf(n));
    }
	
    LinkedList<Facteur> facteurs;
    int multiple;
	
    public String toString(){
	String ans=String.valueOf(multiple);
	for(Facteur f : facteurs)
	    ans=ans+f.toString();
	return ans;
    }

    // cf section 3.3
    // prev={B, C, U, V}
    static BigInteger[][] remonte_HENSEL(BigInteger[][] prev, BigInteger[] prod, BigInteger p){
    	PolyMod[] tmp=new PolyMod[4];
    	for(int i=0; i<4; i++)
    	    tmp[i]=new PolyMod(p.multiply(p), prev[i]);
    	PolyMod multp=PolyMod.soustraction(new PolyMod(p.multiply(p), prod), PolyMod.multiplication(tmp[0], tmp[1]));
    	PolyMod[] div=PolyMod.division(PolyMod.multiplication(multp, tmp[2]), tmp[1]);
        BigInteger[][] ans=new BigInteger[4][];
    	PolyMod ans0=PolyMod.somme(PolyMod.somme(tmp[0], PolyMod.multiplication(tmp[3], multp)), PolyMod.multiplication(div[0], tmp[0]));
    	PolyMod ans1=PolyMod.somme(tmp[1], div[1]);
        PolyMod aux=PolyMod.somme(PolyMod.multiplication(tmp[2], ans0), PolyMod.multiplication(tmp[3], ans1));
        aux.coeff[0]=aux.coeff[0].subtract(ONE);
        PolyMod[] div2=PolyMod.division(PolyMod.multiplication(aux, tmp[2]), ans1);
    	ans[0]=ans0.coeff;
    	ans[1]=ans1.coeff;
        ans[2]=PolyMod.soustraction(tmp[2], div2[1]).coeff;
        aux.multiplie(-1);
    	aux.coeff[0]=aux.coeff[0].add(ONE);
        ans[3]=PolyMod.soustraction(PolyMod.multiplication(tmp[3], aux), PolyMod.multiplication(div2[0], ans0)).coeff;
        return ans;
    }
    
    /* poly est sans facteurs multiples. */
    static LinkedList<BigInteger[]> factorise_HENSEL(BigInteger[] poly){
    	BigInteger res=PolyZ.resultant(poly, PolyZ.derivation(poly));
    	int p=3;
    	while(res.mod(toBI(p)).equals(ZERO))
	    p=Arithmetique.nextPrime(p);
	return factorise_HENSEL_rec(poly, p);
    }

    static LinkedList<BigInteger[]> factorise_HENSEL_rec(BigInteger[] poly, int p){
	LinkedList<BigInteger[]> ans=new LinkedList<BigInteger[]>();
	PolyMod pol=new PolyMod(toBI(p), poly);
	ArrayList<PolyMod> factP=Berlekamp(pol);
	int len=factP.size();
	if(len==1){ //poly est irreductible
	    ans.add(poly);
	    return ans;
	}
	    
	for(PolyMod elt : factP)
	    elt.multiplie(Arithmetique.inverse(elt.p, elt.coeff[elt.coeff.length-1]).intValue());
	BigInteger borne=PolyZ.majore_coeff_facteurs(poly, 1+poly.length/2);
	BigInteger current;
	BigInteger[][] data=new BigInteger[4][];
	ArrayList<PolyMod> facteurs_modulaires=new ArrayList<PolyMod>();
	PolyMod tmp;
	
	/* On remonte tous les facteurs modulo p^k */
    	for(PolyMod elt : factP){
    	    data[1]=elt.coeff;
    	    tmp=null;
    	    for(PolyMod elt2 : factP)
		if(elt2!=elt){
		    if(tmp==null)
			tmp=elt2;
		    else
			tmp=PolyMod.multiplication(tmp, elt2);
		}
    	    data[0]=tmp.coeff;
    	    PolyMod[] bezout=PolyMod.trouveBezout(elt, tmp);
    	    data[3]=bezout[0].coeff;
    	    data[2]=bezout[1].coeff;
    	    current=toBI(p);
    	    while(current.compareTo(borne)<0){
		data=remonte_HENSEL(data, poly, current);
		current=current.multiply(current);
    	    }
    	    facteurs_modulaires.add(new PolyMod(current, data[1]));
    	}
    	BigInteger[] diviseur=teste_TOUT(poly, facteurs_modulaires, facteurs_modulaires.size(), 0, null);
    	if(diviseur==null){ // poly est irreductible
    	    ans.add(poly);
    	    return ans;
    	}
    	LinkedList<BigInteger[]> l1=factorise_HENSEL_rec(diviseur, p);
    	LinkedList<BigInteger[]> l2=factorise_HENSEL_rec(PolyZ.division(poly, diviseur)[0], p);
    	for(BigInteger[] tab : l1)
    	    ans.add(tab);
    	for(BigInteger[] tab : l2)
    	    ans.add(tab);
    	return ans;
    }

    /* Methode recursive qui teste tous les sous-produits de facteurs,
     * a la recherche d'un facteur non trivial de poly.
     * index est l'indice du facteur de facteurs qu'on est en train de manipuler.
     * candidat est le sous produit qu'on est en train d'examiner et qu'on construit en traversant facteurs */
    static BigInteger[] teste_TOUT(BigInteger[] poly, ArrayList<PolyMod> facteurs, int len, int index, PolyMod candidat){
	if(len==index){
	    if(candidat==null)
		return null;
	    candidat.reduit();
	    for(int i=0; i<candidat.coeff.length; i++)
		if(candidat.coeff[i].compareTo(candidat.p.divide(toBI(2)))>0)
		    candidat.coeff[i]=candidat.coeff[i].subtract(candidat.p);
	    BigInteger[][] div=PolyZ.division(poly, candidat.coeff);
	    if(div[1].length==0&&div[0].length!=1) // candidat est un diviseur && il est non trivial.
		return candidat.coeff;
	    return null;
	}
	BigInteger[] ans=teste_TOUT(poly, facteurs, len, index+1, candidat);
	if(ans!=null)
	    return ans;
	if(candidat==null)
	    candidat=facteurs.get(index);
	else
	    candidat=PolyMod.multiplication(candidat, facteurs.get(index));
	return teste_TOUT(poly, facteurs, len, index+1, candidat);
    }
    
    static LinkedList<Facteur> factorise_unitaire(BigInteger[] poly){
        BigInteger[] polyDeriv = PolyZ.derivation(poly);
        BigInteger[] pgcd=Rationnel.pgcd(poly, polyDeriv);
        if(pgcd.length==1){
	    LinkedList<Facteur> ans= new LinkedList<Facteur>();
	    LinkedList<BigInteger[]> facts=factorise_HENSEL(poly);
	    for(BigInteger[] fact : facts)
		ans.add(new Facteur(fact, 1));
	    return ans;
        }
        LinkedList<Facteur> l1=factorise_unitaire(pgcd);
        BigInteger[] poly2=new BigInteger[]{ONE};
        for(Facteur f : l1){
	    f.multiplicite++;
	    poly2=PolyZ.multiplication(poly2, PolyZ.exponentiation(f.pol, f.multiplicite));
        }
        BigInteger[] quotient=PolyZ.division(poly, poly2)[0];
        if(quotient.length<2)
	    return l1;
        LinkedList<BigInteger[]> l2=factorise_HENSEL(quotient);
        for(BigInteger[] f : l2)
	    l1.add(new Facteur(f, 1));
        return l1;
    }
    
    static Factorisation factorise_quelconque(int [] poly){
    	Factorisation ans=new Factorisation();
    	int contenu=PolyZ.contenu(poly);
    	for(int i =0; i< poly.length; i++)
	    poly[i] /= contenu;
    	ans.multiple=contenu;
    	ans.facteurs=factorise_primitif(poly);
    	return ans;
    }
    
    static LinkedList<Facteur> factorise_primitif(int[] poly){
    	BigInteger[] polyPrim = new BigInteger[poly.length];
    	BigInteger an = toBI(poly[poly.length-1]);
    	int i,j;
    	BigInteger ann_1 = ONE;
    	for(i=poly.length-2;i>-1;i--){
	    polyPrim[i] = toBI(poly[i]).multiply(ann_1);
	    ann_1=ann_1.multiply(an);
    	}
    	polyPrim[poly.length-1] = ONE;
    	
    	LinkedList<Facteur> listfact = factorise_unitaire(polyPrim);
    	BigInteger contFact;
    	for(Facteur f : listfact){
	    BigInteger anN = ONE;
	    for(i = 0;i<f.pol.length;i++){
		f.pol[i] =f.pol[i].multiply(anN);
		anN=anN.multiply(an);
	    }
	    contFact = PolyZ.contenu(f.pol);
	    for(i=0;i<f.pol.length;i++){
		f.pol[i] =f.pol[i].divide(contFact);
	    }
    	}
        return listfact;
    }
    
    /* Calcule la matrice de l'application lineaire
     * Q->Q^p-Q, dans la base des puissances de X dans l'ordre decroissant.
     * Echelonne et reduit cette matrice avec le pivot de Gauss. */
    static int[][] aux(PolyMod poly){
    	int tailleMat = poly.coeff.length-1;
    	int [][]mat = new int [tailleMat][tailleMat];
    	int p=poly.p.intValue();
    	for(int i=0;i<tailleMat-1;i++){
	    mat[i][tailleMat-1] = 0;
    	}
    	mat[tailleMat-1][tailleMat-1] =1;
    	int count = 2;
    	int deg;
    	int [] tab;
    	//permet d'avoir la matrice
    	for(int i = 1; i<tailleMat;i++){
	    deg = i*p;
	    tab =new int [deg+1]; //represente X^i a  faire diviser
	    for(int j=0;j<deg;j++){
		tab[j]=0;
	    }
	    tab[deg] = 1;
	    PolyMod divid = new PolyMod(poly.p,tab);
	    PolyMod reste = PolyMod.division(divid,poly)[1];
	    for(int j=reste.coeff.length-1; j>-1;j--)
		mat[tailleMat-1-j][tailleMat-i-1] = reste.coeff[j].mod(poly.p).intValue();
    	}
    	for(int i=0;i<tailleMat;i++){
	    mat[i][i]--;
    	}
    	Arithmetique.pivot(p, mat);
    	return mat;
    }
    
    //Retourne le nombre de facteur d'un polynome
    static int nbFact(int[][] mat){
    	int ans=0;
    	int len=mat.length;
    	boolean quit=false;
    	while(!quit&&ans<len){
	    for(int i=0; i<len; i++)
		if(mat[len-ans-1][i]!=0){
		    quit=true;
		    break;
		}
	    if(!quit)
		ans++;
    	}
    	return ans;
    }
    
    /* Renvoie un couple de deux polynomes, non associes a poly ou 1,
       et dont le produit vaut poly (une factorisation non triviale de poly). */
    static PolyMod[] coupe(int[][] mat, PolyMod poly, int nbFact){
        PolyMod pol;
        /* il faut que pol soit un polynome non constant du noyau de mat.
	 * comment faire pour le trouver ?
	 * trouver les parametres (les colonnes qui ne contiennent pas de pivots)
	 * chaque parametre correspond a un vecteur d'une base du noyau
	 * un de ces vecteurs est forcement non constant (ie l'un de ses coeffs de degre !=0 est !=0 modulo p)
	 * poly va etre le polynome qui correspond a ce vecteur.
	 * rappel : on avait dit que la derniere ligne (et colonne) de mat correspond a 1
	 * et les autres aux X^i */
        
        int[] pivots=new int[mat.length-nbFact];
	int p=poly.p.intValue();
        for(int i=0; i<mat.length-nbFact; i++)
	    for(int j=0; true; j++)
		if(mat[i][j]%p!=0){
		    pivots[i]=j;
		    break;
		}
	int param;
	int rank=0;
	for(param=0; true; param++){
	    if(rank<pivots.length&&pivots[rank]==param){
		rank++;
		continue;
	    }else
		break;
	}
	/* param induit alors un vecteur non constant du noyau */
	int[] tmp=new int[mat.length];
	int deg=0;
	rank=0;
        	
	for(int i=0; i<mat.length; i++){
    	    if(rank<pivots.length&&pivots[rank]==i){
		rank++;
    		if(mat[rank-1][param]%p==0)
    		    tmp[tmp.length-1-i]=0;
    		else
    		    tmp[tmp.length-1-i]=mat[rank-1][param]%p;
    	    }
    	    else if(i==param)
		tmp[tmp.length-1-i]=-1;
    	    else
		tmp[tmp.length-1-i]=0;
    	    if(tmp[tmp.length-1-i]!=0&&deg==0)
		deg=tmp.length-1-i;
	}
        	
	BigInteger []coeff=new BigInteger[deg+1];
	for(int i=0; i<deg+1; i++)
	    coeff[i]=toBI(tmp[i]);
	pol=new PolyMod(poly.p, coeff);
	PolyMod[] ans=new PolyMod[2];
	for(int i=0; i<p; i++){
	    ans[0]=PolyMod.pgcd(poly, pol);
	    if(ans[0].coeff.length>1 && ans[0].coeff.length!=poly.coeff.length)
		break;
	    pol.coeff[0]=pol.coeff[0].add(ONE);
	}
        	
	ans[1]=PolyMod.division(poly, ans[0])[0];
	return ans;
    }
    
    
    //Applique Berlekamp sur un polynome
    static ArrayList<PolyMod> Berlekamp(PolyMod poly){
    	ArrayList<PolyMod> list=new ArrayList<PolyMod>();
    	int [][] mat = aux(poly);
    	int nbFact  = nbFact(mat);
    	if(nbFact==1){
	    list.add(poly);
	    return list;
    	}
    	
    	PolyMod[] tab = coupe(mat,poly,nbFact);
    	ArrayList<PolyMod> list1 = Berlekamp(tab[0]);
    	ArrayList<PolyMod> list2 = Berlekamp(tab[1]);
    	
    	for(PolyMod p : list1)
	    list.add(p);
    	
    	for(PolyMod p : list2)
	    list.add(p);
    	
    	return list;
    }
}
