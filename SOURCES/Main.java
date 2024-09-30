class Main{
    public static void main(String[] args){
	int[] pol=new int[args.length];
	for(int i=0; i<pol.length; i++)
	    pol[i]=Integer.parseInt(args[i]);
	if(pol.length==0){
	    System.out.println("0");
	    return;
	}
	if(pol[pol.length-1]==0){
	    System.err.println("Erreur : coefficient dominant nul.");
	    return;
	}
	System.out.println(Factorisation.factorise_quelconque(pol).toString());
    }
}
