/*********************************************
 * OPL 12.8.0.0 Model
 * Author: thibaut
 * Creation Date: 26 janv. 2020 at 09:50:47
 *********************************************/
int n=...; //10;
range rangeN=1..n;   
int L=...;  //2;
int W=...;  //3;
int K=...;  //9;
int B=...;  //48;
int w_v[rangeN]=...;//w_v = [4, 14, 3, 2, 20, 5, 14, 10, 10, 18];
float W_v[rangeN] =...; //[1.47239, 0.82384, 0.147042, 1.11452, 1.35826, 1.95929, 1.98969, 1.54618, 2.14362, 1.6074];
int lh[rangeN] = ...; //[11, 8, 8, 10, 8, 9, 10, 10, 8, 10];
float coordinates [rangeN][1..2]=...;
float lij[rangeN][rangeN];// la matrice est symétrique donc c'est optimisable
execute{
for (var i in rangeN){for(var j in rangeN){
	lij[i][j]=Math.pow(  Math.pow(coordinates[i][1]-coordinates[j][1],2)+Math.pow(coordinates[i][2]-coordinates[j][2],2) ,0.5  );
		}
	}
}
 
 tuple Edge{
 	int i;
 	int j;}

{Edge} edges = {<i,j> | i in rangeN, j in rangeN : i<j};



// variables de décision, x,y,alpha,beta
dvar boolean x[e in edges];
dvar boolean yki[1..K][rangeN];
dvar float+ alpha;
dvar float+ beta[e in edges];
dvar float+ gammak[1..K];
dvar float+ zetakv[1..K][rangeN];

dvar float objectif;


// objectif
minimize objectif;

// contraintes

subject to{
objectif>=L*alpha + sum(e in edges) ( lij[e.i][e.j]*x[e] +3*beta[e] );// sommer juste les j>i


forall(k in 1..K, e in edges){x[e]>=yki[k][e.i]+yki[k][e.j]-1;}
/*
forall(k in 1..K, e in edges){x[e]+yki[k][e.i]-yki[k][e.j]<=1;}
forall(k in 1..K, e in edges){x[e]-yki[k][e.i]+yki[k][e.j]<=1;}
forall(k in 1..K, e in edges){-x[e]+yki[k][e.i]+yki[k][e.j]<=1;}
*/
//forall(i in rangeN,j in rangeN){ alpha+beta[i][j]>=(lh[i]+lh[j])*x[i][j]; }
forall(e in edges){ alpha+beta[e]>=(lh[e.i]+lh[e.j])*x[e]; }

forall(k in 1..K){ W*gammak[k]+ sum(v in rangeN)(w_v[v]*yki[k][v]+W_v[v]*zetakv[k][v]) <= B; }
forall(k in 1..K, v in rangeN){gammak[k]+zetakv[k][v]>=w_v[v]*yki[k][v];} //quel sens pour l'inégalité

forall(i in rangeN){sum(k in 1..K)(yki[k][i])==1;}      //j'ai enlevé le pour tout k


/*alpha>=0;
forall(e in edges){beta[e]>=0;}
forall(k in 1..K){gammak[k]>=0;}
forall(v in rangeN,k in 1..K){zetakv[k][v]>=0;}*/


}

main{
thisOplModel.generate();
if (cplex.solve()) {
var obj=cplex.getObjValue();
thisOplModel.postProcess();
writeln("Integer Model"); 
writeln("OBJECTIVE: ",cplex.getObjValue());
writeln(thisOplModel.yki);
}
}


