/*********************************************
 * OPL 12.8.0.0 Model
 * Author: thibaut
 * Creation Date: 26 janv. 2020 at 09:50:47
 *********************************************/
 
int n=...;
range rangeN=1..n;   
int L=...;
int W=...;
int K=...;
int B=...;
int w_v[rangeN]=...;
float W_v[rangeN] =...;
int lh[rangeN] = ...;
float coordinates [rangeN][1..2]=...;
float lij[rangeN][rangeN];

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

float temp;
execute {
	var before = new Date();
	temp = before.getTime();
}

// variables de décision, x,y,alpha,beta
dvar float+ x[e in edges];
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
objectif>=L*alpha + sum(e in edges) ( lij[e.i][e.j]*x[e] +3*beta[e] );


forall(k in 1..K, e in edges){x[e]>=yki[k][e.i]+yki[k][e.j]-1;}

forall(e in edges){ alpha+beta[e]>=(lh[e.i]+lh[e.j])*x[e]; }

forall(k in 1..K){ W*gammak[k]+ sum(v in rangeN)(w_v[v]*yki[k][v]+W_v[v]*zetakv[k][v]) <= B; }
forall(k in 1..K, v in rangeN){gammak[k]+zetakv[k][v]>=w_v[v]*yki[k][v];}

forall(i in rangeN){sum(k in 1..K)(yki[k][i])==1;}

forall (i in 1..K-1, k in i+1..K) { yki[k][i] == 0; }
forall (i in 2..n, k1 in 2..K-1, k2 in k1+1..K : k2>k1) { yki[k2][i] <= sum(j in rangeN : j<i) yki[k1][j]; }


}


main{
	thisOplModel.generate();
	cplex.tilim = 2*60;
	cplex.solve();
	thisOplModel.postProcess();
}


execute{
	var after = new Date();
	var solvingTime = (after.getTime() - temp)/1000;
	
	var gap = cplex.getMIPRelativeGap();
	var bestInteger = cplex.getObjValue();
	var infBound = cplex.getBestObjValue();
	
	var output = new IloOplOutputFile("output.dat");
	var stat = cplex.status;
	writeln("status = " + stat);
	writeln("solving time = " + solvingTime);
	writeln("best integer solution = " + bestInteger);
	writeln("best inf bound = " + infBound);
	writeln("gap = " + gap);
	output.writeln("status = " + cplex.status + ";");
	output.writeln("solving time = " + solvingTime);
	output.writeln("best integer solution = " + bestInteger);
	output.writeln("best inf bound = " + infBound);
	output.writeln("gap = " + gap);
	output.close();
}
