/*********************************************
 * OPL 12.8.0.0 Model
 * Author: thibaut
 * Creation Date: 25 janv. 2020 at 19:43:06
 *********************************************/
 include "robusteMod.mod";
main{
var duree; var before = new Date(); duree = before.getTime();

//thisOplModel.addDataSource(IloOplDataSource)
thisOplModel.generate();
if (cplex.solve()) {
var obj=cplex.getObjValue();
thisOplModel.postProcess();
writeln("Integer Model"); 
writeln("OBJECTIVE: ",cplex.getObjValue());
writeln(thisOplModel.yki);
}
/*
  var source = new IloOplModelSource("robusteMod.mod");
  var cplex = new IloCplex();
  var def = new IloOplModelDefinition(source);
  var opl = new IloOplModel(def,cplex);
  var data = new IloOplDataSource("10_ulysses_3.dat");
  opl.addDataSource(data);
  opl.generate();
  if (cplex.solve()) {
     writeln("OBJ = " + cplex.getObjValue());
  }  else {
     writeln("No solution");
  }
  var opl2 = new IloOplModel(def,cplex);
  var data2= new IloOplDataElements();
  data2.maxOfx=11;
  opl2.addDataSource("10_ulysses_6.dat");
  opl2.generate();

  if (cplex.solve()) {
     writeln("OBJ = " + cplex.getObjValue());
  } else {
     writeln("No solution");
  }

  opl.end();
  opl2.end();
  data.end(); 
  def.end(); 
  cplex.end(); 
  source.end();
*/



var after = new Date(); writeln("solving time ~= ",after.getTime()-duree);

}
 