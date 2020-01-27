#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
ILOSTLBEGIN
using namespace std;

#include <vector>
#include <cmath>
#include <ctime>

typedef IloArray<IloBoolVarArray> VarBoolMatrix;
typedef IloArray<IloNumVarArray> VarNumMatrix;
typedef IloArray<IloNumArray> NumMatrix;

// --- Data of the model
IloInt n;
IloInt L;
IloInt W;
IloInt K;
IloInt B;

IloNumArray w_v;
IloNumArray W_v;
IloNumArray lh;

std::vector<float> coord_x;
std::vector<float> coord_y;
std::vector<std::vector<float>> l;


void getData(string instanceFileName) {
  std::ifstream instanceFile;
  std::string line;
  std::stringstream strLine;
  std::string entity;

  instanceFile = std::ifstream(instanceFileName.c_str());
  if (!instanceFile.is_open()) throw std::runtime_error("PROBLEM : No instance file found");

  std::getline(instanceFile, line);
  strLine = std::stringstream(line);
  for (int i=0; i<3; i++) {
    std::getline(strLine, entity, ' ');
  }
  n = std::stoi(entity);

  std::getline(instanceFile, line);
  strLine = std::stringstream(line);
  for (int i=0; i<3; i++) {
    std::getline(strLine, entity, ' ');
  }
  L = std::stoi(entity);

  std::getline(instanceFile, line);
  strLine = std::stringstream(line);
  for (int i=0; i<3; i++) {
    std::getline(strLine, entity, ' ');
  }
  W = std::stoi(entity);

  std::getline(instanceFile, line);
  strLine = std::stringstream(line);
  for (int i=0; i<3; i++) {
    std::getline(strLine, entity, ' ');
  }
  K = std::stoi(entity);

  std::getline(instanceFile, line);
  strLine = std::stringstream(line);
  for (int i=0; i<3; i++) {
    std::getline(strLine, entity, ' ');
  }
  B = std::stoi(entity);

  std::getline(instanceFile, line);
  strLine = std::stringstream(line);
  for (int i=0; i<3; i++) {
    std::getline(strLine, entity, ' ');
  }
  entity.erase(0,1);
  entity.pop_back();
  w_v.add(std::stof(entity));
  for (int i=0; i<n-1; i++) {
    std::getline(strLine, entity, ' ');
    entity.pop_back();
    w_v.add(std::stof(entity));
  }

  std::getline(instanceFile, line);
  strLine = std::stringstream(line);
  for (int i=0; i<3; i++) {
    std::getline(strLine, entity, ' ');
  }
  entity.erase(0,1);
  entity.pop_back();
  W_v.add(std::stof(entity));
  for (int i=0; i<n-1; i++) {
    std::getline(strLine, entity, ' ');
    entity.pop_back();
    W_v.add(std::stof(entity));
  }

  std::getline(instanceFile, line);
  strLine = std::stringstream(line);
  for (int i=0; i<3; i++) {
    std::getline(strLine, entity, ' ');
  }
  entity.erase(0,1);
  entity.pop_back();
  lh.add(std::stof(entity));
  for (int i=0; i<n-1; i++) {
    std::getline(strLine, entity, ' ');
    entity.pop_back();
    lh.add(std::stof(entity));
  }

  std::getline(instanceFile, line);
  for (int i=0; i<n; i++) {
    std::getline(instanceFile, line);
    strLine = std::stringstream(line);
    std::getline(strLine, entity, ' ');
    coord_x.push_back(std::stof(entity));
    std::getline(strLine, entity, ' ');
    coord_y.push_back(std::stof(entity));
  }

  cout << "n = " << n << endl;
  cout << "L = " << L << endl;
  cout << "W = " << W << endl;
  cout << "K = " << K << endl;
  cout << "B = " << B << endl;
  cout << "w_v = [";
  for (int i=0; i<n-1; i++) {
    cout << w_v[i] <<";";
  }
  cout << w_v[n-1] << "]" << endl;
  cout << "lh = [";
  for (int i=0; i<n-1; i++) {
    cout << lh[i] <<";";
  }
  cout << lh[n-1] << "]" << endl;
  cout << "coord = [" << endl;
  for (int i=0; i<n; i++) {
    cout << "  " << coord_x[i] << " " << coord_y[i] << endl;
  }
  cout << "]" << endl;
}

void initLength(std::vector<std::vector<float>>& l, std::vector<float>& coord_x, std::vector<float>& coord_y) {
  l.clear();
  float d;
  for (int i=0; i<n-1; i++) {
    std::vector<float> row;
    for (int j=i+1; j<n; j++) {
      d = sqrt(pow(coord_x[i] - coord_x[j], 2) + pow(coord_y[i] - coord_y[j], 2));
      row.push_back(d);
    }
    l.push_back(row);
  }
}

void initU1_s(std::vector<std::vector<std::vector<float>>>& U1_s) {
  std::vector<std::vector<float>> mat;
  for (int i=0; i<n-1; i++) {
    std::vector<float> row;
    for (int j=i+1; j<n; j++) {
      row.push_back(l[i][j]);
    }
    mat.push_back(row);
  }
  U1_s.push_back(mat);
}

void initU2_s(std::vector<std::vector<float>>& U2_s) {
  std::vector<float> mat;
  for (int i=0; i<n; i++) {
    mat.push_back(w_v[i]);
  }
  U2_s.push_back(mat);
}

void setModel(IloEnv& env, IloModel& model, VarBoolMatrix& x, VarBoolMatrix& y, IloNumVar& z,
              std::vector<std::vector<std::vector<float>>>& U1_s, std::vector<std::vector<float>>& U2_s) {
  // variables
  for (int i=0; i<n; i++) {
    x[i] = IloBoolVarArray(env, K);
  }
  for (int i=0; i<n-1; i++) {
    y[i] = IloBoolVarArray(env, n-i-1);
  }

  // objective function
  IloObjective obj(env, z, IloObjective::Minimize);
  model.add(obj);

  // Constraints
  std::vector<std::vector<float>> l1;
  for (unsigned int p=0; p<U1_s.size(); p++) {
    IloExpr exprCt1(env);
    l1 = U1_s[p];
    for (int i=0; i<n-1; i++) {
      for (int j=i+1; j<n; j++) {
        exprCt1 += l1[i][j-i-1]*y[i][j-i-1];
      }
    }
    model.add(z >= exprCt1);
    exprCt1.end();
  }

  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
      for (int k=0; k<K; k++) {
        model.add(y[i][j-i-1] >= x[i][k] + x[j][k] - 1);
      }
    }
  }

  std::vector<float> w2;
  for (unsigned int p=0; p<U2_s.size(); p++) {
    w2 = U2_s[p];
    for (int k=0; k<K; k++) {
      IloExpr exprCt3(env);
      for (int i=0; i<n; i++) {
        exprCt3 += w2[i]*x[i][k];
      }
      model.add(exprCt3 <= B);
      exprCt3.end();
    }
  }

  for (int i=0; i<n; i++) {
    IloExpr exprCt4(env);
    for (int k=0; k<K; k++) {
      exprCt4 += x[i][k];
    }
    model.add(exprCt4 == 1);
    exprCt4.end();
  }
}

void getSolution(IloNum& valueSolution, IloNum& bestInfBound, NumMatrix& xSolution, NumMatrix& ySolution,
                 VarBoolMatrix& x, VarBoolMatrix& y, IloCplex& cplex, IloEnv& env) {
  valueSolution = cplex.getObjValue();
  bestInfBound = cplex.getBestObjValue();
  for (int i=0; i<n; i++) {
    IloNumArray sol_x(env, K);
    cplex.getValues(x[i], sol_x);
    xSolution[i] = sol_x;
  }
  for (int i=0; i<n-1; i++) {
    IloNumArray sol_y(env, n-i-1);
    cplex.getValues(y[i], sol_y);
    ySolution[i] = sol_y;
  }
}

void displaySolution(IloNum& valueSolution, IloNum& bestInfBound, NumMatrix& xSolution,
                     NumMatrix& ySolution) {
  cout << "objective Master: " << valueSolution << endl;
  cout << "best inf bound : " << bestInfBound << endl << endl;

  cout << "Vertex > cluster" << endl;
  for (int i=0; i<n; i++) {
    for (int k=0; k<K; k++) {
      if (xSolution[i][k] == 1) {
        cout << "vertex " << i << " : " << k << endl;
        break;
      }
    }
  }

  cout << endl << "Variables y_ij" << endl;
  for (int i=0; i<n-1; i++) {
    cout << i << "  -> ";
    for (int j=i+1; j<n; j++) {
      cout << ySolution[i][j-i-1] << " ";
    }
    cout << endl;
  }
}

void setModel1(IloEnv& env1, IloModel& model1, VarNumMatrix& delta1, NumMatrix& ySolution) {
  // variable
  for (int i=0; i<n-1; i++) {
    delta1[i] = IloNumVarArray(env1, n-i-1, 0, 3);
  }

  // objective function
  IloExpr exprObj(env1);
  for (int i=0; i<n-1; i++){
    for (int j=i+1; j<n; j++) {
      exprObj += ySolution[i][j-i-1]*(l[i][j-i-1] + delta1[i][j-i-1]*(lh[i] + lh[j]));
    }
  }
  IloObjective obj(env1, exprObj, IloObjective::Maximize);
  model1.add(obj);
  exprObj.end();

  // constraints
  IloExpr exprCtsum(env1);
  for (int i=0; i<n-1; i++){
    for (int j=i+1; j<n; j++) {
      exprCtsum += delta1[i][j-i-1];
    }
  }
  model1.add(exprCtsum <= L);
  exprCtsum.end();
}

void getSolution1(IloNum& value1Solution, NumMatrix& delta1Solution, VarNumMatrix& delta1, IloCplex& cplex) {
  value1Solution = cplex.getObjValue();
  for (int i=0; i<n-1; i++) {
    cplex.getValues(delta1[i], delta1Solution[i]);
  }
}

void setModel2(int k, IloEnv& env2, IloModel& model2, NumMatrix& xSolution, IloNumVarArray delta2) {
  // objective function
  IloExpr exprObj(env2);
  for (int i=0; i<n; i++){
    exprObj += xSolution[i][k]*w_v[i]*(1 + delta2[i]);
  }
  IloObjective obj(env2, exprObj, IloObjective::Maximize);
  model2.add(obj);
  exprObj.end();

  // constraints
  for (int i=0; i<n; i++) {
    model2.add(delta2[i] <= W_v[i]);
  }

  IloExpr exprCtsum(env2);
  for (int i=0; i<n; i++){
    exprCtsum += delta2[i];
  }
  model2.add(exprCtsum <= W);
  exprCtsum.end();
}

void getSolution2(IloNum& value2Solution, IloNumArray& delta2Solution, IloNumVarArray& delta2, IloCplex& cplex) {
  value2Solution = cplex.getObjValue();
  cplex.getValues(delta2, delta2Solution);
}

void addCutSP1(IloEnv& env, IloModel& model, VarBoolMatrix& y, IloNumVar& z, NumMatrix& delta1Solution,
               std::vector<std::vector<std::vector<float>>>& U1_s) {
  IloExpr exprCtsum(env);
  float coef;
  std::vector<std::vector<float>> l1;
  for (int i=0; i<n-1; i++) {
    std::vector<float> row;
    for (int j=i+1; j<n; j++) {
      coef = l[i][j-i-1] + delta1Solution[i][j-i-1]*(lh[i] + lh[j]);
      row.push_back(coef);
      exprCtsum += coef*y[i][j-i-1];
    }
    l1.push_back(row);
  }
  U1_s.push_back(l1);
  model.add(z >= exprCtsum);
  exprCtsum.end();
} // to remove

void addCutSP2(IloEnv& env, IloModel& model, VarBoolMatrix& x, IloNumArray& delta2Solution,
               std::vector<std::vector<float>>& U2_s) {
  std::vector<float> w2;
  for (int i=0; i<n; i++) {
    w2.push_back(w_v[i]*(1 + delta2Solution[i]));
  }
  U2_s.push_back(w2);
  for (int k=0; k<K; k++) {
    IloExpr exprCtsum(env);
    for (int i=0; i<n; i++) {
      exprCtsum += w2[i]*x[i][k];
    }
    model.add(exprCtsum <= B);
    exprCtsum.end();
  }
}

void saveResults(IloNum valueSolution, IloNum bestInfBound, double calculationTime, string instanceName,
                 string outputFileName) {
  std::ofstream outputFile;
  outputFile.open(outputFileName.c_str());
  outputFile << "instance: " << instanceName << endl;
  outputFile << "calculation_time: " << calculationTime << endl;
  outputFile << "best integer : " << valueSolution << endl;
  outputFile << "best inf bound : " << bestInfBound << endl;
  outputFile.close();
}

struct Args {
  IloEnv env;
  IloModel model;
  VarBoolMatrix x;
  VarBoolMatrix y;
  IloNumVar z;
  std::vector<std::vector<std::vector<float>>> U1_s;
  std::vector<std::vector<float>> U2_s;

  Args(IloEnv& env, IloModel& model, VarBoolMatrix& x, VarBoolMatrix& y, IloNumVar& z,
       std::vector<std::vector<std::vector<float>>>& U1_s, std::vector<std::vector<float>>& U2_s) :
       env(env), model(model), x(x), y(y), z(z), U1_s(U1_s), U2_s(U2_s) {}
};

/*
ILOLAZYCONSTRAINTCALLBACK6(SP1, IloCplex&, cplex, IloEnv&, env, IloModel&, model, VarBoolMatrix&, y, IloNumVar&, z,
                           std::vector<std::vector<std::vector<float>>>&, U1_s) {
  IloEnv env1;
  IloModel model1(env1);

  NumMatrix ySolution(env1, n);
  for (int i=0; i<n-1; i++) {
    IloNumArray sol_y(env1, n-i-1);
    cplex.getValues(y[i], sol_y);
    ySolution[i] = sol_y;
  }
  IloNum valueSolution = cplex.getObjValue();

  IloNum value1Solution;
  NumMatrix delta1Solution(env1, n-1);
  for (int i=0; i<n-1; i++) {
    delta1Solution[i] = IloNumArray(env1, n-i-1);
  }

  // Variables
  VarNumMatrix delta1(env1, n-1);

  setModel1(env1, model1, delta1, ySolution);

  // Resolution
  IloCplex cplex1(model1);
  cplex1.solve();

  // Results
  getSolution1(value1Solution, delta1Solution, delta1, cplex1);

  // Add lazy constraint if necessary
  if (valueSolution < value1Solution) {
    // add cut
    addCutSP1(env, model, y, z, delta1Solution, U1_s);
  }

  env1.end();
}*/

/*
ILOLAZYCONSTRAINTCALLBACK5(SP2, IloCplex&, cplex, IloEnv&, env, IloModel&, model, VarBoolMatrix&, x,
                           std::vector<std::vector<float>>&, U2_s) {
  for (int k=0; k<K; k++) {
    IloEnv env2;
    IloModel model2(env2);

    NumMatrix xSolution(env2, n);
    for (int i=0; i<n; i++) {
      IloNumArray sol_x(env2, K);
      cplex.getValues(x[i], sol_x);
      xSolution[i] = sol_x;
    }

    IloNum value2Solution;
    IloNumArray delta2Solution(env2, n);

    // Variables
    IloNumVarArray delta2(env2, n, 0, W);

    setModel2(k, env2, model2, xSolution, delta2);

    // Resolution
    IloCplex cplex2(model2);
    cplex2.solve();

    // Results
    getSolution2(value2Solution, delta2Solution, delta2, cplex2);
    if (value2Solution > B) {
      // add cuts
      addCutSP2(env, model, x, delta2Solution, U2_s);
      break;
    }

    env2.end();
  }
}*/

ILOLAZYCONSTRAINTCALLBACK2(SP1, IloCplex, cplex, Args, args) {
  IloEnv env1;
  IloModel model1(env1);

  NumMatrix ySolution(env1, n);
  for (int i=0; i<n-1; i++) {
    IloNumArray sol_y(env1, n-i-1);
    cplex.getValues(args.y[i], sol_y);
    ySolution[i] = sol_y;
  }
  IloNum valueSolution = cplex.getObjValue();

  IloNum value1Solution;
  NumMatrix delta1Solution(env1, n-1);
  for (int i=0; i<n-1; i++) {
    delta1Solution[i] = IloNumArray(env1, n-i-1);
  }

  // Variables
  VarNumMatrix delta1(env1, n-1);

  setModel1(env1, model1, delta1, ySolution);

  // Resolution
  IloCplex cplex1(model1);
  cplex1.solve();

  // Results
  getSolution1(value1Solution, delta1Solution, delta1, cplex1);

  // Add lazy constraint if necessary
  if (valueSolution < value1Solution) {
    // add cut
    addCutSP1(args.env, args.model, args.y, args.z, delta1Solution, args.U1_s);
  }

  env1.end();
}

ILOLAZYCONSTRAINTCALLBACK2(SP2, IloCplex, cplex, Args, args) {
  for (int k=0; k<K; k++) {
    IloEnv env2;
    IloModel model2(env2);

    NumMatrix xSolution(env2, n);
    for (int i=0; i<n; i++) {
      IloNumArray sol_x(env2, K);
      cplex.getValues(args.x[i], sol_x);
      xSolution[i] = sol_x;
    }

    IloNum value2Solution;
    IloNumArray delta2Solution(env2, n);

    // Variables
    IloNumVarArray delta2(env2, n, 0, W);

    setModel2(k, env2, model2, xSolution, delta2);

    // Resolution
    IloCplex cplex2(model2);
    cplex2.solve();

    // Results
    getSolution2(value2Solution, delta2Solution, delta2, cplex2);
    if (value2Solution > B) {
      // add cuts
      addCutSP2(args.env, args.model, args.x, delta2Solution, args.U2_s);
      break;
    }

    env2.end();
  }
}



int main(int argc, char* argv[]){
  string instanceName;
  time_t timeBegin = time(NULL);
  double maxTime = 60;
  string outputFileName = "defaultSave.txt";
  for (int i = 0; i < argc; i++){
      if (string(argv[i]).compare("-instanceName") == 0)
          instanceName = argv[i + 1];
      if (string(argv[i]).compare("-maxTime") == 0)
          maxTime = std::stoi(argv[i + 1]);
      if (string(argv[i]).compare("-outputFileName") == 0)
          outputFileName = argv[i + 1];
  }

  IloEnv env;
  w_v=IloNumArray(env);
  W_v=IloNumArray(env);
  lh=IloNumArray(env);

  // read data
  getData(instanceName);
  initLength(l, coord_x, coord_y);

  // sets of current cuts
  std::vector<std::vector<std::vector<float>>> U1_s;
  std::vector<std::vector<float>> U2_s;
  initU1_s(U1_s);
  initU2_s(U2_s);

  // sets to track evolutions

  IloModel model(env);

  // Variables
  VarBoolMatrix x(env, n);
  VarBoolMatrix y(env, n-1);
  IloNumVar z(env, 0);

  setModel(env, model, x, y, z, U1_s, U2_s);

  // Resolution
  IloNum valueSolution;
  IloNum bestInfBound;
  NumMatrix xSolution(env, n);
  NumMatrix ySolution(env, n);

  // Master problem
  IloCplex cplex(model);
  time_t timer = time(NULL);
  double dt = maxTime - difftime(timer, timeBegin);
  try {
    cplex.setParam(IloCplex::TiLim, dt);
    cplex.use(SP1(env, cplex, Args(env, model, x, y, z, U1_s, U2_s)));
    cplex.use(SP2(env, cplex, Args(env, model, x, y, z, U1_s, U2_s)));
    cplex.solve();
    getSolution(valueSolution, bestInfBound, xSolution, ySolution, x, y, cplex, env);
    displaySolution(valueSolution, bestInfBound, xSolution, ySolution);
  } catch (const IloException& e){
    cerr << e;
    throw;
  }


  env.end();
  timer = time(NULL);
  saveResults(valueSolution, bestInfBound, difftime(timer, timeBegin), instanceName, outputFileName);
  cout << ">> " << difftime(timer, timeBegin) << " seconds" << endl;
  cout << ">> objective value : " << valueSolution << endl;
  cout << ">> best inf bound : " << bestInfBound << endl;
  return 0;
}
