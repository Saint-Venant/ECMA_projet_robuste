#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
ILOSTLBEGIN
using namespace std;

#include <vector>
#include <cmath>
#include <ctime>
#include <cassert>

typedef IloArray<IloBoolVarArray> VarBoolMatrix;
typedef IloArray<IloNumVarArray> VarNumMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloBoolArray> BoolMatrix;

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

  // break symetry
  assert(K <= n);
  for (int i=0; i<K-1; i++) {
    for (int k=i+1; k<K; k++) {
      model.add(x[i][k] == 0);
    }
  }
  for (int i=1; i<n; i++) {
    for (int k=1; k<K-1; k++) {
      for (int kk=k+1; kk<K; kk++) {
        IloExpr exprCtSymetry(env);
        for (int j=0; j<i; j++) {
          exprCtSymetry += x[j][k];
        }
        model.add(x[i][kk] <= exprCtSymetry);
        exprCtSymetry.end();
      }
    }
  }
}

void setLpModel(IloEnv& lpEnv, IloModel& lpModel, VarNumMatrix& lpx, VarNumMatrix& lpy, IloNumVar& lpz,
                std::vector<std::vector<std::vector<float>>>& U1_s, std::vector<std::vector<float>>& U2_s) {
  // variables
  for (int i=0; i<n; i++) {
    lpx[i] = IloNumVarArray(lpEnv, K, 0, 1);
  }
  for (int i=0; i<n-1; i++) {
    lpy[i] = IloNumVarArray(lpEnv, n-i-1, 0, 1);
  }

  // objective function
  IloObjective obj(lpEnv, lpz, IloObjective::Minimize);
  lpModel.add(obj);

  // Constraints
  std::vector<std::vector<float>> l1;
  for (unsigned int p=0; p<U1_s.size(); p++) {
    IloExpr exprCt1(lpEnv);
    l1 = U1_s[p];
    for (int i=0; i<n-1; i++) {
      for (int j=i+1; j<n; j++) {
        exprCt1 += l1[i][j-i-1]*lpy[i][j-i-1];
      }
    }
    lpModel.add(lpz >= exprCt1);
    exprCt1.end();
  }

  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
      for (int k=0; k<K; k++) {
        lpModel.add(lpy[i][j-i-1] >= lpx[i][k] + lpx[j][k] - 1);
      }
    }
  }

  std::vector<float> w2;
  for (unsigned int p=0; p<U2_s.size(); p++) {
    w2 = U2_s[p];
    for (int k=0; k<K; k++) {
      IloExpr exprCt3(lpEnv);
      for (int i=0; i<n; i++) {
        exprCt3 += w2[i]*lpx[i][k];
      }
      lpModel.add(exprCt3 <= B);
      exprCt3.end();
    }
  }

  for (int i=0; i<n; i++) {
    IloExpr exprCt4(lpEnv);
    for (int k=0; k<K; k++) {
      exprCt4 += lpx[i][k];
    }
    lpModel.add(exprCt4 == 1);
    exprCt4.end();
  }

  // break symetry
  assert(K <= n);
  for (int i=0; i<K-1; i++) {
    for (int k=i+1; k<K; k++) {
      lpModel.add(lpx[i][k] == 0);
    }
  }
  for (int i=1; i<n; i++) {
    for (int k=1; k<K-1; k++) {
      for (int kk=k+1; kk<K; kk++) {
        IloExpr exprCtSymetry(lpEnv);
        for (int j=0; j<i; j++) {
          exprCtSymetry += lpx[j][k];
        }
        lpModel.add(lpx[i][kk] <= exprCtSymetry);
        exprCtSymetry.end();
      }
    }
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

void computeCutSP1(std::vector<std::vector<float>>& l1, VarBoolMatrix& y, NumMatrix& delta1Solution,
                   std::vector<std::vector<std::vector<float>>>& U1_s) {
  float coef;
  for (int i=0; i<n-1; i++) {
    std::vector<float> row;
    for (int j=i+1; j<n; j++) {
      coef = l[i][j-i-1] + delta1Solution[i][j-i-1]*(lh[i] + lh[j]);
      row.push_back(coef);
    }
    l1.push_back(row);
  }
  U1_s.push_back(l1);
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

bool improveModel_SP1(IloCplex& lpCplex, std::vector<std::vector<float>>& l1, VarBoolMatrix& y, VarNumMatrix& lpy,
                      std::vector<std::vector<std::vector<float>>>& U1_s) {
  lpCplex.solve();
  IloNum lpValueSolution = lpCplex.getObjValue();
  bool cutFound1 = false;
  float eps = 0.1;

  // --- Solve SP1
  IloEnv env1;
  IloModel model1(env1);

  NumMatrix lpySolution(env1, n);
  for (int i=0; i<n-1; i++) {
    IloNumArray sol_y(env1, n-i-1);
    lpCplex.getValues(sol_y, lpy[i]);
    lpySolution[i] = sol_y;
  }
  IloNum value1Solution;
  NumMatrix delta1Solution(env1, n-1);
  for (int i=0; i<n-1; i++) {
    delta1Solution[i] = IloNumArray(env1, n-i-1);
  }
  VarNumMatrix delta1(env1, n-1);

  setModel1(env1, model1, delta1, lpySolution);

  // Resolution
  IloCplex cplex1(model1);
  cplex1.setOut(env1.getNullStream());
  cplex1.solve();

  // Results
  getSolution1(value1Solution, delta1Solution, delta1, cplex1);

  // Add cut if necessary
  if (lpValueSolution < value1Solution - eps) {
    cutFound1 = true;
    float coef;
    for (int i=0; i<n-1; i++) {
      std::vector<float> row;
      for (int j=i+1; j<n; j++) {
        coef = l[i][j-i-1] + delta1Solution[i][j-i-1]*(lh[i] + lh[j]);
        row.push_back(coef);
      }
      l1.push_back(row);
    }
    U1_s.push_back(l1);
  }

  env1.end();
  return cutFound1;
}

bool improveModel_SP2(IloCplex& lpCplex, std::vector<float>& w2, VarBoolMatrix& x, VarNumMatrix& lpx,
                      std::vector<std::vector<float>>& U2_s) {
  lpCplex.solve();
  bool cutFound2 = false;
  float eps = 0.5;

  // --- Solve SP2-k
  for (int k=0; k<K; k++) {
    IloEnv env2;
    IloModel model2(env2);

    NumMatrix xSolution(env2, n);
    for (int i=0; i<n; i++) {
      IloNumArray sol_x(env2, K);
      lpCplex.getValues(sol_x, lpx[i]);
      xSolution[i] = sol_x;
    }
    IloNum value2Solution;
    IloNumArray delta2Solution(env2, n);
    IloNumVarArray delta2(env2, n, 0, W);

    setModel2(k, env2, model2, xSolution, delta2);

    // Resolution
    IloCplex cplex2(model2);
    cplex2.setOut(env2.getNullStream());
    cplex2.solve();

    // Results
    getSolution2(value2Solution, delta2Solution, delta2, cplex2);

    // Add lazy constraint if necessary
    if (value2Solution > B + eps) {
      cutFound2 = true;
      for (int i=0; i<n; i++) {
        w2.push_back(w_v[i]*(1 + delta2Solution[i]));
      }
      U2_s.push_back(w2);
      break;
    }

    env2.end();
  }
  return cutFound2;
}

struct DataInfo {
  int nbCutsSP1 = 0;
  int nbCutsSP2 = 0;
  std::vector<float> valuesBestInteger;
  std::vector<float> valuesInfBound;

  DataInfo() : nbCutsSP1(0), nbCutsSP2(0) {}
};

struct LPrelaxation {
  IloEnv lpEnv;
  IloModel lpModel;
  IloCplex lpCplex;
  VarNumMatrix lpx;
  VarNumMatrix lpy;
  IloNumVar lpz;

  LPrelaxation(IloEnv& env, IloModel& model, IloCplex& cplex, VarNumMatrix& x, VarNumMatrix& y, IloNumVar& z) : lpEnv(env), lpModel(model),
               lpCplex(cplex), lpx(x), lpy(y), lpz(z) {}
};

ILOLAZYCONSTRAINTCALLBACK7(myLazyConstraintCallback, VarBoolMatrix, x, VarBoolMatrix, y, IloNumVar, z,
                           std::vector<std::vector<std::vector<float>>>&, U1_s,
                           std::vector<std::vector<float>>&, U2_s, DataInfo&, info, LPrelaxation&, relaxation) {
  cout << endl << "BEGIN LAZY CONSTRAINT CALLBACK" << endl;
  IloEnv masterEnv = getEnv();
  IloNum valueSolution = getObjValue();
  IloNum bestInfBound = getBestObjValue();
  info.valuesBestInteger.push_back(valueSolution);
  info.valuesInfBound.push_back(bestInfBound);
  bool integerCut = false;

  // --- Solve SP1
  IloEnv env1;
  IloModel model1(env1);

  NumMatrix ySolution(env1, n);
  for (int i=0; i<n-1; i++) {
    IloNumArray sol_y(env1, n-i-1);
    getValues(sol_y, y[i]);
    ySolution[i] = sol_y;
  }
  IloNum value1Solution;
  NumMatrix delta1Solution(env1, n-1);
  for (int i=0; i<n-1; i++) {
    delta1Solution[i] = IloNumArray(env1, n-i-1);
  }
  VarNumMatrix delta1(env1, n-1);

  setModel1(env1, model1, delta1, ySolution);

  // Resolution
  IloCplex cplex1(model1);
  cplex1.setOut(env1.getNullStream());
  cplex1.solve();

  // Results
  getSolution1(value1Solution, delta1Solution, delta1, cplex1);

  // Add lazy constraint if necessary
  if (valueSolution < value1Solution) {
    integerCut = true;
    IloExpr exprCtsum1(masterEnv);
    IloExpr lpExprCtsum1(relaxation.lpEnv);
    std::vector<std::vector<float>> l1;
    computeCutSP1(l1, y, delta1Solution, U1_s);
    for (int i=0; i<n-1; i++) {
      for (int j=i+1; j<n; j++) {
        exprCtsum1 += l1[i][j-i-1]*y[i][j-i-1];
        lpExprCtsum1 += l1[i][j-i-1]*relaxation.lpy[i][j-i-1];
      }
    }
    add(z >= exprCtsum1);
    relaxation.lpModel.add(relaxation.lpz >= lpExprCtsum1);
    exprCtsum1.end();
    lpExprCtsum1.end();
    info.nbCutsSP1 += 1;
    cout << "Added cut SP1 !" << endl;
  }

  env1.end();

  // --- Solve SP2-k
  for (int k=0; k<K; k++) {
    IloEnv env2;
    IloModel model2(env2);

    NumMatrix xSolution(env2, n);
    for (int i=0; i<n; i++) {
      IloNumArray sol_x(env2, K);
      getValues(sol_x, x[i]);
      xSolution[i] = sol_x;
    }
    IloNum value2Solution;
    IloNumArray delta2Solution(env2, n);
    IloNumVarArray delta2(env2, n, 0, W);

    setModel2(k, env2, model2, xSolution, delta2);

    // Resolution
    IloCplex cplex2(model2);
    cplex2.setOut(env2.getNullStream());
    cplex2.solve();

    // Results
    getSolution2(value2Solution, delta2Solution, delta2, cplex2);

    // Add lazy constraint if necessary
    if (value2Solution > B) {
      integerCut = true;
      std::vector<float> w2;
      for (int i=0; i<n; i++) {
        w2.push_back(w_v[i]*(1 + delta2Solution[i]));
      }
      U2_s.push_back(w2);
      for (int kk=0; kk<K; kk++) {
        IloExpr exprCtsum2(masterEnv);
        IloExpr lpExprCtsum2(relaxation.lpEnv);
        for (int i=0; i<n; i++) {
          exprCtsum2 += w2[i]*x[i][kk];
          lpExprCtsum2 += w2[i]*relaxation.lpx[i][kk];
        }
        add(exprCtsum2 <= B);
        relaxation.lpModel.add(lpExprCtsum2 <= B);
        exprCtsum2.end();
        lpExprCtsum2.end();
      }
      info.nbCutsSP2 += 1;
      cout << "Added cut SP2 !" << endl;
      break;
    }

    env2.end();
  }

  // --- Attempt to improve the LP relaxation
  bool cutFound1 = true;
  bool cutFound2 = true;
  int improveCuts= 0;
  while ((integerCut) && (cutFound1 || cutFound2) && (improveCuts < 5)) {
    improveCuts += 1;
    std::vector<std::vector<float>> l1;
    cutFound1 = improveModel_SP1(relaxation.lpCplex, l1, y, relaxation.lpy, U1_s);
    if (cutFound1) {
      IloExpr masterExpr1(masterEnv);
      IloExpr lpExpr1(relaxation.lpEnv);
      for (int i=0; i<n-1; i++) {
        for (int j=i+1; j<n; j++) {
          masterExpr1 += l1[i][j-i-1]*y[i][j-i-1];
          lpExpr1 += l1[i][j-i-1]*relaxation.lpy[i][j-i-1];
        }
      }
      add(z >= masterExpr1);
      relaxation.lpModel.add(relaxation.lpz >= lpExpr1);
      masterExpr1.end();
      lpExpr1.end();
      cout << "Added lp cut 1 !" << endl;
    }

    std::vector<float> w2;
    cutFound2 = improveModel_SP2(relaxation.lpCplex, w2, x, relaxation.lpx, U2_s);
    if (cutFound2) {
      for (int k=0; k<K; k++) {
        IloExpr masterExpr2(masterEnv);
        IloExpr lpExpr2(relaxation.lpEnv);
        for (int i=0; i<n; i++) {
          masterExpr2 += w2[i]*x[i][k];
          lpExpr2 += w2[i]*relaxation.lpx[i][k];
        }
        add(masterExpr2 <= B);
        relaxation.lpModel.add(lpExpr2 <= B);
        masterExpr2.end();
        lpExpr2.end();
      }
      cout << "Added lp cut 2 !" << endl;
    }
  }

  cout << "END LAZY CONSTRAINT CALLBACK" << endl << endl;
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

  IloEnv lpEnv;

  // read data
  getData(instanceName);
  initLength(l, coord_x, coord_y);

  // sets of current cuts
  std::vector<std::vector<std::vector<float>>> U1_s;
  std::vector<std::vector<float>> U2_s;
  initU1_s(U1_s);
  initU2_s(U2_s);

  // sets to track evolutions
  DataInfo info = DataInfo();

  IloModel model(env);
  IloModel lpModel(lpEnv);

  // Variables
  VarBoolMatrix x(env, n);
  VarBoolMatrix y(env, n-1);
  IloNumVar z(env, 0);

  VarNumMatrix lpx(lpEnv, n);
  VarNumMatrix lpy(lpEnv, n-1);
  IloNumVar lpz(lpEnv, 0);

  setModel(env, model, x, y, z, U1_s, U2_s);
  setLpModel(lpEnv, lpModel, lpx, lpy, lpz, U1_s, U2_s);

  // Resolution
  IloNum valueSolution;
  IloNum bestInfBound;
  NumMatrix xSolution(env, n);
  NumMatrix ySolution(env, n);

  IloCplex lpCplex(lpModel);
  lpCplex.setOut(lpEnv.getNullStream());
  LPrelaxation relaxation = LPrelaxation(lpEnv, lpModel, lpCplex, lpx, lpy, lpz);

  cout << "run improve model SP1" << endl;
  bool cutFound1 = true;
  bool cutFound2 = true;
  int initialCuts = 0;
  while ((cutFound1 || cutFound2) && (initialCuts < 20)) {
    initialCuts += 1;
    std::vector<std::vector<float>> l1;
    cutFound1 = improveModel_SP1(lpCplex, l1, y, lpy, U1_s);
    if (cutFound1) {
      IloExpr masterExpr1(env);
      IloExpr lpExpr1(lpEnv);
      for (int i=0; i<n-1; i++) {
        for (int j=i+1; j<n; j++) {
          masterExpr1 += l1[i][j-i-1]*y[i][j-i-1];
          lpExpr1 += l1[i][j-i-1]*lpy[i][j-i-1];
        }
      }
      model.add(z >= masterExpr1);
      lpModel.add(lpz >= lpExpr1);
      masterExpr1.end();
      lpExpr1.end();
      cout << "Added lp cut 1 !" << endl;
    }

    std::vector<float> w2;
    cutFound2 = improveModel_SP2(lpCplex, w2, x, lpx, U2_s);
    if (cutFound2) {
      for (int k=0; k<K; k++) {
        IloExpr masterExpr2(env);
        IloExpr lpExpr2(lpEnv);
        for (int i=0; i<n; i++) {
          masterExpr2 += w2[i]*x[i][k];
          lpExpr2 += w2[i]*lpx[i][k];
        }
        model.add(masterExpr2 <= B);
        lpModel.add(lpExpr2 <= B);
        masterExpr2.end();
        lpExpr2.end();
      }
      cout << "Added lp cut 2 !" << endl;
    }
  }

  // Master problem
  IloCplex cplex(model);
  time_t timer = time(NULL);
  double dt = maxTime - difftime(timer, timeBegin);
  try {
    cplex.setParam(IloCplex::TiLim, dt);
    cplex.setParam(IloCplex::Param::Strategy::NodeSelect, 0);
    cplex.use(myLazyConstraintCallback(env, x, y, z, U1_s, U2_s, info, relaxation));
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
  cout << ">> nb cuts SP1 : " << info.nbCutsSP1 << endl;
  cout << ">> nb cuts SP2 : " << info.nbCutsSP2 << endl;
  cout << ">> values best integer : ";
  for (unsigned int p=0; p<info.valuesBestInteger.size(); p++) {
    cout << info.valuesBestInteger[p] << " ";
  }
  cout << endl << ">> values inf bound : ";
  for (unsigned int p=0; p<info.valuesInfBound.size(); p++) {
    cout << info.valuesInfBound[p] << " ";
  }
  cout << endl;
  return 0;
}
