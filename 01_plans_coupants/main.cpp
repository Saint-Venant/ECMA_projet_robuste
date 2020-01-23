#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
ILOSTLBEGIN
using namespace std;

#include <vector>
#include <cmath>

typedef IloArray<IloBoolVarArray> VarMatrix;

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

void initU1_s(std::vector<std::vector<std::vector<float>>>& U1_s, std::vector<std::vector<float>>& l) {
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

void initU2_s(std::vector<std::vector<float>>& U2_s, IloNumArray& w_v) {
  std::vector<float> mat;
  for (int i=0; i<n; i++) {
    mat.push_back(w_v[i]);
  }
  U2_s.push_back(mat);
}


int main(int argc, char* argv[]){
  string instanceName;
  for (int i = 0; i < argc; i++){
      if (string(argv[i]).compare("-instanceName") == 0)
          instanceName = argv[i + 1];
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
  initU1_s(U1_s, l);
  initU2_s(U2_s, w_v);

  IloModel model(env);

  // Variables
  VarMatrix x(env, n);
  for (int i=0; i<n; i++) {
    x[i] = IloBoolVarArray(env, K);
  }
  VarMatrix y(env, n-1);
  for (int i=0; i<n-1; i++) {
    y[i] = IloBoolVarArray(env, n-i-1);
  }
  IloNumVar z(env, 0);

  // Objective function
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

  // Resolution
  IloCplex cplex(model);
  cplex.solve();

  // Results
  cout << "objective: " << cplex.getObjValue() << endl << endl;
  IloArray<IloNumArray> xSolution(env, n);
  for (int i=0; i<n; i++) {
    IloNumArray sol_x(env, K);
    cplex.getValues(x[i], sol_x);
    xSolution[i] = sol_x;
  }
  for (int i=0; i<n; i++) {
    for (int k=0; k<K; k++) {
      if (xSolution[i][k] == 1) {
        cout << "vertex " << i << " : " << k << endl;
        break;
      }
    }
  }

  IloArray<IloNumArray> ySolution(env, n);
  for (int i=0; i<n-1; i++) {
    IloNumArray sol_y(env, n-i-1);
    cplex.getValues(y[i], sol_y);
    ySolution[i] = sol_y;
  }
  cout << endl << endl;
  for (int i=0; i<n-1; i++) {
    cout << i << "  -> ";
    for (int j=i+1; j<n; j++) {
      cout << ySolution[i][j-i-1] << " ";
    }
    cout << endl;
  }

  env.end();
  return 0;
}
