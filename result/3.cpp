#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>


using std::cin;
using std::cout;
using std::endl;

int main() {
  std::ifstream afile,cfile;
  //char unuse[255];
  std::vector<double> data[2];
  double A,AA;
  std::string file1 = "a.txt";
  std::string file2 = "detector.txt";

  afile.open(file1);
  while (1) {
      afile >> A;
      //if (A == 0) continue;
      if (!afile.good()) break;
      data[0].push_back(A);
    }
  afile.close();

  cfile.open(file2);
  while(1){
    cfile >> AA;
    //double B = AA*1000;
    if (!cfile.good()) break;
    data[1].push_back(AA);
  }
  cfile.close();
  int size = data[0].size();
  std::ofstream bfile("a_compare.txt",std::ios::app);
  cout << size<<"?"<<data[1].size()<<endl;

  for (int i; i<size; i++){
    if (data[0][i] != 0 ) {
      double C = data[0][i] + data[1][i]*1000;
      if (C < 510) bfile<<i<<"  "<<C<<endl;
    }
  }
  
    //cout << data[ii][i] << endl;
    //std::ofstream bfile("numall.dat",std::ios::app);
    //bfile << data[ii][i] << endl; 
    
  /*printf("found 1 %d times\n", n1);
  printf("found 2 %d times\n", n2);
  printf("found 3 %d times\n", n3);
  printf("found 4 %d times\n", n4);*/
}
