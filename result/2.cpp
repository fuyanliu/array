#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <cassert>


using std::cin;
using std::cout;
using std::endl;

int main(int argc, char* argv[]) {
  std::ifstream afile;
  //char unuse[255];
  std::vector<double> data;
  double A;
  //std::string file = "";
  char file[10];
  strcpy(file,argv[1]);
  //std::string filetxt = file ;
  char c[] = {'.','t','x','t','\0'};
  char cc[] = {'.','d','a','t','\0'};
  int b = sizeof(cc);
  std::string filedat = file ;
  std::string filetxt = "";
  filetxt = filetxt + file[0];
  for(int a = 0; a < b; a++){
      filedat = filedat + cc[a];
      filetxt = filetxt + c[a];
    }
  cout << filedat << endl;
  afile.open(filedat);
  while (1) {
      afile >> A;
      if (A == 0) continue;
      if (!afile.good()) break;
      data.push_back(A);
    }
  afile.close();
  int size = data.size();
  std::ofstream bfile(filetxt,std::ios::app);
  bfile<<"0"<<endl;
  for (int i=1; i<100000; i++){
    double edep = 0;
    for (int j=0; j<size; j++){
      if (data[j] == (double)i) {
        int start = j;
        for (int stop = start+1;stop < start+50; stop++){
          if (data[stop]==(double)i){
            for(int num = start+1;num<stop;num++){
              edep = edep + data[num];
              //cout<<i<<"   "<<edep<<endl;
            }
            //cout<<i<<"   "<<edep<<endl;
            //bfile<<edep<<endl;
            continue;
          }
        }
      }
    }
    bfile<<edep<<endl;
  }
    //cout << data[ii][i] << endl;
    //std::ofstream bfile("numall.dat",std::ios::app);
    //bfile << data[ii][i] << endl; 
    
  /*printf("found 1 %d times\n", n1);
  printf("found 2 %d times\n", n2);
  printf("found 3 %d times\n", n3);
  printf("found 4 %d times\n", n4);*/
}
