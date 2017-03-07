#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using std::cin;
using std::cout;
using std::endl;

int main() {
  std::ifstream afile;
  //char unuse[255];
  std::vector<double> data;
  double A;
  std::string file = "detector";
  //std::string filetxt = file ;
  char c[] = {'.','t','x','t','\0'};
  char cc[] = {'.','d','a','t','\0'};
  int b = sizeof(cc);
  std::string filedat = file ;
  std::string filetxt = file ;
  for(int a = 0; a < b; a++){
      filedat = filedat + cc[a];
      filetxt = filetxt + c[a];
    }
  cout << filedat << endl;
  afile.open(filedat);
  while (1) {
      afile >> A;
      if (!afile.good()) break;
      data.push_back(A);
    }
  afile.close();
  int size = data.size();
  std::ofstream bfile(filetxt,std::ios::app);
  //bfile << data[ii][i] << endl;
  for (int i=0; i<100000; i++){
    for (int j=0; j<size; j++){
      if (data[j] == (double)i){
        if(data[j+1] == data[j]) bfile<<"0"<<endl;
        if(data[j+1] != data[j] && data[j+2] == (double)i) bfile<<data[j+1]<<endl;
        if(data[j+1] != data[j] && data[j+2] != (double)i && data[j+3] == (double)i) bfile<<data[j+2]<<endl; 
      }
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
