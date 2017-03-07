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
  std::vector<double> data[4];
  double A;
  std::string file = "";
  //std::string filetxt = file ;
  
  char c[] = {'a', 'b','c','d','\0'};
  char cc[] = {'.','d','a','t','\0'};
  int b = sizeof(cc);
  //cout << b << endl;
  for(int ii = 0; ii<4; ii++){
    //filetxt = filetxt + c[ii];
    std::string filedat = file ;
    filedat = filedat + c[ii];
    for(int a = 0; a < b; a++){
      filedat = filedat + cc[a];
    }
    cout << filedat << endl;
    afile.open(filedat);
    // file.getline(unuse, 256, '\n');
    // file.getline(unuse, 256, '\n');
    // file.getline(unuse, 256, '\n');
    while (1) {
      afile >> A;
      if (!afile.good()) break;
      data[ii].push_back(A);
    }
    afile.close();
    //cout << data[ii].size() << endl;
    for (int i = 0; i < data[ii].size(); i++) {
    //cout << data[ii][i] << endl;
    //std::ofstream bfile("numall.dat",std::ios::app);
    //bfile << data[ii][i] << endl; 
    }
  }
    
  //cout << data[1].size() << endl;
  int n1=0,n2=0,n3=0,n4=0;
  for (int i = 0; i < data[1].size(); i++) {
   double num = data[0][i] + data[1][i] + data[2][i]+ data[3][i];
   if(num == 1) n1++;
   if(num == 2) n2++;
   if(num == 3) n3++;
   if(num == 4) n4++;
   //std::ofstream bfile("num_135.dat",std::ios::app);
   //bfile << num << endl; 
  }
  printf("found 1 %d times\n", n1);
  printf("found 2 %d times\n", n2);
  printf("found 3 %d times\n", n3);
  printf("found 4 %d times\n", n4);
}
