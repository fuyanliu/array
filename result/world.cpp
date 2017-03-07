#include <cassert>
#include <cmath>
#include "TCanvas.h"


Float_t a;
Float_t b;
void world_tree() {
  ifstream in;
  std::string file = "d";
  //char file[10];
  //strcpy(file,argv[1]);
  std::string filetxt = file ;
  std::string filedat = file ;
  char c[] = {'.', 't', 'x', 't','\0'};
  char cc[] = {'.','d','a','t','\0'};
  for(int ii = 0; ii<4; ii++){
  	filetxt = filetxt + c[ii];
    filedat = filedat + cc[ii];
  }
  cout<< filetxt << endl;
  in.open(filetxt);
  char unuse[255];
  //in.getline(unuse, 256, '\n');
  //in.getline(unuse, 256, '\n');
  //in.getline(unuse, 256, '\n');
  Float_t x;
  Float_t px;
  Int_t nlines = 0;
  Int_t lines = 0;
  a = 0;
  b = 2;
  TFile *f = new TFile("world.root","RECREATE");
  TTree *mytree = new TTree("tree","distrebution");
  mytree->Branch("px",&px,"px/F");
  while (1) {
    in >> x; //从文件中读取一行,分别赋值给x,y,z。
    //if (nlines < 5) printf("x=%8f\n",x);
    px = x;
    a=max(a,px);
    b=min(b,px);
    if (!in.good()) break;
    //mytree->Fill();//填充TNtuple
    Float_t UL = 0;
    if (px ==0){
      px = 0;
    }
    else{
      mytree->Fill();
      px = 1;
      lines++;
    }
    std::ofstream bfile(filedat,std::ios::app);
    bfile<<px<<endl;
    //h1->Fill(px);
    nlines++;
    
  }
  printf("the max number is %f\n",a);
  printf("the min number is %f\n",b );
  printf(" found %d points\n",nlines);
  printf(" Found %d signals\n",lines);

  //TH1F *h1 = new TH1F("h1","x distribution",200,b,a);
  in.close();
  f->Write();
}

void world_hist(){
  //TFile* f = new TFile("world.root");
  //TTree* mytree = (TTree*)f->Get("mytree");
  TChain* chain = new TChain("tree");
  chain->Add("world.root");
  static Float_t px;
  Double_t hx;
  //TBranch* hx = mytree->GetBranch("px");
  //hx->SetAddress(&px);
  chain->SetBranchAddress("px",&px);
  Long64_t nentries = chain->GetEntries();
  for(Long64_t n =0;n<nentries;n++){
    chain->GetEntry(n);
    hx=px;
    //if (n < 5) printf("hx=%8f\n",hx);
  } 
  TH1F *h1 = new TH1F("h1","The  distribution of the target",1200,0,a);


  for (Long64_t i=0;i<nentries;i++){
    chain->GetEntry(i);
    //if (i < 5) printf("hx=%8f\n",hx);
    hx=px;
    h1->Fill(px);
  }
  Double_t mean=h1->GetMean();
  printf("the mean number is %f\n",mean);
  
  /*TCanvas* c = new TCanvas("c","Distribution roofit",1024,800);
  h1->Draw();
  gPad->SetLogy();*/
}
void world(){
  world_tree();
  world_hist();
}
