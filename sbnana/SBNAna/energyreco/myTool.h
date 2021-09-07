#ifndef myTool_h
#define myTool_h

vector<string> removeZombie(const vector<string>& vec_ori){

  cout << "[removeZombie] Input number of files = " << vec_ori.size() << endl;

  vector<string> vec_out;
  for(unsigned int i=0; i<vec_ori.size(); i++){
    TString this_filename = vec_ori.at(i);
    TFile *file = new TFile(this_filename);
    if(file->IsZombie()){
      cout << "[removeZombie] " << vec_ori.at(i) << " is a zombie" << endl;
    }
    else{
      vec_out.push_back( vec_ori.at(i) );
    }
    file->Close();
  }
  return vec_out;

}

#endif
