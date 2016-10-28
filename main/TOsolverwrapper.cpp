
int TOSolver(unsigned int NmaxTEBD,unsigned int chiDMRG,unsigned int maxstepsDMRG, unsigned int maxstepsSSDMRG, char *sname, unsigned int step, unsigned int chimax, double dt, double tr_weight, int spin, int branch, int orbital, bool verbose,double deltaDMRG,int Numfiles);

int main(const int argc,const char **argv){
  TOSolver(200,400,10, 10,"test_label",1,50,-0.05,1e-10,1,1,1,true,1e-10,10);
}


