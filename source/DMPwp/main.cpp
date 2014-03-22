
#include "Dmplwr.h"

int main()
{
  Dmplwr hdmp(1001,1,0.001);                
  hdmp.Load("../../data/");
  hdmp.Learnlwr();                
  hdmp.Reprolwr(1,0.2,"../../data/lwr/");

  return 1;
}
