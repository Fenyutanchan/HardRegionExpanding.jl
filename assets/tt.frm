#-
#include hardregion.frm
#include contractor.frm
Local expr=Den(q1,mmH,0)*Den(q1-k2,0,0)*Den(q1-k1,0,0)*Den(q1,mmH,0);
*Local expr=den(72,60);
#call propagatoridea
#include kin_relation.frm
id SP(-k1,-k1)=0;
#call tensorreduction
#include kin_relation.frm
id SP(-k1,-k1)=0;
.sort
print;
#call partialfractioning(1,1)
#call IBPtechinque
#call Paxevaluate
#include kin_relation.frm
id Log(m1?)=0;
id inverEpson1=0;
print;
.end

