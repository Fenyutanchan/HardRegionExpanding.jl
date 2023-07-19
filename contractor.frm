
symbol diim, epsilon;
dimension diim;

vector q1,q2,q3,q4;
vector k1,...,k10;
vector K1,...,K10;
vector mom,mom0,mom1,...,mom10;
vector unity;

Set Kn: K1,...,K10;
Set kn: k1,...,k10;

set LOOP: q1,q2,q3,q4;
set NonLOOP: k1,...,k10,K1,...,K10;
set MOM: k1,...,k10,K1,...,K10,q1,q2,q3,q4;

symbol pi, im, sqrt2, shat;
auto symbol ver, gc;

***Dirac indices
auto symbol spa, spb, spv spw;

***Lorentz indices
auto index mua, mub, muv, muw, rho, dum, epMU;

Set RHO: rho,rho0,...,rho200;
Set EPMU: epMU1,...,epMU100;
Set DUM: dum1,...,dum100;
Set NonEPMU: mua,mua0,...,mua100,mub,mub0,...,mub100,muv, muv0,...,muv10000,muw, muw0,...,muw10000,dum1,...,dum100;

Set LOR: mua,mua0,...,mua100,mub,mub0,...,mub100,
         muv,muv0,...,muv10000,muw,muw0,...,muw10000,
         rho,rho0,...,rho100,dum,dum0,...,dum100,
         epMU1,...,epMU100;


***For matching both indices and momenta
auto index var;
***For matching numbers
auto symbol int, setint;

symbol mass, mass1,...,mass4, width;

CFunction Den, PowDen, FV, FermionLoopPow, GhostLoopPow;
CFunction VecEpsilon, VecEpsilon1,...,VecEpsilon4, VecEp, VecEpC;
CFunction Spinor, Spinor1,...,Spinor4, FermionChain;

CFunction SpUB, SpVB, SpU, SpV;
Set SPSET: SpUB, SpVB, SpU, SpV;
Set LSPSET: SpUB, SpVB;
Set RSPSET: SpU, SpV;

CFunction UB, VB, U, V;
Set INSPSET: UB, VB, U, V;
Set ILSPSET: UB, VB;
Set IRSPSET: U, V;

CFunction SP(symmetric), LMT(symmetric), Levi(antisymmetric);
CFunction SPC(symmetric);
CFunction GAij, PLij, PRij, ONEij(symmetric), Trace, Trace5, Trace5sym;
CFunction GA, GA5, PL, PR;

CFunction JJ(antisymmetric), FF(antisymmetric);

symbol dZ3x1, dZ3x2; 
symbol dtZ3x1, dtZ3x2; 
symbol dZ2tx1, dZ2tx2; 
symbol dZmtx1, dZmtx2; 
symbol dZ2x1, dZ2x2; 
symbol dZ1x1, dZ1x2; 
symbol dZ4x1, dZ4x2; 
symbol dtZ1x1, dtZ1x2; 
symbol dZ1Fx1, dZ1Fx2; 
symbol dZ1Ftx1, dZ1Ftx2; 
symbol dZgx1, dZgx2;
symbol d2s1, d2s2, d2s3, d2s4;
symbol f2s1, f2s2, f2s3, f2s4;




*----------------------------------------
*** without expansion, we do not expand momentum in FV or GA.
#procedure SimplificationNoExpand()

id im^2 = -1;
id sqrt2^2 = 2;
id sqrt2^(-2) = 1/2;

* FeynCalc has definition GA[6]=(1+GA[5])/2, GA[7]=(1-GA[5])/2
* FORM has a little different g_(6)=1+g_(5), g_(7)=1-g_(5)
* According to textbooks (Peskin's QFT etc.) PL=GA(7), PR=GA(6)
* PR-->tagJ, (1+GA5)/2, GA(6)
* PL-->tagF, (1-GA5)/2, GA(7)
* Same as definition in FeynCalc
*
***id PR(spa1?,spa2?) = GA(spa1,spa2,6);
***id PL(spa1?,spa2?) = GA(spa1,spa2,7);

*
* trivial Dirac indices contraction
*
repeat;
  id ONEij(spa?,spa?) = 4;
  id ONEij(spa1?,spa2?)*ONEij(spa2?,spa3?) = ONEij(spa1,spa3);

  id GAij(spa1?,spa2?,rho?LOR) * GAij(spa2?,spa3?,rho?LOR) = diim*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,mom?MOM) * GAij(spa2?,spa3?,mom?MOM) = SP(mom,mom)*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,rho?)*ONEij(spa3?,spa2?) = GAij(spa1,spa3,rho);
  id GAij(spa2?,spa1?,rho?)*ONEij(spa3?,spa2?) = GAij(spa3,spa1,rho);

  id Spinor?SPSET(int?,spa1?,mom?,mass?)*ONEij(spa1?,spa2?) = Spinor(int,spa2,mom,mass);
endrepeat;
.sort
*
* trivial Lorentz indices contraction
*
repeat;
  id LMT(rho1?,rho1?) = diim;
  id LMT(rho1?,rho2?)^2 = diim;
  id LMT(rho1?,rho2?)*LMT(rho2?,rho3?) = LMT(rho1,rho3);

  id LMT(rho1?,rho2?)*FV(mom1?,rho1?) = FV(mom1,rho2);
  id LMT(rho1?,rho2?)*GAij(spa1?,spa2?,rho1?) = GAij(spa1,spa2,rho2);
  id LMT(rho1?,rho2?)*Levi(rho1?,rho3?,rho4?,rho5?) = Levi(rho2,rho3,rho4,rho5);
  id LMT(rho1?,rho2?)*VecEpsilon?{VecEp,VecEpC}(int1?,rho2?,mom1?,mass1?) = VecEpsilon(int1,rho1,mom1,mass1);
endrepeat;

id GAij(spa1?,spa2?,rho?)*FV(mom?,rho?) = GAij(spa1,spa2,mom);
id FV(mom1?,rho?)*FV(mom2?,rho?) = SP(mom1,mom2);

repeat;
  id FV(mom?,rho?)*Levi(var1?,var2?,var3?,rho?) = Levi(var1,var2,var3,mom);
  id LMT(rho1?,rho2?)*FermionChain(?vars1,GA(rho2?),?vars2) = FermionChain(?vars1,GA(rho1),?vars2);
endrepeat;

repeat id LMT(rho?NonEPMU,rho0?)*VecEpsilon?{VecEp,VecEpC}(int?,rho?NonEPMU,?vars)
  = LMT(EPMU[int],rho0)*VecEpsilon(int,EPMU[int],?vars);
id GAij(spa1?,spa2?,rho?NonEPMU)*VecEpsilon?{VecEp,VecEpC}(int?,rho?NonEPMU,?vars)
  = GAij(spa1,spa2,EPMU[int])*VecEpsilon(int,EPMU[int],?vars);
.sort

*
* vanishing 
*
id FV(mom?,rho?) * VecEpsilon?{VecEp,VecEpC}(int?,rho?,mom?,mass?) = 0;
.sort

*
* use EPMU indices
*
id SP( FV(mom1?,0), VecEpsilon?{VecEp,VecEpC}(int?,0,mom2?,mass?) )
  = FV(mom1,EPMU[int]) * VecEpsilon(int,EPMU[int],mom2,mass);
id FV(mom1?,rho?) * VecEpsilon?{VecEp,VecEpC}(int?,rho?,mom2?,mass?)
  = FV(mom1,EPMU[int]) * VecEpsilon(int,EPMU[int],mom2,mass);
id VecEpsilon1?{VecEp,VecEpC}(int1?, rho?, ?vars1) * VecEpsilon2?{VecEp,VecEpC}(int2?, rho?, ?vars2)
  = SP( VecEpsilon1(int1, 0, ?vars1), VecEpsilon2(int2, 0, ?vars2) );
.sort

*
*Linearly expand momentum polynomial in FV and GA
*
***id FV(var?,rho?) = FV(var,rho);
***id GAij(spa1?,spa2?,var?) = GAij(spa1,spa2,var);
***id SP(rho1?,rho2?) = SP(rho1,rho2);
***id Levi(rho1?,rho2?,rho3?,rho4?) = Levi(rho1,rho2,rho3,rho4);
***.sort

*
* vanishing momentum scalar product also vanishes
*
id SP(mom?kn,mom?kn) = 0;

*
* Explain FermionLoopPow and GhostLoopPow
*
id FermionLoopPow(-1,int?) = (-1)^int;
id GhostLoopPow(-1,int?) = (-1)^int;
.sort


#endprocedure








*----------------------------------------
*** without expansion, the commutative relation between PL/PR and GA is uncertain, since GA[mom] may have mass*unity here.
#procedure contractDiracIndicesNoExpand()

* FeynCalc has definition GA[6]=(1+GA[5])/2, GA[7]=(1-GA[5])/2
* FORM has a little different g_(6)=1+g_(5), g_(7)=1-g_(5)
* According to textbooks (Peskin's QFT etc.) PL=GA(7), PR=GA(6)
* PR-->tagJ, (1+GA5)/2, GA(6)
* PL-->tagF, (1-GA5)/2, GA(7)
* Same as definition in FeynCalc
*
***id PR(spa1?,spa2?) = GA(spa1,spa2,6);
***id PL(spa1?,spa2?) = GA(spa1,spa2,7);

*
* trivial Dirac indices contraction
*
repeat;
  id ONEij(spa?,spa?) = 4;
  id ONEij(spa1?,spa2?)*ONEij(spa2?,spa3?) = ONEij(spa1,spa3);

  id GAij(spa1?,spa2?,rho?LOR) * GAij(spa2?,spa3?,rho?LOR) = diim*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,mom?MOM) * GAij(spa2?,spa3?,mom?MOM) = SP(mom,mom)*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,rho?)*ONEij(spa3?,spa2?) = GAij(spa1,spa3,rho);
  id GAij(spa2?,spa1?,rho?)*ONEij(spa3?,spa2?) = GAij(spa3,spa1,rho);

  id PLij(spa1?,spa2?)*ONEij(spa2?,spa3?) = PLij(spa1,spa3);
  id PRij(spa1?,spa2?)*ONEij(spa2?,spa3?) = PRij(spa1,spa3);

  id PLij(spa1?,spa2?)*PLij(spa2?,spa3?) = PLij(spa1,spa3);
  id PRij(spa1?,spa2?)*PRij(spa2?,spa3?) = PRij(spa1,spa3);

  id PLij(spa1?,spa2?)*PRij(spa2?,spa3?) = 0;
  id PRij(spa1?,spa2?)*PLij(spa2?,spa3?) = 0;


  id Spinor?SPSET(int?,spa1?,mom?,mass?)*ONEij(spa1?,spa2?) = Spinor(int,spa2,mom,mass);
endrepeat;
.sort

***id GAij(spa1?,spa2?,var?) = GAij(spa1,spa2,var);
***.sort



*
* now chainin the Dirac objects in FermionChain according to Dirac indices
*
***Here var? could be mom? or rho?
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa1?,spa2?,mom?) = FermionChain( ILSPSET[setint](int,?vars), GA(mom), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa1?,spa2?,rho?LOR) = FermionChain( ILSPSET[setint](int,?vars), GA(rho), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PLij(spa1?,spa2?) = FermionChain( ILSPSET[setint](int,?vars), PL, spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PRij(spa1?,spa2?) = FermionChain( ILSPSET[setint](int,?vars), PR, spa2 );

***flip
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa2?,spa1?,mom?) = -FermionChain( ILSPSET[setint](int,?vars), GA(mom), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa2?,spa1?,rho?LOR) = -FermionChain( ILSPSET[setint](int,?vars), GA(rho), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PLij(spa2?,spa1?) = FermionChain( ILSPSET[setint](int,?vars), PL, spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PRij(spa2?,spa1?) = FermionChain( ILSPSET[setint](int,?vars), PR, spa2 );


id Spinor1?LSPSET[setint1](int1?,spa?,?var1)*Spinor2?RSPSET[setint2](int2?,spa?,?var2) = FermionChain( ILSPSET[setint1](int1,?var1), IRSPSET[setint2](int2,?var2) );
.sort

repeat;
***Here var? could be mom? or rho?
  id FermionChain(?vars,spa1?)*GAij(spa1?,spa2?,?expr) = FermionChain( ?vars, GA(?expr), spa2 );
  id FermionChain(?vars,spa1?)*GAij(spa1?,spa2?,rho?LOR) = FermionChain( ?vars, GA(rho), spa2 );

  id FermionChain(?vars,spa1?)*PLij(spa1?,spa2?) = FermionChain( ?vars, PL, spa2 );
  id FermionChain(?vars,spa1?)*PRij(spa1?,spa2?) = FermionChain( ?vars, PR, spa2 );

***flip
  id FermionChain(?vars,spa1?)*GAij(spa2?,spa1?,?expr) = -FermionChain( ?vars, GA(?expr), spa2 );
  id FermionChain(?vars,spa1?)*GAij(spa2?,spa1?,rho?LOR) = -FermionChain( ?vars, GA(rho), spa2 );

  id FermionChain(?vars,spa1?)*PLij(spa2?,spa1?) = FermionChain( ?vars, PL, spa2 );
  id FermionChain(?vars,spa1?)*PRij(spa2?,spa1?) = FermionChain( ?vars, PR, spa2 );
endrepeat;

id FermionChain(?vars1,spa?)*Spinor?RSPSET[setint](int?,spa?,?vars2) = FermionChain( ?vars1, IRSPSET[setint](int,?vars2) );
.sort


*
* Look for Trace
*
repeat;
  id once, GAij(spa1?,spa2?,?expr) = Trace(GA(?expr),spa1,spa2);
  repeat;
    id Trace(?vars,spa1?,spa2?)*GAij(spa2?,spa3?,?expr) = Trace(?vars,GA(?expr),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*GAij(spa2?,spa3?,rho?LOR) = Trace(?vars,GA(rho),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PLij(spa2?,spa3?) = Trace(?vars,PL,spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PRij(spa2?,spa3?) = Trace(?vars,PR,spa1,spa3);

    id Trace(?vars,spa1?,spa2?)*GAij(spa3?,spa2?,?expr) = -Trace(?vars,GA(?expr),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*GAij(spa3?,spa2?,rho?LOR) = -Trace(?vars,GA(rho),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PLij(spa3?,spa2?) = Trace(?vars,PL,spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PRij(spa3?,spa2?) = Trace(?vars,PR,spa1,spa3);
  endrepeat;
  id Trace(?vars,spa?symbol_,spa?symbol_) = Trace(?vars);
  id PLij(spa?symbol_,spa?symbol_) = 2;
  id PRij(spa?symbol_,spa?symbol_) = 2; 
  id ONEij(spa?symbol_,spa?symbol_) = 4;
endrepeat;
.sort

*** #call ArrangeTrace();

***move PL and PR to left-side of FermionChain, right-side of spinor
repeat;
  id FermionChain( ?vars1, PL, PR, ?vars2 ) = 0;
  id FermionChain( ?vars1, PR, PL, ?vars2 ) = 0;
  id FermionChain( ?vars1, PL, PL, ?vars2 ) = FermionChain( ?vars1, PL, ?vars2 );
  id FermionChain( ?vars1, PR, PR, ?vars2 ) = FermionChain( ?vars1, PR, ?vars2 );

**id FermionChain( ?vars1, GA(mom?), PL, ?vars2 ) = FermionChain( ?vars1, PR, GA(mom), ?vars2 );
**id FermionChain( ?vars1, GA(mom?), PR, ?vars2 ) = FermionChain( ?vars1, PL, GA(mom), ?vars2 );

**id FermionChain( ?vars1, GA(rho?LOR), PL, ?vars2 ) = FermionChain( ?vars1, PR, GA(rho), ?vars2 );
**id FermionChain( ?vars1, GA(rho?LOR), PR, ?vars2 ) = FermionChain( ?vars1, PL, GA(rho), ?vars2 );
endrepeat;
.sort


repeat;
  id FV(mom?,rho?LOR)*FermionChain(?vars1,GA(rho?LOR),?vars2) = FermionChain(?vars1,GA(mom),?vars2);
  id LMT(rho1?LOR,rho2?LOR)*FermionChain(?vars1,GA(rho1?LOR),?vars2) = FermionChain(?vars1,GA(rho2),?vars2);
endrepeat;
.sort

repeat;
  id FermionChain(?vars1,GA(mom?),GA(mom?),?vars2) = SP(mom,mom)*FermionChain(?vars1,?vars2);
  id SP(mom?kn,mom?kn) = 0;
  id FermionChain(?vars1,GA(rho?LOR),GA(rho?LOR),?vars2) = diim*FermionChain(?vars1,?vars2);

*** Dirac equation for U and V (UB and VB) 
  id FermionChain( ?vars, GA(mom?), U(int?,mom?,0) ) = 0;
  id FermionChain( ?vars, GA(mom?), V(int?,mom?,0) ) = 0;

  id FermionChain( UB(int?,mom?,0), GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,0), GA(mom?), ?vars ) = 0;

  id FermionChain( UB(int?,mom?,0), PL, GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,0), PL, GA(mom?), ?vars ) = 0;

  id FermionChain( UB(int?,mom?,0), PR, GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,0), PR, GA(mom?), ?vars ) = 0;
endrepeat;
.sort

#endprocedure




