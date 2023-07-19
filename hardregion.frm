#-
***我们需要一个脚本来把质量区分轻重
***set setlightmasses: mmH,mH;(针对4W试运行)
#define Nmax "4"
***该项是为了进行硬动量展开预先标记振幅分子上的power
CF invDen,den,Mom,unitary,DEN,Log;
s lin,lilin,Epson1,m,m1,m2,mmH,mH,num1,num2,D,I,FourPi,inverEpson1;
set setheavymasses: mmH,mH;
set setlightmasses: 0;
#define n "4"
Index itrick;
***为了做tensor reduction预放置lortenz定义
i <LLor1=D>, ... , <LLor200=D>;
set LLori: LLor1,...,LLor200;
i <splor1=D>, ... , <splor200=D>;
set splori: splor1,...,splor200;
***Form程序会自动的用id把Dot（p,k1+k2）展开,因此我们需要进行些预处理
#procedure expandmom
id SP(k1?,itrick?)=SP(k1,itrick);
id SP(itrick?,k1?)=SP(itrick,k1);
#do jf=1,200;
id,once FermionChain( ?vars1, GA(k1?), ?vars2 )= Mom(k1,LLor'jf')*FermionChain( ?vars1, GA(LLor'jf'), ?vars2 );
#enddo;
id Mom(itrick?,mua1?)=Mom(itrick,mua1);
#endprocedure


#procedure preparenumerator
id Mom(k1?NonLOOP,mua1?)=lin*Mom(k1,mua1);
id SP(k1?NonLOOP,k2?LOOP)=lin*SP(k1,k2);
id SP(k1?LOOP,k2?NonLOOP)=lin*SP(k1,k2);
id SP(k1?NonLOOP,k2?NonLOOP)=lin*lin*SP(k1,k2);
id lin^('Nmax'+1)=0;
#endprocedure
***这样以来我们已经对分子进行了预备处理,这样先从分子上挑出dimension低于Nmax的term了
#procedure propagatoridea
id Den(q1?,m?,0)=Den(q1,q1,m,0);
argument Den,2;
id q1=0;
endargument;
argument Den,1;
#do vari=1,`n';
id k'vari'=0;
#enddo
endargument;
***我们将Den变为Den（LOOP动量,外腿动量,质量,维数展开项目）,现在我们需要进行处理,进行硬动量展开,分为两类,***重质量和轻质量,SMEFT一般设置轻质量为0
id Den(-q1,k1?,m?,0) = Den(q1,-k1,m,0);
id Den(-q1,0,m?,0) = Den(q1,0,m,0);
***上式是为了识别
id m?setlightmasses=lin*m;
***这一步的实现原因是这样的因为轻质量粒子的硬动量展开分子上会带有m^2,需要在计算中进行区分检查
repeat;
* First for the heavy masses
id,once Den(q1,k1?,m?setheavymasses,0) = Den(q1,0,m,0)*(1-(lin^2*SP(k1,k1)+2*lin*SP(k1,q1))*Den(q1,k1,m,0));
id,once Den(0,k1?,m?setheavymasses,0) = Den(0,0,m,0)*(1-(lin^2*SP(k1,k1))*Den(0,k1,m,0));
* Now for the light masses
id,once Den(q1,k1?,m?!setheavymasses,0) = Den(q1,0,0,0)*(1-(lin^2*SP(k1,k1)+2*lin*SP(k1,q1)-lin^2*m^2)*Den(q1,k1,m,0));
id,once Den(k1,0,m?!setheavymasses,0) = Den(k1,0,0,0)*(1-(-lin^2*m^2)*Den(k1,0,m,0));
id,once lin^('Nmax'+1)=0;
endrepeat;
***究其本质这一项之所以那么进行改造实际上只是为了未来我们适配任何dim4+Nmax算符对应振幅的架构,要是觉得麻烦,实***际上也不会影响计算速度,你可以注释掉写入自己的替换。
id lin*m?setlightmasses=m;
id Den(q1?,0,m?,0)=Den(q1,m,0);
id Den(0,0,m?,0)=Den(0,m,0);
.sort
***我们已经对分子先做了全展开,挑出小于等于lin^Nmax的term,同时在此基础进一步进行了硬动量展开,但是我们注意到***写出Qgraf文件的时候传播子动量非LOOPterm有可能是k1+k2,再做一遍展开就行了。
id SP(k1?,k2?)=SP(k1,k2);
.sort
id Den(0,0,0)=0;
id Den(0,m,0)=-1/m^2;
.sort
#endprocedure

***这一步进行完毕之后我们需要进行的计算很明白就是把分子上因为硬动量展开带来的包含LOOP momentum的多项缩并,这一***步同样李钊老师的程序里面没有需要我们自己做,因此我们需要考虑到,form的顺序执行问题,所以从8-6-4-2个LOOP张***量reduction程序并且撰写保持,并在最后结果写回来
***警告：我们前面对于完整振幅的分子带有gamma矩阵的同样没有做为这样方便的展开,因此我们需要同样做预备处理,最后回***来,当然因为前面已经做了展开,我们只需要考虑LOOPterm

***gamma矩阵,准确的来说我们采取非thooft方案
#procedure gamma5
repeat;
id FermionChain( ?vars1, GA(mom?), PL, ?vars2 ) = FermionChain( ?vars1, PR, GA(mom), ?vars2 );
id FermionChain( ?vars1, GA(mom?), PR, ?vars2 ) = FermionChain( ?vars1, PL, GA(mom), ?vars2 );
id FermionChain( ?vars1, PR, PR, ?vars2 ) = FermionChain( ?vars1, PR, ?vars2 );
id FermionChain( ?vars1, PL, PL, ?vars2 ) = FermionChain( ?vars1, PL, ?vars2 );
id FermionChain( ?vars1, PL,PR, ?vars2 ) = 0;
endrepeat;
***所有PL和PR都挪到左边
#endprocedure
#procedure tensorreduction
*#do jf=1,200;
*id,once GA(k1?) = Mom(k1,LLor'jf')*GA(LLor'jf');
*#enddo;
#do jf=1,200;
id,once SP(k1?,k2?)=Mom(k1,splor'jf')*Mom(k2,splor'jf');
#enddo;
id Mom(q1,mua1?) = lilin*Mom(q1,mua1);
repeat;
id lilin^2=1;
endrepeat;
id lilin=0;
***扔掉奇数项
.sort

id Mom(q1,mua1?)^2=SP(q1,q1);
.sort


id Mom(q1, mua1?)*Mom(q1, mua2?)*
   Mom(q1, mua3?)*Mom(q1, mua4?)*
   Mom(q1, mua5?)*Mom(q1, mua6?)*
   Mom(q1, mua7?)*Mom(q1, mua8?)=
   dd_(mua1, mua2, mua3, mua4, mua5, mua6, mua7, mua8)*
   SP(q1,q1)*SP(q1,q1)*SP(q1,q1)*SP(q1,q1)*105*den(201600,258720);
id Mom(q1, mua1?)*Mom(q1, mua2?)*
   Mom(q1, mua3?)*Mom(q1, mua4?)*
   Mom(q1, mua5?)*Mom(q1, mua6?)=
   dd_(mua1, mua2, mua3, mua4, mua5, mua6)*
   SP(q1,q1)*SP(q1,q1)*SP(q1,q1)*15*den(2880,3120);
id Mom(q1, mua1?)*Mom(q1, mua2?)*
   Mom(q1, mua3?)*Mom(q1, mua4?)=
   dd_(mua1, mua2, mua3, mua4)*
  SP(q1,q1)*SP(q1,q1)*3*den(72,60);
id Mom(q1, mua1?)*Mom(q1, mua2?)= dd_(mua1, mua2)*SP(q1,q1)*den(4,2);

***现在我们已经做完去了全部的拆解,现在需要做一件事就是把所有的转回李钊老师的结构做缩并
repeat;
id Mom(k1?NonLOOP,mua1?)*FermionChain(?vars1,GA(mua1?),?vars2)= FermionChain(?vars1,GA(k1),?vars2);
endrepeat;

repeat;
id Mom(k1?,mua1?)*Mom(k2?,mua1?)=SP(k1,k2);
id Mom(k1?,mua1?)*Mom(k1?,mua1?)=SP(k1,k1);
endrepeat;
***涉及到gamma矩阵有关的缩并,李钊老师没有提供,我们需要手动提供
repeat;
id FermionChain(?vars1,GA(LLor1?LLori),GA(LLor1?LLori),?vars2)=D*FermionChain(?vars1,?vars2);
id  FermionChain(?vars1,GA(LLor1?LLori), GA(k1?NonLOOP),GA(LLor1?LLori),?vars2)=(2-D)*FermionChain(?vars1,GA(k1),?vars2);
id  FermionChain(?vars1,GA(LLor1?LLori), GA(k1?NonLOOP),GA(k2?NonLOOP),GA(LLor1?LLori),?vars2)=4*SP(k1,k2)*FermionChain(?vars1,?vars2)*unitary+(D-4)* FermionChain(?vars1,GA(k1),GA(k2),?vars2);
id  FermionChain(?vars1,GA(LLor1?LLori), GA(k1?NonLOOP),GA(k2?NonLOOP),GA(k3?NonLOOP),GA(LLor1?LLori),?vars2)=-2*FermionChain(vars1,GA(k3),GA(k2),GA(k1),vars2)-(D-4)*FermionChain(?vars1,GA(k1),GA(k2),GA(k3),?vars2);
id  FermionChain(?vars1,GA(LLor1?LLori)* GA(k1?NonLOOP)*GA(k2?NonLOOP)*GA(k3?NonLOOP)*GA(k4?NonLOOP)*GA(LLor1?LLori),?vars2)=2*FermionChain(?vars1,GA(k3),GA(k2),GA(k1),GA(k4),?vars2)+2*FermionChain(?vars1,GA(k4),GA(k1),GA(k2),GA(k3),?vars2)+(D-4)*FermionChain(?vars1,GA(k1),GA(k2),GA(k3),GA(k4),?vars2);
endrepeat;
repeat;
id SP(-k1?NonLOOP,k2?NonLOOP)=-SP(k1,k2);
id SP(k1?NonLOOP,-k2?NonLOOP)=-SP(k1,k2);
id SP(-k1?NonLOOP,-k2?NonLOOP)=SP(k1,k2);
id FermionChain(?vars1,GA(-k1?NonLOOP),?vars2)=-FermionChain(?vars1,GA(k1),?vars2);
endrepeat;
#endprocedure
***李钊老师程序,直接调用费米子链缩并

***on shell条件
***#include kin_relation.frm
***圈图计算,这里在函数中会引用宏定义而非固定语法,是为了预留RGE的接口

#procedure partialfractioning(isnoteft,LOOPorder)
id SP(q1,q1) = invDen(q1, 0,0);
id Den(q1, 0,0)*invDen(q1, 0,0)=1;


repeat;
id,once Den(q1,m1?,0)*Den(q1,m2?!{m1?},0)= DEN(m1,m2)*(Den(q1,m1,0)- Den(q1,m2,0));
endrepeat;


id DEN(m1?,m2?)=1/(m1^2-m2^2);
*id DEN(m1?,0)=1/m^2;


*前面区分了展开,接下来我们需要考虑分子
repeat;
id invDen(q1,0,0)*Den(q1,m1?,0)=1+m1^2*Den(q1,m1,0);
endrepeat;

* 这个地方是我们调用宏命令修饰整个函数的,为了之后rge准备,其中isnoteft这个需要用脚本写,这里还做不到,当然还是预留接口
#if ((isnoteft==0) && (LOOPorder==1))
id only Den(q1,0,0)^2=I/(FourPi^2*Epson1);
* everything else is 0
multiply Epson1;
id Epson1=0;
multiply 1/Epson1;
#endif
*-----  zero of dimreg.
id Den(q1, 0,0)=0;
id invDen(q1, 0,0)=0;

.sort


#endprocedure
***我们做了什么事情？其实就是把所有的硬动量展开后首先是分母分离,然后再把分子上的q2/（q2-m2）项分离开来
***费米子链的contract需要调用李钊老师contract.frm的程序,已经在上面写出。
***最后一项是IBP term+计算有限项目
#procedure IBPtechinque
#do jc=0,8
id,once Den(q1,m?,0)^(10-'jc') =
Den(q1, m,0)^(9-'jc')/m^2*(D/2-(9-'jc'))/(9-'jc');
#enddo;
#endprocedure

***注意下procedure已经进行了1/k^4在非rge下的调用区分（已经做了define数值）,以及圈积分的IBP,下面是evaluate有限项的值,我们的IBP根据前面的工作仅仅进行到1/k^2-m^2为止。
#procedure Paxevaluate
id D =4-2*Epson1;
***注意到表达式的一般格式为：poly1（D）／poly2（D）。poly（D）来源自我们还未定义具体值的den函数,因此只有poly1会显式出值,考虑到poly2（D）值提供的term,可以定义Epson1^2为0先
repeat;
id Epson1^2=0;
id den(num1?,num2?)=1/num1+num2/num1*Epson1*den(num1, num2);
endrepeat;
***在RGEmod下我们必须预留一个区分k^2-m^2和k^2的预置代码，这里还没有添加
id Den(q1,m1?,0)=I*m1^2/FourPi^2*(1+inverEpson1-Log(m1^2));
id inverEpson1*Epson1=1;
id Epson1*inverEpson1=0;
id Epson1=0;

#endprocedure
***上述表达式其实就是做一个展开而已,但是是用恒等式的形式展开的,自己写写看就知道了,到最后所有的den项都会和epson1^2有关系进而为0
***包括收集洛伦兹结构这些表达式已经写在了李钊老师的其他文件里,直接调用就好。
***我们需要定义一个总程序给这个文件进行收尾
#procedure hardregion(isnoteft,LOOPorder)
#call expandmom
#call preparenumerator
#call propagatoridea
#call tensorreduction
***#include kin_relation.frm
#call partialfractioning(isnoteft,LOOPorder)
#call IBPtechinque
#call Paxevaluate
#include color.frm
#endprocedure








