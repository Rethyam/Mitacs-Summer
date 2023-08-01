clear
clc
%%%%%%%%%%%%% Quantum circuit Fontanella et al Appendix B
tic
load BSM.mat
%theta0=2*pi*rand(25,1); % Random initialisation 
theta0=[3.142/2; 4.173/2; 1.392/2; 3.713/2; 2.399/2; .935/2; 2.196/2; 5.014/2; 2.736/2; 1.477/2; 4.472/2; 3.415/2; 6.283/2; 4.244/2; 4.711/2; .717/2; 1.741/2; 1.158/2; 2.531/2; 5.705/2; 3.525/2; 4.582/2; 2.465/2; .098/2; 5.018/2];
%Initialisation with Fontanella's solution
%theta0=[3.142; 4.173; 1.392; 3.713; 2.399; .935; 2.196; 5.014; 2.736; 1.477; 4.472; 3.415; 6.283; 4.244; 4.711; .717; 1.741; 1.158; 2.531; 5.705; 3.525; 4.582; 2.465; .098; 5.018];
%thetaT=[3.203; 4.149; 2.278; 3.512; 2.4; 0.512; 2.578; 5.024; 2.405; 2.406; 4.271; 3.375; 6.331; 4.224; 4.347; 0.510; 2.357; 0.956; 2.443; 5.232; 3.184; 4.716; 2.706;-1.154; 4.817];
% construction of psiInit, in fact psi0
xmin=log(50);
xmax=log(150);
incr=(xmax-xmin)/16;
range=xmin:incr:xmax;
range=range';
range=range(1:16);
a=.5; % We suppose r=0
K=100;
BSMPayoff=max(exp(range)-K,0);
psiInit=exp(-a*range).*max(exp(range)-K,0);
gamma0=norm(psiInit);
psiInit=psiInit/gamma0;
%norm(psiInit) =1
load Hamiltonian.mat
%%%% Ansatz: Minimization on angles of circuits for production of psiInit 
fun0=@(theta)(norm(psiInit-Phi(theta).Amplitudes));
[thetaM fval]=fminunc(fun0,theta0);
%thetaM=theta0;
Payoff0=Phi(thetaM).Amplitudes*gamma0;
Payoff0=exp(a*range).*Payoff0;
norm(Payoff0-BSMPayoff)
deltatau=1/500;
thetaMin=zeros(25,502);
thetaMin(:,1)=thetaM;
for t=1:500
    thetaMinCurr=thetaMin(:,t);
    %Construction of A(t)
    v=zeros(25,16);
    v(1,:)=Phi1(thetaMinCurr).Amplitudes;
    v(2,:)=Phi2(thetaMinCurr).Amplitudes;
    v(3,:)=Phi3(thetaMinCurr).Amplitudes;
    v(4,:)=Phi4(thetaMinCurr).Amplitudes;
    v(5,:)=Phi5(thetaMinCurr);
    v(6,:)=Phi6(thetaMinCurr);
    v(7,:)=Phi7(thetaMinCurr);
    v(8,:)=Phi8(thetaMinCurr).Amplitudes;
    v(9,:)=Phi9(thetaMinCurr).Amplitudes;
    v(10,:)=Phi10(thetaMinCurr).Amplitudes;
    v(11,:)=Phi11(thetaMinCurr).Amplitudes;
    v(12,:)=Phi12(thetaMinCurr);
    v(13,:)=Phi13(thetaMinCurr);
    v(14,:)=Phi14(thetaMinCurr);
    v(15,:)=Phi15(thetaMinCurr).Amplitudes;
    v(16,:)=Phi16(thetaMinCurr).Amplitudes;
    v(17,:)=Phi17(thetaMinCurr).Amplitudes;
    v(18,:)=Phi18(thetaMinCurr).Amplitudes;
    v(19,:)=Phi19(thetaMinCurr);
    v(20,:)=Phi20(thetaMinCurr);
    v(21,:)=Phi21(thetaMinCurr);
    v(22,:)=Phi22(thetaMinCurr).Amplitudes;
    v(23,:)=Phi23(thetaMinCurr).Amplitudes;
    v(24,:)=Phi24(thetaMinCurr).Amplitudes;
    v(25,:)=Phi25(thetaMinCurr).Amplitudes;
    for i=1:25
        for j=1:25
        A(i,j)=real(v(i,:)*v(j,:)');
        end
    end
    v0=Phi(thetaMinCurr).Amplitudes;
    Hv0=H*v0;
    for i=1:25
    C(i,1)=real(v(i,:)*Hv0);
    end
    % Advance 1 step
    %thetaPMin=pinv(A)*C;
    if t==1 thetaP=2*pi*rand(25,1);
    else thetaP=thetaPMin;
    end
    fun=@(thetaP)(norm(A*thetaP-C));
    
    % Trust-Region Reflective optimization
    %[thetaPMin fval(t)]=fminunc(fun,thetaP);

    % Nelder-Mead optimization
    % [thetaPMin, fval(t)] = fminsearch(fun, thetaP);

    % Genetic Algorithm options
    %options = optimoptions('ga', 'PopulationSize', 100, 'MaxGenerations', 100, 'FunctionTolerance', 1e-6);
    % Genetic Algorithm optimization
    %[thetaPMin, fval(t)] = ga(fun, 1, [], [], [], [], [], [], [], options);

    % Simulated Annealing options
    %options = optimoptions('simulannealbnd', 'MaxIterations', 100, 'FunctionTolerance', 1e-6);
    % Simulated Annealing optimization
    %[thetaPMin, fval(t)] = simulannealbnd(fun, thetaP, [], [], options);

    % Pattern Search options without bounds
    %options = optimoptions('patternsearch', 'MaxIterations', 100, 'FunctionTolerance', 1e-6);
    % Pattern Search optimization without bounds
    %[thetaPMin, fval(t)] = patternsearch(fun, thetaP, [], [], [], [], [], [], [], options);

    % PSO options without bounds
    %options = optimoptions('particleswarm', 'SwarmSize', 100, 'MaxIterations', 100, 'FunctionTolerance', 1e-6);
    % PSO optimization without bounds
    %[thetaPMin, fval(t)] = particleswarm(fun, 1, [], [], options);




    %[thetaPMin fval(t)]=fmincon(fun,thetaP,[],[],[],[],-5,5);
    % thetaPMin=pinv(A)*C;
    thetaPMinTemp(:,t)=thetaPMin;
    
    thetaMin(:,t+1)=thetaMinCurr+deltatau*thetaPMin;

    PayoffT=Phi(thetaMin(:,t)).Amplitudes*gamma0;
    PayoffT=exp(a*range).*PayoffT;
    figure(1)
    plot(PayoffT)
end

    hold
    plot(BSMPayoff)
toc
function outState=Phi(theta)

level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi1(theta) %The derivative w.r. to theta1,

% XZ=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);zGate(1);xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi2(theta)

level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);zGate(2);xGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi3(theta)

level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);zGate(3);xGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi4(theta)

level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);zGate(4);xGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function vectorOut=Phi5(theta) %The derivative w.r. to theta5,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);
%Preparation for derivative w.r. 2*theta5
vector2=inState2.Amplitudes;
Z=[1 0; 0 -1];
XZ=[0 -1;1 0];
Operator=-1/2*(kron(-eye(2),XZ)+kron(Z,XZ));
temp=kron(Operator,eye(4));
vectorDer=temp*vector2;

level2Gates=[idGate(1);ryGate(2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
M2=getMatrix(level2Circuit);
vector3=M2*vectorDer;

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
M3=getMatrix(level3Circuit);
vector4=M3*vector3;

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
M4=getMatrix(level4Circuit);
vector5=M4*vector4;

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
M5=getMatrix(level5Circuit);
vector6=M5*vector5;

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
M6=getMatrix(level6Circuit);
vector7=M6*vector6;

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
M7=getMatrix(level7Circuit);
vectorOut=M7*vector7;
end

function vectorOut=Phi6(theta) %The derivative w.r. to theta6,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);
levelXGates=[cryGate(1,2,2*theta(5));idGate(3);idGate(4)];
levelXCircuit=quantumCircuit(levelXGates);
inStateX=simulate(levelXCircuit,inState2);
vector2=inStateX.Amplitudes;
% Prep. for derivative w.r. to 2*theta6
Z=[1 0; 0 -1];
XZ=[0 -1;1 0];
Operator=-1/2*(kron(-eye(2),XZ)+kron(Z,XZ));
temp1=kron(Operator,eye(2));
temp2=kron(eye(2), temp1);
vectorDer=temp2*vector2;
level2Gates=[idGate(1);idGate(2);ryGate(3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
M2=getMatrix(level2Circuit);
vector3=M2*vectorDer;

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
M3=getMatrix(level3Circuit);
vector4=M3*vector3;

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
M4=getMatrix(level4Circuit);
vector5=M4*vector4;

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
M5=getMatrix(level5Circuit);
vector6=M5*vector5;

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
M6=getMatrix(level6Circuit);
vector7=M6*vector6;

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
M7=getMatrix(level7Circuit);
vectorOut=M7*vector7;
end

function vectorOut=Phi7(theta) %The derivative w.r. to theta7,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);
level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));idGate(4)];
level2Circuit=quantumCircuit(level2Gates);
inStateX=simulate(level2Circuit,inState2);
vectorX=inStateX.Amplitudes;
% Prep. for derivative w.r. to 2*theta7
Z=[1 0; 0 -1];
XZ=[0 -1;1 0];
Operator=-1/2*(kron(-eye(2),XZ)+kron(Z,XZ));
temp=kron(eye(4),Operator);
vectorDer=temp*vectorX;

level2Gates=[idGate(1); idGate(2);idGate(3);ryGate(4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
M2=getMatrix(level2Circuit);
vector3=M2*vectorDer;

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
M3=getMatrix(level3Circuit);
vector4=M3*vector3;

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
M4=getMatrix(level4Circuit);
vector5=M4*vector4;

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
M5=getMatrix(level5Circuit);
vector6=M5*vector5;

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
M6=getMatrix(level6Circuit);
vector7=M6*vector6;

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
M7=getMatrix(level7Circuit);
vectorOut=M7*vector7;
end

function outState=Phi8(theta) %The derivative w.r. to theta8,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[zGate(1);xGate(1);ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi9(theta) %The derivative w.r. to theta9,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));zGate(2);xGate(2);ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi10(theta) %The derivative w.r. to theta10,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);
% figure(2)
% plot(level2Circuit)
level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));zGate(3);xGate(3);ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi11(theta) %The derivative w.r. to theta11,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));zGate(4);xGate(4);ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function vectorOut=Phi12(theta) %The derivative w.r. to theta12,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);
vector4=inState4.Amplitudes;
% Prep. for derivative w.r. to theta12
Z=[1 0; 0 -1];
XZ=[0 -1;1 0];
Operator=-1/2*(kron(-eye(2),XZ)+kron(Z,XZ));
temp=kron(Operator,eye(4));
vectorDer=temp*vector4;

level4Gates=[idGate(1);ryGate(2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
M4=getMatrix(level4Circuit);
vector5=M4*vectorDer;

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
M5=getMatrix(level5Circuit);
vector6=M5*vector5;

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
M6=getMatrix(level6Circuit);
vector7=M6*vector6;

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
M7=getMatrix(level7Circuit);
vectorOut=M7*vector7;
end

function vectorOut=Phi13(theta) %The derivative w.r. to theta13,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);
% figure(2)
% plot(level2Circuit)
level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);
levelXGates=[cryGate(1,2,2*theta(12));idGate(3);idGate(4)];
levelXCircuit=quantumCircuit(levelXGates);
inStateX=simulate(levelXCircuit,inState4);
vectorX=inStateX.Amplitudes;
% Prep. for derivative w.r. to theta13
Z=[1 0; 0 -1];
XZ=[0 -1;1 0];
Operator=-1/2*(kron(-eye(2),XZ)+kron(Z,XZ));
temp1=kron(Operator,eye(2));
temp2=kron(eye(2),temp1);
vectorDer=temp2*vectorX;
level4Gates=[idGate(1); idGate(2);ryGate(3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
M4=getMatrix(level4Circuit);
vector5=M4*vectorDer;

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
M5=getMatrix(level5Circuit);
vector6=M5*vector5;

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
M6=getMatrix(level6Circuit);
vector7=M6*vector6;

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
M7=getMatrix(level7Circuit);
vectorOut=M7*vector7;
end

function vectorOut=Phi14(theta) %The derivative w.r. to theta14,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);
% figure(1)
% plot(level1Circuit)
level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);
% figure(2)
% plot(level2Circuit)
level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);
levelXGates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));idGate(4)];
levelXCircuit=quantumCircuit(levelXGates);
inStateX=simulate(levelXCircuit,inState4);
vectorX=inStateX.Amplitudes;
% Prep. for derivative w.r. to theta14
Z=[1 0; 0 -1];
XZ=[0 -1;1 0];
Operator=-1/2*(kron(-eye(2),XZ)+kron(Z,XZ));
temp=kron(eye(4),Operator);
vectorDer=temp*vectorX;
level4Gates=[idGate(1); idGate(2);idGate(3);ryGate(4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
M4=getMatrix(level4Circuit);
vector5=M4*vectorDer;

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
M5=getMatrix(level5Circuit);
vector6=M5*vector5;

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
M6=getMatrix(level6Circuit);
vector7=M6*vector6;

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
M7=getMatrix(level7Circuit);
vectorOut=M7*vector7;
end

function outState=Phi15(theta) %The derivative w.r. to theta15,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[zGate(1);xGate(1);ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi16(theta) %The derivative w.r. to theta16,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));zGate(2);xGate(2);ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi17(theta) %The derivative w.r. to theta17,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));zGate(3);xGate(3);ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi18(theta) %The derivative w.r. to theta18,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));zGate(4);xGate(4);ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);
% figure(6)
% plot(level6Circuit)
level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
% figure(7)
% plot(level7Circuit)
end
function vectorOut=Phi19(theta) %The derivative w.r. to theta19,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);
% figure(1)
% plot(level1Circuit)
level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);
% figure(2)
% plot(level2Circuit)
level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);
% figure(3)
% plot(level3Circuit)
level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);
% figure(4)
% plot(level4Circuit)
level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);
vector6=inState2.Amplitudes;
% Prep. for derivative w.r. to theta19
Z=[1 0; 0 -1];
XZ=[0 -1;1 0];
Operator=-1/2*(kron(-eye(2),XZ)+kron(Z,XZ));
temp=kron(Operator,eye(4));
vectorDer=temp*vector6;
level6Gates=[idGate(1);ryGate(2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
M6=getMatrix(level6Circuit);
vector7=M6*vectorDer;

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
M7=getMatrix(level7Circuit);
vectorOut=M7*vector7;
end

function vectorOut=Phi20(theta) %The derivative w.r. to theta20,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);
levelXGates=[cryGate(1,2,2*theta(19));idGate(3);idGate(4)];
levelXCircuit=quantumCircuit(levelXGates);
inStateX=simulate(levelXCircuit,inState6);
vectorX=inStateX.Amplitudes;
% Prep. for derivative w.r. to theta20
Z=[1 0; 0 -1];
XZ=[0 -1;1 0];
Operator=-1/2*(kron(-eye(2),XZ)+kron(Z,XZ));
temp1=kron(Operator,eye(2));
temp2=kron(eye(2),temp1);
vectorDer=temp2*vectorX;
level6Gates=[idGate(1);idGate(2);ryGate(3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
M6=getMatrix(level6Circuit);
vector7=M6*vectorDer;

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
M7=getMatrix(level7Circuit);
vectorOut=M7*vector7;
end

function vectorOut=Phi21(theta) %The derivative w.r. to theta21,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);
levelXGates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));idGate(4)];
levelXCircuit=quantumCircuit(levelXGates);
inStateX=simulate(levelXCircuit,inState6);
vectorX=inStateX.Amplitudes;
% Prep. for derivative w.r. to theta21
Z=[1 0; 0 -1];
XZ=[0 -1;1 0];
Operator=-1/2*(kron(-eye(2),XZ)+kron(Z,XZ));
temp=kron(eye(4),Operator);
vectorDer=temp*vectorX;
level6Gates=[idGate(1);idGate(2);idGate(3);ryGate(4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
M6=getMatrix(level6Circuit);
vector7=M6*vectorDer;

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
M7=getMatrix(level7Circuit);
vectorOut=M7*vector7;
end

function outState=Phi22(theta) %The derivative w.r. to theta22
% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[zGate(1);xGate(1);ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi23(theta) %The derivative w.r. to theta23,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));zGate(2);xGate(2);ryGate(2,2*theta(23));ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi24(theta) %The derivative w.r. to theta24,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));zGate(3);xGate(3);ryGate(3,2*theta(24));ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end

function outState=Phi25(theta) %The derivative w.r. to theta25,

% XY=-iY, change if I implement Ry(theta) as rotation of theta/2
level1Gates=[xGate(1);ryGate(1,2*theta(1));hGate(2);ryGate(2,2*theta(2));hGate(3);ryGate(3,2*theta(3));hGate(4);ryGate(4,2*theta(4))];
level1Circuit=quantumCircuit(level1Gates);
inState2=simulate(level1Circuit);

level2Gates=[cryGate(1,2,2*theta(5));cryGate(2,3,2*theta(6));cryGate(3,4,2*theta(7))];
level2Circuit=quantumCircuit(level2Gates);
inState3=simulate(level2Circuit,inState2);

level3Gates=[ryGate(1,2*theta(8));ryGate(2,2*theta(9));ryGate(3,2*theta(10));ryGate(4,2*theta(11))];
level3Circuit=quantumCircuit(level3Gates);
inState4=simulate(level3Circuit,inState3);

level4Gates=[cryGate(1,2,2*theta(12));cryGate(2,3,2*theta(13));cryGate(3,4,2*theta(14))];
level4Circuit=quantumCircuit(level4Gates);
inState5=simulate(level4Circuit,inState4);

level5Gates=[ryGate(1,2*theta(15));ryGate(2,2*theta(16));ryGate(3,2*theta(17));ryGate(4,2*theta(18))];
level5Circuit=quantumCircuit(level5Gates);
inState6=simulate(level5Circuit,inState5);

level6Gates=[cryGate(1,2,2*theta(19));cryGate(2,3,2*theta(20));cryGate(3,4,2*theta(21))];
level6Circuit=quantumCircuit(level6Gates);
inState7=simulate(level6Circuit,inState6);

level7Gates=[ryGate(1,2*theta(22));ryGate(2,2*theta(23));ryGate(3,2*theta(24));zGate(4);xGate(4);ryGate(4,2*theta(25))];
level7Circuit=quantumCircuit(level7Gates);
outState=simulate(level7Circuit,inState7);
end