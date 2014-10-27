%MAIN RUN: Initialize and call local ode
function MainRun(InputFilename,LoadDataFlag,handles)
global F Rmu Nspecies Temp Kw muH muOH epsilon visc DeltaP% uCST
global PrepareGridStage PrevTotalCost AChannel
global met2lit

disp('Using Spresso Without Ionic Strength');

set(handles.WhatToPlotPopup, 'UserData',1) %sets WhatToPlotPopup->UserData=1
%means that now MainRun is working and in function
%SpressoGUI->WhatToPlotPopup_Callback now dont
%call Spresso('ChangePlotVariable',handles);

tic
%profile on;

% Evaluate input file
fid=fopen(InputFilename,'r');  str=[];
while ~feof(fid)
    str=[str,fgets(fid)];
end
%read input file and evaluate the strings
eval(str);

% Additional parameters used in time stepping etc
normcontrol=0;
t0=0;
rtol=1E-3;
atol=1E-20;
dtMin=0; dtMax=1E5;
dtInit=0.00001;                  % Initial Time step [sec]
RKOrder=23;
DisableChemFlag=0;


if( strcmp(SpatialDiscFlag, 'SLIP'))
    AdaptGrid.Nconv=10;
    SmoothEps=1;
else
  AdaptGrid.Nconv=10;
  SmoothEps=1;
end

AdaptGrid.Power= 1;
uCST=0;

DeltatNextPlot=0;
PrevTotalCost=inf;
Res=[];


%--------------------------------------------------------------------------
%                                                                         -
%--------------------------------------------------------------------------

for ij=1:size(InputTable,1)
    INP{ij,1}=str2num(InputTable{ij,2});
end %for ij

% cMat contains concentrations
cMat = cMat*met2lit;                % Convert from mol/lit to mol/m^3

%Evaluate the area according to the channel shape
if ChannelShape==1
    AChannel = pi*DChannel^2/4; % Channel Area - circular channel
elseif ChannelShape==2
    AChannel = hChannel * (DChannel - 2 * hChannel) + 0.5 * pi * hChannel^2;  % Channel Area - D Shaped
end
Cur = Current/AChannel;
zetaPot = cell2mat(InputTable(:,3))'*0;    % zeta potential for each zone [mV]

% Load input file if it has been saved from the last run
if LoadDataFlag
    try
        %load([strtok(InputFilename,'.'),'.mat']);
       load([InputFilename,'at']);
    catch
        msgbox('No datafile found');
    end % try
     PrepareGridStage=0;
end % if

if AdaptGrid.Coeff==0
    PrepareGridStage=0;
end

% Calcluate pressure difference in Pa
WaterDensity = 1000; %[Kg/m^3]
g = 9.81; % Gravity [m/s^2]
%pressure is used in case we apply pressure gradient to study taylor
%dispersion etc
DeltaP = WaterDensity*g*(Pressurehead*1E-3); % Pressure difference [N/m^2]

% Resistance coefficient cannot be zero. Set DeltaP to zero to prevent
% singularity.
if bPressure==0
    bPressure=1;
    DeltaP=0;
end
%-------------------------
% SPECIES
Nspecies=size(INP,1); % get no. of species from rows in input table
%----------------------------

% Create derivative matrices
% First order
[A1.full,B1.full]=PrepareDerivMatrices(zVec,1,SpatialDiscFlag);
% Second order
[A2.full,B2.full]=PrepareDerivMatrices(zVec,2,SpatialDiscFlag);
A2.full=sparse(A2.full); B2.full=sparse(B2.full);
A1.full=sparse(A1.full); B1.full=sparse(B1.full);



if strcmp(SpatialDiscFlag,'Upwind')
    A1.UW_Left  = sparse(1/((phiVec(end) - phiVec(1))/(N-1)) * (diag([0,ones(1,N-1)],0)  + diag(-ones(1,N-1),-1)));
    A1.UW_Right = sparse(1/((phiVec(end) - phiVec(1))/(N-1)) * (diag([-ones(1,N-1),0],0) + diag(ones(1,N-1),1)));
end


% Runge Kutta Dormand???Prince method (RK45)
if RKOrder==45
    RK_pow = 1/5;
    RK_A = [1/5, 3/10, 4/5, 8/9, 1, 1];
    RK_B = [
        1/5         3/40    44/45   19372/6561      9017/3168       35/384
        0           9/40    -56/15  -25360/2187     -355/33         0
        0           0       32/9    64448/6561      46732/5247      500/1113
        0           0       0       -212/729        49/176          125/192
        0           0       0       0               -5103/18656     -2187/6784
        0           0       0       0               0               11/84
        0           0       0       0               0               0
        ];
    RK_E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];
elseif RKOrder==23 % Runge Kutta Bogacki-Shampine method (RK23)
    RK_pow = 1/3;
    RK_A = [1/2, 3/4, 1];
    RK_B = [
        1/2         0               2/9
        0           3/4             1/3
        0           0               4/9
        0           0               0
        ];
    RK_E = [-5/72; 1/12; 1/9; -1/8];
end

% Smoothing function parameters
M.full=(A2.full-SmoothEps*B2.full);
M.full(1,:)=zeros(1,N); M.full(end,:)=zeros(1,N); M.full(1,1)=1; M.full(end,end)=1;
M.full=sparse(M.full);

%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

A1.a=diag(A1.full,0); A1.b=diag(A1.full,-1); A1.c=diag(A1.full,+1);
A2.a=diag(A2.full,0); A2.b=diag(A2.full,-1); A2.c=diag(A2.full,+1);
M.a=diag(M.full,0);   M.b=diag(M.full,-1);   M.c=diag(M.full,+1);

%--------------------------------------------------------------------------
%  MARCH IN TIME
%--------------------------------------------------------------------------

if( strcmp(SpatialDiscFlag, 'SLIP'))

    dz=zVec(2)-zVec(1);
    cSize=size(cMat);
    ARGS={cSize,DeltatNextPlot,DeltaCounterNextPlot,dtInit,dtMin,dtMax,rtol, ...
    atol,normcontrol,RK_A,RK_B,RK_E,RK_pow,uCST,A1,B1,A2,B2,M,Cur,INP, ...
    AdaptGrid,RKOrder,SteadyStateFlag,L1,L2,zetaPot,hChannel,bPressure,DeltaP,...
    betaDispersion,Res,DisableChemFlag,SpatialDiscFlag,handles, dz};

    %hard code LocalODE45FV for finite volume scheme with conservation variable J*A*c
    disp('Running Finite Volume SLIP Scheme');
    AreaRatioCellVec = AreaRatio; %getAreaRatio(phiVec);
    AreaRatioCellMat = AreaRatioCellVec*ones(1,Nspecies);

OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*phiVec,N);
dphidzVec=OUT(:,1); dphidzMat=dphidzVec*ones(Nspecies,1)'; %J=dx/dz


    qMat=cMat.*AreaRatioCellMat.*dphidzMat; %here J=dx/dz=1
    cSize=size(cMat);

    IC=[reshape(qMat,numel(cMat),1);phiVec; AreaRatioCellVec.*dphidzVec];

    LocalOde45FV(IC,[t0 tEnd],ARGS);%pass initial condition in IC in y0
else
    cSize=size(cMat);
    ARGS={cSize,DeltatNextPlot,DeltaCounterNextPlot,dtInit,dtMin,dtMax,rtol, ...
    atol,normcontrol,RK_A,RK_B,RK_E,RK_pow,uCST,A1,B1,A2,B2,M,Cur,INP, ...
    AdaptGrid,RKOrder,SteadyStateFlag,L1,L2,zetaPot,hChannel,bPressure,DeltaP,...
    betaDispersion,Res,DisableChemFlag,SpatialDiscFlag,handles};

    IC=[reshape(cMat,numel(cMat),1);phiVec];

    if (max(abs(AreaRatio-1))>0)
    fprintf('Using constant channel cross-section\n');
    fprintf('Variable cross-section is avialable only with SLIP scheme\n');
    end

    LocalOde45(IC,[t0 tEnd],ARGS);%pass initial condition in IC in y0
end


set(handles.WhatToPlotPopup, 'UserData',0) %sets WhatToPlotPopup->UserData=1
%profile viewer;

toc

% Time stepping based on ode45
%probably make another function LocalOde23, and remove all if conditions
%for RK45 or RK23



%-----------------------------
