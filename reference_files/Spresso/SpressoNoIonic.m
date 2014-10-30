function SpressoNoIonic(varargin)
global F Rmu Nspecies Temp Kw muH muOH epsilon visc DeltaP uCST Ngrid
global PrepareGridStage PrevTotalCost
global met2lit
global ColorList
global EquiPLength PolDeg Ngrid%EquilibriumPLength and Ploynomial Degree

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',1,...
    'defaultlinelinewidth',1,'defaultpatchlinewidth',2,'defaultaxesfontname','Times');

ColorList=([0.1 0.5 0.5; 1 0 0; 0 0.5 0;  0 0.75 0.75; 0.75 0.75 0; ...
    0.6 0.2 0; 0.25 0.25 0.25; 0.08 0.17 0.55; 0.4 0.52 0.78; ...
    0.68 0.47 0; 0.8 0.7 1; 0.7 0.78 1; 0 1 1; 1 1 0; 1 0.6 0.9; 0.1 0.7 0.1; 0.5 1 1; 0 0.7 0.1]);

% Constants
F=9.65E4;       % Faraday's const.[C/mol]
Rmu=8.31;       % Universal gas const. [J/mol*K]
Temp=298;       % Temperature [K]
Kw=1E-14;       % Water equilibrium constant
muH=362E-9/F;           % Mobility of Hydronium   % [m^2/s*V]/F --> [mol*s/Kg]
muOH=205E-9/F;          % Mobility of Hydroxide   % [m^2/s*V]/F --> [mol*s/Kg]
met2lit=1000;           % m^3 to liters
visc=1E-3;              % Dynamic viscosity (water) [Pa s]

% If no input arguments - run graphical user interface

if isempty(varargin)
    clc;
    warning off all
    handles = SpressoGUI;
    return;
elseif nargin && ischar(varargin{1})
    feval(str2func(varargin{1}),varargin{2:end});
end


%--------------------------------------------------------------------------
% Calculates the Flux For Finite Volume SLIP Scheme With Area Change
% This function is called in the LocalOde45 function to compute flux 
% at each time step of integration
% at Start of every RK time steps CalcChemEqFlag = 1
% for intermediate RK time steps CalcChemEqFlag = 0
%--------------------------------------------------------------------------

%Area Change is in terms of Area Ratio based on area at leftmost grid point.
%For area change finite volume, the conservation variable is J*c*A

function  [dydt,cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec]= ...
    CalcFluxFVSlip(y0,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz)
global PlotTime
global PrepareGridStage PrevTotalCost Rmu Temp F met2lit
global muH muOH Kw epsilon visc hChannel bPressure DeltaP betaDispersion
global uCST Ngrid PolDeg AChannel

[A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,PolDeg,...
    EquiPLength,LMat,LCube,muCube,ValCube,DCube,N,Nspecies,...
    zetaPot,h,DisableChemFlag,SpatialDiscFlag]=ARGS{:};
cSize=[N,Nspecies];
qMat=reshape(y0(1:cSize(1)*cSize(2)),cSize); %N x Nspecies, q=J*c*A; 
phiVec=y0(cSize(1)*cSize(2)+1:cSize(1)*cSize(2)+N); %x(z) spatial coordinate based on grid adaptation

OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*phiVec,N);
dphidzVec=OUT(:,1); dphidzMat=dphidzVec*ones(Nspecies,1)'; %J=dx/dz 

JAreaRatioCellVec=y0(cSize(1)*cSize(2)+N+1:cSize(1)*cSize(2)+2*N);
AreaRatioCellVec=JAreaRatioCellVec./dphidzVec;
AreaRatioCellMat=AreaRatioCellVec*ones(1,Nspecies); %area ratio at cell centers

OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*AreaRatioCellVec,N);
dAdzVec=OUT(:,1); dAdzMat=dAdzVec*ones(Nspecies,1)';

cMat=qMat./(AreaRatioCellMat.*dphidzMat);

%-----------------------------------------
% Arificial Diffusion for Stabilitizing (size N-1  x Nspecies)
%------------------------------------------
LimiterFlag='ELED';        %'minmod' / 'superbee' / 'VanLeer' /'none'/'ELED'

%cMat=qMat./dphidzMat;

% Derivatives requried for adaptive grid  (but also for next steps)
%gives out 1st derivatives of cH', cMat
%A1.a,b,c are diagonal terms, in matrices Af'=Bf

% OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[cH',cMat],N);
% dcHdz=OUT(:,1); dcdzMat=OUT(:,2:Nspecies+1);   

% 
% % OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[cH',cMat, AreaRatioCellVec.*cH', AreaRatioCellMat.*cMat],N);
% % dcHdz=OUT(:,1); dcdzMat=OUT(:,2:Nspecies+1); dAcHdz=OUT(:,Nspecies+2); dAcdzMat=OUT(:,Nspecies+3:end);
% % 
% % OUT=LinSolve(A2.a,A2.b,A2.c,B2.full*[cMat],N);
% % d2Dcdz2Mat=OUT(:,1:Nspecies);

%--------------------
% Solve Equilibrium
%DisableChemFlag is only for testing some generic cases It is not actually
%inputed from the GUI
if DisableChemFlag 

    muMat=F*sum(ValCube.*muCube,3)';
    DMat=DCube(:,:,1)';
    alphaMat=F^2*sum(ValCube.^2.*muCube,3)';
    betaMat=sum(F*ValCube.*DCube,3)';

    %ValTmp=sum(ValCube,3)';
  %  cMat(:,end)=-sum(ValTmp(:,1:end-1).*cMat(:,1:end-1),2)./ValTmp(:,end);
else
    if CalcChemEqFlag
        [cizCube,cH,cHCubePower,gizCube]=LzCalcEquilibrium(cH,LCube,cMat, ValCube);
        muMat=F*sum(ValCube.*muCube.*gizCube,3)';
        DMat=sum(DCube.*gizCube,3)';
        alphaMat=F^2*sum(ValCube.^2.*muCube.*gizCube,3)';
        betaMat=F*sum(ValCube.*DCube.*gizCube,3)';
    end
end % end of if CalcChemEqFlag

    % Calculate velocity (pressure driven)
    L=phiVec(end)-phiVec(1);
    uBulk = hChannel^2/(bPressure*visc*L)*(-DeltaP); %uBulk is net flow rate 
						     %constant along axis (Q)
  
    
if CalcChemEqFlag
    % Dispersion
    PeMat=uBulk*(hChannel)./DMat;                   %h=diameter  h/2=radius
    DMat=DMat.*(1+betaDispersion*PeMat.^2);
end % if CalcChemEqFlag

if DisableChemFlag
    SigVec = sum(alphaMat'.*cMat',1)';
    SVec   = sum(betaMat'.*cMat',1)';
else
    SigVec = sum(alphaMat'.*cMat',1)' + F^2*(muH*cH*met2lit+muOH*(Kw./cH)*met2lit)';
    SVec   = sum(betaMat'.*cMat',1)'  + Rmu*Temp*(muH*cH'*met2lit-muOH*(Kw./cH)'*met2lit);
end


ConstVoltFlag=0;
if (ConstVoltFlag) %for constant voltage case
    %calculate the current
    Voltage = -350; 
    NumCur=2*sum(diff(SVec)./(SigVec(1:N-1)+ SigVec(2:N) ),1);
    DenCur=trapz(phiVec,1./(SigVec.*AreaRatioCellVec));
    Cur=-(Voltage+NumCur)./DenCur;
%    fprintf('Current=%g\n',Cur);
end

if SteadyStateFlag %&& FirstTimeFlag   % Automatic velocity, based on LE (first chemical on list)
   uCST = -Cur/SigVec(end)*muMat(end,1)-uBulk;    
%    uCST = -Cur/SigVec(1)*muMat(1,1)-uBulk;  
end

%FV Scheme
%|--------------------|-------------------|
%|		      |			  |
%|	  j	      |		j+1	  |
%|		      |			  |	
%|--------------------|-------------------|
%                   j+1/2
% f_{j+1/2}=(f_j+f_{j+1})/2 + aplha_{j+1/2}(diffusion)+D_{j+1/2}(c_j+1-c_j)/dz
%  dju/dt+f_{j+1/2}-f_{j-1/2}=D_{j+1/2}(c_j+1-c_j)/dz-D_{j+1/2}(c_j+1-c_j)/dz

%implementation
%cells 1 and N are "ghost cells": no diffusion, apply BC using incoming/outgoing waves 
%define flux on 1 and N neglecting diffusion terms in electric field etc.
%for cells 2:N-1, solve using SLIP scheme including diffusion terms

%numerical diffusivity without diffusion constant. use this to adapt grid
cMatExtended=[cMat(1,:); cMat;cMat(end,:)];  % length: N+2
dcMat=diff(cMatExtended,[],1);              % length: N+1
LimitMat=LimiterFunc(dcMat(3:end,:),dcMat(1:end-2,:),LimiterFlag, dz); %(N-1)
CorrectedDiffEdgeMat = dcMat(2:end-1,:)-LimitMat;  % Corrected numerical diffusivity size (N-1) x Nspecies

if AdaptGrid.Coeff~=0
    mu_c_Char=max(max(abs(muMat.*cMat)));  % characteristic mobility*concentration
    L=phiVec(end)-phiVec(1); %length of domain
%     AdaptGrid.wtConst=((N-AdaptGrid.PointsAtInterface)/(AdaptGrid.PointsAtInterface))* ...
%         Rmu*Temp*mu_c_Char/(abs(Cur)*L);
    
    AdaptGrid.wtConst=(N/(AdaptGrid.PointsAtInterface))*Rmu*Temp*mu_c_Char/(abs(Cur)*L);
       
% %     
% % % % %----------------------------
% OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[cH',cMat],N);
%  dcHdz=OUT(:,1); dcdzMat=OUT(:,2:Nspecies+1); 
%    wttt=abs([dcdzMat,dcHdz]);
%     wtt=wttt./(ones(N,1)*max([cMat,cH'],[],1));
%     wt_gradient=max(wtt,[],2).^(AdaptGrid.Power);
 
   %weight= wt_gradient;
  %%%  fprintf('Max Weight = %g, Max J = %g\n', max(weight), max(dphidzVec) );
% % % % %---------------------------   

%%%%grid adaptation based on numerical dissipation and area variation

    wttt=abs(dAdzMat);
    wtt=wttt./AreaRatioCellMat; %  (ones(N,1)*max([AreaRatioCellMat.*cMat,AreaRatioCellMat(:,1).*cH'],[],1));
    wt=wtt/(max(wtt)+1);
    
  wttt_ndif=[0; max(abs(CorrectedDiffEdgeMat(2:N-1,:)-CorrectedDiffEdgeMat(1:N-2,:)),[],2); 0];
  wtt_ndif=wttt_ndif/max(wttt_ndif);
  wt_ndif =wtt_ndif.^(AdaptGrid.Power);
  

weight= 0.02*wt + wt_ndif/(max(wt_ndif)) ; 
cost = dphidzVec + 0.75*weight/(sum(weight)*AdaptGrid.wtConst/N); %Here sum w/N is the characteristic vaule of w
%fprintf('Max wt = %g, Max Weight = %g, Max J = %g, Max Cost =%g\n', max(wt), max(weight), max(dphidzVec), max(cost));

    if AdaptGrid.Nconv~=0
        cost=conv(cost,[1:AdaptGrid.Nconv,AdaptGrid.Nconv:-1:1]/AdaptGrid.Nconv^2);
        cost=cost(AdaptGrid.Nconv:end-AdaptGrid.Nconv);
    end
    costTakeDeriv=cost.*dphidzVec;
    costDeriv=LinSolve(A1.a,A1.b,A1.c,B1.full*costTakeDeriv,N);

    gradVec=costDeriv;
    RHS=(A2.full(2:end-1,2:end-1)*gradVec(2:end-1));
    gradSmoothVec=[0;LinSolve(M.a(2:end-1),M.b(2:end-1),M.c(2:end-1),RHS,N-2);0];

    dphidtVec=0.002*AdaptGrid.Coeff/AdaptGrid.wtConst*gradSmoothVec;
    dphidtVec(1)=0; dphidtVec(end)=0;
else
    % PrepareGridStage=0;
    dphidtVec=phiVec*0;
    cost=dphidtVec+1;
    dphidtMat=zeros(N,Nspecies);
end
dphidtMat=dphidtVec*ones(1,Nspecies); %size N x Nspecies

% If in prepare grid stage - calculate derivaties and return
% No need to solve equations
if PrepareGridStage
        
AreaRatioVecExtended=[AreaRatioCellVec(1); AreaRatioCellVec ; AreaRatioCellVec(end)];  % length: N+2
dARatioVec=diff(AreaRatioVecExtended);              % length: N+1
ALimitVec=LimiterFunc(dARatioVec(3:end,1), dARatioVec(1:end-2,1),LimiterFlag, dz); %(N-1)
ACorrectedDiffEdgeVec = dARatioVec(2:end-1,:)-ALimitVec;  % Corrected numerical diffusivity (size (N-1) x Nspecies
AalpVec=0.26*abs(dphidtVec(2:N)+dphidtVec(1:N-1));

AGridFluxCellVec = -AreaRatioCellVec.*dphidtVec;
AGridFluxEdgeVec = (AGridFluxCellVec(1:N-1)  + AGridFluxCellVec(2:N))/2 - AalpVec.*ACorrectedDiffEdgeVec;
dARatiodtVec(2:N-1,1) = -diff(AGridFluxEdgeVec)./dz;
dARatiodtVec(1,1)= -AGridFluxCellVec(2,:)/dz;
dARatiodtVec(N,1)= AGridFluxCellVec(N-1)/dz;

alpTempEdgeVec = 0.5*(-dphidtVec(1:N-1,1).*AreaRatioCellVec(1:N-1,1) - dphidtVec(2:N,1).*AreaRatioCellVec(2:N,1)); 
alpVecEdge = 0.52*abs(alpTempEdgeVec);
alpMatEdge=alpVecEdge*ones(1,Nspecies);

CorrectedAdFluxEdgeMat = - alpMatEdge.*CorrectedDiffEdgeMat;
GridFluxCellMat = -AreaRatioCellMat.*dphidtMat.*cMat;
    GridFluxEdgeMat = (GridFluxCellMat(1:N-1,:)  + GridFluxCellMat(2:N,:))/2- ...
        (cMat(1:N-1,:) + cMat(2:N,:))/2.*(AalpVec.*ACorrectedDiffEdgeVec*ones(1,Nspecies));

dqdtMat(2:N-1,:) =  diff(-GridFluxEdgeMat -CorrectedAdFluxEdgeMat)./(dz);
dqdtMat(1,:)= -GridFluxCellMat(2,:)/dz;
dqdtMat(N,:)= GridFluxCellMat(N-1,:)/dz;    

dydt=[reshape(dqdtMat,numel(dqdtMat),1);dphidtVec; dARatiodtVec];
    %tStepMat=DMat*0;
    dPotendxMat=DMat*0;
    SigVec=dphidtVec*0;
    return;
end

%Grid Flux = flux due to moving grid -dphidt*c
% GridFluxCellMat = -dphidtMat.*cMat;
% GridFluxEdgeMat = (GridFluxCellMat(1:N-1,:)  + GridFluxCellMat(2:N,:))/2; 

%-----------------------
% calculate advective flux on cell centers (size N x Nspecies)
% advective flux = Q*c  + mu_i*ci*j/sigvec
% defined on j=1:N
%-----------------------

SigMat=SigVec*ones(1, Nspecies) ;%matrix of conductivity size N x Nspecies 
AdFluxCellMat=(uBulk+uCST)*cMat + muMat.*cMat*Cur./SigMat; %Cur is the total current

%---------------------------------------
% calculate Electrodiffusion Flux
% Electrodiffusion flux = mu_i*c_i*A/sigma*(sum z_i*D_i*F*dc/dx)
% Electrodiffusion flux = mu_i*c_i*A/sigma*(d/dx SVec) = 1/J*mu_i*c_i*A/sigma*(d/dz SVec)
% Defined for (N-1) cell edges (j+1/2)
%---------------------------------------
SMat=SVec*ones(1, Nspecies); %matrix of S size N x Nspecies 
ElecDiffFactorCellMat=muMat.*cMat./SigMat;

%Case 1: D_eff = (D(i+1)*A^2(i+1) + D(i)*A^2(i))/(A(i+1)*dxdz(i+1)+A(i)*dxdz(i))
ElecDiffFactorEdgeMat=(ElecDiffFactorCellMat(1:N-1,:).*AreaRatioCellMat(1:N-1,:)./dphidzMat(1:N-1,:) + ...
    ElecDiffFactorCellMat(2:N,:).*AreaRatioCellMat(2:N,:)./dphidzMat(2:N,:))/2;
   
% %Case 2: D_eff = (D(i+1)*A^2(i+1) + D(i)*A^2(i))/(A(i+1)*dxdz(i+1)+A(i)*dxdz(i))
% ElecDiffFactorEdgeMat=(ElecDiffFactorCellMat(1:N-1,:).*AreaRatioCellMat(1:N-1,:).^2 + ElecDiffFactorCellMat(2:N,:).*AreaRatioCellMat(2:N,:).^2)...
%         ./(AreaRatioCellMat(1:N-1,:).*dphidzMat(1:N-1,:) + AreaRatioCellMat(2:N,:).*dphidzMat(2:N,:));

ElecDiffFluxEdgeMat=ElecDiffFactorEdgeMat.*diff(SMat)/dz; %size (N-1) x Nspecies

%---------------------------------------
% calculate Molecular Diffusion Flux
% MolecDiffusion flux = 1/J*D*A*(dc/dx)
% MolecDiffusion flux = mu_i*c_i*A/sigma*(d/dx SVec) = 1/J*mu_i*c_i*A/sigma*(d/dz SVec)
% Defined for (N-1) cell edges (j+1/2)
%---------------------------------------

MolecDiffFactorEdgeMat=(AreaRatioCellMat(1:N-1,:)./dphidzMat(1:N-1,:) + ...
                            AreaRatioCellMat(2:N,:)./dphidzMat(2:N,:))/2; 
MolecDiffFluxEdgeMat = MolecDiffFactorEdgeMat.*diff(DMat.*cMat,1)./dz; %size (N-1) x Nspecies
                
% Extend cMat to get numerical diffusivity @ boundaries
% %--------------------------------
% cMatExtended=[cMat(1,:); cMat;cMat(end,:)];  % length: N+2
% dcMat=diff(cMatExtended,[],1);              % length: N+1
% LimitMat=LimiterFunc(dcMat(3:end,:),dcMat(1:end-2,:),LimiterFlag, dz); %(N-1)
% CorrectedDiffEdgeMat = dcMat(2:end-1,:)-LimitMat;  % Corrected numerical diffusivity (size (N-1) x Nspecies
% %-------------

alpTempEdge = (uBulk+uCST) +  0.5*(muMat(1:N-1,:).*Cur./SigMat(1:N-1,:)+muMat(2:N,:).*Cur./SigMat(2:N,:) -...
                        dphidtMat(1:N-1,:).*AreaRatioCellMat(1:N-1,:) - dphidtMat(2:N,:).*AreaRatioCellMat(2:N,:)); 


alpVecEdge= 0.52*(max(abs(alpTempEdge),[],2));
alpMatEdge=alpVecEdge*ones(1,Nspecies);
CorrectedAdFluxEdgeMat =  (AdFluxCellMat(1:N-1, : ) + AdFluxCellMat(2:N, : ))/2 ...
  								- alpMatEdge.*CorrectedDiffEdgeMat;

%---------------------------------
%dA/dt = 0

AreaRatioVecExtended=[AreaRatioCellVec(1); AreaRatioCellVec ; AreaRatioCellVec(end)];  % length: N+2
dARatioVec=diff(AreaRatioVecExtended);              % length: N+1
ALimitVec=LimiterFunc(dARatioVec(3:end,1), dARatioVec(1:end-2,1),LimiterFlag, dz); %(N-1)
ACorrectedDiffEdgeVec = dARatioVec(2:end-1,:)-ALimitVec;  % Corrected numerical diffusivity (size (N-1) x Nspecies
AalpVec=0.26*abs(dphidtVec(2:N)+dphidtVec(1:N-1));
AGridFluxCellVec = -AreaRatioCellVec.*dphidtVec;
AGridFluxEdgeVec = (AGridFluxCellVec(1:N-1)  + AGridFluxCellVec(2:N))/2 - AalpVec.*ACorrectedDiffEdgeVec;

dARatiodtVec(2:N-1,1) = -diff(AGridFluxEdgeVec)./dz;
dARatiodtVec(1,1)= -AGridFluxCellVec(2,:)/dz;
dARatiodtVec(N,1)= AGridFluxCellVec(N-1)/dz;

%----------------------------------  

    GridFluxCellMat = -AreaRatioCellMat.*dphidtMat.*cMat;
    GridFluxEdgeMat = (GridFluxCellMat(1:N-1,:)  + GridFluxCellMat(2:N,:))/2- ...
        (cMat(1:N-1,:) + cMat(2:N,:))/2.*(AalpVec.*ACorrectedDiffEdgeVec*ones(1,Nspecies));
    
dqdtMat(2:N-1,:) =  diff(-GridFluxEdgeMat -CorrectedAdFluxEdgeMat +  MolecDiffFluxEdgeMat-ElecDiffFluxEdgeMat)./(dz);

%------------ zero order extrapolation at boundaries
%%Assume C_0 = C_1 and C_(N+1) = C_N
dqdtMat(1,:) =  (-GridFluxCellMat(2,:) -CorrectedAdFluxEdgeMat(1,:) + AdFluxCellMat(1,: ) ...
    + MolecDiffFluxEdgeMat(1,:)-ElecDiffFluxEdgeMat(1,:))/dz;

dqdtMat(N,:) =  (GridFluxCellMat(N-1,:) + CorrectedAdFluxEdgeMat(N-1,:) - AdFluxCellMat(N, : ) ...
    - MolecDiffFluxEdgeMat(N-1,:) + ElecDiffFluxEdgeMat(N-1,:))/dz;
%------------

% % %Apply boundary conditions-----------%works
% %    %     Left boundary
%     IL=1;
%     V1L=uBulk + uCST + muMat(IL,:).*Cur./SigVec(IL); %- alpVecEdge(1); %Q + mu*I/sigma
%     V2L=-Cur./SigVec(IL)^2.*(cMat(IL,:).*muMat(IL,:));
%     V3L=alphaMat(IL,:);
% 
%     AL=diag(V1L)+V2L'*V3L;      [VAL,DAL] = eig(AL); %AL=VAL*DAL*VAL^{-1}
%     SEigAL=real(DAL);
%     DAL_negative=(SEigAL<0).*DAL; %diag matrix with only -ve eigenvalues
%         
%     %for characteristics leaving the domain, compute derivatives
%     %for characteristics enterting the domain set derivatives = 0
%     
%     %dq/dt= V*Lambda^{-}*V^{-1}*dc/dx + dphi/dt*dc/dx
%     dcdzL=(cMat(2,:)-cMat(1,:))/dz; 
%     %FluxLeftBC =VAL*(DAL_negative)*(VAL\dcdzL');
%     FluxLeftBC =VAL*(DAL_negative)*(VAL\dcdzL');
%     dqdtMat(1,:)= -GridFluxCellMat(2,:)/dz - FluxLeftBC' + alpVecEdge(1)*dcdzL;
   
    
    %     Right boundary
    IR=N;
    V1R = uBulk + uCST + muMat(IR,:).*Cur./SigVec(IR); %Q + mu*I/sigma
    V2R=-Cur./SigVec(IR)^2.*(cMat(IR,:).*muMat(IR,:));
    V3R=alphaMat(IR,:);

    AR=diag(V1R)+V2R'*V3R;  [VAR,DAR] = eig(AR); %AL=VAL*DAL*VAL^{-1}
    SEigAR=real(DAR);
    DAR_positive=(SEigAR>0).*DAR;  
    %for characteristics leaving the domain, compute derivatives
    %for characteristics enterting the domain set derivatives = 0
    %dq/dt= V*Lambda^{-}*V^{-1}*dc/dx + dphi/dt*dc/dx
    
    dcdzR=(cMat(N,:)-cMat(N-1,:))/dz; 
    FluxRightBC =VAR*(DAR_positive/VAR)*dcdzR';
    dqdtMat(N,:)= GridFluxCellMat(N-1,:)/dz - FluxRightBC' - alpVecEdge(N-1)*dcdzR;
% % %------------------ end of boundary conditions------


% % % %Apply boundary conditions-----------%works
% % %    %     Left boundary
%     IL=1;
%     V1L=uBulk + uCST + muMat(IL,:).*Cur./SigVec(IL); %- alpVecEdge(1); %Q + mu*I/sigma
%     V2L=-Cur./SigVec(IL)^2.*(cMat(IL,:).*muMat(IL,:));
%     V3L=alphaMat(IL,:);
% 
%     AL=diag(V1L)+V2L'*V3L;      [VAL,DAL] = eig(AL); %AL=VAL*DAL*VAL^{-1}
%     SEigAL=real(DAL);
%     DAL_negative=(SEigAL<0).*DAL; %diag matrix with only -ve eigenvalues
%         
%     %for characteristics leaving the domain, compute derivatives
%     %for characteristics enterting the domain set derivatives = 0
%     
%     %dq/dt= V*Lambda^{-}*V^{-1}*dc/dx + dphi/dt*dc/dx
%     dcdzL=(cMat(2,:)-cMat(1,:))/dz; 
%     %FluxLeftBC =VAL*(DAL_negative)*(VAL\dcdzL');
%     FluxLeftBC =VAL*(DAL_negative)*(VAL\dcdzL');
%     dqdtMat(1,:)= -GridFluxCellMat(2,:)/dz - FluxLeftBC' + alpVecEdge(1)*dcdzL;
%    
%     
%     %     Right boundary
%     IR=N;
%     V1R = uBulk + uCST + muMat(IR,:).*Cur./SigVec(IR); %Q + mu*I/sigma
%     V2R=-Cur./SigVec(IR)^2.*(cMat(IR,:).*muMat(IR,:));
%     V3R=alphaMat(IR,:);
% 
%     AR=diag(V1R)+V2R'*V3R;  [VAR,DAR] = eig(AR); %AL=VAL*DAL*VAL^{-1}
%     SEigAR=real(DAR);
%     DAR_positive=(SEigAR>0).*DAR;  
%     %for characteristics leaving the domain, compute derivatives
%     %for characteristics enterting the domain set derivatives = 0
%     %dq/dt= V*Lambda^{-}*V^{-1}*dc/dx + dphi/dt*dc/dx
%     
%     dcdzR=(cMat(N,:)-cMat(N-1,:))/dz; 
%     FluxRightBC =VAR*(DAR_positive/VAR)*dcdzR';
%     dqdtMat(N,:)= GridFluxCellMat(N-1,:)/dz - FluxRightBC' - alpVecEdge(N-1)*dcdzR;
% % % %------------------ end of boundary conditions------
    
OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[SVec],N);
dSdzMat=OUT(:,1)*ones(1, Nspecies);
dPotendxMat=-(Cur./AreaRatioCellMat +  dSdzMat./dphidzMat)./SigMat; %Cur is current density of ARatio=1

dydt=[reshape(dqdtMat,numel(dqdtMat),1);dphidtVec; dARatiodtVec];

%-------------------------------------------------

%ode45 time stepping for Finite Volume
%the conservation variable is q=J*c*A, where J=dx/dz
function LocalOde45FV(y0,tspan,ARGS)
global Nspecies betaDispersion bPressure hChannel Ngrid
global PlotTime errVec PrepareGridStage AChannel
global uCST
[cSize,DeltatNextPlot,DeltaCounterNextPlot,dtInit,dtMin,dtMax,rtol,atol,...
    normcontrol,RK_A,RK_B,RK_E,RK_pow,uCST,A1,B1,A2,B2,M,Cur,INP,...
    AdaptGrid,RKOrder,SteadyStateFlag,L1,L2,zetaPot,hChannel,bPressure, ...
    DeltaP,betaDispersion,Res,DisableChemFlag,SpatialDiscFlag,handles, dz]=ARGS{:};
N=cSize(1);

% Additional parameters
qMat=reshape(y0(1:cSize(1)*cSize(2)),cSize);
%phi refers to the physical domain
phiVec=y0(cSize(1)*cSize(2)+1:cSize(1)*cSize(2)+N);
JAreaRatioCellVec=y0(cSize(1)*cSize(2)+N+1:end); 

% Additional parameters
t0=tspan(1); tfinal=tspan(end);
tdir = sign(tfinal - t0);
threshold = atol / rtol;
htspan=dtInit;                  % Initial time step
hmin=dtMin;                 % Minimum time step
hmax=dtMax;
neq=length(y0);
normy=norm(y0);
tNextPlot=0;
Nspecies=cSize(2);

t=t0;  y=y0;
FirstTimeFlag=1; %Flag to indicate the first time step
ODEcounter=0;
CounterNextPlot=1;
PlotCounter=0;
cMatAllTimes=[]; phiVecAllTimes=[]; tVecOut=[];

JAreaRatioCellMat=JAreaRatioCellVec*ones(1,Nspecies); %area ratio at cell centers
cMat = qMat./JAreaRatioCellMat; %initially dx/dz=1

% Create Polynomials for Equilibrium reactions
[PMat,PPrimeMat,Q,QPrime,LMat,ValMat,muMat,KaMat,DMat]=EquilibriumPolynomials(INP);
EquiPLength=max(size(PMat,2),size(Q,2));
PolDeg=size(LMat,2); % Polynomial degree

%Cube here refers to 3-d array
%cols: grid points, rows: no. of species depth: polynomial degree
muCube=repmat(reshape(muMat,[Nspecies,1,PolDeg]),[1,N,1]); %mobilities
DCube=repmat(reshape(DMat,[Nspecies,1,PolDeg]),[1,N,1]); %diffusivities
ValCube=repmat(reshape(ValMat,[Nspecies,1,PolDeg]),[1,N,1]); %valence
LCube=FastRepmatColumns(reshape(LMat,[Nspecies,1,PolDeg]),N); %???
cH=zeros(1,size(qMat,1)); %set initial [H+] = 0 
Ngrid=length(cH);

% THE MAIN LOOP
done = false;
%time stepping starts here
while ~done
    ODEcounter=ODEcounter+1;


    if FirstTimeFlag==1 %initialy it is set = 1 for first run
        CalcChemEqFlag=1; %initially evaluate chemical equilibrium
        errVec=ones(N,1); [muMat,DMat,alphaMat,betaMat]=deal([]);
        h=ones(1,N*(Nspecies+1));  %Just for initial flux calculation
        GridCost=zeros(N,1);
        
        %run CalculateEquilibrium when FirstTimeFlag=1
        %otherwise run LzCalcEquilibrium inside CalcFlux
        [cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
    cH,PMat,Q,PPrimeMat,QPrime,PolDeg,EquiPLength,N,LMat,LCube,cMat,Nspecies);
        
        %Define arguments to call the CalcFlux function
        ARGS={A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag, PolDeg,EquiPLength,...
            LMat,LCube,muCube,ValCube,...
            DCube,N,Nspecies,zetaPot,h,DisableChemFlag,SpatialDiscFlag};
        
        %call CalcFlux function
        [f0,cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec]= ...
            CalcFluxFVSlip(y0,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz);
     
        

        % Initialize method parameters.
        if RKOrder==45,       f = zeros(neq,7);
        elseif RKOrder==23,   f = zeros(neq,4);
        end
        f(:,1) = f0;

        % Compute an initial step size h using y'(t).
        absh = min(hmax, htspan);
        rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^RK_pow);
        if absh * rh > 1
            absh = 1 / rh;
        end
        absh = max(absh, hmin);
        absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
        h = tdir * absh;

        FirstTimeFlag=0; %1st step apready done so set this flag = 0
        CalcChemEqFlag=0; %I dont want to compute chemical equilibrim for 
                          %intermediate steps
                          
    end %end of the First time's if condition
    %could probably move this out of the loop.. 


    h = min(hmax, max(hmin, h));
    if h>(tfinal-t)
        h=tfinal-t;
    end

    
    %NOTE: Chemical equilibrium is not computed in intermediate steps
    % LOOP FOR ADVANCING ONE STEP.
    nofailed = true;                      % no failed attempts
    while true
        
        ARGS={A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,PolDeg,EquiPLength,...
            LMat,LCube,muCube,ValCube,...
            DCube,N,Nspecies,zetaPot,h,DisableChemFlag,SpatialDiscFlag};
        
        [f(:,2),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec]= ...
            CalcFluxFVSlip(y+h.*(f*RK_B(:,1)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz);
        
        [f(:,3),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec]=...
            CalcFluxFVSlip(y+h.*(f*RK_B(:,2)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz);

        if RKOrder==45
            [f(:,4),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec]=...
                CalcFluxFVSlip(y+h.*(f*RK_B(:,3)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz);
            
            [f(:,5),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec]=...
                CalcFluxFVSlip(y+h.*(f*RK_B(:,4)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz);
            
            [f(:,6),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec]=...
                CalcFluxFVSlip(y+h.*(f*RK_B(:,5)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz);
        end

        if     RKOrder==45,  tnew = t + min(h*RK_A(6));
        elseif RKOrder==23,  tnew = t + min(h*RK_A(3));
        end

        if RKOrder==45,
            deltay=h.*(f*RK_B(:,6));
            ynew = y + deltay;
            [f(:,7),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec] = ...
                CalcFluxFVSlip(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz);
            
        elseif RKOrder==23,
            deltay=h.*(f*RK_B(:,3));
            ynew = y + deltay;
            [f(:,4),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec] = ...
                CalcFluxFVSlip(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz);
        end

        % Estimate the error.
            errAll= h.*abs((f * RK_E) ./ max(max(abs(y),abs(ynew)),threshold));
            errMat=reshape(errAll,[N Nspecies+2]);
            errVec=max(errMat,[],2);           
            err=norm(errVec);            
            
        if err < rtol                       % Successful step - go to next time step
            break;
            
        else                                % Failed step - calculate new time step
            if nofailed
                nofailed = false;
                h = max(hmin, h.*max(0.1, 0.8*(rtol./err).^RK_pow));
            else
                h = max(hmin, 0.5 * h);
            end
            done = false;
        end % if err > rtol
    end % while true. End of while loop over the time steps after 1st time step

    % If there were no failures compute a new h.
    if (nofailed && err < 0.9*rtol)
        h=h./max(0.2,1.25*(err./rtol).^RK_pow);
    end
    
    % Advance the integration one step.
    t = tnew;     y = ynew;

    % This flag = 1 because, now I have completed one full time step. 
    % Now to start another time step I need to compute chemical equilibrium
    % In other words, Chemical Equilibrium is computed only at the
    % beginning of time step and NOT for the intermediate time steps of
    % RK45/RK23
    CalcChemEqFlag=1;
    [f(:,1),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec] = ...
        CalcFluxFVSlip(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, dz);

    %Having calculated the flux at start of timestep we don't want to 
    %compute the equilibrium in intermediate time steps. So, set this 
    %flag = 0 again.
    CalcChemEqFlag=0;

        
   
    phiVec=y(cSize(1)*cSize(2)+1:cSize(1)*cSize(2)+N);
    
    OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*phiVec,N);
    dphidzVec=OUT(:,1);  

    JAreaRatioCellVec=y(cSize(1)*cSize(2)+N+1:end);
    AreaRatioCellVec=JAreaRatioCellVec./dphidzVec;
    
    qMat=reshape(y(1:cSize(1)*cSize(2)),cSize);
    cMat=qMat./(JAreaRatioCellVec*ones(1,Nspecies));
    
    if PrepareGridStage
        t=0;
    end

    % Plot solution
    PlotTime=0;
    Cond1=((t>tNextPlot && DeltatNextPlot~=0) || (t==0 && ODEcounter>CounterNextPlot));
    Cond2=(ODEcounter>CounterNextPlot);% && DeltatNextPlot~=0);

    if ~exist('Res','var')
        Cond=1;
    elseif DeltatNextPlot~=0
        Cond=Cond1;
    else
        Cond=Cond2;
    end

    if t>=tfinal
        Cond=1;
    end

    if Cond
        PlotCounter=PlotCounter+1;
        if exist('Res','var')
            Indtmp=size(Res,1)+1;
        else
            Indtmp=1;
        end

        deltayReshaped=reshape(deltay./h,[N,Nspecies+2]);

        dphidzMat=dphidzVec*ones(Nspecies,1)';
        dxMat=dz*dphidzMat;
        
        % Res(Indtmp,:)=sqrt(sum((deltayReshaped).^2,1)/N);
         Res(Indtmp,:)=max(abs(deltayReshaped(:,1:Nspecies+1)),[],1);        
        % disp (['Res=',num2str(Res(end,:)./Res(1,:))]);

        tVecOut(PlotCounter)=t;
        cMatAllTimes(:,:,PlotCounter)=[cMat,cH'];
        %phiVecAllTimes(:,PlotCounter)=phiVec;
        uAllTimes(PlotCounter)=uCST;
        phiVecAllTimes(:,1,PlotCounter)=phiVec;
        phiVecAllTimes(:,2,PlotCounter)=phiVec-uCST*t;
        AreaRatioVecAllTimes(:,PlotCounter)=AreaRatioCellVec;
        JAreaRatioVecAllTimes(:,PlotCounter)=JAreaRatioCellVec;
        SigVecAllTime(:,PlotCounter)=SigVec;
        pHVecAllTimes(:,PlotCounter)=-log10(cH);
        
        keyIn = get(gcf, 'CurrentCharacter');
        if (Res(end,end)/Res(1,end)<0.2 && PrepareGridStage) || strcmp(keyIn,'m')
            %tic
            disp('Migration begins');
            PrepareGridStage=0;
            set(gcf, 'CurrentCharacter','t');
            clear Res
        end

        set(handles.DataHolder,'UserData',[phiVec,cMat]);
        PlotCurrentSelection(handles,phiVec,cMat,muMat,DMat,dPotendxMat,cH,t,h,err,SigVec,dphidzVec, AreaRatioCellVec);

        tNextPlot=t+DeltatNextPlot;
        CounterNextPlot=ODEcounter+DeltaCounterNextPlot;

        PauseButton=get(handles.PauseButton,'value');
        if PauseButton==1
            while PauseButton
                pause(0.1);
                PauseButton=get(handles.PauseButton,'value');
            end %while
        end % if

        StopButton = get(handles.StopSaveButton,'value');
        if StopButton==1 || t>=tfinal
            SimulationStoppedFlag=1;
            set(handles.StopSaveButton,'value',0);
            set(handles.MainAxes,'userdata',cMat);
            break;
        end % if

        StopButton = get(handles.StopDiscardButton,'value');
        if StopButton==1
            SimulationStoppedFlag=0;
            set(handles.StopDiscardButton,'value',0);
            set(handles.MainAxes,'userdata',cMat);
            break;
        end % if

    end
end %while over all time steps
%comments for improvements: 
%bring the first time step if() condition out of this while~done loop

if SimulationStoppedFlag
    dtInit=h;
    save ([(get(handles.FilenameEdit,'String')),'at'],'Nspecies','cMat', ...
        'phiVec','L1','L2','N','t','dtInit','cMatAllTimes','tVecOut','phiVecAllTimes',....
        'uAllTimes','cH','muMat','Res','dPotendxMat', 'AreaRatioVecAllTimes', 'JAreaRatioVecAllTimes', 'SigVecAllTime', 'pHVecAllTimes');
  
% if SimulationStoppedFlag
%     dtInit=h;
%     [(get(handles.FilenameEdit,'String')),'at']
%     
% %     save ([strtok(get(handles.FilenameEdit,'String'),'.'),'.mat'],'Nspecies','cMat', ...
% %         'phiVec','L1','L2','N','t','dtInit','cMatAllTimes','tVecOut','phiVecAllTimes',....
% %         'uAllTimes','cH','muMat','Res','dPotendxMat', 'AreaRatioVecAllTimes', 'JAreaRatioVecAllTimes', 'SigVec'); 
    
end

%--------------------------------------------------------------------------
% Calculates the Flux
% This function is called in the LocalOde45 function to compute flux 
% at each time step of integration
% at Start of every RK time steps CalcChemEqFlag = 1
% for intermediate RK time steps CalcChemEqFlag = 0
%--------------------------------------------------------------------------
function  [dydt,cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec]= ...
    CalcFlux(y0,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag)
global PlotTime
global PrepareGridStage PrevTotalCost Rmu Temp F met2lit
global muH muOH Kw epsilon visc hChannel bPressure DeltaP betaDispersion
global uCST Ngrid PolDeg

[A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,PolDeg,...
    EquiPLength,LMat,LCube,muCube,ValCube,DCube,N,Nspecies,...
    zetaPot,h,DisableChemFlag,SpatialDiscFlag]=ARGS{:};
cSize=[N,Nspecies];
cMat=reshape(y0(1:cSize(1)*cSize(2)),cSize);
phiVec=y0(cSize(1)*cSize(2)+1:end);

%--------------------
% Solve Equilibrium
%DisableChemFlag is only for testing some generic cases.
%It is not actually inputed from the GUI
if DisableChemFlag 

    muMat=F*sum(ValCube.*muCube,3)';
    DMat=DCube(:,:,1)';
    alphaMat=F^2*sum(ValCube.^2.*muCube,3)';
    betaMat=sum(F*ValCube.*DCube,3)';

    %ValTmp=sum(ValCube,3)';
  %  cMat(:,end)=-sum(ValTmp(:,1:end-1).*cMat(:,1:end-1),2)./ValTmp(:,end);
else
    if CalcChemEqFlag
        [cizCube,cH,cHCubePower,gizCube]=LzCalcEquilibrium(cH,LCube,cMat, ValCube);
        
        muMat=F*sum(ValCube.*muCube.*gizCube,3)';
        DMat=sum(DCube.*gizCube,3)'; 
        
%         %sieving matrix start----------------------
%         SievingRatio1=2.5;
%         SievingRatio2=4;
%         SievingVec1=(1+1/SievingRatio1)/2-(1-1/SievingRatio1)/2*tanh(500*(phiVec-0.012)); 
%         SievingVec2=(1+1/SievingRatio2)/2-(1-1/SievingRatio2)/2*tanh(500*(phiVec-0.012)); 
%         muMat(:,4)=muMat( :, 4).*SievingVec1;
%         muMat(:,5)=muMat( :, 5).*SievingVec2;
%         DMat(:,4)=DMat( :, 4).*SievingVec1;
%         DMat(:,5)=DMat( :, 5).*SievingVec2;
%         %sieving matrix end----------------------
%         
        
 
        alphaMat=F^2*sum(ValCube.^2.*muCube.*gizCube,3)';
        betaMat=F*sum(ValCube.*DCube.*gizCube,3)';
    end
end % if CalcChemEqFlag

    % Calculate velocity (pressure driven)
    L=phiVec(end)-phiVec(1);
    uBulk = hChannel^2/(bPressure*visc*L)*(-DeltaP); %assuming constant area
    
if CalcChemEqFlag
    % Dispersion
    PeMat=uBulk*(hChannel)./DMat;                   %h=diameter  h/2=radius
    DMat=DMat.*(1+betaDispersion*PeMat.^2);
end % if CalcChemEqFlag



if DisableChemFlag
    SigVec = sum(alphaMat'.*cMat',1)';
    SVec   = sum(betaMat'.*cMat',1)';
else
    SigVec = sum(alphaMat'.*cMat',1)' + F^2*(muH*cH*met2lit+muOH*(Kw./cH)*met2lit)';
    SVec   = sum(betaMat'.*cMat',1)'  + Rmu*Temp*(muH*cH'*met2lit-muOH*(Kw./cH)'*met2lit);
end


if SteadyStateFlag %&& FirstTimeFlag   % Automatic velocity, based on LE (first chemical on list)
    uCST = -Cur/SigVec(end)*muMat(end,1)-uBulk; %uBluk at LE location??
end

% Derivatives requried for adaptive grid  (but also for next steps)
%gives out 1st derivatives of cH', phiVec, cMat
%A1.a,b,c are diagonal terms, in matrices Af'=Bf
OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[cH',phiVec,cMat],N);
dcHdz=OUT(:,1);   dphidzVec=OUT(:,2);  dcdzMat=OUT(:,3:end);  
dphidzMat=dphidzVec*ones(Nspecies,1)';

if AdaptGrid.Coeff~=0
    mu_c_Char=max(max(abs(muMat.*cMat)));  % characteristic mobility*concentration
    L=phiVec(end)-phiVec(1); %length of domain
    AdaptGrid.wtConst=((N-AdaptGrid.PointsAtInterface)/AdaptGrid.PointsAtInterface)* ...
        Rmu*Temp*mu_c_Char/(abs(Cur)*L);
    
    wttt=abs([dcdzMat,dcHdz]);
    wtt=wttt./(ones(N,1)*max([cMat,cH'],[],1));
    wt=max(wtt,[],2).^(AdaptGrid.Power);

    wt=wt+dphidzVec;
    cost=AdaptGrid.wtConst+wt/(sum(wt)/N);

    if AdaptGrid.Nconv~=0
        cost=conv(cost,[1:AdaptGrid.Nconv,AdaptGrid.Nconv:-1:1]/AdaptGrid.Nconv^2);
        cost=cost(AdaptGrid.Nconv:end-AdaptGrid.Nconv);
    end
    costTakeDeriv=cost.*dphidzVec;
    costDeriv=LinSolve(A1.a,A1.b,A1.c,B1.full*costTakeDeriv,N);

    gradVec=costDeriv;
    RHS=(A2.full(2:end-1,2:end-1)*gradVec(2:end-1));
    gradSmoothVec=[0;LinSolve(M.a(2:end-1),M.b(2:end-1),M.c(2:end-1),RHS,N-2);0];

    dphidtVec=0.5*AdaptGrid.Coeff/AdaptGrid.wtConst*gradSmoothVec;
    dphidtVec(1)=0; dphidtVec(end)=0;
else
    % PrepareGridStage=0;
    dphidtVec=phiVec*0;
    cost=dphidtVec+1;
end
dcdxMat=dcdzMat./dphidzMat;
dphidtMat=dphidtVec*ones(Nspecies,1)';

% If in prepare grid stage - calculate derivaties and return
% No need to solve equations
if PrepareGridStage
    dphidtVec = 10*dphidtVec;
    dcdtMat  = 10*dphidtMat.*dcdxMat;
    dydt=[reshape(dcdtMat,numel(dcdtMat),1);dphidtVec];
    tStepMat=DMat*0;
    dPotendxMat=DMat*0;
    SigVec=dphidtVec*0;
    return;
end


% First deriv - the rest of the derives, not calculated before
% d/dz(cH), d/dz(x), d/dz(c) were computed before
OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[SVec,SigVec,muMat.*cMat,DMat.*cMat],N);
dSdzVec=OUT(:,1);    dSigdzVec=OUT(:,2);
dmucdzMat=OUT(:,3:2+Nspecies);
dDcdzMat=OUT(:,3+Nspecies:2+2*Nspecies);

%Second deriv
OUT=LinSolve(A2.a,A2.b,A2.c,B2.full*[SVec,phiVec,DMat.*cMat],N);
d2Sdz2Vec=OUT(:,1);  d2phidz2Vec=OUT(:,2);  d2Dcdz2Mat=OUT(:,3:2+Nspecies);
d2phidz2Mat=d2phidz2Vec*ones(Nspecies,1)';
    
% Electric potential
dPotendzVec=-1./SigVec.*(Cur*dphidzVec+dSdzVec);
dPotendzMat=dPotendzVec*ones(Nspecies,1)';

dDcdxMat = dDcdzMat./dphidzMat;
dPotendxMat = dPotendzMat./dphidzMat;


d2Dcdx2Mat = (dphidzMat.^-2).*(d2Dcdz2Mat-d2phidz2Mat.*dDcdxMat);
d2Potendx2Vec = -1./(SigVec.*dphidzVec.^2).*(d2Sdz2Vec-dSdzVec.*d2phidz2Vec./dphidzVec+dSigdzVec.*dPotendzVec);
d2Potendx2Mat=d2Potendx2Vec*ones(Nspecies,1)';
dmucdxMat=dmucdzMat./dphidzMat;


switch (SpatialDiscFlag)
    case ('Upwind')
        A1.UW_Left(1,:)=A1.UW_Left(2,:);
        A1.UW_Right(end,:)=A1.UW_Right(end-1,:);

        Flux_x  = (muMat.*dPotendxMat);
        Flux_x2 = -(uBulk+uCST);
        dfdz_UW_Mat=0.5*A1.UW_Left*((Flux_x-abs(Flux_x)).*cMat)+0.5*A1.UW_Right*((Flux_x+abs(Flux_x)).*cMat);

        dfdz_UW_Mat=dfdz_UW_Mat+A1.UW_Right*(Flux_x2.*cMat);
        dfdx_UW_Mat=dfdz_UW_Mat./dphidzMat;
        %dc/dt = d^2/dx^2(D*c) + d/dx(Flux_electric) + x_t*d/dx(c)
        dcdtMat = d2Dcdx2Mat + dfdx_UW_Mat + dphidtMat.*dcdxMat; %equation (4.5) in submitted paper
        
    case ('SLIP')
        LimiterFlag='minmod';        %'minmod' / 'superbee' / 'none'
              
        aMat=muMat.*dPotendxMat-uCST-uBulk;
        fluxMat = aMat.*cMat;                   % size N    (with backward differencing at ends)

        aMat1=muMat.*dPotendxMat;
        aMat2=muMat.*dPotendxMat-uCST-uBulk+dphidtMat;

        alpMat1=-0.5*0.5*(abs(aMat1(1:end-1,:))+abs(aMat1(2:end,:)));      % size: N-2
        alpMat2=-0.5*0.5*(abs(aMat2(1:end-1,:))+abs(aMat2(2:end,:)));      % size: N-2

        CC(:,:,1)=alpMat1;
        CC(:,:,2)=alpMat2;

        alpMat=min(CC,[],3);
        %alpMat=repmat(min(alpMat1,[],2),1,Nspecies);

        dcMat=diff(cMat,[],1);     dcMat=[dcMat(1,:);dcMat;dcMat(end,:)];
        LimitMat=LimiterFunc(dcMat(3:end,:),dcMat(1:end-2,:),LimiterFlag);
        dMat=alpMat.*(dcMat(2:end-1,:)-LimitMat);
        NumFluxMat=0.5*(fluxMat(1:end-1,:)+fluxMat(2:end,:)) - dMat;
        dz=(phiVec(end)-phiVec(1))/N;
               
        dcdtMat=[(fluxMat(2,:)-fluxMat(1,:))/dz;diff(NumFluxMat,[],1)/dz; ...
            (fluxMat(end,:)-fluxMat(end-1,:))/dz]./dphidzMat + dphidtMat.*dcdxMat + d2Dcdx2Mat;
               
    otherwise
        dcdtMat  = d2Dcdx2Mat+dmucdxMat.*dPotendxMat+muMat.*d2Potendx2Mat.*cMat+(dphidtMat-uBulk-uCST).*dcdxMat;
end % switch

    %     Left boundary
    IL=1;
    V1L=muMat(IL,:).*dPotendxMat(IL,:);
    V2L=Cur./SigVec(IL)^2.*(cMat(IL,:).*muMat(IL,:));
    V3L=alphaMat(IL,:)';

    AL=diag(V1L)+V2L'*V3L';      [VAL,DAL] = eig(AL);
    SEigAL=sign(real(diag(DAL)));
    dRLdt=VAL\dcdtMat(IL,:)';
    dRLdt(SEigAL<0)=0;
    dcdtMat(1,:)=real((VAL*dRLdt)') +  (dphidtMat(IL,:)-uBulk-uCST).*dcdxMat(IL,:);

    %     Right boundary
    IR=N;
    V1R=muMat(IR,:).*dPotendxMat(IR,:)-(uBulk+uCST); %changed
    V2R=Cur./SigVec(IR)^2.*(cMat(IR,:).*muMat(IR,:));
    V3R=alphaMat(IR,:)';

    AR=diag(V1R)+V2R'*V3R';
    [VAR,DAR] = eig(AR);
    SEigAR=sign(real(diag(DAR)));
    dRRdt=VAR\dcdtMat(IR,:)';
    dRRdt(SEigAR>0)=0;
    dcdtMat(end,:)=real((VAR*dRRdt)');

dydt=[reshape(dcdtMat,numel(dcdtMat),1);dphidtVec];



function [cH,SigVec,muMat,DMat] = CalculateSpatialProperties(INP,N,FirstTimeFlag,cMat)
global F Rmu  Temp Kw muH muOH met2lit

Nspecies=size(INP,1);
[PMat,PPrimeMat,Q,QPrime,LMat,ValMat,muMat,KaMat,DMat]=EquilibriumPolynomials(INP);
EquiPLength=max(size(PMat,2),size(Q,2));
PolDeg=size(LMat,2); % Polynomial degree
muCube=repmat(reshape(muMat,[Nspecies,1,PolDeg]),[1,N,1]);
DCube=repmat(reshape(DMat,[Nspecies,1,PolDeg]),[1,N,1]);
ValCube=repmat(reshape(ValMat,[Nspecies,1,PolDeg]),[1,N,1]);
LCube=FastRepmatColumns(reshape(LMat,[Nspecies,1,PolDeg]),N);
cH=zeros(1,size(cMat,1));

[cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,cH,PMat,Q,PPrimeMat,QPrime,PolDeg,EquiPLength,N,LMat,LCube,cMat,Nspecies);

muMat=F*sum(ValCube.*muCube.*gizCube,3)';
DMat=sum(DCube.*gizCube,3)';
alphaMat=F^2*sum(ValCube.^2.*muCube.*gizCube,3)';
betaMat=F*sum(ValCube.*DCube.*gizCube,3)';
SigVec = sum(alphaMat'.*cMat',1)' + F^2*(muH*cH*met2lit+muOH*(Kw./cH)*met2lit)';
SVec   = sum(betaMat'.*cMat',1)'  + Rmu*Temp*(muH*cH'*met2lit-muOH*(Kw./cH)'*met2lit);


%--------------------------------------------------------------------------
%  Equilibrium polonomial
%  Called in the LocalOde45 function before the first time step
%  COMMENTS: I would require to call this again and again, if I update
%  pKa's for ionic strength stuff.
%  1) In that case, This would have to be called in CalculateEquilibrium function
%  2) Could remove its call in LocalOde function then
%  3) INP would need to be updated everytime I calculate new kPa
%  4) INP is actually part of ARGS used in MainRun to call LocalOde45
%  5) Would have to take INP as input of CalculateEquilibrium if I need to 
%     call EquilibriumPolynomials(INP) in that function.
%--------------------------------------------------------------------------
function [PMat,PPrimeMat,Q,QPrime,M,zMat,muMat,KaMat,DMat]=EquilibriumPolynomials(INP)
global F Rmu Temp Kw

% PREPARE MATRICES
%------------------
MaxCol=-Inf;
for j=1:size(INP,1)
    MaxCol=max([MaxCol,max(INP{j}(1:3:end))-min(INP{j}(1:3:end))+1]);
end

M=zeros(size(INP,1),MaxCol);

zMat=zeros(size(INP,1),MaxCol);
muMat=zMat; KaMat=zMat; DMat=zMat;

for j=1:size(INP,1)
    zList=INP{j}(1:3:end);
    muList=INP{j}(2:3:end)./(F*abs(zList));  % Defines Santiago mobility
    pKaList=INP{j}(3:3:end);
    KaList=10.^(-pKaList);
    DList=Rmu*Temp*muList; %diffusivity


    [zList,Index]=sort(zList);
    KaList=KaList(Index);
    DList=DList(Index);
    muList=muList(Index);

    Ip1=find(zList==1);     Im1=find(zList==-1);

    zList=[zList(1:Im1),0,zList(Ip1:end)];
    muList=[muList(1:Im1),0,muList(Ip1:end)];
    KaList=[KaList(1:Im1),1,KaList(Ip1:end)];
    DList=[DList(1:Im1),mean(DList),DList(Ip1:end)];

    
    zMat(j,1:length(zList))=zList;
    muMat(j,1:length(muList))=muList;
    KaMat(j,1:length(KaList))=KaList;
    DMat(j,1:length(DList))=DList;

    zListArranged{j}=zList;

    nj=min(zList);    pj=max(zList);

    for z=zList
        if z<0
            M(j,z-nj+1)=prod(KaList(z-nj+1:-nj));
        elseif z>0
            M(j,z-nj+1)=1/prod(KaList(-nj+2:z-nj+1));
        elseif z==0
            M(j,z-nj+1)=1;
        end %if
    end % for z

end %for ij

% CONSTRUCT POLYNOMIALS
%--------------------
Q1=1;
for j=1:size(M,1)
    Q1=conv(Q1,M(j,:));
end %for j
Q2=[-Kw 0 1];
Q=conv(Q1,Q2);

for i=1:size(INP,1)
    tmp=zeros(1,size(M,2));
    tmp(1:length(zListArranged{i}))=zListArranged{i};
    Mmod=M;     Mmod(i,:)=Mmod(i,:).*tmp;

    Pi=1;
    for kl=1:size(Mmod,1)
        Pi=conv(Pi,Mmod(kl,:));
    end %for j
    %PMat(i,:)=Pi;
    Pi=conv([0 1],Pi);  % Convolve with P2
    PMat(i,:)=Pi;

    PiPrime=Pi.*([1:length(Pi)]-1);   PiPrime=[PiPrime(2:end),0];
    PPrimeMat(i,:)=PiPrime;
end %for i
%P2=[0 1];

% Calculate polynomial derivatives
QPrime=Q.*([1:length(Q)]-1);   QPrime=[QPrime(2:end),0];

SizeDiff=size(Q,2)-size(PMat,2);
if SizeDiff>0
    PMat=[PMat,repmat(PMat(:,1)*0,1,SizeDiff)];
    PPrimeMat=[PPrimeMat,repmat(PMat(:,1)*0,1,SizeDiff)];
elseif SizeDiff<0
    Q=[Q,repmat(0,1,SizeDiff)];
    QPrime=[QPrime,repmat(0,1,SizeDiff)];
end

PMat=sparse(PMat);
Q=sparse(Q);

% ---------------------------------

%--------------------------------------------------------------------------
%  Calculate Chemical Equilibrium
%  Using the output of EquilibriumPolynomials, i.e. P, Q etc this function
%  a) computes C_H = [H+].
%  b) knowing C_H and LMat, it computes g_iz
%  Returns c_iz, c_h, g_iz 
%--------------------------------------------------------------------------
function [cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
    cH,PMat,Q,PPrimeMat,QPrime,PolDeg,EquiPLength,N,LMat,LCube,cMat,Nspecies)
global met2lit

% EQUILIBRIUM
if FirstTimeFlag==1  % Solve for pH the first time
    %--------------------------------------------------------------------------
    % Setup constant matrices -  THIS PART RUNS ONLY ONCE PER SIMULATION
    %--------------------------------------------------------------------------
    % Get equilibrium polynomials

    for ij=1:size(cMat,1)
        cTot=cMat(ij,:)/met2lit;   % covert back to mol/liter
        cTotRep=cTot'*ones(1,size(PMat,2));
        P=sum(cTotRep.*PMat,1);

        Poly=zeros(1,EquiPLength);

        Poly(1:size(P,2))=Poly(1:size(P,2))+P;
        Poly(1:size(Q,2))=Poly(1:size(Q,2))+Q;
        Poly=fliplr(Poly);

        roo=roots(Poly);
        roo=roo(imag(roo)==0);
        cH(ij)=roo(roo>0);          % NOTE THAT cH is in mol/lit
    end %for ij
else
    % Solve for pH using Newton Raphson
    count=0;        cHPrev=ones(1,N);
    %            while max(abs(cHPrev-cH)) > max([abs(cH),abs(cHPrev)])*1E-3
    while norm((cHPrev-cH)./max([abs(cH);abs(cHPrev)],[],1),inf) > 1E-6
        count=count+1;
        cHPrev=cH;
        EquicHPolMat=[ones(1,N);cumprod(ones(EquiPLength-1,1)*cH,1)];

        fcH=sum(cMat'.*(PMat*EquicHPolMat),1)/met2lit+Q*EquicHPolMat;
        fPrimecH=sum(cMat'.*(PPrimeMat*EquicHPolMat),1)/met2lit+QPrime*EquicHPolMat;
        cH=cH-fcH./fPrimecH;

        if count>5
            count;
        end
        if count>100
            disp('Too many iterations on cH. Returning');
            return;
        end

    end
end %if first time

cHPolMat=[ones(1,N);cumprod(ones(PolDeg-1,1)*cH,1)];
M1Cube=FastRepmatPages(cMat'./(sparse(LMat)*cHPolMat),PolDeg);
cHCubePower=FastRepmatRows(permute(cHPolMat,[3,2,1]),Nspecies);
cizCube=LCube.*cHCubePower.*M1Cube; %equation 5 in J. Chroma A
gizCube=LCube.*cHCubePower.*FastRepmatPages(1./(sparse(LMat)*cHPolMat),PolDeg);

%helper function in plotting data
function ChangePlotVariable(handles)
global met2lit

N = str2num (get(handles.GridPointsEdit,'string'));
InputTable=cell(handles.InputTable.getData);
InputTable=InputTable(:,2:end);

for ij=1:size(InputTable,1)
    INP{ij,1}=str2num(InputTable{ij,2});
end %for ij

FirstTimeFlag=1;
AxesData=get(handles.DataHolder,'UserData');

if isempty(AxesData)
    L = str2num (get(handles.DomainLengthEdit,'string'));
    xVec = linspace(0,L*1e-3,N)';
    cMat = xVec + nan; cH = xVec + nan; muMat = xVec + nan; DMat = xVec + nan; 
    SigVec = xVec + nan;
else
    xVec=AxesData(:,1); cMat=AxesData(:,2:end); %earlier *met2lit;
    [cH,SigVec,muMat,DMat] = CalculateSpatialProperties(INP,N,FirstTimeFlag,cMat);
end

t=NaN; h=NaN; err=NaN; dphidzVec=xVec+nan; dPotendxMat=xVec + nan;

%%% #plotarea
AreaFunctionText =  get(handles.AreaVariation_Edit, 'String');
if isempty(AreaFunctionText)  
    AreaFunctionText='@(x) 1 + 0*x'; 
    set(handles.AreaVariation_Edit, 'String', AreaFunctionText);
end
try
A = feval( eval(AreaFunctionText), xVec);
AreaRatio=A./max(A);
catch err
if exist('err')
    errordlg('Incorrect function for Area Variation. Example input in vectorized format: @(x) 1 + 0.1*x + 0.2*x.^2 + 0.1*tanh(x)');
    AreaRatio = xVec + nan;
end

end

PlotCurrentSelection(handles,xVec,cMat,muMat,DMat,dPotendxMat,cH,t,h,err,SigVec,dphidzVec, AreaRatio);

function B=FastRepmatColumns(A,reps)
sA=size(A);
B=zeros([sA(1),reps,sA(3)]);
for ij=1:reps
    B(:,ij,:)=A;
end

function B=FastRepmatPages(A,reps)
sA=size(A);
B=zeros([sA(1),sA(2),reps]);
for ij=1:reps
    B(:,:,ij)=A;
end


function B=FastRepmatRows(A,reps)
sA=size(A);
B=zeros([reps,sA(2),sA(3)]);
for ij=1:reps
    B(ij,:,:)=A;
end

%------------------------------------------------------------------
function LVec=LimiterFunc(x,y,LimiterFlag,varargin)

switch (LimiterFlag)
    case ('minmod');
        z(:,:,1)=x; z(:,:,2)=y;
        LVec=0.5*(sign(x)+sign(y)).*min(abs(z),[],3);
    case ('superbee');
        z1(:,:,1)=2*x; z1(:,:,2)=y;
        z2(:,:,1)=x;   z2(:,:,2)=2*y;
              
        m1(:,:,1)=min(abs(z1),[],3);
        m1(:,:,2)=min(abs(z2),[],3);
        
        LVec=0.5*(sign(x)+sign(y)).*max(m1,[],3);
    case ('none')
        LVec=x*0;
        
    case ('ELED') %q=2 gives van-leer q: power in expression of D
        q=2;
        z_xy(:,:,1)=abs(x)+abs(y)+1.0e-14; z_xy(:,:,2)=z_xy(:,:,1)*0; %+0.0*varargin{1}.^1.5;
        D=1-(abs((x-y)./max(z_xy,[],3))).^q;
        LVec=0.5*D.*(x+y);
        
end % switch

function RHS=LinSolve(a,b,c,RHS,N)
%a(1)=a(1)+0.0; b(1)=b(1)+0.0; c(1)=c(1)+0.0; RHS(1)=RHS(1)+0.0; N=N+0.0;
% if (min(b==0) && min(c==0) && min(a==1))
%     OUT=RHS;
% else
%     OUT=TriDiagSolve(a,b,c,RHS);

%a(1)=a(1)+0.0; b(1)=b(1)+0.0; c(1)=c(1)+0.0; RHS(1)=RHS(1)+0.0; N=N+0.0;
%b a c 
if (min(b==0) && min(c==0) && min(a==1))
    OUT=RHS;
else
    A=diag(a,0)+diag(b,-1)+diag(c,1);
    RHS=A\RHS;
end
%--------------------------------------------------------------------------





% Time stepping based on ode45
%probably make another function LocalOde23, and remove all if conditions 
%for RK45 or RK23
function LocalOde45(y0,tspan,ARGS)
global Nspecies betaDispersion bPressure hChannel Ngrid
global PlotTime errVec PrepareGridStage AChannel
global uCST
[cSize,DeltatNextPlot,DeltaCounterNextPlot,dtInit,dtMin,dtMax,rtol,atol,...
    normcontrol,RK_A,RK_B,RK_E,RK_pow,uCST,A1,B1,A2,B2,M,Cur,INP,...
    AdaptGrid,RKOrder,SteadyStateFlag,L1,L2,zetaPot,hChannel,bPressure, ...
    DeltaP,betaDispersion,Res,DisableChemFlag,SpatialDiscFlag,handles]=ARGS{:};
N=cSize(1);

% Additional parameters
cMat=reshape(y0(1:cSize(1)*cSize(2)),cSize);
%phi refers to the physical domain
phiVec=y0(cSize(1)*cSize(2)+1:end);

% Additional parameters
t0=tspan(1); tfinal=tspan(end);
tdir = sign(tfinal - t0);
threshold = atol / rtol;
htspan=dtInit;                  % Initial time step
hmin=dtMin;                 % Minimum time step
hmax=dtMax;
neq=length(y0);
normy=norm(y0);
tNextPlot=0;
Nspecies=cSize(2);

t=t0;  y=y0;
FirstTimeFlag=1; %Flag to indicate the first time step
ODEcounter=0;
CounterNextPlot=1;
PlotCounter=0;
cMatAllTimes=[]; phiVecAllTimes=[]; tVecOut=[];


% Create Polynomials for Equilibrium reactions
[PMat,PPrimeMat,Q,QPrime,LMat,ValMat,muMat,KaMat,DMat]=EquilibriumPolynomials(INP);
EquiPLength=max(size(PMat,2),size(Q,2));
PolDeg=size(LMat,2); % Polynomial degree

%Cube here refers to 3-d array
%cols: grid points, rows: no. of species depth: polynomial degree
muCube=repmat(reshape(muMat,[Nspecies,1,PolDeg]),[1,N,1]); %mobilities
DCube=repmat(reshape(DMat,[Nspecies,1,PolDeg]),[1,N,1]); %diffusivities
ValCube=repmat(reshape(ValMat,[Nspecies,1,PolDeg]),[1,N,1]); %valence
LCube=FastRepmatColumns(reshape(LMat,[Nspecies,1,PolDeg]),N); %???
cH=zeros(1,size(cMat,1)); %set initial [H+] = 0 
Ngrid=length(cH);

% THE MAIN LOOP
done = false;
%time stepping starts here
while ~done
    ODEcounter=ODEcounter+1;


    if FirstTimeFlag==1 %initialy it is set = 1 for first run
        CalcChemEqFlag=1; %initially evaluate chemical equilibrium
        errVec=ones(N,1); [muMat,DMat,alphaMat,betaMat]=deal([]);
        h=ones(1,N*(Nspecies+1));  %Just for initial flux calculation
        GridCost=zeros(N,1);
        
        %run CalculateEquilibrium when FirstTimeFlag=1
        %otherwise run LzCalcEquilibrium inside CalcFlux
        [cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
    cH,PMat,Q,PPrimeMat,QPrime,PolDeg,EquiPLength,N,LMat,LCube,cMat,Nspecies);
        
        
        %Define arguments to call the CalcFlux function
        ARGS={A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag, PolDeg,EquiPLength,...
            LMat,LCube,muCube,ValCube,...
            DCube,N,Nspecies,zetaPot,h,DisableChemFlag,SpatialDiscFlag};
        
        %call CalcFlux function
        [f0,cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec]= ...
            CalcFlux(y0,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag);

        % Initialize method parameters.
        if RKOrder==45,       f = zeros(neq,7);
        elseif RKOrder==23,   f = zeros(neq,4);
        end
        f(:,1) = f0;

        % Compute an initial step size h using y'(t).
        absh = min(hmax, htspan);
        rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^RK_pow);
        if absh * rh > 1
            absh = 1 / rh;
        end
        absh = max(absh, hmin);
        absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
        h = tdir * absh;

        FirstTimeFlag=0; %1st step apready done so set this flag = 0
        CalcChemEqFlag=0; %I dont want to compute chemical equilibrim for 
                          %intermediate steps
                          
    end %end of the First time's if condition
    %could probably move this out of the loop.. 
    %why do you want to check if it is the first time step even after 
    %you have computed N steps?????

    h = min(hmax, max(hmin, h));
    if h>(tfinal-t)
        h=tfinal-t;
    end

    
    %NOTE: Chemical equilibrium is not computed in intermediate steps
    % LOOP FOR ADVANCING ONE STEP.
    nofailed = true;                      % no failed attempts
    while true
        
        ARGS={A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,PolDeg,EquiPLength,...
            LMat,LCube,muCube,ValCube,...
            DCube,N,Nspecies,zetaPot,h,DisableChemFlag,SpatialDiscFlag};
        
        [f(:,2),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec]= ...
            CalcFlux(y+h.*(f*RK_B(:,1)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag);
        [f(:,3),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec]=...
            CalcFlux(y+h.*(f*RK_B(:,2)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag);

        if RKOrder==45
            [f(:,4),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec]=...
                CalcFlux(y+h.*(f*RK_B(:,3)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag);
            
            [f(:,5),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec]=...
                CalcFlux(y+h.*(f*RK_B(:,4)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag);
            
            [f(:,6),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec]=...
                CalcFlux(y+h.*(f*RK_B(:,5)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag);
        end

        if     RKOrder==45,  tnew = t + min(h*RK_A(6));
        elseif RKOrder==23,  tnew = t + min(h*RK_A(3));
        end

        if RKOrder==45,
            deltay=h.*(f*RK_B(:,6));
            ynew = y + deltay;
            [f(:,7),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec] = ...
                CalcFlux(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag);
            
        elseif RKOrder==23,
            deltay=h.*(f*RK_B(:,3));
            ynew = y + deltay;
            [f(:,4),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec] = ...
                CalcFlux(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag);
        end

        % Estimate the error.
            errAll= h.*abs((f * RK_E) ./ max(max(abs(y),abs(ynew)),threshold));
            errMat=reshape(errAll,[N Nspecies+1]);
            errVec=max(errMat,[],2);           
            err=norm(errVec);            
            
        if err < rtol                       % Successful step - go to next time step
            break;
            
        else                                % Failed step - calculate new time step
            if nofailed
                nofailed = false;
                h = max(hmin, h.*max(0.1, 0.8*(rtol./err).^RK_pow));
            else
                h = max(hmin, 0.5 * h);
            end
            done = false;
        end % if err > rtol
    end % while true. End of while loop over the time steps after 1st time step

    % If there were no failures compute a new h.
    if (nofailed && err < 0.9*rtol)
        h=h./max(0.2,1.25*(err./rtol).^RK_pow);
    end
    
    % Advance the integration one step.
    t = tnew;     y = ynew;

    % This flag = 1 because, now I have completed one full time step. 
    % Now to start another time step I need to compute chemical equilibrium
    % In other words, Chemical Equilibrium is computed only at the
    % beginning of time step and NOT for the intermediate time steps of
    % RK45/RK23
    CalcChemEqFlag=1;
    [f(:,1),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec] = ...
        CalcFlux(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag);

    %Having calculated the flux at start of timestep we don't want to 
    %compute the equilibrium in intermediate time steps. So, set this 
    %flag = 0 again.
    CalcChemEqFlag=0;

    cMat=reshape(y(1:cSize(1)*cSize(2)),cSize);
    phiVec=y(cSize(1)*cSize(2)+1:end);
      
    if PrepareGridStage
        t=0;
    end

    % Plot solution
    PlotTime=0;
    Cond1=((t>tNextPlot && DeltatNextPlot~=0) || (t==0 && ODEcounter>CounterNextPlot));
    Cond2=(ODEcounter>CounterNextPlot);% && DeltatNextPlot~=0);

    if ~exist('Res','var')
        Cond=1;
    elseif DeltatNextPlot~=0
        Cond=Cond1;
    else
        Cond=Cond2;
    end

    if t>=tfinal
        Cond=1;
    end

    if Cond
        PlotCounter=PlotCounter+1;
        if exist('Res','var')
            Indtmp=size(Res,1)+1;
        else
            Indtmp=1;
        end

        deltayReshaped=reshape(deltay./h,[N,Nspecies+1]);

        dz=(phiVec(end)-phiVec(1))/N;
        dphidzMat=dphidzVec*ones(Nspecies,1)';
        dxMat=dz*dphidzMat;
        
        % Res(Indtmp,:)=sqrt(sum((deltayReshaped).^2,1)/N);
         Res(Indtmp,:)=max(abs(deltayReshaped),[],1);        
        % disp (['Res=',num2str(Res(end,:)./Res(1,:))]);

        tVecOut(PlotCounter)=t;
        cMatAllTimes(:,:,PlotCounter)=[cMat,cH'];
        %phiVecAllTimes(:,PlotCounter)=phiVec;
        uAllTimes(PlotCounter)=uCST;
        phiVecAllTimes(:,1,PlotCounter)=phiVec;
        phiVecAllTimes(:,2,PlotCounter)=phiVec-uCST*t;
        SigVecAllTime(:,PlotCounter)=SigVec;
        pHVecAllTimes(:,PlotCounter)=-log10(cH);
          
        keyIn = get(gcf, 'CurrentCharacter');
        if (Res(end,end)/Res(1,end)<0.2 && PrepareGridStage) || strcmp(keyIn,'m')
            %tic
            disp('Migration begins');
            PrepareGridStage=0;
            set(gcf, 'CurrentCharacter','t');
            clear Res
        end

        set(handles.DataHolder,'UserData',[phiVec,cMat]);
        PlotCurrentSelection(handles,phiVec,cMat,muMat,DMat,dPotendxMat,cH,t,h,err,SigVec,dphidzVec, phiVec*0+AChannel);

        tNextPlot=t+DeltatNextPlot;
        CounterNextPlot=ODEcounter+DeltaCounterNextPlot;

        PauseButton=get(handles.PauseButton,'value');
        if PauseButton==1
            while PauseButton
                pause(0.1);
                PauseButton=get(handles.PauseButton,'value');
            end %while
        end % if

        StopButton = get(handles.StopSaveButton,'value');
        if StopButton==1 || t>=tfinal
            SimulationStoppedFlag=1;
            set(handles.StopSaveButton,'value',0);
            set(handles.MainAxes,'userdata',cMat);
            break;
        end % if

        StopButton = get(handles.StopDiscardButton,'value');
        if StopButton==1
            SimulationStoppedFlag=0;
            set(handles.StopDiscardButton,'value',0);
            set(handles.MainAxes,'userdata',cMat);
            break;
        end % if

    end
end %while over all time steps
%comments for improvements: 
%bring the first time step if() condition out of this while~done loop

if SimulationStoppedFlag
    dtInit=h;
    save ([(get(handles.FilenameEdit,'String')),'at'],'Nspecies','cMat', ...
        'phiVec','L1','L2','N','t','dtInit','cMatAllTimes','tVecOut','phiVecAllTimes',....
        'uAllTimes','cH','muMat','dPotendxMat', 'SigVecAllTime', 'pHVecAllTimes');
  
    
%     save ([strtok(get(handles.FilenameEdit,'String'),'.'),'.mat'],'Nspecies','cMat', ...
%         'phiVec','L1','L2','N','t','dtInit','cMatAllTimes','tVecOut','phiVecAllTimes',....
%         'uAllTimes','cH','muMat','Res','dPotendxMat');
    
    
end

%%original.. uses matrix operations and solves for the cH using vector
%this function uses iterative method on eqn (6) of J. Chroma A paper.
%does not compute the equilibrium polynomials
%--------------------------------------------------------------------------
%  Calculate Chemical Equilibrium
%  Using the output of EquilibriumPolynomials, i.e. P, Q etc this function
%  a) computes C_H = [H+].
%  b) knowing C_H and LMat, it computes g_iz
%  Returns c_iz, c_h, g_iz 
%--------------------------------------------------------------------------
function [cizCube,cH,cHCubePower,gizCube]=LzCalcEquilibrium(cH,LCube,cMat, ValCube)
global Ngrid Nspecies PolDeg met2lit

% options = optimset('Display','iter','TolFun',1e-8);
% [cH, fval, exitflag] = fsolve(@(cH) MyLzFuncNEW(cH,LCube,cMat, ValCube), cH, options);
% %look at the case where function does not converge
cMat=cMat/met2lit;
cHPrev=ones(1,Ngrid);
count=0;
while norm((cHPrev-cH)./max([abs(cH);abs(cHPrev)],[],1),inf) > 1E-6
    count=count+1;
    cHPrev=cH;
    cH=cH-MyLzFuncNEW(cH,LCube,cMat, ValCube);
end
%fprintf('no. of iterations = %g \n', count);

cHPolMat=[ones(1,Ngrid);cumprod(ones(PolDeg-1,1)*cH,1)];
cHCubePower=FastRepmatRows(permute(cHPolMat,[3,2,1]),Nspecies);
Temp=sum((LCube.*cHCubePower),3);
M1Cube=FastRepmatPages(cMat'./Temp,PolDeg);
cizCube=LCube.*cHCubePower.*M1Cube; %equation 5 in J. Chroma A
gizCube=LCube.*cHCubePower.*FastRepmatPages(1./Temp,PolDeg);


%%original.. uses matrix operations and solves for the cH using vector
%this function uses iterative method on eqn (6) of J. Chroma A paper.
%does not compute the equilibrium polynomials
%--------------------------------------------------------------------------
%  Calculate Chemical Equilibrium
%  Using the output of EquilibriumPolynomials, i.e. P, Q etc this function
%  a) computes C_H = [H+].
%  b) knowing C_H and LMat, it computes g_iz
%  Returns c_iz, c_h, g_iz 
%--------------------------------------------------------------------------
function F=MyLzFuncNEW(cH,LCube,cMat, ValCube)
global Ngrid Nspecies PolDeg Kw N met2lit


%step 1: from the vector cH and LCube, construct cizCube
% step 2
% loop over grid points k=1:Ngrid
%     loop over species j=1:Nspecies
%         zList=zListArranged{j};
%         nj=min(zList);    pj=max(zList);
%     
%                     loop z=zList
%                      RHS(k)=RHS(k)-z*cizcube(j,k,z-nj+1);    
%                     loop over z ends
%       loop of species ends
% loop over grid points end
% cH=RHS=RHS+Kw./cH
%Step 3: Repeat step 2 again.. till you converge

cHPolMat=[ones(1,Ngrid);cumprod(ones(PolDeg-1,1)*cH,1)];
cHCubePower=FastRepmatRows(permute(cHPolMat,[3,2,1]),Nspecies);
%Temp=multiprod(TempL, cHPolMat, [1 2], [1]); %faster multiplication
Temp=sum((LCube.*cHCubePower),3);
TempCube=FastRepmatPages(Temp,PolDeg);

M1Cube=FastRepmatPages(cMat'./Temp,PolDeg);
cizCube=LCube.*cHCubePower.*M1Cube; %equation 5 in J. Chroma A

Temp_z = sum((ValCube.*LCube.*cHCubePower),3);
Temp_z_Cube=FastRepmatPages(Temp_z,PolDeg);
RHS_den=sum(sum((cizCube.*ValCube.^2-cizCube.*ValCube.*Temp_z_Cube./TempCube),3),1);

%M1Cube is denominator of eqn 5
%RHS_den=sum(sum((cizCube.*ValCube.^2-cizCube.*ValCube.^2.*LCube.*cHCubePower./TempCube),3),1);
RHS_num=sum(sum((cizCube.*ValCube),3),1);
F_num=RHS_num+cH-Kw./cH;
F_den=RHS_den./cH+1+Kw./(cH.^2);
F=F_num./F_den;

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



function PlotCurrentSelection(handles,xVec,cMat,muMat,DMat,dPotendx,cH,t,h,err,SigVec,dphidzVec, CrossArea)
global met2lit DeltaP
global ColorList

    WhatToPlot=get(handles.WhatToPlotPopup,'value');
    axes(handles.MainAxes); set(handles.MainAxes,'nextplot','replacechildren');

    switch WhatToPlot
        case (1)
        plot(handles.MainAxes,xVec*1E3,cMat/met2lit,'.-');              % concentration is in mol/m^3. Convert to mol/lit
            xlabel(handles.MainAxes,'x [mm]'); ylabel(handles.MainAxes,'c [mol/lit]');
            %ylim([-0.001 1.2*max(max(cMat))/met2lit]);
        case(2)
            plot(handles.MainAxes,xVec*1E3,-log10(cH),'.-');
            xlabel(handles.MainAxes,'x [mm]'); ylabel(handles.MainAxes,'pH');
        case (3)
            plot(handles.MainAxes,xVec*1E3,muMat,'.-');
            xlabel(handles.MainAxes,'x [mm]'); ylabel(handles.MainAxes,'\mu [m^2/Vs]');
        case (4)
            plot(handles.MainAxes,xVec*1E3,dPotendx(:,1),'.-');
            xlabel(handles.MainAxes,'x [mm]'); ylabel(handles.MainAxes,'E [V/m]');
        case (5)
            plot(handles.MainAxes,xVec*1E3,SigVec,'.-');
            xlabel(handles.MainAxes,'x [mm]'); ylabel(handles.MainAxes,'\sigma [S/m]');
        case (6)
            plot(handles.MainAxes,xVec*1E3,1./dphidzVec,'.-');
            xlabel(handles.MainAxes,'x [mm]'); ylabel(handles.MainAxes,'Grid density');
            [MaxDensity,MaxDensityInd]=min(dphidzVec);
            text(xVec(MaxDensityInd)*1E3, 1./MaxDensity,['Max grid density = ',num2str(1./MaxDensity,3)]);
        case (7) %%% #plotarea
            plot(handles.MainAxes,xVec*1E3, CrossArea/CrossArea(1),'.-');
            xlabel(handles.MainAxes,'x [mm]'); ylabel(handles.MainAxes,'Area Ratio [A/A(x=0)]');
                        
    end %switch
      Voltage=trapz(xVec,dPotendx(:,1));
   
  
    title (handles.MainAxes,['t=',num2str(t),' s  \Deltat=',num2str(h),'s  \DeltaV = ',num2str(Voltage,4),' V  \DeltaP = ',num2str(DeltaP,4),' Pa']);
 
    drawnow;

    
    % prepares derivatives depending on the different numerical schemes
function [A,B]=PrepareDerivMatrices(x,Order,Flag)

% Implicit, 6th order compact calculation of first derivative
% A*f'=B*f

h=x(2)-x(1);
% Assuming uniform grid spacing
N=size(x,1);
A=sparse(N,N);    B=sparse(N,N);

switch (Flag)
    case ('Compact')
        if Order==1
            % Forward derivative for first(/last) grid point - 5th order
            A_First=[1 4];
            B_First=(1/h)*[-37/12    2/3    3   -2/3    1/12];

            % Skewed derivative for second/(N-1) grid point - 5th order
            A_Second=[1/6 1 1/2];
            B_Second=(1/h)*[-10/18   -1/2    1    1/18];


            % Central derivative for internal grid points
            alphaI=1/3;
            A_Internal=[alphaI 1 alphaI];
            B_Internal=[-1/3*(4*alphaI-1)*(1/(4*h))  -2/3*(alphaI+2)*(1/(2*h))  0  2/3*(alphaI+2)*(1/(2*h))  1/3*(4*alphaI-1)*(1/(4*h))];


            A(1,1:length(A_First))      = A_First;        B(1,1:length(B_First)) = B_First;
            A(2,1:length(A_Second))     = A_Second;       B(2,1:length(B_Second))= B_Second;

            A(end,end:-1:end-length(A_First)+1)    = A_First;       B(end,end:-1:end-length(B_First)+1) = - B_First;
            A(end-1,end:-1:end-length(A_Second)+1) = A_Second;      B(end-1,end:-1:end-length(B_Second)+1) = - B_Second;


            for ij=3:N-2  % run on all rows
                A(ij,ij-1:ij+1)=A_Internal;
                B(ij,ij-2:ij+2)=B_Internal;
            end

        elseif Order==2
            % Forward derivative for first(/last) grid point - 5th order
            A_First=[1 137/13];
            B_First=(1/h^2)*[1955/156 -4057/156 1117/78 -55/78 -29/156 7/156];

            % Skewed derivative for second/(N-1) grid point - 5th order
            A_Second=[1/10 1 -7/20];
            B_Second=(1/h^2)*[99/80   -3    93/40   -3/5    3/80];


            % Central derivative for internal grid points
            alphaI=2/11;
            aI=4/3*(1-alphaI)*(1/h^2);
            bI=1/3*(-1+10*alphaI)*(1/(4*h^2));
            A_Internal=[alphaI 1 alphaI];
            B_Internal=[bI   aI  -2*(aI+bI)  aI  bI];


            A(1,1:length(A_First))      = A_First;        B(1,1:length(B_First)) = B_First;
            A(2,1:length(A_Second))     = A_Second;       B(2,1:length(B_Second))= B_Second;

            A(end,end:-1:end-length(A_First)+1)    = A_First;       B(end,end:-1:end-length(B_First)+1) = B_First;
            A(end-1,end:-1:end-length(A_Second)+1) = A_Second;      B(end-1,end:-1:end-length(B_Second)+1) = B_Second;


            for ij=3:N-2  % run on all rows
                A(ij,ij-1:ij+1)=A_Internal;
                B(ij,ij-2:ij+2)=B_Internal;
            end

        end % if

    case {'2nd','Upwind','SLIP'}
        if Order==1
            A=eye(N);
            B=diag(ones(1,N-1)/(2*h),1)+diag(-ones(1,N-1)/(2*h),-1);
            B(1,1:2)=[-1/h,1/h];
            B(end,end-1:end)=[-1/h,1/h];
        elseif Order==2
            A=eye(N);
            B=diag(ones(1,N-1)/h^2,1)+diag(ones(1,N-1)/h^2,-1)-2*diag(ones(1,N)/h^2,0);
            B(1,1:3)=[1,-2,1]/h^2;
            B(end,end-2:end)=[1,-2,1]/h^2;
        end
end % switch


%--------------------------------------------------------------------------
% Prepares inputs etc...
%--------------------------------------------------------------------------
function PrepareInput(handles)
global met2lit
FirstTimeFlag=1;
N = str2num (get(handles.GridPointsEdit,'string'));             if isempty(N), set(handles.GridPointsEdit,'string','100'); N=100; end
L = str2num (get(handles.DomainLengthEdit,'string'));           if isempty(L), set(handles.DomainLengthEdit,'string','10'); L=10; end
InjLen = str2num (get(handles.InjectionWidthEdit,'string'));    if isempty(InjLen), set(handles.InjectionWidthEdit,'string','1'); InjLen=1; end
InjLoc = str2num (get(handles.InjectionPointEdit,'string')); if isempty(InjLoc), set(handles.InjectionPointEdit,'string',num2str(L/2)); InjLoc=L/2; end
InputTable = cell(handles.InputTable.getData);
InputTable=InputTable(:,2:end);

if isempty(InputTable)
    return;
end

xVec=linspace(0,L*1E-3,N)';
AnalyteVec=(1+erf((xVec/1E-3-(InjLoc-InjLen/2))))/2-(1+erf((xVec/1E-3-(InjLoc+InjLen/2))))/2;   AnalyteVec(abs(AnalyteVec)<1E-3)=0;
TEVec=(1-erf((xVec/1E-3-(InjLoc-InjLen/2))))/2;    TEVec(abs(TEVec)<1E-3)=0;
LEVec=(1+erf((xVec/1E-3-(InjLoc+InjLen/2))))/2;    LEVec(abs(LEVec)<1E-3)=0;
BackgroundVec=xVec*0+1;

for ij=1:size(InputTable,1)
    INP{ij,1}=str2num(InputTable{ij,2});
    
    Conc=InputTable{ij,3};
    if ~isnumeric(Conc)
        Conc=str2num(Conc);
    end
    
    
    switch(InputTable{ij,4})
        case('LE')
            cMat(:,ij)=Conc*LEVec;
        case('TE')
            cMat(:,ij)=Conc*TEVec;
        case('Analyte')
            cMat(:,ij)=Conc*AnalyteVec;
        case('Background')
            cMat(:,ij)=Conc*BackgroundVec;
    end
end %for ij
[cH,SigVec,muMat,DMat] = CalculateSpatialProperties(INP,N,FirstTimeFlag,cMat);
t=NaN; h=NaN; err=NaN; dphidzVec=xVec+nan; dPotendxMat=xVec+nan; 


%%% #plotarea
AreaFunctionText =  get(handles.AreaVariation_Edit, 'String');  
if isempty(AreaFunctionText)  AreaFunctionText='@(x) 1 + 0*x'; end
A = feval( eval(AreaFunctionText), xVec); 
AreaRatio=A./max(A);

cMat=cMat*met2lit;                % Convert from mol/lit to mol/m^3  
PlotCurrentSelection(handles,xVec,cMat,muMat,DMat,dPotendxMat,cH,t,h,err,SigVec,dphidzVec, AreaRatio);
set(handles.DataHolder,'UserData',[xVec,cMat]);
set(handles.StartButton,'enable','off');

%-------------------------------------------------------




