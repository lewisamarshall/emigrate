function SpressoIonic(varargin)
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
visc=1E-3;              % Dynamic viscosity (water) [Pa s]F

% If no input arguments - run graphical user interface


if isempty(varargin)
    clc;
    warning off all
    % close all;
    %remove paths from previous runs just to make sure
    %rmpath(['SubroutinesIonic']);
    %rmpath(['SubroutinesNoIonic']);
    %now add path
    %addpath(['SubroutinesIonic']);
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

function  [dydt,cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec, ...
    IonicStrength, LCube, Kw_new]= ...
    CalcFluxFVSlip(y0,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
    IonicCalcFlag, IonicStrength, IonicEffectFlag, dz)

global PlotTime
global PrepareGridStage PrevTotalCost Rmu Temp F met2lit
global muH muOH Kw epsilon visc hChannel bPressure DeltaP betaDispersion
global uCST Ngrid PolDeg

[A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,PolDeg, ...
    EquiPLength,LCube, Kw_new, muCube,ValCube,DCube,N,Nspecies,zetaPot,h, ...
    DisableChemFlag,SpatialDiscFlag, zListArranged, KaListCube]=ARGS{:};

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
%dphidzAreaMat=AreaRatioCellMat.*dphidzMat;
cMat=qMat./(AreaRatioCellMat.*dphidzMat);
%cMat=qMat./dphidzMat;

% Derivatives requried for adaptive grid  (but also for next steps)
%gives out 1st derivatives of cH', cMat
%A1.a,b,c are diagonal terms, in matrices Af'=Bf
% %
% % OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[cH',cMat, AreaRatioCellVec.*cH', AreaRatioCellMat.*cMat],N);
% % dcHdz=OUT(:,1); dcdzMat=OUT(:,2:Nspecies+1); dAcHdz=OUT(:,Nspecies+2); dAcdzMat=OUT(:,Nspecies+3:end);


%-----------------------------------------
% Arificial Diffusion for Stabilitizing (size N-1  x Nspecies)
%------------------------------------------
LimiterFlag='ELED';        %'minmod' / 'superbee' / 'VanLeer' /'none'/'ELED'


%--------------------
% Solve Equilibrium
%DisableChemFlag is only for testing some generic cases.
%It is not actually inputed from the GUI
if DisableChemFlag

    muMat=F*sum(ValCube.*muCube,3)';
    DMat=DCube(:,:,1)';
    alphaMat=F^2*sum(ValCube.^2.*muCube,3)';
    betaMat=sum(F*ValCube.*DCube,3)';
    ValTmp=sum(ValCube,3)';

    %  cMat(:,end)=-sum(ValTmp(:,1:end-1).*cMat(:,1:end-1),2)./ValTmp(:,end);
else
    if CalcChemEqFlag
        [cizCube,cH,cHCubePower,gizCube]=LzCalcEquilibrium(cH,LCube,cMat, ValCube, Kw_new);

        if (IonicCalcFlag==1)
            %%%%must call CalculateEquilibrium before CalculateIonicEffects
            [cizCube,cH,cHCubePower,gizCube, IonicStrength, LCube, Kw_new]= ...
                CalculateIonicEffects(IonicEffectFlag, cH,LCube,cMat, ValCube, cizCube, zListArranged, KaListCube, cHCubePower, gizCube);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        IonicStrength=(0.5*sum(sum(ValCube.^2.*cizCube,3),1)/met2lit+0.5*(cH+Kw_new./cH))*IonicEffectFlag;
        %ionic strength includes [H+] and [OH-]
        %%%Now that we have Ionic strength, we can calculate the new MobilityCube
        ICube=repmat(IonicStrength, [Nspecies,1,PolDeg]);
        muIonicCube = OnsagerFuoss(IonicEffectFlag, ValCube, cizCube, muCube);
        muMat=sum(ValCube.*muIonicCube.*gizCube,3)';
        alphaMat=F*sum(ValCube.^2.*muIonicCube.*gizCube,3)';
        DTemp=DCube.*(1-0.2297*sqrt(ICube));
        DMat=sum(DTemp.*gizCube,3)';
        betaMat=F*sum(ValCube.*DTemp.*gizCube,3)';
        %%% CalculateIonicEffects.
    end % if CalcChemEqFlag
end

% Calculate velocity (pressure driven)
L=phiVec(end)-phiVec(1);
uBulk = hChannel^2/(bPressure*visc*L)*(-DeltaP); %uBulk is net flow rate
%constant along axis (Q)


if CalcChemEqFlag
    % Dispersion
    PeMat=uBulk*(hChannel)./DMat;                 %h=diameter  h/2=radius
    DMat=DMat.*(1+betaDispersion*PeMat.^2);
end % if CalcChemEqFlag

if DisableChemFlag
    SigVec = sum(alphaMat'.*cMat',1)';
    SVec   = sum(betaMat'.*cMat',1)';
else
    muH_new=muH;
    muOH_new=muOH;
    %change the Kw since ionic effects are to be considered
    Kw_new=Kw*10.^(-2*getLogActivity(IonicStrength));
    SigVec = sum(alphaMat'.*cMat',1)' + F^2*(muH_new.*cH*met2lit+muOH_new.*(Kw./cH)*met2lit)';
    SVec   = sum(betaMat'.*cMat',1)'  + Rmu*Temp*(muH_new.*cH*met2lit-muOH_new.*(Kw./cH)*met2lit)';
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
CorrectedDiffEdgeMat = dcMat(2:end-1,:)-LimitMat;  % Corrected numerical diffusivity (size (N-1) x Nspecies


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



%%Case 3:
MolecDiffFactorEdgeMat=(AreaRatioCellMat(1:N-1,:)./dphidzMat(1:N-1,:) + ...
                            AreaRatioCellMat(2:N,:)./dphidzMat(2:N,:))/2;
MolecDiffFluxEdgeMat = MolecDiffFactorEdgeMat.*diff(DMat.*cMat,1)./dz; %size (N-1) x Nspecies


% % %--------------------------------
% % cMatExtended=[cMat(1,:); cMat;cMat(end,:)];  % length: N+2
% % dcMat=diff(cMatExtended,[],1);              % length: N+1
% % LimitMat=LimiterFunc(dcMat(3:end,:),dcMat(1:end-2,:),LimiterFlag, dz); %(N-1)
% % CorrectedDiffEdgeMat = dcMat(2:end-1,:)-LimitMat;  % Corrected numerical diffusivity (size (N-1) x Nspecies
% % %-------------

alpTempEdge = (uBulk+uCST) +  0.5*(muMat(1:N-1,:).*Cur./SigMat(1:N-1,:)+muMat(2:N,:).*Cur./SigMat(2:N,:) -...
                        dphidtMat(1:N-1,:).*AreaRatioCellMat(1:N-1,:) - dphidtMat(2:N,:).*AreaRatioCellMat(2:N,:));


alpVecEdge= 0.52*(max(abs(alpTempEdge),[],2));
alpMatEdge=alpVecEdge*ones(1,Nspecies);
CorrectedAdFluxEdgeMat =  (AdFluxCellMat(1:N-1, : ) + AdFluxCellMat(2:N, : ))/2 ...
  								- alpMatEdge.*CorrectedDiffEdgeMat;
%---------------------------------

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

% % % %Apply boundary conditions-----------
% % %    %     Left boundary
%     IL=1;
%     V1L=uBulk + uCST + muMat(IL,:).*Cur./SigVec(IL); %Q + mu*I/sigma
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
%     FluxLeftBC =VAL*(DAL_negative)*(VAL\dcdzL');
%     dqdtMat(1,:)= -GridFluxCellMat(2,:)/dz - FluxLeftBC' + alpVecEdge(1)*dcdzL;
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

%----------------------------------------------------------------
% % % OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*dphidtMat,N); %Fixed BC
% % %     dDphidtDzMat=OUT(:,1:Nspecies); %J=dx/dz
% % % dqdtMat(1,:)=(dphidtMat(1,:).*dcdzMat(1,:) + cMat(1,:).*dDphidtDzMat(1,:));
% % % dqdtMat(N,:)=(dphidtMat(N,:).*dcdzMat(N,:) + cMat(N,:).*dDphidtDzMat(N,:));
%-----------------------------------------------------------------
OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[SVec],N);
dSdzMat=OUT(:,1)*ones(1, Nspecies);
dPotendxMat=-(Cur./AreaRatioCellMat +  dSdzMat./dphidzMat)./SigMat; %Cur is current density of ARatio=1

dydt=[reshape(dqdtMat,numel(dqdtMat),1);dphidtVec; dARatiodtVec];

%-------------------------------------------------

%ode45 time stepping for Finite Volume
%the conservation variable is q=J*c*A, where J=dx/dz
function LocalOde45FV(y0,tspan,ARGS)
global Nspecies betaDispersion bPressure hChannel Ngrid PolDeg EquiPLength
global PlotTime errVec PrepareGridStage AChannel met2lit
global uCST
global Kw

[cSize,DeltatNextPlot,DeltaCounterNextPlot,dtInit,dtMin,dtMax,rtol, ...
    atol,normcontrol,RK_A,RK_B,RK_E,RK_pow,uCST,A1,B1,A2,B2,M,Cur,INP, ...
    AdaptGrid,RKOrder,SteadyStateFlag,L1,L2,zetaPot,hChannel,bPressure,DeltaP,...
    betaDispersion,Res,DisableChemFlag,SpatialDiscFlag,IonicEffectFlag, ...
    PcentIonicChange, handles, dz]=ARGS{:};

N=cSize(1);

% Additional parameters
qMat=reshape(y0(1:cSize(1)*cSize(2)),cSize);
%phi refers to the physical domain
phiVec=y0(cSize(1)*cSize(2)+1:cSize(1)*cSize(2)+N);
AreaRatioCellVec=y0(cSize(1)*cSize(2)+N+1:end);

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

AreaRatioCellMat=AreaRatioCellVec*ones(1,Nspecies); %area ratio at cell centers
cMat = qMat./AreaRatioCellMat; %initially dx/dz=1

% Create Polynomials for Equilibrium reactions
[PCube,PPrimeCube,QMat,QPrimeMat,LCube,ValCube,muCube,DCube,zListArranged, KaListCube] ...
    =EquilibriumPolynomials(INP);

FirstTimeFlag=1; %initialy it is set = 1 for first run
CalcChemEqFlag=1;
cH=zeros(1,size(cMat,1));
[cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
    cH,PCube,QMat,PPrimeCube,QPrimeMat,LCube,cMat);


[cizCube,cH,cHCubePower,gizCube, IonicStrength, LCube, Kw_new]= ...
    CalculateIonicEffects(IonicEffectFlag, cH,LCube,cMat/met2lit, ValCube, cizCube, zListArranged, KaListCube, cHCubePower, gizCube);
cH=zeros(1,size(qMat,1)); %set initial [H+] = 0

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

        [cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
            cH,PCube,QMat,PPrimeCube,QPrimeMat,LCube,cMat);

        %Define arguments to call the CalcFlux function
        ARGS={A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,...
            PolDeg,EquiPLength,LCube, Kw_new, muCube,ValCube, DCube, N, Nspecies, ...
            zetaPot,h,DisableChemFlag,SpatialDiscFlag, zListArranged, KaListCube};

        IonicCalcFlag=1*IonicEffectFlag; %compute only IonicEffectFlag=1 :)

        %call CalcFlux function
        [f0,cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec, ...
    IonicStrength, LCube, Kw_new]= ...
            CalcFluxFVSlip(y0,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
            IonicCalcFlag, IonicStrength, IonicEffectFlag, dz);

        IonicStrengthOld=IonicStrength;

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

    IonicCalcFlag=0;
    %NOTE: Chemical equilibrium is not computed in intermediate steps
    % LOOP FOR ADVANCING ONE STEP.
    nofailed = true;                      % no failed attempts
    while true


        %Define arguments to call the CalcFlux function
        ARGS={A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,...
            PolDeg,EquiPLength,LCube, Kw_new, muCube,ValCube, DCube, N, Nspecies, ...
            zetaPot,h,DisableChemFlag,SpatialDiscFlag,  zListArranged, KaListCube};

        [f(:,2),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec, ...
    IonicStrength, LCube, Kw_new]= ...
            CalcFluxFVSlip(y+h.*(f*RK_B(:,1)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
            IonicCalcFlag, IonicStrength, IonicEffectFlag, dz);

    [f(:,3),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec, ...
    IonicStrength, LCube, Kw_new]= ...
            CalcFluxFVSlip(y+h.*(f*RK_B(:,2)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
            IonicCalcFlag, IonicStrength, IonicEffectFlag, dz);


        if RKOrder==45
          [f(:,4),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec, ...
                IonicStrength, LCube, Kw_new]= ...
                CalcFluxFVSlip(y+h.*(f*RK_B(:,3)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
                IonicCalcFlag, IonicStrength, IonicEffectFlag, dz);

           [f(:,5),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec, ...
                IonicStrength, LCube, Kw_new]= ...
                CalcFluxFVSlip(y+h.*(f*RK_B(:,4)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
                IonicCalcFlag, IonicStrength, IonicEffectFlag, dz);

            [f(:,6),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec, ...
                IonicStrength, LCube, Kw_new]= ...
                CalcFluxFVSlip(y+h.*(f*RK_B(:,5)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
                IonicCalcFlag, IonicStrength, IonicEffectFlag, dz);
        end

        if     RKOrder==45,  tnew = t + min(h*RK_A(6));
        elseif RKOrder==23,  tnew = t + min(h*RK_A(3));
        end

        if RKOrder==45,
            deltay=h.*(f*RK_B(:,6));
            ynew = y + deltay;

              [f(:,7),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec, ...
    IonicStrength, LCube, Kw_new]= ...
                CalcFluxFVSlip(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
                IonicCalcFlag, IonicStrength, IonicEffectFlag, dz);




        elseif RKOrder==23,
            deltay=h.*(f*RK_B(:,3));
            ynew = y + deltay;
               [f(:,4),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec, ...
    IonicStrength, LCube, Kw_new]= ...
                CalcFluxFVSlip(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
                IonicCalcFlag, IonicStrength, IonicEffectFlag, dz);
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


    if PrepareGridStage==0 && max(abs(IonicStrength-IonicStrengthOld))>...
            1.0e-2*PcentIonicChange*max(IonicStrengthOld)
        IonicCalcFlag=1*IonicEffectFlag; %if IonicEffectFlag=0 then no need to calculate IonicEquilibrium
    end


    CalcChemEqFlag=1;
    [f(:,1),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec, dphidzAreaMat]= ...
                CalcFluxFVSlip(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, ...
                IonicCalcFlag, IonicStrength, IonicEffectFlag, dz);

            if IonicCalcFlag==1 && CalcChemEqFlag==1
                %plot(IonicStrength-IonicStrengthOld, 'r-');
                %legend(['ode counter =', num2str(ODEcounter)]);
                IonicStrengthOld=IonicStrength;
            end


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

        keyIn = get(gcf, 'CurrentCharacter');
        if (Res(end,end)/Res(1,end)<0.2 && PrepareGridStage) || strcmp(keyIn,'m')
            %tic
            disp('Migration begins');
            PrepareGridStage=0;
            IonicCalcFlag=1*IonicEffectFlag;
            set(gcf, 'CurrentCharacter','t');
            clear Res
        end

        set(handles.DataHolder,'UserData',[phiVec,cMat]);
        if IonicEffectFlag==1
           pH=-log10(cH.*10.^(-0.50850*(sqrt(IonicStrength)./(1+sqrt(IonicStrength))-0.3*IonicStrength)));
           pHVecAllTimes(:,PlotCounter)=pH;
        else
           pH=-log10(cH);
           pHVecAllTimes(:,PlotCounter)=pH;
        end

        PlotCurrentSelection(handles,phiVec,cMat,muMat,DMat,dPotendxMat,pH,t,h,err,SigVec,dphidzVec, AreaRatioCellVec);

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
function  [dydt,cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
    IonicStrength, LCube, Kw_new]= ...
    CalcFlux(y0,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, IonicCalcFlag, ...
    IonicStrength, IonicEffectFlag)
global PrepareGridStage Rmu Temp F met2lit
global muH muOH Kw visc hChannel bPressure DeltaP betaDispersion
global uCST
global Nspecies PolDeg EquiPLength

[A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,PolDeg, ...
    EquiPLength,LCube, Kw_new, muCube,ValCube,DCube,N,Nspecies,zetaPot,h, ...
    DisableChemFlag,SpatialDiscFlag, zListArranged, KaListCube]=ARGS{:};

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
    ValTmp=sum(ValCube,3)';


    %  cMat(:,end)=-sum(ValTmp(:,1:end-1).*cMat(:,1:end-1),2)./ValTmp(:,end);
else
    if CalcChemEqFlag
        [cizCube,cH,cHCubePower,gizCube]=LzCalcEquilibrium(cH,LCube,cMat, ValCube, Kw_new);

        if (IonicCalcFlag==1)
            %%%%must call CalculateEquilibrium before CalculateIonicEffects
            [cizCube,cH,cHCubePower,gizCube, IonicStrength, LCube, Kw_new]= ...
                CalculateIonicEffects(IonicEffectFlag, cH,LCube,cMat, ValCube, cizCube, zListArranged, KaListCube, cHCubePower, gizCube);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        IonicStrength=0.5*sum(sum(ValCube.^2.*cizCube,3),1)/met2lit+0.5*(cH+Kw_new./cH)*IonicEffectFlag;
        %ionic strength includes [H+] and [OH-]
        %%%Now that we have Ionic strength, we can calculate the new MobilityCube
        ICube=repmat(IonicStrength, [Nspecies,1,PolDeg]);
        muIonicCube = OnsagerFuoss(IonicEffectFlag, ValCube, cizCube, muCube);
        muMat=sum(ValCube.*muIonicCube.*gizCube,3)';
        alphaMat=F*sum(ValCube.^2.*muIonicCube.*gizCube,3)';
        DTemp=DCube.*(1-0.2297*sqrt(ICube));
        DMat=sum(DTemp.*gizCube,3)';
        betaMat=F*sum(ValCube.*DTemp.*gizCube,3)';
        %%% CalculateIonicEffects.
    end % if CalcChemEqFlag
end

% Calculate velocity (pressure driven)
L=phiVec(end)-phiVec(1);
uBulk = hChannel^2/(bPressure*visc*L)*(-DeltaP);

if CalcChemEqFlag
    % Dispersion
    PeMat=uBulk*(hChannel)./DMat;                   %h=diameter  h/2=radius
    DMat=DMat.*(1+betaDispersion*PeMat.^2);
end % if CalcChemEqFlag
if DisableChemFlag
    SigVec = sum(alphaMat'.*cMat',1)';
    SVec   = sum(betaMat'.*cMat',1)';
else
    %muH_new=muH-(0.2297*muH+31.410e-9)*sqrt(IonicStrength)./(1+0.329*4.0*sqrt(
    %IonicStrength));
    %muOH_new=muOH-(0.2297*muOH+31.410e-9)*sqrt(IonicStrength)./(1+0.329*4.0*sqrt(IonicStrength));
    muH_new=muH;
    muOH_new=muOH;
    %change the Kw since ionic effects are to be considered
    %Kw_new=Kw*10.^(2*0.50850*(sqrt(IonicStrength)./(1+sqrt(IonicStrength))-0.3*IonicStrength));

    Kw_new=Kw*10.^(-2*getLogActivity(IonicStrength));
    SigVec = sum(alphaMat'.*cMat',1)' + F^2*(muH_new.*cH*met2lit+muOH_new.*(Kw./cH)*met2lit)';
    SVec   = sum(betaMat'.*cMat',1)'  + Rmu*Temp*(muH_new.*cH*met2lit-muOH_new.*(Kw./cH)*met2lit)';
end


%calculate uBulk
%OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[cH',phiVec,cMat],N); for derivative
%dcHdz=OUT(:,1)
%constant volatge formulation


%constant voltage formulation
%calculate dS/dx

% % %comment for constant current formulation
% % Voltage=-1200;
% % OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[SVec,phiVec],N);
% % dSVecdz=OUT(:,1);   dphidzVec=OUT(:,2);
% % dSdxVec=dSVecdz./dphidzVec;
% % Den1=trapz(phiVec,1./SigVec) + 10.0e-3/0.025; %+L_right/sigma_{right}+L_left/sigma_{left}
% % Num1=trapz(phiVec,dSdxVec./SigVec);
% % Cur=-(Voltage+Num1)./Den1;
% % %constant voltage formulation ends
% % %fprintf('current = %g\n', Cur);

if SteadyStateFlag %&& FirstTimeFlag   % Automatic velocity, baed on LE (first chemical on list)
    uCST = -Cur/SigVec(end)*muMat(end,1)-uBulk;
end

% Derivatives requried for adaptive grid  (but also for next steps)
OUT=LinSolve(A1.a,A1.b,A1.c,B1.full*[cH',phiVec,cMat],N);
dcHdz=OUT(:,1);   dphidzVec=OUT(:,2);  dcdzMat=OUT(:,3:end);
dphidzMat=dphidzVec*ones(Nspecies,1)';

if AdaptGrid.Coeff~=0
    mu_c_Char=max(max(abs(muMat.*cMat)));  % characteristic mobility*concentration
    L=phiVec(end)-phiVec(1);
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
        dcdtMat = d2Dcdx2Mat + dfdx_UW_Mat + dphidtMat.*dcdxMat;

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
V1L=muMat(IL,:).*dPotendxMat(IL,:)-0*(uBulk+uCST);
V2L=Cur./SigVec(IL)^2.*(cMat(IL,:).*muMat(IL,:));
V3L=alphaMat(IL,:)';

AL=diag(V1L)+V2L'*V3L';      [VAL,DAL] = eig(AL);
SEigAL=sign(real(diag(DAL)));
dRLdt=VAL\dcdtMat(IL,:)';
dRLdt(SEigAL<0)=0;
dcdtMat(1,:)=real((VAL*dRLdt)') +  (dphidtMat(IL,:)-uBulk-uCST).*dcdxMat(IL,:);

%     Right boundary
IR=N;
V1R=muMat(IR,:).*dPotendxMat(IR,:)-(uBulk+uCST);
V2R=Cur./SigVec(IR)^2.*(cMat(IR,:).*muMat(IR,:));
V3R=alphaMat(IR,:)';

AR=diag(V1R)+V2R'*V3R';
[VAR,DAR] = eig(AR);
SEigAR=sign(real(diag(DAR)));
dRRdt=VAR\dcdtMat(IR,:)';
dRRdt(SEigAR>0)=0;
dcdtMat(end,:)=real((VAR*dRRdt)');

dydt=[reshape(dcdtMat,numel(dcdtMat),1);dphidtVec];

%%original.. uses matrix operations and solves for the cH using vector

%--------------------------------------------------------------------------
%  Calculate Chemical Equilibrium
%  Using the output of EquilibriumPolynomials, i.e. P, Q etc this function
%  a) computes C_H = [H+].
%  b) knowing C_H and LMat, it computes g_iz
%  Returns c_iz, c_h, g_iz
%--------------------------------------------------------------------------
function [cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
    cH,PCube,QMat,PPrimeCube,QPrimeMat,LCube,cMat)
global met2lit Ngrid Nspecies EquiPLength PolDeg

% EQUILIBRIUM
if FirstTimeFlag==1  % Solve for pH the first time
    %----------------------------------------------------------------------
    %----
    % Setup constant matrices -  THIS PART RUNS ONLY ONCE PER SIMULATION
    %--------------------------------------------------------------------------
    % Get equilibrium polynomials

    for ij=1:size(cMat,1) %iterate over all grids
        cTot=cMat(ij,:)/met2lit;   % convert back to mol/liter
        cTotRep=cTot'*ones(1,size(PCube(:,:,ij),2));
        P=sum(cTotRep.*PCube(:,:,ij),1);

        Poly=zeros(1,EquiPLength);

        Poly(1:size(P,2))=Poly(1:size(P,2))+P;
        Poly(1:size(QMat,2))=Poly(1:size(QMat,2))+QMat(ij,:); %from QMat
        Poly=fliplr(Poly);

        roo=roots(Poly);
        roo=roo(imag(roo)==0);
        cH(ij)=roo(roo>0);          % NOTE THAT cH is in mol/lit
    end %for ij
else
    % Solve for pH using Newton Raphson
    count=0;        cHPrev=ones(1,Ngrid);
    %%%            while max(abs(cHPrev-cH)) > max([abs(cH),abs(cHPrev)])*1E-3
    %     while norm((cHPrev-cH)./max([abs(cH);abs(cHPrev)],[],1),inf) > 1E-6
    %         count=count+1;
    %         cHPrev=cH;
    %         EquicHPolMat=[ones(1,Ngrid);cumprod(ones(EquiPLength-1,1)*cH,1)];
    %         %% this makes a makes a matrix of powers of cH which is multiplied
    %         %by the companian matrix
    %         fcH=sum(cMat'.*(PMat*EquicHPolMat),1)/met2lit+Q*EquicHPolMat;
    %         fPrimecH=sum(cMat'.*(PPrimeMat*EquicHPolMat),1)/met2lit+QPrime*EquicHPolMat;
    %         cH=cH-fcH./fPrimecH;
    %
    %         if count>5
    %             count;
    %         end
    %         if count>100
    %             disp('Too many iterations on cH. Returning');
    %             return;
    %         end
    %
    %     end
    %can speed up, this using multiprod

    fcH=zeros(1,Ngrid);
    fPrimecH=zeros(1,Ngrid);
    while norm((cHPrev-cH)./max([abs(cH);abs(cHPrev)],[],1),inf) > 1E-6
        count=count+1;
        cHPrev=cH;
        %EquicHPolMat=[ones(1,Ngrid);cumprod(ones(EquiPLength-1,1)*cH,1)];
        for k=1:Ngrid
            %        EquicHPolMat=[1;cumprod(ones(EquiPLength-1,1)*cH(k),1)];
            EquicHPolMat=cH(k).^[0:1:EquiPLength-1]'; %[1;cumprod(ones(EquiPLength-1,1)*cH(k),1)];
            %% this makes a makes a matrix of powers of cH which is multiplied
            %by the companian matrix
            fcH(1,k)=sum(cMat(k,:)'.*(PCube(:,:,k)*EquicHPolMat),1)/met2lit+QMat(k,:)*EquicHPolMat;
            fPrimecH(1,k)=sum(cMat(k,:)'.*(PPrimeCube(:,:,k)*EquicHPolMat),1)/met2lit+QPrimeMat(k,:)*EquicHPolMat;
        end

        cH=cH-fcH./fPrimecH;


        if count>5
            count;
        end
        if count>1000
            disp('Too many iterations on cH. Returning');
            return;
        end

    end



    %===testing multiprod for equilibrium calculation

    %       while norm((cHPrev-cH)./max([abs(cH);abs(cHPrev)],[],1),inf) > 1E-6
    %         count=count+1;
    %         cHPrev=cH;
    %         EquicHPolMat=[ones(1,Ngrid);cumprod(ones(EquiPLength-1,1)*cH,1)];
    %         TempPEMult1=multiprod(PCube, EquicHPolMat, [1 2], [1]);
    %         TempPEMult2=multiprod(PPrimeCube, EquicHPolMat, [1 2], [1]);
    %                 %% this makes a makes a matrix of powers of cH which is multiplied
    %         %by the companian matrix
    %         for k=1:Ngrid
    %         fcH=sum(cMat(k,:)'.*TempPEMult1(:,:,k),1)/met2lit +QMat(k,:)*EquicHPolMat;
    %         fPrimecH=sum(cMat(k,:)'.*TempPEMult2(:,:,k),1)/met2lit +QPrimeMat(k,:)*EquicHPolMat;
    %         cH(k)=cH(k)-fcH/fPrimecH;
    %         end
    %         %fcH=fcH+QMat*EquicHPolMat;
    %         %fPrimecH=fPrimecH+QPrimeMat*EquicHPolMat;
    %         %cH=cH-fcH./fPrimecH;
    %    end
    %         if count>5
    %             count;
    %         end
    %         if count>1000
    %             disp('Too many iterations on cH. Returning');
    %             return;
    %         end

    %======= testing for mulitprod for equilibrium calculation ends here



end %if first time
%cH=cH*met2lit;
cHPolMat=[ones(1,Ngrid);cumprod(ones(PolDeg-1,1)*cH,1)];

Temp=zeros(Nspecies,Ngrid);
TempL=permute(LCube,[1,3,2]); %8x2x150
Temp=multiprod(TempL, cHPolMat, [1 2], [1]); %faster multiplication
%========testing
%for k=1:1:Ngrid
%    Temp(:,k)=TempL(:,:,k)*cHPolMat(:,k);
%end
%=====testing

M1Cube=FastRepmatPages(cMat'./Temp,PolDeg);
% M1Cube=FastRepmatPages(cMat'./(sparse(LMat)*cHPolMat),PolDeg);
cHCubePower=FastRepmatRows(permute(cHPolMat,[3,2,1]),Nspecies);
cizCube=LCube.*cHCubePower.*M1Cube; %equation 5 in J. Chroma A
% %LCube should be properly computed
gizCube=LCube.*cHCubePower.*FastRepmatPages(1./Temp,PolDeg);

% THIS FUNCTION IS REDUNDANT AND IS NOT USED ANYWHERE
%
%
%--------------------------------------------------------------------------
%  Calculate Chemical Equilibrium when Ionic Effects are ON
%  Using the output of EquilibriumPolynomials, i.e. P, Q etc this function
%  a) computes C_H = [H+].
%  b) knowing C_H and LMat, it computes g_iz
%  Returns c_iz, c_h, g_iz
%--------------------------------------------------------------------------
function [cizCube,cH,cHCubePower,gizCube]=CalculateEquilibriumIonic(FirstTimeFlag,...
    cH,PMat,Q,PPrimeMat,QPrime,PolDeg,EquiPLength,N,LMat,LCube,cMat,Nspecies)
global met2lit Ngrid

% EQUILIBRIUM
if FirstTimeFlag==1  % Solve for pH the first time
    %----------------------------------------------------------------------
    %----
    % Setup constant matrices -  THIS PART RUNS ONLY ONCE PER SIMULATION
    %--------------------------------------------------------------------------
    % Get equilibrium polynomials

    for ij=1:Ngrid %iterate over all grids
        cTot=cMat(ij,:)/met2lit;   % convert back to mol/liter
        cTotRep=cTot'*ones(1,size(PMat,2)); %probably put PCube(ij)
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
    for k=1:N
        while abs(cHPrev(k)-cH(k))/max(abs(cH(k)),abs(cHPrev(k))) > 1E-6
            count=count+1;
            cHPrev(k)=cH(k);
            EquicHPolMat=[1;cumprod(ones(EquiPLength-1,1)*cH(k),1)];
            %% this makes a makes a matrix of powers of cH which is multiplied
            %by the companian matrix
            fcH=sum(cMat(k,:)'.*(PMat*EquicHPolMat),1)/met2lit+Q*EquicHPolMat;
            fPrimecH=sum(cMat(k,:)'.*(PPrimeMat*EquicHPolMat),1)/met2lit+QPrime*EquicHPolMat;
            cH(k)=cH(k)-fcH/fPrimecH;

            if count>5
                count;
            end
            if count>100
                disp('Too many iterations on cH. Returning');
                return;
            end

        end
    end
end %if first time

cHPolMat=[ones(1,N);cumprod(ones(PolDeg-1,1)*cH,1)]; %will have to make again
M1Cube=FastRepmatPages(cMat'./(sparse(LMat)*cHPolMat),PolDeg);
cHCubePower=FastRepmatRows(permute(cHPolMat,[3,2,1]),Nspecies);
cizCube=LCube.*cHCubePower.*M1Cube;

gizCube=LCube.*cHCubePower.*FastRepmatPages(1./(sparse(LMat)*cHPolMat),PolDeg);

%----------------------------------------------------------------------
% Calculates Ionic Effects.
% a) Returns same values as CalculateEquilibrium()
% b) Inputs: ARG, KaList, zArrangedList
% c) Called within CalcFLux()
% d) This function calls CalculateEquilibrium() and
% RecomputeEquilibriumPolynomials() to iterate for ionic effects

%call CalculateEquilibrium() prior to call this


%updates: Feb 27, 2009
%corrected bug in definition of NewKaListCube{k}{i}, pKaFactor
%----------------------------------------------------------------------
function [cizCube,cH,cHCubePower,gizCube, IonicStrength, LCube, Kw_new]= ...
    CalculateIonicEffects(IonicEffectFlag, cH,LCube, cMat, ValCube, cizCube, zListArranged, KaListCube, cHCubePower, gizCube)
global Ngrid Nspecies met2lit Kw PolDeg

Kw_new=Kw*ones(1,Ngrid);
if (IonicEffectFlag==0)
    Temp=sum((LCube.*cHCubePower),3);
    gizCube=LCube.*cHCubePower.*FastRepmatPages(1./Temp,PolDeg);
    %I know this is redundant.. let it be here for a while
    IonicStrength=zeros(1,Ngrid);
    return;
end



% see how zListArranged, KaList can be taken as input in CalcFlux(), as this
% function is called in CalcFlux()


convergence=false;
count=0;

while ~convergence

    if(count>40)
        disp('no convergence in Ionic Effects');
        break;
    end
    count=count+1;
    cHPrev=cH;
    %START LOOP
    %step 2: Compute I=0.5 sum_i sum_z c_{iz}*z_{iz}^2
    IonicStrength=0.5*sum(sum(ValCube.^2.*cizCube,3),1)/met2lit+cH+Kw_new./cH;
    %cizCube should be computed properly

    % see if you can iterate over the grid points
    %but will have to changed the way LCube, gizCube etc needs to be
    %computed

    %step 3: Using I, zArrangedList and KaList compute delta_pKa_Cube
    %step 4: NewKaListCube = KaListCube + delta_pKa_Cube
    %pKaFactor = 0.50850*(sqrt(IonicStrength)./(1+sqrt(IonicStrength))-0.3*IonicStrength);
    pKaFactor = -getLogActivity(IonicStrength);
    NewKaListCube=KaListCube; %initialize so that it doesn't grow inside the loop

    for i=1:Nspecies
        zListHeavy{i}=zListArranged{i}-Myheaviside(zListArranged{i}-0.1);
        %precomputing these so that it does not go into loop over grid points
    end


    for k=1:Ngrid
        for i=1:Nspecies
            %     NewKaListCube{k}{i}=NewKaListCube{k}{i}.*10.^(-pKaFactor(k)*2.0*(zListArranged{i}-...
            %         heaviside(zListArranged{i}-0.1)));
            NewKaListCube{k}{i}=NewKaListCube{k}{i}.*10.^(-pKaFactor(k)*2.0*zListHeavy{i});
            %i have taken abs because for base, delta_pKb=delta_pKa
        end
    end
    Kw_new=Kw*10.^(2*pKaFactor); %Kw/gamma^2 this goes into ReComputeEquilibriumPolynomials
    %Kw_new=Kw;
    %step 5: Call RecomputeEquilibriumPolynomials using "NewKaList" to get new LCube
    LCube=RecomputeEquilibriumPolynomials(zListArranged, NewKaListCube);
    %step 6: Call LzCalcEquilibrium()
    [cizCube,cH,cHCubePower,gizCube]=LzCalcEquilibrium(cH,LCube,cMat, ValCube, Kw_new);

    %step 7: Check for convergence
    %convergence of cH values
    if (norm((cHPrev-cH)./max([abs(cH);abs(cHPrev)],[],1),inf) < 1E-6)
        convergence=true;
        break;
    end



    %GO UP, if NO congergence, else END of LOOK. Return Values
end
%  disp('iteration took counts = ');
%  disp(count);

%step 8: After convergence, compute new IonicStrength
IonicStrength=0.5*sum(sum(ValCube.^2.*cizCube,3),1)/met2lit+0.5*(cH+Kw_new./cH);
%step 9: Compute New muCube and using that compute DCube

%updates: Feb 27, 2009: coorected value for muMat
%Kw_new is now computed based on the ionic strength


function [pH,SigVec,muMat,DMat] = CalculateSpatialProperties(INP,N,FirstTimeFlag,cMat, IonicEffectFlag)
global F Rmu  Temp Kw muH muOH met2lit Nspecies PolDeg Ngrid EquiPLength

%cMat=cMat/met2lit;
Nspecies=size(INP,1);
%[PMat,PPrimeMat,Q,QPrime,LMat,ValMat,muMat,KaMat,DMat]=EquilibriumPolynomials(INP);

Ngrid=size(cMat,1); %grid size, global variable

[PCube,PPrimeCube,QMat,QPrimeMat,LCube,ValCube,muCube,DCube,...
    zListArranged, KaListCube]=EquilibriumPolynomials(INP);

EquiPLength=max(size(PCube,2),size(QMat,2));
PolDeg=size(LCube,3); % Polynomial degree

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cH=zeros(1,size(cMat,1));

FirstTimeFlag=1;
[cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
    cH,PCube,QMat,PPrimeCube,QPrimeMat,LCube,cMat);



[cizCube,cH,cHCubePower,gizCube, IonicStrength, LCube, Kw_new]= ...
    CalculateIonicEffects(IonicEffectFlag, cH,LCube,cMat, ValCube, cizCube, zListArranged, KaListCube, cHCubePower, gizCube);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Now that we have Ionic strength, we can calculate the new MobilityCube
IonicStrength=IonicStrength*IonicEffectFlag;
ICube=repmat(IonicStrength, [Nspecies,1,PolDeg]);
%ICube will be zero if IonicEffectFlag=0
%in that case mu D etc will be caluclated without Ionic Strength dependence
%original
%muMat=sum(gizCube.*(F*ValCube.*muCube.*(1-0.2297*sqrt(ICube)./(1+0.329*4.0*sqrt(ICube))) ...
%        -(sign(ValCube)).*31.410e-9.*sqrt(ICube)./(1+0.329*4.0*sqrt(ICube))),3)';

% muMat=sum(gizCube.*(F*ValCube.*muCube.*(1-0.2297*abs(ValCube).*sqrt(ICube)./(1+1.5*sqrt(ICube))) ...
%       -(ValCube).*31.410e-9.*sqrt(ICube)./(1+1.5*sqrt(ICube))),3)';

muIonicCube = OnsagerFuoss(IonicEffectFlag, ValCube, cizCube, muCube);
muMat=sum(ValCube.*muIonicCube.*gizCube,3)';
%%%use correct units
%%%I have used value 4A = atomic radii, but actually its depends on the
%%%problem

DTemp=DCube.*(1-0.2297*sqrt(ICube));

%muMat=F*sum(ValCube.*muCube.*gizCube,3)';
%DMat=sum(DCube.*gizCube,3)'; %original without ionic strength
DMat=sum(DTemp.*gizCube,3)';
%alphaMat=F^2*sum(ValCube.^2.*muCube.*gizCube,3)';
alphaMat=F*sum(ValCube.^2.*muIonicCube.*gizCube,3)';
%  alphaMat=F*sum(ValCube.*(gizCube.*(F*ValCube.*muCube.*(1-0.2297*abs(ValCube).*sqrt(ICube)./(1+1.5*sqrt(ICube))) ...
%        -(ValCube).*31.410e-9.*sqrt(ICube)./(1+1.5*sqrt(ICube)))),3)';


%betaMat=F*sum(ValCube.*DCube.*gizCube,3)'; %original without ionic strength
%betaMat=F*sum(ValCube.*DTemp.*gizCube,3)';

DTemp=DCube.*(1-0.2297*sqrt(ICube));
DMat=sum(DTemp.*gizCube,3)';
%        betaMat=F*sum(ValCube.*DCube.*gizCube,3)';
betaMat=F*sum(ValCube.*DTemp.*gizCube,3)';



muH_new=muH;
muOH_new=muOH;
Kw_new=Kw*10.^(-2*getLogActivity(IonicStrength*IonicEffectFlag));

%compute new Kw because of change in ionic strength
SVec   = sum(betaMat'.*cMat',1)'  + Rmu*Temp*(muH_new.*cH*met2lit-muOH_new.*(Kw_new./cH)*met2lit)';
SigVec = sum(alphaMat'.*cMat',1)' + F^2*(muH_new.*cH*met2lit+muOH_new.*(Kw_new./cH)*met2lit)';

if IonicEffectFlag==1
    pH=-log10(cH.*10.^getLogActivity(IonicStrength));
else
    pH=-log10(cH);
end


%helper function in plotting data
function ChangePlotVariable(handles)
global met2lit

%disp('ChangePlotVariable was called');

N = str2num (get(handles.GridPointsEdit,'string'));
InputTable=cell(handles.InputTable.getData);
InputTable=InputTable(:,2:end);

for ij=1:size(InputTable,1)
    INP{ij,1}=str2num(InputTable{ij,2});
end %for ij

FirstTimeFlag=1;
%IonicEffectFlag=1; %later get this from handles
IonicEffectFlag = 2-get(handles.IonicEffectPopup,'value');
AxesData=get(handles.DataHolder,'UserData');
if isempty(AxesData)
 L = str2num (get(handles.DomainLengthEdit,'string'));
 xVec = linspace(0,L*1e-3,N)';
 cMat = xVec + NaN; pH = xVec + NaN; muMat = xVec + NaN; DMat = xVec + NaN;
 SigVec = xVec + nan;
 else
xVec=AxesData(:,1); cMat=AxesData(:,2:end); %earlier *met2lit;
[pH,SigVec,muMat,DMat] = CalculateSpatialProperties(INP,N,FirstTimeFlag,cMat, IonicEffectFlag);
end

%disp('ChangePlotVariable->ChangeSpatialProperties ended');
t=NaN; h=NaN; err=NaN; dphidzVec=xVec+nan; dPotendxMat=xVec+nan;

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


PlotCurrentSelection(handles,xVec,cMat,muMat,DMat,dPotendxMat,pH,t,h,err,SigVec,dphidzVec, AreaRatio);



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
function [PCube,PPrimeCube,QMat,QPrimeMat,LCube,ValCube,muCube,DCube, ...
    zListArranged, KaListCube]=EquilibriumPolynomials(INP)
global F Rmu Temp Kw Ngrid EquiPLength PolDeg

%Ngrid=150; %testing

% PREPARE MATRICES
%------------------
MaxCol=-Inf;
for j=1:size(INP,1) %iterate on species
    MaxCol=max([MaxCol,max(INP{j}(1:3:end))-min(INP{j}(1:3:end))+1]);
end

LMat=zeros(size(INP,1),MaxCol);

ValMat=zeros(size(INP,1),MaxCol);
muMat=ValMat; KaMat=ValMat; DMat=ValMat;

for j=1:size(INP,1) % for every specie
    zList=INP{j}(1:3:end);
    muList=INP{j}(2:3:end)./(F*abs(zList));  % Defines Santiago mobility
    pKaList=INP{j}(3:3:end);
    KaList=10.^(-pKaList);

    %%make 2-D-Kalist
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


    ValMat(j,1:length(zList))=zList;
    muMat(j,1:length(muList))=muList;
    KaMat(j,1:length(KaList))=KaList;
    DMat(j,1:length(DList))=DList;

    zListArranged{j}=zList;

    for k=1:Ngrid
        KaListCube{k}{j}=KaList; %k=grid point, j=species
    end

    nj=min(zList);    pj=max(zList);


    for z=zList
        if z<0
            %LMat is LMat
            LMat(j,z-nj+1)=prod(KaList(z-nj+1:-nj));
        elseif z>0
            LMat(j,z-nj+1)=1/prod(KaList(-nj+2:z-nj+1));
        elseif z==0
            LMat(j,z-nj+1)=1;
        end %if
    end % for z

end %for ij

% CONSTRUCT POLYNOMIALS
%--------------------
Q1=1;
for j=1:size(LMat,1)
    Q1=conv(Q1,LMat(j,:));
end %for j
Q2=[-Kw 0 1];
Q=conv(Q1,Q2);

for i=1:size(INP,1)
    tmp=zeros(1,size(LMat,2));
    tmp(1:length(zListArranged{i}))=zListArranged{i};
    Mmod=LMat;     Mmod(i,:)=Mmod(i,:).*tmp;

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

% PMat=sparse(PMat); %these matrices are not sparse!
% Q=sparse(Q);

EquiPLength=max(size(PMat,2),size(Q,2)); %make these using Pcube
PolDeg=size(LMat,2); % Polynomial degree %make this using LCube

%Cube here refers to 3-d array
%cols: grid points, rows: no. of species depth: polynomial degree
muCube=repmat(reshape(muMat,[size(INP,1),1,PolDeg]),[1,Ngrid,1]); %mobilities
DCube=repmat(reshape(DMat,[size(INP,1),1,PolDeg]),[1,Ngrid,1]); %diffusivities
ValCube=repmat(reshape(ValMat,[size(INP,1),1,PolDeg]),[1,Ngrid,1]); %valence
LCube=FastRepmatColumns(reshape(LMat,[size(INP,1),1,PolDeg]),Ngrid); %L
PCube=FastRepmatColumns(reshape(PMat,[size(INP,1),1,EquiPLength]),Ngrid);
PPrimeCube=FastRepmatColumns(reshape(PPrimeMat,[size(INP,1),1, EquiPLength]),Ngrid); %L
PCube=permute(PCube,[1, 3, 2]);
PPrimeCube=permute(PPrimeCube,[1, 3, 2]);

QMat=repmat(Q,Ngrid,1);
QPrimeMat=repmat(QPrime,Ngrid,1);
%First time LCube is just made by repeating LMat as pKas are same on all
%grid points. But later for Ionic strengths it should be made differently

% ---------------------------------

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

function y=getLogActivity(IonicStrength)
%y=-0.50850*(sqrt(IonicStrength)./(1+sqrt(IonicStrength))-0.3*IonicStrength);
%used in peakmaster.. check once again
%lc
y=-0.50850*(sqrt(IonicStrength)./(1+1.5*sqrt(IonicStrength))-0.1/0.5085*IonicStrength);

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

    case ('ELED')
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

if (min(b==0) && min(c==0) && min(a==1))
    OUT=RHS;
else
    A=diag(a,0)+diag(b,-1)+diag(c,1);
    RHS=A\RHS;
end

% Time stepping based on ode45
%probably make another function LocalOde23, and remove all if conditions
%for RK45 or RK23
function LocalOde45(y0,tspan,ARGS)
global Nspecies betaDispersion bPressure hChannel Ngrid PolDeg EquiPLength
global PlotTime errVec PrepareGridStage met2lit AChannel
global uCST
global Kw
[cSize,DeltatNextPlot,DeltaCounterNextPlot,dtInit,dtMin,dtMax,rtol,atol,...
    normcontrol,RK_A,RK_B,RK_E,RK_pow,uCST,A1,B1,A2,B2,M,Cur,INP,...
    AdaptGrid,RKOrder,SteadyStateFlag,L1,L2,zetaPot,hChannel,bPressure, ...
    DeltaP,betaDispersion,Res,DisableChemFlag,SpatialDiscFlag,IonicEffectFlag, PcentIonicChange, handles]=ARGS{:};
N=cSize(1);
Ngrid=N; %global variable

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

%  figure(1);
%  set(0,'CurrentFigure',1);


% Create Polynomials for Equilibrium reactions
[PCube,PPrimeCube,QMat,QPrimeMat,LCube,ValCube,muCube,DCube,zListArranged, KaListCube] ...
    =EquilibriumPolynomials(INP);

FirstTimeFlag=1; %initialy it is set = 1 for first run
CalcChemEqFlag=1;
cH=zeros(1,size(cMat,1));
[cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
    cH,PCube,QMat,PPrimeCube,QPrimeMat,LCube,cMat);


[cizCube,cH,cHCubePower,gizCube, IonicStrength, LCube, Kw_new]= ...
    CalculateIonicEffects(IonicEffectFlag, cH,LCube,cMat/met2lit, ValCube, cizCube, zListArranged, KaListCube, cHCubePower, gizCube);


% % %Cube here refers to 3-d array
% % %cols: grid points, rows: no. of species depth: polynomial degree
% % muCube=repmat(reshape(muMat,[Nspecies,1,PolDeg]),[1,N,1]); %mobilities
% % DCube=repmat(reshape(DMat,[Nspecies,1,PolDeg]),[1,N,1]); %diffusivities
% % ValCube=repmat(reshape(ValMat,[Nspecies,1,PolDeg]),[1,N,1]); %valence
% % LCube=FastRepmatColumns(reshape(LMat,[Nspecies,1,PolDeg]),N); %L,
% % %First time LCube is just made by repeating LMat as pKas are same on all
% % %grid points. But later for Ionic strengths it should be made differently

cH=zeros(1,size(cMat,1)); %set initial [H+] = 0

%---------------------
% Make a datastructure containing KaListComplete for all grid points
% That would be used in the ARGS to be used in CalculateEquilibrium and
% CalcFlux and CaculateEquilibriumIonic and CalculateIonicEffects
%---------------------

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

        %         [cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
        %     cH,PCube,QMat,PPrimeCube,QPrimeMat,LCube,cMat)

        %Calculate Equilibrium using Moran's function for first time flag only
        %Otherwise in CalcFlux Bahga's equilibrium function is used
        [cizCube,cH,cHCubePower,gizCube]=CalculateEquilibrium(FirstTimeFlag,...
            cH,PCube,QMat,PPrimeCube,QPrimeMat,LCube,cMat);

        %Define arguments to call the CalcFlux function
        ARGS={A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,...
            PolDeg,EquiPLength,LCube, Kw_new, muCube,ValCube, DCube, N, Nspecies, ...
            zetaPot,h,DisableChemFlag,SpatialDiscFlag,  zListArranged, KaListCube};

        IonicCalcFlag=1*IonicEffectFlag; %compute only IonicEffectFlag=1 :)

        %call CalcFlux function
        [f0,cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
            IonicStrength, LCube, Kw_new]= ...
            CalcFlux(y0,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag,IonicCalcFlag, IonicStrength, IonicEffectFlag);

        IonicStrengthOld=IonicStrength;

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
    %could probably move this out of the loop.
    %why do you want to check if it is the first time step even after
    %you have computed N steps?????

    h = min(hmax, max(hmin, h));
    if h>(tfinal-t)
        h=tfinal-t;
    end

    IonicCalcFlag=0;
    %NOTE: Chemical equilibrium is not computed in intermediate steps
    % LOOP FOR ADVANCING ONE STEP.
    nofailed = true;                      % no failed attempts
    while true

        ARGS={A1,B1,A2,B2,M,Cur,AdaptGrid,FirstTimeFlag,SteadyStateFlag,...
            PolDeg,EquiPLength,LCube, Kw_new, muCube,ValCube, DCube, N, Nspecies, ...
            zetaPot,h,DisableChemFlag,SpatialDiscFlag,  zListArranged, KaListCube};

        [f(:,2),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
            IonicStrength, LCube, Kw_new]= ...
            CalcFlux(y+h.*(f*RK_B(:,1)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, IonicCalcFlag, IonicStrength, IonicEffectFlag);
        [f(:,3),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
            IonicStrength, LCube, Kw_new]=...
            CalcFlux(y+h.*(f*RK_B(:,2)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, IonicCalcFlag, IonicStrength, IonicEffectFlag);

        if RKOrder==45
            [f(:,4),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
                IonicStrength, LCube, Kw_new]=...
                CalcFlux(y+h.*(f*RK_B(:,3)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, IonicCalcFlag, IonicStrength, IonicEffectFlag);

            [f(:,5),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
                IonicStrength, LCube, Kw_new]=...
                CalcFlux(y+h.*(f*RK_B(:,4)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, IonicCalcFlag, IonicStrength, IonicEffectFlag);

            [f(:,6),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
                IonicStrength, LCube, Kw_new]=...
                CalcFlux(y+h.*(f*RK_B(:,5)),t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, IonicCalcFlag, IonicStrength, IonicEffectFlag);
        end

        if     RKOrder==45,  tnew = t + min(h*RK_A(6));
        elseif RKOrder==23,  tnew = t + min(h*RK_A(3));
        end

        if RKOrder==45,
            deltay=h.*(f*RK_B(:,6));
            ynew = y + deltay;
            [f(:,7),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
                IonicStrength, LCube, Kw_new] = ...
                CalcFlux(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, IonicCalcFlag, IonicStrength, IonicEffectFlag);

        elseif RKOrder==23,
            deltay=h.*(f*RK_B(:,3));
            ynew = y + deltay;
            [f(:,4),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
                IonicStrength, LCube, Kw_new] = ...
                CalcFlux(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, IonicCalcFlag, IonicStrength, IonicEffectFlag);
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
    %     disp('time is = ');
    %     disp(t);
    %     disp('ode counter =');
    %     disp(ODEcounter);

    % This flag = 1 because, now I have completed one full time step.
    % Now to start another time step I need to compute chemical equilibrium
    % In other words, Chemical Equilibrium is computed only at the
    % beginning of time step and NOT for the intermediate time steps of
    % RK45/RK23

    %
    %         if PrepareGridStage==0 && mod(ODEcounter,10)==1
    %             IonicCalcFlag=1;
    %            % plot(muMat);
    %         end
    %


    if PrepareGridStage==0 && max(abs(IonicStrength-IonicStrengthOld))>...
            1.0e-2*PcentIonicChange*max(IonicStrengthOld)
        IonicCalcFlag=1*IonicEffectFlag; %if IonicEffectFlag=0 then no need to calculate IonicEquilibrium
        %disp('Using Ionic Strength');
        %fprintf('max(IonicStrength) = %f\n', max(IonicStrength));
        % plot(muMat);
    end
    %when IonicCalcFlag=1, then new LMat etc are computed. This is passed out
    %through the CalcFlux function.
    %Now, when the IonicCalcFlag=0, we do not compute new pH. But the LMat etc
    %are based on the changes in pKa over last time, which were passed through
    %the CalcFlux


    %    CalcChemEqFlag==1, means that we will compute the chemical
    %    equilibrium using the LMat etc which were updated last time
    %    IonicCalcFlag was 1.
    %    CalcChemEqFlag==1 with IonicCalcFlag==0 means, old values of pKas are
    %    used, but new mobility is used based not newly computed
    %    ChemicalEquilibrium

    %    CalcChemEqFlag==1 with IonicCalcFlag==1 means, new pKas are computed
    %    and new IonicStrength is used to update the mobilities.



    CalcChemEqFlag=1;
    [f(:,1),cH,GridCost,muMat,DMat,alphaMat,betaMat,dPotendxMat,SigVec,dphidzVec,...
        IonicStrength, LCube, Kw_new] = ...
        CalcFlux(ynew,t,cH,ARGS,GridCost,muMat,DMat,alphaMat,betaMat,CalcChemEqFlag, IonicCalcFlag, IonicStrength, IonicEffectFlag);
    %===================testing=============
    if IonicCalcFlag==1 && CalcChemEqFlag==1
        %plot(IonicStrength-IonicStrengthOld, 'r-');
        %legend(['ode counter =', num2str(ODEcounter)]);
        IonicStrengthOld=IonicStrength;
    end
    %fprintf('maximum difference in Ionic Strength = %f \n',...
    %    max(abs(IonicStrength-IonicStrengthOld)));

    %=======================================

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

        keyIn = get(gcf, 'CurrentCharacter');
        if (Res(end,end)/Res(1,end)<0.2 && PrepareGridStage) || strcmp(keyIn,'m')
            %tic
            disp('Migration begins');
            PrepareGridStage=0;
            IonicCalcFlag=1*IonicEffectFlag;
            set(gcf, 'CurrentCharacter','t');
            clear Res
        end

        set(handles.DataHolder,'UserData',[phiVec,cMat]);

        if IonicEffectFlag==1
           pH=-log10(cH.*10.^(-0.50850*(sqrt(IonicStrength)./(1+sqrt(IonicStrength))-0.3*IonicStrength)));
           pHVecAllTimes(:,PlotCounter)=pH;
        else
           pH=-log10(cH);
           pHVecAllTimes(:,PlotCounter)=pH;
        end


        PlotCurrentSelection(handles,phiVec,cMat,muMat,DMat,dPotendxMat,pH,t,h,err,SigVec,dphidzVec, phiVec*0 + AChannel);

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
    save ([get(handles.FilenameEdit,'String'), 'at'],'Nspecies','cMat', ...
        'phiVec','L1','L2','N','t','dtInit','cMatAllTimes','tVecOut','phiVecAllTimes',....
        'uAllTimes','cH','muMat','Res','dPotendxMat', 'SigVecAllTime', 'pHVecAllTimes');
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
function [cizCube,cH,cHCubePower,gizCube]=LzCalcEquilibrium(cH,LCube,cMat, ValCube, Kw_new)
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
    cH=cH-MyLzFuncNEW(cH,LCube,cMat, ValCube, Kw_new);
    if count>200
        disp('Too many iterations on cH. Returning');
        return;
    end
end
cMat=cMat*met2lit;
%fprintf('no. of iterations = %g \n', count);
% else
%     disp('initial condition is already the solution');
%     fprintf('Max absolute Residue = %g \n', max(abs(F)));
% end


cHPolMat=[ones(1,Ngrid);cumprod(ones(PolDeg-1,1)*cH,1)];
cHCubePower=FastRepmatRows(permute(cHPolMat,[3,2,1]),Nspecies);
Temp=sum((LCube.*cHCubePower),3);
M1Cube=FastRepmatPages(cMat'./Temp,PolDeg);
cizCube=LCube.*cHCubePower.*M1Cube; %equation 5 in J. Chroma A
gizCube=LCube.*cHCubePower.*FastRepmatPages(1./Temp,PolDeg);

%variable list
%cMat Ngrid x Nspecies
%cH   1 x Ngrid
%
%
%set(handles.Whattoplotpopup, 'UserData',1)
function MainRun(InputFilename,LoadDataFlag,handles)
global F Rmu Nspecies Temp Kw muH muOH epsilon visc DeltaP% uCST
global PrepareGridStage PrevTotalCost AChannel
global met2lit
global PolDeg Ngrid

disp('Using Spresso With Ionic Strength');

set(handles.WhatToPlotPopup, 'UserData',1) %sets WhatToPlotPopup->UserData=1
%means that now MainRun is working and in function
%SpressoGUI->WhatToPlotPopup_Callback now dont
%call Spresso('ChangePlotVariable',handles);

%profile on;
tic

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

atol=1E-20;
dtMin=0; dtMax=1E5;
dtInit=0.00001;                  % Initial Time step [sec]

RKOrder=23;

DisableChemFlag=0;

if( strcmp(SpatialDiscFlag, 'SLIP'))
    rtol=1E-3;
    AdaptGrid.Nconv=10;
    SmoothEps=1;
else
  rtol=1E-3;
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
Ngrid=size(cMat,1); %grid size, global variable
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

%profile on;
%--------------------------------------------------------------------------
%  MARCH IN TIME
%--------------------------------------------------------------------------



if( strcmp(SpatialDiscFlag, 'SLIP'))

    dz=zVec(2)-zVec(1);
    cSize=size(cMat);

    ARGS={cSize,DeltatNextPlot,DeltaCounterNextPlot,dtInit,dtMin,dtMax,rtol, ...
        atol,normcontrol,RK_A,RK_B,RK_E,RK_pow,uCST,A1,B1,A2,B2,M,Cur,INP, ...
        AdaptGrid,RKOrder,SteadyStateFlag,L1,L2,zetaPot,hChannel,bPressure,DeltaP,...
        betaDispersion,Res,DisableChemFlag,SpatialDiscFlag,IonicEffectFlag, ...
        PcentIonicChange, handles, dz};

    %hard code LocalODE45FV for finite volume scheme with conservation variable J*A*c
    disp('Running Finite Volume SLIP Scheme');
    AreaRatioCellVec = AreaRatio; %getAreaRatio(phiVec);
    AreaRatioCellMat = AreaRatioCellVec*ones(1,Nspecies);

    qMat=cMat.*AreaRatioCellMat; %here J=dx/dz=1
    cSize=size(cMat);

    IC=[reshape(qMat,numel(cMat),1);phiVec; AreaRatioCellVec];

    LocalOde45FV(IC,[t0 tEnd],ARGS);%pass initial condition in IC in y0
else

    cSize=size(cMat);
    IC=[reshape(cMat,numel(cMat),1);phiVec];
    %PcentIonicChange=1; %now input is from handles
    %IonicEffectFlag=1; %make it zero if ionic effects are off
    %later this should be loaded from input file
    ARGS={cSize,DeltatNextPlot,DeltaCounterNextPlot,dtInit,dtMin,dtMax,rtol, ...
        atol,normcontrol,RK_A,RK_B,RK_E,RK_pow,uCST,A1,B1,A2,B2,M,Cur,INP, ...
        AdaptGrid,RKOrder,SteadyStateFlag,L1,L2,zetaPot,hChannel,bPressure,DeltaP,...
        betaDispersion,Res,DisableChemFlag,SpatialDiscFlag,IonicEffectFlag, PcentIonicChange, handles};

    if (max(abs(AreaRatio-1))>0)
    fprintf('Using constant channel cross-section\n');
    fprintf('Variable cross-section is avialable only with SLIP scheme\n');
    end


    LocalOde45(IC,[t0 tEnd],ARGS);

end
set(handles.WhatToPlotPopup, 'UserData',0) %sets WhatToPlotPopup->UserData=0
toc
%profile viewer;

function c = multiprod(a, b, idA, idB)
%MULTIPROD  Multiplying 1-D or 2-D subarrays contained in two N-D arrays.
%   C = MULTIPROD(A,B) is equivalent  to C = MULTIPROD(A,B,[1 2],[1 2])
%   C = MULTIPROD(A,B,[D1 D2]) is eq. to C = MULTIPROD(A,B,[D1 D2],[D1 D2])
%   C = MULTIPROD(A,B,D1) is equival. to C = MULTIPROD(A,B,D1,D1)
%
%   MULTIPROD performs multiple matrix products, with array expansion (AX)
%   enabled. Its first two arguments A and B are "block arrays" of any
%   size, containing one or more 1-D or 2-D subarrays, called "blocks" (*).
%   For instance, a 5???6???3 array may be viewed as an array containing five
%   6???3 blocks. In this case, its size is denoted by 5???(6???3). The 1 or 2
%   adjacent dimensions along which the blocks are contained are called the
%   "internal dimensions" (IDs) of the array (???).
%
%   1) 2-D by 2-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, [DA1 DA2], [DB1 DB2]) contains the products
%         of the P???Q matrices in A by the R???S matrices in B. [DA1 DA2] are
%         the IDs of A; [DB1 DB2] are the IDs of B.
%
%   2) 2-D by 1-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, [DA1 DA2], DB1) contains the products of the
%         P???Q matrices in A by the R-element vectors in B. The latter are
%         considered to be R???1 matrices. [DA1 DA2] are the IDs of A; DB1 is
%         the ID of B.
%
%   3) 1-D by 2-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, DA1, [DB1 DB2]) contains the products of the
%         Q-element vectors in A by the R???S matrices in B. The vectors in A
%         are considered to be 1???Q matrices. DA1 is the ID of A; [DB1 DB2]
%         are the IDs of B.
%
%   4) 1-D BY 1-D BLOCK(S) (*)
%      (a) If either SIZE(A, DA1) == 1 or SIZE(B, DB1) == 1, or both,
%             C = MULTIPROD(A, B, DA1, DB1) returns products of scalars by
%             vectors, or vectors by scalars or scalars by scalars.
%      (b) If SIZE(A, DA1) == SIZE(B, DB1),
%             C = MULTIPROD(A, B, [0 DA1], [DB1 0]) or
%             C = MULTIPROD(A, B, DA1, DB1) virtually turns the vectors
%             contained in A and B into 1???P and P???1 matrices, respectively,
%             then returns their products, similar to scalar products.
%             Namely, C = DOT2(A, B, DA1, DB1) is equivalent to
%             C = MULTIPROD(CONJ(A), B, [0 DA1], [DB1 0]).
%      (c) Without limitations on the length of the vectors in A and B,
%             C = MULTIPROD(A, B, [DA1 0], [0 DB1]) turns the vectors
%             contained in A and B into P???1 and 1???Q matrices, respectively,
%             then returns their products, similar to outer products.
%             Namely, C = OUTER(A, B, DA1, DB1) is equivalent to
%             C = MULTIPROD(CONJ(A), B, [DA1 0], [0 DB1]).
%
%   Common constraints for all syntaxes:
%      The external dimensions of A and B must either be identical or
%      compatible with AX rules. The internal dimensions of each block
%      array must be adjacent (DA2 == DA1 + 1 and DB2 == DB1 + 1 are
%      required). DA1 and DB1 are allowed to be larger than NDIMS(A) and
%      NDIMS(B). In syntaxes 1, 2, and 3, Q == R is required, unless the
%      blocks in A or B are scalars.
%
%   Array expansion (AX):
%      AX is a powerful generalization to N-D of the concept of scalar
%      expansion. Indeed, A and B may be scalars, vectors, matrices or
%      multi-dimensional arrays. Scalar expansion is the virtual
%      replication or annihilation of a scalar which allows you to combine
%      it, element by element, with an array X of any size (e.g. X+10,
%      X*10, or []-10). Similarly, in MULTIPROD, the purpose of AX is to
%      automatically match the size of the external dimensions (EDs) of A
%      and B, so that block-by-block products can be performed. ED matching
%      is achieved by means of a dimension shift followed by a singleton
%      expansion:
%      1) Dimension shift (see SHIFTDIM).
%            Whenever DA1 ~= DB1, a shift is applied to impose DA1 == DB1.
%            If DA1 > DB1, B is shifted to the right by DA1 - DB1 steps.
%            If DB1 > DA1, A is shifted to the right by DB1 - DA1 steps.
%      2) Singleton expansion (SX).
%            Whenever an ED of either A or B is singleton and the
%            corresponding ED of the other array is not, the mismatch is
%            fixed by virtually replicating the array (or diminishing it to
%            length 0) along that dimension.
%
%   MULTIPROD is a generalization for N-D arrays of the matrix
%   multiplication function MTIMES, with AX enabled. Vector inner, outer,
%   and cross products generalized for N-D arrays and with AX enabled are
%   performed by DOT2, OUTER, and CROSS2 (MATLAB Central, file #8782).
%   Elementwise multiplications (see TIMES) and other elementwise binary
%   operations with AX enabled are performed by BAXFUN (MATLAB Central,
%   file #23084). Together, these functions make up the ???ARRAYLAB toolbox???.
%
%   Input and output format:
%      The size of the EDs of C is determined by AX. Block size is
%      determined as follows, for each of the above-listed syntaxes:
%      1) C contains P???S matrices along IDs MAX([DA1 DA2], [DB1 DB2]).
%      2) Array     Block size     ID(s)
%         ----------------------------------------------------
%         A         P???Q  (2-D)     [DA1 DA2]
%         B         R    (1-D)     DB1
%         C (a)     P    (1-D)     MAX(DA1, DB1)
%         C (b)     P???Q  (2-D)     MAX([DA1 DA2], [DB1 DB1+1])
%         ----------------------------------------------------
%         (a) The 1-D blocks in B are not scalars (R > 1).
%         (b) The 1-D blocks in B are scalars (R = 1).
%      3) Array     Block size     ID(s)
%         ----------------------------------------------------
%         A           Q  (1-D)     DA1
%         B         R???S  (2-D)     [DB1 DB2]
%         C (a)       S  (1-D)     MAX(DA1, DB1)
%         C (b)     R???S  (2-D)     MAX([DA1 DA1+1], [DB1 DB2])
%         ----------------------------------------------------
%         (a) The 1-D blocks in A are not scalars (Q > 1).
%         (b) The 1-D blocks in A are scalars (Q = 1).
%      4)     Array     Block size         ID(s)
%         --------------------------------------------------------------
%         (a) A         P        (1-D)     DA1
%             B         Q        (1-D)     DB1
%             C         MAX(P,Q) (1-D)     MAX(DA1, DB1)
%         --------------------------------------------------------------
%         (b) A         P        (1-D)     DA1
%             B         P        (1-D)     DB1
%             C         1        (1-D)     MAX(DA1, DB1)
%         --------------------------------------------------------------
%         (c) A         P        (1-D)     DA1
%             B         Q        (1-D)     DB1
%             C         P???Q      (2-D)     MAX([DA1 DA1+1], [DB1 DB1+1])
%         --------------------------------------------------------------
%
%   Terminological notes:
%   (*) 1-D and 2-D blocks are generically referred to as "vectors" and
%       "matrices", respectively. However, both may be also called
%       ???scalars??? if they have a single element. Moreover, matrices with a
%       single row or column (e.g. 1???3 or 3???1) may be also called ???row
%       vectors??? or ???column vectors???.
%   (???) Not to be confused with the "inner dimensions" of the two matrices
%       involved in a product X * Y, defined as the 2nd dimension of X and
%       the 1st of Y (DA2 and DB1 in syntaxes 1, 2, 3).
%
%   Examples:
%    1) If  A is .................... a 5???(6???3)???2 array,
%       and B is .................... a 5???(3???4)???2 array,
%       C = MULTIPROD(A, B, [2 3]) is a 5???(6???4)???2 array.
%
%       A single matrix A pre-multiplies each matrix in B
%       If  A is ........................... a (1???3)    single matrix,
%       and B is ........................... a 10???(3???4) 3-D array,
%       C = MULTIPROD(A, B, [1 2], [3 4]) is a 10???(1???4) 3-D array.
%
%       Each matrix in A pre-multiplies each matrix in B (all possible
%       combinations)
%       If  A is .................... a (6???3)???5   array,
%       and B is .................... a (3???4)???1???2 array,
%       C = MULTIPROD(A, B, [1 2]) is a (6???4)???5???2 array.
%
%   2a) If  A is ........................... a 5???(6???3)???2 4-D array,
%       and B is ........................... a 5???(3)???2   3-D array,
%       C = MULTIPROD(A, B, [2 3], [2]) is   a 5???(6)???2   3-D array.
%
%   2b) If  A is ........................... a 5???(6???3)???2 4-D array,
%       and B is ........................... a 5???(1)???2   3-D array,
%       C = MULTIPROD(A, B, [2 3], [2]) is   a 5???(6???3)???2 4-D array.
%
%   4a) If both A and B are .................. 5???(6)???2   3-D arrays,
%       C = MULTIPROD(A, B, 2) is .......... a 5???(1)???2   3-D array, while
%   4b) C = MULTIPROD(A, B, [2 0], [0 2]) is a 5???(6???6)???2 4-D array
%
%   See also DOT2, OUTER, CROSS2, BAXFUN, MULTITRANSP.

% $ Version: 2.1 $
% CODE      by:            Paolo de Leva
%                          (Univ. of Rome, Foro Italico, IT)    2009 Jan 24
%           optimized by:  Paolo de Leva
%                          Jinhui Bai (Georgetown Univ., D.C.)  2009 Jan 24
% COMMENTS  by:            Paolo de Leva                        2009 Feb 24
% OUTPUT    tested by:     Paolo de Leva                        2009 Feb 24
% -------------------------------------------------------------------------

error( nargchk(2, 4, nargin) ); % Allow 2 to 4 input arguments
switch nargin % Setting IDA and/or IDB
    case 2, idA = [1 2]; idB = [1 2];
    case 3, idB = idA;
end

% ESC 1 - Special simple case (both A and B are 2D), solved using C = A * B

if ndims(a)==2 && ndims(b)==2 && ...
        isequal(idA,[1 2]) && isequal(idB,[1 2])
    c = a * b; return
end

% MAIN 0 - Checking and evaluating array size, block size, and IDs

sizeA0 = size(a);
sizeB0 = size(b);
[sizeA, sizeB, shiftC, delC, sizeisnew, idA, idB, ...
    squashOK, sxtimesOK, timesOK, mtimesOK, sumOK] = ...
    sizeval(idA,idB, sizeA0,sizeB0);

% MAIN 1 - Applying dimension shift (first step of AX) and
%          turning both A and B into arrays of either 1-D or 2-D blocks

if sizeisnew(1), a = reshape(a, sizeA); end
if sizeisnew(2), b = reshape(b, sizeB); end

% MAIN 2 - Performing products with or without SX (second step of AX)

if squashOK % SQUASH + MTIMES (fastest engine)
    c = squash2D_mtimes(a,b, idA,idB, sizeA,sizeB, squashOK);
elseif timesOK % TIMES (preferred w.r. to SX + TIMES)
    if sumOK, c = sum(a .* b, sumOK);
    else      c =     a .* b; end
elseif sxtimesOK % SX + TIMES
    if sumOK, c = sum(bsxfun(@times, a, b), sumOK);
    else      c =     bsxfun(@times, a, b); end
elseif mtimesOK % MTIMES (rarely used)
    c = a * b;
end

% MAIN 3 - Reshaping C (by inserting or removing singleton dimensions)

[sizeC sizeCisnew] = adjustsize(size(c), shiftC, false, delC, false);
if sizeCisnew, c = reshape(c, sizeC); end


function c = squash2D_mtimes(a, b, idA, idB, sizeA, sizeB, squashOK)
% SQUASH2D_MTIMES  Multiproduct with single-block expansion (SBX).
%    Actually, no expansion is performed. The multi-block array is
%    rearranged from N-D to 2-D, then MTIMES is applied, and eventually the
%    result is rearranged back to N-D. No additional memory is required.
%    One and only one of the two arrays must be single-block, and its IDs
%    must be [1 2] (MAIN 1 removes leading singletons). Both arrays
%    must contain 2-D blocks (MAIN 1 expands 1-D blocks to 2-D).

if squashOK == 1 % A is multi-block, B is single-block (squashing A)

    % STEP 1 - Moving IDA(2) to last dimension
    nd = length(sizeA);
    d2 = idA(2);
    order = [1:(d2-1) (d2+1):nd d2]; % Partial shifting
    a = permute(a, order); % ...???Q

    % STEP 2 - Squashing A from N-D to 2-D
    q = sizeB(1);
    s = sizeB(2);
    lengthorder = length(order);
    collapsedsize = sizeA(order(1:lengthorder-1));
    n = prod(collapsedsize);
    a = reshape(a, [n, q]); % N???Q
    fullsize = [collapsedsize s]; % Size to reshape C back to N-D

else % B is multi-block, A is single-block (squashing B)

    % STEP 1 - Moving IDB(1) to first dimension
    nd = length(sizeB);
    d1 = idB(1);
    order = [d1 1:(d1-1) (d1+1):nd]; % Partial shifting
    b = permute(b, order); % Q???...

    % STEP 2 - Squashing B from N-D to 2-D
    p = sizeA(1);
    q = sizeA(2);
    lengthorder = length(order);
    collapsedsize = sizeB(order(2:lengthorder));
    n = prod(collapsedsize);
    b = reshape(b, [q, n]); % Q???N
    fullsize = [p collapsedsize]; % Size to reshape C back to N-D

end

% FINAL STEPS - Multiplication, reshape to N-D, inverse permutation
invorder(order) = 1 : lengthorder;
c = permute (reshape(a*b, fullsize), invorder);


function [sizeA, sizeB, shiftC, delC, sizeisnew, idA, idB, ...
    squashOK, sxtimesOK, timesOK, mtimesOK, sumOK] = ...
    sizeval(idA0,idB0, sizeA0,sizeB0)
%SIZEVAL   Evaluation of array size, block size, and IDs
%    Possible values for IDA and IDB:
%        [DA1 DA2], [DB1 DB2]
%        [DA1 DA2], [DB1]
%        [DA1],     [DB1 DB2]
%        [DA1],     [DB1]
%        [DA1 0],   [0 DB1]
%        [0 DA1],   [DB1 0]
%
%    sizeA/B     Equal to sizeA0/B0 if RESHAPE is not needed in MAIN 1
%    shiftC, delC    Variables controlling MAIN 3.
%    sizeisnew   1x2 logical array; activates reshaping of A and B.
%    idA/B       May change only if squashOK ~= 0
%    squashOK    If only A or B is a multi-block array (M-B) and the other
%                is single-block (1-B), it will be rearranged from N-D to
%                2-D. If both A and B are 1-B or M-B arrays, squashOK = 0.
%                If only A (or B) is a M-B array, squashOK = 1 (or 2).
%    sxtimesOK, timesOK, mtimesOK    Flags controlling MAIN 2 (TRUE/FALSE).
%    sumOK       Dimension along which SUM is performed. If SUM is not
%                needed, sumOK = 0.

% Initializing output arguments

idA = idA0;
idB = idB0;
squashOK = 0;
sxtimesOK = false;
timesOK = false;
mtimesOK = false;
sumOK = 0;
shiftC = 0;
delC = 0;

% Checking for gross input errors

NidA = numel(idA);
NidB = numel(idB);
idA1 = idA(1);
idB1 = idB(1);
if  NidA>2 || NidB>2 || NidA==0 || NidB==0 || ...
        ~isreal(idA1) ||    ~isreal(idB1)   || ...
        ~isnumeric(idA1) || ~isnumeric(idB1)   || ...
        0>idA1  ||          0>idB1    || ... % negative
        idA1~=fix(idA1) ||  idB1~=fix(idB1)   || ... % non-integer
        ~isfinite(idA1) ||  ~isfinite(idB1) % Inf or NaN
    error('MULTIPROD:InvalidDimensionArgument', ...
        ['Internal-dimension arguments (e.g., [IDA1 IDA2]) must\n', ...
        'contain only one or two non-negative finite integers']);
end

% Checking Syntaxes containing zeros (4b/c)

declared_outer = false;
idA2 = idA(NidA); % It may be IDA1 = IDA2 (1-D block)
idB2 = idB(NidB);

if any(idA==0) || any(idB==0)

    % "Inner products": C = MULTIPROD(A, B, [0 DA1], [DB1 0])
    if idA1==0 && idA2>0 && idB1>0 && idB2==0
        idA1 = idA2;
        idB2 = idB1;
        % "Outer products": C = MULTIPROD(A, B, [DA1 0], [0 DB1])
    elseif idA1>0 && idA2==0 && idB1==0 && idB2>0
        declared_outer = true;
        idA2 = idA1;
        idB1 = idB2;
    else
        error('MULTIPROD:InvalidDimensionArgument', ...
            ['Misused zeros in the internal-dimension arguments\n', ...
            '(see help heads 4b and 4c)']);
    end
    NidA = 1;
    NidB = 1;
    idA = idA1;
    idB = idB1;

elseif (NidA==2 && idA2~=idA1+1) || ...  % Non-adjacent IDs
        (NidB==2 && idB2~=idB1+1)
    error('MULTIPROD:InvalidDimensionArgument', ...
        ['If an array contains 2-D blocks, its two internal dimensions', ...
        'must be adjacent (e.g. IDA2 == IDA1+1)']);
end

% ESC - Case for which no reshaping is needed (both A and B are scalars)

scalarA = isequal(sizeA0, [1 1]);
scalarB = isequal(sizeB0, [1 1]);
if scalarA && scalarB
    sizeA = sizeA0;
    sizeB = sizeB0;
    sizeisnew = [false false];
    timesOK = true; return
end

% Computing and checking adjusted sizes
% The lengths of ADJSIZEA and ADJSIZEB must be >= IDA(END) and IDB(END)

NsA = idA2 - length(sizeA0); % Number of added trailing singletons
NsB = idB2 - length(sizeB0);
adjsizeA = [sizeA0 ones(1,NsA)];
adjsizeB = [sizeB0 ones(1,NsB)];
extsizeA = adjsizeA([1:idA1-1, idA2+1:end]); % Size of EDs
extsizeB = adjsizeB([1:idB1-1, idB2+1:end]);
p = adjsizeA(idA1);
q = adjsizeA(idA2);
r = adjsizeB(idB1);
s = adjsizeB(idB2);
scalarsinA = (p==1 && q==1);
scalarsinB = (r==1 && s==1);
singleA = all(extsizeA==1);
singleB = all(extsizeB==1);
if q~=r && ~scalarsinA && ~scalarsinB && ~declared_outer
    error('MULTIPROD:InnerDimensionsMismatch', ...
        'Inner matrix dimensions must agree.');
end

% STEP 1/3 - DIMENSION SHIFTING (FIRST STEP OF AX)
%   Pipeline 1 (using TIMES) never needs left, and may need right shifting.
%   Pipeline 2 (using MTIMES) may need left shifting of A and right of B.

shiftA = 0;
shiftB = 0;
diffBA = idB1 - idA1;
if scalarA % Do nothing
elseif singleA && ~scalarsinB, shiftA = -idA1 + 1; %  Left shifting A
elseif idB1 > idA1,            shiftA = diffBA;    % Right shifting A
end
if scalarB % Do nothing
elseif singleB && ~scalarsinA, shiftB = -idB1 + 1; %  Left shifting B
elseif idA1 > idB1,            shiftB = -diffBA;   % Right shifting B
end

% STEP 2/3 - SELECTION OF PROPER ENGINE AND BLOCK SIZE ADJUSTMENTS

addA  = 0; addB  = 0;
delA  = 0; delB  = 0;
swapA = 0; swapB = 0;
idC1 = max(idA1, idB1);
idC2 = idC1 + 1;
checktimes = false;

if (singleA||singleB) &&~scalarsinA &&~scalarsinB % Engine using MTIMES

    if singleA && singleB
        mtimesOK = true;
        shiftC=idC1-1; % Right shifting C
        idC1=1; idC2=2;
    elseif singleA
        squashOK = 2;
        idB = [idB1, idB1+1] + shiftB;
    else % singleB
        squashOK = 1;
        idA = [idA1, idA1+1] + shiftA;
    end

    if NidA==2 && NidB==2 % 1) 2-D BLOCKS BY 2-D BLOCKS
        % OK
    elseif NidA==2        % 2) 2-D BLOCKS BY 1-D BLOCKS
        addB=idB1+1; delC=idC2;
    elseif NidB==2        % 3) 1-D BLOCKS BY 2-D BLOCKS
        addA=idA1; delC=idC1;
    else                  % 4) 1-D BLOCKS BY 1-D BLOCKS
        if declared_outer
            addA=idA1+1; addB=idB1;
        else
            addA=idA1; addB=idB1+1; delC=idC2;
        end
    end

else % Engine using TIMES (also used if SCALARA || SCALARB)

    sxtimesOK = true;

    if NidA==2 && NidB==2 % 1) 2-D BLOCKS BY 2-D BLOCKS

        if scalarA || scalarB
            timesOK=true;
        elseif scalarsinA && scalarsinB % scal-by-scal
            checktimes=true;
        elseif scalarsinA || scalarsinB || ... % scal-by-mat
                (q==1 && r==1)  % vec-by-vec ("outer")
        elseif p==1 && s==1 % vec-by-vec ("inner")
            swapA=idA1; sumOK=idC1; checktimes=true;
        elseif s==1 % mat-by-vec
            swapB=idB1; sumOK=idC2;
        elseif p==1 % vec-by-mat
            swapA=idA1; sumOK=idC1;
        else % mat-by-mat
            addA=idA2+1; addB=idB1; sumOK=idC2; delC=idC2;
        end

    elseif NidA==2 % 2) 2-D BLOCKS BY 1-D BLOCKS

        if scalarA || scalarB
            timesOK=true;
        elseif scalarsinA && scalarsinB % scal-by-scal
            addB=idB1; checktimes=true;
        elseif scalarsinA % scal-by-vec
            delA=idA1;
        elseif scalarsinB % mat-by-scal
            addB=idB1;
        elseif p==1 % vec-by-vec ("inner")
            delA=idA1; sumOK=idC1; checktimes=true;
        else % mat-by-vec
            addB=idB1; sumOK=idC2; delC=idC2;
        end

    elseif NidB==2 % 3) 1-D BLOCKS BY 2-D BLOCKS

        if scalarA || scalarB
            timesOK=true;
        elseif scalarsinA && scalarsinB % scal-by-scal
            addA=idA1+1; checktimes=true;
        elseif scalarsinB % vec-by-scal
            delB=idB2;
        elseif scalarsinA % scal-by-mat
            addA=idA1+1;
        elseif s==1 % vec-by-vec ("inner")
            delB=idB2; sumOK=idC1; checktimes=true;
        else % vec-by-mat
            addA=idA1+1; sumOK=idC1; delC=idC1;
        end

    else % 4) 1-D BLOCKS BY 1-D BLOCKS

        if scalarA || scalarB
            timesOK=true;
        elseif declared_outer % vec-by-vec ("outer")
            addA=idA1+1; addB=idB1;
        elseif scalarsinA && scalarsinB % scal-by-scal
            checktimes=true;
        elseif scalarsinA || scalarsinB % vec-by-scal
        else % vec-by-vec
            sumOK=idC1; checktimes=true;
        end
    end
end

% STEP 3/3 - Adjusting the size of A and B. The size of C is adjusted
%            later, because it is not known yet.

[sizeA, sizeisnew(1)] = adjustsize(sizeA0, shiftA, addA, delA, swapA);
[sizeB, sizeisnew(2)] = adjustsize(sizeB0, shiftB, addB, delB, swapB);

if checktimes % Faster than calling BBXFUN
    diff = length(sizeB) - length(sizeA);
    if isequal([sizeA ones(1,diff)], [sizeB ones(1,-diff)])
        timesOK = true;
    end
end


function [sizeA, sizeisnew] = adjustsize(sizeA0, shiftA, addA, delA, swapA)
% ADJUSTSIZE  Adjusting size of a block array.

% Dimension shifting (by adding or deleting trailing singleton dim.)
if     shiftA>0, [sizeA,newA1] = addsing(sizeA0, 1, shiftA);
elseif shiftA<0, [sizeA,newA1] = delsing(sizeA0, 1,-shiftA);
else   sizeA = sizeA0;  newA1  = false;
end
% Modifying block size (by adding, deleting, or moving singleton dim.)
if      addA, [sizeA,newA2] = addsing(sizeA, addA+shiftA, 1); % 1D-->2D
elseif  delA, [sizeA,newA2] = delsing(sizeA, delA+shiftA, 1); % 2D-->1D
elseif swapA, [sizeA,newA2] = swapdim(sizeA,swapA+shiftA); % ID Swapping
else                 newA2  = false;
end
sizeisnew = newA1 || newA2;


function [newsize, flag] = addsing(size0, dim, ns)
%ADDSING   Adding NS singleton dimensions to the size of an array.
%   Warning: NS is assumed to be a positive integer.
%   Example: If the size of A is ..... SIZE0 = [5 9 3]
%            NEWSIZE = ADDSING(SIZE0, 3, 2) is [5 9 1 1 3]

if dim > length(size0)
    newsize = size0;
    flag = false;
else
    newsize = [size0(1:dim-1), ones(1,ns), size0(dim:end)];
    flag = true;
end


function [newsize, flag] = delsing(size0, dim, ns)
%DELSING   Removing NS singleton dimensions from the size of an array.
%   Warning: Trailing singletons are not removed
%   Example: If the size of A is SIZE0 = [1 1 1 5 9 3]
%            NEWSIZE = DELSING(SIZE, 1, 3) is  [5 9 3]

if dim > length(size0)-ns % Trailing singletons are not removed
    newsize = size0;
    flag = false;
else % Trailing singl. added, so NEWSIZE is guaranteed to be 2D or more
    newsize = size0([1:dim-1, dim+ns:end, dim]);
    flag = true;
end


function [newsize, flag] = swapdim(size0, dim)
%SWAPDIM   Swapping two adjacent dimensions of an array (DIM and DIM+1).
%   Used only when both A and B are multi-block arrays with 2-D blocks.
%   Example: If the size of A is .......... 5???(6???3)
%            NEWSIZE = SWAPIDS(SIZE0, 2) is 5???(3???6)

newsize = [size0 1]; % Guarantees that dimension DIM+1 exists.
newsize = newsize([1:dim-1, dim+1, dim, dim+2:end]);
flag = true;

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
function F=MyLzFuncNEW(cH,LCube,cMat, ValCube, Kw_new)
global Ngrid Nspecies PolDeg Kw met2lit


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
M1Cube=FastRepmatPages(cMat'./Temp,PolDeg);
TempCube=FastRepmatPages(Temp,PolDeg);
cizCube=LCube.*cHCubePower.*M1Cube; %equation 5 in J. Chroma A

%M1Cube is denominator of eqn 5


Temp_z = sum((ValCube.*LCube.*cHCubePower),3);
Temp_z_Cube=FastRepmatPages(Temp_z,PolDeg);
RHS_den=sum(sum((cizCube.*ValCube.^2-cizCube.*ValCube.*Temp_z_Cube./TempCube),3),1);

%M1Cube is denominator of eqn 5
%RHS_den=sum(sum((cizCube.*ValCube.^2-cizCube.*ValCube.^2.*LCube.*cHCubePower./TempCube),3),1);
RHS_num=sum(sum((cizCube.*ValCube),3),1);

F_num=RHS_num+cH-Kw_new./cH;
F_den=RHS_den./cH+1+Kw_new./(cH.^2);
F=F_num./F_den; %for Newton raphson

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
function F=MyLzFuncNEWTest(cH,LCube,cMat, ValCube, Kw_new)
global Ngrid Nspecies PolDeg Kw


cMat=cMat/1000;

cHPolMat=[ones(1,Ngrid);cumprod(ones(PolDeg-1,1)*cH,1)];
cHCubePower=FastRepmatRows(permute(cHPolMat,[3,2,1]),Nspecies);
%Temp=multiprod(TempL, cHPolMat, [1 2], [1]); %faster multiplication
Temp=sum((LCube.*cHCubePower),3);
M1Cube=FastRepmatPages(cMat'./Temp,PolDeg);
TempCube=FastRepmatPages(Temp,PolDeg);
cizCube=LCube.*cHCubePower.*M1Cube; %equation 5 in J. Chroma A

%M1Cube is denominator of eqn 5

RHS_den=sum(sum((cizCube.*ValCube.^2-cizCube.*ValCube.^2.*LCube.*cHCubePower./TempCube),3),1);
RHS_num=sum(sum((cizCube.*ValCube),3),1);
F_num=RHS_num+cH-Kw_new./cH;

% F_den=RHS_den./cH+1+Kw_new./(cH.^2);
%F_den=RHS_den./cH+1+Kw_new./(cH.^2);
%F=RHS_num;
F=RHS_den./cH;
%F=F_num./F_den; %for Newton raphson


%this function computes the mobility based on Onsager-Fuoss theory.
%the algo is:
%a) extract the non zero entries of cizCube, ValCube, muCube using find
%function.. (later give this as input)
%construct: omega,
%Follow the following for each grid point....
%b) once you have row matrices of mobility, valence and concentrations.
%construct: mu i.e. ratio c*z^2/IonicStrength
%c) now construct matrix H
%d) using H construct vectors r0, r1, ... , r5
%e) compute required factors for applying Onsager-Fuoss correction
%Repeat this for all grid points

%f) assemble back the matrices using the result of find
function muIonicCube = OnsagerFuoss(IonicEffectFlag, ValCube, cizCube, muCube)
global Nspecies PolDeg Ngrid F

muIonicCube=F*muCube;
if (IonicEffectFlag==0)
    return;

elseif (IonicEffectFlag==1)

    %extract the non zero elements
    [X, Y]=find(ValCube(:,1,:)~=0); %need to do this for one matrix only
    %same for all grid points

    %muCube is actually Omega in Onsage-Fuoss paper and is same for
    %all grid points
    c=[0.2929 -0.3536 0.0884 -0.0442 0.0276 -0.0193];
    %coefficients in onsager-fuoss paper
    for i=1:1:length(X)
        omega(1,i)=muCube(X(i),1,Y(i));
        z(1,i)=ValCube(X(i),1,Y(i));
    end
    mob=omega.*abs(z);
    XLen=length(X);
    mob_new=zeros(1,XLen);
    conc=zeros(1,length(X)); %concentrations

    II=eye(XLen,XLen);

    for k=1:1:Ngrid
        for i=1:1:XLen
            conc(1,i)=cizCube(X(i),k,Y(i));
        end
        %extract concentrations

        %now that we have concentrations at k_th grid point
        %make matrices mu
        IonicStr=sum(conc.*z.^2);
        mu=conc.*z.^2./IonicStr; %total potential

        %now make matrix H
        for j=1:XLen
            for i=1:XLen
                h(j,i)=mu(i)*omega(i)/(omega(i)+omega(j));
            end
        end

        d=sum(h,2);
        d1=diag(d);
        h=h+d1; %makes matrix H.. adds the delta_(i,j) terms

        B=2*h-II;

        r(:,1)=(z-sum(z.*mu)/sum(mu.*abs(z)./mob)*(abs(z)./mob))'; %check for absolute signs
        for i=2:6
            r(:,i)=B*r(:,i-1);
        end
        factor=c*r'; %factors computed in table 3, page 2755 of Onsager-Fuoss

        %Now compute the ionic strength dependence
        mob_new=F*omega-(F*0.78420*z.*factor.*omega+31.410e-9).*sqrt(IonicStr/2000)./(1+1.5*sqrt(IonicStr/2000));
        %mob_new=F*omega-(F*0.78420*z.*factor.*omega+31.410e-9).*sqrt(IonicStr/2000);
        %Robinson-Stokes, just checking dont use this.
        %mob_new=F*omega-(F*0.2297*z.*omega+31.410e-9).*sqrt(IonicStr/2000)./(1+1.5*sqrt(IonicStr/2000));
        %Assemble matrix back

        for i=1:1:XLen
            muIonicCube(X(i),k,Y(i))=mob_new(1,i);
        end

    end
end




function PlotCurrentSelection(handles,xVec,cMat,muMat,DMat,dPotendx,pH,t,h,err,SigVec,dphidzVec, CrossArea)
global met2lit DeltaP
global ColorList

WhatToPlot=get(handles.WhatToPlotPopup,'value');
axes(handles.MainAxes); set(handles.MainAxes,'nextplot','replacechildren');

switch WhatToPlot
    case (1)
        plot(handles.MainAxes,xVec*1E3,cMat/met2lit,'.-');              % concentration is in mol/m^3. Convert to mol/lit
        xlabel(handles.MainAxes,'x [mm]'); ylabel(handles.MainAxes,'c [mol/lit]');
    case(2)
        plot(handles.MainAxes,xVec*1E3,pH,'.-');
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
    case (7)
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
%IonicEffectFlag=1; %later get this from handles
IonicEffectFlag = 2-get(handles.IonicEffectPopup,'value');
[pH,SigVec,muMat,DMat] = CalculateSpatialProperties(INP,N,FirstTimeFlag,cMat, IonicEffectFlag);
t=NaN; h=NaN; err=NaN; dphidzVec=xVec+nan; dPotendxMat=xVec+nan; CrossArea=xVec+nan;

%%% #plotarea
AreaFunctionText =  get(handles.AreaVariation_Edit, 'String');
if isempty(AreaFunctionText)  AreaFunctionText='@(x) 1 + 0*x'; end
A = feval( eval(AreaFunctionText), xVec);
AreaRatio=A./max(A);

cMat=cMat*met2lit;                % Convert from mol/lit to mol/m^3
PlotCurrentSelection(handles,xVec,cMat,muMat,DMat,dPotendxMat,pH,t,h,err,SigVec,dphidzVec, AreaRatio);
set(handles.DataHolder,'UserData',[xVec,cMat]);
set(handles.StartButton,'enable','off');
%-------------------------------------------------------

%--------------------------------------------------------------------------
%  Recompute Equilibrium polonomials
%  Will be called in a new function called IonicStrengthEffects()
%  Instead of EquilibriumPolynomials() this will use inputs as KaMat etc,
%  As they are already in ordered form, instead of INP
%--------------------------------------------------------------------------

%updated Feb 27, 2009: Included Kw_new to include ionic strengths

%function [PMat,PPrimeMat,Q,QPrime,LMat]=RecomputeEquilibriumPolynomials(zListArranged, KaListCube)

function LCube=RecomputeEquilibriumPolynomials(zListArranged, KaListCube)
global Nspecies Ngrid PolDeg

LCube=zeros(Nspecies, Ngrid, PolDeg);

for k=1:1:Ngrid
    for j=1:1:Nspecies
        zList=zListArranged{j};
        KaList=KaListCube{k}{j};
        nj=min(zList);    %pj=max(zList);

        for z=zList
            if z<0
                LCube(j,k,z-nj+1)=prod(KaList(z-nj+1:-nj));
            elseif z>0
                LCube(j,k,z-nj+1)=1/prod(KaList(-nj+2:z-nj+1));
            elseif z==0
                LCube(j,k,z-nj+1)=1;
            end %if
        end % for z

    end %for j
end %for loop k


function Y = Myheaviside(X)
%HEAVISIDE    Step function.
%    HEAVISIDE(X) is 0 for X < 0, 1 for X > 0, and NaN for X == 0.

Y = zeros(size(X));
Y(X > 0) = 1;
Y(X == 0) = NaN;
