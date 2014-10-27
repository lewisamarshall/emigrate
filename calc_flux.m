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
