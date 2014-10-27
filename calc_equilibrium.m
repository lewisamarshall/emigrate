function [PMat,PPrimeMat,Q,QPrime,M,zMat,muMat,KaMat,DMat]=EquilibriumPolynomials(INP)
global F Rmu Temp Kw

% PREPARE MATRICES
%------------------
MaxCol = -Inf;
for j = 1:size(INP,1)
    MaxCol=max([MaxCol,max(INP{j}(1:3:end))-min(INP{j}(1:3:end))+1]);
end

M = zeros(size(INP,1),MaxCol);

zMat = zeros(size(INP, 1),MaxCol);
muMat = zMat; KaMat=zMat; DMat=zMat;

for j = 1:size(INP,1)
    zList = INP{j}(1:3:end)
    muList = INP{j}(2:3:end)./(F*abs(zList));  % Defines Santiago mobility
    pKaList = INP{j}(3:3:end)
    KaList = 10.^(-pKaList)
    DList = Rmu*Temp*muList #diffusivity


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
