def LocalOde45(y0,tspan,ARGS):

N=cSize(1);

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
