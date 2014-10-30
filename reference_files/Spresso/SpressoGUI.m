function varargout = SpressoGUI(varargin)
% SpressoGUI M-file for SpressoGUI.fig
%      SpressoGUI, by itself, creates a new SpressoGUI or raises the existing
%      singleton*.
%
%      H = SpressoGUI returns the handle to a new SpressoGUI or the handle to
%      the existing singleton*.
%
%      SpressoGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SpressoGUI.M with the given input
%      arguments.
%
%      SpressoGUI('Property','Value',...) creates a new SpressoGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpressoGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpressoGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpressoGUI

% Last Modified by GUIDE v2.5 09-Aug-2012 00:27:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpressoGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SpressoGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before SpressoGUI is made visible.
function SpressoGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpressoGUI (see VARARGIN)

% Choose default command line output for SpressoGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Create the database table
mtable=DatabaseTable_CreateFcn([], eventdata, handles);
handles.DatabaseTable=mtable;
guidata(handles.figure1,handles);

% Create the input table
mtable=InputTable_CreateFcn(hObject, eventdata, handles);
handles.InputTable=mtable;
guidata(handles.figure1,handles);
setappdata(handles.InputTablePanel,'mtable',mtable);

position=get(handles.Panel1,'position');
set(handles.Panel2,'position',position,'visible','off');
SpressoNoIonic('ChangePlotVariable',handles);

% UIWAIT makes SpressoGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = SpressoGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function mtable=DatabaseTable_CreateFcn(hObject, eventdata, handles)
warning('off','MATLAB:xlsread:Mode')
[num, txt]= xlsread('Database/SpressoDatabase.xls',1,'','basic');
%[num, txt]= xlsread(fullfile(ctfroot, 'Database','SpressoDatabase.xls'),1,'','basic');

index=1;
for ij=2:size(num,1)
    cell_data{index,1}=txt{ij,2};
    a=([num(index,5:end)]);  a=a(~isnan(a));

    str=[];
    for kl=1:length(a)/3
        str=[str,mat2str(a(3*(kl-1)+1:3*kl)),'   '];
    end

    cell_data{index,2}=str;
    index=index+1;
end
Column.Titles={'Name','Valence, Mobility, pKa'};

mtable = createTable(handles.DatabasePanel,Column.Titles,cell_data,0);
mtable.setEditable(1,false);
mtable.setEditable(2,false);
mtable.getTable.setAutoResizeMode(false);
mtable.getTable.getColumnModel.getColumn(0).setPreferredWidth(150);
mtable.getTable.getColumnModel.getColumn(1).setPreferredWidth(400);

function mtable=InputTable_CreateFcn(hObject, eventdata, handles)
Column.Titles={' ','Name','Valence, Mobility, pKa','c [M]','Location'};
Column.Editable=[false,true, true, true, true];
Column.Format={'text','text', 'text', 'numeric', 'choice'};%{'LE' 'TE' 'Analyte'}};
Column.Width={20,100,160,56,50};
cell_data={};


mtable = createTable(handles.InputTablePanel,Column.Titles,[],0);
mtable.setEditable(1,false);  % Disable editing first column (see note
jtable = mtable.getTable;
cb = javax.swing.JComboBox({'LE','TE','Analyte','Background'}); cb.setEditable(true);  % prepare an editable drop-down CellEditor
editor = javax.swing.DefaultCellEditor(cb);
jtable.getColumnModel.getColumn(4).setCellEditor(editor);  % assign this editor to second column (see Note 2)
setappdata(mtable,'handles',handles);
mtable.DataChangedCallback = @InputTableDataChangeCallback;
mtable.getTable.setAutoResizeMode(false);
mtable.getTable.getColumnModel.getColumn(0).setPreferredWidth(15);
mtable.getTable.getColumnModel.getColumn(1).setPreferredWidth(100);
mtable.getTable.getColumnModel.getColumn(2).setPreferredWidth(200);
mtable.getTable.getColumnModel.getColumn(3).setPreferredWidth(50);
mtable.setEditable(1,false);

function InputTableDataChangeCallback(varargin)
PartialHandles=getappdata(varargin{1},'handles');
mtable=getappdata(PartialHandles.InputTablePanel,'mtable');
handles=guihandles(PartialHandles.figure1);
handles.InputTable=mtable;

% Make sure concentration values are numeric
InputData=cell(handles.InputTable.getData);
for ij=1:size(InputData,1)
    Conc=InputData{ij,4};
    if ~isnumeric(Conc)
        Conc=str2num(Conc);
        handles.InputTable.getTable.setValueAt(Conc,ij-1,4-1);
    end
end % for ij

if isempty(handles)
    return;
end

if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('PrepareInput',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('PrepareInput',handles);
end

set(handles.StartButton,'enable','off');

function AddButton_Callback(hObject, eventdata, handles)
global ColorList

PreviousInput = cell(handles.InputTable.getData);

RowNumber = handles.DatabaseTable.getTable.getSelectedRows + 1;
Database  = cell(handles.DatabaseTable.getData);
if isempty(RowNumber)
    return;
end

CurrentRow=size(PreviousInput,1)+1;
CurrentColorString=num2str(round(255*ColorList(CurrentRow,:)),'%d,');

NewLine{1}=['<HTML><FONT color=rgb(',CurrentColorString(1:end-1),')>@@@</FONT></HTML>'];
NewLine{2}=Database{RowNumber,1};
NewLine{3}=Database{RowNumber,2};
NewLine{4}='0';
NewLine{5}='LE';
handles.InputTable.getTable.getModel.addRow(NewLine);  % appends a row to the bottom of the table

% Rearrange colors
for ij=1:CurrentRow
    CurrentColorString=num2str(round(255*ColorList(ij,:)),'%d,');
    handles.InputTable.getTable.setValueAt(['<HTML><FONT color=rgb(',CurrentColorString(1:end-1),')>@@@</FONT></HTML>'],ij-1,0);
end

if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('PrepareInput',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('PrepareInput',handles);
end

set(handles.StartButton,'enable','off');

function RemoveButton_Callback(hObject, eventdata, handles)

global ColorList
RowNumber = handles.InputTable.getTable.getSelectedRows;
if ~isempty(RowNumber)
    handles.InputTable.getTable.getModel.removeRow(RowNumber);
end

Input = cell(handles.InputTable.getData);
for ij=1:size(Input,1)
    CurrentColorString=num2str(round(255*ColorList(ij,:)),'%d,');
    handles.InputTable.getTable.setValueAt(['<HTML><FONT color=rgb(',CurrentColorString(1:end-1),')>@@@</FONT></HTML>'],ij-1,0);
end
if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('PrepareInput',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('PrepareInput',handles);
end

set(handles.StartButton,'enable','off');


% % --- Executes when selected cell(s) is changed in InputTable.
% function InputTable_CellSelectionCallback(hObject, eventdata, handles)
% if ~isempty(eventdata.Indices)
% RowNumber=eventdata.Indices(1);
% else
%     RowNumber=1;
% end
% set(hObject,'UserData',RowNumber);


function SteadyStatePopup_Callback(hObject, eventdata, handles)
msg={['Setting the ''Moving Frame'' option to ''Yes'' will cause the program to solve the equation in a moving frame of reference. ',...
    'This may significantly reduce the computational cost by reducing the length of the channel and the number of grid points.'],...
    '',...
    ['The moving frame speed is determined by right boundary of the first constituent in the input list. For example, for an ITP simulation ',...
    'the leading electrolyte should be set on the right side of the channel and appear first in the input list. This way, the solution will be obtained ',...
    'in a frame of reference moving with the plug']};
set(handles.MessageBoxEdit,'String',msg);set(handles.StartButton,'enable','off');

function SteadyStatePopup_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Yes|No');

function SpatialDiscPopup_Callback(hObject, eventdata, handles)
msg={['SLIP (and 1st order): these are dissipative schemes that are non oscillatory, very robust, and can tolerate much higher grid spacing and remain stable.',...
    ' These schemes can even be used without the adaptive grid. However these schemes may result in diffused interfaces.'],...
    ['SLIP is always preferable compared to 1st order. Using SLIP together with the adaptive grid can yield quite accurate results even for very high current densities.',...
    ' SLIP scheme is recommended for fast and robust simulations for routine optimization of electrophoresis assays.'],...
    'The current implementation of SLIP scheme also allows axial variation in channel cross-sectional area.',...
    '',...
    ['6th order compact (and 2nd order centered): these are non-dissipative schems that are highly accurate but may exhibit oscillations',...
    ' if there is an insufficient number of grid points. 6th order is always preferable compared to the 2nd order, as it will be able to sustain higher current densities without developing oscillations.'],...
    '6th order compact scheme is recommended for simulations where the accuracy of zone boundaries is important.', ...
    'If 6th order compact scheme yields ocsillatory solutions, either increase grid density or switch to the SLIP scheme.'};
set(handles.MessageBoxEdit,'String',msg);

set(handles.StartButton,'enable','off');

function SpatialDiscPopup_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','6th order compact|2nd order centered|1st order Upwind|SLIP');

function AdaptiveGridSpeedEdit_Callback(hObject, eventdata, handles)
msg={['The ''Adaptive Grid Speed'' field determines the rate at which grid points cluster in regions of high gradients.',...
     ' The higher the value, the faster the grid clustering occures, which could be useful in cases where sharp gradients rapidly evolve.',...
     ' However, values that are too high may cause the simulation to run slower (when the adaptive grid equation constrains the time step, instead of the advection-diffusion equations.'],...
     '',...
     'Typical values are between 0.1 and 1.  0.1 is recommened for CE type simulations (shallow gradients), and 1 is recommended for ITP type simulations (sharp gradients).',...
     'Set the adaptive grid speed to 0 to disable the adaptive grid.'};
set(handles.MessageBoxEdit,'String',msg);

set(handles.StartButton,'enable','off');

function AdaptiveGridSpeedEdit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','1');

function AdaptiveGridInterfaceNodesEdit_Callback(hObject, eventdata, handles)
msg={['The ''Clustering Level'' field determines the grid density at regions of high gradients.',...
     'This will influence the adaptive grid only when the grid speed is non zero. The higher this value is, ',...
     'the higher the grid point density would be. However, high grid density results in smaller time steps and a slower simulation.',...
     ' Increase this value only when required to overcome oscillations.'],...
     '',...
     ['Typical values are between 0.1 and 5. 0.1 is typically used for CE simulation, 1 for ITP simulations and 5 for high current density ITP simulations. ',...
     'If it is not possible to obtain a smooth solution using this range,',...
     ' it is usallly more efficient to set it back to 1 and increase the number of grid points.']};
set(handles.MessageBoxEdit,'String',msg);

set(handles.StartButton,'enable','off');

function AdaptiveGridInterfaceNodesEdit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','1');

function PlotIntervalEdit_Callback(hObject, eventdata, handles)
msg={['The ''Plot interval'' field determines the number of time steps between updates of figure.',...
     ' It also determines which times steps are saved to the mat file (every display to the figure is also saved to the file)']};
set(handles.MessageBoxEdit,'String',msg);
set(handles.StartButton,'enable','off');

function PlotIntervalEdit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String', '100');

function CurrentEdit_Callback(hObject, eventdata, handles)
msg={'The ''Current'' field sets the electric current running through the channel.',...
     'This value, together with the cross section definitions, is used to calculated the current density.',...
     '',...
     'Higher current densities lead to sharper concentration gradients which requires more aggresive grid adaptation or an increase in the number of grid points.',...
     'Using a dissipative scheme (1st order or SLIP) makes the solution less sensitive to current density, at the cost of boundary thinkness inaccuracies.'};
set(handles.MessageBoxEdit,'String',msg);
set(handles.StartButton,'enable','off');

function CurrentEdit_CreateFcn(hObject, eventdata, handles)

function DomainLengthEdit_Callback(hObject, eventdata, handles)
msg={'''Domain Length'' sets the length of the simulated channel.',...
    'Generally, the longer the channel the more grid points and therefore computational time are required.',...
    'When possible, use the ''Moving Frame'' option to allow a solution in a frame of reference moving the the area of interest.'};
set(handles.MessageBoxEdit,'String',msg);
if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('PrepareInput',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('PrepareInput',handles);
end

set(handles.StartButton,'enable','off');

function DomainLengthEdit_CreateFcn(hObject, eventdata, handles)

function ChannelShapePopup_Callback(hObject, eventdata, handles)
msg={'For a cirucular cross section, Dim 1 should be set as the diameter.',...
    'For a D shaped cross section, Dim 1 is the depth.',...
    ['Cross section definitions are used to calculate the cross section area and the resulting current density.'...
    'Dim 1. is also used as characteristic scales for calculations of pressure driven flow and dispersion.']};
set(handles.MessageBoxEdit,'String',msg);
set(handles.StartButton,'enable','off');

function ChannelShapePopup_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Circular|D Shape');

function ChannelDimension1Edit_Callback(hObject, eventdata, handles)
msg={'For a cirucular cross section, Dim 1 should be set as the diameter.',...
    'For a D shaped cross section, Dim 1 is the depth.',...
    ['Cross section definitions are used to calculate the cross section area and the resulting current density.'...
    'Dim 1. is also used as characteristic scales for calculations of pressure driven flow and dispersion.']};
set(handles.MessageBoxEdit,'String',msg);

set(handles.StartButton,'enable','off');


function ChannelDimension1Edit_CreateFcn(hObject, eventdata, handles)

function ChannelDimension2Edit_Callback(hObject, eventdata, handles)
msg={'For a cirucular cross section, Dim 1 should be set as the diameter.',...
    'For a D shaped cross section, Dim 1 is the depth.',...
    ['Cross section definitions are used to calculate the cross section area and the resulting current density.'...
    'Dim 1. is also used as characteristic scales for calculations of pressure driven flow and dispersion.']};
set(handles.MessageBoxEdit,'String',msg);
set(handles.StartButton,'enable','off');

function ChannelDimension2Edit_CreateFcn(hObject, eventdata, handles)

function DatabaseFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to DatabaseFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function DatabaseFileEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BrowseDatabasePush_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseDatabasePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function InjectionPointEdit_Callback(hObject, eventdata, handles)
if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('PrepareInput',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('PrepareInput',handles);
end

set(handles.StartButton,'enable','off');

function InjectionPointEdit_CreateFcn(hObject, eventdata, handles)

function InjectionWidthEdit_Callback(hObject, eventdata, handles)
if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('PrepareInput',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('PrepareInput',handles);
end

set(handles.StartButton,'enable','off');

function InjectionWidthEdit_CreateFcn(hObject, eventdata, handles)

function GridPointsEdit_Callback(hObject, eventdata, handles)
msg={'Use the ''Grid Points'' field to set the number of grid points used in the simulation.',...
    'Typical numbers range between 100 and 500 grid points, depending on the length of the channel, the steepness of the gradients and the chosen spatial discretization scheme.',...
    'When using dissipative schemes (1st order, SLIP) you may use may use less grid points, at the cost of accuracy, as these schemes will remain non-oscillatory.',...
    'When using the accurate high resolution 6th order scheme you must ensure there are engough grid points to avoid oscillations.',...
    'See more information on clicking the adaptive grid paramters.'};
set(handles.MessageBoxEdit,'String',msg);

if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('PrepareInput',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('PrepareInput',handles);
end

set(handles.StartButton,'enable','off');




function GridPointsEdit_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in PauseButton.

function PauseButton_Callback(hObject, eventdata, handles)
% Callback not used except for message. Button is read from main program
msg={'The ''Pause'' toggle button pauses the simulation until the button is released.',...
    ' The simulation does not stop, and no data file is written',...
    '',...
    'The option is useful for examining the current state of the simulation, and for readjusting the axes without contant updates from the simulation.'};
set(handles.MessageBoxEdit,'String',msg);

function StopSaveButton_Callback(hObject, eventdata, handles)
% Callback not used except for message. Button is read from main program
msg={'Pressing the ''Stop/Save'' button will stop the simulation and will save the simulation data as a .mat file with the same name as the input file.',...
    '',...
    'The mat file includes the grid spacing and the concentration profiles for all previous times and may be used for postprocessing of the data.',...
    ''...
    'The mat file can also be used to resume the simulation form the time it was last stopped by checking the ''Load mat'' box before pressing ''Start''.'};
set(handles.MessageBoxEdit,'String',msg);



function WhatToPlotPopup_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Concentration|pH|Effective Mobility|Electric Field|Conductivity|Grid Density|Area Variation');

function FilenameEdit_Callback(hObject, eventdata, handles)

function FilenameEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function InputTable_CellEditCallback(hObject, eventdata, handles)
if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('PrepareInput',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('PrepareInput',handles);
end

set(handles.StartButton,'enable','off');

function WhatToPlotPopup_Callback(hObject, eventdata, handles)

if get(hObject, 'UserData')==1
    return;
end

if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('ChangePlotVariable',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('ChangePlotVariable',handles);
end


% if ~get(hObject, 'UserData')
%      disp('I was called');
%      Spresso('ChangePlotVariable',handles);
%  end
%  disp('WhatToPlotPopup was called');
%Spresso('ChangePlotVariable',handles);
%if u r using ChangePlotVariable then modify it

function PrepareGridPopup_Callback(hObject, eventdata, handles)
msg={['Setting ''Prepare grid'' to ''Yes'' will cause the program to start the simulation by clustering grid points at regions of high gradient',...
    'before activating any physical mechanisms such as electromigration or diffusion.'],...
    '',...
    ['This is useful for cases where high gradients are expected. Clustering grid points in advance usualy enables to',...
    'successfully simulate higher current densities without running into oscillations'],...
    '',...
    ['For cases with shallower gradients the grid adaptation during the simulation is usualy sufficient, and you ',...
    'may turn this option off to save some computational time at the begining of the simulation.']};
set(handles.MessageBoxEdit,'String',msg);

set(handles.StartButton,'enable','off');

function PrepareGridPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','Yes|No');

function SaveInputFile(handles,Filename)
msgboxFlag=0;

IonicEffectFlag=2-get(handles.IonicEffectPopup,'value');            if isempty(IonicEffectFlag), msgboxFlag=1;  end
PcentIonicChange = str2num (get(handles.PcentIonicChangeEdit,'string'));    if isempty(PcentIonicChange), msgboxFlag=1;  end
tEnd = str2num (get(handles.SimulationTimeEdit,'string'));        if isempty(tEnd), msgboxFlag=1;  end
N = str2num (get(handles.GridPointsEdit,'string'));               if isempty(N), msgboxFlag=1;  end
L = str2num (get(handles.DomainLengthEdit,'string'));             if isempty(L), msgboxFlag=1;  end
InjLen = str2num (get(handles.InjectionWidthEdit,'string'));      if isempty(InjLen), msgboxFlag=1;  end
InjLoc = str2num (get(handles.InjectionPointEdit,'string'));      if isempty(InjLoc), msgboxFlag=1;  end
DChannel = str2num (get(handles.ChannelDimension2Edit,'string')); if isempty(DChannel), msgboxFlag=1;  end
hChannel = str2num (get(handles.ChannelDimension1Edit,'string')); if isempty(hChannel), msgboxFlag=1;  end
ChannelShape=get(handles.ChannelShapePopup,'Value');              if isempty(ChannelShape), msgboxFlag=1;  end
Current = str2num (get(handles.CurrentEdit,'string'));            if isempty(Current), msgboxFlag=1;  end

SteadyStateFlag=2-get(handles.SteadyStatePopup,'value');            if isempty(SteadyStateFlag), msgboxFlag=1;  end
SpatialDiscFlag=get(handles.SpatialDiscPopup,'value');            if isempty(SpatialDiscFlag), msgboxFlag=1;  end
PrepareGridStage=2-get(handles.PrepareGridPopup,'value');           if isempty(PrepareGridStage), msgboxFlag=1;  end
AdaptGrid.Coeff=str2num(get(handles.AdaptiveGridSpeedEdit,'string'));       if isempty(AdaptGrid.Coeff), msgboxFlag=1;  end
AdaptGrid.PointsAtInterface=str2num(get(handles.AdaptiveGridInterfaceNodesEdit,'string'));  if isempty(AdaptGrid.PointsAtInterface), msgboxFlag=1;  end
DeltaCounterNextPlot=str2num(get(handles.PlotIntervalEdit,'string'));       if isempty(DeltaCounterNextPlot), msgboxFlag=1;  end

bPressure = str2num (get(handles.bPressureEdit,'string'));            if isempty(bPressure), msgboxFlag=1;  end
Pressurehead = str2num (get(handles.PressureheadEdit,'string'));            if isempty(Pressurehead), msgboxFlag=1;  end
betaDispersion = str2num (get(handles.betaDispersionEdit,'string'));            if isempty(betaDispersion), msgboxFlag=1;  end

InputTable=cell(handles.InputTable.getData);
InputTable=InputTable(:,2:end);

AreaFunctionText = get(handles.AreaVariation_Edit, 'String'); if isempty(AreaFunctionText), msgboxFlag=1;  end

if msgboxFlag==1
    msgbox('Please input all parameters');
    return;
end

fid=fopen(Filename,'w');
fprintf(fid,['%%INPUT FILE \n\n']);
fprintf(fid,'IonicEffectFlag = %d;  %%IonicStrengthFlag on/off \n', IonicEffectFlag);
fprintf(fid,'PcentIonicChange = %f;  %% change in IonicStrength before next evaluation \n', PcentIonicChange);
fprintf(fid,'L1 = 0;   \n');
fprintf(fid,'L2 = %e;  %% Domain length [m]\n',L*1E-3);
fprintf(fid,'N = %d;  %% Number of grid points\n',N);
fprintf(fid,'DChannel = %e;  %% Channel width\n',DChannel*1E-6);
fprintf(fid,'hChannel = %e;  %% Channel depth\n',hChannel*1E-6);
fprintf(fid,'ChannelShape = %d;  %% Channel Shape 1-circular 2-D shaped\n',ChannelShape);


fprintf(fid,'Current = %fE-6;  %% Current [A]\n',Current);
fprintf(fid,'tEnd = %f; %% [sec] Stop simulation at this time \n',tEnd);

fprintf(fid,'\n\nSteadyStateFlag = %d;  %% Solve in frame of reference moving with LE velocity\n',SteadyStateFlag);
fprintf(fid,'PrepareGridStage = %d; %% Enable/disable grid adaptation stage\n',PrepareGridStage);
fprintf(fid,'AdaptGrid.Coeff = %f; %% Adaptive grid speed\n',AdaptGrid.Coeff);
fprintf(fid,'AdaptGrid.PointsAtInterface = %d; %% Clustering level of grid points\n',AdaptGrid.PointsAtInterface);
fprintf(fid,'DeltaCounterNextPlot = %d; %% Number of iterations before plotting\n',DeltaCounterNextPlot);

switch(SpatialDiscFlag)
    case(1)
        fprintf(fid,'SpatialDiscFlag = ''%s''; %% Spatial discretization scheme\n','Compact');
    case(2)
        fprintf(fid,'SpatialDiscFlag = ''%s''; %% Spatial discretization scheme\n','2nd');
    case(3)
        fprintf(fid,'SpatialDiscFlag = ''%s''; %% Spatial discretization scheme\n','Upwind');
    case(4)
        fprintf(fid,'SpatialDiscFlag = ''%s''; %% Spatial discretization scheme\n','SLIP');
end

fprintf(fid,'\nInjLen = %f; %% Injection width [mm]\n',InjLen);
fprintf(fid,'InjLoc = %f; %% Injection location [mm]',InjLoc);

fprintf(fid,'\nPressurehead = %f; %% Pressure head [mm]\n',Pressurehead);
fprintf(fid,'bPressure = %f; %% Hydraulic resistance coefficient\n',bPressure);
fprintf(fid,'betaDispersion = %f; %% Taylor-Aris Dispersion coefficient\n',betaDispersion);


fprintf(fid,'\n\n%% Vectors defining distribution of analytes\n');
fprintf(fid,'zVec=linspace(L1,L2,N)'';  %% Define grid\n');
fprintf(fid,'phiVec=zVec; xVec=phiVec;  %% Start with uniform grid\n');
fprintf(fid,'AnalyteVec=(1+erf((xVec/1E-3-(InjLoc-InjLen/2))))/2-(1+erf((xVec/1E-3-(InjLoc+InjLen/2))))/2;   AnalyteVec(abs(AnalyteVec)<1E-3)=0;\n');
fprintf(fid,'TEVec=(1-erf((xVec/1E-3-(InjLoc-InjLen/2))))/2;    TEVec(abs(TEVec)<1E-3)=0;\n');
fprintf(fid,'LEVec=(1+erf((xVec/1E-3-(InjLoc+InjLen/2))))/2;    LEVec(abs(LEVec)<1E-3)=0;\n');
fprintf(fid,'BackgroundVec=xVec*0+1;\n');

fprintf(fid,'\n\n%% Variation in Cross-section: Area normalized by left most Area\n');


fprintf(fid,'AreaFunctionText = ''%s''; %% Function defining area variation\n', AreaFunctionText);
fprintf(fid,'A = feval( eval(AreaFunctionText), phiVec);\n');
fprintf(fid,'AreaRatio=A./max(A);\n');

fprintf(fid,'\n\n%%Analytes chemistry\nInputTable={');
for kl=1:size(InputTable,1)
    str=[];
    for ij=1:4
        object=InputTable{kl,ij};
        if isnumeric(object),
            str=[str,num2str(object)];
        else
            str=[str,'''',object,''''];
        end
        if ij~=4, str=[str,','];
        else
            str=[str,';\n\t\t\t'];
        end
    end
    fprintf(fid,str);
end % for kl
fprintf(fid,'};\n');

fprintf(fid,'\n%% Analyte distributions\n');
for ij=1:size(InputTable,1)
    fprintf(fid,['cMat(:,',num2str(ij),') = InputTable{',num2str(ij),',3} * ',InputTable{ij,4},'Vec;\n']);
end % for ij
fclose(fid);

function handles=LoadFilePush_Callback(hObject, eventdata, handles)
msg={'''Load file'' reads an existing input file and fills all the GUI boxes with that information.',...
    '',...
    'If a corresponding .mat file is found, you will be able to choose whether to load it or not.',...
    'Choosing ''Yes'' will automatically check the ''Loat mat'' box so that the data will be loaded after hitting ''Start''.',...
    'You may still ignore the data by unchecking the box, in which case only the input file will be processed after hitting ''Start''.',...
    '',...
    'You may also make manual changes to the input file (for example, for defining complex initial concentration profiles'};

set(handles.MessageBoxEdit,'String',msg);

[filename, pathname] = uigetfile({'*.m', 'MATLAB Files (*.m)';'*.*','All Files (*.*)'},'Load Input File');
if filename==0
    return;
end
LoadFile([pathname,filename],handles);
set(handles.FilenameEdit,'String',[pathname,filename]);

function ReloadFilePush_Callback(hObject, eventdata, handles)
msg={'The ''Reload'' button simply reloads the input file specified in the ''Filename'' field.',...
    'This is useful when making frequent manual changes to the input file (for example, for defining complex intitial concentration profiles)',...
    'Note that the ''Save'' button should not be used when manualy modifying the input file, or the file will be overwritten by the GUI content.' };
set(handles.MessageBoxEdit,'String',msg);

Filename=get(handles.FilenameEdit,'String');
LoadFile(Filename,handles);

function LoadFile(Filename,handles)
global ColorList
% Evaluate input file
fid=fopen(Filename,'r');  str=[];
while ~feof(fid)
    str=[str,fgets(fid)];
end
eval(str);
%DatFileExist=exist([strtok(Filename,'.'),'.mat']);
DatFileExist=exist([Filename,'at']);

if DatFileExist==2
    LoadDataFlag=questdlg('A corresponding data file has been found. Load data file?','Data file found');
    switch (LoadDataFlag)
        case('Yes')
            pause(0.1); % insert pause to fix questdlg bug
            %load([strtok(Filename,'.'),'.mat']);
            load([Filename,'at']);
            xVec=phiVec;
            set(handles.LoadDataCheckbox,'value',1);
        case('No')
           pause(0.1); % insert pause to fix questdlg bug
           set(handles.LoadDataCheckbox,'value',0);
        case('Cancel')
            return;
    end
end % if

ErrorFlag=0;

set(handles.IonicEffectPopup,'value',2-IonicEffectFlag);
set(handles.PcentIonicChangeEdit,'string',num2str(PcentIonicChange));

if (IonicEffectFlag == 1)
    set(handles.PcentIonicChangeEdit,'enable','on');
elseif (IonicEffectFlag == 2)
    set(handles.PcentIonicChangeEdit,'enable','off');
end

set(handles.DomainLengthEdit,'string',num2str(L2*1E3));
set(handles.GridPointsEdit,'string',num2str(N));
set(handles.ChannelDimension2Edit,'string',num2str(DChannel*1E6));
set(handles.ChannelDimension1Edit,'string',num2str(hChannel*1E6));
set(handles.CurrentEdit,'string',num2str(Current));
set(handles.SteadyStatePopup,'value',2-SteadyStateFlag);
set(handles.PrepareGridPopup,'value',2-PrepareGridStage);
set(handles.AdaptiveGridSpeedEdit,'string',num2str(AdaptGrid.Coeff));
set(handles.AdaptiveGridInterfaceNodesEdit,'string',num2str(AdaptGrid.PointsAtInterface));
set(handles.InjectionWidthEdit,'string',num2str(InjLen));
set(handles.InjectionPointEdit,'string',num2str(InjLoc));
set(handles.CurrentEdit,'string',num2str(Current*1E6));
set(handles.ChannelShapePopup,'value',ChannelShape);
set(handles.PlotIntervalEdit, 'string', num2str(DeltaCounterNextPlot));

if exist('AreaFunctionText', 'var')
set(handles.AreaVariation_Edit, 'String', AreaFunctionText);
else
 set(handles.AreaVariation_Edit, 'String', '@(x) 1 + 0*x');
 ErrorFlag = 1;
 errordlg('Area Variation not defined. Setting automatically to unifom cross-section. Save input file to proceed.');
end

try
set(handles.SimulationTimeEdit,'string',tEnd);
set(handles.PressureheadEdit,'string',Pressurehead);
set(handles.bPressureEdit,'string',bPressure);
set(handles.betaDispersionEdit,'string',betaDispersion);
catch
end

TmpTable=cell(size(InputTable,1),size(InputTable,2)+1);
TmpTable(:,2:end)=InputTable;
InputTable=TmpTable;

for ij=1:size(InputTable,1)
    CurrentColorString=num2str(round(255*ColorList(ij,:)),'%d,');
    InputTable{ij,1}=['<HTML><FONT color=rgb(',CurrentColorString(1:end-1),')>@@@</FONT></HTML>'];
end
% Get current table
CurrentInput=cell(handles.InputTable.getData);
%Remove all rows
for ij=1:size(CurrentInput,1)
    handles.InputTable.getTable.getModel.removeRow(0);
end
% Add new rows
 NewLine=InputTable(1,:);
  for ij=1:size(InputTable,1)
     handles.InputTable.getTable.getModel.addRow(InputTable(ij,:));  % appends a row to the bottom of the table
 end

switch(SpatialDiscFlag)
    case ('Compact');
        val=1;
    case ('2nd');
        val=2;
    case ('Upwind')
        val=3;
    case ('SLIP')
        val=4;
end
set(handles.SpatialDiscPopup,'value',val);
set(handles.DataHolder,'UserData',[xVec,cMat]);
if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('ChangePlotVariable',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('ChangePlotVariable',handles);
end
if ~ErrorFlag
    set(handles.StartButton,'enable','on');
end

function StartButton_Callback(hObject, eventdata, handles)
msg={['When hitting the ''Start'' button the program will process the input file specified in the ''Filename'' field and NOT the content of the GUI',...
     ' Therefore, any changes in the GUI fields disables the ''Start''button until the ''Save'' or ''Save as'' buttons are pressed'],...
     '',...
     'If the ''Load mat'' box is checked, the corresponding data file (.mat) will be loaded after processing the input file'};

set(handles.MessageBoxEdit,'String',msg);

LoadDataFlag=get(handles.LoadDataCheckbox,'value');
%if Ionic strength is on then run Spresso normally. Since directory
%SubroutinesIonic is added to path
%Otherwise remove the directory SubroutinesIonic from path
%and add SubroutinesNoIonic to path
%
if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
% rmpath(['SubroutinesNoIonic']);
% addpath(['SubroutinesIonic']);
SpressoIonic('MainRun',get(handles.FilenameEdit,'String'),LoadDataFlag,handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
% rmpath(['SubroutinesIonic']);
% addpath(['SubroutinesNoIonic']);
SpressoNoIonic('MainRun',get(handles.FilenameEdit,'String'),LoadDataFlag,handles);
end

function SaveButton_Callback(hObject, eventdata, handles)
msg={'Clicking the ''Save'' button saves the content of the GUI to file name specified in the ''Filename'' field.',...
     '',...
     ['When hitting the ''Start'' button the program will read this input file and NOT the content of the GUI',...
     ' Therefore, any changes in the GUI fields disables the ''Start''button until the ''Save'' or ''Save as'' buttons are pressed']};

set(handles.MessageBoxEdit,'String',msg);

Filename=get(handles.FilenameEdit,'String');
if isempty(Filename)
    msgbox({'Please enter a file name in the''Filename'' field','or use the ''Save As'' option'},'No file name specificed','error');
    return;
end
SaveInputFile(handles,Filename);
set(handles.StartButton,'enable','on');

function SaveAsButton_Callback(hObject, eventdata, handles)
msg={'Clicking the ''Save as'' allows you to specify the name and location of a new input file. Following this, the ''Filename'' field will be updated with the full path to the file.',...
     '',...
     ['When hitting the ''Start'' button the program will read this input file and NOT the content of the GUI',...
     ' Therefore, any changes in the GUI fields disables the ''Start''button until the ''Save'' or ''Save as'' buttons are pressed']};

set(handles.MessageBoxEdit,'String',msg);

[filename, pathname] = uiputfile({'*.m', 'MATLAB Files (*.m)';'*.*','All Files (*.*)'},'Load Input File');
if filename==0
    return;
end
set(handles.FilenameEdit,'String',[pathname,filename]);
SaveInputFile(handles,[pathname,filename]);
set(handles.StartButton,'enable','on');



function StartButton_CreateFcn(hObject, eventdata, handles)
set(hObject,'enable','off');

function LoadDataCheckbox_Callback(hObject, eventdata, handles)
msg={['If the ''load mat'' box is checked, after hitting ''Start'' the program will look for an existing .mat file',...
    ' with the same name as the input file and will load it.'],...
    ' ',...
    'This is useful when you wish to resume the simulation from the point it stopped:',...
    '- click the ''Stop/Save'' button when the simulation is running to save all the data',...
    '- Next time you run the simulation check the ''load mat'' box and press Start',...
    ' ',...
    'If a .mat file exists but you do not wish to load, just uncheck the ''load mat'' box and only the input file will be processed'};

set(handles.MessageBoxEdit,'String',msg);

function StopDiscardButton_Callback(hObject, eventdata, handles)

function LogoAxes_CreateFcn(hObject, eventdata, handles)
axes(hObject);
Img = importdata('SpressoLogo.png');
Vec=size(Img);
set(hObject,'units','pixels');
Position=get(hObject,'Position'); Position(3:4)=[Vec(2),Vec(1)];
set(hObject,'Position',Position);
imagesc(Img);axis equal;
col=get(gcf,'color');
set(gca,'yticklabel',[],'xticklabel',[],'xtick',[],'ytick',[],'ycolor',col,'xcolor',col,'box','on');

function MainAxes_CreateFcn(hObject, eventdata, handles)
global ColorList
set(hObject,'ColorOrder',ColorList);


function Panel1Radiobutton_Callback(hObject, eventdata, handles)
set(handles.Panel2,'visible','off');
set(handles.Panel1,'visible','on');

function Panel2Radiobutton_Callback(hObject, eventdata, handles)
set(handles.Panel1,'visible','off');
set(handles.Panel2,'visible','on');

function SimulationTimeEdit_Callback(hObject, eventdata, handles)
set(handles.StartButton,'enable','off');

function SimulationTimeEdit_CreateFcn(hObject, eventdata, handles)

function betaDispersionEdit_Callback(hObject, eventdata, handles)
msg={'The ''beta'' parameter determines the Taylor-Aris dispersion coeffieicent in the form',...
    'Deff = D * (1 + beta * Pe^2) ',...
    'where Pe is based on the parameter Dim1 (diameter for circular channel, depth for D shaped channels).'};
set(handles.MessageBoxEdit,'String',msg);set(handles.StartButton,'enable','off');
set(handles.StartButton,'enable','off');

function betaDispersionEdit_CreateFcn(hObject, eventdata, handles)

function PressureheadEdit_Callback(hObject, eventdata, handles)
msg={['The ''Pressure head'' parameter determines the pressure (in mm H2O) applied to the channel. ',...
    'This can be used, for example, for simulating sationary ITP. A positive value corresponds to pressure ',...
    'applied to the right well (right boundary).']};
set(handles.MessageBoxEdit,'String',msg);set(handles.StartButton,'enable','off');
set(handles.StartButton,'enable','off');

function PressureheadEdit_CreateFcn(hObject, eventdata, handles)

function bPressureEdit_Callback(hObject, eventdata, handles)
msg={'The ''b'' parameter hydraulic resistance coeffieicent for calculating pressure driven flow. The velocity is then:',...
    'u = - DP * d^2 / (b * mu * L)  ',...
    ['where d is the parameter Dim1 (diameter for circular channel, depth for D shaped channels).',...
    ' DP is the pressure difference across the channel, L is the length of the channel, and mu is the dynamic viscosity']};
set(handles.MessageBoxEdit,'String',msg);set(handles.StartButton,'enable','off');
set(handles.StartButton,'enable','off');

function bPressureEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MessageBoxEdit_CreateFcn(hObject, eventdata, handles)
msg={'Float over the different input field to get fast tooltips.',...
    '',...
    ['Click and change to values of the field to get more ellaborated information in this box.',...
    ' To change the value of an edit box after entering a new value, you must either hit "Enter" or click outside the box.']};
set(hObject,'String',msg);


function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);


% --------------------------------------------------------------------
function AboutSpressoPush_ClickedCallback(hObject, eventdata, handles)
Data=get(hObject,'CData');
Data(isnan(Data))=1;
MsgString = {'\bfSpresso v0.1',...
             'Stanford''s public release electrophoretic separation solver',...
             '',...
             '\rmCode by Moran Bercovici, Supreet Singh Bahga, Sanjiva K. Lele and Juan G. Santiago',...
             'Stanford Microlfuidics Laboratory - http://microfluidics.stanford.edu',...
             '',...
             'Database largely based on publications by T. Hirokawa, Journal of Chromatography',...
             '',...
             'Java GUI tables by Yair Altman, http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=14225&objectType=FILE',...
             '',...
             'Free for use for all purposes',...
             '---------------------------------------------',...
             'The code is provided as an open source so that you may modify, add, remove, embed, or link part or all of the codes with additional components and other programs.',...
             '',...
             'Disclaimer',...
             '-----------------',...
             'Use at your own risk - the code is provided on an "as is" basis without warranty of any kind, express or implied.' };
%MsgString=textwrap(MsgString, 45);
CreateStruct.Interpreter='tex';
CreateStruct.WindowStyle='non-modal';
h=msgbox(MsgString,'About Spresso','custom',Data,jet(64),CreateStruct);


% --- Executes on selection change in IonicEffectPopup.
function IonicEffectPopup_Callback(hObject, eventdata, handles)
% hObject    handle to IonicEffectPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns IonicEffectPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IonicEffectPopup
msg={['Setting the ''Ionic Strength Dependence'' option to ''Yes'' will cause the program to used ionic strength effects on mobility, pKa and diffusivity.',...
    '   This may reduce the speed by factor of 2-3 due to additonal computations for ionic strength']};
set(handles.MessageBoxEdit,'String',msg);
set(handles.StartButton,'enable','off');
if get(hObject, 'UserData')==1
    return;
end

%change the paths accordingly depending upon choice of ionic strength
%on/off
    if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
    SpressoIonic('ChangePlotVariable',handles); %testing
    set(handles.PcentIonicChangeEdit,'enable','on');
        %         rmpath(['SubroutinesNoIonic']);
%         addpath(['SubroutinesIonic']);
    elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
        SpressoNoIonic('ChangePlotVariable',handles); %testing
        set(handles.PcentIonicChangeEdit,'enable','off');
%         rmpath(['SubroutinesIonic']);
%        addpath(['SubroutinesNoIonic']);
    end



% --- Executes during object creation, after setting all properties.
function IonicEffectPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IonicEffectPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','Yes|No');
set(hObject,'value',2);



function PcentIonicChangeEdit_Callback(hObject, eventdata, handles)
msg={'Use the ''% change of Ionic Strength'' field to set the ionic strength difference after which the program solves for ionic strength',...
    'This would range between 0-100',...
    'Choice 0 means that ionic strength calculation is done everytime step',...
    'Typically 1% would work well.'};
set(handles.MessageBoxEdit,'String',msg);

if get(handles.IonicEffectPopup,'value')==1 %ionic strength is On
SpressoIonic('PrepareInput',handles);
elseif get(handles.IonicEffectPopup,'value')==2 %ionic strength is Off
SpressoNoIonic('PrepareInput',handles);
end
set(handles.StartButton,'enable','off');


% --- Executes during object creation, after setting all properties.
function PcentIonicChangeEdit_CreateFcn(hObject, eventdata, handles)
 set(hObject,'enable','off');


function AreaVariation_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to AreaVariation_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AreaVariation_Edit as text
%        str2double(get(hObject,'String')) returns contents of AreaVariation_Edit as a double

msg={'Use the Area Variation field to set axial variation of channel cross-section',...
    'Available with SLIP scheme only', ...
    'E.g.: Uniform cross-section: @(x) 1+0*x',...
    'For 1:10 area reduction: @(x) (1+1/10)/2-(1-1/10)/2*tanh(100*(x-20e-3))',...
    'Use only vectorized notation, e.g. 1+x.^2',...
    'Cross-sectional variation is automatically normalized to the area', ...
    'of inlet (left side) of the domain ', ...
    'defined by fields Dim1 and Dim2'};

set(handles.MessageBoxEdit,'String',msg);
set(handles.StartButton,'enable','off');
SpressoNoIonic('ChangePlotVariable',handles); %cross-section is same with or without ionic strength

% --- Executes during object creation, after setting all properties.
function AreaVariation_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AreaVariation_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

set(hObject, 'String', '@(x) 1 + 0*x');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
