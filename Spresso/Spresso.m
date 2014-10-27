function Spresso(varargin)
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
muH=362E-9/F;           % Mobility of Hydronium   %% [m^2/s*V]/F --> [mol*s/Kg]
muOH=205E-9/F;          % Mobility of Hydroxide    [m^2/s*V]/F --> [mol*s/Kg]
met2lit=1000;           % m^3 to liters
visc=1E-3;              % Dynamic viscosity (water) [Pa s]

% If no input arguments - run graphical user interface


if isempty(varargin)
   clc;
    warning off all
    close all;
    handles = SpressoGUI;
    return;
elseif nargin && ischar(varargin{1})
    feval(str2func(varargin{1}),varargin{2:end});
end
end
