function [Kd1,k1_,Kd2,k2_,kA,Kdlux,klux_,kR,kI,        ...
          pR,pI,alphaR,alphaI,dA,dC2,dC,dR,dI,dmR,dmI, ...
          D,dilution,growth,                           ...
          diffusiveLossRate,colonyD]= modelCoeffs_weber(y,growthOn,diffusionOn)
%% Constants
% The kinetic constants are in (nM, 1/min)
% Space is in mm
Kd1=    100;   %[nM]
k1_=    10;    %[1/min]
Kd2=    20;    %[nM]
k2_=    1;     %[1/min]
kA=     .04;   %[1/min]
Kdlux=  200;   %[nM]
klux_=  10;    %[1/min]
b=      20;
kR=     200/b; %[1/min]     % transcription rates "fitted" (revise?)
kI=     50/b;  %[1/min]
dmR=    0.347; %[1/min]
dmI=    0.347; %[1/min]
pR=     b*dmR; %[1/min]
pI=     b*dmI; %[1/min]
alphaR= .001;
alphaI= .01;
dA=     .001;  %[1/min]
dC2=    .002;  %[1/min]
dC=     .002;  %[1/min]
dR=     .002;  %[1/min]
dI=     .01;   %[1/min]
D=      10;    % membrane permeability
taf=    45;    %[min]
%lamda=  .8;
%V0=     1.5*1e-9;   %[mm3]
%Vtot=   200*1e-9;   %[mm3]

% Perhaps the translation rate should be scaled by <r_current/r_paper>, as ribosome count
% is directly related to division rate? (metabolic capacity model)
r_paper= log(2)/taf;
V0_mean= 2.164e-9;    %[mm3]   average cell size across its cell cycle
% how closely together the bacteria are packed (1.05 -> separation: 5% of length per bacterium)
packingFactor= 1.1;
colonyHalo= 2;

% Extra parameters, not from QS paper
cellDNA= 1.5347;   % amount of DNA per cell

% Growth model constants (logistic++)
growth.r= 2.20/60;              % /60: time multiplier (hour->min)
growth.m= 0.52; growth.n= 3.5;
growth.Nmax= 10^(6);
%growth.Nmin= N0*(1-1e-6);
growth.Nmin= 0;   % disregard lag phase

%% Coefficients dependent on y
% Diffusion-related
D= diffusionOn*D;     % enable/disable internal diffusion, to accomodate pde
Vcelltot= V0_mean*y(11);
colonyD= D;%*y(11);
% total colony size
Vcolony= Vcelltot*packingFactor^3*colonyHalo^3;
colonyArea= Vcolony^(2/3)*5;
diffusiveLossRate= 0.1*colonyArea; % [1/min]
% dilution: <r> in paper
dilution= Vcelltot/(Vcolony-Vcelltot);

%% Growth model
growth.p1= 1-(y(11)./growth.Nmax).^growth.m;
growth.p2= 1-(growth.Nmin./y(11)).^growth.n;
growth.dN= growthOn*growth.r*y(11)*growth.p1*growth.p2;
growth.divRate= growth.dN ./ y(11);
growth.rateScaling= growth.divRate ./ r_paper;

%% Metabolic capacity model
% Simple, hypothetical model
%metabolicCapacity= max(0.6,growth.rateScaling);
%pI= pI * metabolicCapacity; pR= pR * metabolicCapacity;
