function [p,growth]= modelCoeffs_weber(y,N0,growthOn,diffusionOn)
  %[Kd1,k1_,Kd2,k2_,kA,Kdlux,klux_,kR,kI,        ...
  %        pR,pI,alphaR,alphaI,dA,dC2,dC,dR,dI,dmR,dmI, ...
  %        D,dilution,growth,                           ...
  %        colonyD]= modelCoeffs_weber(y,growthOn,diffusionOn)
%% Constants
% The kinetic constants are in (nM, 1/min)
% Space is in mm
p.Kd1=    100;   %[nM]
p.k1_=    10;    %[1/min]
p.Kd2=    20;    %[nM]
p.k2_=    1;     %[1/min]
p.kA=     .04;   %[1/min]
p.Kdlux=  200;   %[nM]
p.klux_=  10;    %[1/min]
b=        20;
p.kR=     200/b; %[1/min]     % transcription rates "fitted" (revise?)
p.kI=     50/b;  %[1/min]
p.dmR=    0.347; %[1/min]
p.dmI=    0.347; %[1/min]
p.pR=     b*p.dmR; %[1/min]
p.pI=     b*p.dmI; %[1/min]
p.alphaR= .001;
p.alphaI= .01;
p.dA=     .001;  %[1/min]
p.dC2=    .002;  %[1/min]
p.dC=     .002;  %[1/min]
p.dR=     .002;  %[1/min]
p.dI=     .01;   %[1/min]
p.D=      10;    % membrane permeability
taf=      45;    %[min]
%lamda=  .8;
%V0=     1.5*1e-9;   %[mm3]
Vtot=     5000;   %[mm3]

% Perhaps the translation rate should be scaled by <r_current/r_paper>, as ribosome count
% is directly related to division rate? (metabolic capacity model)
r_paper= log(2)/taf;
V0_mean= 1.6996e-9;    %[mm3]   average cell size across its cell cycle
% how closely together the bacteria are packed (1.05 -> separation: 5% of length per bacterium)
%packingFactor= 1.1;
%colonyHalo= 2;

% Extra parameters, not from QS paper
cellDNA= 1.5347;   % amount of DNA per cell

% Growth model constants (logistic++)
%growth.r= r_paper;              % /60: time multiplier (hour->min)
growth.r= 1.5/60;
growth.m= 0.52; growth.n= 3.5;
growth.Nmax= 1*10^(8.9);
%growth.Nmin= N0*(1-1e-6);
growth.Nmin= 0;   % disregard lag phase

%% Coefficients dependent on y
% Diffusion-related
p.D= diffusionOn*p.D;     % enable/disable internal diffusion, to accomodate pde
%Vcelltot= V0_mean*y(11);
%p.colonyD= p.D*y(11);

%Vcelltot= V0_mean*growth.Nmax*Vtot/1e3;
%p.colonyD= p.D*growth.Nmax;

Vcelltot= V0_mean*y(11);
p.colonyD= p.D;
%Vcolony= 4.5889e-8;
Vcolony= Vtot;

% total colony size
%Vcolony= Vcelltot*packingFactor^3*colonyHalo^3;
%Vcolony= Vtot;
% The exterior volume doesn't correspond between the local and the distributed models in terms
% of [AHL]
% dist: 2.94e-06 -> local: 1.2e-06
%Vcolony= 1.2e-06;
% dilution: <r> in paper
p.dilution= Vcelltot/(Vcolony-Vcelltot);

%% Growth model
growth.p1= 1-(y(11)./growth.Nmax).^growth.m;
growth.p2= 1-(growth.Nmin./y(11)).^growth.n;
growth.dN= growthOn*growth.r*y(11)*growth.p1*growth.p2;
if p.dilution > 9
  p.dilution= 9;
  growth.dN= 0;
end
growth.divRate= growth.dN ./ y(11);
growth.rateScaling= growth.divRate ./ r_paper;

  %growth.dN= 0;
%growth.divRate= growth.r;
%% Metabolic capacity model
% Simple, hypothetical model
%metabolicCapacity= max(0.6,growth.rateScaling);
%pI= pI * metabolicCapacity; pR= pR * metabolicCapacity;


