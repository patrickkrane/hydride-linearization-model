%DefineReactorProperties.m
%

%This function sets the values for the dimensions of the metal hydride
%reactor and the constant properties of the coolant and hydrogen.

%DefineInitalReactorProperties.m requires:=
% The Purdue Metal Hydride Toolbox

%@author: Patrick Krane      email: pkrane@purdue.edu
function [Reactors,InProps] = DefineReactorProperties(Reactors)
%% Variable Definition
%Defining hydride materials
Reactors.hydA      = MetalHydride('Mm[1]Ni[4.5]Cr[0.5].mh'); %Low-temp side hydride Nd[1]Ni[4.8]Sn[0.2]
Reactors.hydB      = MetalHydride('La[1]Ni[5].mh'); %High-temp side hydride
%Defining dimensions of the hydride reactors based on the predefined number of tubes
Reactors.LA        = 1.7727; %Length of reactor A
Reactors.LB        = 1.5363; %Length of reactor B (difference in values is based on different storage capacities of the materials)
Reactors.eps       = 0.2; %Porosity of the hydride bed
Reactors.Ri        = 0.004; %Inner Radius of the reactor/ radius of coolant tube (m)
Reactors.Ro        = 0.007; %Outer Radius of the reactor (m)
Reactors.Asec      = pi * Reactors.Ro ^ 2; %Cross-sectional area of the control volume (one tube and surrounding shell) used in the model
Reactors.Atot      = Reactors.ntube * Reactors.Asec; %Cross-sectional area of the reactor
Asec_hyd           = pi * (Reactors.Ro ^ 2 - Reactors.Ri ^2); %Cross-sectional area of the metal hydride shell in the control volume
Reactors.Rtot      = sqrt(Reactors.Atot / pi); %Reactor radius
Reactors.VsecA     = Asec_hyd * Reactors.LA; %Total volume of the control volume in reactor A
Reactors.VsecB     = Asec_hyd * Reactors.LB; %Total volume of the control volume in reactor B
Reactors.V_rA      = Reactors.VsecA * Reactors.ntube; %Total volume of  hydride reactor A (m^3)
Reactors.V_rB      = Reactors.VsecB * Reactors.ntube; %Total volume of  hydride reactor B (m^3)
Reactors.mA        = Reactors.V_rA * (Reactors.hydA.rho * (1- Reactors.eps)); %Mass of hydride in reactor A
Reactors.mB        = Reactors.V_rB * (Reactors.hydB.rho * (1- Reactors.eps)); %Mass of hydride in reactor B
Reactors.mA_sec    = Reactors.VsecA * (Reactors.hydA.rho * (1- Reactors.eps)); %Mass of hydride in the control volume in reactor A
Reactors.mB_sec    = Reactors.VsecB * (Reactors.hydB.rho * (1- Reactors.eps)); %Mass of hydride in the control volume in reactor A
%Defining coolant properties for 47 % Propylenglycol water
InProps.cp_c       = 3550; %Specific Heat (J/kg K)- could make this 3620
InProps.rho_c      = 1044; %Density (kg/m^3)
InProps.kcool      = 0.394; %Thermal conductivity (W/mK)
%Defining material properties of hydrogen
MMH                = 0.00202; %Molar mass of H2, in kg/mol
InProps.Rbar       = 8.314; %J/mol*K, universal gas constant
InProps.cp_H	   = 14307; %J/(kg * K)
InProps.R_H	       = InProps.Rbar / MMH; %J/kg*K, gas constant for hydrogen
InProps.mu_H       = 0.88e-5; %Viscosity of hydrogen, approximate for this temp range (Pa*s)
InProps.rho_H_1b   = 0.085; %Density of H at 1 bar, approximate value for this tempertaure range(kg/m^3)
InProps.Patm       = 1.01325e5; %Atmospheric pressure
%Defining properties of the hydrogen line between the reactors
Reactors.D_Hline   = 0.01; %Diameter of hydrogen pipe
Reactors.Kloss_H   = 40; %Loss coefficient
Reactors.comp_eff  = 0.8;%Compressor efficiency