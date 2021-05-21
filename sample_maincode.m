%sample_maincode.m
%

%This script shows how the linearize_model function can be used inside a
%larger script.

%DefineInitalReactorProperties.m requires:=
% The Purdue Metal Hydride Toolbox- can be downloaded at https://github.com/PurdueH2Lab/MetalHydrideToolbox
% DefineInitialReactorProperties.m
% linearize_model.m
% CalculateEffectiveness.m

%@author: Patrick Krane      email: pkrane@purdue.edu

addpath '\\nas01.itap.purdue.edu\puhome\My Documents\GitHub\CHPB-34\chpb-34-metal-hydrides\MetalHydrideToolbox-master'
op_mode        = 1; %1-charge    (H moves from reactor B to reactor A)
                    %0-discharge (H moves from reactor A to reactor B)
Reactors.ntube = 400; %Number of tubes in the shell-and-tube heat exchangers
[Reactors,InProps] = DefineReactorProperties(Reactors); %Defining the dimensions of the reactors and relevant material properties
%Defining the operating point that the model will be linearized at
if op_mode == 1
    ReactorSF.TcA= 275; %Cold-side coolant temperature (K)
    ReactorSF.TcB= 315; %Hot-side coolant temperature (K)
    TiA= 280; %Temperature in reactor A
    TiB= 310; %Temperature in reactor B
    wiA= 0.006; %Weight fraction in reactor A
    wiB= 0.007; %Weight fraction in reactor B
    PiA= 4.8e5; %Pressure in reactor A
    PiB= 2.8e5; %Pressure in reactor B
    mdotA= 0.2; %Cold-side coolant mass flow rate
    mdotB= 0.2; %Hot-side coolant mass flow rate
    dP  = 2.1e5; %Compressor pressure difference (Pa)
    x0= [TiA,PiA,wiA,TiB,PiB,wiB]'; %State variables at the linearization point
    u0= [mdotA,mdotB,dP,ReactorSF.TcA,ReactorSF.TcB]'; %Input and distrubance variables at the linearization point
else
    ReactorSF.TcA= 285; %Cold-side coolant temperature (K)
    ReactorSF.TcB= 303; %Hot-side coolant temperature (K)
    TiA= 280; %Temperature in reactor A
    TiB= 308; %Temperature in reactor B
    wiA= 0.006; %Weight fraction in reactor A
    wiB= 0.007; %Weight fraction in reactor B
    PiA= 2.9e5; %Pressure in reactor A
    PiB= 3.6e5; %Pressure in reactor B
    mdotA= 0.2; %Cold-side coolant mass flow rate
    mdotB= 0.2; %Hot-side coolant mass flow rate
    dP  = -8e4; %Compressor pressure difference (Pa)
    x0= [TiB,PiB,wiB,TiA,PiA,wiA]'; %State variables at the linearization point
    u0= [mdotB,mdotA,-dP,ReactorSF.TcB,ReactorSF.TcA]'; %Input and distrubance variables at the linearization point
end
%Creating the linear state-space matrices for the linear model using the
%linearization point defined above
[A,B,Bd,F_x0,C,D,Dd,G_x0] = linearize_model(x0,u0,op_mode,ReactorSF,Reactors,InProps);