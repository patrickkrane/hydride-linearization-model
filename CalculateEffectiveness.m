%CalculateEffectiveness.m
%

%This function calculates heat exchanger effectiveness for the metal
%hydride reactors for input coolant temperatures and mass flow rates

%CalculateEffectiveness.m requires:=
% N/A

%@author: Patrick Krane      email: pkrane@purdue.edu
function ReactorSF = CalculateEffectiveness(mdotA,mdotB,op_mode,ReactorSF,Reactors,InProps)
ReactorSF.mdotA_tube= mdotA / Reactors.ntube; %Mass flow through one tube in reactor A
ReactorSF.mdotB_tube= mdotB / Reactors.ntube; %Mass flow through one tube in reactor B
if (ReactorSF.mdotA_tube > 0)
    %Calculating Prandtl number as a function of temperature
    if (ReactorSF.TcA > 240) && (ReactorSF.TcA < 273)
        Pr_A  = 1.55632892*ReactorSF.TcA^2 - 843.81295994*ReactorSF.TcA + 114544.44198525;
    else
        Pr_A  = (7.99696268e-6) * (ReactorSF.TcA^4) - (0.0108427886 * ReactorSF.TcA^3) + (5.508419166 * ReactorSF.TcA^2) - (1243.0847 * ReactorSF.TcA) + 105189.5;
    end
    mu_A      = InProps.kcool * Pr_A / InProps.cp_c; %Kinematic viscosity- used to calculate Reynolds number
    %Determining effectiveness of the reactors from the coolant properties
    ReA       = 2 * ReactorSF.mdotA_tube / (pi * Reactors.Ri * mu_A); %Coolant flow Reynolds number
    if op_mode == 1
        Nu_A  = 0.023 * ReA ^ 0.8 * Pr_A ^ 0.4; %Turbulent Nusselt number
    else
        Nu_A  = 0.023 * ReA ^ 0.8 * Pr_A ^ 0.3; %Turbulent Nusselt number
    end
    ReactorSF.h_coolA   = Nu_A * InProps.kcool / (2 * Reactors.Ri); %Heat transfer coefficient for reactor A (W/ m^2 K)
    NTU_A     = (2 * pi * Reactors.Ri * Reactors.LA * ReactorSF.h_coolA) / (ReactorSF.mdotA_tube * InProps.cp_c); %Number of transfer units
    ReactorSF.eff_A     = 1 - exp(-NTU_A); %Effectiveness for reactor A
else
    %If no flow in the reactor
    ReactorSF.h_coolA = 0;
    ReactorSF.eff_A   = 0;
end
if (ReactorSF.mdotB_tube > 0)
    %Calculating Prandtl number as a function of temperature
    if ReactorSF.TcB < 273
        Pr_B  = 1.55632892*ReactorSF.TcB^2 - 843.81295994*ReactorSF.TcB + 114544.44198525;
    else
        Pr_B  = (7.99696268e-6) * (ReactorSF.TcB^4) - (0.0108427886 * ReactorSF.TcB^3) + (5.508419166 * ReactorSF.TcB^2) - (1243.0847 * ReactorSF.TcB) + 105189.5;
    end
    mu_B      = InProps.kcool * Pr_B / InProps.cp_c; %Kinematic viscosity- used to calculate Reynolds number
    ReB       = 2* ReactorSF.mdotB_tube / (pi * Reactors.Ri * mu_B); %Coolant flow Reynolds number
    if op_mode == 1
        Nu_B  = 0.023 * ReB ^ 0.8 * Pr_B ^ 0.3; %Turbulent Nusselt number
    else
        Nu_B  = 0.023 * ReB ^ 0.8 * Pr_B ^ 0.4; %Turbulent Nusselt number
    end
    ReactorSF.h_coolB   = Nu_B * InProps.kcool / (2 * Reactors.Ri); %Heat transfer coefficient for reactor B (W/ m^2 K)- seperate for if/when seperate Nusselt #s
    NTU_B     = (2 * pi * Reactors.Ri * Reactors.LB * ReactorSF.h_coolB) / (ReactorSF.mdotB_tube * InProps.cp_c); %Number of transfer units
    ReactorSF.eff_B     = 1 - exp(-NTU_B); %Effectiveness for reactor B
else
    %If no flow in the reactor
    ReactorSF.h_coolB = 0;
    ReactorSF.eff_B = 0;
end