%linearize_model.m
%

%This function takes as an input the operating mode of a metal hydride
%system and a operating point to linearize around, and produces a
%linearized version of the governing equations, of the form f(x)= A(x-x0) +
%B(u-u0) + f(x0), by outputting the matrices A, B and f(x0)

%linearize_model.m requires:
%  Purdue Metal Hydride Toolbox
%  CalculateEffectiveness.m

%@author: Patrick Krane      email: pkrane@purdue.edu

function [A,B,Bd,F_x0,C,D,Dd,G_x0] = linearize_model(x0,u0,op_mode,ReactorSF,Reactors,InProps)
%addpath 'W:\MetalHydrideToolbox-master'
%% Define initial conditions and variables
%The order of variables in x and u varies depending on system state, so
%that each x_n or u_n never changes whether it is absorbing or desorbing.
%The equations work for a absorbing tank, m, and desorbing tank, n
%In charge mode    (op_mode= 1), m = A and n = B
%In discharge mode (op_mode= 0), m = B and n = A
%Commented-out section: variables are defined outside this function
if op_mode == 0
    mdotA= u0(2);
    mdotB= u0(1);
else
    mdotA= u0(1);
    mdotB= u0(2);
end

K_loss = Reactors.Kloss_H;
ReactorSF      = CalculateEffectiveness(mdotA,mdotB,op_mode,ReactorSF,Reactors,InProps);

%Defining hydride properties on a molar basis (values taken from Toolbox)
%Properties are for MmNi(4.5)Cr(0.5) and LaNi(5); change if materials
%change
Ea_Aa     = 30867;
Ea_Bd     = 16420;
Ca_Aa     = 620;
Ca_Bd     = 9.57;
Ea_Ad     = 28250;
Ea_Ba     = 20000;
Ca_Ad     = 90;
Ca_Ba     = 50;
A_A       = 1065.4;
A_B       = 322.08;
dH_Aa     = -23534;
dH_Ba     = -31168;
dH_Ad     = -25500;
dH_Bd     = -32151;
dS_Aa     = -96.675;
dS_Ba     = -111.4;
dS_Ad     = -100.3;
dS_Bd     = -112.8;

%% Linearization Equation Constants
%Defining values for constants in the linearization equations.
%As mentioned earlier, for op_mode= 1, m = A and n = B, so x0 = 
%[TA0,PA0,wA0,TB0,PB0,wB0] and u0 = [mA0,mB0,dP0,TcA,TcB]
if op_mode ~= 0
    K_1       = InProps.cp_c / Reactors.ntube;
    K_2       = - Ea_Aa / InProps.Rbar;
    K_3       = Ca_Aa * Reactors.hydA.dH('abs') * Reactors.mA * Reactors.hydA.wMax / Reactors.ntube;
    K_4       = Ca_Aa * Reactors.hydA.dH('abs') * Reactors.mA  / Reactors.ntube;
    K_5       = InProps.Patm;
    K_6       = InProps.Rbar;
    K_7       = (pi * Reactors.D_Hline ^ 2) / (4 * sqrt(K_loss) * Reactors.ntube);
    K_8       = InProps.rho_H_1b / 1e5;
    K_9       = Reactors.hydA.Cp * Reactors.mA / Reactors.ntube;
    K_10      = Reactors.eps * InProps.cp_H * Reactors.V_rA / (InProps.R_H * Reactors.ntube);
    K_11      = Reactors.hydA.Cp * Reactors.mA / Reactors.ntube;
    K_12      = -Ca_Aa * Reactors.mA * Reactors.hydA.wMax * InProps.R_H /  (Reactors.eps * Reactors.V_rA);
    K_13      = -Ca_Aa * Reactors.mA * InProps.R_H /  (Reactors.eps * Reactors.V_rA);
    K_14      = -(pi * Reactors.D_Hline ^ 2 * InProps.R_H) / (4 * sqrt(K_loss) * Reactors.eps * Reactors.V_rA);
    K_15      = Ca_Aa * Reactors.hydA.wMax;
    K_16      = Ca_Aa;
    K_17      = Reactors.hydB.Cp * Reactors.mB/ Reactors.ntube;
    K_18      = Reactors.eps * InProps.cp_H * Reactors.V_rB / (InProps.R_H * Reactors.ntube);
    %There is no K_19; it was used in a previous version of this code but
    %is no longer neededl
    K_20      = Ca_Bd * Reactors.hydB.dH('des') * Reactors.mB / Reactors.ntube;
    K_21      = - Ea_Bd / InProps.Rbar;
    K_22      = Reactors.hydB.Cp * Reactors.mB / Reactors.ntube;
    K_23      = -Ca_Bd * Reactors.mB * InProps.R_H /  (Reactors.eps * Reactors.V_rB);
    K_24      = -(pi * Reactors.D_Hline ^ 2 * InProps.R_H) / (4 * sqrt(K_loss) * Reactors.eps * Reactors.V_rB);
    K_25      = Ca_Bd;
    K_26      = InProps.R_H / InProps.cp_H;
    K_27      = InProps.cp_H;
    K_28      = 1 / Reactors.comp_eff;
    K_29      = InProps.R_H / Reactors.comp_eff;
    K_30      = Reactors.hydA.wMax;
    K_31      = (4 * InProps.Rbar) / (A_A + 4* InProps.Rbar * Reactors.hydA.Tc);
    K_33      = 2 * InProps.Rbar * Reactors.hydA.Tc;
    K_34      = dH_Aa;
    K_35      = -dS_Aa;
    K_36      = A_A;
    K_38      = Reactors.hydB.wMax;
    K_39      = (4 * InProps.Rbar) / (A_B + 4* InProps.Rbar * Reactors.hydB.Tc);
    K_41      = 2 * InProps.Rbar * Reactors.hydB.Tc;
    K_42      = dH_Bd;
    K_43      = -dS_Bd;
    K_44      = A_B;
    K_46      = -2 * pi * Reactors.Ri * Reactors.LA * Reactors.ntube / InProps.cp_c;
    K_47      = -2 * pi * Reactors.Ri * Reactors.LB * Reactors.ntube/ InProps.cp_c;
    x_alphA   = ((K_30 / 2) * (1 - sqrt(1 - K_31 * x0(1))));
    x_betaA   = ((K_30 / 2) * (1 + sqrt(1 - K_31 * x0(1))));
    x_alphB   = ((K_38 / 2) * (1 - sqrt(1 - K_39 * x0(4))));
    x_betaB   = ((K_38 / 2) * (1 + sqrt(1 - K_39 * x0(4))));
    mu_alph0A = K_34 + K_35 * x0(1) + K_36 * (x_alphA / K_30 - 0.5) - (K_33 * (1 - 2 * x_alphA / K_30) + K_6 * x0(1) * log(x_alphA / (K_30 - x_alphA)));
    mu_beta0A = K_34 + K_35 * x0(1) + K_36 * (x_betaA / K_30 - 0.5) - (K_33 * (1 - 2 * x_betaA / K_30) + K_6 * x0(1) * log(x_betaA / (K_30 - x_betaA)));
    mu_alph0B = K_42 + K_43 * x0(4) + K_44 * (x_alphB / K_38 - 0.5) - (K_41 * (1 - 2 * x_alphB / K_30) + K_6 * x0(4) * log(x_alphB / (K_38 - x_alphB)));
    mu_beta0B = K_42 + K_43 * x0(4) + K_44 * (x_betaB / K_38 - 0.5) - (K_41 * (1 - 2 * x_betaB / K_30) + K_6 * x0(4) * log(x_betaB / (K_38 - x_betaB)));
    K_32      = mu_alph0A;
    K_37      = mu_beta0A;
    K_40      = mu_alph0B;
    K_45      = mu_beta0B;
    eff_m     = ReactorSF.eff_A;
    qhm       = ReactorSF.h_coolA / mdotA ^ 0.8;
    if (ReactorSF.TcA > 240) && (ReactorSF.TcA < 273)
        Pr_m      = 1.55632892 * ReactorSF.TcA^2 - 843.81295994 * ReactorSF.TcA + 114544.44198525;
        dPrm_du4  = 2 * 1.55632892 * ReactorSF.TcA - 843.81295994;
    else
        Pr_m      = (7.99696268e-6) * (ReactorSF.TcA^4) - (0.0108427886 * ReactorSF.TcA^3) + (5.508419166 * ReactorSF.TcA^2) - (1243.0847 * ReactorSF.TcA) + 105189.5;
        dPrm_du4  = 4 * (7.99696268e-6) * (ReactorSF.TcA^3) - (3 * 0.0108427886 * ReactorSF.TcA^2) + (2 * 5.508419166 * ReactorSF.TcA) - 1243.0847;
    end
    if (ReactorSF.TcB > 240) && (ReactorSF.TcB < 273)
        Pr_n      = 1.55632892 * ReactorSF.TcB^2 - 843.81295994 * ReactorSF.TcB + 114544.44198525;
        dPrn_du5  = 2 * 1.55632892 * ReactorSF.TcB - 843.81295994;
    else
        Pr_n      = (7.99696268e-6) * (ReactorSF.TcB^4) - (0.0108427886 * ReactorSF.TcB^3) + (5.508419166 * ReactorSF.TcB^2) - (1243.0847 * ReactorSF.TcB) + 105189.5;
        dPrn_du5  = 4 * (7.99696268e-6) * (ReactorSF.TcB^3) - (3 * 0.0108427886 * ReactorSF.TcB^2) + (2 * 5.508419166 * ReactorSF.TcB) - 1243.0847;
    end
    dqhmdu4   = 0.023 * (2 * InProps.cp_c / (pi * Reactors.ntube)) ^ 0.8 * (InProps.kcool ^ 0.2 / ( 2 * Reactors.Ri ^ 1.8)) * (-0.4 * Pr_m ^ -1.4) * dPrm_du4;
    eff_n     = ReactorSF.eff_B;
    qhn       = ReactorSF.h_coolB / mdotB ^ 0.8;
    dqhndu4   = 0.023 * (2 * InProps.cp_c / (pi * Reactors.ntube)) ^ 0.8 * (InProps.kcool ^ 0.2 / ( 2 * Reactors.Ri ^ 1.8)) * (-0.5 * Pr_n ^ -1.5) * dPrn_du5;
%For op_mode= 0, m = B and n = A, so x0 = 
%[TB0,PB0,wB0,TA0,PA0,wA0] and u0 = [mB0,mA0,dP0,TcB,TcA]
else
    K_1       = InProps.cp_c / Reactors.ntube;
    K_2       = - Ea_Ba / InProps.Rbar;
    K_3       = Ca_Ba * Reactors.hydB.dH('abs') * Reactors.mB * Reactors.hydB.wMax / Reactors.ntube;
    K_4       = Ca_Ba * Reactors.hydB.dH('abs') * Reactors.mB  / Reactors.ntube;
    K_5       = InProps.Patm;
    K_6       = InProps.Rbar;
    K_7       = (pi * Reactors.D_Hline ^ 2) / (4 * sqrt(K_loss) * Reactors.ntube);
    K_8       = InProps.rho_H_1b / 1e5;
    K_9       = Reactors.hydB.Cp * Reactors.mB / Reactors.ntube;
    K_10      = Reactors.eps * InProps.cp_H * Reactors.V_rB / (InProps.R_H * Reactors.ntube);
    K_11      = Reactors.hydB.Cp * Reactors.mB / Reactors.ntube;
    K_12      = -Ca_Ba * Reactors.mB * Reactors.hydB.wMax * InProps.R_H /  (Reactors.eps * Reactors.V_rB);
    K_13      = -Ca_Ba * Reactors.mB * InProps.R_H /  (Reactors.eps * Reactors.V_rB);
    K_14      = -(pi * Reactors.D_Hline ^ 2 * InProps.R_H) / (4 * sqrt(K_loss) * Reactors.eps * Reactors.V_rB);
    K_15      = Ca_Ba * Reactors.hydB.wMax;
    K_16      = Ca_Ba;
    K_17      = Reactors.hydA.Cp * Reactors.mA / Reactors.ntube;
    K_18      = Reactors.eps * InProps.cp_H * Reactors.V_rA / (InProps.R_H * Reactors.ntube);
    K_20      = Ca_Ad * Reactors.hydA.dH('des') *Reactors. mA / Reactors.ntube;
    K_21      = - Ea_Ad / InProps.Rbar;
    K_22      = Reactors.hydA.Cp * Reactors.mA / Reactors.ntube;
    K_23      = -Ca_Ad * Reactors.mA * InProps.R_H /  (Reactors.eps * Reactors.V_rA);
    K_24      = -(pi * Reactors.D_Hline ^ 2 * InProps.R_H) / (4 * sqrt(K_loss) * Reactors.eps * Reactors.V_rA);
    K_25      = Ca_Ad;
    K_26      = InProps.R_H / InProps.cp_H;
    K_27      = InProps.cp_H;
    K_28      = 1 / Reactors.comp_eff;
    K_29      = InProps.R_H / Reactors.comp_eff;
    K_30      = Reactors.hydB.wMax;
    K_31      = (4 * InProps.Rbar) / (A_B + 4* InProps.Rbar * Reactors.hydB.Tc);
    K_33      = 2 * InProps.Rbar * Reactors.hydB.Tc;
    K_34      = dH_Ba;
    K_35      = -dS_Ba;
    K_36      = A_B;
    K_38      = Reactors.hydA.wMax;
    K_39      = (4 * InProps.Rbar) / (A_A + 4* InProps.Rbar * Reactors.hydA.Tc);
    K_41      = 2 * InProps.Rbar * Reactors.hydA.Tc;
    K_42      = dH_Ad;
    K_43      = -dS_Ad;
    K_44      = A_A;
    K_46      = -2 * pi * Reactors.Ri * Reactors.LB * Reactors.ntube / InProps.cp_c;
    K_47      = -2 * pi * Reactors.Ri * Reactors.LA * Reactors.ntube / InProps.cp_c;
    x_alphB   = ((K_30 / 2) * (1 - sqrt(1 - K_31 * x0(1))));
    x_betaB   = ((K_30 / 2) * (1 + sqrt(1 - K_31 * x0(1))));
    x_alphA   = ((K_38 / 2) * (1 - sqrt(1 - K_39 * x0(4))));
    x_betaA   = ((K_38 / 2) * (1 + sqrt(1 - K_39 * x0(4))));
    mu_alph0B = K_34 + K_35 * x0(1) + K_36 * (x_alphB / K_30 - 0.5) - (K_33 * (1 - 2 * x_alphB / K_30) + K_6 * x0(1) * log(x_alphB / (K_38 - x_alphB)));
    mu_beta0B = K_34 + K_35 * x0(1) + K_36 * (x_betaB / K_30 - 0.5) - (K_33 * (1 - 2 * x_betaB / K_30) + K_6 * x0(1) * log(x_betaB / (K_38 - x_betaB)));
    mu_alph0A = K_42 + K_43 * x0(4) + K_44 * (x_alphA / K_38 - 0.5) - (K_41 * (1 - 2 * x_alphA / K_30) + K_6 * x0(4) * log(x_alphA / (K_30 - x_alphA)));
    mu_beta0A = K_42 + K_43 * x0(4) + K_44 * (x_betaA / K_38 - 0.5) - (K_41 * (1 - 2 * x_betaA / K_30) + K_6 * x0(4) * log(x_betaA / (K_30 - x_betaA)));
    K_32      = mu_alph0B;
    K_37      = mu_beta0B;
    K_40      = mu_alph0A;
    K_45      = mu_beta0A;
    eff_m     = ReactorSF.eff_B;
    if (ReactorSF.TcA > 240) && (ReactorSF.TcA < 273)
        Pr_n      = 1.55632892 * ReactorSF.TcA^2 - 843.81295994 * ReactorSF.TcA + 114544.44198525;
        dPrn_du5  = 2 * 1.55632892 * ReactorSF.TcA - 843.81295994;
    else
        Pr_n      = (7.99696268e-6) * (ReactorSF.TcA^4) - (0.0108427886 * ReactorSF.TcA^3) + (5.508419166 * ReactorSF.TcA^2) - (1243.0847 * ReactorSF.TcA) + 105189.5;
        dPrn_du5  = 4 * (7.99696268e-6) * (ReactorSF.TcA^3) - (3 * 0.0108427886 * ReactorSF.TcA^2) + (2 * 5.508419166 * ReactorSF.TcA) - 1243.0847;
    end
    if (ReactorSF.TcB > 240) && (ReactorSF.TcB < 273)
        Pr_m      = 1.55632892 * ReactorSF.TcB^2 - 843.81295994 * ReactorSF.TcB + 114544.44198525;
        dPrm_du4  = 2 * 1.55632892 * ReactorSF.TcB - 843.81295994;
    else
        Pr_m      = (7.99696268e-6) * (ReactorSF.TcB^4) - (0.0108427886 * ReactorSF.TcB^3) + (5.508419166 * ReactorSF.TcB^2) - (1243.0847 * ReactorSF.TcB) + 105189.5;
        dPrm_du4  = 4 * (7.99696268e-6) * (ReactorSF.TcB^3) - (3 * 0.0108427886 * ReactorSF.TcB^2) + (2 * 5.508419166 * ReactorSF.TcB) - 1243.0847;
    end
    qhm       = ReactorSF.h_coolB / mdotB ^ 0.8;
    dqhmdu4   = 0.023 * (2 * InProps.cp_c / (pi * Reactors.ntube)) ^ 0.8 * (InProps.kcool ^ 0.2 / Reactors.Ri ^ 1.8) * (-0.4 * Pr_m ^ -1.4) * dPrm_du4;
    eff_n     = ReactorSF.eff_A;
    qhn       = ReactorSF.h_coolA / mdotA ^ 0.8;
    dqhndu4   = 0.023 * (2 * InProps.cp_c / (pi * Reactors.ntube)) ^ 0.8 * (InProps.kcool ^ 0.2 / Reactors.Ri ^ 1.8) * (-0.5 * Pr_n ^ -1.5) * dPrn_du5;
end
%Derivatives of the heat exchanger effectiveness: this depends on mdot and
%Tc, so for heat exchanger m, it depends on u1 and u4; for heat exchanger
%n, it depends on u2 and u5
deffmdu1    = 0.2 * K_46 * qhm * exp(K_46 * qhm / u0(1) ^ 0.2) / (u0(1) ^ 1.2);
deffmdu4    = -K_46 * dqhmdu4 * exp(K_46 * qhm / u0(1) ^ 0.2) / (u0(1) ^ 0.2);
deffndu2    = 0.2 * K_47 * qhn * exp(K_47 * qhn / u0(2) ^ 0.2) / (u0(2) ^ 1.2);
deffndu5    = -K_47 * dqhndu4 * exp(K_47 * qhn / u0(2) ^ 0.2) / (u0(2) ^ 0.2);
%Chemical potential, mu, for hydrides m and n
if x0(3) < ((K_30 / 2) * (1 - sqrt(1 - K_31 * x0(1))))
    mum     = K_32 + K_33 * (1 - 2 * x0(3) / K_30) + K_6 * x0(1) * log(x0(3) / (K_30 - x0(3)));
    dmumdx1 = K_6 * log(x0(3) / (K_30 - x0(3)));
    dmumdx3 = (2 * K_33 / K_30) + (K_6 * x0(1) / x0(3)) * (1 + x0(3) / (K_30 - x0(3))); 
elseif x0(3) > ((K_30 / 2) * (1 + sqrt(1 - K_31 * x0(1))))
    mum     = K_37 + K_33 * (1 - 2 * x0(3) / K_30) + K_6 * x0(1) * log(x0(3) / (K_30 - x0(3)));
    dmumdx1 = K_6 * log(x0(3) / (K_30 - x0(3)));
    dmumdx3 = (2 * K_33 / K_30) + (K_6 * x0(1) / x0(3)) * (1 + x0(3) / (K_30 - x0(3)));
else
    mum     = K_34 + K_35 * x0(1) + K_36 * (x0(3) / K_30 - 0.5);
    dmumdx1 = K_35;
    dmumdx3 = K_36 / K_30;
end
if x0(6) < ((K_38 / 2) * (1 - sqrt(1 - K_39 * x0(4))))
    mun     = K_40 + K_41 * (1 - 2 * x0(6) / K_38) + K_6 * x0(4) * log(x0(6) / (K_38 - x0(6)));
    dmundx4 = K_6 * log(x0(6) / (K_38 - x0(6)));
    dmundx6 = (2 * K_41 / K_38) + (K_6 * x0(4) / x0(6)) * (1 + x0(6) / (K_38 - x0(6)));
elseif x0(6) > ((K_38 / 2) * (1 + sqrt(1 - K_39 * x0(4))))
    mun     = K_45 + K_41 * (1 - 2 * x0(6) / K_38) + K_6 * x0(4) * log(x0(6) / (K_38 - x0(6)));
    dmundx4 = K_6 * log(x0(6) / (K_38 - x0(6)));
    dmundx6 =  (2 * K_41 / K_38) + (K_6 * x0(4) / x0(6)) * (1 + x0(6) / (K_38 - x0(6)));
else
    mun     = K_42 + K_43 * x0(4) + K_44 * (x0(6) / K_38 - 0.5);
    dmundx4 = K_43;
    dmundx6 = K_44 / K_38;
end
%Inlet hydrogen enthalpy for reactor m, and its derivatives with respect to
%x1, x4, x5, and u3; this equals zero when hydrogen flows from the absorbing
%to the desorbing reactor. In that case, we instead have inlet enthalpy for
%reactor n, and its derivatives with respect to x1 and x4 (since flow is
%usually from the desorbing to the absorbing reactor, it is assumed flow
%only goes the other way when the compressor is not active)
if sign(x0(5)+(u0(3))-x0(2)) == 1
    hm          = K_27 * (x0(4)* (1 + K_28 *((((x0(5) + (u0(3)))/x0(5))^ K_26)-1))-x0(1));
    dhmdx1      = -K_27;
    dhmdx4      = K_27 * (1 + K_28 * (((x0(5) + (u0(3)))/x0(5))^ K_26 - 1));
    dhmdx5      = -K_29 * x0(4) * ((u0(3)) * (x0(5) + (u0(3))) ^ (K_26 - 1)) / (x0(5) ^ (K_26+1));
    dhmdu3      = K_29 * x0(4) * (x0(5) + (u0(3))) ^ (K_26 - 1) / (x0(5) ^ K_26);
    hn          = 0;
    dhndx1      = 0;
    dhndx4      = 0;
else
    hm          = 0;
    dhmdx1      = 0;
    dhmdx4      = 0;
    dhmdx5      = 0;
    dhmdu3      = 0;
    hn          = K_27 * (x0(1) - x0(4));
    dhndx1      = K_27;
    dhndx4      = -K_27;
end

%% Defining A, B, and F(x0)
%f(x)= [dT_m/dt, dP_m/dt, dw_m/dt,dT_n/dt, dP_n/dt, dw_n/dt]
%A= df/dx and B= df/du
% Define the state matrix, A, of the linearized model
A(1,1)      = (-K_1*u0(1)*eff_m+exp(K_2/x0(1))*(K_3-K_4*x0(3))*((-K_2/(x0(1)^2))*(log(x0(2)/K_5)-mum/(K_6*x0(1)))+mum/(K_6*x0(1)^2)-dmumdx1/(K_6*x0(1)))...
            +K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*dhmdx1)/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1))...
            +((K_10*x0(2)/x0(1)^2)*(K_1*eff_m*u0(1)*(u0(4)-x0(1))+exp(K_2/x0(1))*(K_3-K_4*x0(3))*(log(x0(2)/K_5)-mum/(K_6*x0(1)))...
            +K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*hm))/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1))^2; %d(f_1)/d(T_m)
A(1,2)      = ((exp(K_2/x0(1))/x0(2))*(K_3-K_4*x0(3))+K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8/((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2))))...
            *((u0(3))-2*x0(2))*(hm/2))/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1))-(K_10/x0(1))*(K_1*eff_m*u0(1)*(u0(4)-x0(1))+exp(K_2/x0(1))*...
            (K_3-K_4*x0(3))*(log(x0(2)/K_5)-mum/(K_6*x0(1)))+K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*hm)/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1))^2; %d(f_1)/d(P_m)
A(1,3)      = (exp(K_2/x0(1))*(-K_4*(log(x0(2)/K_5)-mum/(K_6*x0(1)))-((K_3-K_4*x0(3))*dmumdx3/(K_6*x0(1)))))/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1))...
            -K_11*(K_1*eff_m*u0(1)*(u0(4)-x0(1))+exp(K_2/x0(1))*(K_3-K_4*x0(3))*(log(x0(2)/K_5)-mum/(K_6*x0(1)))...
            +K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*hm)/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1))^2; %d(f_1)/d(w_m)
A(1,4)      = K_7 * sign(x0(5)+(u0(3))-x0(2))* sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*dhmdx4/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1)); %d(f_1)/d(T_n)
A(1,5)      = ((K_7*sqrt(K_8))/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1)))*(((((u0(3))+2*x0(5))*hm)/(2*sqrt((x0(5)+x0(2))...
            *abs(x0(5)+(u0(3))-x0(2)))))+sign(x0(5)+(u0(3))-x0(2))*sqrt((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*dhmdx5); %d(f_1)/d(P_n)
A(1,6)      = 0; %d(f_1)/d(w_n)
A(2,1)      = exp(K_2/x0(1))*(K_12-K_13*x0(3))*((1-K_2/x0(1))*(log(x0(2)/K_5)-mum/(K_6*x0(1)))+mum/(K_6*x0(1))-dmumdx1/K_6)-K_14*sign(x0(5)+(u0(3))-x0(2))...
            *sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2))); %d(f_2)/d(T_m)
A(2,2)      = x0(1)*((K_12-K_13*x0(3))*exp(K_2/x0(1))/x0(2)-K_14*((-2*x0(2)+(u0(3)))/2)*sqrt(K_8/((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2))))); %d(f_2)/d(P_m)
A(2,3)      = x0(1)*exp(K_2/x0(1))*(-K_13*(log(x0(2)/K_5)-mum/(K_6*x0(1)))-(K_12-K_13*x0(3))*(dmumdx3/(K_6*x0(1)))); %d(f_2)/d(w_m)
A(2,4)      = 0; %d(f_2)/d(T_n)
A(2,5)      = -(K_14/2)*x0(1)*(2*x0(5)+(u0(3)))*sqrt(K_8/((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))); %d(f_2)/d(P_n)
A(2,6)      = 0; %d(f_2)/d(w_n)
A(3,1)      = (K_15-K_16*x0(3))*exp(K_2/x0(1))*((-K_2/(x0(1)^2))*(log(x0(2)/K_5)-mum/(K_6*x0(1)))+mum/(K_6*x0(1)^2)-dmumdx1/(K_6*x0(1))); %d(f_3)/d(T_m)
A(3,2)      = (K_15-K_16*x0(3)) * exp(K_2/x0(1)) / x0(2); %d(f_3)/d(P_m)
A(3,3)      = exp(K_2/x0(1))*(-K_16*(log(x0(2)/K_5)-mum/(K_6*x0(1)))-(K_15-K_16*x0(3))*(dmumdx3/(K_6*x0(1)))); %d(f_3)/d(w_m)
A(3,4)      = 0; %d(f_3)/d(T_n)
A(3,5)      = 0; %d(f_3)/d(P_n)
A(3,6)      = 0; %d(f_3)/d(w_n)
A(4,1)      = -K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*dhndx1/(K_17+K_18*x0(5)/x0(4)+K_22*x0(6)); %d(f_4)/d(T_m)
A(4,2)      = (K_7*sqrt(K_8)/(K_17+K_18*x0(5)/x0(4)+K_22*x0(6)))*((2*x0(5)-(u0(3)))*hn)/(2*sqrt((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))); %d(f_4)/d(P_m)     %Change if using hn
A(4,3)      = 0; %d(f_4)/d(w_m)
A(4,4)      = (-K_1*u0(2)*eff_n+K_20*x0(6)*exp(K_21/x0(4))*((-K_21/(x0(4)^2))*(log(x0(5)/K_5)-mun/(K_6*x0(4)))+mun/(K_6*x0(4)^2)-dmundx4/(K_6*x0(4)))-K_7*sign(x0(5)+(u0(3))-x0(2))...
            *sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*dhndx4)/((K_17+K_18*x0(5)/x0(4)+K_22*x0(6)))+(K_18*x0(5)/(x0(4)^2))*(K_1*eff_n*u0(2)*(u0(5)-x0(4))+K_20*x0(6)*exp(K_21/x0(4))...
            *(log(x0(5)/K_5)-mun/(K_6*x0(4)))-K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*hn)/(K_17+K_18*x0(5)/x0(4)+K_22*x0(6))^2; %d(f_4)/d(T_n)
A(4,5)      = ((K_20*x0(6)*exp(K_21/x0(4))/x0(5))-K_7*sqrt(K_8)*sign(x0(5)+(u0(3))-x0(2))*((u0(3))+2*x0(5))*hn/(2*sqrt((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))))/...
            ((K_17+K_18*x0(5)/x0(4)+K_22*x0(6)))-((K_18/x0(4))*(K_1*eff_n*u0(2)*(u0(5)-x0(4))+K_20*x0(6)*exp(K_21/x0(4))*(log(x0(5)/K_5)-mun/(K_6*x0(4))))-K_7*sign(x0(5)+(u0(3))-x0(2))*...
            sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*hn)/(K_17+K_18*x0(5)/x0(4)+K_22*x0(6))^2; %d(f_4)/d(P_n)
A(4,6)      = (K_20*exp(K_21/x0(4))*(log(x0(5)/K_5)-mun/(K_6*x0(4))-x0(6)*dmundx6/(K_6*x0(4))))/(K_17+K_18*x0(5)/x0(4)+K_22*x0(6))...
            -K_22*(K_1*eff_n*u0(2)*(u0(5)-x0(4))+K_20*x0(6)*exp(K_21/x0(4))*(log(x0(5)/K_5)-mun/(K_6*x0(4)))-K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*hn)...
            /(K_17+K_18*x0(5)/x0(4)+K_22*x0(6))^2; %d(f_4)/d(w_n)
A(5,1)      = 0; %d(f_5)/d(T_m)
A(5,2)      = (K_24*sqrt(K_8)*x0(4)*(-2*x0(2)+(u0(3))))/(2*sqrt(((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2))))); %d(f_5)/d(P_m)
A(5,3)      = 0; %d(f_5)/d(w_m)
A(5,4)      = exp(K_21/x0(4))*(K_23*x0(6)*((1-K_21/x0(4))*(log(x0(5)/K_5)-mun/(K_6*x0(4)))+mun/(K_6*x0(4))-dmundx4/K_6))+K_24*sign(x0(5)+(u0(3))-x0(2))*...
            (sqrt(K_8*((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2))))); %d(f_5)/d(T_n)
A(5,5)      = x0(4)*((K_23*x0(6)*exp(K_21/x0(4))/x0(5))+(K_24*(2*x0(5)+(u0(3)))/2)*sqrt(K_8/((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2))))); %d(f_5)/d(P_n)
A(5,6)      = K_23*x0(4)*exp(K_21/x0(4))*(log(x0(5)/K_5)-mun/(K_6*x0(4))+(-dmundx6*x0(6))/(K_6*x0(4))); %d(f_5)/d(w_n)
A(6,1)      = 0; %d(f_6)/d(T_m)
A(6,2)      = 0; %d(f_6)/d(P_m)
A(6,3)      = 0; %d(f_6)/d(w_m)
A(6,4)      = K_25*x0(6)*exp(K_21/x0(4))*((-K_21/(x0(4)^2))*(log(x0(5)/K_5)-mun/(K_6*x0(4)))+mun/(K_6*x0(4)^2)-dmundx4/(K_6*x0(4))); %d(f_6)/d(T_n)
A(6,5)      = K_25 * x0(6) * exp(K_21/x0(4)) / x0(5); %d(f_6)/d(P_n)
A(6,6)      = K_25*exp(K_21/x0(4))*(log(x0(5)/K_5)-mun/(K_6*x0(4))-x0(6)*dmundx6/(K_6*x0(4))); %d(f_6)/d(w_n)

% Define the input matrix, B, of the linearized model
B(1,1)      = K_1 * (u0(4) - x0(1)) * (u0(1) * deffmdu1 + eff_m) / (K_9+K_11*x0(3)+K_10*x0(2)/x0(1)); %d(f_1)/d(mdot_m)
B(1,2)      = 0; %d(f_1)/d(mdot_n)
B(1,3)      = (K_7*sqrt(K_8)/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1)))*(((x0(5)+x0(2))*hm)/(2*sqrt(((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))))+dhmdu3*sign(x0(5)+(u0(3))-x0(2))*sqrt(((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2))))); %d(f_1)/d(delP_comp)
B(1,4)      = K_1 * u0(1) * ((u0(4)-x0(1))*deffmdu4 + eff_m) / (K_9+K_11*x0(3)+K_10*x0(2)/x0(1)); %d(f_1)/d(TK_m)
B(1,5)      = 0; %d(f_1)/d(TK_n)
B(2,1)      = 0; %d(f_2)/d(mdot_m)
B(2,2)      = 0; %d(f_2)/d(mdot_n)
B(2,3)      = -0.5 * K_14 * x0(1) * (x0(5) + x0(2)) * sqrt(K_8/((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))); %d(f_2)/d(delP_comp)
%-0.5 * K_14 * x0(1) * (x0(5) + x0(2)) * sqrt(K_8/((x0(5)+x0(2))*abs(x0(5)+u0(3)-x0(2))));
B(2,4)      = 0; %d(f_2)/d(TK_m)
B(2,5)      = 0; %d(f_2)/d(TK_n)
B(3,1)      = 0; %d(f_3)/d(mdot_m)
B(3,2)      = 0; %d(f_3)/d(mdot_n)
B(3,3)      = 0; %d(f_3)/d(delP_comp)
B(3,4)      = 0; %d(f_3)/d(TK_m)
B(3,5)      = 0; %d(f_3)/d(TK_n)
B(4,1)      = 0; %d(f_4)/d(mdot_m)
B(4,2)      = K_1 * (u0(5) - x0(4)) * (u0(2) * deffndu2 + eff_n) / (K_17+K_18*x0(5)/x0(4)+K_22*x0(6)); %d(f_4)/d(mdot_n)
B(4,3)      = ((-K_7 * (x0(5) + x0(2)) * hn) / (2 * (K_17 + K_18 * (x0(5) / x0(4)) + K_22 * x0(6)))) * sqrt(K_8 / ((x0(5) + x0(2)) * abs(x0(5) + u0(3) - x0(2)))); %d(f_4)/d(delP_comp)
B(4,4)      = 0; %d(f_4)/d(TK_m)
B(4,5)      = K_1 * u0(2) * ((u0(5)-x0(4))*deffndu5 + eff_n) / (K_17+K_18*x0(5)/x0(4)+K_22*x0(6)); %d(f_4)/d(TK_n)
B(5,1)      = 0; %d(f_5)/d(mdot_m)
B(5,2)      = 0; %d(f_5)/d(mdot_n)
B(5,3)      = 0.5 * K_24 * x0(4) * (x0(5) + x0(2)) * sqrt(K_8/((x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))); %d(f_5)/d(delP_comp)
B(5,4)      = 0; %d(f_5)/d(TK_m)
B(5,5)      = 0; %d(f_5)/d(TK_n)
B(6,1)      = 0; %d(f_6)/d(mdot_m)
B(6,2)      = 0; %d(f_6)/d(mdot_n)
B(6,3)      = 0; %d(f_6)/d(delP_comp)
B(6,4)      = 0; %d(f_6)/d(TK_m)
B(6,5)      = 0; %d(f_6)/d(TK_n)

Bmat = B; clear B;

B = Bmat(:,1:3);
Bd = Bmat(:,4:5);

% Evaluate the system derivatives at the nominal operating condition
F_x0(1)     = (K_1*eff_m*u0(1)*(u0(4)-x0(1))+exp(K_2/x0(1))*(K_3-K_4*x0(3))*(log(x0(2)/K_5)-mum/(K_6*x0(1)))...
            +K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*hm)/(K_9+K_11*x0(3)+K_10*x0(2)/x0(1));
F_x0(2)     = x0(1)*(exp(K_2/x0(1))*(K_12-K_13*x0(3))*(log(x0(2)/K_5)-mum/(K_6*x0(1)))-K_14*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2))));
F_x0(3)     = (K_15-K_16*x0(3))*exp(K_2/x0(1))*(log(x0(2)/K_5)-mum/(K_6*x0(1)));
F_x0(4)     = (K_1*eff_n*u0(2)*(u0(5)-x0(4))+K_20*x0(6)*exp(K_21/x0(4))*(log(x0(5)/K_5)-mun/(K_6*x0(4)))...
            -K_7*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2)))*hn)/(K_17+K_18*x0(5)/x0(4)+K_22*x0(6));
F_x0(5)     = x0(4)*(K_23 * x0(6) * exp(K_21/x0(4)) * (log(x0(5)/K_5)-mun/(K_6*x0(4)))+K_24*sign(x0(5)+(u0(3))-x0(2))*sqrt(K_8*(x0(5)+x0(2))*abs(x0(5)+(u0(3))-x0(2))));
F_x0(6)     = K_25*x0(6)*exp(K_21/x0(4))*(log(x0(5)/K_5)-mun/(K_6*x0(4)));
F_x0        = F_x0';

%% %% Defining G_x0, C and D
G_x0(1)= InProps.cp_c * eff_m * u0(1) * (x0(1) - u0(4));
G_x0(2)= InProps.cp_c * eff_n * u0(2) * (x0(4) - u0(5));

C= zeros(2,6);
C(1,1)= InProps.cp_c * eff_m * u0(1);
C(2,4)= InProps.cp_c * eff_n * u0(2);

D= zeros(2,5);
D(1,1)= InProps.cp_c * (x0(1) - u0(4)) * (eff_m + deffmdu1 * (u0(1)));
D(2,2)= InProps.cp_c * (x0(4) - u0(5)) * (eff_n + deffndu2 * (u0(2)));
D(1,4)= InProps.cp_c * u0(1) * (-eff_m + deffmdu4 * (x0(1) - u0(4)));
D(2,5)= InProps.cp_c * u0(2) * (-eff_n + deffndu5 * (x0(4) - u0(5)));

Dmat = D; clear D;

D = Dmat(:,1:3);
Dd = Dmat(:,4:5);