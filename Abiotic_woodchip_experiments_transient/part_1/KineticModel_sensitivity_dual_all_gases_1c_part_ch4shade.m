clear

%first-type BC;
%kinetic model;
%see Van Genuchten pg. 16 and 28:



%data = csvread('ObjC_1_gas.csv',1,0);
data = csvread('Obj1c_part1_gas.csv',1,0);
data_time = data(:,1).*60; %time in sec
ethane = data(:,2);
helium = data(:,3);
n2o = data(:,4);
sf6 = data(:,5);
ch4 = data(:,6);



data_br = csvread('Obj1c_part1_br.csv',1,0);
%data_br = csvread('Obj1b_br.csv',1,0);
Br_time = data_br(:,1).*60
Br = data_br(:,2);


%If Vg/Vw = .01, then
%R_He = 2
%R_SF6 = 2.5
%R_ethane = 1.2
%R_N2O = 1.015

%time units are in seconds;
%length units are in m;

C_o = 1; %dimensionless
t = 0:50:50000;   %default delta t is 25
D = 2.1875e-6; 


%2.69e-6, dispersion coefficients, m2/s, default 5e-7; determined from bromide only.
    %%%%%%3.1e-6 for bromide
v_m = 7.875e-5;  

%velocity, 7.09e-5,  m/s, default 1e-4 m/s
    %%%%%%%3e-5 for bromide
    %try 1.19e-4 to 1.43e-4, expected values based on reactor geometry +
    %porosity
L = .3;
theta_m = 0.586; %mobile phase volume (dynamic, macropores)
theta_im = 0.31; %immobile phase volume (stagnant, micropores)
phi_m = theta_m/(theta_m + theta_im)  %fraction of liquid phase that is mobile

P = v_m*L/D; %Peclet number
PV = v_m.*t.*phi_m./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
data_PV_br = v_m.*Br_time.*phi_m./L;

data_PV = v_m.*data_time.*phi_m./L

data_PV_br_plot = v_m.*Br_time.*phi_m./L


figure(11)
clf
h1 = plot(data_PV_br_plot, Br, '*k','MarkerSize', 18)
hold on
h2 = plot(data_PV, n2o, 'dm','MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'MarkerSize', 18)
hold on
h3 = plot(data_PV, sf6, 'vg','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'MarkerSize', 18)
hold on
h4 = plot(data_PV, ethane, 'sb','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'MarkerSize', 18)
hold on
h5 = plot(data_PV, ch4, 'oc','MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', 'MarkerSize', 18)
hold on
h6 = plot(data_PV, helium, '^r','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 18)
hold on





measured = Br; %define which solute you are modeling

R = 1.000001; %ethane 2.7
%R = [2.0:.5:6.0];
R_m = 1.000001;
%R_m = [1.0:.5:3.0];
alpha = 0.0001;  %ethane 1.3e-4
%alpha = [.9e-4:.1e-4:2.5e-4];
loopcnt = 0;

results_matrix = zeros(length(v_m)*length(D),4);
%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is alpha
%column 4 of results_matrix is sum of squared residuals



beta = phi_m*R_m/R;


for j = 1:length(v_m)
    
    for k = 1:length(D)
        
        
w = alpha*L*(R-1)/(v_m(j)*theta_m);
P = v_m(j)*L/D(k); %Peclet number
PV = v_m(j).*t.*phi_m./L;     %pore volumes, tau in equations in Table 2 in Van Genuchten
data_PV_br = v_m(j).*Br_time.*phi_m./L;


%G term in Van Genuchten, equation SI-1 in Table 2:
G_T = .5 *  erfc( sqrt(P./(4.*beta.*R.*PV)).*(beta*R-PV)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R.*PV) ).*(beta*R + PV) );
%IMPORTANT NOTE: Equation SI-1 in Table 2 has the final term as beta*R-tau.
%I think this is an error and that it should be beta*R+tau.  At the top of page 18 they note: 
%"Expressions for G(T) in Table 2 follow from those given in Table 1 by replacing T by
%tau and R by (beta*R).  In Table 1, the final term is R+T, so based on
%their statement the final term should be B*R+tau.
%This is how it is coded above, and this equation gives an effluent C
%ranging from 0 to 1, while the original equation from Van Genuchten has an
%effluent C ranging from 0 to a number much greater than 1.

%exponential term multiplied by G
term = exp(-w.*PV/(beta*R));


integrand_vec = zeros(1,length(PV));

delta_tau = 0.001;   %tau time increment, also d_tau dummy variable
                    %default .001

for i = 2:length(PV);  %2:length(PV)
    
    tau = [0:delta_tau:PV(i)];
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R.*tau)).*(beta*R-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R.*tau) ).*(beta*R + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* ( I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


prod = G_tau.*H;

integrand = delta_tau*trapz(prod(1:length(prod)-1));   %last element of product is NaN, since the last element of H is NaN ...
                                                    %due to division by 0
                                                    %in b term
integrand_vec(i) = integrand;


end


BrC = G_T.*term + w/R .* integrand_vec;


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV_br, edges);

RSS = sum((BrC(I) - measured').^2)


loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = v_m(j);
     results_matrix(loopcnt,3) = D(k);
     results_matrix(loopcnt,4) = RSS;

    end
    
end

plot(PV, BrC,'-k','LineWidth',1.5, 'handlevisibility','off')
hold on
%End of Bromide










measured = n2o; %define which solute you are modeling

%fit parameters:
%R = 2;      %MUST BE > 1
%alpha = 1e-3;    %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha
%for determination of R-squared:



%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter




R = 1.35; %ethane 2.7
%R = [1.1:0.02:1.2];
%R_m = [1.1:0.1:2.8];
R_m = 1.2;
alpha = 9.4e-5;  %ethane 1.3e-4
%alpha = [1e-5 1e-6 1e-7 1e-8 1e-9];%8e-6:4e-6:8e-5];

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),5);
%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is R_m
%column 4 of results_matrix is alpha
%column 5 of results_matrix is sum of squared residuals




for j = 1:length(R)
    
    for k = 1:length(alpha)
        
        for l = 1:length(R_m)
        
        
w = alpha(k)*L*(R(j)-1)/(v_m*theta_m);
beta = phi_m*R_m(l)/R(j);



%G term in Van Genuchten, equation SI-1 in Table 2:
G_T = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*PV)).*(beta*R(j)-PV)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*PV) ).*(beta*R(j) + PV) );
%IMPORTANT NOTE: Equation SI-1 in Table 2 has the final term as beta*R-tau.
%I think this is an error and that it should be beta*R+tau.  At the top of page 18 they note: 
%"Expressions for G(T) in Table 2 follow from those given in Table 1 by replacing T by
%tau and R by (beta*R).  In Table 1, the final term is R+T, so based on
%their statement the final term should be B*R+tau.
%This is how it is coded above, and this equation gives an effluent C
%ranging from 0 to 1, while the original equation from Van Genuchten has an
%effluent C ranging from 0 to a number much greater than 1.

%exponential term multiplied by G
term = exp(-w.*PV/(beta*R(j)));

%some figures to visualize how each of these terms are affected by
%different fit parameters:
%use these to troubleshoot problems if the model acts weird with certain
%paremeters.
%figure(5)
%plot(PV,term)
%title('Term')

integrand_vec = zeros(1,length(PV));

delta_tau = 0.001;   %tau time increment, also d_tau dummy variable
                    %default .001

for i = 2:length(PV);  %2:length(PV)
    
    tau = [0:delta_tau:PV(i)];
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* ( I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


prod = G_tau.*H;

integrand = delta_tau*trapz(prod(1:length(prod)-1));   %last element of product is NaN, since the last element of H is NaN ...
                                                    %due to division by 0
                                                    %in b term
integrand_vec(i) = integrand;


end


%first_term = G_T.*term;
%figure(14)
%plot(PV,first_term)
%title('G*exp term')

%C = G(1:length(G)-1).*term(1:length(term)-1) + (w/R).*ConvInt;   %equation 35 in van genuchten
C = G_T.*term + w/R(j) .* integrand_vec;


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV, edges);

RSS = sum((C(I) - measured').^2)



loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = R(j);
     results_matrix(loopcnt,3) = R_m(l);
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;

    end
    
    end

end

C(C > 1) = 1;

plot(PV, C,'-m','LineWidth',1.5, 'handlevisibility','off')
hold on



%%%End of N2O

measured = sf6; %define which solute you are modeling

%fit parameters:
%R = 2;      %MUST BE > 1
%alpha = 1e-3;    %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha
%for determination of R-squared:



%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter




R = 48.75; %ethane 2.7
%R = [1.1:0.02:1.2];
%R_m = [1.1:0.1:2.8];
R_m = 1.175;
alpha = 1.35e-6;  %ethane 1.3e-4
%alpha = [1e-5 1e-6 1e-7 1e-8 1e-9];%8e-6:4e-6:8e-5];

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),5);
%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is R_m
%column 4 of results_matrix is alpha
%column 5 of results_matrix is sum of squared residuals




for j = 1:length(R)
    
    for k = 1:length(alpha)
        
        for l = 1:length(R_m)
        
        
w = alpha(k)*L*(R(j)-1)/(v_m*theta_m);
beta = phi_m*R_m(l)/R(j);



%G term in Van Genuchten, equation SI-1 in Table 2:
G_T = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*PV)).*(beta*R(j)-PV)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*PV) ).*(beta*R(j) + PV) );
%IMPORTANT NOTE: Equation SI-1 in Table 2 has the final term as beta*R-tau.
%I think this is an error and that it should be beta*R+tau.  At the top of page 18 they note: 
%"Expressions for G(T) in Table 2 follow from those given in Table 1 by replacing T by
%tau and R by (beta*R).  In Table 1, the final term is R+T, so based on
%their statement the final term should be B*R+tau.
%This is how it is coded above, and this equation gives an effluent C
%ranging from 0 to 1, while the original equation from Van Genuchten has an
%effluent C ranging from 0 to a number much greater than 1.

%exponential term multiplied by G
term = exp(-w.*PV/(beta*R(j)));

%some figures to visualize how each of these terms are affected by
%different fit parameters:
%use these to troubleshoot problems if the model acts weird with certain
%paremeters.
%figure(5)
%plot(PV,term)
%title('Term')

integrand_vec = zeros(1,length(PV));

delta_tau = 0.001;   %tau time increment, also d_tau dummy variable
                    %default .001

for i = 2:length(PV);  %2:length(PV)
    
    tau = [0:delta_tau:PV(i)];
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* ( I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


prod = G_tau.*H;

integrand = delta_tau*trapz(prod(1:length(prod)-1));   %last element of product is NaN, since the last element of H is NaN ...
                                                    %due to division by 0
                                                    %in b term
integrand_vec(i) = integrand;


end


%first_term = G_T.*term;
%figure(14)
%plot(PV,first_term)
%title('G*exp term')

%C = G(1:length(G)-1).*term(1:length(term)-1) + (w/R).*ConvInt;   %equation 35 in van genuchten
C = G_T.*term + w/R(j) .* integrand_vec;


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV, edges);

RSS = sum((C(I) - measured').^2)


loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = R(j);
     results_matrix(loopcnt,3) = R_m(l);
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;

    end
    
    end

end


plot(PV, C,'-g','LineWidth',1.5, 'handlevisibility','off')
hold on



%%%End of SF6








measured = ethane; %define which solute you are modeling

%fit parameters:
%R = 2;      %MUST BE > 1
%alpha = 1e-3;    %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha
%for determination of R-squared:



%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter




R = 6.1; %ethane 2.7
%R = [1.1:0.02:1.2];
%R_m = [1.1:0.1:2.8];
R_m = 1.225;
alpha = 1.68e-5;  %ethane 1.3e-4
%alpha = [1e-5 1e-6 1e-7 1e-8 1e-9];%8e-6:4e-6:8e-5];

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),5);
%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is R_m
%column 4 of results_matrix is alpha
%column 5 of results_matrix is sum of squared residuals




for j = 1:length(R)
    
    for k = 1:length(alpha)
        
        for l = 1:length(R_m)
        
        
w = alpha(k)*L*(R(j)-1)/(v_m*theta_m);
beta = phi_m*R_m(l)/R(j);



%G term in Van Genuchten, equation SI-1 in Table 2:
G_T = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*PV)).*(beta*R(j)-PV)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*PV) ).*(beta*R(j) + PV) );
%IMPORTANT NOTE: Equation SI-1 in Table 2 has the final term as beta*R-tau.
%I think this is an error and that it should be beta*R+tau.  At the top of page 18 they note: 
%"Expressions for G(T) in Table 2 follow from those given in Table 1 by replacing T by
%tau and R by (beta*R).  In Table 1, the final term is R+T, so based on
%their statement the final term should be B*R+tau.
%This is how it is coded above, and this equation gives an effluent C
%ranging from 0 to 1, while the original equation from Van Genuchten has an
%effluent C ranging from 0 to a number much greater than 1.

%exponential term multiplied by G
term = exp(-w.*PV/(beta*R(j)));

%some figures to visualize how each of these terms are affected by
%different fit parameters:
%use these to troubleshoot problems if the model acts weird with certain
%paremeters.
%figure(5)
%plot(PV,term)
%title('Term')

integrand_vec = zeros(1,length(PV));

delta_tau = 0.001;   %tau time increment, also d_tau dummy variable
                    %default .001

for i = 2:length(PV);  %2:length(PV)
    
    tau = [0:delta_tau:PV(i)];
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* ( I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


prod = G_tau.*H;

integrand = delta_tau*trapz(prod(1:length(prod)-1));   %last element of product is NaN, since the last element of H is NaN ...
                                                    %due to division by 0
                                                    %in b term
integrand_vec(i) = integrand;


end


%first_term = G_T.*term;
%figure(14)
%plot(PV,first_term)
%title('G*exp term')

%C = G(1:length(G)-1).*term(1:length(term)-1) + (w/R).*ConvInt;   %equation 35 in van genuchten
C = G_T.*term + w/R(j) .* integrand_vec;


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV, edges);

RSS = sum((C(I) - measured').^2)


loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = R(j);
     results_matrix(loopcnt,3) = R_m(l);
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;

    end
    
    end

end


plot(PV, C,'-b','LineWidth',1.5, 'handlevisibility','off')
hold on



%%%End of Ethane





measured = ch4; %define which solute you are modeling

%fit parameters:
%R = 2;      %MUST BE > 1
%alpha = 1e-3;    %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha
%for determination of R-squared:



%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter




R = 7; %ethane 2.7
%R = [1.1:0.02:1.2];
%R_m = [1.1:0.1:2.8];
R_m = 1.185;
alpha = 1.6e-5;  %ethane 1.3e-4
%alpha = [1e-5 1e-6 1e-7 1e-8 1e-9];%8e-6:4e-6:8e-5];

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),5);
%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is R_m
%column 4 of results_matrix is alpha
%column 5 of results_matrix is sum of squared residuals




for j = 1:length(R)
    
    for k = 1:length(alpha)
        
        for l = 1:length(R_m)
        
        
w = alpha(k)*L*(R(j)-1)/(v_m*theta_m);
beta = phi_m*R_m(l)/R(j);



%G term in Van Genuchten, equation SI-1 in Table 2:
G_T = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*PV)).*(beta*R(j)-PV)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*PV) ).*(beta*R(j) + PV) );
%IMPORTANT NOTE: Equation SI-1 in Table 2 has the final term as beta*R-tau.
%I think this is an error and that it should be beta*R+tau.  At the top of page 18 they note: 
%"Expressions for G(T) in Table 2 follow from those given in Table 1 by replacing T by
%tau and R by (beta*R).  In Table 1, the final term is R+T, so based on
%their statement the final term should be B*R+tau.
%This is how it is coded above, and this equation gives an effluent C
%ranging from 0 to 1, while the original equation from Van Genuchten has an
%effluent C ranging from 0 to a number much greater than 1.

%exponential term multiplied by G
term = exp(-w.*PV/(beta*R(j)));

%some figures to visualize how each of these terms are affected by
%different fit parameters:
%use these to troubleshoot problems if the model acts weird with certain
%paremeters.
%figure(5)
%plot(PV,term)
%title('Term')

integrand_vec = zeros(1,length(PV));

delta_tau = 0.001;   %tau time increment, also d_tau dummy variable
                    %default .001

for i = 2:length(PV);  %2:length(PV)
    
    tau = [0:delta_tau:PV(i)];
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* ( I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


prod = G_tau.*H;

integrand = delta_tau*trapz(prod(1:length(prod)-1));   %last element of product is NaN, since the last element of H is NaN ...
                                                    %due to division by 0
                                                    %in b term
integrand_vec(i) = integrand;


end


%first_term = G_T.*term;
%figure(14)
%plot(PV,first_term)
%title('G*exp term')

%C = G(1:length(G)-1).*term(1:length(term)-1) + (w/R).*ConvInt;   %equation 35 in van genuchten
ch4C = G_T.*term + w/R(j) .* integrand_vec;


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV, edges);

RSS = sum((ch4C(I) - measured').^2)


loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = R(j);
     results_matrix(loopcnt,3) = R_m(l);
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;

    end
    
    end

end


plot(PV, ch4C,'-c','LineWidth',1.5, 'handlevisibility','off')
hold on



%%%End of CH4






measured = helium; %define which solute you are modeling


%fit parameters:
%R = 2;      %MUST BE > 1
%alpha = 1e-3;    %rate constant, try starting in the range 5e-4 to 5e-3

                    %1/s, 100 d^-1 from Geistlinger (2005) = 1e-3 1/s
                    %vulava has 5e-4 - 1e-3  1/s for alpha for SF6, so
                    %probably faster for Helium, N2O
                    %velocity in bioreactor on the order of .003 m/s, much slower than in Vulava, so 
                    %expect lower alpha
%for determination of R-squared:



%dimensionless parameters
%defined in Van Genuchten pg. 15, with the exception of w, beta, and R
%For w and beta, defined on pg. 28 for one-site kinetic
%non-equilibrium adsorption.  For R, comes from Fry et al. (1994) and Vulava et al. (2002): 
%R = 1 + KH*(Vg/Vw)   %use this with specific K_H values to determine Vg/Vw
%based on fit R parameter




R = 21; %ethane 2.7
%R = [1.1:0.02:1.2];
%R_m = [1.1:0.1:2.8];
R_m = 1.01;
alpha = 2.4e-5;  %ethane 1.3e-4
%alpha = [1e-5 1e-6 1e-7 1e-8 1e-9];%8e-6:4e-6:8e-5];

loopcnt = 0;

results_matrix = zeros(length(R)*length(alpha),5);
%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is R_m
%column 4 of results_matrix is alpha
%column 5 of results_matrix is sum of squared residuals




for j = 1:length(R)
    
    for k = 1:length(alpha)
        
        for l = 1:length(R_m)
        
        
w = alpha(k)*L*(R(j)-1)/(v_m*theta_m);
beta = phi_m*R_m(l)/R(j);



%G term in Van Genuchten, equation SI-1 in Table 2:
G_T = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*PV)).*(beta*R(j)-PV)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*PV) ).*(beta*R(j) + PV) );
%IMPORTANT NOTE: Equation SI-1 in Table 2 has the final term as beta*R-tau.
%I think this is an error and that it should be beta*R+tau.  At the top of page 18 they note: 
%"Expressions for G(T) in Table 2 follow from those given in Table 1 by replacing T by
%tau and R by (beta*R).  In Table 1, the final term is R+T, so based on
%their statement the final term should be B*R+tau.
%This is how it is coded above, and this equation gives an effluent C
%ranging from 0 to 1, while the original equation from Van Genuchten has an
%effluent C ranging from 0 to a number much greater than 1.

%exponential term multiplied by G
term = exp(-w.*PV/(beta*R(j)));

%some figures to visualize how each of these terms are affected by
%different fit parameters:
%use these to troubleshoot problems if the model acts weird with certain
%paremeters.
%figure(5)
%plot(PV,term)
%title('Term')

integrand_vec = zeros(1,length(PV));

delta_tau = 0.001;   %tau time increment, also d_tau dummy variable
                    %default .001

for i = 2:length(PV);  %2:length(PV)
    
    tau = [0:delta_tau:PV(i)];
    
G_tau = .5 *  erfc( sqrt(P./(4.*beta.*R(j).*tau)).*(beta*R(j)-tau)) + .5*exp(P).*erfc( sqrt (P./(4.*beta.*R(j).*tau) ).*(beta*R(j) + tau) );
    
    a = w.*tau;
    b = w.*(PV(i) - tau)./( (1 - beta)*R(j) );
    
    squig = 2*sqrt(a.*b);
    
%In equation 36, the I_0 and I_1 are modified bessel functions, see
%notation on pg. 50-51
I_0 = besseli(0,squig);  %modified bessel function of the first kind, order 0
I_1 = besseli(1,squig);  %modified bessel function of the first kind, order 1
    
%H term defined by Equation 36:
H = exp(-a-b).* ( I_0/beta + (squig.*I_1)./(2.*b.*(1-beta)) );


prod = G_tau.*H;

integrand = delta_tau*trapz(prod(1:length(prod)-1));   %last element of product is NaN, since the last element of H is NaN ...
                                                    %due to division by 0
                                                    %in b term
integrand_vec(i) = integrand;


end


%first_term = G_T.*term;
%figure(14)
%plot(PV,first_term)
%title('G*exp term')

%C = G(1:length(G)-1).*term(1:length(term)-1) + (w/R).*ConvInt;   %equation 35 in van genuchten
C = G_T.*term + w/R(j) .* integrand_vec;


%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV, edges);

RSS = sum((C(I) - measured').^2)


loopcnt = loopcnt + 1

     results_matrix(loopcnt,1) = loopcnt;
     results_matrix(loopcnt,2) = R(j);
     results_matrix(loopcnt,3) = R_m(l);
     results_matrix(loopcnt,4) = alpha(k);
     results_matrix(loopcnt,5) = RSS;

    end
    
    end

end

plot(PV, C,'-r','LineWidth',1.5, 'handlevisibility','off')
hold on
set(gca,'FontSize',20)
ylabel('C/C_o')
xlabel('PV')

%%%End of helium





% figure(111)
% clf
% %plot(PV, C,'-')
% ylim([-0.1 1.1])
% xlim([0 2])
% ylabel('C/C_o')
% xlabel('PV')



%column 1 of results_matrix is loopcnt
%column 2 of results_matrix is R
%column 3 of results_matrix is R_m
%column 4 of results_matrix is alpha
%column 5 of results_matrix is sum of squared residuals

y = [results_matrix(:,1)'; results_matrix(:,2)'; results_matrix(:,3)'; results_matrix(:,4)'; results_matrix(:,5)']
fileID = fopen('results_matrix.txt','w');
%fprintf(fileID, 'Exponential Function\n\n');
fprintf(fileID,'%0.3f %0.3e %0.3e %0.3e %0.4f\n',y);
fclose(fileID);




y = [PV; C];
fileID = fopen('exptable.txt','w');
%fprintf(fileID, 'Exponential Function\n\n');
fprintf(fileID,'%f %f\n',y);
fclose(fileID);



%find indices in PV that are closest PV to data: 
%edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV, edges);

%edges = [-Inf, mean([Br_time(2:end); Br_time(1:end-1)]), +Inf];
%I = discretize(Br_time, edges);



%find indices in PV that are closest PV to data: 
edges = [-Inf, mean([PV(2:end); PV(1:end-1)]), +Inf];
%I = discretize(data_PV_br, edges);
I = discretize(data_PV, edges);

RSS = sum((C(I) - measured').^2)
TSS = sum((measured - mean(measured)).^2)
R2 = 1 - (RSS/TSS)



figure(11)
%plot(PV, C,'--r')
ylim([-0.1 1.2])
xlim([0 8])
ylabel('C/C_o')
xlabel('PV')
%columnlegendpadded(2, {'Br', 'N_2O', 'SF_6', 'C_2H_6', 'CH_4', 'He'}, 'location', 'southeast')
% text(.4, .8, sprintf('R = %0.2f', R))
% text(.4, .9, sprintf('alpha = %0.2e s^{-1}', alpha))
% text(.4, .7, sprintf('R^2 = %0.3f', R2))
saveas(gcf,'output.eps','epsc')
%legend('Bromide','N_2O', 'SF_6', 'Ethane', 'Methane', 'Helium', 'Location', 'southeast'); %
%title('Abiotic, Transient, Woodchip')
set(gca, 'FontSize', 34);

transparency = 0.5
color = [0.5 0.5 0.5];
ch4fill=jbfill(PV,BrC,ch4C,'c', 'c', 0,transparency)
uistack(ch4fill,'bottom');
legend('Entrapped CH_4','Bromide','N_2O', 'SF_6', 'C_2H_6', 'CH_4', 'He',  'Location', 'southeast'); %








