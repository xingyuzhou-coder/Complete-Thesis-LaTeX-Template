clear all;
% main function

%parameters
global m0        % electron mass (Unit:kg)
global h         % Planck's constant (unit:J*s)
global h1        % Planck's constant (unit:J*s)
global Kb        % Boltzmann constant (Unit:eV/K)
global pi        % pi
global eps0      % vacuum dielectric constant (Unit:F/m)
global e         % electron charge (Unit:C)
global ev        % ev to J (unit:J/ev)
global T         % room temperature (Unit:K)
global N         % Mesh number
global L         % Well structure (unit:m)
global dz        % Step lenth (unit:m) 
global me        % Electron effecive mass distribution
global mhh       % hole effecive mass
global Ed        % Donor ionization energy distribution (Unit:eV)
global eps       % dielectric constant (Unit:F/m)
global Esf3      % Material-3 surface potential (Unit:eV)
global miu;      % Sheet polarizition charge(unit:C/m^2)
global A         % Schrodinger Matrix
global Psi       % Wavefunction
global Psi_0     % Origin Wavefunction
global Ei        % eigenergy (Unit:eV)
global Psi1      % 1st Wavefunction
global Ei1       % 1st eigenergy (Unit:eV)
global na        % Free electron sheet concentration (Unit: 1/m^2)
global n         % Free electron distribution (unit:1/m^3)
global n_j_sum   % Total Free electron concentration (Unit: 1/m^2)
global Nd        % ionized donor concentration (unit:1/m^3)
global Nd_sum    % Total ionized donor concentration (unit:1/m^2)
global Nd0_sum   % Total donor concentration (unit:1/m^2)
global Esd3      % Material-3 surface donor state ionization energy (Unit:eV)
global Esa3      % Material-3 acceptor donor state ionization energy (Unit:eV)
global Esf3      % Material-3 surface potential (Unit:eV)
global Nsd0      % surface donor state concentration (unit:1/m^2)
global Nsa0      % surface acceptor state concentration (unit:1/m^2)
global Nsd       % surface donor state concentration (unit:1/m^3)
global Nsa       % surface donor state concentration (unit:1/m^3)
global ro        % Total charge concentration (unit:1/m^2)
global sum_ro    % Square Total charge concentration (unit:1/m^3)
global Charge    % Total charge distribution (unit:C^2/m/ev)
global aB        % Effective Bohr radius
global rs        % the ratio of the mean electron distance to the effective Bohr radius
global C         % Poisson Matrix
global Vex       % local exchange potential
global V         % total potential
global Nd0
global na_j
global n_j

% load parameter
Parameter;                  
global Interlayer  % Interlayer width (unit:m)
global Rms         % surface RMS (unit:m)
global ii
for ii=3:3;                               
    Interlayer(ii) = (ii-1)*4.9816e-10;
    Structure;                         
    results = num2str((ii-1)/2);  
    HEFT = strcat('HEFT', results, '.txt'); % output file name
    
% original iteration value
global V_0       
global V_1        
global Fermi;
global Standard

% Conduction band
global Ec         % Conduction band
global Ec0        % Oringin conduction band
global Ei_c       % eigenergy (Unit:eV)
global Psi1_c     % 1st Wavefunction
global Ei1_c      % 1st eigenergy (Unit:eV)
global Psi_c      % Free electron distribution
global VH         % Hatree pontential
global Iteration  % Iteration time
global V_dif      % Force difference (eV)
 
V_0 = Ec0-Ec0;
V_1 = Ec0;
Iteration = 0; 
while max(abs(V_1-V_0)) > Standard   
    Iteration = Iteration + 1;
    V_1 = 0.01*V_1 + (1-0.01)*V_0;        %0.01 relaxation factor
    V_0 = V_1;
    V = V_0;
    m = me;
    Psi = Schrodinger(m,V);  
    Fermi = V(1) - Esf3; 
    Fermi_level(Fermi);
    VH = Poisson;
    V_1 = Ec0 - VH + Vex;  
    V_dif(Iteration) = max(abs(V_1-V_0));
    if Iteration >= 800;
        break;
    end
end
Ec = V_1 - 0.5*Fermi;               
Ei_c = Ei- Fermi;          
Psi_c = Psi;          
Psi1_c = Psi1;     

% electron concentration at each layer
global Interface1  % interface-1 mesh number    
global Interface2  % interface-2 mesh number    
global dz          % Step lenth (unit:m)        
global n           % Free electron concentration (Unit: 1/m^3)
global n_per       % Different layer electron concentration percent (Unit: 1/m^3)
global FWHM        % FWHM of electron  (Unit: m)

n_AlGaN(ii) = sum(n(1:Interface1-1-Rms/dz))/sum(n)*100; 
n_AlN(ii)   = sum(n(Interface1-Rms/dz:Interface2+Rms/dz))/sum(n)*100;
n_GaN(ii)   = sum(n(Interface2+1+Rms/dz:N+1))/sum(n)*100;

%FWHM of electron distribution
for i=1:N+1;    
  if n(i) >= mean(n);
        FWHM_1(i) = 1;
  else
        FWHM_1(i) = 0;
  end
end

FWHM(ii) = sum(FWHM_1)*dz*1e10; %unit :A

% Graph output
figure (1);
plot(1:N+1,Ec);
hold on;                            
title ('energy bend')           

figure (2);
plot(1:N+1,Psi1_c);
grid on;
title ('electron wavefunction');    

figure (3);                          
subplot (4,1,1);            
plot(1:N+1,VH);
grid on;                      
title('hartree force');           
subplot (4,1,2);           
plot(1:N+1,Nd);
grid on;                       
title('donor concentration');       
subplot (4,1,3);                   
plot(1:N+1,n);
grid on;
title('carrier concentration');     
subplot (4,1,4);                     
plot(1:N+1,ro);
grid on;                             
title('sheet charge concentration')  

figure(4)
plot(1:N+1,ro);
grid on;
hold on;
title('sheet charge concentration')

%txt output
fid=fopen(HEFT,'at');
fprintf(fid,'%s\t %f\t %s\t \n','ALN depth(A)', Interlayer(ii)*1e10);
fprintf(fid,'%s\t %g\t %s\t %f\t \n','2DEG sheet carrier density:',sum(n*dz), 'FWHM(A):', FWHM(ii));
fprintf(fid,'%s\t %f\t %s\t %f\t %s\t %f\t \n', 'AlGaN-electron:',n_AlGaN,  'AlN-electron:', n_AlN, 'GaN-electron:', n_GaN);
fprintf(fid,'%s\t %f\t  \n',  'Iteration:', Iteration);
fprintf(fid,'%s\t %f\t  \n',  'Fermi level:', Fermi);
fprintf(fid,'%s\t %s\t %s\t \n', 'depth(A):','Ec:','carrier concentration:');
iii = 1:N+1;
Out = [iii',Ec',n]';  
fprintf(fid,'%3.0f\t %4.5f\t %g\t \n', Out);
fclose(fid);
end