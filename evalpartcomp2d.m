%% Mechanical Properties of DGEBA/DETA Epoxy
% Properties taken from Garcia et.al reference:
% Mech. Prop. of Epoxy Networks Based on DGEBA and Aliphatic Amines

E1=2.55e9; % Young's Modulus in GPa
nu1=0.35; % Poisson's Ratio 

%% Mechanical Properties of Carbonyl Iron Particles
% Properties found at: 
% https://www.americanelements.com/carbonyl-iron-powder-7439-89-6

E2=211e9; % Young's Modulus in GPa
nu2=0.29; % Poisson's Ratio 

%% Lame Inputs
lambda = [E1*nu1/((1+nu1)*(1-2*nu1)), E2*nu2/((1+nu2)*(1-2*nu2))];
mu = [E1/(2*(1+nu1)),E2/(2*(1+nu2))];
lambda = 2*mu.*lambda./(lambda + 2*mu); %adjust for plane stress

%% Generate microstructure
tvf = 0.15; % target particle volume fraction
res = 300; % microstructure resolution
plotflag = true;
fprintf('Creating microstructure...')
Z = partcomp2d(tvf,res,plotflag); % material map
fprintf('done.\n')

%% Compute stiffness
lx = 1; % width of RVE
ly = 1; % height of RVE
phi = 90; % angle of RVE corner (just leave this at 90)
fprintf('Computing stiffness...')
CH = homogenize(lx, ly, lambda, mu, phi, Z); % stiffness tensor
fprintf('done.\n')
disp(CH)