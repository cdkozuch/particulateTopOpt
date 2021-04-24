clear

%% define flags
drawflag = 0;

%% Inclusion inputs
% creates a binary field for homogenization
% places an ellipse inclusion in a square

lx = 1; % width of RVE
ly = 1; % height of RVE
phi = 90; % angle between x and y dimension of RVE

num_theta = 5; % number of inclusion angles
m = 100; % m x m grid
vf = 0.15; % inclusion area fraction

num_ar = 25; % number of aspect ratios
max_a = 0.95*lx/2;
ar = linspace(1,(max_a^2)*pi/vf,num_ar); % aspect ratio fractions
a = sqrt(vf*ar/pi); % size of major radius
b = a./ar; % size of minor radius
theta = linspace(0,pi,num_theta); % angle between major axis and x axis

%% Homogenization inputs
E1 = 1e9*[10,0.1]'; % elastic modulus for matrix
num_stiff = numel(E1); %number of stiffness ratios
E2 = 1e9*ones(num_stiff,1); % elastic modulus for inclusion
nu1 = 0.33; % poisson's ratio for matrix
nu2 = 0.33; % poisson's ratio for inclusion

% Lame parameters
lambda = [E1*nu1/((1+nu1)*(1-2*nu1)), E2*nu2/((1+nu2)*(1-2*nu2))];
mu = [E1/(2*(1+nu1)),E2/(2*(1+nu2))];
lambda = 2*mu.*lambda./(lambda + 2*mu); %adjust for plane stress

%% Topology Optimization inputs
nelx = 200; % number of elements in x
nely = 100; % number of elements in y
volfrac = 0.5; % target volume fraction
penal = 3; % penalty parameter
rmin = 1.5; % filter radius
ft = 1; % filter type (1 = sensativity filtering, 2 = density filtering)

%% Loop

% preallocate
emat = zeros(num_stiff,num_ar,num_theta);
rhomat = zeros(nely,nelx,num_stiff,num_ar,num_theta);

for prob=1:2 %type of 2D problem (1 = Michell Beam, 2 = half MBB Beam)
    for i=2:num_stiff
        for j=num_ar:num_ar
            x = inclusion(a(j),b(j),0,m,drawflag);
            C = homogenize(lx, ly, lambda(i,:), mu(i,:), phi, x);

            parfor k=1:num_theta
                r = cos(theta(k)); s = sin(theta(k));
                T1 = [r^2 s^2 2*r*s; s^2 r^2 -2*r*s; -r*s r*s (r^2)-(s^2)];
                T2 = [r^2 s^2 r*s; s^2 r^2 -r*s; -2*r*s 2*r*s (r^2)-(s^2)];
                Crot = (T1\C)*T2;

                [emat(i,j,k),rhomat(:,:,i,j,k)] = anisoTopOpt(prob,...
                    nelx,nely,Crot,volfrac,penal,rmin,ft,drawflag);
            end
        end
    end
    
    switch prob
        case 1
            filename = sprintf('data_michell4.mat');
            save(filename)
        case 2
            filename = sprintf('data_mbb4.mat');
            save(filename)
    end
end





