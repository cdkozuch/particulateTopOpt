function Z = partcomp2d(tvf,res,plotflag)
%% PARTCOMP generates 2D microstructure for particle composite
% Inputs:
%   tvf (scalar): target particle volume fraction
%   res (scalar): resolution of material map
%   plotflag (bool): true for plotting
% Outputs:
%   Z (res x res): map of material (1=matrix, 2=particle) 

%% define radii for particles
% Current assumptions for particle size are ~6.5 microns on average
% with a normal distribution between 4 and 9 microns.
% Waiting on confirmation from Sigma

davg = 6.5e-6; % average diameter
dmin = 4e-6; % diameter lower bound
dmax = 9e-6; % diameter upper bound (not needed since assumed symmetric)
sigma_count = 2; % number of standard deviations between davg and dmin
sigma = (davg - dmin)/sigma_count; % standard deviation
np = 1e3; % number of particles
Rtot = 0.5*normrnd(davg,sigma,np,1); % list of possible radii

%% Generate centers for particles
lx = 3.5e1*davg; % width of RVE
ly = lx; % height of RVE
ng = 1e3; % number of grid points

cxv = linspace(-lx/2,lx/2,ng);
cyv = linspace(-ly/2,ly/2,ng);
iv = randperm(ng,np);
jv = randperm(ng,np);
Ctot = [cxv(iv)',cyv(jv)']; % list of possible centers

%% Choose radii and centers to satisfy volume fraction & no intersection
% this method is extremely inefficient

fnp = 0;
pvf = 0;

for i=1:np
    % check for collision
    if fnp >= 2
        s = Ctot(i,:);
        S = repmat(s,fnp,1);
        if any(sum((S-C).^2,2) <= (Rtot(i)+R).^2)
            continue
        end

        % check for wall collisions
        if abs(s(1))+Rtot(i) > lx/2 || abs(s(2))+Rtot(i) > ly/2
            continue
        end
    end
    
    % add particle
    fnp = fnp + 1;
    R(fnp,1) = Rtot(i);
    C(fnp,:) = Ctot(i,:);
    
    pvf = pvf + pi*(R(fnp)^2)/(lx*ly);
    
    % stop when target volume fraction is reached
    if pvf > tvf
        break
    end
end

%% Generate microstructure
xgv = linspace(-lx/2,lx/2,res); % gridpoints in x direction
ygv = linspace(-lx/2,ly/2,res); % gridpoints in y direction
[X,Y] = meshgrid(xgv,ygv);
Z = ones(size(X)); % 1 for matrix, 2 for particle
for i=1:res
    for j=1:res
        p = [X(i,j),Y(i,j)];
        P = repmat(p,fnp,1);
        if any(sum((P-C).^2,2) <= R.^2);
            Z(i,j) = 2;
        end
    end
end

if plotflag
    histogram(R*2)
    title('Particle Diameter Distribution')
    figure
    imagesc(Z); axis off
    title('Particle Placement')
end
end
