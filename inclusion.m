function V = inclusion(a,b,theta,n,drawflag)
%INCLUSION creates a 2D binary field for input into homogenize function
%Input
%   a (scalar): length of major axis
%   b (scalar): length of minor axis
%   theta (scalar): angle between major axis and x axis
%   n (integer): n x n field
%   drawflag (boolean): 0 = draw figure, 1 = don't draw figure

%% create n x n node grid with size of 1 x 1
xgv = linspace(-0.5,0.5,n);
ygv = xgv;
[X,Y] = meshgrid(xgv,ygv);
V = zeros(size(X));

%% assign 1 for outside ellipse, 2 for inside ellipse
theta = -theta; %counter-clockwise
for i=1:n
    for j=1:n
        x = X(i,j); y = Y(i,j);
        val = ((x*cos(theta) + y*sin(theta))/a)^2 + ((x*sin(theta) - y*cos(theta))/b)^2;
        if val > 1
            V(i,j) = 1;
        else
            V(i,j) = 2;
        end
    end
end

%% draw figure
if drawflag
    figure
    colormap(gray)
    imagesc(V)
end
end
