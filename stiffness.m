%This code is based on code written by Siva Srinivas Kolukula
%https://www.mathworks.com/matlabcentral/fileexchange/31788-the-plane-stress-problem?focused=5205477&tab=function

function k = stiffness(D)
coordinates = 0.5*[1,1;-1,1;-1,-1;1,-1]; %nodal coordinates
xx = coordinates(:,1); yy = coordinates(:,2);
nnel=4; %number of nodes per element
ndof=2; %number of dofs per node (UX,UY)
edof=nnel*ndof; %degrees of freedom per element     
%--------------------------------------------------------------------------
%  initialization of matrices and vectors
%--------------------------------------------------------------------------
B = zeros(3,edof); %kinematic matrix for bending
dhdx = zeros(1,2);
dhdy = zeros(1,2);
jacobian=zeros(2,2);
%--------------------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%--------------------------------------------------------------------------
k = zeros(edof,edof); %initialization of stiffness matrix
%--------------------------------------------------------------------------
%  numerical integration for stiffness matrix
%--------------------------------------------------------------------------
G = [-0.577350269189626; 0.577350269189626]; %gausspoints
for intx=1:2
    xi = G(intx,1); %sampling point in x-axis
    for inty=1:2
        eta = G(inty,1); %sampling point in y-axis         
        %shape functions
        dhdr = 0.25*[eta-1, 1-eta, 1+eta, -eta-1];
        dhds = 0.25*[xi-1, -xi-1, xi+1, 1-xi];
        %Jacobian
        for i=1:nnel
            jacobian(1,1)=jacobian(1,1)+dhdr(i)*xx(i);
            jacobian(1,2)=jacobian(1,2)+dhdr(i)*yy(i);
            jacobian(2,1)=jacobian(2,1)+dhds(i)*xx(i);
            jacobian(2,2)=jacobian(2,2)+dhds(i)*yy(i);
        end        
        detjacob=det(jacobian); %determinant of Jacobian
        invjacob=inv(jacobian); %inverse of Jacobian matrix        
        %shape function derivatives
        for i=1:nnel
            dhdx(i)=invjacob(1,1)*dhdr(i)+invjacob(1,2)*dhds(i);
            dhdy(i)=invjacob(2,1)*dhdr(i)+invjacob(2,2)*dhds(i);
        end
        %kinematic stiffness
        for i=1:nnel
            i1=(i-1)*2+1;  
            i2=i1+1;
            B(1,i1)=dhdx(i);
            B(2,i2)=dhdy(i);
            B(3,i1)=dhdy(i);
            B(3,i2)=dhdx(i);
        end
        k = k+B'*D*B*detjacob; %stiffness matrix
    end
end
end
