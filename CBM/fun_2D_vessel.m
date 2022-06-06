function vesselness = fun_2D_vessel(I, sigmas, spacing, tau, brightondark)


verbose = 1;

if nargin<5
    brightondark = false; % default mode for 2D is dark vessels compared to the background
end

I = single(I);

for j = 1:length(sigmas)
    
%     if verbose
%         disp(['Current filter scale (sigma): ' num2str(sigmas(j)) ]);
%     end
    
    [~, Lambda2] = imageEigenvalues(I,sigmas(j),spacing,brightondark);
    if brightondark == true
        Lambda2 = -Lambda2;
    end  
    
    % proposed filter at current scale
    Lambda3 = Lambda2;
    
    Lambda_rho = Lambda3;
    Lambda_rho(Lambda3 > 0 & Lambda3 <= tau .* max(Lambda3(:))) = tau .* max(Lambda3(:));
    Lambda_rho(Lambda3 <= 0) = 0;
    response = Lambda2.*Lambda2.*(Lambda_rho-Lambda2).* 27 ./ (Lambda2 + Lambda_rho).^3;    
    
    response(Lambda2 >= Lambda_rho./2 & Lambda_rho > 0) = 1;    
    response(Lambda2 <= 0 | Lambda_rho <= 0) = 0;
    response(~isfinite(response)) = 0;   
    
    %max response over multiple scales
    if(j==1)
        vesselness = response;
    else        
        vesselness = max(vesselness,response);
    end
        
    clear response Lambda2 Lambda3
end

vesselness = vesselness ./ max(vesselness(:)); % should not be really needed   
vesselness(vesselness < 1e-2) = 0;

function [Lambda1, Lambda2] = imageEigenvalues(I,sigma,spacing,brightondark)
% calculates the two eigenvalues for each voxel in a volume

% Calculate the 2D hessian
[Hxx, Hyy, Hxy] = Hessian2D(I,sigma,spacing);

% Correct for scaling
c=sigma.^2;
Hxx = c*Hxx; 
Hxy = c*Hxy;
Hyy = c*Hyy;

B1 = - (Hxx+Hyy);
B2 = Hxx .* Hyy - Hxy.^2;

T = ones(size(B1));

if brightondark == true
    T(B1<0) = 0;
    T(B2==0 & B1 == 0) = 0;
else
    T(B1>0) = 0;
    T(B2==0 & B1 == 0) = 0;
end

clear B1 B2;

indeces = find(T==1);

Hxx = Hxx(indeces);
Hyy = Hyy(indeces);
Hxy = Hxy(indeces);

% Calculate eigen values
[Lambda1i,Lambda2i]=eigvalOfHessian2D(Hxx,Hxy,Hyy);

clear Hxx Hyy Hxy;

Lambda1 = zeros(size(T));
Lambda2 = zeros(size(T));

Lambda1(indeces) = Lambda1i;
Lambda2(indeces) = Lambda2i;

% some noise removal
Lambda1(~isfinite(Lambda1)) = 0;
Lambda2(~isfinite(Lambda2)) = 0;

Lambda1(abs(Lambda1) < 1e-4) = 0;
Lambda2(abs(Lambda2) < 1e-4) = 0;


function [Dxx, Dyy, Dxy] = Hessian2D(I,Sigma,spacing)


if nargin < 3, Sigma = 1; end

if(Sigma>0)
    F=imgaussian(I,Sigma,spacing);
else
    F=I;
end

% Create first and second order diferentiations
Dy=gradient2(F,'y');
Dyy=(gradient2(Dy,'y'));
clear Dy;

Dx=gradient2(F,'x');
Dxx=(gradient2(Dx,'x'));
Dxy=(gradient2(Dx,'y'));
clear Dx;

function D = gradient2(F,option)

% Fx = gradient2(F,'x');

[k,l] = size(F);
D  = zeros(size(F),class(F)); 

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:) = (F(2,:) - F(1,:));
    D(k,:) = (F(k,:) - F(k-1,:));
    % Take centered differences on interior points
    D(2:k-1,:) = (F(3:k,:)-F(1:k-2,:))/2;
case 'y'
    D(:,1) = (F(:,2) - F(:,1));
    D(:,l) = (F(:,l) - F(:,l-1));
    D(:,2:l-1) = (F(:,3:l)-F(:,1:l-2))/2;
otherwise
    disp('Unknown option')
end
        
function I=imgaussian(I,sigma,spacing,siz)


if(~exist('siz','var')), siz=sigma*6; end

if(sigma>0)

    % Filter each dimension with the 1D Gaussian kernels\
    x=-ceil(siz/spacing(1)/2):ceil(siz/spacing(1)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(1))^2)));
    H = H/sum(H(:));    
    Hx=reshape(H,[length(H) 1]);
    
    x=-ceil(siz/spacing(2)/2):ceil(siz/spacing(2)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(2))^2)));
    H = H/sum(H(:));    
    Hy=reshape(H,[1 length(H)]);
    
    I=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
end

function [Lambda1,Lambda2]=eigvalOfHessian2D(Dxx,Dxy,Dyy)
% This function calculates the eigen values from the
% hessian matrix, sorted by abs value

% Compute the eigenvectors of J, v1 and v2
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);

% Compute the eigenvalues
mu1 = 0.5*(Dxx + Dyy + tmp);
mu2 = 0.5*(Dxx + Dyy - tmp);

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
check=abs(mu1)>abs(mu2);

Lambda1=mu1; Lambda1(check)=mu2(check);
Lambda2=mu2; Lambda2(check)=mu1(check);

