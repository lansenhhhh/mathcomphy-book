%Click "Run Section"
load('assigment3.mat') % asymmetric 17 by 17
x_true=double(imread('Hubble.jpg'));
center=ceil(size(PSF)/2);
radius=ceil(max(size(PSF))/2);
sz=size(x_true);
A=Aclass(PSF,center,'reflexive',sz(1),sz(2));
b=A*x_true;
x_true=x_true(radius+1:end-radius,radius+1:end-radius);
b=b(radius+1:end-radius,radius+1:end-radius);
bc='reflexive';
[n,m]=size(b);
A=Aclass(PSF,center,bc,n,m);

% Use the generated A and b to test your GMRES algorithm.

% Here, you should initialize x0 (initial guess for x) if needed. For simplicity, we'll use zeros
x0 = zeros(n*m, 1); % Initial guess for GMRES, ensure it's a column vector

% Solving the system Ax = b using your GMRES solver. Make sure casestudies_gmres can handle Aclass
[xk, Res] = casestudies_gmres(A, b, 0.03, 50); % Adjusted call to include initial guess x0

% Store the relative residual norm of each iteration in a vector called Res
% The k-th element of Res corresponds to the k-th residual norm.
% Store the estimate when the algorithm is terminated as xk.
%% plotting residual norm
figure
plot(Res,'LineWidth',2,'Color',[0 0 0])
hold on
plot(1e-1*ones(maxdim,1),'k--')
plot(1e-2*ones(maxdim,1),'k--')
plot(1e-3*ones(maxdim,1),'k--')
plot(1e-4*ones(maxdim,1),'k--')
hold off
xlabel('iteration')
ylabel('relative residual norm')
%% plotting restored image
figure
imshow(reshape(xk,n,n),[])