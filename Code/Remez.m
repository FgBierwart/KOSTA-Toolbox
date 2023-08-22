
clear 
close   
clc 

addpath(genpath('C:\Program Files\YALMIP'))
addpath(genpath('C:\Users\fgbierwa\OneDrive - Université de Namur\Documents\Doctorat - Recherche\SOSTOOLS'))
addpath(genpath('C:\Users\fgbierwa\OneDrive - Université de Namur\Documents\Doctorat - Recherche\SeDuMi'))

% Construction of the Psi_i (n functions) 
syms x1 x2; x = [x1;x2]; 

% 1. Monomials 

order = 12; 
z = monomials([x1;x2],1:order);

% 2. Chebyshev polynomials (here for 2D)

% order = 8;  
% z = chebyshevT(0:order,x(1)); 
% for i = 2:2
%     z = kron(z,chebyshevT(0:order,x(i))); 
% end

n = length(z); 
a = 1;  

% Function we wnat to approximate 

% K = 0.2; 
% f = [K*sin(x1-x2)-sin(x1),K*sin(x2-x1)-sin(x2)]; 

% f = 0.2*sin(x1-x2)-sin(x1); 

f = [-(x1+x2)/sqrt(1+(x1+x2)^2)]; 

w = 4; 
f = (1/w)*subs(f,x,w*x);

c = (a-(-a))*rand(2,100000)-a;

for i = 1:length(f)
    [alpha(i,:),~] = Remez_alg(f(i),n,z,a);
    tmp = matlabFunction(abs(f(i)-value(alpha(i,:))*z)); 
    % tmp = matlabFunction(abs(f(i)-value(alpha(i,:))*transpose(z))); 
    h(i) = max(reshape(tmp(c(1,:),c(2,:)),1,[])); 
end

%%

% F_sym = value(alpha)*z;

F_sym(1) = x2; 
F_sym(2) = value(alpha)*z;

%%
clc
close;       
xr{1} = linspace(-a,a,1000); 
xr{2} = linspace(-a,a,1000); 
c = Grid(xr); 
c = flip(c,1); 
tmp = matlabFunction(abs(f-value(alpha)*z));  
disp(max(reshape(tmp(c(1,:),c(2,:)),1,[]))); 
disp(value(h)); 

%%

close all
[X,Y] = meshgrid(xr{1},xr{2});

f_app = matlabFunction(value(alpha)*transpose(z)); 
fct = matlabFunction(f); 

%%

b = surf(X,Y,f_app(X,Y)); 
set(b,'FaceColor',[0 0 1], ...
      'FaceAlpha',0.5,'EdgeColor','none')

hold on

b = surf(X,Y,fct(X,Y)); 
set(b,'FaceColor',[1 0 0], ...
      'FaceAlpha',0.5,'EdgeColor','none')   

%%

% b = surf(X,Y,abs(f_app(X,Y)-fct(X,Y))); 
% set(b,'FaceColor',[0 0 1], ...
%       'FaceAlpha',0.5,'EdgeColor','none')

hold on

b = surf(X,Y,value(h)*(X.^2+Y.^2)); 
set(b,'FaceColor',[1 0 0], ...
      'FaceAlpha',0.5,'EdgeColor','none')
