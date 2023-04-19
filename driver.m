%%Driver code. We assume that n=2, inner radius = 5h.
clear; close all;
clc
format long e

nx = 50;% number of nodes along real axis
h = .1;% node spacing

bounds = [ -4,4; -4,4; 0,4]; % Lower and upper in the three display
al = 0.1; %Order of the fractional derivative, 0<alpha<1

bx = [-nx*h/2,nx*h/2,-nx*h/2,nx*h/2]; % Domain to be displayed
nxy=sign(bx).*ceil(abs(bx)/h); % Number of nodes inside the domain

f = @(x) sinh(x); %Function to differentiate and its true derivative
true_D_alp_f = @(x,al) x.^(1-al).*hypergeom([1],[1-al/2,(3-al)/2],x.^2/4)/gamma(2-al);

Df = compute_frac_der(f,al,nxy,h);
display_function(Df,true_D_alp_f,al,nxy,h,bounds); % Create the three subplots

