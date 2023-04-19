function W = stencil(al)
%%Computes the endpoint correction stencils (n is set to 2)
%INPUTS: al: Order of the fractrional derivative
%OUTPUTS: W: weights in the stencil (If al=-1, the output is the regular 
% (non-singular) endpoint stencil with 1/2 weight at the central node)
n=2;

x_ep = (-n:n);
za = x_ep-1i*x_ep.'; zk = za(:).'; % Stencil nodes

%Calculates stencil weights
A = zk.^((0:(2*n+1)^2-1).');IA = inv(A);
vec = (0:(2*n+1)^2-1); zet = zeta(al+1-vec(:));
vS = ([1 -1 1i -1i].^vec(:)).*zet(:,[1 1 1 1]);
cS = IA*vS;
W = sign(al)*reshape(cS,2*n+1,2*n+1,4);
end