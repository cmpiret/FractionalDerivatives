function W = weights_midpoint(al,z0,zv)
%Computes the central weights using the midpoint expansion method
%INPUTS: al: Order of the fractrional derivative
%        z0: evaluation point
%        zv: surrounding nodes at which we want to compute the weights
%OUTPUTS: W: weights

num_data = length(zv); kvec = (0:(num_data-1))'; %Initializes values
Z = repmat(zv,1,num_data); Z(:,1) = 1; Z = cumprod(Z,2);
b = z0(:)/2; Ones_d = ones(1,num_data);

dk(1) = 1./(al-1);
for k=1:num_data-1
    dk(k+1) = (k*dk(k)+1)/(al-k-1);
end

C=-((2*b(:,Ones_d)).^(1-al)).*((-b(:,Ones_d)).^(kvec(:,1)')).*(dk(1,:));
D_al_Z0 = (1/gamma(1-al))*([0,C(:,1:end-1)]*diag(kvec));

W = D_al_Z0/Z;
