function Df = compute_frac_der(f,al,nxy,h)
%%INPUT: f: function to differentiate
%%al : order of the fractional derivative (0<al<1)
%%nxy : number of nodes in all directions.
%%h: interval in x and in y
%%OUTPUT: Df: alpha^th derivative (Caputo) of f(z) with base point at z=0
x = h*(nxy(1)-2:nxy(2)+2); y = h*(nxy(3)-2:nxy(4)+2);
[xr,xi] = meshgrid(x,y(end:-1:1)); z = xr+1i*xi; fval = f(z);

ind_z_list = find(real(z(:))>=h*nxy(1) & real(z(:))<=h*nxy(2) & ...
    imag(z(:))>=h*nxy(3) & imag(z(:))<=h*nxy(4));%Nodes inside (display) domain
[row,col] = ind2sub(size(z),ind_z_list(:));%coord shift to display domain

D = sparse(length(ind_z_list),length(z(:))); %Initializes global diff matrix
Df = zeros(length(ind_z_list(:)),1);%Initializes derivative approx
w_reg = stencil(-1); w_sing = stencil(al);%Corr stencils (reg and sing)
cntr_x=3-nxy(1); cntr_y=3+nxy(4);%Position of the origin in the matrix

for ind_z=1:length(ind_z_list) %For each evaluation point z0 = x0 + i y0
    x0 = xr(row(ind_z),col(ind_z)); y0 = xi(row(ind_z),col(ind_z));

    if (sqrt(x0.^2+y0.^2)<=10*h) %Taylor Expansion using midpoint method
        t=linspace(0,2*pi,26)'; t(end)=[];
        ind = ismember(z(:),h*(round(x0/(2*h)+5*cos(t))+1i*round(y0/(2*h)+5*sin(t))));
        D(ind_z_list(ind_z),ind)=weights_midpoint(al,x0 + 1i*y0,z(ind)-(x0 + 1i*y0)/2);

    else %End Correction method
        W = zeros(size(z)); Ws = zeros(size(z));%reg and sing weights matrices
        sing_vals = (1./(x0 + 1i*y0-z).^(al+1)); %singular term
        sing_vals(row(ind_z),col(ind_z))=0;%Zeroes out the NaN at eval point

        if x0<0 & y0>=0 %shifts position of branch cut on left side of C
            sing_vals(imag(z)>y0) = exp(-1i*2*pi*al)*sing_vals(imag(z)>y0);
        elseif x0<0 & y0<0
            sing_vals(imag(z)<=y0) = exp(1i*2*pi*al)*sing_vals(imag(z)<=y0);
        end

        %Chooses a path for the particular evaluation point x0 + 1i y0
        ind_x0 = round(x0/h);ind_y0 = round(y0/h);%Eval point coordinates
        if ind_x0==0 | ind_y0==0
            path = [0 0; ind_x0 ind_y0];
        elseif abs(ind_x0)<abs(ind_y0)
            path = [0 0; ind_x0 0; ind_x0 ind_y0];
        elseif abs(ind_y0)>0
            path = [0 0; 0 ind_y0; ind_x0 ind_y0];
        end

        W(cntr_y,cntr_x) =(x0 + 1i*y0)/al;

        for idx_pth = 2:size(path,1) %For each path segment
            %Gets path direction in complex plane, indices ofapth nodes, etc.
            stpx = sign(path(idx_pth,1)-path(idx_pth-1,1));dx = stpx+(stpx==0);
            stpy = sign(path(idx_pth,2)-path(idx_pth-1,2));dy = stpy+(stpy==0);
            path_yi=cntr_y-(path(idx_pth-1,2):dy:path(idx_pth,2));
            path_xi=cntr_x+(path(idx_pth-1,1):dx:path(idx_pth,1));
            dir = stpx+1i*stpy;  x_ep = (-2:2);

            tr_path =  h*dir*ones(size(z(path_yi,path_xi)));%TR weights
            tr_path(1)=0;tr_path(end)=0; %TR endpoints
            W(path_yi,path_xi)=W(path_yi,path_xi)+tr_path;%Records TR weights

            dir_ind_strt = find([1 -1 1i -1i]==dir);%dir of start and end of
            dir_ind_end = find([-1 1 -1i 1i]==dir);%path for correct stencil

            W(path_yi(1)+x_ep,path_xi(1)+x_ep)=...% Gets reg start corrections
                W(path_yi(1)+x_ep,path_xi(1)+x_ep)+h*dir*w_reg(:,:,dir_ind_strt);

            if idx_pth == size(path,1) % Get final (singular) end correction
                Ws(path_yi(end)+x_ep,path_xi(end)+x_ep)=Ws(path_yi(end)+...
                    x_ep,path_xi(end)+x_ep)-(h^-al)*exp(-1i*angle(dir)*...
                    (1-2*(angle(dir)==pi & y0<0))*al)*(w_sing(:,:,dir_ind_end));

            else % Gets regular end corrections
                W(path_yi(end)+x_ep,path_xi(end)+x_ep)=...
                    W(path_yi(end)+x_ep,path_xi(end)+x_ep)+h*dir*w_reg(:,:,dir_ind_end);
            end
        end
        % Records matrices W and Ws into the j^th row of global Diff Mat D
        [~,Dpos,Dval] = find((-al/gamma(1-al))*(W(:).*sing_vals(:)).');
        [~,Dpos_s,Dval_s] = find((-al/gamma(1-al))*(Ws(:)).');
        D(ind_z_list(ind_z),Dpos)=Dval;
        D(ind_z_list(ind_z),Dpos_s)=D(ind_z_list(ind_z),Dpos_s)+Dval_s;
    end
end
Df = D(ind_z_list,:)*fval(:);
end