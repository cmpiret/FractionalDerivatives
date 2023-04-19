function display_function(Df,D_alp_f,al,nxy,h,bounds)
%%Displays the fractional derivative of function f(z)
%%Modified version to display real, imaginary, magnitude of Df, and error
lw = 4; %line width
bx = h*nxy; %scaled domain bounds
x = h*(nxy(1):nxy(2)); y = h*(nxy(3):nxy(4)); %computational domain
[xr,xi] = meshgrid(x,y(end:-1:1)); z = complex(xr,xi); Z=z(:);
cntr_x=-nxy(1)+1; cntr_y=nxy(4)+1; %central node

Df = reshape(Df,-nxy(3)+nxy(4)+1,-nxy(1)+nxy(2)+1);%reshapes approx derivative
Df_below = Df(cntr_y:end,:); Df_below(1,:)=conj(Df_below(1,:));
Df_top = Df(1:cntr_y,:);

true = D_alp_f(z,al);%true derivative

error=abs(Df-true)./abs(true);%Relative error
error(error==0)=eps; error(isinf(error))=NaN; %gets rid of empty error pixels
max_error = max(error(:));

ax=figure(1); hold on

subplot(2,2,1) % Plots the real part using mesh
mesh(xr(cntr_y:end,:),xi(cntr_y:end,:),real(Df_below), 'edgecolor', 'k');
hold on;
mesh(xr(1:cntr_y,:),xi(1:cntr_y,:),real(Df_top), 'edgecolor', 'k');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
set( gca,'FontSize', 16);
xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(1,:));
plot3(xr(cntr_y,:),xi(cntr_y,:),real(Df(cntr_y,:)),'r','LineWidth',lw); % Highlights real axis
plot3(xr(cntr_y,:),xi(cntr_y,:),real(Df_below(1,:)),'r','LineWidth',lw);
title(['$$Re(D^{\alpha} f)$$'],'interpreter','latex')

subplot(2,2,2) % Plots the imaginary part using mesh
mesh(xr(cntr_y:end,:),xi(cntr_y:end,:),imag(Df_below), 'edgecolor', 'k');
hold on;
mesh(xr(1:cntr_y,:),xi(1:cntr_y,:),imag(Df_top), 'edgecolor', 'k');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
set( gca,'FontSize', 16);
xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(2,:));
plot3(xr(cntr_y,:),xi(cntr_y,:),imag(Df(cntr_y,:)),'r','LineWidth',lw); % Highlight real axis
plot3(xr(cntr_y,:),xi(cntr_y,:),imag(Df_below(1,:)),'r','LineWidth',lw);
title(['$$Im(D^{\alpha} f)$$'],'interpreter','latex')


subplot(2,2,3) % Plot the magnitude as a surface
p = surf(xr,xi,abs(Df),angle(-Df)); % Display the surface
set (p,'EdgeColor','none'); colormap hsv(600); hold on;
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
title(['$$|D^{\alpha} f|$$'],'interpreter','latex')
set( gca,'FontSize', 16);
xlim(bx(1:2)); ylim(bx(3:4)); zlim(bounds(3,:));
plot3(xr(cntr_y,:),xi(cntr_y,:),abs(Df(cntr_y,:)),'k','LineWidth',lw); % Highlight real axis
axes('Position',[0.05 0.05 .17 .17]) % Add color wheel for phase information
[th,r] = meshgrid(linspace(-pi,pi),linspace(0,1));
[X,Y] = pol2cart(th+pi,r);

contourf(X,Y,th,100,'linestyle','none'); hold on % Show colors wheel
plot([-1 1],[0,0],'k'); plot([0 0],[-1,1],'k'); % Show Re and Im axes
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k'); colormap hsv(600);
axis equal; axis off

ax4=subplot(2,2,4)
[M,c] = contourf(real(z),imag(z),log10(abs(error)));
c.LineWidth = 1;
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
set( gca,'FontSize', 16);
colorbar;
title({['$$Relative \ error$$'];['$$\alpha = $$' num2str(al)]},'interpreter','latex')
hold on
plot(real(z),imag(z),'k.')
rectangle('Position',[-10*h -10*h 20*h 20*h],'Curvature',[1 1],'EdgeColor','r','LineWidth',2)
clim([-16,log10(max_error)])
colormap(ax4,"parula")
axis equal

ax.Position = [100 100 800 800];