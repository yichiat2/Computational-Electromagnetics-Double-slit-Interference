format long;
clear all;

% Define Physical Constants
e0 = 8.85E-12; % F/m
mu0 = 4*pi*1E-7; % H/m
c = 1/sqrt(e0*mu0); % m/s

% Flags for activating different functionalities
flagPML = 1; % PML Layer Activation
flagPMC = 0; % PMC Layer Activation
flagGS  = 1; % Source type: modulated Gaussian pulse (0) or single-frequency source (1) 

% Device and Mesh Parameters
xLen = 200E-6; % Box length in x: 200 um
yLen = 200E-6; % Box length in y: 200 um
NGX = 200;     % number of x grids  
NGY = 200;     % number of y grids
NGT = 200;     % number of t grids
NGPML = 20;    % number of grids assigned to PML layers
sorx = NGX/2;  % The x coordinate of the source 
sory = NGY/2;  % The y coordinate of the source 
NGXS2PMC = 20; % number of grids between the source and the PMC layer along the x direction
NGXPMC = 10;   % number of grids of the PMC layer along the x direction
NGYBAR = 4;    % length of the PMC bar defined in the number of grids along the y direciton
NGYSLT = 1;    % width of each slit defined in the number of grids along the y direciton 
NGYPMC = 1/2 * (NGY - (2*NGYSLT+NGYBAR)); % length of each PMC fence defined in the number of grids along the y direciton
dx = xLen/ NGX;  % individual x-grid length in m
dy = yLen/ NGY;  % individual y-grid length in m 
dt =  1/(c*sqrt(  (1/dx^2)  + (1/dy^2) )); % individual t-grid length in s 
lenPMLx = NGPML * dx;  % PML thickness in m along the x direction
lenPMLy = NGPML * dy;  % PML thickness in m along the y direction
 
 
% Define the PML conductivity
mPML = 3;   % m-th polynomial for PML conductivity
rhomaxx = -(mPML+1)/(2* lenPMLx *sqrt(mu0/e0) ) * log(1E-11);
rhomaxy = -(mPML+1)/(2* lenPMLy *sqrt(mu0/e0) ) * log(1E-11);

% Matrices Initialization
Ezx = zeros(NGT,NGX,NGY);   % the x contribution of Ez field
Ezy = zeros(NGT,NGX,NGY);   % the y contribution of Ez field
Ez = zeros(NGT,NGX,NGY);    % the total contribution of Ez field
Hx = zeros(NGT,NGX,NGY);    % the total contribution of Hx field
Hy = zeros(NGT,NGX,NGY);    % the total contribution of Hy field
S = zeros(NGT,NGX,NGY);		% source stored in 3D (2D in space + 1D in time )
PMCAREA = zeros(NGX,NGY);   % 2D array of flags used to mark the boundary conditions of PMC regions

e = zeros(NGX,NGY);         % e stands for epsilon (permitivity) 
mu = zeros(NGX,NGY);        % mu stands for permeability
rhox = zeros(NGX,NGY);      % conductivity in the x direction
rhoy = zeros(NGX,NGY);      % conductivity in the y direction

e(:,:)  = e0;				% initialization (vacuum)
mu(:,:) = mu0;              % initialization (vacuum)


% Setup the conductivity of the PML layers
if flagPML == 1
    for i=1:NGPML
        rhox(i,:) = rhomaxx * ( (NGPML+1-i)/NGPML )^mPML;
        rhoy(:,i) = rhomaxy * ( (NGPML+1-i)/NGPML )^mPML;
    end
    for i=NGX-NGPML+1:NGX
        rhox(i,:) = rhomaxx * ( (i-NGX+NGPML)/NGPML )^mPML;
    end
    for i=NGY-NGPML+1:NGY
        rhoy(:,i) = rhomaxy * ( (i-NGY+NGPML)/NGPML )^mPML;
    end
    %contourf(rhox); 
    %contourf(rhoy);
end

% Mark the boundary condition that conductivity of the PML layers
if flagPMC == 1
    
    xini = sorx + NGXS2PMC;        % define the x indices of the PMC layers
    xend = xini + NGXPMC;    
    yslt1ini = NGYPMC;             % define the y indices of the first slit
    yslt1end = yslt1ini + NGYSLT; 
    yslt2ini = yslt1end + NGYBAR;  % define the y indices of the second slit
    yslt2end = yslt2ini + NGYSLT;
     
	% The surfaces perpendicular to the x direction is marked as 1, meaning that Hy (i,j) = 0
    PMCAREA(xini, 1: yslt1ini) = 1; 
    PMCAREA(xini, yslt1end+1: yslt2ini) = 1; 
    PMCAREA(xini, yslt2end+1: yslt2end + NGYPMC) = 1;
    PMCAREA(xend, 1: yslt1ini) = 1;
    PMCAREA(xend, yslt1end+1: yslt2ini) = 1;
    PMCAREA(xend, yslt2end+1: yslt2end + NGYPMC) = 1;
    
	% The surfaces perpendicular to the y direction is marked as -1, meaning that Hx (i,j) = 0
    PMCAREA(xini+1 : xend, 1) = -1; 
    PMCAREA(xini+1 : xend, yslt1ini) = -1;
    PMCAREA(xini+1 : xend, yslt1end) = -1;
    PMCAREA(xini+1 : xend, yslt2ini) = -1;
    PMCAREA(xini+1 : xend, yslt2end) = -1;
    PMCAREA(xini+1 : xend, yslt2end+NGYPMC) = -1;
	
	% The cornors of slits should obey Hx (i,j) = Hy (i,j) = 0, marked as 2 
	% Failed to specify these boundary conditions result in the leakages of electric fields in the PMC layers
    PMCAREA(xend, yslt1ini) = 2; 
    PMCAREA(xend, yslt2ini) = 2;
 %   contourf(PMCAREA(:,:)');
end

    

% Setup Pulse Source
tmpt0 = 30*dt; % a delay 
if flagGS == 0 
    amp =   1E6;           % pulse amplitude
    omega0 =  2*pi*3E12;   % modulation angular frequency 
    taup =   5*dt;         % dispersion parameter
elseif flagGS == 1
    amp =   1E5; 
    omega0 =  2*pi*3E13; 
    taup =   10*dt; 
end


for l=1:NGT
    tmpt = (l-1) * dt; 

    if flagGS == 0
        S(l, sorx, sory) = -amp*exp(-1/2* ((tmpt0-tmpt)/taup)^2) * sin(omega0 * (tmpt0-tmpt));
    elseif flagGS == 1
        S(l, sorx, sory) = -amp*(1-exp(-tmpt/taup))*sin(omega0 * tmpt);   
    end
end
% Plot the source as a function of time 
figure(1); 
tmp = 1:NGT; 
tmp = (tmp-1) * dt * 1E15;
plot(tmp, S(:, sorx, sory));
%title(['\fontsize{10} Modulated Gaussian Current Pulse']); 
xlabel('Time (fs)','FontSize',10);
ylabel('Amplitude (A/m^2)','FontSize',10);
    
    
    
% The main FDTD (time-stepping) algorithm (Eq. 4 in the report)
for l = 2:NGT  % index of time
    
    % Update Ez in the 2D space (Eq. 4 in the report)
    for i = 2:NGX-1
        for j = 2:NGY-1
    
        bx = e(i,j)/dt + rhox(i,j)/2; 
        by = e(i,j)/dt + rhoy(i,j)/2;
        ax = e(i,j)/dt - rhox(i,j)/2;
        ay = e(i,j)/dt - rhoy(i,j)/2;
        
        Ezx(l,i,j) = 1/bx * (ax * Ezx(l-1,i,j) + 1/dx * (Hy(l-1,i,j)-Hy(l-1,i-1,j) ));
        Ezy(l,i,j) = 1/by * (ay * Ezy(l-1,i,j) - 1/dy * (Hx(l-1,i,j)-Hx(l-1,i,j-1) ));
        Ez(l,i,j) =  Ezx(l,i,j) + Ezy(l,i,j); 
        
        if (i==sorx && j==sory ) 
            Ez(l,i,j) =  Ez(l,i,j)- 1/bx * S(l-1,i,j);
        end

        end
    end
    
    % Update Hx and Hy in the 2D space
    for i = 2:NGX-1
        for j = 2:NGY-1
    
        bx = e(i,j)/dt + rhox(i,j)/2;
        by = e(i,j)/dt + rhoy(i,j)/2;
        ax = e(i,j)/dt - rhox(i,j)/2;
        ay = e(i,j)/dt - rhoy(i,j)/2;
        
        Hx(l,i,j) = 1/by * (ay * Hx(l-1,i,j) - e(i,j)/(mu(i,j)*dy) * (Ez(l,i,j+1) - Ez(l,i,j)));
        Hy(l,i,j) = 1/bx * (ax * Hy(l-1,i,j) + e(i,j)/(mu(i,j)*dx) * (Ez(l,i+1,j) - Ez(l,i,j)));
        
		% Apply the PMC boundary conditions
        if (flagPMC == 1)
            if PMCAREA(i,j) == 1
                Hy(l,i,j) = 0;
            end
            if PMCAREA(i,j) == -1
                Hx(l,i,j) = 0;
            end
            if PMCAREA(i,j) == 2
                Hx(l,i,j) = 0;
                Hy(l,i,j) = 0;
            end
        end
        
        end 
    end
    
  %  A way to display the evolution of E and H during each time step
  
  %  %Movie type colour scaled image plot of Ez
  %  imagesc(dx*(1:NGX-1)*1e+6,(1e+6*dy*(1:NGY-1))',squeeze(Ez(l,:,:)),[-1,1]);colorbar;
  %  title(['\fontsize{20}Colour-scaled image plot of Ez in a spatial domain with PEC boundary and at time = ',num2str(round(l*dt*1e+15)),' fs']); 
  %  xlabel('x (in um)','FontSize',20);
  %  ylabel('y (in um)','FontSize',20);
  %  set(gca,'FontSize',20);
  %  getframe;

end


% Generic plot functions to display E and H at specified time steps
figure(2)
for l=1:NGT
    M = squeeze(Ez(l,:,:));
    colormap jet;
    
    surf(M','FaceAlpha',1,'edgecolor','none');
    axis([1 NGX 1 NGY -1.5 1.5]);
    caxis([-1.5,1.5]); 
    colorbar;
    view(2);
    title(['\fontsize{15} E_z (V/m) at time = ',num2str(round(l*dt*1e+15)),' fs']); 
    xlabel('x (um)','FontSize',15);
    ylabel('y (um)','FontSize',15);
    pause(0.02);
end


% correspond to Fig. 7(a) 
figure(3)
plot(squeeze(Ez(170,:,101)));
xlabel('x (um)','FontSize',15);
ylabel('E_z (V/m)','FontSize',15);

% correspond to Fig. 7(b) 
figure(4)
plot(squeeze(Ez(170,159,:)));
xlabel('y (um)','FontSize',15);
ylabel('E_z (V/m)','FontSize',15);



