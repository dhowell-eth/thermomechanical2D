% Exercise 11:
% Thermomechanical solver for variable-viscosity, buoyancy driven flows using a staggered
% grid with markers approach (now including heating terms)

% Model domain is defined as a rectangular grid with a circular zone in the
% middle which may have a diffent density and viscosity than the
% surrounding medium and an overlying layer of "sticky air".
%
% Mechanical solver assumes free-slip boundary conditions. Thermal solver
% assumes constant T at top/bottom and insulation along the left and right
% boundaries
%
% Written for Numerical Modeling with Taras, Fall 2017
% Dorran Howell
% dhowell@student.ethz.ch

%clf
clear

%% --------------------DEFINE MODEL PARAMETERS-------------------------%
Nx = 41;                        % # Grid nodes x direction 41x31 for turnin
Ny = 31;                        % # Grid nodes y direction
nxm = (Nx-1)*4;                 % # Markers X direction
nym = (Ny-1)*4;                 % # Markers Y direction
useRandomMarkerSpacing = false; % Flag for whether to introduce random noise to marker spacing

xsize = 100000;                 % Model length x-direction [m]
ysize = 100000;                 % Model length y-direction [m]

Gy = 9.81 ;                     % Magnitude of Gravity y-component (m/s^2)

rInside = 20000;                % Interior Radius [m] 20000
rhoInside = 3200;               % Interior Density [kg/m^3] 3200
rhoOutside = 3300;              % Exterior Density [kg/m^3] 3300
rhoAir = 1;                     % Sticky Air Density [kg/m^3] 3300

etaInside = 1e19;               % Interior Viscocity [kg/(s*m)]
etaOutside = 1e21;              % Exterior Viscosity [kg/(s*m)]
etaAir = 1e17;                  % Sticky Air Viscosity [kg/(s*m)]

kInside = 2;                    % Interior Thermal conductivity [W/(m*K)]
kOutside = 3;                   % Exterior Thermal conductivity [W/(m*K)]
kAir = 3000;                    % Sticky Air Thermal conductivity [W/(m*K)]

pcpInside = 3200000;            % Interior Volumetric Heat Capacity [J/(kg*K)]
pcpOutside = 3300000;           % Exterior Volumetric Heat Capacity [J/(kg*K)]
pcpAir = 3300000;               % Sticky Air Volumetric Heat Capacity [J/(kg*K)]

tempInside = 1700;              % Interior Initial Temperature [K]
tempOutside = 1500;             % Exterior Initial Temperature [K]
tempAir = 273;                  % Sticky Air Initial Temperature [K]

alphaInside = 2e-5;             % Interior Initial Temperature [K]
alphaOutside = 3e-5;            % Exterior Initial Temperature [K]
alphaAir = 0;                   % Sticky Air Initial Temperature [K]

hRInside = 3e-8;                % Interior Thermal Expansion Coefficient [1/K]
hROutside = 2e-8;               % Exterior Thermal Expansion Coefficient [1/K]
hRAir = 0;                      % Sticky Air Thermal Expansion Coefficient [1/K]

bcPressure = 1e9;               % BC Pressure @ index (i=2,j=3) [Pa];
bcMode = 'freeSlip';            % BC Mode: 'freeSlip' or 'noSlip'
bcTempTop = 273;                % Temperature BC @ Top Boundary [K]
bcTempBottom = 1500;            % Temperature BC @ Bottom Boundary [K]

nSteps = 10;                    % # of Iterations
dt = 0;                         % Initial timestep [s] 0e12
dTMax = 50;                     % Max T change per timestep [K]
dxMax = 0.5;                    % Maximum change in x allowed for a timestep (fraction of grid spacing)
dyMax = 0.5;                    % Maximum change in y allowed for a timestep (fraction of grid spacing)
%----------------------------------------------------------------------%

%% Model Setup
% Check BC input and set coefficients
if strcmp(bcMode,'freeSlip')
    bcLeft = [1,-1];      % 1*vy1 + -1*vy2 = 0
    bcRight = [1,-1];     % 1*vy1 + -1*vy2 = 0
    bcTop = [1,-1];       % 1*vx1 + -1*vx2 = 0
    bcBottom = [1,-1];    % 1*vx1 + -1*vx2 = 0
    bcVelocities = 1;
elseif strcmp(bcMode,'noSlip')
    bcLeft = [1,-1/3];    % 1*vy1 + -1*vy2 = 0
    bcRight = [1,-1/3];   % 1*vy1 + -1*vy2 = 0
    bcTop = [1,-1/3];     % 1*vx1 + -1*vx2 = 0
    bcBottom = [1,-1/3];  % 1*vx1 + -1*vx2 = 0
    bcVelocities = -1;
else
    fprintf('ERROR: Unrecognized bcMode value. Must be one of {''freeSlip'', ''noSlip''}\n');
    return
end
% Determine grid spacing
dx = xsize/(Nx-1);
dy = ysize/(Ny-1);
% Determine marker spacing
dxm = xsize/nxm;
dym = ysize/nym;
% Populate nodal geometrical arrays
x = 0:dx:xsize;
y = 0:dy:ysize;
% Positions of staggered grid nodes
xVx = 0:dx:xsize;
yVx = -dy/2:dy:ysize-dy/2;
xVy = -dx/2:dx:xsize-dx/2;
yVy = 0:dy:ysize;
xP = -dx/2:dx:xsize+dx/2;
yP = -dy/2:dy:ysize+dy/2;
% Initialize global matrices
N = Nx*Ny*3; % Number of unknowns
Nm = nxm*nym; % Number of markers

L = sparse(N,N);
R = zeros(N,1);
LT = sparse(N,N);
RT = zeros(N,1);

% Intialize geometrical solution arrays for nodes and markers
vxNodes = zeros(Ny,Nx);
vyNodes = zeros(Ny,Nx);
pressureNodes = zeros(Ny,Nx);
tempNodes = zeros(Ny,Nx);
dTempNodes = zeros(Ny,Nx);
kNodes = zeros(Ny,Nx);
pcpNodes = zeros(Ny,Nx);

alphaNodes = zeros(Ny,Nx); % NEW - HW11
hRNodes = zeros(Ny,Nx);  % NEW - HW11
hANodes = zeros(Ny,Nx);  % NEW - HW11
hSNodes = zeros(Ny,Nx);  % NEW - HW11

vxMarkers = zeros(Nm,1);
vyMarkers = zeros(Nm,1);
rhoMarkers = zeros(Nm,1);
etaMarkers = zeros(Nm,1);
tempMarkers = zeros(Nm,1);
dTempMarkers = zeros(Nm,1);
kMarkers = zeros(Nm,1);
pcpMarkers = zeros(Nm,1);

alphaMarkers = zeros(Nm,1); % NEW - HW11
hRMarkers = zeros(Nm,1);  % NEW - HW11

% Initialize Arrays for Velocities @ Pressure Nodes
vxP = zeros(length(yP),length(xP));
vyP = zeros(length(yP),length(xP));

% Initialize arrays for various stress/strain terms to be used in adiabatic/shear
% heating calculations and later:
% Terms @ P nodes
sigmaXXP = zeros(Ny,Nx);
sigmaXYP = zeros(Ny,Nx);
sigmaYYP = zeros(Ny,Nx);
epsilonXXP = zeros(Ny,Nx);
epsilonXYP = zeros(Ny,Nx);
epsilonYYP = zeros(Ny,Nx);
% Only storing product terms at basic grid nodes (for shear heating calculations)
sigmaEpsilonXXnodes = zeros(Ny,Nx);
sigmaEpsilonXYnodes = zeros(Ny,Nx);
sigmaEpsilonYYnodes = zeros(Ny,Nx);

%% Populate initial marker positions and material properties
xm = zeros(Nm,1);
ym = zeros(Nm,1);
m = 1; % marker index
for j=1:1:nxm
    for i=1:1:nym
        
        if useRandomMarkerSpacing
            % !!! Set position (introducing some random noise to the locations)
            xm(m) = (j-0.5) * dxm + (rand-0.5)*dxm;
            ym(m) = (i-0.5) * dym + (rand-0.5)*dym;
        else
            % !!! Set position (introducing some random noise to the locations)
            xm(m) = (j-0.5) * dxm;
            ym(m) = (i-0.5) * dym;
        end
        
        % Upper 20% of Model is in "Air"
        if (ym(m) < 0.2*ysize)
           rhoMarkers(m) = rhoAir;
           etaMarkers(m) = etaAir;
           kMarkers(m) = kAir;
           pcpMarkers(m) = pcpAir;
           tempMarkers(m) = tempAir;
           alphaMarkers(m) = alphaAir; % NEW - HW11
           hRMarkers(m) = hRAir; % NEW - HW11
        % Rest follows the planet geometry of HW7
        else
            % Set material properties based on location within model
            rCurrent = sqrt((xsize/2 - xm(m))^2.0 + (ysize/2 - ym(m))^2.0);
            if (rCurrent > rInside)
                rhoMarkers(m) = rhoOutside;
                etaMarkers(m) = etaOutside;
                kMarkers(m) = kOutside;
                pcpMarkers(m) = pcpOutside;
                tempMarkers(m) = tempOutside;
                alphaMarkers(m) = alphaOutside; % NEW - HW11
                hRMarkers(m) = hROutside; % NEW - HW11
                
            else
                rhoMarkers(m) = rhoInside;
                etaMarkers(m) = etaInside;
                kMarkers(m) = kInside;
                pcpMarkers(m) = pcpInside;
                tempMarkers(m) = tempInside;
                alphaMarkers(m) = alphaInside; % NEW - HW11
                hRMarkers(m) = hRInside; % NEW - HW11
            end
        end
        % Move to next marker
        m = m + 1;
    end 
end

%% Additional setup
% Initialize pressure scaling factor (to try to eliminate precision error in inverse func.)
pscale = etaOutside/dx;

%% ------------------ BEGIN TIME LOOP -----------------------
try
currentTime = 0;                                % Current time [s]
testVals = [];                                  % Diagnostic values for HW grading

for tIter =  1:1:nSteps

% Interpolate Material Properties from markers to nodes
RHO = interpMarkersToNodes2D(rhoMarkers,xm,ym,x,y,dx,dy);
ETA = interpMarkersToNodes2D(etaMarkers,xm,ym,x,y,dx,dy);
kNodes = interpMarkersToNodes2D(kMarkers,xm,ym,x,y,dx,dy);
pcpNodes = interpMarkersToNodes2D(pcpMarkers,xm,ym,x,y,dx,dy);
% TODO [done] - add new properties from HW 11
alphaNodes = interpMarkersToNodes2D(alphaMarkers,xm,ym,x,y,dx,dy);
hRNodes = interpMarkersToNodes2D(hRMarkers,xm,ym,x,y,dx,dy);

% Note: Using a seperate routine for interpolating temperature (to ensure energy
% is conserved)
tempNodes = heatInterpMarkersToNodes2D(tempMarkers,pcpMarkers,xm,ym,x,y,dx,dy);
% Apply BC to temperature nodes:
% Constant T at Top and Bottom
tempNodes(1,:) = bcTempTop;
tempNodes(end,:) = bcTempBottom;
% dT/dx = 0 for left and right sides
tempNodes(:,1) = tempNodes(:,2);
tempNodes(:,end) = tempNodes(:,end-1);

% Populate Global Matrices
for j=1:1:Nx
    for i=1:1:Ny
        % Get global indexing
        kvx =  ((j-1)*Ny + i - 1)*3 + 1;
        kvy = kvx + 1;
        kpm = kvx + 2;
       
        % -------------- X-Stokes Equations --------------
        % Boundary Conditions:
        if (i == 1) || (i==2) || (j==1) || (j==Nx) || (i==Ny)
            if i==1 % external pts
                L(kvx,kvx) = 1; 
                R(kvx) = 0;
            elseif j==1 % left boundary
                L(kvx,kvx) = 1; 
                R(kvx) = 0;
            elseif j==Nx % right boundary
                L(kvx,kvx) = 1; 
                R(kvx) = 0;
            elseif i==2 % top boundary
                L(kvx,kvx) = bcTop(1);
                L(kvx,kvx+3) = bcTop(2);
                R(kvx) = 0;
            elseif i==Ny % bottom boundary
                L(kvx,kvx) = bcBottom(1);
                L(kvx,kvx-3) = bcBottom(2);
                R(kvx) = 0;
            end
        % Interior Points:
        else
            % -------------- X-Stokes Stencil --------------- %
            % vx, vy - velocities
            % Na/Nb/N1/N2 - viscosity
            % P - pressure
            % p1/p2 - density
            % ** See p.90 of Numerical  Modeling Script for additional details
            % 0------*-----0-------*-----0
            % |            |             |
            % |            |             |
            % |            +vx2          |
            % |            |             |
            % |      vy1   |p1,N1 vy3    |
            % 0------*-----0-------*-----0
            % |            |             |
            % | P/Na #     +vx3    # P/Nb| 
            % |            |             |
            % |      vy2   |      vy4    |
            % 0------*----(i,j)----*-----0
            % |            | p2,N2       |
            % |            |             |
            % |            +vx4          |
            % |            |             |
            % |            |             |
            % 0------------0-------------0

            % Get viscosity values

            etaA = 1 / ((1/ETA(i,j) + 1/ETA(i-1,j) +  1/ETA(i-1,j-1) + 1/ETA(i,j-1)) * (1/4));
            etaB = 1 / ( ( 1/ETA(i,j) + 1/ETA(i-1,j) +  1/ETA(i-1,j+1) + 1/ETA(i,j+1) ) * (1/4) );           
            eta1 = ETA(i-1,j); 
            eta2 = ETA(i,j);
            
            % Populate vx terms
            L(kvx,kvx - Ny*3) = 2*etaA / (dx^2.0);               %vx1
            L(kvx,kvx - 3) = eta1 / (dy^2.0);                    %vx2
            L(kvx,kvx) =  (-2*etaA / (dx^2.0)) + (-2*etaB / (dx^2.0)) + (-1*eta2/(dy^2.0)) + ...
                          (-1*eta1/(dy^2.0));                    %vx3
            L(kvx,kvx + 3) = eta2 / (dy^2.0)  ;                  %vx4
            L(kvx,kvx + Ny*3) = 2*etaB / (dx^2.0);               %vx5
            
            % Populate vy terms
            L(kvx,kvy-3) = eta1 / (dx*dy);                       %vy1
            L(kvx,kvy) = -1*eta2 / (dx*dy);                      %vy2
            L(kvx,kvy + Ny*3 - 3) = -1*eta1 / (dx*dy);           %vy3
            L(kvx,kvy + Ny*3) = eta2 / (dx*dy);                  %vy4
            
            % Populate P terms
            L(kvx,kpm) = pscale/dx;                              %P1
            L(kvx,kpm + Ny*3) = -pscale/dx;                      %P2
            % RHS
            R(kvx) = 0;
        end
       % -----------------END X-Stokes--------------------
       
       % ------------------ Y-Stokes Equations ----------------------------
        % -------------- Y-Stokes Stencil --------------- %
        % vx, vy - velocities
        % Na/Nb/N1/N2 - viscosity
        % P - pressure
        % p1/p2 - density
        % ** See p.91 of Numerical  Modeling Script for additional details
        % 0------------0-------*-----0-------------0
        % |            |      vy2    |             |
        % |            |             |             |
        % |            +vx1   #P/Na  +vx3          |
        % |            |             |             |
        % |      vy1   |      vy3    |      vy5    |
        % 0------*-----0-------*-----0-------*-----0
        % |            |        (i,j)|             |
        % |            |             |             | 
        % |            +vx2   #P/Nb  +vx4          |
        % |            |             |             |
        % |            |      vy4    |             |
        % 0------------0-------*-----0-------------0

       % Boundary Conditions:
        if (i == 1) || (j==2) || (j==1) || (j==Nx) || (i==Ny)
            if j==1 % external pts
                L(kvy,kvy) = 1; 
                R(kvy) = 0;
            elseif i==1 % top boundary
                L(kvy,kvy) = 1;
                R(kvy) = 0;
            elseif i==Ny % bottom boundary
                L(kvy,kvy) = 1;
                R(kvy) = 0;
            elseif j==2 % left boundary
                L(kvy,kvy) = bcLeft(1); 
                L(kvy,kvy+Ny*3) = bcLeft(2);
                R(kvy) = 0;
            elseif j==Nx % right boundary
                L(kvy,kvy) = bcRight(1); 
                L(kvy,kvy-Ny*3) = bcRight(2);
                R(kvy) = 0;
            end
        % Interior Points:
        else
            % Get viscosity values
            etaA = 1 / ( ( 1/ETA(i,j) + 1/ETA(i-1,j) +  1/ETA(i-1,j-1) + 1/ETA(i,j-1) ) * (1/4) );
            etaB = 1 / ( ( 1/ETA(i,j) + 1/ETA(i+1,j) +  1/ETA(i+1,j-1) + 1/ETA(i,j-1) ) * (1/4) ); 
            eta1 = ETA(i,j-1);
            eta2 = ETA(i,j);
            
            % Find dRHO/dx and dRHO/dy for diffusion correction
            % Note: Averaging the two neighboring dRho/dy values since no density
            % nodes directly line up with the vy3 node (stencil above)
            dRhoDx = (RHO(i,j) - RHO(i,j-1)) / dx;
            dRhoDy = (((RHO(i+1,j) - RHO(i-1,j))/(2*dy)) + ((RHO(i+1,j-1) - RHO(i-1,j-1))/(2*dy))) / 2  ;
            
            % Populate vx terms
            L(kvy,kvx - Ny*3) = eta1 / (dx*dy) - (dRhoDx*Gy*dt)/4;                         %vx1 - updated for HW8
            L(kvy,kvx + 3 - Ny*3) = -1*eta1 / (dx*dy) - (dRhoDx*Gy*dt)/4;                  %vx2 - updated for HW8
            L(kvy,kvx) = -1*eta2 / (dx*dy) - (dRhoDx*Gy*dt)/4;                             %vx3 - updated for HW8
            L(kvy,kvx+3) = eta2 / (dx*dy) - (dRhoDx*Gy*dt)/4;                              %vx4 - updated for HW8
   
            % Populate vy terms
            L(kvy,kvy - Ny*3) = eta1 / (dx^2.0);                        %vy1
            L(kvy,kvy - 3) = 2*etaA / (dy^2.0);                         %vy2
            L(kvy,kvy) = (-1*eta1 / (dx^2.0)) + (-1*eta2 / (dx^2.0)) + (-2*etaA / (dy^2.0)) ...
                       + (-2*etaB / (dy^2.0))  - dRhoDy*Gy*dt;                             %vy3 - updated for HW8
            L(kvy,kvy + 3) = 2*etaB / (dy^2.0);                         %vy4            
            L(kvy,kvy + Ny*3) = eta2 / (dx^2.0);                        %vy5 
            
            % P terms
            L(kvy,kpm) = pscale/dy;                                     %P1
            L(kvy,kpm+3) = -pscale/dy;                                  %P2
            % RHS
            R(kvy) = -1*Gy*( (RHO(i,j) + RHO(i,j-1)) / 2 );             
        end
        % ------------------------------------------
        
        % -------------- Continuity Equations --------------
        % Boundary Conditions:
        if (i==1) || (j==1) || (i==2 && j==2) || (i==Ny && j==Nx) || (i==2 && j==Nx) || (i==Ny && j==2) || (i==2 && j==3)
            if (i==1) || (j==1)         % External pts
                L(kpm,kpm) = 1; 
                R(kpm) = 0;
            elseif (i==2 && j==2)       % Upper left node
                L(kpm,kpm) = 1;
                L(kpm,kpm+3*Ny) = -1;
                R(kpm) = 0;
            elseif (i==Ny && j==Nx)     % Bottom right node
                L(kpm,kpm) = 1;
                L(kpm,kpm-3*Ny) = -1;
                R(kpm) = 0;
            elseif (i==2 && j==Nx)      % Upper right node
                L(kpm,kpm) = 1;
                L(kpm,kpm-3*Ny) = -1;
                R(kpm) = 0;
            elseif (i==Ny && j==2)      % Bottom Left node
                L(kpm,kpm) = 1;
                L(kpm,kpm+3*Ny) = -1;
                R(kpm) = 0;
            elseif (i==2 && j==3)       % BC Node, defined as pt. (2,3)
                L(kpm,kpm) = 1*pscale;
                R(kpm) = bcPressure;
            end
        else
            % P terms
            L(kpm,kvx - Ny*3) = -1/dx;                 %vx1
            L(kpm,kvx) = 1/dx;                         %vx2
            L(kpm,kvy-3) = -1/dy;                      %vy1
            L(kpm,kvy) = 1/dy;                         %vy2
            % RHS
            R(kpm) = 0; 
        end
        % ------------------------------------------
    end 
end

% Solve Problem
S = L\R;

% Reload into geometric arrays
for j=1:1:Nx
    for i=1:1:Ny
        % Get global indexing
        kvx =  ((j-1)*Ny + i - 1)*3 + 1;
        kvy = kvx + 1;
        kpm = kvx + 2;
        % Transfer values
        vxNodes(i,j) = S(kvx);
        vyNodes(i,j) = S(kvy);
        pressureNodes(i,j) = S(kpm)*pscale; % note: converting to real pressure
    end
end 
aaa(1,1)=vxNodes(10,10);
aaa(2,1)=vyNodes(10,10);
aaa(3,1)=pressureNodes(10,10);

% Calculate vx, vy @ P-nodes then interpolate back to markers
% Stencil:              
% 0-------*-----0
% |      vy1    |             
% |             |             
% +vx1   #P    +vx2          
% |             |            
% |      vy2    |             
% 0-------*-----0 (i,j)
%             
% vxP = avg(vx1,vx3)
% vyP = avg(vy1,vy2)

% 1 - Calculate Velocities @ interior P nodes
for j=2:1:Nx
    for i=2:1:Ny
        vxP(i,j) = (vxNodes(i,j) + vxNodes(i,j-1)) / 2;
        vyP(i,j) = (vyNodes(i,j) + vyNodes(i-1,j)) / 2;
    end
end

% 2 - Set boundary condition values for external pts.
% top/bottom
vxP(1,:) = bcVelocities * vxP(2,:);         % Vx Top
vxP(end,:) = bcVelocities * vxP(end-1,:);   % Vx Bottom
vyP(1,:) = -1 * vyP(2,:);                   % Vy Top
vyP(end,:) = -1 * vyP(end-1,:);             % Vy Bottom
% left/right
vxP(:,1) = -1 * vxP(:,2);                   % Vx Left
vxP(:,end) = -1 * vxP(:,end-1);             % Vx Right
vyP(:,1) = bcVelocities * vyP(:,2);         % Vy Left
vyP(:,end) = bcVelocities * vyP(:,end-1);   % Vy Right

% 3 - Interpolate vx, vy from P nodes to markers
vxMarkers = interpNodesToMarkers2D(vxP,xm,ym,xP,yP,dx,dy,[1,Nx,1,Ny]);
vyMarkers = interpNodesToMarkers2D(vyP,xm,ym,xP,yP,dx,dy,[1,Nx,1,Ny]);

% Determine timestep based on velocities and constraints for maximum change
% in x and y
maxVx = max(max(abs(vxNodes)));
maxVy = max(max(abs(vyNodes)));
dt = min(dxMax*dx/maxVx,dyMax*dy/maxVy);

% TODO [done] Compute various stress/strain terms (only applies to internal
% points
% Stencil:              
% 0-------*-----0 eta
% | eta  vy1    |             
% |             |       # = P node (Position we are finding sigma, epsilon values at)        
% +vx1   #P    +vx2          
% |             |            
% |             |  
% | eta  vy2    |
% 0-------*-----0 eta (i,j)
for j=2:1:Nx-1
    for i=2:1:Ny-1
        % 1. ETA() is defined at basic nodes, compute it at the P node
        % (using a harmonic average)
        etaP = 4 / (1/ETA(i,j) + 1/ETA(i-1,j) + 1/ETA(i,j-1) + 1/ETA(i-1,j-1));
        % 2. Compute values @ P nodes
        epsilonXXP(i,j) = (vxNodes(i,j) - vxNodes(i,j-1)) / dx;
        epsilonYYP(i,j) = (vyNodes(i,j) - vyNodes(i-1,j)) / dy;
        sigmaXXP(i,j) = 2 * etaP * epsilonXXP(i,j);       
        sigmaYYP(i,j) = 2 * etaP * epsilonYYP(i,j);
        % NOTE: The 2 terms below are not used or needed for HW 11 (these
        % are directly computed at the basic nodes rather than averaging)
        %epsilonXYP(i,j) = 0;        
        %sigmaXYP(i,j) = 0;          
        % 
    end 
end
% TODO [done] - Compute shear heating terms @ Basic Nodal Points
for j=2:1:Nx-1
    for i=2:1:Ny-1
        % 3. Compute values @ basic nodes
        % XY 
        thisEpsilonXY = 0.5 * ((vxNodes(i+1,j)-vxNodes(i,j))/dy + (vxNodes(i,j+1)-vxNodes(i,j))/dx);
        thisSigmaXY = 2 * ETA(i,j) * thisEpsilonXY;
        sigmaEpsilonXYnodes(i,j) =  thisSigmaXY * thisEpsilonXY;
        % XX & YY
        % First, compute sum of products for surrounding P nodes
        thisXXProductSum = sigmaXXP(i,j)*epsilonXXP(i,j) + sigmaXXP(i,j+1)*epsilonXXP(i,j+1) + ...
                           sigmaXXP(i+1,j)*epsilonXXP(i+1,j) + sigmaXXP(i+1,j+1)*epsilonXXP(i+1,j+1);
        thisYYProductSum = sigmaYYP(i,j)*epsilonYYP(i,j) + sigmaYYP(i,j+1)*epsilonYYP(i,j+1) + ...
                           sigmaYYP(i+1,j)*epsilonYYP(i+1,j) + sigmaYYP(i+1,j+1)*epsilonYYP(i+1,j+1); 
        % Then, average:
        sigmaEpsilonXXnodes(i,j) = thisXXProductSum / 4;
        sigmaEpsilonYYnodes(i,j) = thisYYProductSum / 4;

        % Finally, Compute shear heating terms
        hSNodes(i,j) = sigmaEpsilonXXnodes(i,j) + sigmaEpsilonYYnodes(i,j) + 2*sigmaEpsilonXYnodes(i,j);
    end 
end

% TODO [done] - Compute Adiabatic Heating terms
for j=2:1:Nx-1
    for i=2:1:Ny-1
        hANodes(i,j) = tempNodes(i,j) * alphaNodes(i,j) * vyNodes(i,j) * RHO(i,j) * Gy;
    end
end

for r=1:1:2
    % Prior to advecting markers, solve T equation
    for j=1:1:Nx
        for i=1:1:Ny
            % Get global index
            gkT = ((j-1)*Ny + i - 1) + 1;
            % Boundary Conditions for Temp Eq.
            if (i == 1) || (j==1) || (j==Nx) || (i==Ny)
                if (i==1) % Top
                    LT(gkT,gkT) = 1;
                    RT(gkT) = bcTempTop;
                elseif (i==Ny) % Bottom
                    LT(gkT,gkT) = 1;
                    RT(gkT) = bcTempBottom;
                elseif (j==1) % Left
                    LT(gkT,gkT) = 1; 
                    LT(gkT,gkT+Ny) = -1; 
                    RT(gkT) = 0;
                elseif (j==Nx) % Right
                    LT(gkT,gkT) = 1; 
                    LT(gkT,gkT-Ny) = -1; 
                    RT(gkT) = 0;
                end
            % Interior Points:
            else
                % STENCIL:
                %         0 T2
                %         |
                %         |
                % T1  T3_0|T3_dt
                % 0-------0-------0 T5
                %         |pcp(i,j)
                %         |
                %         |
                %         0 T4

                % Find k values
                kx1 = (kNodes(i,j) + kNodes(i,j-1)) / 2;
                kx2 = (kNodes(i,j) + kNodes(i,j+1)) / 2;
                ky1 = (kNodes(i,j) + kNodes(i-1,j)) / 2;
                ky2 = (kNodes(i,j) + kNodes(i+1,j)) / 2;

                %LHS
                LT(gkT,gkT-Ny) = -kx1 / (dx^2.0);                         % T1_dt
                LT(gkT,gkT-1) = -ky1 / (dy^2.0);                          % T2_dt
                LT(gkT,gkT) = (pcpNodes(i,j) / dt) + (kx1+kx2)/dx^2.0 + ... 
                            (ky1+ky2)/dy^2.0 ;                            % T3_dt  
                LT(gkT,gkT+1) = -ky2 / (dy^2.0);                          % T4_dt
                LT(gkT,gkT+Ny) = -kx2 / (dx^2.0);                         % T5_dt

                % RHS
                % RT(gkT) = tempNodes(i,j) / dt * pcpNodes(i,j);
                RT(gkT) = tempNodes(i,j) / dt * pcpNodes(i,j) - hRNodes(i,j) -  ...
                    hANodes(i,j) - hSNodes(i,j); % NEW - HW11
            end
        end
    end

    % Solve For New Temperature
    ST = LT\RT;

    % Reload into geometric arrays
    for j=1:1:Nx
        for i=1:1:Ny
            % Get global indexing
            gkT =  ((j-1)*Ny + i - 1) + 1;
            % Transfer dT,T values into geometric arrays
            % dT must be transfered first with how this is coded
            dTempNodes(i,j) = ST(gkT) - tempNodes(i,j);
            tempNodes(i,j) = ST(gkT);
        end
    end 

    % Apply dT constraints
    thisMaxDT = max(max(abs(dTempNodes)));
    if (r==1) && (thisMaxDT > dTMax)
        dt = dt / thisMaxDT * dTMax;
    else
        break
    end
end

% Interpolate new T from nodes to Markers
if (tIter == 1)
    % Use nodal T in interpolation for first iteration
    tempMarkers = interpNodesToMarkers2D(tempNodes,xm,ym,x,y,dx,dy,[1,Nx,1,Ny]);
else
    % Use nodal dT in interpolation for the following interations
    dTempMarkers = interpNodesToMarkers2D(dTempNodes,xm,ym,x,y,dx,dy,[1,Nx,1,Ny]);
    % Update marker T
    tempMarkers = tempMarkers + dTempMarkers;
end
%% Advect markers using 4th order Runga Kutta:
% 1.) Compute effective velocities, vxEff, vyEff.
vxEff = zeros(Nm,1);
vyEff = zeros(Ny,1);
for m=1:1:Nm
   % Pt. A
   xA = xm(m);
   yA = ym(m);
   vxA = vxMarkers(m);
   vyA = vyMarkers(m);
   % Pt. B
   xB = xA + vxA*(dt/2);
   yB = yA + vyA*(dt/2);
   vxB = interpNodesToMarkers2D(vxP,(xB),(yB),xP,yP,dx,dy,[1,Nx,1,Ny]);
   vyB = interpNodesToMarkers2D(vyP,(xB),(yB),xP,yP,dx,dy,[1,Nx,1,Ny]);
   % Pt. C
   xC = xA + vxB*(dt/2);
   yC = yA + vyB*(dt/2);
   vxC = interpNodesToMarkers2D(vxP,(xC),(yC),xP,yP,dx,dy,[1,Nx,1,Ny]);
   vyC = interpNodesToMarkers2D(vyP,(xC),(yC),xP,yP,dx,dy,[1,Nx,1,Ny]);
   % Pt. D
   xD = xA + vxC*(dt);
   yD = yA + vyC*(dt);
   vxD = interpNodesToMarkers2D(vxP,(xD),(yD),xP,yP,dx,dy,[1,Nx,1,Ny]);
   vyD = interpNodesToMarkers2D(vyP,(xD),(yD),xP,yP,dx,dy,[1,Nx,1,Ny]);
   % Effective Velocities
   vxEff(m) = (1/6) * (vxA + 2*vxB + 2*vxC + vxD);
   vyEff(m) = (1/6) * (vyA + 2*vyB + 2*vyC + vyD);
end

% 2.) Advect markers (preventing markers from leaving model domain)
for m=1:1:Nm 
    xNew = xm(m) + dt*vxEff(m);
    yNew = ym(m) + dt*vyEff(m);
    if (xNew > 0 && xNew < xsize) && (yNew > 0 && yNew < ysize)
        xm(m) = xNew;
        ym(m) = yNew;
    end
end

% Update time time
currentTime = currentTime + dt;

% Plot results For this timestep
figure(1);clf
% RHO
sp1 = subplot(3,3,1);
colormap('jet')
p1 = contourf(x,y,RHO,3);
colormap(jet)
if tIter == 1
    lims1 = caxis;
else
    caxis(lims1); 
end
colorbar();
axis ij image
hold on
quiver(x(1:5:Nx),y(1:5:Ny),vxNodes(1:5:Ny,1:5:Nx),vyNodes(1:5:Ny,1:5:Nx))
title('rho')
xlabel('x [m]')
ylabel('y [m]')

% VX
subplot(3,3,2)
p2 = pcolor(x,y,vxNodes);
set(p2,'edgecolor','none')
if tIter == 1
    lims2 = caxis;
else
    caxis(lims2); 
end
colorbar();
axis ij image
hold on
quiver(x(1:5:Nx),y(1:5:Ny),vxNodes(1:5:Ny,1:5:Nx),vyNodes(1:5:Ny,1:5:Nx))
title('vx')
xlabel('x [m]')
ylabel('y [m]')

% VY
subplot(3,3,3);
p3 = pcolor(x,y,vyNodes);
set(p3,'edgecolor','none')
if tIter == 1
    lims3 = caxis;
else
    caxis(lims3); 
end
colorbar();
axis ij image
hold on
quiver(x(1:5:Nx),y(1:5:Ny),vxNodes(1:5:Ny,1:5:Nx),vyNodes(1:5:Ny,1:5:Nx))
title('vy')
xlabel('x [m]')
ylabel('y [m]')

% Pressure
subplot(3,3,4);
p4 = pcolor(x,y,pressureNodes);
set(p4,'edgecolor','none')
if tIter == 1
    lims4 = caxis;
else
    caxis(lims4); 
end
colorbar();
axis ij image
hold on
quiver(x(1:5:Nx),y(1:5:Ny),vxNodes(1:5:Ny,1:5:Nx),vyNodes(1:5:Ny,1:5:Nx))
title('pressure')
xlabel('x [m]')
ylabel('y [m]')

% ETA
subplot(3,3,5)
p5 = contourf(x,y,ETA,3);
if tIter == 1
    lims5 = caxis;
else
    caxis(lims5); 
end
colorbar()
axis ij image
hold on
quiver(x(1:5:Nx),y(1:5:Ny),vxNodes(1:5:Ny,1:5:Nx),vyNodes(1:5:Ny,1:5:Nx))
title('eta')
xlabel('x [m]')
ylabel('y [m]')

% Temperature
subplot(3,3,6)
p6 = contourf(x,y,tempNodes);
if tIter == 1
    lims6 = caxis;
else
    caxis(lims6); 
end
colorbar()
axis ij image
title('T')
xlabel('x [m]')
ylabel('y [m]')

% Radioactive Heating
subplot(3,3,7)
p7 = pcolor(x,y,hRNodes);
shading interp;
if tIter == 1
    lims7 = caxis;
else
    caxis(lims7); 
end
colorbar()
axis ij image
title('H_r')
xlabel('x [m]')
ylabel('y [m]')

% Adiabatic Heating
subplot(3,3,8)
p8 = pcolor(x,y,hANodes);
shading interp;
if tIter == 1
    lims8 = caxis;
else
    caxis(lims8); 
end
colorbar()
axis ij image
title('H_a')
xlabel('x [m]')
ylabel('y [m]')

% Shear Heating
subplot(3,3,9)
p9 = pcolor(x,y,hSNodes);
shading interp;
if tIter == 1
    lims9 = caxis;
else
    caxis(lims9); 
end
colorbar()
axis ij image
title('H_s')
xlabel('x [m]')
ylabel('y [m]')

drawnow

fprintf('Current Iteration: %d \n',tIter)

end
% ------------------- END TIME LOOP ----------------------

catch ERR
    fprintf('!!! SIMULATION CRASHED AFTER %d iterations. See below error message for more details\n\n',tIter)
    rethrow(ERR)
end
% ------------------- END TRY-CATCH BLOCK-----------------



