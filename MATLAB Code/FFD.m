classdef FFD < handle
    properties
        % list all variables that will be passed between methods in the
        % model here. think of them as global variables within FFD. if a
        % value is included in the definition then it is treated as the
        % default value.
                
        % dependent variables
        tauEnd = 1      % simulation end time
        rbar            % r-coordinates for calculations
        zbar            % z-coordinates for calculations
        tau             % non-dimensional time vector
        tauNow = 0      % current iteration non-dimensional time
        
        % Updated by Jinghu on Apr 25 16:00
        rMaxIndex       % node num in r dir including two end points
        zMaxIndex       % node num in z dir including two end points
        Urbar           % r direction velocity storage 
        Uzbar           % z direction velocity storage 
        Ubar            % full velocity vector (KP-04/25)
        Ustar           % intermediate velocity vector (KP-04/25)
        Pbar            % pressure storage 
        
        testAr
        
        % Staggered grid real node number:
        % pressure data is (zMaxIndex-1)-by-(rMaxIndex-1)
        % z velocity data is (zMaxIndex)-by-(rMaxIndex-1)
        % r velocity data is (zMaxindex-1)-by-(rMaxIndex)
        
        %For ghost point
        % vz has 2 columns of (M-2) ghost points
        % vr has 2 rows    of (N-2) ghost points
        % P  has 2 columns of (M-1) and 2 rows of (N-1) ghost points
        
        
        % Updated by Jinghu on Apr 25 16:00
        rbarOutlet = 0.3                % outlet radius
        UzInlet = 1.0                   % inlet z dir velocity
        
        % static model parameters with default values
        H = 1                           % (m) initial height of particles
        ztop                            % current nondimensional height of 
                                        % particle region
        a0 = 0.04                       % non-dimensional radius of 
                                        % opening at bottom
        b = 0.4                         % outer nondimensional radius
        g = 9.80665                     % (m/s2)
        rhoPack = 2000                  % (kg/m3) particle packed bulk 
                                        % density
        rhoLoose = 1810                 % (kg/m3) particle loose bulk 
                                        % density
        mu = 1 %2.5*1.81e-5                % (kg/ms) viscosity of particles 
                                        % moving in air(Bicerano, Douglas 
                                        % and Brune, 1999)
        nu                              % (m2/s) dynamic viscosity of 
                                        % particles moving in air  
        Uinf = 0.01                     % (m/s) normalization velocity
        
        % state matrices for discretized equations
        ArStar              % intermediate r-momentum state matrix (KP-04/25)
        AzStar              % intermediate z-momentum state matrix (KP-04/25)
        ArN                 % intermediate r-momentum advection operator 
                            % state matrix (KP-04-25)                            
        BrN                 % "" ""
        AzN                 % intermediate z-momentum advection operator 
                            % state matrix (KP-04-25)
        BzN                 % "" ""
        Nr                  % intermediate r-momentum advection operator
        Nz                  % intermediate z-momentum advection operator
        
        
        % nondimensional terms
        Re                  % Reynolds number w.r.t. Uinf
        Fr                  % Freud number w.r.t Uinf
        
        % discretization parameters
        dtau = 0.1     % nondimensional time-step
        dt             % (s) time-step
        drbar = 0.1    % nondimensional radius large mesh size
        dzbar = 0.1    % nondimensional height large mesh size 
        drhat = 0.01   % nondimensional radius small mesh size
        dzhat = 0.01   % nondimensional height small mesh size
        
        % figures and tables
        vtbl            % table containing all static variables
        ubarfig         % figure showing velocity in top boundary
        wbarfig         % figure showing velocity in center boundary
        computeBUfig    % figure showing approximated bulk velocities
    end      
    methods 
        % methods are defined to operate on class properties. the entire
        % framework for the numerical simulation will be constructed here.
        function obj = FFD()
            % include any initializations for properties here. It will be
            % ran whenever a new class instantiation is performed.
            obj.ztop = obj.H;
            obj.rbar = 1e-6:obj.drbar:obj.b;          
            obj.zbar = 0:obj.dzbar:1;
            obj.rMaxIndex = size(obj.rbar,2);
            obj.zMaxIndex = size(obj.zbar,2);

            % Staggered grid real node number:
            % pressure data is (zMaxIndex-1)-by-(rMaxIndex-1)
            % z velocity data is (zMaxIndex)-by-(rMaxIndex-1)
            % r velocity data is (zMaxindex-1)-by-(rMaxIndex)          
            
            
            % Updated 05-25 JH staggered grid mod
            %obj.Urbar = ones(length(obj.rbar), length(obj.zbar));
            %obj.Uzbar = zeros(length(obj.rbar), length(obj.zbar));
            obj.Urbar = zeros(obj.zMaxIndex-1,obj.rMaxIndex);
            obj.Uzbar = zeros(obj.zMaxIndex,obj.rMaxIndex-1);
            obj.Ubar = [obj.Urbar(:); obj.Uzbar(:)];          % (KP-04/25)
            obj.Pbar  = zeros(obj.zMaxIndex-1,obj.rMaxIndex-1);
            
            
            obj.tau = 0:obj.dtau:obj.tauEnd;         
            obj.dt = obj.dtau*obj.H/obj.Uinf;
            obj.nu = obj.mu/obj.rhoLoose;  
            obj.Re = obj.H*obj.Uinf/obj.nu;              
            obj.Fr = obj.Uinf/sqrt(obj.g*obj.H);
        end
        function reInitObj(obj)
            % recomputes static variables. should be ran if any of the
            % system properties are changed externally to ensure that 
            % everything is consistent.
            obj.ztop = obj.H;
            obj.rbar = 1e-6:obj.drbar:obj.b;          
            obj.zbar = 0:obj.dzbar:1;
            
            obj.rMaxIndex = size(obj.rbar,2);
            obj.zMaxIndex = size(obj.zbar,2);
            
            
            %obj.Urbar = ones(length(obj.rbar), length(obj.zbar));
            %obj.Uzbar = zeros(length(obj.rbar), length(obj.zbar));
            obj.Urbar = zeros(obj.zMaxIndex-1,obj.rMaxIndex);
            obj.Uzbar = zeros(obj.zMaxIndex,obj.rMaxIndex-1);
            %obj.Ubar = [obj.Urbar(:); obj.Uzbar(:)];          % (KP-04/25)
            obj.Pbar  = zeros(obj.zMaxIndex-1,obj.rMaxIndex-1);
            
            obj.tau = 0:obj.dtau:obj.tauEnd;         
            obj.dt = obj.dtau*obj.H/obj.Uinf;
            obj.nu = obj.mu/obj.rhoLoose;  
            obj.Re = obj.H*obj.Uinf/obj.nu;              
            obj.Fr = obj.Uinf/sqrt(obj.g*obj.H);
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add any additional functions here
        % as a general rule, try to keep the number of operations for 
        % each method at a minimum. if a function is more than ~10 lines
        % then it should probably be split into multiple methods.
        % in general, the methods can be treated as a typical matlab
        % function with the acception of the method for passing in
        % properties. the general form for functions that have an output is
        function x = exampleFunction(obj, input1, input2)
            x = input1*obj.Re + input2*obj.Fr;
        end
        % and the general form for operations is
        function exampleOperation(obj, input1, input2)
            % updates meshed domain
            obj.rbar = 1e-6:obj.drbar:input1;
            obj.zbar = 0:obj.dzbar:input2;
            reInitObj(obj);
        end
        % where the inputs are just arbitrary variables that aren't defined
        % as properties. it is typically a good practice to limit the usage
        % of these and instead just operate with class properties. if a
        % non-property variable is being passed into multiple functions 
        % then it should instead be defined as a property.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function computeUStar(obj)
            % computes the intermediate non-divergence free r velocity
            % component by solving the decoupled r-momentum equation with
            % an implicit-explicit technique.
            n = length(obj.zbar); m = length(obj.rbar); nm = n*m;
            
            % compute intermediate state matrices
            computeArStar(obj);
            computeAzStar(obj);
                                        
            % compute intermediate advection operators
            computeNr(obj);
            computeNz(obj);
            
            % compute intermediate r-velocity
            obj.Ustar = [obj.ArStar, zeros(nm); zeros(nm), obj.AzStar] ...
                        \[obj.Nr; obj.Nz];
                    
            % set boundary conditions
            applyUrBoundaries(obj);
            applyUzBoundaries(obj);                       
        end
        function computeArStar(obj)
            % evaluates state matrix for intermediate r-velocity comp.
            
            % n the number of z nodes; m the number of r nodes
            n = length(obj.zbar); m = length(obj.rbar); 
            
            % the number of real points m*(n-1)
            n_1m = (n-1)*m;
            
            % two rows of ghost points added
            % the number of total points m * (n-1) + 2m = m * (n+1)
            % the number of ghost points 2m
            TotalNum = m*(n+1);
            
            % obj.a0 the outlet radius
            %find out how many nodes belong to outlet
            OutletNodeNum = sum((obj.rbar)<=obj.a0);
            
            % Omega1*u_ij + Omega2*u_i(j+1) + Omega3*u_(i+1)j
            %  +Omega4*u_i(j-1) + Omega5*u_(i-1)j = N
          
            % self coefficient Omega1; have to be 1
            Omega1_GhostPoint = ones(m,1);
            % real internal points Omega1
            Omega1_main = -1/obj.dtau - 2/(obj.Re*obj.drbar^2) ...
                     - 1./(obj.Re*repmat(obj.rbar', n-1, 1).^2) ...
                     - 2/(obj.Re*obj.dzbar^2);  
            
            % left boundary condition: Ur_ij = 0 (symmetry) Omega1=1.0
            Omega1_main(1:m:n_1m)=1.0;
            % right boundary condition Ur_ij = 0 (wall) Omega1 = 1.0
            Omega1_main(m:m:n_1m)=1.0;
            
            % Omega2 for ghost points: all should be zero
            % the 1st upper subdiagonal; first element popped off
            % first element will be cut-off
            Omega2_UpperGhostPoint = zeros(m,1);
            % last elment of lower ghost points won't have Omega2
            Omega2_LowerGhostPoint = zeros(m-1,1);
            
            % real internal points Omega2
            Omega2_main = 1/(2*obj.Re*repmat(obj.rbar',n-1,1)*obj.drbar)...
                    + 1/(obj.Re*obj.drbar^2);
            Omega2_main = Omega2_main';%Omega2 now ((n-1)xm)-by-1
            
            % left boundary condition: Ur_ij = 0 (symmetry) Omega2 = 0
            Omega2_main(1:m:n_1m)=0.0;
            % right boundary condition Ur_ij = 0 (wall) Omega2 = 0
            Omega2_main(m:m:n_1m)=0.0;
          
            % upper boundary condition Ur_ij + Ur_(i+1)j = 0; Omega3 = 1.0
            Omega3_UpperGhostPoint = ones(m,1);
            
            % two upper corner ghost points won't have Omega3;
            % U_ij=0
            Omega3_UpperGhostPoint([1,m])=0.0;
            
            % real internal points Omega3
            Omega3_main = (1/(obj.Re*obj.dzbar^2))*ones(n_1m, 1);
            
            % left and right boundary points: Omega3 = 0
            Omega3_main(1:m:n_1m) = 0.0;
            Omega3_main(m:m:n_1m) = 0.0;
            
            % Omega4 should be zero for all ghost points
            % first upper ghost points doesn't have Omega4
            Omega4_UpperGhostPoint = zeros(m-1,1);
            Omega4_LowerGhostPoint = zeros(m,1);
           
            % real internal points Omega4
            Omega4_main = -1/(2*obj.Re*repmat(obj.rbar',n-1,1)*obj.drbar) ...
                    + 1/(obj.Re*obj.drbar^2);
            Omega4_main = Omega4_main';
            % left and right boundary points: Omega4 = 0
            Omega4_main(1:m:n_1m)=0.0;
            Omega4_main(m:m:n_1m) = 0.0;
       
            % lower Boundary conditions: Ur_ij + Ur_(i-1)j = 0 Omega5 = 1.0
            Omega5_LowerGhostPoint = ones(m,1);
            
            % for outlet ghost points: Ur_ij - Ur_(i-1)j = 0 (no gradient)
            % Omega5 = -1.0
            Omega5_LowerGhostPoint([1:OutletNodeNum])=-1.0;
            % two lower corner ghost points don't have meaningful eqn
            Omega5_LowerGhostPoint([1,m])=0.0;
            
            % real internal points Omega5
            Omega5_main = (1/(obj.Re*obj.dzbar^2))*ones(n_1m, 1);
            
            % left and right boundary points: Omega5 = 0
            Omega5_main(1:m:n_1m)=0.0;
            Omega5_main(m:m:n_1m)=0.0;
           
            % total Omega vectors to be put in matrix as sub-diagnals
            
            TOmega1 = [Omega1_GhostPoint;Omega1_main;Omega1_GhostPoint];
            
            % first 0.0 is to be cut off by spdiags function
            TOmega2 = [0.0;Omega2_UpperGhostPoint;Omega2_main;...
                       Omega2_LowerGhostPoint];
            
            % first mx1 vectors are to be cut off by spdiags function
            TOmega3 = [zeros(m,1);Omega3_UpperGhostPoint;Omega3_main];
            
            % last 0.0 is to be cut off by spdiags function
            TOmega4 = [Omega4_UpperGhostPoint;Omega4_main;...
                       Omega4_LowerGhostPoint;0.0];
            
            % last mx1 vectors are to be cut off by spdiags function
            TOmega5 = [Omega5_main;Omega5_LowerGhostPoint;zeros(m,1)];
            
            % build up Ar*
            obj.ArStar = spdiags([TOmega1, TOmega2, TOmega3,...
                         TOmega4, TOmega5], ...
                         [0, 1, m, -1, -m], TotalNum, TotalNum);     
            
              
            
            
            % compute diagonal elements of state matrix
            %Omega1 = -1/obj.dtau - 2/(obj.Re*obj.drbar^2) ...
            %         - 1./(obj.Re*repmat(obj.rbar', n, 1).^2) ...
            %         - 2/(obj.Re*obj.dzbar^2);
            %Omega2 = 1/(2*obj.Re*repmat(obj.rbar', n, 1)*obj.drbar) ...
            %        + 1/(obj.Re*obj.drbar^2);
            %Omega3 = 1/(obj.Re*obj.dzbar^2);
            %Omega4 = -1/(2*obj.Re*repmat(obj.rbar', n, 1)*obj.drbar) ...
            %        + 1/(obj.Re*obj.drbar^2);
            %Omega5 = 1/(obj.Re*obj.dzbar^2);
            %
            % compile sparse diagonal state matrix                      
            %obj.testAr = spdiags([Omega1, Omega2', Omega3*ones(nm, 1), ...
            %             Omega4', Omega5*ones(nm, 1)], ...
            %             [0, 1, m, -1, -m], nm, nm);  
        end
        
        function computeNr(obj)
            % computes advection operator for intermediate
            % non-divergence-free r-momentum equation
            n = length(obj.zbar); m = length(obj.rbar); nm = n*m;
            
            OutletNodeNum = sum((obj.rbar)<=obj.a0);
            
            TotalNum = m*(n+1);
            n_1m = m*(n-1);
            
            
            
            
            
            % compute diagonal elements of state matrices
            Lambda1 = 1/(2*obj.drbar); Lambda2 = -1/(2*obj.drbar);
            Pi1 = 1/(2*obj.dzbar); Pi2 = -1/(2*obj.dzbar); 
            
            % compile sparse diagonal state matrices
            obj.ArN = spdiags([Lambda1*ones(nm, 1), ...
                               Lambda2*ones(nm, 1)], [1, -1], nm, nm);              
            obj.BrN = spdiags([Pi1*ones(nm, 1), Pi2*ones(nm, 1)], ...
                              [m, -m], nm, nm);             
            Z = spdiags(0.25*ones(nm, 4), [0, m, m-1, -1], nm, nm);

            % compute advection operator
            ur = obj.Ubar(1:nm);
            uz = obj.Ubar(nm+1:end);
            obj.Nr = diag(ur)*(obj.ArN*ur - ones(nm, 1)./obj.dtau) ...
                + diag(Z*uz)*obj.BrN*ur;                                 
        end
        
        function computeAzStar(obj)
            % evaluates state matrix for intermediate r-velocity comp.
            
            % n the number of z nodes; m the number of r nodes
            % the number of real points (m-1)*n
            n = length(obj.zbar); m = length(obj.rbar); m_1n = n*(m-1);
            
            % two cols of ghost points added
            % the number of total points (m-1) * n + 2n = (m+1) * n
            % the number of ghost points 2n
            
            OutletNodeNum = sum((obj.rbar)<=obj.a0);
            
            %total points = real points (real internal + real boundary)
            %               + ghost points
            TotalNum = n*(m+1);
            
            %       Omega1*u_ij + Omega2*u_i(j+1) + Omega3*u_(i+1)j
            %   + Omega4*u_i(j-1) + Omega5*u_(i-1)j = B
            
            % internal real points Omega1
            Omega1_main = (-1/obj.dtau - 2/(obj.Re*obj.drbar^2) ...
                     - 2/(obj.Re*obj.dzbar^2))*ones(TotalNum,1);
            
            % left boundary Omega1
            Omega1_main(1:m+1:TotalNum)=1.0;
            % right boundary Omega1
            Omega1_main(m+1:m+1:TotalNum)=1.0;
            % upper boundary Omega1
            Omega1_main(1:m+1)=1.0;
            % lower boundary Omega1
            Omega1_main((m+1)*(n-1)+1:end)=1.0;
            
            % internal real points Omega2
            Omega2_main = 1/(2*obj.Re*(repmat(([-obj.drbar,...
                     obj.rbar])',n, 1) ...
                     + obj.drbar/2)*obj.drbar) + 1/(obj.Re*obj.drbar^2);
            
            Omega2_main=Omega2_main';
            % left boundary except for two corner ghost points
            % left boundary condition: zero gradient 
            % Uz,ij = Uz,i(j+1)
            % Omega2 = -1.0
            Omega2_main(m+2:m+1:(m+1)*(n-1))=-1.0;
            % Omega2 for two left corner ghost points = 0
            Omega2_main(1)=0.0;
            Omega2_main(1:m+1)=0.0;
            
            % Lower boundary condition Omega2 = 0
            Omega2_main([(m+1)*(n-1)+1:end])=0.0;
            
            % Right boundary condition Omega2 = 0
            Omega2_main(m+1:m+1:TotalNum)=0.0;
            
            % Last element doesn't have Omega2
            Omega2_main(end)=[];
            % first zero is to be cut off
            Omega2_main = [0.0;Omega2_main];
            
            
            
            % real internal points Omega3
            Omega3_main = (1/(obj.Re*obj.dzbar^2))*ones(TotalNum,1);
            % upper boundary: Uz,ij = 1
            Omega3_main(1:m+1)=0.0;
            % left boundary Omega3 = 0
            Omega3_main(1:m+1:TotalNum)=0.0;
            
            % right boundary Omega3 = 0.0
            Omega3_main(m+1:m+1:TotalNum)=0.0;
            
            % last row of nodes (lower boundary) don't have Omega3
            Omega3_main((m+1)*(n-1)+1:end)=[];
            % zeros are to be cut off
            Omega3_main = [zeros(m+1,1);Omega3_main];
            
            % real internal points Omega4
            Omega4_main = -1/(2*obj.Re*(repmat(([-obj.drbar,...
                     obj.rbar])',n, 1) ...
                     + obj.drbar/2)*obj.drbar) + 1/(obj.Re*obj.drbar^2);
            Omega4_main=Omega4_main';
            % left boundary Omega4 = 0
            Omega4_main(1:m+1:end)=0.0;
            % lower boundary Omega4 = 0
            Omega4_main([(m+1)*(n-1)+1:end])=0.0;
            
            % right boundary Omega4 = 0.0
            Omega4_main(m+1:m+1:(m+1)*(n-1))=1.0;
            % upper boundary Omega4 = 0.0
            Omega4_main(1:m+1)=0.0;
            
            % first element doesn't have Omega4
            Omega4_main(1)=[];
            %zeros to be cut off
            Omega4_main = [Omega4_main;0.0];
            
            % real internal points Omega5
            Omega5_main = 1/(obj.Re*obj.dzbar^2)*ones(TotalNum,1);
            % left boundary Omega5 = 0.0
            Omega5_main(1:m+1:end)=0.0;
            % right boundary
            Omega5_main(m+1:m+1:end)=0.0;
            % lower boundary outlet condition: zero gradient
            % U_ij = U_(i-1)j  Omega5 = -1.0
            Omega5_main((m+1)*(n-1)+1:(m+1)*(n-1)+OutletNodeNum)=-1.0;
            % lower boundary outlet condition: wall
            % U_ij = 0; Omega5 = 0.0
            Omega5_main((m+1)*(n-1)+OutletNodeNum:end)=0.0;
            % first (m+1) elements don't have Omega5
            Omega5_main(1:m+1)=[];
            Omega5_main = [Omega5_main;zeros(m+1,1)];
            obj.AzStar = spdiags([Omega1_main,Omega2_main,...
                         Omega3_main, Omega4_main, ...
                         Omega5_main],[0, 1, m+1, -1, -m-1],TotalNum, TotalNum); 
            
            
            % compute diagonal elements of state matrix
            Omega1 = -1/obj.dtau - 2/(obj.Re*obj.drbar^2) ...
                     - 2/(obj.Re*obj.dzbar^2);
            Omega2 = 1/(2*obj.Re*(repmat(obj.rbar', n, 1) ...
                     + obj.drbar/2)*obj.drbar) + 1/(obj.Re*obj.drbar^2);
            Omega3 = 1/(obj.Re*obj.dzbar^2);
            Omega4 = -1/(2*obj.Re*(repmat(obj.rbar', n, 1) ...
                     + obj.drbar/2)*obj.drbar) + 1/(obj.Re*obj.drbar^2);
            Omega5 = 1/(obj.Re*obj.dzbar^2);
            
            % compile sparse diagonal state matrix
            %obj.AzStar = spdiags([Omega1*ones(nm, 1), ...
            %             Omega2', Omega3*ones(nm, 1), Omega4', ...
            %             Omega5*ones(nm, 1)], [0, 1, m, -1, -m], nm, nm);                      
                                   
        end
        
        function computeNz(obj)
            % computes advection operator for intermediate
            % non-divergence-free r-momentum equation
            n = length(obj.zbar); m = length(obj.rbar); nm = n*m;
            
            % compute diagonal elements of state matrices
            Lambda1 = 1/(2*obj.drbar); Lambda2 = -1/(2*obj.drbar);
            Pi1 = 1/(2*obj.dzbar); Pi2 = -1/(2*obj.dzbar); 
            
            % compile sparse diagonal state matrices
            obj.AzN = spdiags([Lambda1*ones(nm, 1), ...
                               Lambda2*ones(nm, 1)], [1, -1], nm, nm);              
            obj.BzN = spdiags([Pi1*ones(nm, 1), Pi2*ones(nm, 1)], ...
                              [m, -m], nm, nm);                           
            R = spdiags(0.25*ones(nm, 4), [0, -m, -m+1, 1], nm, nm);

            % compute advection operator
            ur = obj.Ubar(1:nm);
            uz = obj.Ubar(nm+1:end);
            obj.Nz = diag(R*ur)*(obj.ArN*ur - ones(nm, 1)./obj.dtau) ...
                   + diag(uz)*(obj.BrN*ur - ones(nm, 1)./obj.dtau) ...
                   - ones(nm, 1)./obj.Fr^2;                                 
        end
        
        function applyUrBoundaries(obj)
            % sets boundary conditions for r-velocity at intermediate
            % computation step.
            n = length(obj.zbar); m = length(obj.rbar); nm = n*m;
            ur = obj.Ustar(1:nm);

            % boundary 1 at z = 1 (free surface)
            ur(1:m-1) = 0;
            
            % boundary 2 at r = 0 (centerline)
            ur(1:m:end-m) = 0;
            
            % boundary 3 at z = 0 (bin opening)
            [~, ia] = min(abs(obj.a0 - obj.rbar));
            ur(end-m:end-m+ia) = ur(end-2*m:end-2*m+ia);
            
            % boundary 4 at z = 0 (bottom wall)
            ur(end-m+ia+1:end) = 0;
            
            % boundary 5 at r = b (outer radius wall)
            ur(m:m:end) = 0;
            
            % reasign velocity
            obj.Ustar(1:nm) = ur;                      
        end
        
        function applyUzBoundaries(obj)
            % sets boundary conditions for r-velocity at intermediate
            % computation step.
            n = length(obj.zbar); m = length(obj.rbar); nm = n*m;
            uz = obj.Ustar(nm+1:end);
            
            % boundary 1 at z = 1 (free surface)
            uz(1:m-1) = obj.Uinf;
            
            % boundary 2 at r = 0 (centerline)
            uz(1:m:end-m) = uz(2:m:end-m+1);
            
            % boundary 3 at z = 0 (bin opening)
            [~, ia] = min(abs(obj.a0 - obj.rbar));
            uz(end-m:end-m+ia) = uz(end-2*m:end-2*m+ia);
            
            % boundary 4 at z = 0 (bottom wall)
            uz(end-m+ia+1:end) = 0;
            
            % boundary 5 at r = b (outer radius wall)
            uz(m:m:end) = 0;

            % reasign velocity
            obj.Ustar(nm+1:end) = uz;                      
        end
       
    end
             
end
