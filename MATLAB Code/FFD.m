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
        
        % pressure data is (zMaxIndex-1)-by-(rMaxIndex-1)
        % z velocity data is (zMaxIndex)-by-(rMaxIndex-1)
        % r velocity data is (zMaxindex-1)-by-(rMaxIndex)
        % the difference is due to staggered grid
        
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

            % Updated by Jinghu on Apr 25 16:00
            obj.Urbar = ones(length(obj.rbar), length(obj.zbar));
            obj.Uzbar = zeros(length(obj.rbar), length(obj.zbar));
%             obj.Urbar = zeros(obj.zMaxIndex-1,obj.rMaxIndex);
%             obj.Uzbar = zeros(obj.zMaxIndex,obj.rMaxIndex-1);
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
            
            % Updated by Jinghu on Apr 25 16:00
            obj.Urbar = ones(length(obj.rbar), length(obj.zbar));
            obj.Uzbar = zeros(length(obj.rbar), length(obj.zbar));
%             obj.Urbar = zeros(obj.zMaxIndex-1,obj.rMaxIndex);
%             obj.Uzbar = zeros(obj.zMaxIndex,obj.rMaxIndex-1);
            obj.Ubar = [obj.Urbar(:); obj.Uzbar(:)];          % (KP-04/25)
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
        
        % Updated by Jinghu on Apr 25 16:00
        function grad_inside = InsideNodeGrd(obj,cur_zindex,cur_rindex)
        % calculate the gradient values of the pressure nodes
        % "pressure nodes" are the nodes at the finite volume center
        % Messed up; need to think more
        grad_inside = 0;
        end
        
        
        function iterateUrStar(obj) % (KP-04/25)
            % computes the intermediate non-divergence free r velocity
            % component by solving the decoupled r-momentum equation with
            % an implicit-explicit technique.
            n = length(obj.zbar); m = length(obj.rbar); nm = n*m;
            
            % compute diagonal elements of state matrix
            Omega1 = 1/obj.dtau + 2/(obj.Re*obj.drbar^2) ...
                     + 1./(obj.Re*repmat(obj.rbar', n, 1)) ...
                     + 2/(obj.Re*obj.dzbar^2);
            Omega2 = -1/(2*obj.Re*repmat(obj.rbar', n, 1)*obj.drbar) ...
                    - 1/(obj.Re*obj.drbar^2);
            Omega3 = -1/(obj.Re*obj.dzbar^2);
            Omega4 = 1/(2*obj.Re*repmat(obj.rbar', n, 1)*obj.drbar) ...
                    - 1/(obj.Re*obj.drbar^2);
            Omega5 = -1/(obj.Re*obj.dzbar^2);
            
            % compile sparse diagonal state matrix
            obj.ArStar = sparse(diag(Omega1) ...
                   + diag(Omega2(1:end-1), 1) ...
                   + diag(Omega3*ones(nm-m, 1), m) ...
                   + diag(Omega4(2:end), -1) ...
                   + diag(Omega5*ones(nm-m, 1), -m));               
               
            % compute intermediate r-momentum advection operator
            Nrhs = computeNrhs(obj);
            
            % compute intermediate r-velocity
            obj.Ustar(1:nm) = obj.ArStar\Nrhs;
                        
        end
        
        function N = computeNrhs(obj) % (KP-04/25)
            % computes advection operator for intermediate
            % non-divergence-free r-momentum equation
            n = length(obj.zbar); m = length(obj.rbar); nm = n*m;
            
            % compute diagonal elements of state matrices
            Lambda1 = -1/(2*obj.drbar); Lambda2 = 1/(2*obj.drbar);
            Pi1 = -1/(2*obj.dzbar); Pi2 = 1/(2*obj.dzbar); 
            
            % compile sparse diagonal state matrices
            obj.ArN = sparse(diag(Lambda1*ones(nm-1, 1), 1) ...
                    + diag(Lambda2*ones(nm-1, 1), -1));
               
            obj.BrN = sparse(diag(Pi1*ones(nm-m, 1), m) ...
                    + diag(Pi2*ones(nm-m, 1), -m));
            
            % compute advection operator
            ur = obj.Ubar(1:nm);
            uz = obj.Ubar(nm+1:end);
            N = diag(ur)*(ones(nm, 1)./obj.dtau - obj.ArN*ur) ...
                + diag(uz)*obj.BrN*ur;                                 
        end
       
    end
             
end
