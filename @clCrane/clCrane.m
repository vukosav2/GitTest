% All physical properties are in standard metric units (kg, N, m, s, etc.)

classdef clCrane < handle
    
    properties (Constant, GetAccess=private)
        szMeshFileRootName = 'Mesh';
        szPlotFolder = 'PlotDump/';
        nVersion = 7;
    end
    
    properties (GetAccess=public, SetAccess=public)
        
        nVerLoaded                      % File loaded version
        
        triOpt = struct(...
            'method', 1, ...
            'meshIdx', 11, ...          % Only applies to method 1
            'S', [], ...                % 3D matrix storing simplex vertices
            'V', [], ...                % XY coordinates of all vertices
            'tri', [], ...              % Connectivity of vertices for simplices
            'tr', []);                  % MATLAB barycentric rep. of V and tri
                     
        simPar = struct(...
            'end_time', 5, ...          % Simulation stopping conditions
            'end_distance', 10.8, ...   % Simulation stopping conditions
            'tSpanCheck', 0.01, ...     % Time spanned between each simplex check
            'methodODE', 5, ...         % See note below
            'tResolution', 1/500);      % See note below
        
        % methodODE:
        % Method 1 updates the controller much more often, which may be
        % unrealistic, and note that ode45 interpolates, which messes with the
        % embedded control logic
        % 1 uses fnLcranemodel, 2 uses fnLcranemodel2, 3 uses the regulator
        % simplex approach, 4 uses fnNcranemodel with adaptable pole
        % placement (using different linearized systems - caution, this is
        % a problem with this), and 5 uses fnNcranemodel with just pole
        % placement (system linearized about zero)
 
        % tResolution:
        % For ODE solver, gives the time steps for solved solution sampling
        % points returned by the solver (internally it still
        % uses adaptive time intervals though)
        
        % Note for pLambda: it indirectly sets the pendulum length
        % Select pLambda between [0,1] to linearly scale the pendulum length
        % 0 => length is equal to the (room height - obstacle height), meaning the
        % obstacle is trivial to clear
        % 1 => the length is equal to the room height, the hardest case
        param = struct(...
            'hr', 0.33, ...             % height of room
            'wr', 2.2, ...              % width of room
            'ho', 0.099, ...            % height of obstacle
            'wo', 0.3, ...              % width of obstacle
            'M', 1.0731, ...            % mass of cart
            'm', 0.23, ...              % mass of pendulum bob
            'g', 9.8, ...               % gravitational acceleration
            'pLambda', 1, ...           % (Sets L, see note above)
            'L', [], ...                % length of pendulum
            'A', [], ...                % Linearized A matrix about zeros
            'B', [], ...                % Linearized B matrix about zeros
            'a', []);                   % Affine part: zero for this case
              
        contr = struct(...
            'U', zeros(12,3), ...           % Initial Forces on Triangulation
            'P4D', [-3,-4,-5,-6], ...       % Poles
            'poleContr', 5, ...             % See note below
            'val_SwitchToPoleContr', 9);  % Used with poleContr

        % Pole placement type of controller and conditions for use
        %  4 for 4D pole placement based on z threshold ("val_SwitchToPoleContr")
        %  5 for 4D pole placement based on reaching simplex
        %   "val_SwitchToPoleContr"      
        
        % Simulation processing variables
        g_polesActive
        g_nCountToActivate
        
        % Simulation results
        soln = struct(...
            'x0', [-0.35; 0; 0; 0], ...     	% Initial condition
            'xf', [0.7; 0; 0; 0], ...          	% Desired equilibrium location
            't', [], ...
            'x', [], ...
            'u', [], ...
            'sIdx', [], ...
            'pOn', [], ...
            'bExitedPolytope', 0, ...           % true if exited/failure
            'testID', []);                      % string containing testID       
        
        animPar = struct(...
            'nFramesPerSecond', 60, ...
            'fTimeFlowFactor', 1, ...       % < 1 => slower; > 1 => faster
            'graphLineStyle', '-', ...      % The usual MATLAB styles, like 
            ...                                 '.', '-', '--', etc.
            'trackLineStyle', '--', ...
            'bDrawTrackCtrlChangePts', true, ...
            'drawSchPos', [], ...           % [-0.5, pi/6]; % z and theta position
            'movieSaveName', []);           % extension and FPS are set in name later 
            
    end
    
    methods
        
        % Constructor
        function obj = clCrane(loadName)
            obj.nVerLoaded = obj.nVersion;
            fnParseLogFile(obj, loadName);
            fnInitializeScenario(obj);
        end
        
%         fnInitializeScenario(obj)
%         [A, B, a] = fnLinearizeCraneSys(obj, xLin, uLin)
%         fnParseLogFile(obj, fName)
%         fnTriangulation(obj)
%         fnDefinePolytopeVertices(obj)
%         simplexIdx = fnCheckSimplex(obj, x)
%         fnCraneControl(obj)
%         [dx, u, sIdx, pOn] = fnLcranemodel(~, x, obj, K, g, K4D, g4D)
%         [dx] = fnLcranemodel2(~, x, obj, K, g)
%         [dx, u, sIdx, pOn] = fnNcranemodel(~, x, obj, K, g, K4D, g4D)
%         [K_, g_, sIdx, pOn] = fnControlLogic(obj, x, K, g, K4D, g4D)
%         fnCraneAnimation(obj)
%         fnCreateLogFileOfTest(obj, resultsFolderName)
        
        function fnInitializeScenario(obj)
            
            % Compute pendulum length
            obj.param.L = obj.param.hr - (1 - obj.param.pLambda) * obj.param.ho;
            
            % Store the system matrices, linearized about the most obvious eqb point
            [obj.param.A, obj.param.B, obj.param.a] = fnLinearizeCraneSys(obj, ...
                [0; 0; 0; 0], 0);
            
            % Create the triangulation
            fnTriangulation(obj);
            
            % Turn warning off and on, as it will complain, in general, that not all
            % the points in V are part of the triangulation, tri
            warning('off','all');
            obj.triOpt.tr = TriRep(obj.triOpt.tri, ...
                obj.triOpt.V(:, 1), obj.triOpt.V(:, 2));
            warning('on','all');
            
            % Verify if inside polytope
            if fnCheckSimplex(obj, obj.soln.x0) == 0
                error('The state initial vector falls outside the polytope');
            end
            
            if fnCheckSimplex(obj, obj.soln.xf) == 0
                error('The state equilibrium vector falls outside the polytope');
            end
            
            % Post process the control values at the simplex vertices
            if size(obj.contr.U, 1) ~= size(obj.triOpt.S, 3)
                nS = size(obj.triOpt.S, 3);
                if size(obj.contr.U, 1) > nS
                    obj.contr.U(nS+1:end, :) = []; % delete these
                else
                    obj.contr.U = [obj.contr.U; ...
                        zeros(nS - size(obj.contr.U, 1), 3)]; % pad with extra zeros
                end
            end
        end
        
        function [A, B, a] = fnLinearizeCraneSys(obj, xLin, uLin)
            
            % Unpack parameters for convenience
            M = obj.param.M;
            m = obj.param.m;
            L = obj.param.L;
            g = obj.param.g;
            
            % xLin is the state and uLin is the input about which to linearize
            % Note that xLin(1) and xLin(2) should usually be zero because
            % these terms do not show up anywhere in the linearized system matrices
            
            % Define more terms for convenience
            t1 = M + m*sin(xLin(3))^2;
            
            % A matrix
            a23 = ((m*g*cos(2*xLin(3)) + m*L*cos(xLin(3))*xLin(4)^2) * (t1) - ...
                (0.5*m*g*sin(2*xLin(3)) + m*L*sin(xLin(3))*xLin(4)^2 + uLin) * ...
                (m*sin(2*xLin(3)))) / (t1)^2;
            a24 = (2 * xLin(4) * m * L * sin(xLin(3))) / t1;
            a43 = (-((M+m)*g*cos(xLin(3)) + m*L*xLin(4)^2*cos(2*xLin(3)) - ...
                uLin*sin(xLin(3))) * (L * t1) + (m * L * sin(xLin(3))) * ...
                ((M+m)*g*sin(xLin(3)) + 0.5*m*L*xLin(4)^2*sin(2*xLin(3)) + ...
                uLin*cos(xLin(3)))) / (L * t1)^2;
            a44 = -(m * L * xLin(4) * sin(2 * xLin(3))) / (L * t1);
            A = [0 1 0 0;
                0 0 a23 a24;
                0 0 0 1;
                0 0 a43 a44];
            
            % B matrix
            b2 = 1 / t1;
            b4 = -cos(xLin(3)) / (L * t1);
            B = [0; b2; 0; b4];
            
            % a matrix (affine part)
            % This part vanishes if the linearization point is an equilibrium point
            a1 = xLin(2);
            a2 = (0.5*m*g*sin(2*xLin(3)) + m*L*sin(xLin(3))*xLin(4)^2 + uLin) / t1;
            a3 = xLin(4);
            a4 = -((M+m)*g*sin(xLin(3)) + 0.5*m*L*xLin(4)^2*sin(2*xLin(3)) ...
                + uLin*cos(xLin(3))) / (L * t1);
            a = [a1; a2; a3; a4];
            
            obj.param.A = A;
            obj.param.B = B;
            obj.param.a = a;
        end
        
        function fnTriangulation(obj)
            
            %%  Obtain the hard coded vertex coordinates (scaled by problem parameters)
            
            fnDefinePolytopeVertices(obj);
            
            %% Three options for the triangulation mesh
            
            triangFolder = strrep(pwd, '\', '/');
            if (triangFolder(end) ~= '/')
                triangFolder = [triangFolder '/'];
            end
            triangFolder = [triangFolder 'TriangulationMeshes'];
            if (~exist(triangFolder, 'dir'))
                mkdir(triangFolder); % create the folder if it does not exist
            end
            
            switch obj.triOpt.method
                case 1
                    % Load selected file
                    fileLoad = sprintf('%s/%s%d.mat', triangFolder, ...
                        obj.szMeshFileRootName, obj.triOpt.meshIdx);
                    load(fileLoad, 'SimplexConnec');
                case 2
                    % Create a new triangulation
                    SimplexConnec = clCrane.fnsSelectNewTriangulation(obj.triOpt.V);
                    fileIndices = clCrane.fnsTriangMeshFileNamingIdx(...
                        obj.szMeshFileRootName, triangFolder);
                    if ~isempty(fileIndices)
                        newIdx = max(fileIndices) + 1;
                    else
                        newIdx = 1;
                    end
                    % Save the mesh
                    fileSave = sprintf('%s/%s%d.mat', triangFolder, ...
                        obj.szMeshFileRootName, newIdx);
                    save(fileSave, 'SimplexConnec');
                    obj.triOpt.meshIdx = newIdx; % save the actual mesh ID
                case 3
                    % Display existing meshes and allow user to select
                    fileIndices = clCrane.fnsTriangMeshFileNamingIdx(...
                        obj.szMeshFileRootName, triangFolder);
                    obj.triOpt.meshIdx = clCrane.fnsDisplayAllMeshes(fileIndices, ...
                        triangFolder, obj.triOpt.V, obj.szMeshFileRootName);
                    close all; % anything else inadvertently closed?
                    
                    % Load the mesh corresponding to the selection
                    fileLoad = sprintf('%s/%s%d.mat', triangFolder, ...
                        obj.szMeshFileRootName, obj.triOpt.meshIdx);
                    load(fileLoad, 'SimplexConnec');
                    
                    % Warn the user that the control values at the indices should be
                    % updated!
                    disp('Warning: U matrix with control values should be updated');
                otherwise
                    strerr = ['Invalid triangulation option selected\n'...
                        '1 - pre-specified by name\n'...
                        '2 - create new one\n'...
                        '3 - pick an existing one with visual aid'];
                    error(strerr);
            end
            
            %%  Store the triangulation
            
            nS = size(SimplexConnec, 1);
            obj.triOpt.S = zeros(2, 3, nS);
            for i = 1:nS
                c = SimplexConnec(i, :);
                obj.triOpt.S(:, :, i) = [obj.triOpt.V(c(1), :); ...
                    obj.triOpt.V(c(2), :); obj.triOpt.V(c(3), :)]';
            end
            
            obj.triOpt.tri = SimplexConnec;
        end
        
        function fnDefinePolytopeVertices(obj)
            
            % Unpack parameters
            hr = obj.param.hr;
            wr = obj.param.wr;
            ho = obj.param.ho;
            wo = obj.param.wo;
            L = obj.param.L;
            
            %% Core Outer Polytope Vertices
            
            % Avoid collision with walls and crane
            min_angle_out = -90 * pi / 180;
            min_dist_out = - wr / 2;
            mid_dist_out = -(wr / 2 - L);
            
            obj.triOpt.V(1,:) = [min_dist_out, 0];
            obj.triOpt.V(2,:) = [mid_dist_out, min_angle_out];
            obj.triOpt.V(3,:) = [-min_dist_out, min_angle_out];
            obj.triOpt.V(4,:) = [-min_dist_out, 0];
            obj.triOpt.V(5,:) = [-mid_dist_out, -min_angle_out];
            obj.triOpt.V(6,:) = [min_dist_out, -min_angle_out];
            
            % Make the adjustment to the two points on the right-top side to factor in
            % the non-linearity of the pendulum approaching the wall. This is an
            % optimized solution, cutting off as little of the original polytope
            [obj.triOpt.V(5,:), obj.triOpt.V(4,:)] = clCrane.fnsFixPolytope(L, wr / 2);
            
            %% Core Inner Polytope Vertices
            
            % Avoid collision with obstacle
            min_angle_in =  -acos((hr - ho) / L);
            min_dist_in = -(wo / 2 + sqrt(L^2 - (hr - ho)^2));
            mid_dist1_in = -wo / 2;
            mid_dist2_in = -(wo / 2 - sqrt(L^2 - (hr - ho)^2));
            
            obj.triOpt.V(7,:) = [mid_dist1_in, 0];
            obj.triOpt.V(8,:) = [mid_dist2_in, min_angle_in];
            obj.triOpt.V(9,:) = [-min_dist_in, min_angle_in];
            obj.triOpt.V(10,:) = [-mid_dist1_in, 0];
            obj.triOpt.V(11,:) = [-mid_dist2_in, -min_angle_in];
            obj.triOpt.V(12,:) = [min_dist_in, -min_angle_in];
            
            %% Make further corrections
            if obj.nVerLoaded >= 7
                % Left side of the obstacle for theta max and min
                obj.triOpt.V(1,:) = -obj.triOpt.V(4,:);
                obj.triOpt.V(2,:) = -obj.triOpt.V(5,:);
                
                % Similar to the other fix, but now for theta upper and
                % lower along the right side of the obstacle
                [obj.triOpt.V(11,:), obj.triOpt.V(9,:)] = clCrane.fnsFixPolytope2(...
                    L, wo, hr, ho);
                
                % By symmetry, now for theta upper and lower along left
                % side of obstacle
                obj.triOpt.V(8,:) = -obj.triOpt.V(11,:);
                obj.triOpt.V(12,:) = -obj.triOpt.V(9,:);
                
                % Lastly, the ones in between need to be fixed to make a
                % straight line
                obj.triOpt.V(7,:) = [mean(obj.triOpt.V([8;12],1)), 0];
                obj.triOpt.V(10,:) = [mean(obj.triOpt.V([9;11],1)), 0];
            end
            
            %% Additional Vertices
            
            % These may be set anywhere on or inside the polytope defined by the first
            % 12 vertices. In principal, any arbitrary amount of vertices can now be
            % defined, but practically define new vertices sparingly to keep complexity
            % low.
            % New vertices may be added directly to the code, but once it is, NEVER
            % delete it (it doesn't have to be used though in every triangulation)
            % Vertices are indexed by their row position in the global vertex
            % coordinate vector, and should be automatically scaled by problem
            % parameters
            
            obj.triOpt.V(13,:) = [mean(obj.triOpt.V(2:3,1)), min_angle_out];
            obj.triOpt.V(14,:) = [mean(obj.triOpt.V(5:6,1)), -min_angle_out];
        end
        
        function simplexIdx = fnCheckSimplex(obj, x)
            
            eps = 1e-6;
            
            nS = size(obj.triOpt.tri, 1);
            point = ones(nS, 2);
            point(:, 1) = point(:, 1) * x(1);
            point(:, 2) = point(:, 2) * x(3);
            
            % Convert coordinates from cartesian to barycentric coordinates
            B = cartToBary(obj.triOpt.tr, (1:length(obj.triOpt.tri))', point);
            
            % Determine the current simplex based on barycentric properties
            simplexIdx = 0;
            for i = 1:nS
                if abs(abs(B(i, 1)) + abs(B(i, 2)) + abs(B(i, 3)) - 1) < eps
                    simplexIdx = i;
                    break;
                end
            end
            
        end
        
        function fnCraneControl(obj)
            
            % fnCraneControl.m
            % Description: solver function that creates all the controllers, runs the
            % simulation with ODE45, and returns the state variable solution
            % Inputs:
            % S - the triangulation
            % x0 - initial condition
            % end_time - simulation length
            % end_distance - out of bounds check
            % param - structure of param values definining the scenario and model
            % xf - the desired equilibrium
            % U - control values at the verticies to contruct the feedback controllers
            % P3D and P4D - desired poles for 3D and 4D pole placement controllers
            % poleContr and val_SwitchToPoleContr - the method of switching to pole
            %   placement during the simulation
            % tSpanCheck - the length of time per simulation for methodODE = 2
            % methodODE - specifies the option of how to integrate the equations
            %   (infrequent vs frequent update to the controller)
            % Outputs:
            % t - a vector of the time variable
            % x - a vector of the states (as a function of time)
            % u - the control signal (as a function of time)
            % sIdx - the simplex index indicating the current simplex (as a function of
            %   time)
            % pOn - vector that indicates with a 1 when the switch to pole placement
            %   occurred
            % bExitedPolytope - a flag indicating whether the solution trajectory
            %   remained inside the polytope for all time (if yes, the flag is true)
            % U2 - updated control values at verticies using poleContr = 5,
            % automatically calculated by reverse engineering the pole placement
            % controller for the given triangulation (to be saved in the log)
            
            % Initialize
            obj.soln.t = [];
            obj.soln.x = [];
            obj.soln.u = [];
            obj.soln.sIdx = [];
            obj.soln.pOn = [];
            obj.soln.bExitedPolytope = 0;
            
            nS = size(obj.triOpt.S, 3);
            
            %% Create the controller gains
            
            % Find the gain needed to stabilize the system using pole placement
            % acker is used instead of place to support the same pole location
            % placement. acker is safe to use since the problem dimension = 4 < 10
            if obj.simPar.methodODE == 1 || obj.simPar.methodODE == 2 ...
                    || obj.simPar.methodODE == 5
                K4D = -ones(nS, 1) * acker(obj.param.A, obj.param.B, obj.contr.P4D);
                g4D = -K4D(:, 1) * obj.soln.xf(1, 1);
            elseif obj.simPar.methodODE == 4
                % Compute a different pole placement controller for each
                % simplex, using an appropriately linearized system for the
                % non-linear simulation
                K4D = zeros(nS, 4);
                g4D = zeros(nS);
                for i = 1:nS
                    % For the current simplex, compute the average theta
                    % location
                    thAvg = mean(obj.triOpt.S(2, :, nS));
                    [A, B, a] = fnLinearizeCraneSys(obj, [0; 0; thAvg; 0], 0);
                    K4D(i, :) = -acker(A, B, obj.contr.P4D);
                    
                    % See notes for derivation
                    % However, correction made because tilda state
                    % variables should have been used, so this doesn't make
                    % sense anymore!
                    g4D(i) = -K4D(i, 1) * (obj.soln.xf(1, 1)) + ...
                        (a(4) - a(2)) / (B(2) - B(4));
                    
                end
            end
            
            % Compute Triangulation Controller Parameters
            K = zeros(nS, 2);
            g = zeros(nS, 1);
            
            for i = 1:nS
                [K(i,:), g(i)] = clCrane.fnsFeedback(obj.triOpt.S(:,:, i), ...
                    obj.contr.U(i,:));
                
                % Perform a quick check whether the system will be stabilized.
                % Theoretically it is always impossible using this form of K
                evals = eig(obj.param.A + obj.param.B * [K(i, 1), 0, K(i, 2), 0]);
                clCrane.fnsSearchUnstableEigs(evals);
            end
            
            %% ODE solver
            
            % Tracks whether the conditions for switching to the pole placement
            % controller have been met
            obj.g_polesActive = 0; % reset
            
            switch(obj.simPar.methodODE)
                case 1
                    obj.g_nCountToActivate = 5;
                    % In this approach, the function handle contains the 
                    % state space model as well as the logic to choose the 
                    % appropriate feedback controller
                    obj.soln.t = (0:obj.simPar.tResolution:obj.simPar.end_time)';
                    [~, obj.soln.x] = ode45(@(t, x) fnLcranemodel(t, x, obj, ...
                        K, g, K4D, g4D), obj.soln.t, obj.soln.x0, []);
                    
                case 2
                    % In this approach, the controller is changed less often; 
                    % the equations are integrated for a small fixed length 
                    % of time before checking the current simplex and 
                    % updating the (feedback) controller to use. The switch
                    % statement resembles the one in fnLcranemodel.m
                    
                    % Set options for ODE45 function
                    %options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5 1e-5]);
                    
                    % Save the overall time and state values from the entire simulation
                    obj.soln.t = 0;
                    obj.soln.x = obj.soln.x0';
                    
                    obj.g_nCountToActivate = 1;
                    % Simulate the linearized crane model using ode45 for a
                    % set distance and time
                    while obj.soln.t(end) < obj.simPar.end_time && ...
                            obj.soln.x(end, 1) < obj.simPar.end_distance
                        
                        [K_, g_, ~, ~] = fnControlLogic(obj, obj.soln.x(end, :)', ...
                            K, g, K4D, g4D);
                        
                        T = (obj.soln.t(end):obj.simPar.tResolution:...
                            (obj.soln.t(end) + obj.simPar.tSpanCheck))';
                        [~, X] = ode45(@(t, x) fnLcranemodel2(t, x, obj, K_, g_), ...
                            T, obj.soln.x(end, :)', []);
                        
                        % Accumulate the time and state values for the animation
                        % This is inefficient due to constant memory allocations, 
                        % but speed is not critical during the simulation
                        obj.soln.t = [obj.soln.t; T(2:end)];
                        obj.soln.x = [obj.soln.x; X(2:end, :)];
                    end
                case 3
                    % Still under major experimentation
                    % Simulates the control using regulator-simplex design
                    % Only the linearized system is used
                    % The set-up has not been done to stop the cart: if the
                    % parameters below are selected correctly, the cart can
                    % clear the obstacle at best
                    
                    % Order the vertices in each simplex, v0, v1, v2
                    % Ensure o(i, 1) corresponds to the point across from
                    % exit facet for the ith simplex (ie. v0 -> F0) to know
                    % how to push trajectories
                    o = [1,2,3;
                        2,1,3;
                        1,2,3;
                        1,2,3;
                        3,2,1;
                        1,2,3;
                        2,1,3];
                    
                    % Initial conditions
                    xStart = obj.soln.x0;
                    tStart = 0;
                    tSz = obj.simPar.tResolution/10;
                    dt = 1;
                    
                    for i = 1:size(o,1)
                        % Parameters to solve related equations, hard coded for
                        % now until solvability is explored more
                        
                        if i < 6
                            poles = 1*[-1, -2, -3, -4]; % for tracking convergence
                            oG = [1; 5; 6];
                            oq = [2; 3; 4];
                            qG = [2; 7; -4];
                            q0 = [8; -5; 1];
                            % Enforce exact tracking: trouble is finding an
                            % initial guess that allows fsolve to find a
                            % solution...
%                             oG = [1; 2; 3; 5; 6];
%                             oq = [4];
%                             qG = [2; 7; 5; -3; -4];
%                             q0 = [0];
                        else
                            % Exosystem barely behaves, but plant has very
                            % hard time following, strength poles
                            poles = 2*[-1, -2, -3, -4]; % for tracking convergence
                            oG = [1; 2; 5];
                            oq = [6; 3; 4];
                            qG = [0.25; 4; 0.25];
                            q0 = [8; -5; 1];
                        end
                        
                        % Run simulation
                        t_ = (tStart:tSz:(tStart+dt))';
                        [~, x_, u_, idxExit] = clCrane.fnsRegulatorDesign(obj.param.A, ...
                            obj.param.B, t_, xStart, i, obj.triOpt.S(:,o(i,1),i), ...
                            obj.triOpt.S(:,o(i,2),i), obj.triOpt.S(:,o(i,3),i), ...
                            oG, oq, qG, q0, poles);
                        
                        if isempty(idxExit)
                            % Cannot continue because exited polytope
                            obj.soln.t = [obj.soln.t; t_];
                            obj.soln.x = [obj.soln.x; x_];
                            obj.soln.u = [obj.soln.u; u_];
                            obj.soln.sIdx = [obj.soln.sIdx; i * ones(size(t_,1), 1)];
                            break;
                        else
                            % Append solution, just up to point where it
                            % exited the simplex
                            xStart = x_(idxExit, 1:4)';
                            tStart = t_(idxExit);
                            
                            obj.soln.t = [obj.soln.t; t_(1:idxExit)];
                            obj.soln.x = [obj.soln.x; x_(1:idxExit, :)];
                            obj.soln.u = [obj.soln.u; u_(1:idxExit, :)];
                            obj.soln.sIdx = [obj.soln.sIdx; i * ones(idxExit, 1)];
                        end
                        
                    end
                    
                case 4
                    % Exactly like case 1 but using the non-linear model
                    % Here the pole placement controllers have been
                    % linearized locally
                    obj.g_nCountToActivate = 2;
                    % In this approach, the function handle contains the 
                    % state space modelas well as the logic to choose the 
                    % appropriate feedback controller
                    obj.soln.t = (0:obj.simPar.tResolution:obj.simPar.end_time)';
                    [~, obj.soln.x] = ode45(@(t, x) fnNcranemodel(t, x, obj, ...
                        K, g, K4D, g4D), obj.soln.t, obj.soln.x0, []);
                case 5
                    % Exactly like case 1 but using the non-linear model
                    % Here only one pole placement controller is used,
                    % based on linearization about zero
                    obj.g_nCountToActivate = 2;
                    % In this approach, the function handle contains the 
                    % state space modelas well as the logic to choose the 
                    % appropriate feedback controller
                    obj.soln.t = (0:obj.simPar.tResolution:obj.simPar.end_time)';
                    [~, obj.soln.x] = ode45(@(t, x) fnNcranemodel(t, x, obj, ...
                        K, g, K4D, g4D), obj.soln.t, obj.soln.x0, []);
                otherwise
                    error('Only support methodODE = 1, 2, 3, 4, or 5');
            end
            
            %% Obtain control input history over the solution trajectory
            
            nPts = size(obj.soln.t, 1);
            if obj.simPar.methodODE == 1 || obj.simPar.methodODE == 2
                
                obj.soln.u = zeros(nPts, 1);
                obj.soln.sIdx = zeros(nPts, 1);
                obj.soln.pOn = zeros(nPts, 1);
                obj.g_polesActive = 0;  % reset
                obj.g_nCountToActivate = 1;
                for i = 1:nPts
                    [~, obj.soln.u(i), obj.soln.sIdx(i), obj.soln.pOn(i)] = ...
                        fnLcranemodel(obj.soln.t(i), obj.soln.x(i, :)', ...
                        obj, K, g, K4D, g4D);
                end
            elseif obj.simPar.methodODE == 3
                obj.soln.pOn = zeros(nPts, 1);
            elseif obj.simPar.methodODE == 4 || obj.simPar.methodODE == 5
                obj.soln.u = zeros(nPts, 1);
                obj.soln.sIdx = zeros(nPts, 1);
                obj.soln.pOn = zeros(nPts, 1);
                obj.g_polesActive = 0;  % reset
                obj.g_nCountToActivate = 1;
                for i = 1:nPts
                    [~, obj.soln.u(i), obj.soln.sIdx(i), obj.soln.pOn(i)] = ...
                        fnNcranemodel(obj.soln.t(i), obj.soln.x(i, :)', ...
                        obj, K, g, K4D, g4D);
                end
            end
            
            % Check whether the trajectory ever exitted the polytope
            obj.soln.bExitedPolytope = ~isempty(find(obj.soln.sIdx == 0));
            
        end
        
        function [dx, u, sIdx, pOn] = fnLcranemodel(~, x, obj, K, g, K4D, g4D)
            
            [K_, g_, sIdx, pOn] = fnControlLogic(obj, x, K, g, K4D, g4D);
            u = K_ * x + g_;
            dx = obj.param.A * x + obj.param.B * u;
        end
        
        function [dx] = fnLcranemodel2(~, x, obj, K, g)
            
            dx = (obj.param.A + obj.param.B * K) * x + obj.param.B * g;
        end
        
        function [dx, u, sIdx, pOn] = fnNcranemodel(~, x, obj, K, g, K4D, g4D)
            
            % Nonlinear Crane Model
            
            [K_, g_, sIdx, pOn] = fnControlLogic(obj, x, K, g, K4D, g4D);
            u = K_ * x + g_;
            
            dx(1,1) = x(2);
            dx(2,1) = (obj.param.m*obj.param.L*x(4)^2*sin(x(3)) + ...
                obj.param.m*obj.param.g*cos(x(3))*sin(x(3)) + u ) / ...
                ( obj.param.M + obj.param.m * sin(x(3))^2 );
            dx(3,1) = x(4);
            dx(4,1) = -( (obj.param.M+obj.param.m)*obj.param.g*sin(x(3)) + ...
                obj.param.m*obj.param.L*x(4)^2*cos(x(3))*sin(x(3)) + ...
                u*cos(x(3)) ) / ...
                (obj.param.L*( obj.param.M + obj.param.m * sin(x(3))^2 ));
        end
        
        function [K_, g_, sIdx, pOn] = fnControlLogic(obj, x, K, g, K4D, g4D)
            
            % This function constructs the appropriate affine controller to use,
            % u = K_ * x + g_, based on the current simplex (selects the 
            % right K and g from the arrays provided) and handles when to 
            % switch to the pole placement controller, again providing 
            % the appropriate K and g based on the selected poles and 
            % equilibrium position
            
            pOn = 0; % tracks whether the switch to pole placement 
                     %  controller occurred
            
            % Update the controller to use
            sIdx = fnCheckSimplex(obj, x);	% Control using triangulation
            
            if obj.g_polesActive < obj.g_nCountToActivate
                % Use simplex controllers before switching to pole
                % placement controller
                if sIdx == 0
                    K_ = [0 0 0 0];
                    g_ = 0;
                else
                    K_ = [K(sIdx, 1), 0, K(sIdx, 2), 0];
                    g_ = g(sIdx);
                end
            else
                % Use 4D pole placement
                if sIdx == 0
                    K_ = [0 0 0 0];
                    g_ = 0;
                else
                    K_ = K4D(sIdx, :); % Recall place/acker gives the negative
                    g_ = g4D(sIdx);
                end
            end
            
            % Determine whether to turn on pole placement controller
            % (permanently)
            if obj.g_polesActive < obj.g_nCountToActivate
                switch obj.contr.poleContr
                    case 4
                        % Switch to pole placement when desired z position is
                        % reached
                        if x(1) > obj.contr.val_SwitchToPoleContr
                            obj.g_polesActive = obj.g_polesActive + 1;
                            pOn = 1;
                        end
                    case 5
                        % Switch to pole placement when desired simplex is reached
                        if sIdx == round(obj.contr.val_SwitchToPoleContr)
                            obj.g_polesActive = obj.g_polesActive + 1;
                            pOn = 1;
                        end
                    otherwise
                        error('Only support option "4" and "5"');
                end
            end
        end
        
        function fnCraneAnimation(obj)
            
            % Function name: fnCraneAnimation
            % Function m-file: fnCraneAnimation.m
            % Function Description: given a state variable solution along with the time
            %   variable, displays an animation in real-time of the solution along with
            %   some plots of states vs. each other or time.  Also shows the control
            %   signal vs. time and the velocity vs. angular velocity trajectory. The
            %   problem scenario is scaled by the parameters values in param.
            %   Inputs:
            %   t - independent variable (time in seconds)
            %   x(:, 1) - horizontal position of the cart in meters
            %   x(:, 2) - horizontal velocity of the cart in m/s
            %   x(:, 3) - angle of the pendulum arm wrt to a vertical reference line
            %       in radians
            %   x(:, 4) - angular velocity in rad/s.
            %   u - control signal as a function of time
            %   sIdx - the simplex index indicating the current simplex (as a function
            %       of time)
            %   pOn - vector that indicates with a 1 when the switch to pole placement
            %       occurred
            %   S - the triangulation coordinate vertices
            %   param - structure of param values definining the scenario and model
            %   end_time - simulation length
            %   nFramesPerSecond - number of frames or updates to make to the animation
            %       per second to make it seem like a movie
            %   fTimeFlowFactor - speed up or slow down factor for the animation
            %   graphLineStyle and trackLineStyle - MATLAB string for a plotting
            %       setting
            %   bDrawTrackCtrlChangePts - true: will draw circles and star on the
            %       animation trajectory denoting where a new simplex is entered and
            %       when the change to the pole placement controller starts,
            %       respectively. Warning: since the animation time is sampled, several
            %       simplices may be passed at the same frame, so do not count circles
            %       to determine the simplex regions exactly
            
            % Written by K. Wedeward Fall 96
            % Modified by J. Watkins Fall 96
            
            % Source: http://www.ee.nmt.edu/~wedeward/EE443L/FA00/Ballbeam.m
            
            % Modified by N. Karunaharan, H. Kolahdoozan, S. Nagaratnam, P.
            % Thillaichandran 2012.
            
            % Additional source: 
            %  http://www.princeton.edu/~dsavrans/docs/MATLAB_Animation.pdf
            
            % Modified by Mario Vukosavljev, Dec 27 2012
            
            if obj.soln.x(1,1) == obj.soln.x(end,1)
                disp('No movement, so animation will not be displayed');
                return;
            end
            disp('Starting the Animation');
            
            %% Select the (equally time spaced) data points for the animation
            
            movieIndices = clCrane.fnsSelectAnimationFrames(obj.soln.t, ...
                obj.animPar.nFramesPerSecond, obj.animPar.fTimeFlowFactor);
            nFrames = size(movieIndices, 1);
            
            %% Cart and bob positions
            
            % Obtain expressions for the bob and cart position in the
            % world frame, based on the motion of the state variables
            % The bob refers to the center of the bob, and the cart refers 
            % to the bottom center of the cart
            % Obtain values only for the animation frames
            
            bobPos = [obj.soln.x(movieIndices, 1) + ...
                obj.param.L * sin(obj.soln.x(movieIndices, 3)), ...
                -obj.param.L * cos(obj.soln.x(movieIndices, 3))];
            
            cartPosX = obj.soln.x(movieIndices, 1);
            % Y position of the bottom of the cart is always zero
            
            %% Create figures for displaying graphs
            
            scrSz = get(0,'ScreenSize');
            h = figure('Name', 'Animation', 'units', 'normalized', ...
                'outerposition', [0 0 scrSz(4)/scrSz(3) 1]); %, 'visible', 'off');
            
            % First plot is of the polytope in (z, theta) coordinates
            h11 = subplot(4, 4, [1, 2, 5, 6]);
            clCrane.fnsDrawSimplices(obj.triOpt.S, obj.triOpt.V, 0, []);
            h111 = plot(obj.soln.x(1, 1), obj.soln.x(1, 3), ...
                obj.animPar.graphLineStyle, 'Color', 'b', 'LineWidth', 2);
            
            % Second plot contains subplots of each state as a function of time
            
            % Determine max and min values that occur for each state to determine 
            % the plot size (to remain fixed throughout the animation)
            minVals = min(obj.soln.x);
            maxVals = max(obj.soln.x);
            buffSz = 0.1 * (maxVals - minVals);
            maxVals = maxVals + buffSz;
            minVals = minVals - buffSz;
            
            h21 = subplot(4, 4 , 3);
            h211 = plot(obj.soln.t(1), obj.soln.x(1, 1), ...
                obj.animPar.graphLineStyle, 'Color', 'k');
            axis(h21, [0, obj.simPar.end_time, minVals(1), maxVals(1)]);
            title('$z$ (m) vs. $t$ (s)', 'Interpreter', 'Latex', ...
                'LineWidth', 2);
            grid on;
            
            h22 = subplot(4, 4, 4);
            h221 = plot(obj.soln.t(1), obj.soln.x(1, 2), ...
                obj.animPar.graphLineStyle, 'Color', 'k');
            axis(h22, [0, obj.simPar.end_time, minVals(2), maxVals(2)]);
            title('$\dot{z}$ (m/s) vs. $t$ (s)', 'Interpreter', 'Latex', ...
                'LineWidth', 2);
            grid on;
            
            h23 = subplot(4, 4, 7);
            h231 = plot(obj.soln.t(1), obj.soln.x(1, 3), ...
                obj.animPar.graphLineStyle, 'Color', 'k');
            axis(h23, [0, obj.simPar.end_time, minVals(3), maxVals(3)]);
            title('$\theta$ (rad) vs. $t$ (s)', 'Interpreter', 'Latex', ...
                'LineWidth', 2);
            grid on;
            
            h24 = subplot(4, 4, 8);
            h241 = plot(obj.soln.t(1), obj.soln.x(1, 4), ...
                obj.animPar.graphLineStyle, 'Color', 'k');
            axis(h24, [0, obj.simPar.end_time, minVals(4), maxVals(4)]);
            title('$\dot{\theta}$ (rad/s) vs. $t$ (s)', 'Interpreter', 'Latex', ...
                'LineWidth', 2);
            grid on;
            
            %% Set up the animation objects
            
            h31 = subplot(4, 4, 9:16);
            
            % Create circles and a square, with arbitrary scaling that look aesthetic:
            theta = 0 : (pi / 100) : (2 * pi);
            s = (obj.param.hr / 10) * [-1.5, -1.5, 1.5, 1.5, -1.5; 0, 2, 2, 0, 0];
            xscirc = (obj.param.hr / 20) * sin(theta); % for the pendulum bob
            yscirc = (obj.param.hr / 20) * cos(theta);
            wheelrad = obj.param.hr / 28;
            xwcirc = wheelrad * sin(theta); % for the cart wheel
            ywcirc = wheelrad * cos(theta);
            xssquare = s(1, :);
            yssquare = s(2, :);
            
            % Create the walls
            minz = - obj.param.wr / 2;
            maxz = obj.param.wr / 2;
            miny = - obj.param.hr;
            maxy = 0;
            line([minz, minz], [maxy, miny], 'Color', 'k', 'LineWidth', 2), hold on;
            line([maxz, maxz], [maxy, miny], 'Color', 'k', 'LineWidth', 2), hold on;
            line([minz, maxz], [miny, miny], 'Color', 'k', 'LineWidth', 2), hold on;
            line([minz, maxz], [maxy, maxy], 'Color', 'k', 'LineWidth', 2), hold on;
            
            % Create the obstacle
            % vertical space between ceiling and obstacle
            space = -(obj.param.hr - obj.param.ho);
            obsHalfW = obj.param.wo / 2;
            patch(obsHalfW * [-1, 1, 1, -1, -1], [miny miny space space miny], ...
                0.7 * ones(1, 3)), hold on;
            
            % Create the cart and pendulum components
            
            % Connect mass and pivot:
            % Note: the bob trace should be plotted first so that the other objects
            % overlap it if they cross; this is the same reason why the scenario walls
            % and obstacle were drawn first
            trace = plot(bobPos(1, 1), bobPos(1, 2), obj.animPar.trackLineStyle, ...
                'Color', 'b', 'LineWidth', 2);
            arm = line([cartPosX(1, 1), bobPos(1, 1)], [0, bobPos(1, 2)], ...
                'Color', 'r', 'LineWidth', 2);
            cart = patch(xssquare + cartPosX(1, 1), yssquare, 'g');
            bob = patch(xscirc + bobPos(1, 1), yscirc + bobPos(1, 2), 'b');
            % Bonus: create the wheels of the cart
            wheelShift = (max(xssquare) - min(xssquare)) / 2 - wheelrad;
            lwheel_circle = patch(xwcirc + cartPosX(1, 1) - wheelShift, ...
                ywcirc - wheelrad, 'k');
            rwheel_circle = patch(xwcirc + cartPosX(1, 1) + wheelShift, ...
                ywcirc - wheelrad, 'k');
            
            % Set axes to proper values:
            axis('equal');
            grid('on');
            axis([1.1 * minz, 1.1 * maxz, miny, -miny/2]);
            
            xlabel('Cart Position (m)', 'Interpreter', 'Latex');
            ylabel('Vertical Position (m)', 'Interpreter', 'Latex');
            h311 = title(sprintf('\n'), 'Interpreter', 'Latex');
            
            %% Classify the trajectory into pieces over the time interval
            
            % Determine at what time the trajectory switched to another simplex
            controlSw = {1, obj.soln.sIdx(1)}; % start time index, simplex 
                                               %  index (0 means out)
            nRegions = 1;
            for i = 2:size(obj.soln.t, 1)
                % Check for a controller switch
                if (obj.soln.sIdx(i) ~= controlSw{nRegions, 2})
                    nRegions = nRegions + 1;
                    controlSw{nRegions, 1} = i;
                    controlSw{nRegions, 2} = obj.soln.sIdx(i);
                end
            end
            
            % Show where the controllers changed to pole placement
            switchIdx = find(obj.soln.pOn == 1);
            
            movie_switchEvents = zeros(size(movieIndices, 1), 1);
            if obj.animPar.bDrawTrackCtrlChangePts
                % For the animation, determine the appropriate sampled frame 
                % times when these swithces occured
                for i = 2:nRegions
                    [~, idx] = min(abs(controlSw{i, 1} - movieIndices));
                    movie_switchEvents(idx) = 1;
                end
                
                % Pole placement switch in sampled time
                if (~isempty(switchIdx))
                    [~, idx] = min(abs(switchIdx - movieIndices));
                    if movie_switchEvents(idx) == 1
                        movie_switchEvents(idx) = 3; % both simplex region 
                                                     %  and pole pl. change occurred
                    else
                        movie_switchEvents(idx) = 2; % only the pole pl. change occurred
                    end
                end
            end
            
            %% Run the animation
            
            invTimeFactor = (1 / obj.animPar.fTimeFlowFactor);
            actualFPS = round(obj.animPar.fTimeFlowFactor * nFrames / obj.soln.t(end));
            flagNegPause = 0;
            
            if ~isempty(obj.animPar.movieSaveName)
                % Initialize movie
                
                numberStr = sprintf('%.2f', obj.animPar.fTimeFlowFactor);
                % The only dot allowed is the extension
                numberStr = strrep(numberStr, '.', 'dot');
                
                movieFullPath = sprintf('%s_timeFlow_x%s', ...
                    obj.animPar.movieSaveName, numberStr);
                writerObj = VideoWriter(movieFullPath, 'MPEG-4');
                writerObj.FrameRate = actualFPS;
                open(writerObj);
            end
            
            for i = 1 : nFrames - 1
                
                tFrameProc = tic;
                j = movieIndices(i); % obtain the appropriate index for this frame
                range = movieIndices(1:i);
                
                % Update the supporting plots
                % Plot the cart position vs pendulum angle
                set(h111, 'XData', obj.soln.x(range, 1));
                set(h111, 'YData', obj.soln.x(range, 3));
                
                % Plot the states as a function of time
                set(h211, 'XData', obj.soln.t(range));
                set(h211, 'YData', obj.soln.x(range, 1));
                set(h221, 'XData', obj.soln.t(range));
                set(h221, 'YData', obj.soln.x(range, 2));
                set(h231, 'XData', obj.soln.t(range));
                set(h231, 'YData', obj.soln.x(range, 3));
                set(h241, 'XData', obj.soln.t(range));
                set(h241, 'YData', obj.soln.x(range, 4));
                
                % Update moving elements
                set(arm, 'XData', [cartPosX(i, 1), bobPos(i, 1)]);
                set(arm, 'YData', [0, bobPos(i, 2)]);
                set(cart, 'XData', xssquare + cartPosX(i, 1));
                set(bob, 'XData', xscirc + bobPos(i, 1));
                set(bob, 'YData', yscirc + bobPos(i, 2));
                set(trace, 'XData', bobPos(1:i, 1));
                set(trace, 'YData', bobPos(1:i, 2));
                set(lwheel_circle, 'XData', xwcirc + cartPosX(i, 1) - wheelShift);
                set(rwheel_circle, 'XData', xwcirc + cartPosX(i, 1) + wheelShift);
                
                % Draw the simplex controller type and region changes
                switch movie_switchEvents(i)
                    case 1 % simplex region change
                        hold on;
                        htemp = plot(bobPos(i, 1), bobPos(i, 2), 'ko');
                        uistack(htemp, 'bottom');
                    case 2 % change to pole placement
                        hold on;
                        htemp = plot(bobPos(i, 1), bobPos(i, 2), 'k*');
                        uistack(htemp, 'bottom');
                    case 3 % both changes
                        hold on;
                        htemp = plot(bobPos(i, 1), bobPos(i, 2), 'ko', ...
                            bobPos(i, 1), bobPos(i, 2), 'k*');
                        uistack(htemp, 'bottom');
                end
                
                strTitle = sprintf('Crane Problem Animation\n');
                strTitle = sprintf('%sAnimation time = %.1f s', strTitle, ...
                    obj.soln.t(j));
                set(h311, 'String', strTitle, 'FontSize', 12);
                
                if ~isempty(obj.animPar.movieSaveName)
                    writeVideo(writerObj, getframe(h)); % save a movie frame
                end
                
                % Calculate how much to pause to simulate real flow of time
                tBtwFrames = (obj.soln.t(movieIndices(i + 1)) - ...
                    obj.soln.t(j) - toc(tFrameProc)) * invTimeFactor;
                if (tBtwFrames > 0)
                    pause(tBtwFrames);
                else
                    flagNegPause = 1;
                    pause(0.0001);
                end
            end
            
            if (flagNegPause)
                disp('Warning: real flow of time was not shown');
                disp('Try to decrease the number of frames per second');
            end
            
            if ~isempty(obj.animPar.movieSaveName)
                % Complete the movie
                close(writerObj);
            end
            
            % Save the final image (as vector format, for latex)
            set(h,'Units', 'Inches');
            pos = get(h, 'Position');
            set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
                'PaperSize', [pos(3), pos(4)]);
            print(h, '-dpdf', [obj.szPlotFolder 'pdfFig']);
            
            %% Plot control input history over the solution trajectory
            
            minu = min(obj.soln.u);
            maxu = max(obj.soln.u);
            padu = 0.1 * (maxu - minu);
            uline = [minu - padu, maxu + padu];
            
            % Plot the control signal
            hCon = figure;
            plot(obj.soln.t, obj.soln.u, 'b-', 'LineWidth', 2);
            
            % Plot when the switch to pole placement occurred
            if ~isempty(switchIdx)
                hold on, plot(obj.soln.t(switchIdx), obj.soln.u(switchIdx), 'k*');
            end
            
%             % Label the regions by the simplex index
%             for i = 1:nRegions-1
%                 tRange = [obj.soln.t(controlSw{i, 1}), obj.soln.t(controlSw{i+1, 1})];
%                 hold on, plot([tRange(2), tRange(2)], uline, 'k--');
%                 hold on, text(mean(tRange), minu - padu, ...
%                     sprintf('%d', controlSw{i, 2}), 'BackgroundColor', [.7 .9 .7]);
%             end
%             hold on, text(mean([obj.soln.t(controlSw{end, 1}), obj.soln.t(end)]), ...
%                 minu - padu, sprintf('%d', controlSw{end, 2}), ...
%                 'BackgroundColor', [.7 .9 .7]);
            
            title('Control input');
            xlabel('Time (s)');
            ylabel('Control input (N)');
            grid on;
            
            % Save the final image (as vector format, for latex)
            set(hCon,'Units', 'Inches');
            pos = get(hCon, 'Position');
            set(hCon, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
                'PaperSize', [pos(3), pos(4)]);
            print(hCon, '-dpdf', [obj.szPlotFolder 'contrSignal']);
            
            %% Plot the ang_velocity vs. velocity trajectory
            % Show in particular which simplex the trajectory was in, using color
            
            figure;
            
            unique_sIdx = unique(obj.soln.sIdx); % simplex indices that 
                                                 %  were visited, 0 means out
            nIdx = size(unique_sIdx, 1);
            cc = hsv(nIdx); % obtain a set of colors
            
            legendstrs = cell(nIdx, 1);
            for i = 1:nIdx
                legendstrs{i} = sprintf('%d', unique_sIdx(i));
            end
            
            % Plot the trajectory piece-wise for each region encountered, and mark 
            % with a color (already created above, complete with a legend)
            hLegend = -1 * ones(nIdx, 1);
            for i = 1:nRegions-1
                v = obj.soln.x(controlSw{i, 1}:controlSw{i+1, 1}, [2, 4]);
                pos = find(unique_sIdx == controlSw{i, 2});
                hold on;
                % Save the handle while plotting
                % Note: this may overwrite a handle, but ok since they 
                % represent the same simplex region
                hLegend(pos) = plot(v(:, 1), v(:, 2), 'Color', cc(pos, :), ...
                    'LineStyle', '-');
            end
            v = obj.soln.x(controlSw{end, 1}:end, [2, 4]);
            pos = find(unique_sIdx == controlSw{end, 2});
            hold on;
            hLegend(pos) = plot(v(:, 1), v(:, 2), 'Color', cc(pos, :), ...
                'LineStyle', '-');
            legend(hLegend, legendstrs); % display the legend
            
            % Display where the switch to pole placement controller occurred
            if ~isempty(switchIdx)
                hold on, plot(obj.soln.x(switchIdx, 2), obj.soln.x(switchIdx, 4), 'k*');
            end
            
            % Display the start and end velocities (for clarity, although 
            % it is always (0, 0)
            hold on, plot(0, 0, 'ko');
            
            title('Velocity trajectory');
            xlabel('Velocity (m/s)');
            ylabel('Angular velocity (rad/s)');
        end
        
        function fnCreateLogFileOfTest(obj, resultsFolderName)
            %% Source Code Version Number
            
            % Increment by one every time this file is modified
            % VERSION = 1; % baseline
            % VERSION = 2; % added triangulation creating tools
            % VERSION = 3; % added polytope fix, simplified ode45 call,
            %   and added 4.2 option (switch in simplex
            %   "zVal_SwitchToPoleContr")
            % VERSION = 4; % added a second way to integrate the ODEs,
            %   using fnLcranemodel2
            % VERSION = 5; % removed the pass/fail print at the bottom of
            %   the file, renamed "zVal_SwitchToPoleContr" to
            %   "val_SwitchToPoleContr", generalized to an arbitrary size
            %   for U/number of vertices, and corrected to also print
            %   complex part of poles
            % VERSION = 6; % added the tResolution variable.
            % nVersion = 7; % huge update to MATLAB class format
            
            if isempty(resultsFolderName)
                return;
            end
            
            %% Construct the save path and file
            
            % Convert backslashes to forward slashes for portability and ensure 
            % there is a forward slash at the end
            currentFolder = strrep(pwd, '\', '/');
            if (currentFolder(end) ~= '/')
                currentFolder = [currentFolder '/'];
            end
            
            % Verify the results folder name provided
            if (~isempty(strfind(resultsFolderName, ' ')))
                error('No spaces are allowed in the results folder name');
            end
            
            if (~isempty(strfind(resultsFolderName, '\')))
                error('No backslashes are allowed in the results folder name');
            end
            
            % Name of the directory to save the results in:
            dirSave = [currentFolder resultsFolderName];
            if (dirSave(end) ~= '/')
                dirSave = [dirSave '/'];
            end
            
            % Create the directory if it does not exist
            if (exist(dirSave, 'dir') == 0)
                [s, mess, ~] = mkdir(dirSave);
                if (s ~= 1)
                    error(mess);
                end
            end
            
            % Pause at least one second to ensure all files have a unique name, in case
            % the simulation runs extremely fast
            pause(1);
            
            % Create a timestamp of the test
            cl = clock;
            timestampStr = sprintf('%d%s%s-%s%s%s', ...
                cl(1), clCrane.fnsTwoDigitPaddedStr(cl(2)),...
                clCrane.fnsTwoDigitPaddedStr(cl(3)), ...
                clCrane.fnsTwoDigitPaddedStr(cl(4)), ...
                clCrane.fnsTwoDigitPaddedStr(cl(5)), ...
                clCrane.fnsTwoDigitPaddedStr(round(cl(6))));
            
            % String to be included in the file name denoting the version
            verStr = sprintf('Ver-%d', obj.nVersion);
            
            % String denoting pass or fail depending on whether it exited polytope
            if obj.soln.bExitedPolytope
                strTestRes = 'Fail';
            else
                strTestRes = 'Pass';
            end
            
            % Create the file name
            filename = sprintf('%s_%s%d_%s_%s.txt', verStr, obj.szMeshFileRootName, ...
                obj.triOpt.meshIdx, timestampStr, strTestRes);
            
            %% Create file log
            
            % Create the file handle
            fid = fopen([dirSave filename], 'w');
            
            fprintf(fid, 'soln.x0 = [%f; %f; %f; %f];\n', obj.soln.x0);
            fprintf(fid, 'soln.xf = [%f; %f; %f; %f];\n', obj.soln.xf);
            
            fprintf(fid, 'simPar.end_time = %f;\n', obj.simPar.end_time);
            fprintf(fid, 'simPar.end_distance = %f;\n', obj.simPar.end_distance);
            fprintf(fid, 'simPar.tSpanCheck = %f;\n', obj.simPar.tSpanCheck);
            fprintf(fid, 'simPar.methodODE = %d;\n', obj.simPar.methodODE);
            fprintf(fid, 'simPar.tResolution = %f;\n', obj.simPar.tResolution);
            
            fprintf(fid, 'param.hr = %f;\n', obj.param.hr);
            fprintf(fid, 'param.wr = %f;\n', obj.param.wr);
            fprintf(fid, 'param.ho = %f;\n', obj.param.ho);
            fprintf(fid, 'param.wo = %f;\n', obj.param.wo);
            fprintf(fid, 'param.M = %f;\n', obj.param.M);
            fprintf(fid, 'param.m = %f;\n', obj.param.m);
            fprintf(fid, 'param.g = %f;\n', obj.param.g);
            fprintf(fid, 'param.pLambda = %f;\n', obj.param.pLambda);
            
            fprintf(fid, ...
                'contr.P4D = [(%f + %fi), (%f + %fi), (%f + %fi), (%f + %fi)];\n', ...
                real(obj.contr.P4D(1)), imag(obj.contr.P4D(1)), ...
                real(obj.contr.P4D(2)), imag(obj.contr.P4D(2)), ...
                real(obj.contr.P4D(3)), imag(obj.contr.P4D(3)), ...
                real(obj.contr.P4D(4)), imag(obj.contr.P4D(4)));
            fprintf(fid, 'contr.poleContr = %f;\n', obj.contr.poleContr);
            fprintf(fid, 'contr.val_SwitchToPoleContr = %f;\n', ...
                obj.contr.val_SwitchToPoleContr);
            
            for i = 1:size(obj.contr.U, 1)
                fprintf(fid, 'contr.U(%d, :) = [%f, %f, %f];\n', i, obj.contr.U(i, :));
            end
            fprintf(fid, '\n');
            
            fclose(fid);
        end
        
        function fnParseLogFile(obj, fName)
            
            if isempty(fName)
                return;
            end
            
            % Reads in all the variables stored in the file
            % Each line has the form "varname" = value(s)
            
            %% Open tbe file
            
            fName = strrep(fName, '\', '/'); % for consistency
            
            fid = fopen(fName, 'r');
            if fid == -1
                error('Could not open %s\n', fName);
            end
            
            %% Initialize necessary variables
            
            obj.nVerLoaded = 1; % This will be determined from the file name
            obj.soln.testID = '';
            obj.triOpt.method = 1; % This will not be changed
            obj.triOpt.meshIdx = 1; % Will be read from file name
            obj.contr.U = [];
            
            %% Read necessary contents from the filename
            
            lastBackslash = strfind(fName, '/');
            if isempty(lastBackslash)
                fHeader = fName;
            else
                fHeader = fName(lastBackslash(end):end);
            end
            underScores = strfind(fHeader, '_'); % Delimiter between quantities 
                                                 %  in file name
            
            % Version number
            verStart = strfind(fHeader, 'Ver-');
            obj.nVerLoaded = str2num(fHeader(verStart+4:underScores(1)-1));
            
            % Mesh index
            if obj.nVerLoaded > 1
                % Version 1 only used Mesh1 before triangulation tools were created
                
                meshStart = strfind(fHeader, 'Mesh');
                obj.triOpt.meshIdx = str2num(fHeader(meshStart+4:underScores(2)-1));
            end
            
            % Test ID
            obj.soln.testID = fHeader(underScores(end-1)+1:underScores(end)-1);
            
            %% Parse the file
            
            str = fgets(fid);
            
            while ischar(str)
                eq = strfind(str, '=');
                varName = strtrim(str(1:eq-1));
                val = strtrim(str(eq+1:end));
                
                if obj.nVerLoaded < 7
                    % For backwards compatibility
                    while ~isempty(strfind(varName, 'U('))
                        % Read all control values; there are at least 12 lines, but
                        % there may be more (unknown amount)
                        obj.contr.U(end+1, :) = sscanf(val, '[%f, %f, %f]');
                        
                        % Update to next line
                        str = fgets(fid);
                        eq = strfind(str, '=');
                        varName = strtrim(str(1:eq-1));
                        val = strtrim(str(eq+1:end));
                    end
                    
                    switch(varName)
                        case 'IC'
                            obj.soln.x0 = sscanf(val, '[%f; %f; %f; %f]');
                        case 'equil_x'
                            obj.soln.xf = sscanf(val, '[%f; %f; %f; %f]');
                        case 'end_time'
                            obj.simPar.end_time = sscanf(val, '%f');
                        case 'end_distance'
                            obj.simPar.end_distance = sscanf(val, '%f');
                        case 'tSpanCheck'
                            obj.simPar.tSpanCheck = sscanf(val, '%f');
                        case 'methodODE'
                            obj.simPar.methodODE = sscanf(val, '%d');
                        case 'tResolution'
                            obj.simPar.tResolution = sscanf(val, '%f');
                        case 'param.room_height'
                            obj.param.hr = sscanf(val, '%f');
                        case 'param.room_width'
                            obj.param.wr = sscanf(val, '%f');
                        case 'param.obs_height'
                            obj.param.ho = sscanf(val, '%f');
                        case 'param.obs_width'
                            obj.param.wo = sscanf(val, '%f');
                        case 'param.mcart'
                            obj.param.M = sscanf(val, '%f');
                        case 'param.mpend'
                            obj.param.m = sscanf(val, '%f');
                        case 'param.gravity'
                            obj.param.g = sscanf(val, '%f');
                        case 'param.pLambda'
                            obj.param.pLambda = sscanf(val, '%f');
                        case 'P4D'
                            if obj.nVerLoaded < 5
                                % For backwards compatibility, before saving complex part
                                % of the poles
                                obj.contr.P4D = sscanf(val, '[%f, %f, %f, %f]');
                            else
                                tmp = sscanf(val, ...
                                    '[(%f + %fi), (%f + %fi), (%f + %fi), (%f + %fi)]');
                                obj.contr.P4D = [tmp(1) + 1i * tmp(2), ...
                                    tmp(3) + 1i * tmp(4), ...
                                    tmp(5) + 1i * tmp(6), ...
                                    tmp(7) + 1i * tmp(8)];
                            end
                        case 'poleContr'
                            val = sscanf(val, '%f');
                            % Redefined the options after VERSION 6, so fix
                            % it here to be consistent with the new options
                            switch(val)
                                case 3
                                    error('No longer supported, but useless anyway');
                                case 4
                                    obj.contr.poleContr = 4;
                                case 4.1
                                    error('No longer supported, but subset of case 4.2');
                                case 4.2
                                    obj.contr.poleContr = 5;
                                case 5
                                    error('Never used/implemented');
                                otherwise
                                    error('Getting here should not be possible');
                            end
                        case 'val_SwitchToPoleContr'
                            obj.contr.val_SwitchToPoleContr = sscanf(val, '%f');
                        case 'zVal_SwitchToPoleContr'
                            % Backwards compatibility before VERSION 5
                            obj.contr.val_SwitchToPoleContr = sscanf(val, '%f');
                    end
                    
                else
                    while ~isempty(strfind(varName, 'contr.U('))
                        % Read all control values; there are at least 12 lines, but
                        % there may be more (unknown amount)
                        obj.contr.U(end+1, :) = sscanf(val, '[%f, %f, %f]');
                        
                        % Update to next line
                        str = fgets(fid);
                        eq = strfind(str, '=');
                        varName = strtrim(str(1:eq-1));
                        val = strtrim(str(eq+1:end));
                    end
                    
                    switch(varName)
                        case 'soln.x0'
                            obj.soln.x0 = sscanf(val, '[%f; %f; %f; %f]');
                        case 'soln.xf'
                            obj.soln.xf = sscanf(val, '[%f; %f; %f; %f]');
                        case 'simPar.end_time'
                            obj.simPar.end_time = sscanf(val, '%f');
                        case 'simPar.end_distance'
                            obj.simPar.end_distance = sscanf(val, '%f');
                        case 'simPar.tSpanCheck'
                            obj.simPar.tSpanCheck = sscanf(val, '%f');
                        case 'simPar.methodODE'
                            obj.simPar.methodODE = sscanf(val, '%d');
                        case 'simPar.tResolution'
                            obj.simPar.tResolution = sscanf(val, '%f');
                        case 'param.hr'
                            obj.param.hr = sscanf(val, '%f');
                        case 'param.wr'
                            obj.param.wr = sscanf(val, '%f');
                        case 'param.ho'
                            obj.param.ho = sscanf(val, '%f');
                        case 'param.wo'
                            obj.param.wo = sscanf(val, '%f');
                        case 'param.M'
                            obj.param.M = sscanf(val, '%f');
                        case 'param.m'
                            obj.param.m = sscanf(val, '%f');
                        case 'param.g'
                            obj.param.g = sscanf(val, '%f');
                        case 'param.pLambda'
                            obj.param.pLambda = sscanf(val, '%f');
                        case 'contr.P4D'
                            if obj.nVerLoaded < 5
                                % For backwards compatibility, before saving complex part
                                % of the poles
                                obj.contr.P4D = sscanf(val, '[%f, %f, %f, %f]');
                            else
                                tmp = sscanf(val, ...
                                    '[(%f + %fi), (%f + %fi), (%f + %fi), (%f + %fi)]');
                                obj.contr.P4D = [tmp(1) + 1i * tmp(2), ...
                                    tmp(3) + 1i * tmp(4), ...
                                    tmp(5) + 1i * tmp(6), ...
                                    tmp(7) + 1i * tmp(8)];
                                
                            end
                        case 'contr.poleContr'
                            obj.contr.poleContr = sscanf(val, '%f');
                        case 'contr.val_SwitchToPoleContr'
                            obj.contr.val_SwitchToPoleContr = sscanf(val, '%f');
                        case 'contr.zVal_SwitchToPoleContr'
                            % Backwards compatibility before VERSION 5
                            obj.contr.val_SwitchToPoleContr = sscanf(val, '%f');
                    end
                    
                end
                
                str = fgets(fid); % get next line
            end
        end
        
    end
    methods(Static)
        
%         fnsSearchUnstableEigs(lambda)
%         fnsNonLin_ComplexConstraints(hr, wr, ho, wo, L)
%         fnsSaveSchematic(hr, wr, ho, wo, L, ...
%             drawPos_z, drawPos_th, labelFlag)
%         fnsTwoDigitPaddedStr(num)
%         chosenIdx = fnsDisplayAllMeshes(fileIndices, triangFolder, V, ...
%             meshFileRootName)
%         fnsDrawSimplices(S, V, bDetailed, hDraw)
%         [K,g] = fnsFeedback(v,u)
%         [Vleft, Vright] = fnsFixPolytope(pL, halfW)
%         [Vleft, Vright] = fnsFixPolytope2(L, wo, hr, ho)
%         bIntersect = fnsNonLinPolytope_Intersection(x1c, x3c, lx, ly)
%         area = fnsPolygonAreaIntersection(set1, set2)
%         bIntersect = fnsTriangulation_Intersections(X, Y, xObs, yObs, ...
%             chosenVertices, SimplexConnec, nS)
%         SimplexConnec = fnsSelectNewTriangulation(V)
%         movieIndices = fnsSelectAnimationFrames(t, nFramesPerSecond, fTimeFlowFactor)
%         fileIndices = fnsTriangMeshFileNamingIdx(fileRootName, path)
        [t, x, u, idxExit] = fnsRegulatorDesign(A, B, t, x0, sIdx, v0, v1, v2, ...
            oG, oq, qG, q0, poles)
        
        function fileIndices = fnsTriangMeshFileNamingIdx(fileRootName, path)
            
            nChar = size(fileRootName, 2);
            list = dir(path);
            
            % Store the file IDs (can analyse distribution if necessary)
            fileIndices = zeros(size(list, 1), 1); 
            nIDs = 0;
            for i = 1:size(list, 1)
                if (list(i).isdir)
                    continue;
                end
                
                k = strfind(list(i).name, '.mat');
                if (~isempty(k))
                    % locate the period
                    
                    [x, status] = str2num(list(i).name(nChar + 1 : k - 1));
                    if (~status)
                        continue;
                    end
                    nIDs = nIDs + 1;
                    fileIndices(nIDs) = x;
                end
            end
            
            fileIndices = fileIndices(find(fileIndices > 0)); % clean up
            fileIndices = sort(fileIndices); % sort by increasing order for clarity
        end
        
        function movieIndices = fnsSelectAnimationFrames(t, nFramesPerSecond, ...
                fTimeFlowFactor)
            
            % Construct a subset of the given time and state data, spaced to produce
            % the desired number of frames per second
            nDataPoints = size(t, 1); % == size(x, 1)
            tIdeal = (0 : (fTimeFlowFactor / nFramesPerSecond) : t(end));
            nFrames = size(tIdeal, 2);  % total number of frames in the animation
            movieIndices = zeros(nFrames, 1); % Indices referring to data points 
                                              %  selected for the animation
            
            k = 0;
            pos = 1; % as time increases, can track the latest index
            for i = tIdeal
                % Find the closest time to the current time and store the index
                prevDist = abs(i - t(pos));
                for j = pos + 1 : nDataPoints
                    curDist = abs(i - t(j));
                    if (curDist > prevDist)
                        % done; the previous data point was the minimum
                        k = k + 1;
                        movieIndices(k) = j - 1;
                        pos = j - 1; % set the next appropriate starting point
                        break;
                    else
                        % update
                        prevDist = curDist;
                    end
                end
            end
            
            % Post process the indices: eliminate the zeros at the end and repeated
            % indices
            movieIndices = movieIndices(movieIndices > 0);
            movieIndices = unique(movieIndices);
            
            trueFPS = round(fTimeFlowFactor * size(movieIndices, 1) / t(end));
            if (trueFPS ~= nFramesPerSecond)
                fprintf('Unable to display at %d frames per second\n', nFramesPerSecond);
                fprintf('Try decreasing the time interval in t and x\n');
                fprintf('Display will be %d frames per second instead\n', trueFPS);
            end
        end
        
        function SimplexConnec = fnsSelectNewTriangulation(V)
            
            % This function guides the user to select valid simplices to completely
            % triangulate the polytope and returns the vertex indices for each
            % simplex selected. The order of vertices per simplex is sorted based on
            % the vertex indices, but the simplices themselves are not. Hence, two
            % triangulations may be equivalent up to the point of the order of the
            % simplices.
            % The triangulation is determined to be complete when the area of the
            % polytope is filled.
            
            X = V(:, 1);
            Y = V(:, 2);
            
            nVert = size(V, 1);
            
            %% Determine the clicking tolerance
            
            % Compute the pairwise distance between each vertex and find the minimum
            % distance.  Work with the square to avoid computing the square roots
            eps = 1e6;
            for i = 1:nVert - 1
                for j = i + 1:nVert
                    distsq = V(i, :)' - V(j, :)';
                    distsq = distsq(1)^2 + distsq(2)^2; % square of distance
                    if (distsq < eps)
                        eps = distsq;
                    end
                end
            end
            
            % The tolerance can be less than half the minimum distance, or less than
            % quarter of the square of the distance
            eps = eps / 4;
            
            %% Main variables to update
            
            % Area of max outer region - area of obstacle
            areaPolytope = polyarea(X([1:6, 1]), Y([1:6, 1])) - ...
                polyarea(X([7:12, 7]), Y([7:12, 7]));
            
            areaCovered = 0; % area covered so far
            
            nS = 0; % current number of simplicies created
            SimplexConnec = zeros(1, 3);
            
            %% Set up the plot
            
            scrsz = get(0,'ScreenSize');
            h = figure('Name', 'Select new triangulation', ...
                'OuterPosition', [50, 100, scrsz(3) * 0.8, scrsz(4) * 0.8]);
            
            % Fill the obstacle hole with a gray color
            xObs = [X(7:12)', X(7)];
            yObs = [Y(7:12)', Y(7)];
            graycolObs = 0.3;
            graycolSimp = 0.9;
            hold on, patch(xObs, yObs, graycolObs * ones(1, 3));
            
            % Plot the points individually to obtain handles to be able to switch
            % marker colors on and off as points are selected
            hVert = zeros(nVert, 1);
            for i = 1:nVert
                hold on, hVert(i) = plot(X(i), Y(i), 'o');
            end
            
            set(hVert, 'MarkerFaceColor', 'w');
            uistack(hVert, 'top'); % emphasize the vertices
            
            title(sprintf('Total number of simplices created: %d\n', nS));
            
            %% User point selection
            
            epsA = 1e-6; % Hopefully this is a reasonable value
            
            % Interact with the user until the necessary number of simplicies have been
            % created
            while(abs(areaCovered - areaPolytope) > epsA)
                nV = 0;
                chosenVertices = zeros(1, 3);
                while (nV < 3) % always must select three vertices
                    [x, y] = ginput(1);
                    
                    % Check which vertex this is closest to
                    distsqs = zeros(nVert, 1);
                    for i = 1:nVert
                        distsq = V(i, :)' - [x; y];
                        distsqs(i) = distsq(1)^2 + distsq(2)^2; % square of distance
                    end
                    idxPt = find(distsqs < eps);
                    
                    if isempty(idxPt)
                        strInput = sprintf('Total number of simplices created: %d\n', nS);
                        strInput = sprintf('%sPlease click on a valid point', strInput);
                        title(strInput);
                        continue; % invalid point, try again
                    elseif size(idxPt, 2) > 1
                        strInput = sprintf('Total number of simplices created: %d\n', nS);
                        strInput = sprintf('%sAmbiguous selection: click closer to a point', strInput);
                        title(strInput);
                        continue; % invalid point, try again
                    else
                        % Success message to user
                        title(sprintf('Total number of simplices created: %d\n', nS));
                    end
                    
                    % Add to the list of vertices for the current simplex only if the
                    % selection is valid (unique in the triplet)
                    if isempty(find(chosenVertices == idxPt, 1))
                        nV = nV + 1;
                        chosenVertices(nV) = idxPt;
                        set(hVert(idxPt), 'MarkerFaceColor', 'r');
                    else
                        strInput = sprintf('Total number of simplices created: %d\n', nS);
                        strInput = sprintf('%sPlease select unique vertices per simplex', strInput);
                        title(strInput);
                    end
                end
                
                % Sort the selection for convenience and consistency as order does not
                % matter
                chosenVertices = sort(chosenVertices);
                idxLoop = [chosenVertices, chosenVertices(1)];
                
                % First verify that the current simplex has non-zero area (i.e. it is
                % not a degenerate simplex - a line)
                areaNew = polyarea(X(idxLoop'), Y(idxLoop'));
                if areaNew < epsA
                    % User must reselect
                    strInput = sprintf('Total number of simplices created: %d\n', nS);
                    strInput = sprintf('%sPlease create a simplex with non-zero area', strInput);
                    title(strInput);
                else
                    % Then verify that the current simplex is unique
                    bMatch = false; % assume unique
                    for i = 1:nS
                        if isequal(SimplexConnec(nS, :), chosenVertices)
                            bMatch = true; % not unique, user must reselect
                            break;
                        end
                    end
                    
                    if bMatch
                        % User must reselect
                        strInput = sprintf('Total number of simplices created: %d\n', nS);
                        strInput = sprintf('%sPlease create a unique simplex', strInput);
                        title(strInput);
                    else
                        % Now verify that the current simplex is valid in the polytope
                        % Use intersection of polygons
                        % Specifically, compare the selected simplex overlap area with all the
                        % previous simplices and verify that there is no intersection
                        % Also ensure that the simplex does not overlap through the obstacle
                        % region
                        bIntersect = clCrane.fnsTriangulation_Intersections(X, Y, xObs, yObs, ...
                            chosenVertices, SimplexConnec, nS);
                        
                        if bIntersect
                            % User must reselect
                            strInput = sprintf('Total number of simplices created: %d\n', nS);
                            strInput = sprintf('%sInvalid selection: no simplices ', strInput);
                            strInput = sprintf('%smay overlap with each other or ', strInput);
                            strInput = sprintf('%sthe obstacle', strInput);
                            title(strInput);
                        else
                            % All tests passed and the simplex is valid
                            % Draw the simplex to assist the user visually in further selections
                            nS = nS + 1;
                            SimplexConnec(nS, :) = chosenVertices;
                            
                            hold on, patch(X(idxLoop'), Y(idxLoop'), graycolSimp * ones(1, 3));
                            uistack(hVert, 'top'); % emphasize the vertices
                            % Success message to user
                            title(sprintf('Total number of simplices created: %d\n', nS));
                            
                            % Update the amount of area of the polytope covered
                            areaCovered = areaCovered + areaNew;
                        end
                    end
                end
                
                % Clear the marker colors
                set(hVert, 'MarkerFaceColor', 'w');
                
            end
            
            close(h);
        end
        
        function bIntersect = fnsTriangulation_Intersections(X, Y, xObs, yObs, ...
                chosenVertices, SimplexConnec, nS)
            
            % This function determines whether the simplex specified by connectivity in
            % the input chosen vertices intersects with any of the previously selected
            % vertices.  Vertex indices are linked to the physical 2D points in [X, Y]
            % Also ensure that the simplex does not intersect through the obstacle,
            % described as a polygon with coordinates [xObs', yObs']
            % This function is an improvement to the previous function
            % 'Triangulation_Intersections'
            
            bIntersect = false;
            eps = 1e-6; % Hopefully this is a reasonable value
            
            % Set up the currently selected simplex polygon
            chosenCoord = [X([chosenVertices, chosenVertices(1)]), ...
                Y([chosenVertices, chosenVertices(1)])];
            
            % Set up the obstacle polygon
            obsCoord = [xObs', yObs'];
            
            % Check for area intersection
            if fnsPolygonAreaIntersection(chosenCoord, obsCoord) > eps
                bIntersect = true;
                return;
            end
            
            for s = 1:nS
                % Set up a previously selected simplex polygon
                simplexCoord = [X([SimplexConnec(s, :), SimplexConnec(s, 1)]), ...
                    Y([SimplexConnec(s, :), SimplexConnec(s, 1)])];
                
                % Check for area intersection
                if clCrane.fnsPolygonAreaIntersection(chosenCoord, simplexCoord) > eps
                    bIntersect = true;
                    return;
                end
            end
        end
        
        function area = fnsPolygonAreaIntersection(set1, set2)
            
            % This function finds the polygon describing the intersection region 
            % of the two polygons specified in set1 and set2
            % This function then returns the area of this intersecting polygon
            
            % Find all of the intersection points between all pairs of line 
            % segments on each polygon
            % This function does not exist in MATLAB 2011!
            [xi, yi] = polyxpoly(set1(:, 1), set1(:, 2), set2(:, 1), set2(:, 2));
            
            % Find all points of one polygon that lie in or on the other polygon
            in1 = inpolygon(set2(:, 1) , set2(:, 2), set1(:, 1), set1(:, 2));
            in2 = inpolygon(set1(:, 1) , set1(:, 2), set2(:, 1), set2(:, 2));
            
            xP = [xi; set2(in1, 1); set1(in2, 1)];
            yP = [yi; set2(in1, 2); set1(in2, 2)];
            
            if isempty(xP)
                area = 0;
            else
                try
                    c = convhull(xP, yP);
                    area = polyarea(xP(c), yP(c));
                catch
                    % Degenerate case, when not enough points specified or the points
                    % are collinear
                    area = 0;
                end
            end
        end
        
        function [Vleft, Vright] = fnsFixPolytope2(L, wo, hr, ho)
            
            dispFlag = 0;
            dispRes = 0;
            
            h = hr - ho;            % space between ceiling and obstacle, m
            ths = acos(h / L);      % (theta_star) angle such that bob meets the 
                                    %  obstacle tip, rad
            w = sqrt(L^2 - h^2);    % z length projected from pendulum when at 
                                    %  (theta_star) angle, m
            N = 25;
            
            x1_1 = linspace(wo/2 - w, wo/2, N);
            x3_1 = atan((wo/2 - x1_1)/h);
            
            x1_2 = linspace(wo/2, wo/2 + w, N);
            x3_2 = -asin((x1_2 - wo/2)/L);
            
            % Iteratively move (x1,x3) from (wo/2 + w, -ths) to (wo/2 + w + dz, -ths)
            % and vary the angle from the maximum until it intersects the nonlinear
            % curves beyond tangency. The maximum angle at the start position is
            % specified by the slope of -1/h, so all slopes will be tested more
            % positive than this. Tangency is tested by intersects of the test line
            % with the nonlinear curves (there are two of them). Ideally one stops when
            % the intersection is tangent to one or both of the nonlinear curves.
            % However approximating the nonlinear curves as polylines means a tangent
            % intersection will look like two very close intersections; hence an
            % epsilon is specified to determine which intersections are tangent or not.
            % Intersections are tested separately for the test line vs. the two
            % nonlinear curves
            
            dz = 0.01;
            z = wo/2 + w;
            zMax = wo + w;
            dm = 0.1;
            mMin = -1/h;
            mMax = -h/L^2;
            
            bestA = 1e9;
            bestz = z;
            bestm = mMin;
            
            if dispRes
                zs = [];
                counter = [];
                idx = 0;
            end
            
            while z < zMax
                if dispRes
                    idx = idx + 1;
                    zs = [zs; z];
                    counter = [counter; 0];
                end
                
                m = mMin;
                
                while m < mMax
                    if dispRes
                        counter(idx) = counter(idx) + 1;
                    end
                    
                    % Calculate other z position, based current z and m
                    z2 = 2*ths/m + z; 
                    
                    % Check for intersections
                    if clCrane.fnsNonLinPolytope_Intersection(x1_1, x3_1, ...
                            [z, z2], [-ths, ths])
                        break;
                    end
                    
                    if clCrane.fnsNonLinPolytope_Intersection(x1_2, x3_2, ...
                            [z, z2], [-ths, ths])
                        break;
                    end
                    
                    % If valid, check if this is the best area
                    xPoly = [x1_1, x1_2, z, z2];
                    yPoly = [x3_1, x3_2, -ths, ths];
                    A = polyarea(xPoly, yPoly);
                    if A < bestA
                        bestA = A;
                        bestz = z;
                        bestm = m;
                    end
                    
                    if dispFlag
                        h = figure;
                        hold all;
                        plot(x1_1, x3_1, 'b-', x1_2, x3_2, 'b-');
                        patch(xPoly, yPoly, 'c');
                        close(h);
                    end
                    
                    m = m + dm;
                end
                
                z = z + dz;
            end
            
            bestz2 = 2*ths/bestm + bestz;
            Vleft = [bestz2, ths];
            Vright = [bestz, -ths];
            
            if dispRes
                xPoly = [x1_1, x1_2, bestz, bestz2];
                yPoly = [x3_1, x3_2, -ths, ths];
                h = figure;
                hold all;
                plot(x1_1, x3_1, 'b-', x1_2, x3_2, 'b-');
                patch(xPoly, yPoly, 'c');
                figure;
                plot(zs, counter);
            end
        end
        
        function bIntersect = fnsNonLinPolytope_Intersection(x1c, x3c, lx, ly)
            
            % Helper function for fnsFixPolytope2()
            
            eps = 1e-6;
            bIntersect = 0;
            [xi, yi] = polyxpoly(x1c, x3c, lx, ly);
            nInt = size(xi, 1);
            
            switch(nInt)
                case 0
                    % no intersection, ok
                case 1
                    % exact tangency, ok
                    % But more often, it intersects non-tangentially
                    bIntersect = 1;
                case 2
                    dsquared = ((xi(1)-xi(2))^2 + (yi(1)-yi(2))^2);
                    if dsquared > eps
                        bIntersect = 1;
                    end
                otherwise
                    error('A line cannot intersection this curve more than twice?');
            end
            
        end
        
        function [Vleft, Vright] = fnsFixPolytope(pL, halfW)
            
            % Obtain a set of points representing the nonlinear relationship between
            % the distance from the wall (z_max) and the maximum pendulum swing allowed
            % before hitting the wall (theta_max)
            stepSz = 0.001;
            w_l =  halfW - pL;
            z_max = linspace(w_l + stepSz, halfW - stepSz, 100);
            theta_max = asin((halfW - z_max) / pL);
            
            % See Mario's notes. Find the line, tangent to the nonlinear curve found
            % above, such that the area between the curve and the line is minimal.
            % The area is found by constructing a polygon of the shape and using
            % polyarea()
            % An analytic approach is possible, but requires solving a non-linear
            % equation for the best slope. Once the slope is known, the point of
            % tangency (zp,thp) can be easily found using a derived equation
            
            ms = linspace(-pi / pL, -1 / pL, 100); % negative inf slope is unfeasible
            As = zeros(1, size(ms ,2));
            for i = 1:size(ms, 2)
                m = ms(i); % current slope
                zp = halfW - sqrt(pL^2 - 1 / m^2);      % tangent point
                thp = asin(sqrt(1 - 1 / (pL * m)^2));   % tangent point
                z0 = (pi / 2 - thp) / m + zp; % z point where line intersects 
                                              %  theta = pi/2 (max)
                zend = halfW - pL;
                fend2 = m * (halfW - zp) + thp;     % theta point where
                                                    %  the line intersects z = halfW
                poly = [z_max, halfW, halfW, z0, zend; theta_max, 0, ...
                    fend2, pi/2, pi/2];
                As(i) = polyarea(poly(1, :), poly(2, :));  % area
            end
            
            % Find the best area, and therefore best line
            [~, idx] = min(As);
            m_min = ms(idx);
            zp = halfW - sqrt(pL^2 - 1 / m_min^2);
            thp = asin(sqrt(1 - 1 / (pL * m_min)^2));
            z0 = (pi / 2 - thp) / m_min + zp;
            fend2 = m_min * (halfW - zp) + thp;
            
            % Take the two end points of the best line, which are returned 
            % as the replacement points
            Vleft = [z0, pi/2];
            Vright = [halfW, fend2];
        end
        
        function [K,g] = fnsFeedback(v,u)
            
            % Copyright (c) Graeme Ashford
            % feedback Version 1.0
            %
            % Version History
            % 2010-10-08 - v1.0 - original
            %
            % Function Description:
            % feedback takes as input the set of vertices, v, and the control inputs at
            % each vertex, u.  v is organized as follows:
            %
            %      |      |
            % v = [v0 ... vn]
            %      |      |
            %
            % where the v0,...,vn are n-dim column vectors.  u is organized as
            % follows:
            %
            %      |      |
            % u = [u0 ... un]
            %      |      |
            %
            % where the u0,...,un are m-dim column vectors representing the inputs at
            % v0,...,vn respectively.
            %
            % Usage Example:
            %
            % [K,g] = feedback(v,u);
            %
            n = size(v,1);
            V = [v;ones(1,n+1);];
            u = u/V;
            K = u(:,1:n);
            g = u(:,n+1);
            
        end
        
        function fnsDrawSimplices(S, V, bDetailed, hDraw)
            
            X = V(:, 1);
            Y = V(:, 2);
            
            % Fill the obstacle hole with a gray color
            xObs = [X(7:12)', X(7)];
            yObs = [Y(7:12)', Y(7)];
            graycol = 0.7;
            hold on, patch(xObs, yObs, graycol * ones(1, 3));
            
            % Loop through the simplices
            for i = 1:size(S, 3)
                xVertex = [S(1, :, i), S(1, 1, i)];
                yVertex = [S(2, :, i), S(2, 1, i)];
                hold on, plot(xVertex, yVertex, 'k-', 'LineWidth', 2);
            end
            
            % Plot a label in each simplex identifying the associated index value
            for i = 1:size(S, 3)
                hold on, text(mean(S(1, :, i)), mean(S(2, :, i)), ...
                    sprintf('$S_{%d}$', i), 'BackgroundColor', [0.4 0.8039 0.6667], ...
                    'FontSize', 12, 'Interpreter', 'Latex');
            end
            
            if bDetailed
                % Label vertex coordinates
                for i = 1:size(X, 1)
                    hold on, text(X(i), Y(i), sprintf('$v_{%d}$', i), ...
                        'BackgroundColor', [0.4980    1.0000    0.8314], ...
                        'FontSize', 12, 'Interpreter', 'Latex');
                end
                
                if bDetailed > 1
                    % In each simplex, label the order of the points as
                    % i, ii, iii
                    strLab = {'1', '2', '3'};
                    for i = 1:size(S, 3)
                        for j = 1:3
                            xPos = 0.5 * (mean(S(1, :, i)) + S(1, j, i));
                            yPos = 0.5 * (mean(S(2, :, i)) + S(2, j, i));
                            hold on, text(xPos, yPos, sprintf('%s', strLab{j}), ...
                                'BackgroundColor', [0.5961    0.9843    0.5961], ...
                                'FontSize', 12, 'Interpreter', 'Latex');
                        end
                    end
                end
            end
            
            xlabel('Cart Position (m)', 'Interpreter', 'Latex');
            ylabel('Crane Angle (rad)', 'Interpreter', 'Latex');
            title('State Space Polytope', 'FontSize', 12, 'Interpreter', 'Latex');
            axis([min(X) * 1.1, max(X) * 1.1, min(Y) * 1.1, max(Y) * 1.1]);
            
            if ~isempty(hDraw)
                scrSz = get(0,'ScreenSize');
                set(hDraw, 'units', 'normalized', ...
                    'outerposition', [0 0 scrSz(4)/scrSz(3) 1]);
                set(hDraw,'Units', 'Inches');
                pos = get(hDraw, 'Position');
                set(hDraw, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
                    'PaperSize', [pos(3), pos(4)]);              
                print(hDraw, '-dpdf', ['PlotDump/' 'Simplices']);
            end
            
        end
        
        function chosenIdx = fnsDisplayAllMeshes(fileIndices, triangFolder, V, ...
                meshFileRootName)
            
            % This function loops through all the meshes in the library and displays
            % them all as subplots. The user must then select one of the meshes.
            % Note: the plot may get unwieldly as the number of meshes increases
            
            nWindows = size(fileIndices, 1);
            AR = 3/2; % hard-coded aspect ratio of display, # columns/ # rows
            nRows = ceil(sqrt(nWindows / AR));
            nCols = ceil(AR * nRows);
            if (nRows * nCols < nWindows)
                error('Fix the subplot sizing code!!');
            end
            
            scrsz = get(0,'ScreenSize');
            figure('OuterPosition', [150, 150, scrsz(3) * 0.8, scrsz(4) * 0.8]);
            for k = 1:nWindows
                fileLoad = sprintf('%s/%s%d.mat', triangFolder, meshFileRootName, ...
                    fileIndices(k));
                load(fileLoad, 'SimplexConnec');
                nS = size(SimplexConnec, 1);
                S = zeros(2, 3, nS);
                for i = 1:nS
                    c = SimplexConnec(i, :);
                    S(:, :, i) = [V(c(1), :); V(c(2), :); V(c(3), :)]';
                end
                subplot(nRows, nCols, k);
                clCrane.fnsDrawSimplices(S, V, 0, []);
                title(sprintf('%s%d', meshFileRootName, fileIndices(k)));
            end
            
            % Prompt the user to select a mesh (by its index)
            while (1)
                chosenIdx = input('\nPlease select an existing mesh (input an integer).\n');
                if ~isempty(find(fileIndices == chosenIdx))
                    break;
                end
            end
        end
        
        function fnsSearchUnstableEigs(lambda)
            % Evokes an error if all the eigenvalues (in lambda) are stable
            
            eps = 1e-6;
            
            bUnstable = false;
            for i = 1:4
                if real(lambda(i) >= -eps)
                    bUnstable = true;
                end
            end
            if (bUnstable == false)
                error('It should be impossible to be stable!');
            end
        end
        
        function fnsNonLin_ComplexConstraints(hr, wr, ho, wo, L)
            
            %% Plotting script for the derived non-linear constraints
            % Specifies the ranges of theta as non-linear piece-wise functions of z
            % Let x1 = z, x3 = theta
            
            %% Useful derived quantities, constants
            
            h = hr - ho;            % space between ceiling and obstacle, m
            ths = acos(h / L);      % (theta_star) angle such that bob 
                                    %  meets the obstacle tip, rad
            w = sqrt(L^2 - h^2);    % z length projected from pendulum 
                                    %  when at (theta_star) angle, m
            
            % Check for problem well-posedness
            bPass = 0;
            if w^2 > 0 && w^2 < 0.25 * (wr - wo)^2
                bPass = 1;
            end
            
            %% Set up the plot
            
            hDraw = figure;
            
            if bPass
                strVal = 'Well-posed';
            else
                strVal = 'Ill-posed';
            end
            
            title(sprintf('Polytope\n%s', strVal));
            hold all;
            grid on;
            N = 50;
            linW = 2;
            
            %% Plot theta_max
            
            x1 = linspace(-wr/2, wr/2 - L, N);
            x3 = pi/2 * ones(1, N);
            plot(x1, x3, 'b-', 'LineWidth', linW+1);
            
            x1 = linspace(wr/2 - L, wr/2, N);
            x3 = asin((wr/2 - x1)/L);
            plot(x1, x3, 'b-', 'LineWidth', linW+1);
            
            plot([wr/2, wr/2], [0, -pi/2], 'b-', 'LineWidth', linW+1);
            
            
            %% Plot theta_min
            
            plot([-wr/2, -wr/2], [pi/2, 0], 'r-', 'LineWidth', linW+1);
            
            x1 = linspace(-wr/2, -wr/2 + L, N);
            x3 = -asin((x1 + wr/2)/L);
            plot(x1, x3, 'r-', 'LineWidth', linW+1);
            
            x1 = linspace(-wr/2 + L, wr/2, N);
            x3 = -pi/2 * ones(1, N);
            plot(x1, x3, 'r-', 'LineWidth', linW+1);
            
            %% Plot theta_lower
            
            x1 = linspace(-wo/2 - w, -wo/2, N);
            x3 = asin(-(wo/2 + x1)/L);
            plot(x1, x3, 'r-', 'LineWidth', linW);
            
            x1 = linspace(-wo/2, w - wo/2, N);
            x3 = -atan((x1 + wo/2)/h);
            plot(x1, x3, 'r-', 'LineWidth', linW);
            
            x1 = linspace(w - wo/2, wo/2 + w, N);
            x3 = -ths * ones(1, N);
            plot(x1, x3, 'r-', 'LineWidth', linW);
            
            
            
            %% Plot theta_upper
            
            x1 = linspace(-wo/2 - w, wo/2 - w, N);
            x3 = ths * ones(1, N);
            plot(x1, x3, 'b-', 'LineWidth', linW);
            
            x1 = linspace(wo/2 - w, wo/2, N);
            x3 = atan((wo/2 - x1)/h);
            plot(x1, x3, 'b-', 'LineWidth', linW);
            
            x1 = linspace(wo/2, wo/2 + w, N);
            x3 = -asin((x1 - wo/2)/L);
            plot(x1, x3, 'b-', 'LineWidth', linW);
            
            
            
            %% Plot lines that cover up the curves
            
%             [Vleft, Vright] = clCrane.fnsFixPolytope(L, wr/2);
%             plot([Vleft(1), Vright(1)], [Vleft(2), Vright(2)], 'k--', ...
%                 'LineWidth', 2); % optimized theta max fix
%             plot(-[Vleft(1), Vright(1)], -[Vleft(2), Vright(2)], 'k--', ...
%                 'LineWidth', 2); % optimized theta min fix
%             
%             [Vleft2, Vright2] = clCrane.fnsFixPolytope2(L, wo, hr, ho);
%             plot([Vleft2(1), Vright2(1)], [Vleft2(2), Vright2(2)], 'k--', ...
%                 'LineWidth', 2); % optimized theta upper fix
%             plot(-[Vleft2(1), Vright2(1)], -[Vleft2(2), Vright2(2)], 'k--', ...
%                 'LineWidth', 2); % optimized theta upper fix
            
%             % Infeasible theta upper fix
%             plot([wo/2 - w, wo/2 + w], [ths, -ths], 'k-.');
%             
%             % Simple theta min fix
%             plot([-wr/2, L*pi/2 - wr/2], [0, -pi/2], 'k--', 'LineWidth', linW); 
%             
%             % Simple theta max fix
%             plot([wr/2, -L*pi/2 + wr/2], [0, pi/2], 'k--', 'LineWidth', linW); 
%             
%             % Simple theta upper fix
%             plot([wo/2 + w, wo/2 + w - 2*ths*h], [-ths, ths], 'k--', 'LineWidth', linW); 
%             
%             % Simple theta lower fix
%             plot([-(wo/2 + w), -(wo/2 + w) + 2*ths*h], [ths, -ths], 'k--', 'LineWidth', linW); 
%             
%             % Other simple theta upper fix direction
%             plot([wo/2 - w, wo/2 - w + 1], [ths, -(h/L^2)*1 + ths], 'm--'); 
            
%             % Bad version
%             plot([wr/2 - L, wr/2], [pi/2, 0], ...
%                 'k--', 'LineWidth', linW);
%             plot([-wr/2 + L, -wr/2], [-pi/2, 0], ...
%                 'k--', 'LineWidth', linW);
%             plot([-wo/2 - w, w - wo/2], [ths, -ths], ...
%                 'k--', 'LineWidth', linW);
%             plot([wo/2 - w, wo/2 + w], [ths, -ths], ...
%                 'k--', 'LineWidth', linW);
            
            
            xlabel('z');
            ylabel('\theta');
            
            % Save
            set(hDraw,'Units', 'Inches');
            pos = get(hDraw, 'Position');
            set(hDraw, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
                'PaperSize', [pos(3), pos(4)]);
            print(hDraw, '-dpdf', ['PlotDump/' 'NonLinPolytopeWFixes']);
            
%             close(hDraw);
%             hDraw = figure;
%             hold all;
%             plot([wr/2 - L, wr/2], [pi/2, 0], ...
%                 'k-', 'LineWidth', linW);
%             plot([-wr/2 + L, -wr/2], [-pi/2, 0], ...
%                 'k-', 'LineWidth', linW);
%             plot([-wo/2 - w, w - wo/2], [ths, -ths], ...
%                 'k-', 'LineWidth', linW);
%             plot([wo/2 - w, wo/2 + w], [ths, -ths], ...
%                 'k-', 'LineWidth', linW);
%             plot([-wr/2, -wr/2], [0, pi/2], 'k-', 'LineWidth', linW);
%             plot([-wr/2, wr/2 - L], [pi/2, pi/2], 'k-', 'LineWidth', linW);
%             plot([wr/2, wr/2], [0, -pi/2], 'k-', 'LineWidth', linW);
%             plot([wr/2, -wr/2 + L], [-pi/2, -pi/2], 'k-', 'LineWidth', linW);
%             plot([-wo/2 - w, wo/2 - w], [ths, ths], 'k-', 'LineWidth', linW);
%             plot([w - wo/2, wo/2 + w], [-ths, -ths], 'k-', 'LineWidth', linW);
%             
%             xlabel('z');
%             ylabel('\theta');
%             axis off
%             
%             set(hDraw,'Units', 'Inches');
%             pos = get(hDraw, 'Position');
%             set(hDraw, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
%                 'PaperSize', [pos(3), pos(4)]);
%             print(hDraw, '-dpdf', ['PlotDump/' 'GenericPolytope']);
            
         
        end
        
        function fnsSaveSchematic(hr, wr, ho, wo, L, ...
                drawPos_z, drawPos_th, labelFlag, labelFlag2)
            
            % Create circles and a square, with arbitrary scaling that look aesthetic:
            theta = 0 : (pi / 100) : (2 * pi);
            s = (hr / 10) * [-1.5, -1.5, 1.5, 1.5, -1.5; 0, 2, 2, 0, 0];
            xscirc = (hr / 20) * sin(theta); % for the pendulum bob
            yscirc = (hr / 20) * cos(theta);
            wheelrad = hr / 28;
            xwcirc = wheelrad * sin(theta); % for the cart wheel
            ywcirc = wheelrad * cos(theta);
            xssquare = s(1, :);
            yssquare = s(2, :);
            
            hDraw = figure;
            hold all;
            
            % Create the walls
            minz = - wr / 2;
            maxz = wr / 2;
            miny = - hr;
            maxy = 0;
            line([minz, minz], [maxy, miny], 'Color', 'k', 'LineWidth', 2);
            line([maxz, maxz], [maxy, miny], 'Color', 'k', 'LineWidth', 2), hold on;
            line([minz, maxz], [miny, miny], 'Color', 'k', 'LineWidth', 2), hold on;
            line([minz, maxz], [maxy, maxy], 'Color', 'k', 'LineWidth', 2), hold on;
            
            % Create the obstacle
            % vertical space between ceiling and obstacle
            space = -(hr - ho);
            obsHalfW = wo / 2;
            patch(obsHalfW * [-1, 1, 1, -1, -1], [miny miny space space miny], ...
                0.7 * ones(1, 3));
            
            % Typically for use in the report
            % z and theta contain a selected value good for display of the pendulum
            % configuration
            
            axis equal
            axis off
            %set(hDraw, 'visible', 'off');
            %set(hDraw, 'outerposition', [0 0 1 1]);
            
            %% Positions
            
            tracePos = [drawPos_z + L * sin(drawPos_th), -L * cos(drawPos_th)];
            bobPos = [drawPos_z(end) + L * sin(drawPos_th(end)), -L * cos(drawPos_th(end))];
            cartPosX = drawPos_z(end);
            % Y position of the bottom of the cart is always zero
            
            %% Draw and label main features
            
            % Arm
            ha = line([cartPosX(1, 1), bobPos(1, 1)], [0, bobPos(1, 2)], ...
                'Color', 'r', 'LineWidth', 2);
            uistack(ha, 'bottom');
            if labelFlag
                rad = norm([0.5 * abs(cartPosX(1, 1) - bobPos(1, 1)), ...
                    0.5 * bobPos(1, 2)]);
                if drawPos_th(end) < 0
                    dth = -10*(pi/180);
                else
                    dth = 10*(pi/180);
                end
                ang = drawPos_th(end) + dth;
                pos = [cos(ang), -sin(ang); sin(ang), cos(ang)] * [0; -rad];
                text('Interpreter', 'Latex', ...
                    'Position', [pos(1) + cartPosX(1, 1), pos(2)],...
                    'HorizontalAlignment', 'center', ...
                    'String', '$L$');
            end
            
            if labelFlag
                % Zero angle reference line
                line([cartPosX(1, 1), cartPosX(1, 1)], [0, bobPos(1, 2)], ...
                    'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
                text('Interpreter', 'Latex', ...
                    'Position', [cartPosX(1, 1) + tan(drawPos_th(end)) * ...
                    mean([0, abs(bobPos(1, 2))])/2, mean([0, bobPos(1, 2)])], ...
                    'HorizontalAlignment', 'center', ...
                    'String', '$\theta$');
                % Curved line
                rad = L*cos(drawPos_th(end))/3;
                ang = linspace(0, drawPos_th(end), 30);
                curve = zeros(2, size(ang, 2));
                for i = 1:size(ang,2)
                    curve(:, i) = [cos(ang(i)), -sin(ang(i)); sin(ang(i)), ...
                        cos(ang(i))] * [0; -rad];
                    curve(1, i) = curve(1, i) + cartPosX(1,1);
                end
                plot(curve(1,:), curve(2,:), 'Color', 'k', 'LineWidth', 1, ...
                    'LineStyle', '-');
            end
            
            if labelFlag2
                % Zero angle reference line
                line([cartPosX(1, 1), cartPosX(1, 1)], [0, bobPos(1, 2)], ...
                    'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
                text('Interpreter', 'Latex', ...
                    'Position', [cartPosX(1, 1) + tan(drawPos_th(end)) * ...
                    mean([0, abs(bobPos(1, 2))])/2, mean([0, bobPos(1, 2)])], ...
                    'HorizontalAlignment', 'center', ...
                    'String', '$\theta^{\star}$');
                % Curved line
                rad = L*cos(drawPos_th(end))/3;
                ang = linspace(0, drawPos_th(end), 30);
                curve = zeros(2, size(ang, 2));
                for i = 1:size(ang,2)
                    curve(:, i) = [cos(ang(i)), -sin(ang(i)); sin(ang(i)), ...
                        cos(ang(i))] * [0; -rad];
                    curve(1, i) = curve(1, i) + cartPosX(1,1);
                end
                addpath([docroot '/techdoc/creating_plots/examples']);
                plot(curve(1,:), curve(2,:), 'Color', 'k', 'LineWidth', 1, ...
                    'LineStyle', '-');
                
                vPos = -(hr-ho) - 0.05;
                D = sqrt(L^2 - (hr-ho)^2);
                if drawPos_z < 0
                    [arrowx,arrowy] = dsxy2figxy(gca, [-D-wo/2, -wo/2], vPos*[1,1]);
                    avgPos = mean([-D, -wo/2]);
                else
                    [arrowx,arrowy] = dsxy2figxy(gca, [D+wo/2, wo/2], vPos*[1,1]);
                    avgPos = mean([D, wo/2]);
                end
                annotation('doublearrow', arrowx, arrowy);
                text('Interpreter', 'Latex', ...
                    'Position', [avgPos, vPos - 0.05], 'String', 'd');
            end
            
            % Cart
            patch(xssquare + cartPosX(1, 1), yssquare, 'g');
            if labelFlag
                text('Interpreter', 'Latex', ...
                    'Position', [cartPosX(1, 1), mean([min(yssquare),max(yssquare)])], ...
                    'HorizontalAlignment', 'center', ...
                    'String', '$M$');
            end
            
            % Trace
            % Only draw if enough data points are given
            plot(tracePos(:,1), tracePos(:,2), 'b-');
            
            % Bob
            % If multiple data points given, the bob position at the end is
            % drawn
            hb = patch(xscirc + bobPos(1, 1), yscirc + bobPos(1, 2), 'b');
            uistack(hb, 'top');
            if labelFlag
                rad = L;
                if drawPos_th(end) < 0
                    dth = -5*(pi/180);
                else
                    dth = 5*(pi/180);
                end
                ang = drawPos_th(end) + dth;
                pos = [cos(ang), -sin(ang); sin(ang), cos(ang)] * [0; -rad];
                text('Interpreter', 'Latex', ...
                    'Position', [pos(1) + cartPosX(1, 1), pos(2)], ...
                    'HorizontalAlignment', 'left', ...
                    'String', '$m$');
                %'Position', [bobPos(1, 1) + (max(xscirc)-min(xscirc)), bobPos(1, 2)], ...
            end
            
            % Bonus: create the wheels of the cart
            wheelShift = (max(xssquare) - min(xssquare)) / 2 - wheelrad;
            patch(xwcirc + cartPosX(1, 1) - wheelShift, ...
                ywcirc - wheelrad, 'k');
            patch(xwcirc + cartPosX(1, 1) + wheelShift, ...
                ywcirc - wheelrad, 'k');
            
            if labelFlag
                % Label the z axis
                line([0, 0], [min(yssquare), 1.5 * max(yssquare)], ...
                    'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
                text('Interpreter', 'Latex', ...
                    'Position', [0, mean([min(yssquare),max(yssquare)])], ...
                    'HorizontalAlignment', 'left', ...
                    'String', '$\rightarrow z$');
                
                % Label control input
                text('Interpreter', 'Latex', ...
                    'Position', [cartPosX(1, 1) - (max(xssquare)-min(xssquare)), ...
                    mean([min(yssquare),max(yssquare)])], ...
                    'HorizontalAlignment', 'right', ...
                    'String', '$u \rightarrow$');
            end
            
            %% Label room dimensions
            
            if labelFlag
                addpath([docroot '/techdoc/creating_plots/examples']);
                
                lineoffset = 0.15*hr;
                
                % Room width
                roomWx = [-wr/2, wr/2];
                roomWy = -(hr+lineoffset)*[1, 1];
                % Transform from data space to figure space
                [arrowx,arrowy] = dsxy2figxy(gca, roomWx, roomWy);
                annotation('doublearrow', arrowx, arrowy);
                text('Interpreter', 'Latex', ...
                    'Position', [0, -(hr+2*lineoffset)], ...
                    'HorizontalAlignment', 'center', ...
                    'String', '$w_{room}$', 'BackgroundColor', 'w');
                
                % Room height
                roomHx = wr/2+lineoffset * [1, 1];
                roomHy = [0, -hr];
                % Transform from data space to figure space
                [arrowx,arrowy] = dsxy2figxy(gca, roomHx, roomHy);
                annotation('doublearrow', arrowx, arrowy);
                text('Interpreter', 'Latex', ...
                    'Position', [wr/2+1.6*lineoffset, -hr/2],...
                    'HorizontalAlignment', 'left', ...
                    'String', '$h_{room}$', 'BackgroundColor', 'w');
                
                % Obstacle width
                obsWx = [-wo/2, wo/2];
                obsWy = (-(hr-ho)+lineoffset)*[1, 1];
                % Transform from data space to figure space
                [arrowx,arrowy] = dsxy2figxy(gca, obsWx, obsWy);
                annotation('doublearrow', arrowx, arrowy);
                text('Interpreter', 'Latex', ...
                    'Position', [0, -(hr-ho)+1.8*lineoffset],...
                    'HorizontalAlignment', 'center', ...
                    'String', '$w_{obs}$', 'BackgroundColor', 'w');
                
                % Obstacle height
                obsHx = wo/2+lineoffset * [1, 1];
                obsHy = [-(hr-ho), -hr];
                % Transform from data space to figure space
                [arrowx,arrowy] = dsxy2figxy(gca, obsHx, obsHy);
                annotation('doublearrow', arrowx, arrowy);
                text('Interpreter', 'Latex', ...
                    'Position', [wo/2+1.5*lineoffset, -ho/2-(hr-ho)],...
                    'HorizontalAlignment', 'left', ...
                    'String', '$h_{obs}$', 'BackgroundColor', 'w');
                
            end
            
            %% Save and close the drawing
            
            zStr = sprintf('%.2f', drawPos_z(end));
            zStr = strrep(zStr, '.', 'dot');
            thStr = sprintf('%.2f', drawPos_th(end));
            thStr = strrep(thStr, '.', 'dot');
            saveName = sprintf('Crane_z=%s_th=%s', zStr, thStr);
            
            set(hDraw,'Units', 'Inches');
            pos = get(hDraw, 'Position');
            pos2 = axis;
            set(hDraw, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
                'PaperSize', [1.2 * pos(3) * ((maxz - minz)/(pos2(2) - pos2(1))), ...
                1.6 * pos(4) * ((maxy - miny)/(pos2(4) - pos2(3)))]);
            
            print(hDraw, '-dpdf', ['PlotDump/' saveName]);
            %close(hDraw);
            
        end
        
        function strNum = fnsTwoDigitPaddedStr(num)
            
            if num < 10
                strNum = sprintf('0%d', num);
            else
                strNum = sprintf('%d', num);
            end
        end
        
    end
end