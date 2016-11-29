function result = general_direct_collocation(duration, targetangle, Ncolloc, oldresult)
% function result = general_direct_collocation(duration, targetangle, Ncolloc, oldresult)
% Jeremy Wong 2016-11-29
% implements via Hargraves 1987 Spline collocation constraint.
% Thereby replaces transcription method of van den Bogert.
% Allows black-box uses of optimization by replacing 'getstated' with arbitrary
% ode function returning ddt(state) from state. 
% There are 3 solution types:
% 1 Newton (as per the pdf 'DirectCollocation_simple.pdf')
% (github.com/cato469/DC)
% 2 IPOPT
% 3 fmincon
% Fmincon has been used to check the derivatives of the constraint
% equations. So there is a long derivativecheck at the beginning.
% This example:
% Finds optimal motion of a torque-driven pendulum.  Task is to move from one
% static posture to another in a given time.
% Contributing Author: Ton van den Bogert <a.vandenbogert@csuohio.edu>
% This work is licensed under a Creative Commons Attribution 3.0 Unported License.
% http://creativecommons.org/licenses/by/3.0/deed.en_US
% Inputs:
%	duration		duration of the movement (s)
%	targetangle		target angle (deg)
%	Ncolloc				number of collocation nodes to use
% 	oldresult		(optional) initial guess
% Notes:
% 1. This code may be useful as a template for solving other optimal control problems, such
%    as cart-pole upswing.
% 2. IPOPT will be used when it is installed, otherwise Newton's method.  IPOPT is recommended
%    because it is more rowbust.  Newton's method can still solve most problems, but you may
%    need to solve a sequence of problems of increasing difficulty to ensure convergence.
% 3. The solution may be a local optimum, especially for tasks that involve multiple
%    revolutions.  Try different initial guesses (maybe random) to check for this.

% The following examples all converge with IPOPT or Newton:
%	r = pend(1.0, 180, 100);		% swing up in 1 second, 100 collocation nodes
%   r = pend(3.0, 180, 100);		% now do it in 3 seconds, note the countermovement
%   r = pend(10.0, 180, 100);		% now do it in 10 seconds, multiple countermovements are seen
%   r = pend(5.0, 720, 300);		% do two full revolutions in 5 seconds, 300 collocation nodes
%   r2 = pend(..,..,..,r);			% use previous result r as initial guess

% settings
MaxIterations = 300;
if exist('ipopt') == 3
    method = 'ipopt';
elseif 1
    disp('IPOPT is not installed.');
    disp('Newton method will be used, may be more sensitive to initial guess.');
    disp('Hit ENTER to continue...');
    pause
    method = 'newton';
else
    method = 'fmincon';
end

% initializations
close all
tic
h = duration/(Ncolloc-1);		% time interval between nodes
times = h*(0:Ncolloc-1)';		% list of time points
Nstates = 2;
Nfdiff = 1;
Nconstraints = Nstates*(Ncolloc-Nfdiff)+ 4;			% Ncolloc-2 dynamics constraints and 4 task constraints
Nq = Nstates*Ncolloc;
Ncontrols = 1;
% model parameters
L = 1;			%	length of pendulum (m)
m = 1;			%	mass of pendulum (kg)
I = 1;			% 	moment of inertia relative to pivot (kg m^2)
g = 9.81;		%	gravity (m s^-2)

% state variable is x: angle relative to hanging down
% control variable is u: torque applied at joint

% if oldresult was provided, use it as initial guess, otherwise use zero initial guess (pendulum hangs down, no torque)
if (nargin == 4)
    %     oldN = numel(oldresult.t);
    %     oldreltime = (0:oldN-1)'/(oldN-1);			% sample times of old result, from 0 to 1
    %     newreltime = (0:Ncolloc-1)'/(Ncolloc-1);				% sample times for new optimization, from 0 to 1
    %     x = interp1(oldreltime, oldresult.x, newreltime);
    %     u = interp1(oldreltime, oldresult.u, newreltime);
    x = oldresult.x;
    u = oldresult.u;
else
    %    x = randn(Nq,1);
    %    u = randn(Ncolloc*Ncontrols,1);
    x = zeros(Nq,1);
    u = zeros(Ncolloc,1);
end

% encode initial guess of unknowns into a long column vector X
X0 = [x ; u];
ix = (1:Nq);				% index to elements in X where angles x are stored
iphi = ix(1:Ncolloc);
Ndyncon=Nstates * (Ncolloc-Nfdiff);
iu = Nq + (1:Ncolloc*Ncontrols);			% index to elements in X where controls u are stored
NX = size(X0,1);		% number of unknowns
show(X0, confun(X0));	% show the initial guess
drawnow

if nargin ==4
    [init_con,a,b,c,extra] = confun(X0);
    init_conjac = conjac(X0);
    figure;plot(init_con);
    figure;
    plot(extra.eom,'r');hold on;
    plot(extra.finite_difference,'g');
    acc = finite_differences(X0(101:200),h);
    legend({ 'eom','finite-difference'});
    fprintf('checking\n');
end;
% X0 = X0 + 0.001*randn(size(X0));		% perturb the initial guess a little before optimizing

if strcmp(method,'ipopt')
    % solve the NLP with IPOPT
    funcs.objective = @objfun;
    funcs.gradient  = @objgrad;
    funcs.constraints = @confun;
    funcs.jacobian    = @conjac;
    funcs.jacobianstructure = @conjacstructure;
    options.cl = zeros(Nconstraints,1);
    options.cu = zeros(Nconstraints,1);
    options.ipopt.max_iter = MaxIterations;
    options.ipopt.hessian_approximation = 'limited-memory';
    [X, info] = ipopt(X0,funcs,options);
elseif strcmp(method, 'newton')
    % solve the NLP using Newton iteration on the KKT conditions
    X = X0;
    ctol = 1e-2;		% constraint tolerance
    ftol = 1e-4;		% cost function tolerance
    xtol = 1e-4;		% solution tole rance
    F = 1e10;
    
    
    
    for iter=1:MaxIterations
        Fprev = F;
        
        % evaluate objective function F and constraint violations c
        F = objfun(X);
        G = objgrad(X);
        H = objhess(X);
        c = confun(X);
        J = conjac(X);
        
        % form the linearized KKT system K*x = b
        K = [H J'; J sparse(Nconstraints,Nconstraints)];
        b = [-G; -c];
        
        % solve the linear system K*dZ=b
        % Z is a vector containing the unknowns X and the Lagrange multipliers.  dZ is the change in this iteration
        dZ = K\b;
        dX = dZ(1:NX);					% the first NX are the elements of X
        
        % do a half Newton step (converges slower than full Newton step, but more likely to converge)
        % for more robust convergence, we should do a line search here to make sure we always have progress
        X = X + dX/2;
        rmsC = sqrt(mean(c.^2));
        rmsdX = sqrt(mean(dX.^2));
        fprintf('Iter: %3d  F=%10.5e  rms(c)=%10.5e   rms(dX)=%10.5e\ncolloc', iter,F,rmsC,rmsdX);
        
        if (max(abs(c)) < ctol) && (abs(F-Fprev)<ftol) && (mean(abs(dX))<xtol)
            break;
        end
    end
    if iter >= MaxIterations
        disp('Maximum number of iterations exceeded.');
    else
        disp('Optimal solution found');
    end
elseif strcmp(method,'fmincon')
    % Bounds on the optimization parameters
    q_LB(1:Ncolloc) = -180*pi/180;      q_UB(1:Ncolloc) = 180*pi/180;
    qdot_LB(1:Ncolloc) = -180*pi/180*1000;
    qdot_UB(1:Ncolloc) = 180*pi/180*1000;
    % act_LB(1:N) = 0.001;      act_UB(1:Ncoord*N) = 1;
    % lm_LB(1:N) = 0;      lm_UB(1:Ncoord*N) = 0.5;
    
    u_LB(1:Ncolloc) = -5;    u_UB(1:Ncolloc) = 5;
    lb = [q_LB qdot_LB u_LB]';
    ub = [q_UB qdot_UB u_UB]';
    % lb = [q_LB qdot_LB act_LB lm_LB u_LB]';
    % ub = [q_UB qdot_UB act_UB lm_UB u_UB]';
    func.confun = @(x)confun(x);
    func.objfun = @(x)objfun(x);
    X=X0;
    % Set some parameters for fmincon and then run the optimization
    % options = optimoptions('algorithm','interior-point','TolFun',1e-4,'TolX',1e-4, ...
    %                    'TolCon',1e-4,'FinDiffType','central','MaxFunEvals',1e5, ...
    %                    'Hessian','bfgs','display','iter','DerivativeCheck','on',...
    %                    'SpecifyObjectiveGradient','true','SpecifyConstraintGradient',...
     a=ver;
    if isequal(a(1).Release,'(R2015a)')
         options = optimoptions(@fmincon,'Algorithm','interior-point','DerivativeCheck','on','TolFun',1e-4,'TolX',1e-4, ...
        'TolCon',1e-4,'FinDiffType','central','MaxFunEvals',1e5, ...
        'Hessian','bfgs','display','iter-detailed',...
        'FunValCheck','on','Diagnostics','on','GradConstr','on');
    elseif isequal(a(1).Release,'(R2016b')
        options = optimoptions(@fmincon,'Algorithm','interior-point','CheckGradients',true,'TolFun',1e-4,'TolX',1e-4, ...
        'TolCon',1e-4,'FinDiffType','central','MaxFunEvals',1e5, ...
        'Hessian','bfgs','display','iter-detailed',...
        'FunValCheck','on','Diagnostics','on','GradConstr','on');
    end;   
    
    % start a timer
    tic;
    [Xopt,fval,exitflag,output] = fmincon(func.objfun,X0,[],[],[],[],lb,ub,func.confun,options);
    
    
else
    error('method not recognized');
end

% plot results
show(X, confun(X));

% make movie of the solution
disp('Hit ENTER to generate animation...');
pause
fps = Ncolloc/duration;
avi = VideoWriter('pend.avi');
open(avi);
figure(2);
clf;
set(gcf,'Position',[5 100 650 650]);
set(gcf, 'color', 'white');
s = 1.5*L;
for i_vid=1:Ncolloc
    plot([-s s],[0 0],'k','LineWidth',2);
    hold on
    plot([0 L*cos(X(i_vid)-pi/2)], [0 L*sin(X(i_vid)-pi/2)],'b-o','LineWidth',2);
    axis('equal');
    axis('square');
    axis([-s s -s s]);
    title(['t = ' num2str(times(i_vid),'%8.3f')]);
    dovid = 0;
    if dovid
        if (i_vid==1)
            F = getframe(gca);
            frame = [1 1 size(F.cdata,2) size(F.cdata,1)];
        else
            F = getframe(gca,frame);
        end
        writeVideo(avi,F);end;
    drawnow;
    hold off;
end
close(avi);
%	close(2);

% store results
result.t = times;
result.x = X(1:Nq);
result.u = X(Nq+1:end);
[c,ceq,J,gradceq,extra] = confun(X)
fprintf('done collocation.\n');

% start of embedded functions

%=========================================================
    function F = objfun(X)
        % objective function: integral of squared controls
        F = h * sum(X(iu).^2);
    end

%=========================================================
    function G = objgrad(X)
        % gradient of the objective function coded in objfun
        G = zeros(NX,1);
        G(iu) = 2 * h * X(iu);
    end

%=========================================================
    function H = objhess(X)
        % hessian of objective function coded in objfun
        H = spalloc(NX,NX,Ncolloc);
        H(iu,iu) = 2 * h * speye(Ncolloc,Ncolloc);
    end

%=========================================================
    function out = getstated(x,u)
        %this should be replaced by the call to opensim.
        out(1)=x(2);
        phi = x(1);
        out(2) = ( -m * g * L*sin(phi) + u) / I;
        out = out(:);
    end
%=========================================================
%DC_f
% this is a direct transcription of the collocation for constraint
% equations. 
    function out = DC_f(fn,xp0,xp1,up0,up1)
        uphalf = (up0+up1)/2;
        dstate0 = fn(xp0,up0);
        dstate1 = fn(xp1,up1);
        out = getstated(1/2*(xp0+xp1)+h/8*(dstate0-dstate1),uphalf) + 3/(2*h)*(xp0-xp1) + 1/4*(dstate0+dstate1);
    end
%=========================================================
    function [c,ceq,J,gradceq,extra] = confun(X)
        ceq = []; gradceq=[];
        J = conjac(X);
        J = J';
        % constraint function (dynamics constraints and task constraints)
        
        % size of constraint vector
        c = zeros(Ndyncon,1);
        
        % dynamics constraints
        % Note: torques at node 1 and node Ncolloc do not affect movement and will therefore
        % always be zero in a minimal-effort solution.
        [xp0,xp1,xp2,dxdt]=deal(zeros(Nstates,1));
        extra = struct;
        for i=1:Ncolloc-Nfdiff
            for j = 1:Nstates
                xp0(j) = X(i+(j-1)*(Ncolloc));
                xp1(j) = X(i+1+(j-1)*(Ncolloc));
                xp2(j) = X(i+2+(j-1)*(Ncolloc));
                dxdt(j) = (xp2(j)-xp0(j))/(2*h);
                if j ==1
                    extra.fd2(i+j-1) = (xp2(1) - 2*xp1(1) + xp0(1))/h^2;							% three-point formula for angular acceleration

                end;
            end;
            up1=zeros(Ncontrols,1);
            for l=1:Ncontrols
                up1(l) = X(i + 1+ Nq+ (l-1)*(Ncolloc));
            end
            up0=zeros(Ncontrols,1);
            for l=1:Ncontrols
                up0(l) = X(i + 0+ Nq+ (l-1)*(Ncolloc));
            end
            stated = getstated(xp1,up1);
            temp_dc1=getstated(1/2*(xp0+xp1) + h/8*(getstated(xp0,up0) - getstated(xp1,up1)),(up0+up1)/2);
            temp_dc2 = 3/(2*h)*(xp0-xp1) + 1/4*(getstated(xp0,up0) + getstated(xp1,up1));
            
            for j =1:Nstates
%                 c(i+(j-1)*(Ncolloc-Nfdiff)) =  dxdt(j) - stated(j);		% equation of motion must be satisfied
                extra.finite_difference(i+(j-1)*(Ncolloc-Nfdiff)) =  dxdt(j);
                extra.eom(i+(j-1)*(Ncolloc-Nfdiff)) =  stated(j);
                extra.dc1(i+(j-1)*(Ncolloc-Nfdiff)) = temp_dc1(j);
                extra.dc2(i+(j-1)*(Ncolloc-Nfdiff)) = temp_dc2(j);
                c(i+(j-1)*(Ncolloc-Nfdiff)) = temp_dc1(j) + temp_dc2(j);
            end; 
        end
        
        % task constraints
        % initial position must be zero:
        c(Ndyncon+1) 	= X(1);
        % initial velocity must be zero:
        c(Ndyncon+2) 	= X(2)-X(1);
        %c(Ndyncon+2) 	= X(Ncolloc+1);
        % final position must be at target angle:
        c(Ndyncon+3) 	= X(Ncolloc) - targetangle*pi/180;
        % final velocity must be zero: HERE I AM USING POSITION DIFFERENCE.
        % STAY WITH VANDENBOGERT.
        c(Ndyncon+4) 	= X(Ncolloc)-X(Ncolloc-1);
        %c(Ndyncon+4) 	= X(Nq);
        
        % show current iterate, every 0.1 second
        %		if toc > 0.1
        show(X,c);
        tic;
        %		end
        
    end

%=========================================================
    function J = conjac(X)
        
        % size of Jacobian
        J = spalloc(Nconstraints,NX,4*Nconstraints + 6);
        
        % dynamics constraints
        for i=1:Ncolloc-Nfdiff
            % Jacobian matrix: derivatives of c(i) with respect to the elements of X
            %xp1 is the one that matters. since it is in the EOM.
            states0 = zeros(Nstates,1);
            states1 = zeros(Nstates,1);
            for j = 1:Nstates
                states0(j)=X(i+0+(j-1)*(Ncolloc));
                states1(j)=X(i+1+(j-1)*(Ncolloc));
            end;
            states0 = states0(:);
            states1 = states1(:);
            
            %get all controls
            cmds0 = zeros(Ncontrols,1);
            cmds1 = zeros(Ncontrols,1);
            for l = 1:Ncontrols
                cmds0(l) = X(i + Nq + 0 + (l-1)*Ncolloc);
                cmds1(l) = X(i + Nq + 1 + (l-1)*Ncolloc);
            end;
            
            %get stated
            dstate0 = getstated(states0,cmds0);
            dstate1 = getstated(states1,cmds1);
            
            fn=@getstated;
            deltapert = 1e-6;
            for j = 1:Nstates
                i_con = i+(j-1)*(Ncolloc-Nfdiff);
                C = DC_f(fn,states0,states1,cmds0,cmds1);
                for k = 1:Nstates % perturb states.
                    pert     = zeros(Nstates,1);
                    pert(k)  = deltapert;
                    
                    d_state0 = DC_f(fn,states0 + pert,states1,cmds0,cmds1);
                    dC_dstate0 = (d_state0(j)-C(j))/deltapert;
                    
                    d_state1 = DC_f(fn,states0,states1 + pert,cmds0,cmds1);
                    dC_dstate1 = (d_state1(j)-C(j))/deltapert;

                    J(i_con, i + 0 + (k-1)  *Ncolloc) = ...
                        dC_dstate0;
                    J(i_con, i + 1 + (k-1)  *Ncolloc) = ...
                        dC_dstate1;

                end;
                for l = 1:Ncontrols % perturb controls.
                    pert = zeros(Ncontrols,1);
                    pert(l)  = deltapert;
                    d_up0 = DC_f(fn,states0,states1,cmds0+pert,cmds1);
                    J(i_con, i + 0 + Nq +(l-1)*Ncolloc) = ...
                        (d_up0(j)-C(j))/deltapert;
                    
                    d_up1 = DC_f(fn,states0,states1,cmds0,cmds1+pert);
                    J(i_con, i + 1 + Nq +(l-1)*Ncolloc) = ...
                        (d_up1(j)-C(j))/deltapert;
                    
                end;
            end;
        end
        
        % task constraints        
        % initial position must be zero:
        J(Ndyncon+1,	1) = 1;
        % initial velocity must be zero:
        J(Ndyncon+2,2) = 1;
        J(Ndyncon+2,1) = -1;
        %J(Ndyncon+2,Ncolloc+1) = 1;
        % final position must be at target angle:
        J(Ndyncon+3,Ncolloc) = 1;
        % final velocity must be zero:
        J(Ndyncon+4,Ncolloc) = 1;
        J(Ndyncon+4,Ncolloc-1) = -1;
        %J(Ndyncon+4,Nq) = 1;
    end
%=========================================================
    % this function 'conjacstructure' defines the sparsity of the problem.
    % used in Ipopt.
    function J = conjacstructure(X)
        
        % size of Jacobian
        J = spalloc(Nconstraints,NX,4*Nconstraints + 6);
        
        % dynamics constraints
        for i=1:Ncolloc-Nfdiff
            % Jacobian matrix: derivatives of c(i) with respect to the elements of X
            
            for j = 1:Nstates
                ind_eq =  i+(j-1)*(Ncolloc-Nfdiff);
                for k = 1:Nstates
                    J(ind_eq, i + 0 + (k-1)*Ncolloc) = 1;
                    J(ind_eq, i + 1 + (k-1)*Ncolloc) = 1;
                    
                end;
                for l = 1:Ncontrols
                    J(ind_eq, i + Nq + 0 +(l-1)*Ncolloc) 	= 1;
                    J(ind_eq, i + Nq + 1 +(l-1)*Ncolloc) 	= 1;

                end;
            end;
        end
        
        % task constraints
        
        % initial position must be zero:
        J(Ndyncon+1,	1) = 1;
        % initial velocity must be zero:
        %J(Ndyncon+2,Ncolloc+1) = 1;
        J(Ndyncon+2,2) = 1;
        J(Ndyncon+2,1) = 1;
        
        % final position must be at target angle:
        J(Ndyncon+3,Ncolloc) = 1;
        % final velocity must be zero:
        % J(Ndyncon+4,Nq) = 1;
        J(Ndyncon+4,Ncolloc) = 1;
        J(Ndyncon+4,Ncolloc-1) = 1;
        
    end
%============================================================
    function show(X,c)
        % plot the current solution
        x = X(iphi);
        u = X(iu);
        figure(1)
        subplot(3,1,1);plot(times,x*180/pi);title('angle')
        subplot(3,1,2);plot(times,u);title('torque');
        subplot(3,1,3);plot(c);title('constraint violations');
    end
end
