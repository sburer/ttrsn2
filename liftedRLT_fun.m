function [rm_old,rg,rm_new,cuts,rv] = liftedRLT_fun(filename)

load(filename);

tol = 1.0e-14;
% damp = 0.75;
damp = 1;
num_trials = 10^5;
if damp < 1.0
  warning('Using a dampening factor.');
end
maxiter = 10;
myset = sdpsettings('verbose',0);

SOCRLT_or_liftedRLT_or_both = 1; % 1, 2, or 3
  
  n = size(Q,1);
  rm_old = log10(eigY(n+1)/eigY(n));
  da = zeros(n,1);

  %% Derive the sqrt of H

  [V,D] = eig(H); d = diag(D); Hsqrt = V*diag(sqrt(d))*V';

  %% Setup first SDP

  x = sdpvar(n,1);
  X = sdpvar(n,n,'symm');
  Y = [1,x';x,X];

  con = [ Y >= 0 ];
  con = [ con ; trace(X) <= r1^2 ];
  con = [ con ; H(:)'*X(:) - 2*h'*H*x + h'*H*h <= r2^2 ];

  con_SOCRLT    = [];
  con_liftedRLT = [];

  obj = Q(:)'*X(:) + c'*x;

  %% Setup separation problem (objective TBD)

  a = sdpvar(n,1);
  A = sdpvar(n,n,'symm');
  sepcon = [ [1,a';a,A] >= 0 ];
  sepcon = [ sepcon ; trace(A) <= 1 ];

  %% Iterate

  iter = 1;
  stop = 0;
  while stop == 0 & iter <= maxiter

    %% Solve current relaxation

    if SOCRLT_or_liftedRLT_or_both == 1
      solvesdp([con;con_SOCRLT],obj,myset);
    elseif SOCRLT_or_liftedRLT_or_both == 2
      solvesdp([con;con_liftedRLT],obj,myset);
    else
      solvesdp([con;con_SOCRLT;con_liftedRLT],obj,myset);
    end

    %fprintf('iter = %d   relaxed value = %f\n', iter, double(obj));
    rv = double(obj);

    %% Store portion of solution

    dx = double(x);
    dX = double(X);
    eigY = eig(double(Y));

    %% Setup objective for SOC-RLT separation

    tmp1 = r1*(dx-h);
    tmp2 = h*dx' - dX;

    sepobj1 = (r2^2)*(r1^2 - 2*r1*a'*dx + dx'*A*dx);

    sepobj21 = tmp1'*H*tmp1;
    sepobj22 = tmp1'*H*tmp2;
    sepobj23 = tmp2'*H*tmp2;

    sepobj = sepobj1 - (sepobj21 + 2*sepobj22*a + sepobj23(:)'*A(:));

    %% Solve separation subproblem

    solvesdp(sepcon,sepobj,myset);

    %% If no separation, quit

    if double(sepobj)/(norm(Q,'fro') + norm(c)) > tol 
      %fprintf('Terminating because no SOCRLT active or violated.\n');
      stop = 1;
    else

      %% Specify new SOCRLT contraint

      da = double(a);
      new_SOCRLT = [ norm( Hsqrt*(r1*(x-h) - X*da + h*da'*x )) <= r2*(r1 - da'*x) ];

      %% If we want liftedRLT constraints...

      if SOCRLT_or_liftedRLT_or_both > 1

        %% Add new_SOCRLT to latest relaxation and solve

        if SOCRLT_or_liftedRLT_or_both == 1
          solvesdp([con;con_SOCRLT;new_SOCRLT],obj,myset);
        elseif SOCRLT_or_liftedRLT_or_both == 2
          solvesdp([con;con_liftedRLT;new_SOCRLT],obj,myset);
        else
          solvesdp([con;con_SOCRLT;con_liftedRLT;new_SOCRLT],obj,myset);
        end

        % Extract y and z

        y = r1*da;
        y = r1*y/norm(y);
        z = (r1*dx - dX*da)/(r1 - da'*dx);
        z = r2*z/norm(Hsqrt*z);
        aa = y/(r1*norm(y));
        bb = H*z/r2^2;
%         norm(y)
%         norm(z)
%         y'*H*y
%         z'*H*z
%         y'*H*z
        1 - aa'*y; % Should be 0
        1 - bb'*z; % Should be 0

        % Calculate projection matrix

        dx = double(x);
        dX = double(X);
        P = sdpvar(n);
        newcon = [P >= 0 ; trace(P) == n-1 ; P*(y-z) == 0];
        newobj = P(:)'*dX(:) - 2*y'*P*dx + y'*P*y;
        solvesdp(newcon,-newobj,myset);
        P = double(P);

        % Extract positive eigenvector

        [V,D] = eig(P);
        d = diag(D);
        [~,ind] = max(d);
        if numel(ind) > 1
          eig(P)
          error('Unexpected rank(P) > 1');
        end
        v = V(:,ind);

        % Calculate lambda
        
        %lambda = calc_lambda_rand(r1,Hsqrt,r2,y,aa,bb,v,num_trials);
        lambda = lambda_generalhd(eye(n)/r1^2,zeros(n,1),H/r2^2,h,y,z,v);

        % Save new liftedRLT constraint

        con_liftedRLT = [con_liftedRLT;...
           1 - aa'*x - bb'*x + aa'*X*bb - damp*lambda*(v'*X*v - 2*y'*v*v'*x + (y'*v)^2) >= 0];

      end

      %% Add new_SOCRLT to con_SOCRLT
      con_SOCRLT = [ con_SOCRLT ; new_SOCRLT ];

    end
    
    iter = iter + 1;
  end

  dx = double(x);
  dX = double(X);
  eigY = eig(double(Y));
  rg = ((dx'*Q*dx + c'*dx) - double(obj))/abs(dx'*Q*dx + c'*dx);
  rm_new = log10(eigY(n+1)/eigY(n));
  cuts = iter - 2;
