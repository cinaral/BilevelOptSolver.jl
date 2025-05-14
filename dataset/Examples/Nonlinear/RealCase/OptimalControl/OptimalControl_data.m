function data = OptimalControl_data(rate)

geometry = @squareg;
[mesh.p, mesh.e, mesh.t] = initmesh(geometry, 'hmax', rate);

% pdeplot(mesh.p, mesh.e, mesh.t)


% Extract inner points
n  = size(mesh.p, 2);
ip = setdiff(1:n, mesh.e(1:2,:));
ni = length(ip);

% Set up stiffness and mass matrix
[data.K_full, data.M_full, ~] = assema( mesh.p, mesh.t, 1, 1, 0 );

% Set up lumped mass matrix
Ml_full = spdiags(sum(data.M_full)', 0, n, n);

% Reduce to interior points
K       = data.K_full(ip, ip);
data.W  = data.M_full(ip, ip);
Ml      = Ml_full(ip, ip);

% Set up operators for lower-level problem

data.U = Ml;
data.C = speye(size(data.W));
data.A = K;
data.B = Ml;

% Set up bounds
data.ua = -10 * ones(ni,1);
data.ub = 10 * ones(ni,1);
  
% Define the actions P and Q
p1 = @(z)(10 * exp( -((z(1,:)-.7).^2 + (z(2,:)-.3).^2) * 5  ));
p2 = @(z)(10 * exp( -((z(1,:)+.4).^2 + (z(2,:)-.5).^2) * 10  ));


data.P  = [p1(mesh.p(:,ip))', p2(mesh.p(:,ip))'];
data.Q  = sparse(ni, 2);
data.Pt = data.P';
data.Qt = data.Q';
data.Ct = data.C';

%upper level objective set up
data.target = .2*p1(mesh.p)' + .3*p2(mesh.p)';
 
data.k      = [-0.1 ; -0.3 ];
data.ip     = ip;
data.ni     = ni;
data.n      = n;
data.sigma  = 1e-2;

data.fxx    = data.Pt*data.W*data.P + data.sigma*data.Qt*data.U*data.Q;  
data.fxy    = -[data.Pt*data.W*data.C data.sigma*data.Qt*data.U]'; 
data.fyy    = [data.Ct*data.W*data.C zeros(ni);zeros(ni) data.sigma*data.U]; 
            
end

