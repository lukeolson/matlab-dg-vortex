% Driver script for solving the 2D Euler equations
Globals2D;

% Order of polynomials used for approximation
N = 1;

% Read in Mesh
filename = 'vortexA04.neu';
InitialSolution = @IsentropicVortexIC2D;
ExactSolution   = @IsentropicVortexIC2D;
BCSolution      = @IsentropicVortexBC2D;

% read mesh from file
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D(filename);

% set up nodes and basic operations
StartUp2D;
disp(Cyl)

% turn cylinders into walls
ids = find(BCType==Cyl);
BCType(ids) = Wall;

refnum = 0
for i = 1:refnum,
    disp(K)
    refineflag = ones(K,1);
    Hrefine2D(refineflag);
    StartUp2D;
    disp(K)
end;

BuildBCMaps2D;

% compute initial condition
Q = feval(InitialSolution, x, y, 0);

% Solve Problem
FinalTime = 1;
[Q] = Euler2D(Q, FinalTime, BCSolution);

gamma = 1.4;
time = FinalTime;

rho = Q(:,:,1); rhou = Q(:,:,2); rhov = Q(:,:,3); Ener = Q(:,:,4);
u = rhou./rho;
v = rhov./rho;
p = (gamma-1)*(Ener-0.5*(rhou.*u + rhov.*v));
c = sqrt(gamma*p./rho);
magu = sqrt(u.^2 + v.^2);
mach = max(max(magu./c));

QEX = feval(ExactSolution, x, y, time);
rhoEX = QEX(:,:,1); rhouEX = QEX(:,:,2); rhovEX = QEX(:,:,3); EnerEX = QEX(:,:,4);

L2norm = @(u) sqrt(u(:)' * (MassMatrix*(J.*u))(:));
disp(['||rho|| = ', num2str(L2norm(rho-rhoEX))]);
disp(['||rho u|| = ', num2str(L2norm(rhou-rhouEX))]);

%disp(['||x|| = ', num2str(L2norm(x))]);
