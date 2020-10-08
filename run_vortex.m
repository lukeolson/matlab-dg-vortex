Globals2D;
N = 1;
reflevel = 0;

allnorms = zeros(8,3)
for reflevel = 0:2
    for N = 1:8
        disp(['------------- N = ', num2str(N)])
        [Q, QEX, norms] = vortex(N, reflevel);
        allnorms(N, reflevel+1) = norms(1)
    end
end

%disp(['||rho|| = ', num2str(L2norm(rho-rhoEX))]);
% disp(['||rho u|| = ', num2str(L2norm(rhou-rhouEX))]);
