function v = L2norm(u)
    Globals2D;
    v = sqrt(u(:)' * (MassMatrix*(J.*u))(:));
end
