function [u,w,tUsed]=computeNWTVelocity(simRes,xOut,zOut,tOut)

switch simRes.nwtSpec.solver.type
    case {'melhos','melhos_kInt'}
        [u,w,tUsed]=computeMELHOSVelocity(simRes.phiLS,simRes.x,simRes.z,simRes.eta,...
            simRes.alpha,simRes.t,simRes.nwtSpec.sim.depth,xOut,zOut,tOut,simRes.nwtSpec.solver.M);
    otherwise
        error('a:a','Cannot compute velocity field for solver type %s\n',simRes.nwtSpec.solver.type)
end

end