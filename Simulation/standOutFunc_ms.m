%% Output Function: COM X Position
function comX = standOutFunc_ms(x,u)
global bodyParam

m    = bodyParam.m1+bodyParam.m2;
phi  = x(6);  % Floor angle

comLegX   = -bodyParam.h1*sin(x(2)+phi);
comTrunkX = -bodyParam.L1*sin(x(2)+phi)-bodyParam.h2*sin(x(4)+x(2)+phi);
comX = (bodyParam.m1*comLegX + bodyParam.m2*comTrunkX)/m;

end

