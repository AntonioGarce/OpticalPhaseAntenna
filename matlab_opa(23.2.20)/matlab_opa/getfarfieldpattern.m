function au = getfarfieldpattern(phi0,d_phi,phase_error,lambda,z,range,d) % z: distance to far field. range: simulation range on far field
   % d : width of antenna grid
    [m,n] = size(phi0);
   % n: number of antennas
    
    k = 2*pi/lambda;
    for l = 1:m
        phi= phi0(l,:) + d_phi + phase_error;
        x = z*tan(range(l));
        r = (z^2+x.^2).^0.5;
        u(l) = sum(exp(1i*phi))*z*exp(1i*k*r)/(1i*lambda*r.^2)*d;
    end
    
    au = abs(u)/n;
end