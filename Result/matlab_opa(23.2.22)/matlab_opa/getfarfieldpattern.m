function au = getfarfieldpattern(phi0,d_phi,phase_error,lambda,z,range,d) % z: distance to far field. range: simulation range on far field
   % d : width of antenna grid
    [m,n] = size(phi0);
   % n: number of antennas
    
    k = 2*pi/lambda;
    for i = 1:m
        phi= phi0(i,:) + d_phi + phase_error;
        x = z*tan(range(i));
        r = (z^2+x.^2).^0.5;
        u(i) = sum(exp(1i*phi))*z*exp(1i*k*r)/(1i*lambda*r.^2)*d;
    end
    
    au = abs(u)/n;
end