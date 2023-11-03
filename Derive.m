function uPrime = Derive(u,uPrime,v)
% v is the viscous term, v = 0 is inviscid Burgers

    N = length(u);

    for k = -floor(N/2): floor(N/2) - 1
        
        if k < 0
            neg = 1;
        else
            neg = 0;
        end
        % Diffusion term with viscosity: V*u_xx
        uPrime(N*neg+k+1) = v*((1i*k)^2)*u(N*neg+k+1);
       
        for j = -floor(N/2): floor(N/2) - 1
          
            if j < 0
                jneg = 1;
            else 
                jneg = 0;
            end
            l = k - j;
            if l < 0
                lneg = 1;
            else
                lneg = 0;
            end

            if -floor(N/2) <= l && l <= floor(N/2) - 1
                % Convolution term for u*u_x 
                uPrime(N*neg+k+1) = uPrime(N*neg+k+1) - (1i*l)*u(N*jneg+j+1)*u(N*lneg+l+1);
            end
        end

    end
    
    

    end

    