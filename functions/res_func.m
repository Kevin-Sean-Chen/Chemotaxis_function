function [K] = res_func(THETA, t_vec, type)
%%%
% Response function with parameter THETA and time vector t_vec
% Specify the response type: 
            % exponential: A*exp(-t_vec/tau)  (with amplitude and time scale)
            % alpha function:  a^2*t_vec*exp(-a*t_vec)
            % respone function: a*exp(-a*t_vec)*( (a*t_vec)^5/factorial(5) - (a*t_vec)^7/factorial(7))
            % Basis: A*Basis (with Basis being the basis functions)
% Return the response kernel K with the same length as tau
%%%

    if strcmp(type, 'exp')
        Amp = THETA(1);
        tau = THETA(2);
        K = Amp*exp(-t_vec/tau);
        
    elseif strcmp(type, 'alpha')
        a = THETA;
        K = a^2*t_vec.*exp(-a*t_vec);
        
    elseif strcmp(type, 'rf')
        a = THETA;
        K = a*exp(-a*t_vec).*( (a.*t_vec).^5/factorial(5) - (a*t_vec).^7/factorial(7));
        
    elseif strcmp(type, 'basis')
        nB = length(THETA);
        [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);
        K = THETA * cosBasis';
        
    end

end