function [Range,Angle,Veloc] = para_extra_for_VCPD_IMDF(wx,wy,Y,L,P,c0,lambda,delta_f,Ts,Mt,Mr,N1,K1)
    % ----------------------------------------------------------------------------  
    % The function is used to extract target parameters for the VCPD/IMDF method.
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ----------------------------------------------------------------------------
    Range = zeros(1,L);
    Angle = zeros(1,L);
    Veloc = zeros(1,L);
    [I1,I2,I3] = size(Y);
    Y3_T = reshape(permute(Y,[2,1,3]),I1*I2,I3).';

    for i = 1 : L
        Angle(1,i) = asind((1/(2*pi)) * 2 * wx(i));
        if wy(i) > 0
            Range(1,i) = (c0/(2*delta_f)) * ((2*pi - wy(i))/(2*pi));
        else
            Range(1,i) = (c0/(2*delta_f)) * (-wy(i)/(2*pi));
        end
    end

    A = exp(1j * pi * (0 : (Mr - 1)).' * sin(Angle * pi / 180));
    B = exp(-1j * 2 * pi * delta_f * (0 : N1 - 1).' * (2 * Range / c0));
    C = Y3_T*pinv(kr(A,B).'); 
    At = exp(1j * pi * (0 : (Mt - 1)).' * sin(Angle * pi / 180));

    %% The polynomial method for velocity estimation
    for i = 1 : L
        Pequal = diag(P.'*At(:,i));
        W = conj(Pequal)*(norm(C(:,i),'fro')^2*eye(K1)-C(:,i)*C(:,i)')*Pequal.';
        coe = zeros(1, 2*K1-1); % Computing coefficients for polynomial construction
        for j = -(K1-1):K1-1
            coe(-j+K1) = sum(diag(W,j));
        end
        a1 = roots(coe);
        b1 = a1(abs(a1)<1);
        [~,I1] = sort(abs(abs(b1)-1));
        Veloc(1,i) = (0.25/(pi*Ts))*angle(b1(I1(1)))*lambda;
    end
end