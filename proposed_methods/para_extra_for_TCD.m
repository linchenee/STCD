function [Range,Angle,Veloc] = para_extra_for_TCD(A,B,C,P,L,c0,lambda,K,delta_f,Ts,Mt)
    % ---------------------------------------------------------------------------    
    % The function is used to extract target parameters for the UTCD/STCD method.
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ---------------------------------------------------------------------------
    Range = zeros(1,L);
    Angle = zeros(1,L);
    Veloc = zeros(1,L);
    for i = 1 : L
        Angle(i) = asind((1/(2*pi)) * 2 * PUMA(A(:, i)));
        temp = PUMA(B(:, i));
        if temp > 0
            Range(i) = (c0/(2*delta_f)) * ((2*pi - temp)/(2*pi));
        else
            Range(i) = (c0/(2*delta_f)) * (-temp/(2*pi));
        end

        %% The polynomial method for velocity estimation
        Ati = exp(1j * pi * (0 : (Mt - 1)).' * sin(Angle(i) * pi / 180));
        Pequal = diag(P.'*Ati);
        W = conj(Pequal)*(norm(C(:,i),'fro')^2*eye(K)-C(:,i)*C(:,i)')*Pequal.';
        coe = zeros(1, 2*K-1); % Computing coefficients for polynomial construction
        for j = -(K-1):K-1
            coe(-j+K) = sum(diag(W,j));
        end
        a1 = roots(coe);
        b1 = a1(abs(a1)<1);
        [~,I1] = sort(abs(abs(b1)-1));
        Veloc(1,i) = (0.25/(pi*Ts))*angle(b1(I1(1)))*lambda;
    end
end