function [X,A,B,C] = generate_echo(Angle,Delay,Doppl,Gain,W,num_tar,N,K,delta_f,Ts,Mr,Mt)
    % ------------------------------------------------------------    
    % The function is used to generate the echo tensor in Eq. (6).
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ------------------------------------------------------------
    A = exp(1j * pi * (0 : (Mr - 1)).' * sin(Angle.' * pi / 180));
    B = exp(-1j * 2 * pi * delta_f * (0 : N - 1).' * Delay);
    C = zeros(K, num_tar);
    At = exp(1j * pi * (0 : (Mt - 1)).' * sin(Angle.' * pi / 180));
    for i = 1 : num_tar
        D = diag(exp(1j*2*pi*(0:K-1)*Doppl(i)*Ts));
        C(:, i) = Gain(i) * D * W.'*At(:,i);
    end
    X = permute(reshape(kr(A,B)*C.',N,Mr,K),[2,1,3]);
end
