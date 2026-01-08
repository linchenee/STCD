function [CRB1,CRB2,CRB3,Fisher] = cal_CRB(Angle,Delay,Doppl,alpha,N1,N2,N3,sigma_2,mask,delta_f,Ts,pre,L,c0,lambda)
% ------------------------------------------
% The function is used to calculate CRB.
% version 1.0 - 01/07/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ------------------------------------------
A = exp(1j*pi*(0:(N1-1)).'*sin(Angle*pi/180));
B = exp(-1j*2*pi*delta_f*(0:N2-1).'*Delay);
D = exp(1j*2*pi*Ts*(0:N3-1).'*Doppl);
C = (pre.'*A).*D*diag(alpha);

A_tilde = 1j*pi*diag(0:(N1-1))*A*diag(cos(Angle*pi/180));
B_tilde = -1j*2*pi*delta_f*diag(0:(N2-1))*B;
C_tilde = 1j*2*pi*Ts*diag(0:(N3-1))*C;

num = 1;
FIM = cell.empty;
for mode_1 = 1:3
    for mode_2 = 1:3
        if mode_1 == 1
            temp1 = kron(A_tilde.',kr(B,C).');
        elseif mode_1 == 2
            temp1 = kron(B_tilde.',kr(C,A).');    
        elseif mode_1 == 3
            temp1 = kron(C_tilde.',kr(A,B).'); 
        end
        
        if mode_2 == 1
            temp2 = kron(conj(A_tilde),conj(kr(B,C)));
        elseif mode_2 == 2
            temp2 = kron(conj(B_tilde),conj(kr(C,A)));    
        elseif mode_2 == 3
            temp2 = kron(conj(C_tilde),conj(kr(A,B)));
        end
        
        C_N = Noise_Covariance(N1,N2,N3,mode_1,mode_2,sigma_2,mask);
        FIM_temp1 = 2*real((temp1*C_N*temp2)/(sigma_2^2));
        FIM{num} = convert_covari_vec2mat(FIM_temp1,L);
        num = num+1;
    end
end

Fisher = [FIM{1},FIM{2},FIM{3}; ...
         FIM{4},FIM{5},FIM{6};...
         FIM{7},FIM{8},FIM{9}];

invF = inv(Fisher);
invF = diag(invF);
crb1 = real(sum(invF(1:L))); % angle
crb2 = real(sum(c0^2*invF((L+1):2*L)/4)); % delay -> range
crb3 = real(sum(lambda^2*invF((2*L+1):3*L)/4)); % doppler -> velocity
CRB1 = sqrt(crb1/L);
CRB2 = sqrt(crb2/L);
CRB3 = sqrt(crb3/L);

function [C_N] = Noise_Covariance(M,T,K,mode_1,mode_2,sigma_2,mask)
    % The cross expectation matrix of Circular Symmetric Gaussian
    C_N = sparse(M*T*K,M*T*K);
    for m = 1:M
        for t = 1:T
            for k = 1:K
               if mask(m,t,k) == 1 
                if mode_1 == 1
                    index_1 = k+(t-1)*K+(m-1)*T*K;
                elseif mode_1 == 2
                    index_1 = m+(k-1)*M+(t-1)*M*K;
                elseif mode_1 == 3
                    index_1 = t+(m-1)*T+(k-1)*M*T;
                end
                
                if mode_2 == 1
                    index_2 = k+(t-1)*K+(m-1)*T*K;
                elseif mode_2 == 2
                    index_2 = m+(k-1)*M+(t-1)*M*K;
                elseif mode_2 == 3
                    index_2 = t+(m-1)*T+(k-1)*M*T;             
                end
                    
                C_N(index_1,index_2) = sigma_2;
               end
            end
        end
    end
end

function B = convert_covari_vec2mat(A,L)
    % Convert the L^2*L^2 matrix A into L*L matrix B
    B = zeros(L,L);
    for l1 = 1:L
        for l2 = 1:L
            B(l1,l2) = A(L*(l1-1)+l1,L*(l2-1)+l2);
        end
    end
end

end