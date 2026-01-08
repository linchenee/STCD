function [out1,out2,out3] = STCD(Y,L,omega,I,J,xi,eta,thr,A1,B1,C1)
% -------------------------------------------------------------    
% Structured tensor completion and decomposition (STCD) method
% version 1.0 - 01/06/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -------------------------------------------------------------
option.tol = 1e-3;
[I1,I2,I3] = size(Y);
m1 = 2^nextpow2(I1);
m2 = 2^nextpow2(I2);
rank1 = 1;
rank2 = 1;
T1 = construct_T(I1,rank1);
T2 = construct_T(I2,rank2);   
T11 = construct_T(I1,L);
T22 = construct_T(I2,L);

phi = 1e100;
A = A1;
B = B1;
C = C1;
S = permute(reshape(kr(A,B)*C.',I2,I1,I3),[2,1,3]);

for iter = 1:I
    phi_old = phi;
    S1_T = reshape(permute(S,[3,2,1]),I2*I3,I1).'; 
    S2_T = reshape(permute(S,[1,3,2]),I3*I1,I2).';                   
    S3_T = reshape(permute(S,[2,1,3]),I1*I2,I3).'; 

    %% Update the factor matrix A
    G1 = kr(B,C).';
    temp1 = G1*G1';
    temp2 = S1_T*G1';
    %% The first iteration in PGD, i.e., j=1
    A = PGD(S1_T*pinv(G1),T11,T1,m1,rank1,eta,option,temp1,temp2);
    fitA_old = norm(S1_T-A*G1,'fro');
    %% Form the second iteration to the Jth iteration in PGD
    if J > 1
        A1 = A;
        for j = 2:J
            A1 = PGD(A1,T11,T1,m1,rank1,eta,option,temp1,temp2);
            fitA = norm(S1_T-A1*G1,'fro');
            if fitA >= thr*fitA_old
                break;
            else
                A = A1;
                fitA_old = fitA;
            end
        end
        % fprintf('t=%d,fitA=%f\n', t, fitA);
    end

    %% Update the factor matrix B
    G2 = kr(C,A).';
    temp1 = G2*G2';
    temp2 = S2_T*G2';
    B = PGD(S2_T*pinv(G2),T22,T2,m2,rank2,eta,option,temp1,temp2);
    fitB_old = norm(S2_T-B*G2,'fro');
    if J > 1
        B1 = B;
        for j = 2:J
            B1 = PGD(B1,T22,T2,m2,rank2,eta,option,temp1,temp2);
            fitB = norm(S2_T-B1*G2,'fro');
            if fitB >= thr*fitB_old
                break;
            else
                B = B1;
                fitB_old = fitB;
            end
        end
        % fprintf('t=%d,fitB=%f\n', t, fitB);
    end
    
    %% Update the factor matrix C
    G3 = kr(A,B).';
    C = S3_T*pinv(G3); 
    phi = norm(S3_T-C*G3,'fro');    
    if (abs((phi-phi_old)/phi_old)<xi)
     break;
    end

    %% update S
    S = permute(reshape(kr(A,B)*C.',I2,I1,I3),[2,1,3]); % mode-3 matrix -> tensor
    S(omega) = Y(omega); 
end
out1 = A; out2 = B; out3 = C;
end


