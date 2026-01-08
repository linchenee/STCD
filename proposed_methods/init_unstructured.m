function [out1,out2,out3] = init_unstructured(Y,R,omega,Tol)
% ------------------------------------------------------------------------------------    
% The function is used for unstructured initialization of the UTCD and STCD-1 methods.
% version 1.0 - 01/06/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ------------------------------------------------------------------------------------
[I1,I2,I3] = size(Y);
phiprevious = realmax;
nRestart = 5;
for iRe = 1:nRestart
    % rng(iRe);
    A = randn(I1,R);
    B = randn(I2,R);
    C = randn(I3,R);
    temp4 = kr(A,B);
    S = permute(reshape(temp4*C.',I2,I1,I3),[2,1,3]); % mode-3 matrix -> tensor
    phi = 1e100;

    temp3 = B'*B;
    for iter = 1:5
        phi_old = phi;
        S1 = reshape(permute(S,[3,2,1]),I2*I3,I1); 
        S2 = reshape(permute(S,[1,3,2]),I3*I1,I2);                   
        S3 = reshape(permute(S,[2,1,3]),I1*I2,I3); 

        temp1 = C'*C;
        A = ((temp3.*temp1)\ (kr(B,C)'*S1) ).';
        temp2 = A'*A;
        B = ((temp1.*temp2)\ (kr(C,A)'*S2) ).'; 
        temp3 = B'*B;
        temp4 = kr(A,B);
        C = ((temp2.*temp3)\ (temp4'*S3) ).';

        %% update Z
        S = permute(reshape( temp4*C.',I2,I1,I3),[2,1,3]); % mode-3 matrix -> tensor
        S(omega) = Y(omega); 

        phi=norm(S3-temp4*C.','fro');
        if (abs((phi-phi_old)/phi_old)<Tol)
         break;
        end
    end

    if phi<phiprevious
      out1 = A; out2 = B; out3 = C;
      phiprevious = phi;
    end
end
end