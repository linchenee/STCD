function [est_para] = MUSICstep2(Y, T, L, delta_f, c0)
% -------------------------------------------------------------------  
% Multiple signal classification (MUSIC) method for range estimation.
% version 1.0 - 01/06/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -------------------------------------------------------------------
N1 = size(Y,1);
num_grid_Para = 1e4;
grid_Para = linspace(0,300,num_grid_Para);
min_Para = 0;
max_Para = 300;
dic = exp(-1j*2*pi*delta_f*(0:N1-1).'*(2*grid_Para/c0));

R = 1/T * (Y*Y');
[U,D] = eig(R);
[~, ind] = sort(diag(D), 'descend');
U = U(:, ind);
Uz = U(:,L+1:end);
U = Uz*Uz';

%% First stage: Searching on girds
spectrum = zeros(size(dic,2), 1);
for i = 1:size(dic,2)
    spectrum(i,1) = (1/abs( dic(:,i)' * U * dic(:,i) ));
end

interval1 = 1;
est_para = findmax(grid_Para,spectrum,L,interval1);

%% Second stage: Searching off girds
interval = 1;
obj_func = @(x) abs( beam(x,N1,delta_f,c0)' * U * beam(x,N1,delta_f,c0) );
for i = 1:L
    est_para(i) = fminsearchbnd(obj_func,est_para(i), max(est_para(i)-interval,min_Para), min(est_para(i)+interval,max_Para));
end

function [y] = beam(x,N1,delta_f,c0)
    y = exp(-1j*2*pi*delta_f*(0:N1-1).'*(2*x/c0));
end

end
