function [est_para] = MUSICstep1(Y, T, L)
% -------------------------------------------------------------------  
% Multiple signal classification (MUSIC) method for angle estimation.
% version 1.0 - 01/06/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -------------------------------------------------------------------
Mr = size(Y,1);
num_grid = 1e3;
grid = linspace(0,90,num_grid)-45;
dic = exp(1j*pi*(0:(Mr-1)).'*sin(grid*pi/180));

[U,D] = eig(1/T * (Y*Y'));
[~, ind] = sort(diag(D), 'descend');
U = U(:, ind);
Uz = U(:,L+1:end);
U = Uz*Uz';

%% First stage: Searching on girds
spectrum = zeros(size(dic,2), 1);
for i = 1:size(dic,2)
    spectrum(i,1) = (1/abs( dic(:,i)' * U * dic(:,i) ));
end

interval = 0.5;
est_para = findmax(grid,spectrum,L,interval);

%% Second stage: Searching off girds
min_Angle = -45;
max_Angle = 45;
interval = 1;
obj_func = @(x) abs( beam(x,Mr)' * U * beam(x,Mr) );
for i = 1:L
    est_para(i) = fminsearchbnd(obj_func,est_para(i), max(est_para(i)-interval,min_Angle), min(est_para(i)+interval,max_Angle));
end

function [y] = beam(x,Mr)
    y = exp(1j*pi*(0:(Mr-1)).'*sin(x*pi/180));
end

end
