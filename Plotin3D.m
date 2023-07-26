clear
clc

load('14MEPs.mat');
index = 2;
side = 'C';

frequency = 25000;
time = 1/frequency:1/frequency:0.01;
time = time*1000;
z = 0.1:0.1:3.1;
z = transpose(z);

if index == 0
    if side == 'I'
        y = PreInjIpsi(:, 1:250);
        % c = NonPreInjIpsi(:, 1:250);
    elseif side == 'C'
        y = PreInjContra(:, 1:250);
        % c = NonPreInjContra(:, 1:250);
    end
elseif index == 1
    if side == 'I'
        y = PostInjIpsi(:, 1:250);
        % c = NonPostInjIpsi(:, 1:250);
    elseif side == 'C'
        y = PostInjContra(:, 1:250);
        % c = NonPostInjContra(:, 1:250);
    end
elseif index == 2
    if side == 'I'
        y = PostMEPsIpsi(:, 1:250);
        % c = NonPostMEPsIpsi(:, 1:250);
    elseif side == 'C'
        y = PostMEPsContra(:, 1:250);
        % c = NonPostMEPsContra(:, 1:250);
    end
end

[X, Z] = meshgrid(time, z);

% tiledlayout(2, 1)
% nexttile
surf(X, Z, y);
colorbar
caxis([-0.5 0.5])
title('Stimulated Breaths')
xlabel('Time (ms)')
ylabel('Current (mA)')
%zlim([-0.5 0.5])

% nexttile
% surf(X, Z, c);
% colorbar
% caxis([-0.5 0.5])
% title('Unstimulated Breaths')
% xlabel('Time (ms)')
% ylabel('Current (mA)')
%zlim([-0.5 0.5])