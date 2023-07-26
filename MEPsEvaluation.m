clear
clc

%put all three files in same directory
%MUST be saved with the format 'animalnumber' 'preinj/postinj/poststim'

%name of file to save
currentfile = '28MEPs.mat';
F = dir('*.mat');
namelist = {F.name};

%setup all values stim and no stim
PreInjIpsi= zeros(31, 525);
PostInjIpsi = zeros(31, 525);
PostMEPsIpsi = zeros(31, 525);
PreInjContra= zeros(31, 525);
PostInjContra = zeros(31, 525);
PostMEPsContra = zeros(31, 525);

frequency = 25000;

%move through all files
for i = 1:numel(namelist)

    thename = namelist(i);
    fname = char(thename);

    [PreInjIpsi, PostInjIpsi, PostMEPsIpsi, PreInjContra, PostInjContra, PostMEPsContra,  time] = ...
        runMEPs(PreInjIpsi, PostInjIpsi, PostMEPsIpsi, PreInjContra, PostInjContra, PostMEPsContra, fname, frequency);
    
end

p2p = zeros(31, 6);

[p2p(:, 1)] = analyzeMEPs(PreInjIpsi);
[p2p(:, 2)] = analyzeMEPs(PostInjIpsi);
[p2p(:, 3)] = analyzeMEPs(PostMEPsIpsi);

[p2p(:, 4)] = analyzeMEPs(PreInjContra);
[p2p(:, 5)] = analyzeMEPs(PostInjContra);
[p2p(:, 6)] = analyzeMEPs(PostMEPsContra);

save(currentfile, 'PreInjIpsi', 'PreInjContra', 'PostInjIpsi', 'PostInjContra', 'PostMEPsIpsi', 'PostMEPsContra', 'time', 'p2p')     
        
        
        
        
        
        
        
   
