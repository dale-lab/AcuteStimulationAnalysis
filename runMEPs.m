function [PreInjIpsi, PostInjIpsi, PostMEPsIpsi, PreInjContra, PostInjContra, PostMEPsContra, time] = ...
        runMEPs(PreInjIpsi, PostInjIpsi, PostMEPsIpsi, PreInjContra, PostInjContra, PostMEPsContra, fname, frequency)

    load(fname); 
    number = extractBefore(fname, "p");
    type = extractBefore(fname, ".mat");
    type = extractAfter(type, "p");
    
    if number == "17" || number == "14" || number == "21" || number == "16" 
        injury = "L";
    elseif number == "28" || number == "24" || number == "23" || number == "27" || number == "18" || number == "15"
        injury = "R";
    end
    
    %0 = preinj 1 = postinj 2 = poststim
    preinj = strcmp(type, "reinj");
    postinj = strcmp(type, "ostinj");
    poststim = strcmp(type, "oststim");
    
    if preinj == 1
        timepoint = 0;
    elseif postinj == 1
        timepoint = 1;
    elseif poststim == 1
        timepoint = 2;
    end
    
    if injury == 'R'
        ipsiraw = strcat(who('*Ch1'), '.values');
        ipsiRAW = eval(ipsiraw{1});
        contraraw = strcat(who('*Ch3'), '.values');
        contraRAW = eval(contraraw{1});
        ipsismooth = strcat(who('*Ch2'), '.values');
        ipsiSMOOTH = eval(ipsismooth{1});
        contrasmooth = strcat(who('*Ch4'), '.values');
        contraSMOOTH = eval(contrasmooth{1});
    else
        ipsiraw = strcat(who('*Ch3'), '.values');
        ipsiRAW = eval(ipsiraw{1});
        contraraw = strcat(who('*Ch1'), '.values');
        contraRAW = eval(contraraw{1});
        ipsismooth = strcat(who('*Ch4'), '.values');
        ipsiSMOOTH = eval(ipsismooth{1});
        contrasmooth = strcat(who('*Ch2'), '.values');
        contraSMOOTH = eval(contrasmooth{1});
    end
    stim = strcat(who('*Ch8'), '.values');
    STIM = eval(stim{1});
    time = transpose((1/frequency:1/frequency:numel(ipsiRAW)/frequency));
    
    if numel(time) > contraSMOOTH
        time = time(1:end-1);
    end
    
    stimtimes = zeros(1, 1);
    i = 1;
    count = 1;
    contracount = 0;
    mean = 0;
    
    while i < numel(STIM) - 2
        if STIM(i) > 0.05
            while STIM(i+1) > 0.05 && STIM(i + 2) > 0.05
                i = i + 1;
            end
            stimtimes(count) = i;
            count = count + 1;
            mean = mean + contraSMOOTH(i);
            contracount = contracount + 1;
            i = i + 12500;
        end
        i = i + 1;
    end
    threshold = mean/contracount;
    
    nonstims = zeros(1, 1);
    count = 1;
    yes = 1;
    b = 1201;
    
    
    while b < numel(contraSMOOTH)-1250
        if contraSMOOTH(b-1) < threshold && contraSMOOTH(b) > threshold
            for g = b:b+1250
                if contraSMOOTH(g) < threshold
                    yes = 0;
                end
            end
            if yes == 1
                for t = b-1200:b+1200
                    if STIM(t) > 0.05
                        yes = 0;
                    end
                end
                if yes == 1
                    nonstims(count) = b;
                    count = count + 1;
                end
                while b <= numel(contraSMOOTH)&& contraSMOOTH(b) > threshold
                    b = b + 1;
                end
            end
            b = b + 1;
        end
        yes = 1;
        b = b + 1;
    end
    
    pulsenum = 1;
    count = 0;
    v = 1;
    avenostimipsi = zeros(1, 525);
    avestimipsi = zeros(1, 525);
    avenostimcontra = zeros(1, 525);
    avestimcontra = zeros(1, 525);
    thestims = floor(numel(stimtimes)/4);
    max = thestims*4;
    while v < max
        while count < 4
            for i = stimtimes(v)+1:stimtimes(v)+525
                avestimipsi(1, i-stimtimes(v)) = avestimipsi(1, i-stimtimes(v)) + ipsiRAW(i);
                avestimcontra(1, i-stimtimes(v)) = avestimcontra(1, i-stimtimes(v)) + contraRAW(i);
            end
            count = count + 1;
            if v < max
                v = v + 1;
            end
        end

        count = 0;
            avestimipsi = avestimipsi/4;
            avestimcontra = avestimcontra/4;
            if timepoint == 0
                PreInjIpsi(pulsenum, :) = avestimipsi;
                PreInjContra(pulsenum, :) = avestimcontra;
            elseif timepoint == 1
                PostInjIpsi(pulsenum, :) = avestimipsi;
                PostInjContra(pulsenum, :) = avestimcontra; 
            elseif timepoint == 2
                PostMEPsIpsi(pulsenum, :) = avestimipsi;
                PostMEPsContra(pulsenum, :) = avestimcontra; 
            end

        pulsenum = pulsenum + 1;
        avestimipsi = zeros(1, 525);
        avestimcontra = zeros(1, 525);
    end

end
