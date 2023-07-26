%% Load file, set variables

clc
clear

load('AI28postMEPs.mat');


%preinj = 1, postinj = 2, stim/sham = 3, postMEPs = 4, phrenicotomy = 5
timepoint = 4;
injury = 'R'; %L or R for injury side
stiman = 'N';%Y or N for stim vs sham
maxwidth = 500;%maximum width of EKG contaminations to blank on either side of event identification
catchlow = -0.075;%lower limit of ipsilesional noise
catchhigh = 0.075;%upper limit of ipsilesional noise
thresh = 0.045;%threshold to identify breaths of contralesional moving average
high = 4;%slope value indicating EKG contamination
artificialadd = 550; %# of extra data points to blank even if slope is below low (use if negative dip in slope after EKG based on pattern)
contrazero = 0.0225;%moving average noise level above 0

%likely don't change these variables
blank = 1;%1 or 0 for blank EKG or not
low = 0;%slope value to cut off EKG blanking time
frequency = 25000;%sampling frequency of imported data
wiggle = 2000; %how far off +/- expected distance of EKGs should be searched for missing event

%% Format data
%Right raw = channel 1, right moving averaged = channel 2 
%left raw = channel 3, left moving averaged = channel 4
%stim marker = channel 8

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
time = transpose((0:1/frequency:numel(ipsiRAW)/frequency));


%% Initialize variables
pulses = 0;
maxpeak = 0;
AUC = 0; %these three characteristic holders reset each breath below
breath = 1;%results index
i = 1; %moves through contralesional data
w = 0;%how far moved in a single breath; makes sure it's long enough

results = zeros(1, 10);
starts = zeros(1, 1);
ends = zeros(1, 1);

%% Calculate contralateral metrics
%output; results(1) = breath start time; results(2) = breath end time;
%results(3) = breath length; results(4) = maximum moving averaged value
%results(5) = modulus; results(6) = duty cycle
%results(7) = breathing rate; results(8) = number of stimulations in the breath

while i < numel(time) - 1
    if contraSMOOTH(i) > thresh
        starttime = i;%start of breath
        while i <= numel(contraSMOOTH) && contraSMOOTH(i) > thresh
            w = w + 1;
            i = i + 1;
        end
        endtime = i;%end of breath
        if w >= 4000%check breath is long enough, w/frequency seconds minimum
            starts(breath) = ceil(starttime);
            ends(breath) = ceil(endtime);
            for t = starttime:endtime
                if t <= numel(contraRAW)
                    AUC = contraSMOOTH(t)/frequency - contrazero/frequency + AUC;  %calculate modulus, subtract zero value
                    if contraSMOOTH(t)-contrazero > maxpeak
                        maxpeak = contraSMOOTH(t) - contrazero; %find highest point, subtract zero value
                    end
                    if t < numel(STIM) && STIM(t) < 0.05 && STIM(t+1) >= 0.05
                        pulses = pulses + 1;
                    end
                end
            end
            %save results, move onto next breath
            results(breath, 1) = starttime/frequency;
            results(breath, 2) = endtime/frequency;
            results(breath, 3) = (endtime-starttime)/frequency;
            results(breath, 4) = maxpeak;
            results(breath, 5) = AUC;
            results(breath, 8) = pulses;
            breath = 1 + breath;
        end
        pulses = 0;
        maxpeak = 0;
        AUC = 0;
        w = 0;
    end
    i = i + 1;
end

for u = 1:numel(results(:, 1))-1
    total = results(u+1, 1) - results(u, 1); %full length of breathing cycle
    results(u, 6) = results(u, 3)/total; %duty cycle
    results(u, 7) = 60/total; %rate
end


%% Get mean centered and cumulative values on ipsilesional side

meancentered = zeros(1, numel(ipsiRAW));
cumulativesum = zeros(1, numel(ipsiRAW));

%delete noise contribution
for f = 1:numel(ipsiRAW)
    if ipsiRAW(f) > catchhigh
        meancentered(f) = ipsiRAW(f) - catchhigh;
    elseif ipsiRAW(f) < catchlow
        meancentered(f) = ipsiRAW(f) - catchlow;
    else
        meancentered(f) = 0;
    end
end

%calculate cumulative sum
cumulativesum(1) = abs(meancentered(1));
for g = 2:numel(meancentered)
    cumulativesum(g) = cumulativesum(g-1) + abs(meancentered(g));
end

%calculate slope of cumulative sum
sumslope = zeros(1, numel(cumulativesum));
for i = 251:numel(sumslope)
sumslope(i) = (cumulativesum(i)-cumulativesum(i-250));
end


%% First pass blank major EKG contributions

fixedcumulative = zeros(1, numel(cumulativesum)); %holds blanked qsum
allowed = ones(1, numel(cumulativesum)); %holds if EKG blanked or not; 1 = no blank, 0 = blank
o = 1; %moves through sumslope
k = 0; %holds number moved forwards from detection reset each breath
m = 0; %holds number move backwards from detection reset each breath

while o <= numel(sumslope)
    if sumslope(o) > high
        t = o;
        while t <= numel(sumslope) && sumslope(t) > low && m < maxwidth %blank forwards
            allowed(t) = 0;
            t = t + 1;
            m = m + 1;
        end
        for r = t:t + artificialadd %if need to blank further forwards, does the blanking
            if r <= numel(cumulativesum)
                allowed(r) = 0;
            end
        end
        b = o;
        while b > 0 && sumslope(b) > low && k < maxwidth %blank backwards
            allowed(b) = 0;
            b = b - 1;
            k = k + 1;
        end
        o = t + 1;  
        k = 0;
        m = 0;
    end
    o = o + 1;
end

%% Second pass blank additional detections based on EKG rate

count = 1;%index
ekgs = zeros(1,1);%stores time of all EKG events
widths = zeros(1, 1); %stores width of detected EKG events
p = 1;%moves through allowed
l = 0;%stores width of detected EKGs


if blank == 1
    while p <= numel(allowed)
        if allowed(p) < 1 %find time and width of all detections
            ekgs(count) = p; %stores width of detected EKGs
            p = p + 1;
            while p <= numel(allowed) && allowed(p) < 1
                p = p + 1;
                l = l + 1;
            end
            widths(count, 1) = l;
            count = count + 1;
            l = 0;
        end
        p = p + 1;
    end

    diff = zeros(numel(ekgs)-1, 1); %holds distance between each EKG event

    for h = 2:numel(ekgs)
        diff(h-1) = ekgs(h) - ekgs(h-1); %calculate distances
    end

    sorted = sort(diff); %sort distances between by length

    count =0; %holds number of detections for each bin
    bins = zeros(1, 1); %find number of detected EKG in each width - looking for average heartrate
    binnum = 1;%index
    k = 1; %moves through each detection

    while k < numel(sorted)
        while k < numel(sorted) && sorted(k) < binnum*100
            count = count + 1;
            k = k + 1;
        end
        bins(binnum, 1) = count; %for each bin of 100/frequency time width, how many EKG detections?
        count = 0;
        binnum = binnum + 1;
    end

    if numel(bins) > 70 %get rid of all counts with way too long a space between them, find the time between EKG that is most common
        [M, I] = max(bins(1:70));
    else
        [M, I] = max(bins);
    end

    index = I * 100;%most common time 
    highindex = (I+1.5)*100;  %allow for variety of times between EKGs for heart rate changes
    lowindex = (I-1.5)*100;

    sum = 0; 
    num = 0;
    for o = 1:numel(diff)
        if diff(o) < highindex && diff(o) > lowindex %if the distance between EKGs is correct
            sum = sum + widths(o);
            num = num + 1;
        end
    end
    blanks = ceil(sum/num);
    %find average blanking time for correctly identified EKGs

    %if EMG bursts are high enough, EKGs are not as obvious, possible that they are mis-identified, treat differently
    for h = 2:numel(ekgs)-1
        if ekgs(h) - ekgs(h-1) > highindex || ekgs(h) - ekgs(h-1) < lowindex || ekgs(h+1) - ekgs(h) > highindex || ekgs(h+1) - ekgs(h) < lowindex
            %if the distance between EKGs is greater or less than
            %expected, get rid of that whole blanked EKG section as it
            %is likely a breath
            for l = ekgs(h-1):ekgs(h+1)
                allowed(l) = 1;
            end
        end
    end
    

    jump = ceil(blanks/2) + index; %distance that you need to move from one detection to find the next 

    b = ekgs(1)+jump;%starting point index for moving through allowed
    t = 1;%starting point index for moving through EKGs
    h = -1;

    %go back through and identify any missed EKGs
    while b < numel(allowed)
        %reinitialize
        yes = 0;%0 or 1 for if the next event was detected or not
        max = 0;%highest slope value found in range
        newind = 0;%where highest slope value is found in range
        for i = b-wiggle:b+wiggle%search area around where next EKG should be +/- wiggle value
            if i <= numel(sumslope) && sumslope(ceil(i)) > max%find the maximum slope and location in that detection area
                max = sumslope(ceil(i));
                newind = i;
            end
            if i <= numel(sumslope) && allowed(ceil(i)) == 0%check to see if the area at highest slope is already blanked
                yes = 1;
            end
        end
        if yes == 1 && t < numel(ekgs)%if area is already blanked, move on to the next one
            t = t + 1;
            while ekgs(t) + jump < b-wiggle
                t = t + 1;
            end
            b = ekgs(t) + jump;
        else%if area is not blanked, do so
            if newind == 0 && t <numel(ekgs)% check for case that max slope is at very beginning (no blank)
                t = t + 1;
                b = ekgs(t) + jump;
            elseif newind == 0 && t >=numel(ekgs)% check for case that max slope is at the beginning (no blank) plus working past last detected EKG
                b = numel(allowed);
            else
                for k = newind-blanks/2:newind+blanks/2%blank the EKG
                    allowed(ceil(k)) = 0;
                end
                %move index over 
                if b ~= h
                    h = b;
                    b = newind + jump;
                else
                    b = numel(allowed);
                end
            end
        end
    end
end


%% get output characteristic

%get new cumulative values with EKG blanking
fixedcumulative(1) = 0;
for y = 2:numel(sumslope)
   if allowed(y) == 1
       fixedcumulative(y) = fixedcumulative(y-1) + abs(meancentered(y));
   else
       fixedcumulative(y) = fixedcumulative(y-1) + 0;
   end
end


%find the increase in cumulative sum throughout each breath and between
%breaths
for t = 1:breath - 1
    Ti = starts(t);
    Te = ends(t);
    if Te > numel(fixedcumulative)
        Te = numel(fixedcumulative);
    end
    if t < breath - 1
        Ti2 = starts(t+1);
    else
        Ti2 = numel(fixedcumulative);
    end
    during = fixedcumulative(Te) - fixedcumulative(Ti);
    between = fixedcumulative(Ti2) - fixedcumulative(Te);
    results(t, 9) = during;
    results(t, 10) = between;
end

%% Minute by minute and averages overall 

minbymin = zeros(ceil(results(end, 1)/60), 9);
mins = 1;
b = 1;

averages = zeros(3, 8);

%find sum minute by minute values

if stiman == 'N' && timepoint == 3
    while b <= numel(results(:, 1)) && mins <= numel(minbymin(:, 1))
        while results(b) < mins * 60
            minbymin(mins, 1) = minbymin(mins, 1) + 1;
            minbymin(mins, 2:9) = minbymin(mins, 2:9) + results(b, 3:10);
            if b <numel(results(:, 1))
                b = b + 1;
            else
                break
            end
        end
        mins = mins + 1;
    end

    %get average/breath
    for c = 1:numel(minbymin(:, 1))
        minbymin(c, :) = minbymin(c, :)/minbymin(c, 1);
        minbymin(c, 1) = c;
    end

    averages(1, :) = mean(minbymin(1:2, 2:9));
    averages(2, :) = mean(minbymin(3:end-2, 2:9));
    averages(3, :) = mean(minbymin(end-1:end, 2:9));

elseif stiman == 'Y' && timepoint == 3
    while results(b, 8) == 0
        b = b + 1;
    end
    stimstart = results(b, 1);
    c = numel(results(:, 1));
    while results(c, 8) == 0
        c = c - 1;
    end
    stimend = results(c, 1);

    b = 1;

    while results(b, 1) < stimstart
        minbymin(mins, 1) = minbymin(mins, 1) + 1;
        minbymin(mins, 2:9) = minbymin(mins, 2:9) + results(b, 3:10);
        b = b + 1;
    end
    mins = mins + 1;

    while results(b, 1) <= stimend
        while results(b) < mins * 60 + stimstart
            minbymin(mins, 1) = minbymin(mins, 1) + 1;
            minbymin(mins, 2:9) = minbymin(mins, 2:9) + results(b, 3:10);
            b = b + 1;
        end
        mins = mins + 1;
    end

    while b <= numel(results(:, 1))
        minbymin(mins, 1) = minbymin(mins, 1) + 1;
        minbymin(mins, 2:9) = minbymin(mins, 2:9) + results(b, 3:10);
        b = b + 1;
    end

    %get average/breath
     for c = 1:numel(minbymin(:, 1))
        minbymin(c, :) = minbymin(c, :)/minbymin(c, 1);
        minbymin(c, 1) = c;
    end

    if stiman == 'Y'
        stimbegin = 0;
        for r = 1:numel(results(:, 1))
            if results(r, 8) == 0
                stimbegin = stimbegin+1;
            else
                break
            end
        end
    
        r = numel(results(:, 1));
        while r >= 0
            if results(r, 8) == 0
                r = r - 1;
            else
                break
            end
        end
    else
        stimbegin = 0;
        for r = 1:numel(results(:, 1))
            if results(r, 1) <= 120
                stimbegin = stimbegin+1;
            else
                break
            end
        end
    
        r = numel(results(:, 1));
        while r >= 0
            if results(r, 1) >= results(end, 1) - 120
                r = r - 1;
            else
                break
            end
        end
    end

    averages(1, :) = mean(results(1:stimbegin, 3:10));
    averages(2, :) = mean(results(stimbegin+1:r, 3:10));
    averages(3, :) = mean(results(r+1:end, 3:10));

end

if timepoint ~= 3
    averages = zeros(1, 8);
    averages = mean(results(:, 3:10));
end


%characteristics

characteristics = zeros(14, 1);

characteristics(1) = timepoint;
characteristics(2) = injury;
characteristics(3) = stiman;
characteristics(4) = maxwidth;
characteristics(5) = catchlow;
characteristics(6) = catchhigh;
characteristics(7) = thresh;
characteristics(8) = high;
characteristics(9) = artificialadd;
characteristics(10) = contrazero;
characteristics(11) = blank;
characteristics(12) = low;
characteristics(13) = frequency;
characteristics(14) = wiggle;


%% double check output 

hold off
plot(time(1:end-1), allowed)
hold on
plot(time(1:end-1), meancentered)
ylim([-0.2 1.2])

