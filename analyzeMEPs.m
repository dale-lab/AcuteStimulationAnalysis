function [p2p1to9] = analyzeMEPs(MEPs)

    p2p1to9 = zeros(31, 1);
    high = 0;
    low = 0;

    for j = 1:31
        for i = 1:225
            if MEPs(j, i) >= high
                high = MEPs(j, i);
            elseif MEPs(j, i) <= low
                low = MEPs(j, i);
            end
        end
        p2p1to9(j, 1) = high-low;
        high = 0;
        low = 0;
    end

end