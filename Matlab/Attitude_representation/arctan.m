function [r] = arctan(num,den)
%arctan Returns the atan with coordinate plane correction
%   Corrects fepending to the corresponding cuadrant
    
    if (num >= 0 && den >= 0)
        r = atan(num / den);
    elseif (num >= 0 && 0 > den)
        r = atan(num / den) + pi;
    elseif (num < 0 && den < 0)
        r = atan(num / den) - pi;
    else
        r = atan(num / den);
    end
end

