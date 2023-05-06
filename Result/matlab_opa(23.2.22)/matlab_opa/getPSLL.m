function [psll,mainpeak, mainpeakpos,maxsidepeak,maxsidepeakpos] = getPSLL(u)
    numpeak = 0;
    u_diff = diff(u);
    for i = 1:length(u)-2
        if u_diff(i)>0 && u_diff(i+1)<0
            numpeak = numpeak + 1;
            p(numpeak) = i;
            mag(numpeak) = u(i);
            if(u(i)<u(i+1))
                p(numpeak) = i+1;
                mag(numpeak) = u(i+1);
            end
        end
    end
    
%     [mainpeak, mainpeakpos] = max(mag);
    
%     mag(mainpeakpos) = 0;
    lenu = length(u);
    midpos = floor(lenu/2);
    mainpeak = max(u(floor(midpos-5:midpos+5)));
%     [maxsidepeak, maxsidepeakpos] = max(mag);
%     if maxsidepeak == 0
        maxsidepeak = max(u(1:midpos-5));
        maxsidepeakr = max(u(midpos+5:end));
        if maxsidepeak < maxsidepeakr
            maxsidepeak = maxsidepeakr;
        end
%     end
    psll = 20*log10(maxsidepeak/mainpeak);
    if isinf(psll)
        disp('inf');
    end
end