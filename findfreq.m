function fr = findfreq(sig) %in Hz
    delta = 1e-5; %remember that only with dt 0.01 and t 1000 ms
    val = sig;
    h=detrend(val);
    [p,peaks]=findpeaks(h,'MinPeakProminence',(max(h)-min(h))/4);
    fr=1/(mean(diff(peaks))*delta);
return