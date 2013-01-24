function [ts_off]=synchronize_ts(ts1,ts2)    
    tmp1=abs(ts1(1)-ts2);
    [val1 ts_off1]=min(tmp1);
    tmp2=abs(ts2(1)-ts1);
    [val2 ts_off2]=min(tmp2);
    if (val1<val2)
        ts_off=ts_off1;
    else
        ts_off=-ts_off2;
    end
end
