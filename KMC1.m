% KMC under a constant length

clear all;

e0 = 50;
n = 100;

yb = 50;
Es = 46.13404982;
Eh = 40.61457757;

%yb = 80;
%Es = 43.05273078;
%Eh = 128.5724768;

m = 100;
dts = exp(Es)/(n+1)/m;
dth = exp(Eh)/m;

t(1) = 0;
S(1) = 1;
k = 1;
%tmax = 1.e20;

kh = 1;
%th(kh) = 0;
dwth(kh) = 0;

ks = 0;

while ks < 201
    k = k+1;
    S(k) = S(k-1);

    rn = rand(1);

    if S(k) > 0
        dt = dts;
        dwth(kh) = dwth(kh)+dt;

        if rn < 1/m
            S(k) = 0;
            ks = ks+1;
            dwts(ks) = 0;
%            ts(ks) = t(k-1)+dt;
        end
    else
        dt = dth;
        dwts(ks) = dwts(ks)+dt;

        if rn < 1/m
            S(k) = 1;
            kh = kh+1;
            dwth(kh) = 0;
%            th(kh) = t(k-1)+dt;
        end
    end

    t(k) = t(k-1)+dt;
%    if t(k) > tmax
%        break;
%    end
end

sum_dwth = 0;
sum_dwts = 0;

ndwt = ks-1;
for i = 1:ndwt
%    dwth(i) = ts(i)-th(i);
%    dwts(i) = th(i+1)-ts(i);

    sum_dwth = sum_dwth + dwth(i);
    sum_dwts = sum_dwts + dwts(i);
end

mean_twh = sum_dwth/ndwt;
mean_tws = sum_dwts/ndwt;
