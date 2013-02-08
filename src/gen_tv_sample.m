function tv = gen_tv_sample(T,L,p)

F = (T+L)/2 - 1;
B = (T-L)/2;

intercept = F-B;

f = 0;
b = 0;
t = 0;
l = 0;

tv = zeros(T-1,1);

for i=1:T-1
    temp = rand;
    d = sum(tv);
    if d == 0
        tv(i)=1;
    elseif d == intercept
        tv(i)=-1;
    elseif sum(tv==-1)==B
        tv(i)=1;
    elseif temp<p
        tv(i)=1;
    elseif temp>p
        tv(i)=-1;
    end
end


