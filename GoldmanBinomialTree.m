function Price=GoldmanBinomialTree(cp,vol,rate,T,principal,coupon,conversionprice,creditspread,callschedule,mcallschedule,N)
%cp:当前股价
%vol:股票波动率
%rate:无风险利率
%T:转债期限
%principal:本金
%frequency:付息频率
%conversionprice：转股价?
%callschedule：赎回价
%mcallshcedule：到期赎回价
conversionratio=100/conversionprice;
dt=T/N;
up=exp(vol*sqrt(dt));
down=1/up;
s=zeros(N,N);
v=zeros(N,N);
p=zeros(N,N);
%建立二叉树
for n=0:N-1;
    for j=0:n;
        s(n+1,j+1)=cp*(up)^j*(down)^(n-j);
    end;
end;
interest=zeros(1,T);
interest=principal*coupon;
for j=N:-1:1;
    payment=principal+interest;
    if s(N,j)>conversionprice;
        v(N,j)=max(conversionratio*s(N,j),mcallschedule);
    else
        v(N,j)=mcallschedule;
    end;
    if v(N,j)==conversionratio*s(N,j);
        p(N,j)=1;
    else
        p(N,j)=0;
    end;
end;

for i=N-1:-1:1;
    for j=i:-1:1;
        call=callschedule;
        p(i,j)=0.5*(p(i+1,j+1)+p(i+1,j));
        creditadjustedrate=p(i,j)*rate+(1-p(i,j))*creditspread;
        H=0.5*((v(i+1,j+1)+dt*interest(ceil(T*i/N)))/(1+creditadjustedrate*dt)+(v(i+1,j)+dt*interest(ceil(T*i/N)))/(1+creditadjustedrate*dt)); 
        if s(i,j)>=conversionprice;
            v(i,j)=max(conversionratio*s(i,j),min(H,call+T*(N-i)*dt*interest(ceil(T*i/N))));
        else
            v(i,j)=min(H,call+T*(N-i)*dt*interest(ceil(T*i/N)));
        end;
    end;
end;
Price=v(1,1);