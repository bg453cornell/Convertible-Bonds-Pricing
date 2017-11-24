function [ price ] = LSM( cp, X, T, r, coupon, sigma, mcallschedule, Nrepl, Nstep)
%UNTITLED2 Summary of this function goes here
%cp:当前股价
%X：初始转股价
%T：债券期限
%r:无风险利率
%coupon：票面利率（向量形式）
%sigma：年化日波动率
%mcallschedule:到期赎回价
%Nrepl：路径总数量
%Nstep：每条路径的步数
%   Detailed explanation goes here
dt = T/Nstep;
s = cp * ones(Nrepl, Nstep);

%MonteCarloSimulation生成Nrepl条股价序列
for i = 1:Nrepl
    for j = 1:(Nstep-1)
        s(i, j+1) = s(i, j) * exp((r-0.5*sigma^2) * dt+sigma * sqrt(dt)* randn);
    end
end

X = X * ones(Nrepl, Nstep);  %生成每条路径的转股价
p = zeros(Nrepl, 1); %存放转债价格

%为防止回售而调整转股价
for j = 1:Nrepl
    for k = 0:T-1
        for i = (1 + k*Nstep/T + round(0.9/(k+1))*0.5*k*Nstep/T):((k+1)*Nstep/T -30*floor((k+1)/T))
            if size(find(s(j,i:i+29)<X(j,i:i+29)*0.8))>=15
                X(j,i+30:Nstep) = mean(s(j,i:i+29));
                break
            end
        end
    end
end
   
for j=1:Nrepl %统计触发赎回条款的路径
    for i=(0.5*Nstep/T+1):Nstep-30
        if (find(s(j, i:i+29)>= 1.3*X(j,i:i+29)))>=15 
            p(j,1)=((100/X(j,i+30))*s(j,i+30)+sum(coupon(1,1:floor(i*T/Nstep))))*exp(-r*dt*i);
            break
        end
    end
end

%将股票模拟路径触发赎回去除，剩下的路径转股价已经确定，完全不受可转债附加条款的影响
for m = 1:Nrepl
    if p(m,1)>0 %如果该条路径已经提前转股
        s(m,:)=0; %删去该条路径中的所有股价数据
    end
end

A = 100*s(:,Nstep)./X(:,Nstep); %每条路径最终转换价值的列向量A
cashflows = max(A, mcallschedule); %mcallschedule为到期赎回价值，cashflows用来存储每条路径在不同时刻的cb价值
for i=1:Nrepl
    if A(i,1)==0 %由于赎回导致提前转股的路径的股价已经删去，因此该路径的A就为0
        cashflows(i, 1)=0;
    end
end

for step = (Nstep-1):-1:0.5*Nstep/T
    cashflows(:,1) = cashflows(:,1).* exp(-r*dt); %将上一期的现金流贴现到现在
    ConversionPrice = 100*s(:,step)/X(:,step); %计算第step期转换价值
    PathNo = find(s(:, step)>X(:,step)); %记录股价大于转股价的位置
    x = s(PathNo, step); 
    y = cashflows(PathNo, 1);
    RegressionMatrix = [ones(length(y),1),x, x.^2];
    b = regress(y,RegressionMatrix);
    HoldingPrice = [ones(length(s(:,step)),1),s(:,step),s(:,step).^2]*b; %生成step期的持有价格
    cashflows(PathNo,1)=max(ConversionPrice(PathNo,1),HoldingPrice(PathNo,1));
end

cashflows0 = cashflows.*exp(-r*dt*(0.5*Nstep/T));
price = mean(cashflows0(:,1)+p(:,1));

end

