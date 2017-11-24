function [ price ] = MonteCarlo(F0, T, cr, alpha, miu, sigmar, cH, s, sigmaH, td, EndDft, ExoDft, ti,Nrepl, Nstep)
%UNTITLED2 Summary of this function goes here
%F0：贷款本金
%T：贷款期限（年）
%cr：当前利率
%alpha：CIR模型中利率调整速度
%miu：CIR模型中利率均值
%sigmar：CIR模型中利率波动率
%cH：当前房价
%s：租金率
%sigmaH：GB模型中房价波动率
%td：违约成本率（相对房价的百分比）
%EndDft：dt时间内理性违约比例
%ExoDft：dt时间内非理性违约比例
%ti: 房产变现成本率
%Nrepl：路径总数量
%Nstep：每条路径的步数

%   Detailed explanation goes here
dt = T/Nstep;
r = cr * ones(Nrepl, Nstep);
H = cH * ones(Nrepl, Nstep);

A = F0 * (c/12) * (1+c/12)^(12*Nstep)/((1+c/12)^(12*Nstep)-1); %每月支付金额
F = F0 * ones(1, Nstep);
%生成贷款余额序列
for i = 1:(Nstep-1)
    F(1, i+1) = F(1, i) * (1+c/12)-A;
end

%MonteCarloSimulation生成Nrepl条利率和房价序列
for i = 1:Nrepl
    for j = 1:(Nstep-1)
        r(i, j+1) = r(i, j) + alpha*(miu-r(i, j))*dt + sqrt(r(i, j))*sigmar*sqrt(dt)*randn;
        H(i, j+1) = H(i, j)* exp((r(i, j) - s - 0.5*sigmaH^2) * dt+sigmaH * sqrt(dt)* randn);
    end
end

%生成理性违约概率序列
p = zeros(Nrepl, Nstep+1);
k = zeros(Nrepl, Nstep+1);
W = zeros(Nrepl, Nstep+1);
for i = 1:Nrepl
    for j = 1:Nstep
        k(i, j+1) = (log((c/s*F(1,j+1)-r(i,j+1)/s*td)/H(i,j)) - ((r(i,j+1)-s)-0.5*sigmaH^2)*dt)/(sigmaH*sqrt(dt));%违约距离
        pl = normcdf(k(i,j+1),0,1); %出现期权违约最优情况的概率
        p(i, j+1) = pl *  EndDft * dt + (1-pl) * (1 - EndDft) * dt + ExoDft * dt; %违约概率
        W(i, j+1) = H(i,j+1)*(1-ti); %违约现金流
    end
end

%构建贷款价值序列
V = zeros(Nrepl, Nstep + 1);
V(:,Nstep+1) =(ones(Nrepl,1) - p(:, Nstep+1)).* A + p(:,Nstep+1).* W(i, j+1);
for i = 1:Nrepl
    for j = Nstep:-1:1
        V(i,j) = (1 - p(i, j))* V(i, j+1)*exp(-r(i,j)*dt) + p(i,j)* W(i,j);
    end
end

price = mean(V(:,1));

end

