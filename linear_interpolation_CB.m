StockPremium = xlsread('CBs.xlsx',1,'K5:K22')
CBPrice = xlsread('CBs.xlsx',1, 'D5:D22')
ConversionValue = xlsread('CBs.xlsx',1, 'I5:I22')
ConversionPremium = (CBPrice - ConversionValue) ./ ConversionValue
data = [StockPremium ConversionPremium]
data1 = sortrows(data, 1)
interPoint = -0.7507
y1 = interp1(data1(:,1), data1(:,2), interPoint,'linear')
y2 = interp1(data1(:,1), data1(:,2), interPoint,'spline')
y = [y1 y2]
plot(interPoint, y, '*', data1(:,1), data1(:,2), 'go-')
x = xlsread('CBs.xlsx',1, 'I4')
price = y .* x



