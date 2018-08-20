function [y,x] = pamHardThreshold2(input)

input = input(:);
symbols = [-3 -1 3 1];
input = repmat(input,1,length(symbols));
symbols = repmat(symbols,size(input,1),1);

[~,x] = min(abs(input - symbols).^2,[],2);
y = symbols(:,x);
y = y(1,:);
y = y(:);
x = x - 1;
