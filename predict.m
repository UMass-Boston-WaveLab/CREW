%% This is a new function called predict. 
% It allows you to create a new vector of x_hat (predicted future values)
function [x_hat] = predict(a,x)
x_hat = sum(-fliplr(a(2:end)).*x(1:(end-1)));
%new_vec = [x x_hat];
end
%%