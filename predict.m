%% This is a new function called predict. 
% It allows you to create a new vector of x_hat (predicted future values)
function [x_hat] = predict(a,x)
x_hat = dot(-a(2:end),fliplr(x((length(x)-length(a)+2):end)));
%new_vec = [x x_hat];
end
%%