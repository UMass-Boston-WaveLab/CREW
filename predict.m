
function [new_vec] = predict(a,x)
x_hat = dot(a,x);
new_vec = [x(2:end) x_hat];
end