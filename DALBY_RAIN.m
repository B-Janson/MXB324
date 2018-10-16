function []=DALBY_RAIN()

x = 1:18;
% a years worth of data with repeated 3 months prior and after
% proper year needed from point 4-15
y = [58.4 71 94.1 76.1 80 57.4 19.1 36.7 31.2 22.7 23.5 30.4 58.4 71 94.1 76.1 80 57.4];
cs = spline(x,[0 y 0]);
xx = linspace(1,18,101);
plot(x,y,'o',xx,ppval(cs,xx),'-');
% Now print out the coefficients
coefficients = cs.coefs
fprintf('The equation for the different segments are:\n');
for k = 1 : size(coefficients, 1)
	fprintf('y = %7.4f * x^3 + %7.4f * x^2 + %7.4f * x + %7.4f\n', coefficients(k,:));
end


end