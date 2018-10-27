function [RFT]=DALBY_RAIN(PARAMS,t)
% rainfall
%PARAMS.r_f          = [0.00171,0.0171,0.000171];  % rainfall constant [normal,flood,drought]
%PARAMS.r_t          = 1;        % rain type 1=normal, 2=flood, 3=drought
%PARAMS.r_m          = 2;        % rain model 1=constant, 2=cosine, 3=interpol

% t needs to be a value between 4-15..... can be a decimal...
% 4.0 = beginning of jan, 15.9999 = end of december.

% fix t for using in this section
 t=16.234 % test t value
int_t = floor(t)
decimal_t = t-int_t

% drought factor: set this to make it flood (more than 1) or drought (less than 1)
% change this to accept a variable from PARAMS...
drought_factor = 0.654; % 1=100% of rainfall, 1.2=more rainfall, 0.8=less rainfall

x = 1:18;
% a years worth of data with repeated 3 months prior and after
% proper year needed is from point 4-15
y = [58.4 71 94.1 76.1 80 57.4 19.1 36.7 31.2 22.7 23.5 30.4 58.4 71 94.1 76.1 80 57.4];

csp = csape(x,y,'periodic'); %enforces periodicity
xx = linspace(1,18,101);
plot(x,y,'o',xx,ppval(csp,xx),'-');

% Now print out the coefficients
coefficients = csp.coefs
fprintf('The equation for the different segments are:\n');
for k = 1 : size(coefficients, 1)
	fprintf('y = %7.4f * x^3 + %7.4f * x^2 + %7.4f * x + %7.4f\n', coefficients(k,:));
end

% time input t from 1-12 will return the avg rainfall for that month
month = int_t;
for i = 0:3
    month_rain(i+1) = csp.coefs((i+month)+(17*i-i));
end

fprintf('month: %f', month)
fprintf(' Decimal: %f', decimal_t)
month_rain()

RFT = 0; % initialise Rainfall WRT time
RFT = month_rain(1)*decimal_t^3 + month_rain(2)*decimal_t^2 + month_rain(3)*decimal_t^1 + month_rain(4)*drought_factor;

%put in for loop
% for i=0:2
%     avg = avg + month_rain(i+1)*decimal_t^();
% end

fprintf('avg rain for given t: %f', RFT)


end




