function []=DALBY_RAIN()

X = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
Y = [94.1 76.1 80 57.4 19.1 36.7 31.2 22.7 23.5 30.4 58.4 71 94.1 76.1];
[P,R,S] = lagrangepoly(X,Y);

xx = 0.5 : 0.01 : 14.5;
plot(xx,polyval(P,xx),X,Y,'or',R,S,'.b',xx,spline(X,Y,xx),'--g')
grid
axis([0.5 14.5 -50 110])


end