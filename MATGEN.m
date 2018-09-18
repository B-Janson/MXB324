function [J]=MATGEN(DIM,h,T)

%% Warmup

n=DIM.n;
m=DIM.m;
XY=DIM.XY;
NT=DIM.NT;

%Initialise Thingos
J=zeros(n*m,n*m);



%% Create Matrix

%Bottom left Corner
[J(1,1),J(1+n,1),J(1,2)]=V1(1,h(1));

for i=2:n*m-1
    if NT(i) == 2
        [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V2(i,h(i));
    elseif NT(i) == 3
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V3(i,h(i));
    elseif NT(i) == 4
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V4(i,h(i));
    elseif NT(i) == 5
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V5(i,h(i));
    elseif NT(i) == 6
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V6(i,h(i));
    elseif NT(i) == 7
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V7(i,h(i));
    elseif NT(i) == 8
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V8(i,h(i));
    elseif NT(i) == 9
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V9(i,h(i));
    elseif NT(i) == 10
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V10(i,h(i));
    elseif NT(i) == 11
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V11(i,h(i));
    elseif NT(i) == 12
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V12(i,h(i));
    elseif NT(i) == 13
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V13(i,h(i));
    elseif NT(i) == 14
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V14(i,h(i));
    elseif NT(i) == 15
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V15(i,h(i));
    elseif NT(i) == 16
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V16(i,h(i));
    elseif NT(i) == 17
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V17(i,h(i));
    elseif NT(i) == 18
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V18(i,h(i));
    elseif NT(i) == 19
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V19(i,h(i));
    elseif NT(i) == 20
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V20(i,h(i));
    elseif NT(i) == 21
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V21(i,h(i));
    elseif NT(i) == 22
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V22(i,h(i));
    elseif NT(i) == 23
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V23(i,h(i));
    elseif NT(i) == 24
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V24(i,h(i));
    elseif NT(i) == 25
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V25(i,h(i));
    elseif NT(i) == 26
                [J(i,i),J(i+n,1),J(i-n,1),J(i,i-1),J(i,i+1)]=V26(i,h(i));
    end   
end
%Top right corner
        [J(i,i),J(i-n,1),J(i,i-1)]=V27(n*m,h(n*m));

end