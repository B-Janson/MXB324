function dh=Jsolv(J,F,DIM,PARAMS)



if strcmp(PARAMS.method,'Full')
    %Regular n*m x n*m problem
    dh=J\(-F);
    
elseif strcomp(PARAMS.method,'Column')
    
 n=DIM.n;
m=DIM.m;   
    
    
    
    
    
    
    
    
end

end