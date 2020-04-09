function cur(n1,n2,val)
global b;

if (n1~=0) 
    b(n1,1) = b(n1,1)+ val;
end
if (n2~=0)
    b(n2,1) = b(n2,1) + val;
end
 if (n1 ~= 0) && (n2 ~= 0)
    b(n1,1) = b(n1,1) + val;
    b(n2,2) = b(n2,1) - val;
    
 end
 %End