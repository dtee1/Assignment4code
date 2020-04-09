function vol(n1,n2,val)

global G b C;

d = size(G,1); % current size of the MNA
xr = d+1; % new row
b(xr) = 0; % add new row
G(xr,xr) = 0; % add new row/column
C(xr,xr) = 0; % add new row/column

if (n1 ~= 0)
    G(n1,xr) = 1;
    G(xr,n1) = 1;
end

if (n2 ~= 0)
    G(n2,xr) = -1;
    G(xr,n2) = -1;
end
b(xr) = val;
%End
