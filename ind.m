function ind(n1,n2,val)
global G C b;
d = size(G,1); % current size of the MNA
xr = d+1;      % new row/column
b(xr) = 0;     % add new row
% Matlab automatically increases the size of a matrix
% if you use an index that is bigger than the current size.
G(xr,xr) = 0; % add new row/column
C(xr,xr) = 0;
if (n1 ~= 0)
    G(n1,xr) = 1;
    G(xr,n1) = 1;
end
if (n2 ~= 0)
    G(n2,xr) = 1;
    G(xr,n2) = 1;
end
if (n1 ~= 0) && (n2 ~= 0)
    G(n1,xr) = 1;
    G(xr,n1) = 1;
    G(n2,xr) = -1;
    G(xr,n2) = -1; 
end

C(xr,xr) = -val;
%End

