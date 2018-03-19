j = 0;
m = pi;
n = 1;
f = (m-j)/(n+1);
x = linspace(j,m,n+2);
y = linspace(j,m,n+2);

[k,l] = meshgrid(x,y);
a = k';
b = l';

c = 2:n+1;
d = 2:n+1;
a_int = a(c,d);
b_int = b(c,d);

f_n = @(x,y)-2*n*sin(n*x)*cosh(n*y);

right = f_n(a_int,b_int);

t_sol1 = (m-y)*sin(n*x)*sinh(n*y); % true solution

t_sol = t_sol1;

right(:,1) = right(:,1) - t_sol(c,1)/f^2;  % adjust right hand side to include boundary conditions
right(:,n) = right(:,m) - t_sol(c,n+2)/f^2;
right(1,:) = right(1,:) - t_sol(1,d)/f^2;
right(m,:) = right(n,:) - t_sol(n+2,d)/f^2;

g = reshape(rights,n*n,1);

R = speye(n);
E = ones(n,1);
O = spdiags([E-4*E E], [-1 0 1],n,n);
Z = spdiags([E E], [-1 1],n,n);
S = (kron(R,O) + kron(Z,R))/f^2;

sol = S\g; % solution to linear system

t_sol(c,d) = reshape(sol,n,n);

err = max(max(abs(t_sol - t_sol1)));

