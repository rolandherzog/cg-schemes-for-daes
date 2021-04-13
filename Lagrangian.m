function [val,der] = Lagrangian(x,xbar,i)
% This function evaluates the i-th 1D Lagrangian interpolation polynomial
% belonging to the set xbar (column vector) of interpolation points. The
% evaluation is done at the point(s) given in x (row vector). The derivative is
% also evaluated. The size of xbar determines the order of the polynomial.

% Remember and unify the dimensions of the inputs x and xbar
[m,n] = size(x);
x = x(:)';
xbar = xbar(:);

% Evaluate the function values of the Lagrangian
if (length(xbar) > 1)
	num = x - xbar;
	den = xbar(i) - xbar;
	num = [num(1:i-1,:); num(i+1:end,:)];
	den = [den(1:i-1); den(i+1:end)];
	val = prod(num,1) / prod(den);
else
	val = ones(size(x));
end

% Evaluate the derivative of the Lagrangian as well
der = zeros(size(val));
if (length(xbar) > 1)
	for k=1:size(num,1)
		der = der + prod([num(1:k-1,:); num(k+1:end,:)],1);
	end
	der = der / prod(den);
end

% Reshape the output to the dimensions of the input x
val = reshape(val,m,n);
der = reshape(der,m,n);

