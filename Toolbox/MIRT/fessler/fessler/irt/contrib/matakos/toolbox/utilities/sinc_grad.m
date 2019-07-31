function y = sinc_grad(x,w)
%| function y = sinc_grad(x,w)
%|
%| Function that calculated the derivative of a sinc function, where
%| sinc(w*x) = sin(pi*w*x)/(pi*w*x)
%|
%| in:
%|	x	[N 1]	Vector of values where derivative of sinc is evaluated
%|	w	scalar	width of the sinc function (scaling parameter)
%|
%| out:
%|	y	[N 1]	The derivative of the sinc function evaluated at locations
%|				specified by x
%|
%| Copyright 10-06-02, Antonis Matakos, University of Michigan

idx = x==0;

y = cos(pi*w*x)./x - sin(pi*w*x)./(pi*w*x.^2);
y(idx) = 0;
