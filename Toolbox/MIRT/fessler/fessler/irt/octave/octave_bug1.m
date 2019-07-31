function fun1

a = {1, 2};
b = 'foo';

fun2(a{1}, a{2}, b) % works

fun2(a{:}, b) % fails


function fun2(varargin)

nargin % always is 3

disp(argn) % look at the 'a {:}' here the 2nd time!
numel(argn) % always 3

inputname(3)

% fails the second time it is called above with this error
% error: fun2: A(I): index out of bounds; value 3 out of bound 2

