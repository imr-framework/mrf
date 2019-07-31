function y = mrdivide(a,b)
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk

y = a;

y.eigblurmatrix = (a.eigblurmatrix) .^(b);
