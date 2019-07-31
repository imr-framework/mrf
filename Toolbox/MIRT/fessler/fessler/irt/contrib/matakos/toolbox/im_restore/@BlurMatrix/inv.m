function y = inv(a)
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk

y = a;

switch a.boundarycond
    case {'cir','refl'}
        y.eigblurmatrix = (a.eigblurmatrix) .^(-1);
    otherwise
        disp('Unknow method');
end
