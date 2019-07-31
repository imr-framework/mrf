function y = abs(a)
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk

y = a;
switch ob.boundarycond
    case {'cir','refl'}
        y.eigblurmatrix = abs(a.eigblurmatrix);
    otherwise
        disp('Unknow method');
end
