 function ob = minus(a, b)
%function c = plus(a, b)
% "plus" method for this class
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk

if isa(a,'BlurMatrix')
    ob = a;
    if isa(b,'BlurMatrix')
        c = a.eigblurmatrix - b.eigblurmatrix ;
    else
        if numel(b)~=1
            error('Wrong');
        else
            c = a.eigblurmatrix - b;
        end        
    end
else
    if numel(b)~=1
        error('Wrong');
    else
        ob = b;
        c = a - b.eigblurmatrix ;
    end
end
ob.eigblurmatrix = c;