function ob = rdivide(a,b)
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk

if isa(a,'BlurMatrix')
    ob = a;
    if isa(b,'BlurMatrix')
        c = (a.eigblurmatrix) ./ (b.eigblurmatrix) ;
    else
        c = (a.eigblurmatrix) ./ b;
    end
else
    ob = b;
    c = a ./ (b.eigblurmatrix) ;
end
ob.eigblurmatrix = c;
return;
