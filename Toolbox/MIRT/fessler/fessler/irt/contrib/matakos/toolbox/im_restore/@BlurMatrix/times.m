function y = times(ob, x)
%function y = times(ob, x)
% y = G * x	 
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk
if ~isa(ob, 'BlurMatrix')
    if isa(x,'BlurMatrix')
        y = ob .* x.eigblurmatrix;	
    else
        error('Wrong');
    end
else   
    if isa(x,'BlurMatrix')
        error('Tow object, using H1 * H2.');
    end
    y = ob.eigblurmatrix .* x;	
end