function ob = ctranspose(ob)
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk


switch ob.boundarycond
    case {'cir', 'refl'}
        ob.eigblurmatrix = conj(ob.eigblurmatrix);
    case 'imfilter'
        ob.eigblurmatrix = flipud(fliplr(ob.eigblurmatrix));
    otherwise
        disp('Unknow method');
end
        
