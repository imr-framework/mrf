function x = peak_remove(x)

for shot = 1:size(x,2)
   temp_raw = squeeze(x(:,shot)); 
   thresh = abs(2*std(temp_raw));%need to check this for spiral artifact removal
   ind = find(abs(temp_raw) > thresh);
   
   ind = ind(ind>18); %first 10 points might be the ones corresponding to low frequencies 
   
   ind_retain = setdiff(1:length(temp_raw), ind);
   new_raw_r = interp1(ind_retain, real(temp_raw(ind_retain)),1:length(temp_raw)); 
   new_raw_i = interp1(ind_retain, imag(temp_raw(ind_retain)),1:length(temp_raw));
   new_raw = complex(new_raw_r, new_raw_i);
   
   x(:,shot) = new_raw;
    
    
end