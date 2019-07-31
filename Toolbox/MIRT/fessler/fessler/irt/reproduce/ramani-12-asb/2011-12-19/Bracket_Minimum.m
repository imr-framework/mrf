function abcf = Bracket_Minimum(ptemp)
% 2011-11-30, Sathish Ramani, University of Michigan
% Load constants
EPSL = ptemp.Constants.EPSL;
GOLD = ptemp.Constants.GOLD;
GLIMIT = ptemp.Constants.GLIMIT;

CAA = ptemp.CAA;

%% ax and fa
if(isfield(ptemp, 'BRax'))
    BRax = ptemp.BRax;
else
    BRax = 0;
end
BRfa = real( max( CAA(:) ) / min( CAA(:) ) );

%% bx and fb
BRbx = ptemp.BRbx;
t1 = clock;
BRfb = computemeasure(BRbx, ptemp); % Evaluate the meas for the given method at bx

%% cx and fc
BRcx = (BRbx)+GOLD*(BRbx-BRax); % First guess for c. 
BRfc = computemeasure(BRcx, ptemp); % Evaluate the meas for the given method at cx

iBrk = 0;
while (BRfb > BRfc) % Keep returning here until we bracket.
    te = clock-t1;
    te = te(1)*0 + te(2)*0 + te(3)*86400 + te(4)*3600 + te(5)*60 + te(6);
         
    % Show current results
    disp(['             Itr: ', int2str(iBrk), ...
    '; Current Bracket: [', num2str(BRax),', ', num2str(BRbx), ', ', num2str(BRcx), ']; ', ...
    'meas: [', num2str(BRfa), ', ', num2str(BRfb), ', ', num2str(BRfc), ']; ', ...
    'Time = ', num2str(te)]);

    t1 = clock;
    BRr = (BRbx-BRax) * (BRfb-BRfc); % Compute u by parabolic extrapolation from a,b,c.
    BRq = (BRbx-BRcx) * (BRfb-BRfa);
    
    tmp0 = BRq-BRr;
    tmp1 = abs(tmp0);
    
    if(tmp1<EPSL)
        tmp1 = EPSL;
    end
    tmp0 = tmp0/abs(tmp0)*tmp1;

    BRu = (BRbx) - ((BRbx-BRcx) * BRq - (BRbx-BRax) * BRr)/(2.0*tmp0);

    BRulim = BRbx + GLIMIT * (BRcx - BRbx); % We won?t go farther than this. Test various possibilities:

    if((BRbx-BRu)*(BRu-BRcx) > 0.0)  % Parabolic u is between b and c: try it.
        BRfu = computemeasure(BRu, ptemp); % Evaluate the meas for the given method at u

        if(BRfu < BRfc) % Got a minimum between b and c.
            BRax = BRbx;
            BRbx = BRu;
            BRfa = BRfb;
            BRfb = BRfu;
            break;
        elseif(BRfu > BRfb) % Got a minimum between a and u.
            BRcx = BRu;
            BRfc = BRfu;
            break;
        end
        BRu = BRcx + GOLD*(BRcx-BRbx); % Parabolic fit was no use. Use default magnification.

        BRfu = computemeasure(BRu, ptemp); % Evaluate the meas for the given method at u
        
    elseif((BRcx-BRu)*(BRu-BRulim) > 0.0) % Parabolic fit is between c and its allowed limit.
        BRfu = computemeasure(BRu, ptemp); % Evaluate the meas for the given method at u
        if(BRfu < BRfc)
            BRbx = BRcx;
            BRcx = BRu;
            BRu = BRcx + GOLD*(BRcx-BRbx);

            BRfb = BRfc;
            BRfc = BRfu;
            BRfu = computemeasure(BRu, ptemp); % Evaluate the meas for the given method at u
        end
    elseif((BRu-BRulim)*(BRulim-BRcx) >=0.0) % Limit parabolic u to maximum allowed value.
        BRu = BRulim;
        BRfu = computemeasure(BRu, ptemp); % Evaluate the meas for the given method at u
    else % Reject parabolic u, use default magnification.
        BRu = BRcx + GOLD*(BRcx-BRbx);
        BRfu = computemeasure(BRu, ptemp); % Evaluate the meas for the given method at u
    end
    BRax = BRbx;
    BRbx = BRcx;
    BRcx = BRu; % Eliminate oldest point and continue.

    BRfa = BRfb;
    BRfb = BRfc;
    BRfc = BRfu;
    iBrk = iBrk + 1;
end
te = clock-t1;
te = te(1)*0 + te(2)*0 + te(3)*86400 + te(4)*3600 + te(5)*60 + te(6);

abcf(1) = BRax; abcf(2) = BRbx; abcf(3) = BRcx;
abcf(4) = BRfa; abcf(5) = BRfb; abcf(6) = BRfc;

disp(['             Itr: ', int2str(iBrk), ...
    '; Final Bracket: [', num2str(BRax),', ', num2str(BRbx), ', ', num2str(BRcx), ']; ', ...
    '; meas: [', num2str(BRfa), ', ', num2str(BRfb), ', ', num2str(BRfc), ']; ', ...
    'Time = ', num2str(te)]);