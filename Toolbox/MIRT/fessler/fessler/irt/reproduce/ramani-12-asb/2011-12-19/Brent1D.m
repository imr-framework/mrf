function [minmeas_temp lamopt_temp] = Brent1D(abcf, params_temp)
% 2011-11-30, Sathish Ramani, University of Michigan
% abcf contains the variables for the bracketed minimum: a<b<c and fa<fb<fc

%% Load constants
EPSL = params_temp.Constants.EPSL;
CGOLD = params_temp.Constants.CGOLD;
Brenttol = params_temp.Constants.Brenttol;
maxBitr = params_temp.Constants.maxBitr;

Brente = 0.0;

Brenta = abcf(1); % ax
Brentb = abcf(3); % brentmaxlam*nv; //cx

Brentbx = abcf(2); % brentlammid*nv; // bx : ax<bx<cx such that f(ax)>f(bx)<f(cx)
measbx = abcf(5); % meas of the chosen method at bx

Brentx = Brentbx; 
Brentw = Brentbx; 
Brentv = Brentbx;

% Brent initialization (meas at bx, the mid point and the least meas so far)
Brentfw = measbx; 
Brentfv = measbx; 
Brentfx = measbx;

itrb = 0;

disp(['             Itrb=' int2str(itrb) '; stepsize = ' num2str(Brentbx) '; measVAL = ' num2str(measbx)]);

% Start BRENT ITERATIONS
Brentxm = Brentx; 
Brenttol2 = 0.0; 
Brentd = 0.0;
while((itrb<maxBitr) && (abs(Brentx-Brentxm)>=Brenttol2-0.5*(Brentb-Brenta))) % Main loop that tests for DONE here

        Brentxm = 0.5*(Brenta + Brentb);
        Brenttol1 = Brenttol*abs(Brentx) + EPSL;
        Brenttol2 = 2.0*Brenttol1;

        if(abs(Brente) > Brenttol1)
            % Construct a trial parabolic fit
            Brentr = (Brentx - Brentw)*(Brentfx - Brentfv);
            Brentq = (Brentx - Brentv)*(Brentfx - Brentfw);
            Brentp = (Brentx - Brentv)*Brentq - (Brentx - Brentw)*Brentr;

            Brentq = 2.0*(Brentq - Brentr);
            
            if(Brentq > 0.0) 
                Brentp = -Brentp;
            end

            Brentq = abs(Brentq); 
            Brentetemp = Brente; 
            Brente = Brentd;

            % The condition below determine the acceptability of the parabolic fit.
            if((abs(Brentp) >= abs(0.5*Brentq*Brentetemp)) || (Brentp <= Brentq*(Brenta - Brentx)) || Brentp >= Brentq*(Brentb - Brentx)) 
                % Here we take the golden section step  into the larger of the two segments.
                if(Brentx >= Brentxm)	
                    Brente = Brenta - Brentx;
                else
                    Brente = Brentb - Brentx;
                end
                Brentd = CGOLD*Brente;
            else % If good, take the parabolic step.
                Brentd = Brentp/Brentq;
                Brentu = Brentx + Brentd;
                if((Brentu - Brenta < Brenttol2) || (Brentb - Brentu < Brenttol2))
                    Brentd = Brenttol1*(Brentxm-Brentx)/abs(Brentxm-Brentx);
                end
            end
        else % If the trial parabolic fit was not good, use the Golden Ratio method
            if(Brentx >= Brentxm)
                Brente = Brenta - Brentx;
            else
                Brente = Brentb - Brentx;
            end
            Brentd = CGOLD*Brente;
        end

        if(abs(Brentd) >= Brenttol1)	
            Brentu = Brentx + Brentd;
        else
            Brentu = Brentx + Brenttol1*Brentd/abs(Brentd);
        end

        % Compute the meas for the chosen step size Brentu
        t1 = clock;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Brentfu = computemeasure(Brentu, params_temp);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        te = clock-t1;
        te = te(1)*0 + te(2)*0 + te(3)*86400 + te(4)*3600 + te(5)*60 + te(6);
        
%         disp('--------------------------------------------------------------------');
        disp(['             Itrb=' int2str(itrb+1) '; stepsize = ' num2str(Brentu) '; measVAL = ' num2str(Brentfu) '; Time = ' num2str(te)]);
        
        % Now decide what to do with our function evaluation.
        if(Brentfu <= Brentfx) % House keeping...
                if(Brentu >= Brentx)	
                    Brenta = Brentx;
                else
                    Brentb = Brentx;
                end

                Brentv = Brentw;	
                Brentfv = Brentfw;
                
                Brentw = Brentx;
                Brentfw = Brentfx;
                
                Brentx = Brentu;
                Brentfx = Brentfu;
        else
            if(Brentu < Brentx)	
                Brenta = Brentu;
            else
                Brentb = Brentu;
            end

            if((Brentfu <= Brentfw) || (Brentw == Brentx)) 
                Brentv = Brentw;	Brentfv = Brentfw;
                Brentw = Brentu;	Brentfw = Brentfu;
            elseif((Brentfu <= Brentfv) || (Brentv == Brentx) || (Brentv == Brentw))
                    Brentv = Brentu;	Brentfv = Brentfu;
            end
        end
        itrb = itrb + 1;
end
minmeas_temp = Brentfx; % TRUE meas
lamopt_temp = Brentx; % OPTIMUM STEP SIZE