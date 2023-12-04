function grad = gradient_lobe(ga_des,dt,gmax,smax,g_off,verb)
% This scripts calculates gradient lobe with a initial gradient value 
% for desired area. Results in triangluar or trapezoidal gradient shape.
% INPUT
%   grad = gradient_lobe(ga_des,dt,gmax,smax,g_off)
%  ga_des  Desired area of gradient lobe         [s*T/m]
%      dt  Sampling dwell time                   [s]
%    gmax  Maximum gradient strength             [T/m]
%    smax  Maximum slew rate                     [T/m/s]
%   g_off  Offset gradient value (default=0)     [T/m]
%
% OUTPUT
%    grad  Gradient waveform                     [T/m]
%
% Created 7/2018  Rolf Schulte
% Modified 7/2019 Enlin Qian
if (nargin<1), help(mfilename); return; end

if ~exist('g_off','var'), g_off = []; end
if isempty(g_off),        g_off = 0; end; g_off = round(g_off,2);
if ~exist('verb','var'),  verb = []; end
if isempty(verb),         verb = false; end

if gmax<0, error('gmax(=%g)<0',gmax); end
if smax<0, error('smax(=%g)<0',smax); end
if abs(g_off)>gmax, error('abs(g_off)(=%g)>gmax(=%g)',abs(g_off),gmax); end


%% calculate positive lobes; invert when necessary
sgn    = sign(ga_des);
ga_des = abs(ga_des);
g_off  = sgn*g_off;

gm = sqrt(smax*abs(ga_des)+0.5*g_off^2);    % nominal gradient strength
if gm<g_off        % if offset larger, invert again for lobe pointing up
     sgn = -sgn;
     ga_des = -ga_des;
     g_off = -g_off;
     gm = sqrt(-smax*abs(ga_des)+0.5*g_off^2);
end
if verb, fprintf('gm=%g; g_off=%g\n',gm,g_off); end


%% differentiate between trapezoids and triangles
if gm>gmax    % trapezoid: add gradient plateau in middle
    if verb, fprintf('Constructing trapezoid\n'); end
    
    % points for ramp up and down
    n_up   = ceil((gmax-g_off)/smax/dt)+1;
    n_down = ceil(gmax/smax/dt)+1;
    
    % area for ramp up and down
    ga_up   = 0.5*(gmax^2-g_off^2)/smax;
    ga_down = 0.5*gmax^2/smax;
    
    n_plt = ceil((ga_des-(ga_up+ga_down))/gmax/dt);
    % gmax_act = gmax;
else          % triangular gradient
    if verb, fprintf('Constructing triangle\n'); end
    n_up   = ceil((gm-g_off)/smax/dt);
    if n_up==1
        warning('n_up==1; setting to 0');
        n_up = 0;
    else
        n_up = n_up+1;
    end
    n_down = ceil(gm/smax/dt)+1;
    % gmax_act = gm;
    
    if n_up<0, warning('n_up<0'); end
    
    n_plt = 0;
end

%% calculate exact gmax_act to prevent rounding errors
gmax_act = (2*ga_des/dt-n_up*g_off)/(n_up+n_down+2*n_plt);

%% construct actual gradient waveform
grad = sgn*[linspace(g_off,gmax_act,n_up),...
    gmax_act*ones(1,n_plt),...
    linspace(gmax_act,0,n_down)];

%% check waveform
ga_act = sum(grad)*dt;
ga_err = (sgn*ga_des-ga_act)/ga_des;
smax_act = max(abs(diff(grad)))/dt;
if smax_act>smax, warning('smax_act(=%g)>smax(=%g)',smax_act,smax); end
if gmax_act>gmax, warning('gmax_act(=%g)>gmax(=%g)',gmax_act,gmax); end
if abs(ga_err)>1d-6
    warning('abs(ga_err(=%g))>1d-6; ga_des=%g; ga_act=%g',...
        ga_err,ga_des,ga_act); 
end

%% print info
if verb,
    fprintf('gmax=%g; gmax_act=%g [T/m]\n',gmax,gmax_act);
    fprintf('smax=%g; smax_act=%g [T/m]\n',smax,smax_act);    
    fprintf('n_up=%g; n_plt=%g; n_down=%g\n',n_up,n_plt,n_down);
    fprintf('ga_des=%g; ga_act=%g; err=%g[%%]\n',...
        sgn*ga_des,ga_act,ga_err*100);
end

end   % end main function gradient_lobe.m
