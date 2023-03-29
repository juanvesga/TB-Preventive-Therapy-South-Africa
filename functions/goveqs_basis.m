function [out, lam] = goveqs_basis(t, in, M, i, s, r, sel, agg, growth,hivpoints)

invec = in(1:i.nstates);
N = sum(invec(s.nstates));

% -- Get force of infection
lam = (M.lambda*invec/N).*r.cohort + (1-r.cohort)*r.lam;


% -- If year is > 1980 get the HIV acquiaiton rate into the model 
foi=0;
if t>=1980
    foi = get_foi_hiv(t,r,hivpoints);
end


% -- Get all model components together
allmat = M.lin + M.linh +...
    lam(1) * M.nlin.ds  +...
    lam(2) * M.nlin.mdr +...
    foi * M.nlinh;

out = allmat*invec;


% Implement deaths
morts = M.mortvec.*invec;
out = out - morts;


% Implement births
births = sum(morts)*(growth==0) + growth;
out(i.U.hn) = out(i.U.hn)+births*r.cohort ;


% Get the auxiliaries
if r.cohort==1
out(i.aux.inc(1:4)) = agg.inc*(sel.inc.*allmat)*invec;
out(i.aux.inc(5))   = sum((sel.acqu.*allmat)*invec);
else
out(i.aux.inc(1:2)) = agg.inc_c*(sel.inc.*allmat)*invec;    
end

out(i.aux.notif)    = (agg.notif*(sel.notif.*allmat)*invec);
out(i.aux.mort)     = agg.mort*morts;
out(i.aux.pt)       = agg.ipt *(sel.ipt.*(allmat))*invec;
out(i.aux.newart)   = agg.art *(sel.art.*allmat)*invec;
out(i.aux.remote) = agg.remote*(sel.remote.*allmat)*invec;
out(i.aux.recent) = agg.recent*(sel.recent.*allmat)*invec;

%  out(i.aux.newart(1))= agg.art(1,:)*(sel.art.*allmat)*invec;
if (t>2017)
%  out(i.aux.newart(2))= sum((sel.ltbi.*allmat)*invec);
out(i.aux.daly(1:2))  = agg.daly*invec;
out(i.aux.daly(3:4))  = agg.daly*morts;
out(i.aux.ptpmo)      = agg.ptpmo*invec;
out(i.aux.ptcomp)     = agg.ptcomp*(sel.ptcomp.*M.ptcomp)*invec;
out(i.aux.ptforg)     = agg.ptcomp*(sel.ptcomp.*M.ptforg)*invec;
else
% out(i.aux.newart(1))=0;
out(i.aux.daly(1:2))  = 0;
out(i.aux.daly(3:4))  = 0;
out(i.aux.ptpmo)      = 0;
out(i.aux.ptcomp)     = 0;
out(i.aux.ptforg)     = 0;
end    
    
    



