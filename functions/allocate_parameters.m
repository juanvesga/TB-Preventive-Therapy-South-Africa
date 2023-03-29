
function [r,p] = allocate_parameters(x,r,p,xi)

r.beta          =  x(xi.beta);
r.beta_mdr      =  x(xi.beta_mdr);
r.bRRhiv        =  x(xi.bRRhiv);
r.reinfhiv      =  x(xi.reinfhiv);
r.fast_react    =   x(xi.fast);
r.symp_del      =  x(xi.symp_del);
r.careseek      =  x(xi.careseek);
r.ARTrec        =  x(xi.ARTrec);
p.IPThiv        =  x(xi.IPThiv);
p.Dx(1:2)       =  x(xi.Dxnhiv);
p.Dx(3)         =  x(xi.Dxhiv);
p.Tx_init       =  x(xi.Txinit);
r.muTB          =  x(xi.muTB);
r.RRmortTBhiv   =  x(xi.RRmortTBhiv);
p.ntpcov        =  x(xi.ntpcov);
 


