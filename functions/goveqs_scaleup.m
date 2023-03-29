function out = goveqs_scaleup(t, in, M0, M1, times, i, s,r, growth, sel, agg,hivpoints)
   
scale = min((t-times(1))/(times(2)-times(1)),1); 
    if (scale<0),scale=0;end
 
    Mt = M1;
    Mt.lin = M0.lin + scale*(M1.lin-M0.lin);
    Mt.linh = M0.linh + scale*(M1.linh-M0.linh);
    Mt.ptcomp= M0.ptcomp + scale*(M1.ptcomp-M0.ptcomp) ;
    Mt.ptforg= M0.ptforg + scale*(M1.ptforg-M0.ptforg) ;
 
out = goveqs_basis(t, in, Mt, i, s,r, sel, agg, growth,hivpoints);
