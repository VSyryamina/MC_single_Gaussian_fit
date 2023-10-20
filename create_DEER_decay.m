function V_DEER = create_DEER_decay (r, f_r, t, pB, gAB, noise)

f_r = f_r/sum(f_r);
wv = 326 *gAB/ 2 / pi ./ r.^3;

V_DEER  = t*0; 
for i = 1:length(t)
 
    for j = 1:length(r)
                        
        for k = 0:1.5:90           
            st = sind (k); ct2 = (cosd(k))^2;     
            V_interim (i) = (1-pB* (1-(cos (2*pi*( wv(j) *(3*ct2-1) ) * t (i) ))))    * st * f_r(j);     
            V_DEER (i) = V_interim (i) + V_DEER(i);
        end
    end
end

zero_val = find(t==0);
V_DEER = V_DEER / max(V_DEER)+noise * rand (size (V_DEER));
