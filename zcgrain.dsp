import("stdfaust.lib");

lagrange_h(N, x_vals, idx) = par(n, N + 1, prod(k, N + 1, f(n, k)))
    with {
        vals(i) = ba.take(i + 1, x_vals);
        f(n, k) = ((idx - vals(k)) * (n != k) + (n == k)) / 
            ((vals(n) - vals(k)) + (n == k));
    };

lagrangeN(N, x_vals, idx) = lagrange_h(N, x_vals, idx) ,
                            si.bus(N + 1) : si.dot(N + 1);

frwtable(N, S, init, w_idx, x, r_idx) =
    lagrangeN(N, x_vals, f_idx, par(i, N + 1, y_vals(i_idx - int(N / 2) + i)))
    with {
        x_vals = par(i, N + 1, i);
        y_vals(j) = rwtable(S, init, w_idx, x, int(ma.modulo(j, S)));
        f_idx = ma.frac(r_idx) + int(N / 2);
        i_idx = int(r_idx);
    };

ibuffer(r_idx, x) = rwtable(size, .0, index, x, int(ma.modulo(r_idx, size))); 
fbuffer(r_idx, x) = frwtable(5, size, .0, index, x, r_idx);

zc(x) = x * x' < 0; 
up(x) = diff(x) > 0; 
down(x) = diff(x) < 0; 
diff(x) = x - x'; 

CGP(len, pos, pitch, x) =   loop ~ 
                            _ 
    with {
        loop(y) = grain , Lgrain , t
            with {
                grain = fbuffer(offset + line, x);
                t = loop ~ 
                    _ 
                    with {
                        loop(reset) = 
                            (fi.pole(1 - reset, 1) >= 
                                ba.sAndH(1 - 1' + reset, len)) & zc(y);
                    };
                pitch_sah = ba.sAndH(1 - 1' + t, pitch);
                line = fi.pole(1 - t, 1 - t) * pitch_sah; 
                offset = ba.sAndH(t, zc_sel + corr) 
                    with {
                        dir = ma.signum(pitch_sah);
                        zc_sel =    ba.if(  diff(y) * dir > 0, 
                                            zc_up(pos, x), 
                                            zc_down(pos, x));
                        zc_up(read, x) = 
                            ibuffer(read, ba.sAndH(zc(x) & up(x), index)); 
                        zc_down(read, x) = 
                            ibuffer(read, ba.sAndH(zc(x) & down(x), index)); 
                        corr = y_diff / safe_den(x_diff) + (dir - 1) / 2
                            with {
                                y_diff = diff(y);
                                x_diff = ibuffer(zc_sel, diff(x));
                                safe_den(den) = ba.if(  den < 0,
                                                        min(0 - ma.EPSILON, den),
                                                        max(ma.EPSILON, den));
                            };

                    };
                Lgrain = ba.if( lline < L, 
                                lagrangeN(N, x_vals, lline, y_vals), 
                                grain)
                    with {
                        N = 5;
                        L = hslider("L", 16, 4, 64, 1);
                        halfp = (N + 1) / 2;
                        x_vals = par(i, N + 1, (i - halfp) * 
                            (i < halfp) + (i + L - halfp) * (i >= halfp));
                        y_vals =    l_points , 
                                    r_points
                            with {
                                l_points = 
                                    par(i, halfp, 
                                        ba.sAndH(t, y @ (halfp - 1 - i)));
                                r_points = 
                                    par(i, halfp, 
                                        fbuffer(offset + (L + i) * pitch_sah, x));
                            };
                        lline = min(L, (+(1 - t)) ~ 
                                *(1 - t));
                    };
                
            };
    };

size = 192000 * 5; 
index = ba.period(size);
pos = no.noise * size;
pitch = hslider("pitch", 1, -16, 16, .001) + no.noise * 
    hslider("ptc_rand", 0, -16, 16, 1) <: ba.if(<(0), min(-1/16), max(1/16));
len = hslider("len", 100, 0, 192000, 1); 
vol1 = hslider("vol1", 0, 0, 1, .001);
in1 = os.osc(1000); 
mic(x) = +(r * x) ~ (de.delay(size - 1, size - 1) * (1 - r)) * in_sel
    with {
        r = checkbox("rec");
    };
in_sel = checkbox("in_sel");
process = CGP(len, pos, pitch, in1) : *(vol1) , *(vol1) , *(vol1);

