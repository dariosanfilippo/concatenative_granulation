import("stdfaust.lib");

inspect(i, lower, upper) =
      _ <: _ ,
           vbargraph("sig_%i [style:numerical]", lower, upper) : attach;

lagrange_h(N, idx) = par(n, N + 1, prod(k, N + 1, f(n, k)))
    with {
        f(n, k) = ((idx - k) * (n != k) + (n == k)) / ((n - k) + (n == k));
    };

lagrangeN(N, idx) = lagrange_h(N, idx) ,
                    si.bus(N + 1) : si.dot(N + 1);

frwtable(N, S, init, w_idx, x, r_idx) =
    lagrangeN(N, f_idx, par(i, N+1, table(i_idx - int(N / 2) + i)))
    with {
        table(j) = rwtable(S, init, w_idx, x, int(ma.modulo(j, S)));
        f_idx = ma.frac(r_idx) + int(N / 2);
        i_idx = int(r_idx);
    };

size = 192000 * 10; // buffer size in samples
index = ba.period(size); // writing pointer
buffer1(read, x) = rwtable(size, .0, index, x, int(ma.modulo(read, size))); // buffer with wrapped-around reading pointer
buffer2(read, x) = frwtable(N, size, .0, index, x, read)
    with {
        N = 5; // Lagrange interp. order
    };

zc(x) = x * x' < 0; // zero-crossing indicator
up(x) = diff(x) > 0; // positive slope indicator
down(x) = diff(x) < 0; // negative slope indicator
diff(x) = x - x'; // first difference

grain1(len, pos, pitch, x) = loop ~ _
    with {
        loop(y) = buffer2(offset + line, x)
            with {
                t = loop ~ _ // trigger function; condition: grain dur. passed AND output at a ZC
                    with {
                        loop(reset) = (fi.pole(1 - reset, 1) >= ba.sAndH(1 - 1' + reset, len)) & zc(y);
                    };
                pitch_ada = ba.sAndH(t, corr);
                line = fi.pole(1 - t, 1 - t) * ba.if(checkbox("ada"), pitch_ada, pitch); // line function: from 0 to at least len - 1
                offset = ba.sAndH(t, zc_sel + corr); // grain starting position
                //    with {
                        dir = ma.signum(pitch);
                        zc_sel = ba.if(diff(y) * dir > 0, zc_up(pos, x), zc_down(pos, x));
                        //zc_slopes_up(N, read, x) = 
                        zc_up(read, x) = buffer1(read, ba.sAndH(zc(x) & up(x), index)); // buffer containing ZC position of positive slopes
                        zc_down(read, x) = buffer1(read, ba.sAndH(zc(x) & down(x), index)); // buffer containing ZC position of negative slopes
                        corr = y_diff / safe_den(x_diff) + (dir - 1) / 2 : inspect(0, -1000, 1000) // correction for zero-order continuity
                            with {
                                y_diff = diff(y);
                                x_diff = xslope(zc_sel, x);
                                xslope(read, x) = buffer1(read, diff(x)); // buffer containing the first derivative of the input signal
                                safe_den(den) = ba.if(  den < 0, 
                                                        min(0 - ma.EPSILON, den), 
                                                        max(ma.EPSILON, den));
                            };
                  //  };
                
            };
    };

grain2(len, pos, pitch, x) = loop ~ _ : ! , _
    with {
        loop(y) = buffer2(offset + line, x) , crossfade(N, buffer2(offset + line, x) @ (N), buffer2(offset - (N) + line, x)) 
            with {
                N = hslider("xfade", 8, 8, 1024, 1); // number of crossfading samples at the junction corresponding to the lookahead delay
                t = loop ~ _ // trigger function; condition: grain dur. passed AND output at a ZC
                    with {
                        loop(reset) = (fi.pole(1 - reset, 1) >= ba.sAndH(1 - 1' + reset, len)) & zc(y);
                    };
                line = fi.pole(1 - t, 1 - t) * pitch; // line function: from 0 to at least len - 1
                crossfade(N, x1, x2) = x1 * (1 - line) + x2 * line
                    with {
                        line = ((+(1 - t) : min(N - 1)) ~ *(1 - t)) / (N - 1);
                    };
                offset = ba.sAndH(t, zc_sel + corr) // grain starting position
                    with {
                        zc_sel = ba.if(diff(y) > 0, zc_up(pos, x), zc_down(pos, x));
                        zc_up(read, x) = buffer1(read, ba.sAndH(zc(x) & up(x), index)); // buffer containing ZC position of positive slopes
                        zc_down(read, x) = buffer1(read, ba.sAndH(zc(x) & down(x), index)); // buffer containing ZC position of negative slopes
                        //corr = y / safe_den(buffer(zc_sel, x));
                        corr = diff(y) / safe_den(xslope(zc_sel, x)); // possible correction for varying slopes
                        //corr = 1; // basic correction to avoid zero-order hold at the junction
                        xslope(read, x) = buffer1(read, diff(x)); // buffer containing the first derivative of the input signal
                        safe_den(den) = ba.if(den < 0, min(0 - ma.EPSILON, den), max(ma.EPSILON, den));
                    };

            };
    };

grain3(len, pos, pitch, x) = loop ~ _ : ! , _
    with {
        loop(y) = buffer2(offset + line, x) , interpolate(9, 4)
            with {
                t = loop ~ _ // trigger function; condition: grain dur. passed AND output at a ZC
                    with {
                        loop(reset) = (fi.pole(1 - reset, 1) >= ba.sAndH(1 - 1' + reset, len)) & zc(y);
                    };
                line = fi.pole(1 - t, 1 - t) * pitch; // line function: from 0 to at least len - 1
                interpolate(N, P) = 
                    ba.if(  lline < 1, 
                            lagrangeN(N, int(N / 2) + lline, points), 
                            buffer2(offset + line, x) @ hslider("del", 0, 0, P + 1, 1))
                    with {
                        points = l_points , r_points
                            with {
                                l_points = 
                                    par(i, (N + 1) / 2, ba.sAndH(t, y @ ((N + 1) / 2 - 1 - i)));
                                r_points = 
                                    par(i, (N + 1) / 2, buffer2(offset + P + i, x));
                            };
                        lline = ((+(1) : min(P + 1)) ~ *(1 - t)) / (P + 1);
                    };
                offset = ba.sAndH(t, zc_sel + corr) // grain starting position
                    with {
                        zc_sel = ba.if(diff(y) > 0, zc_up(pos, x), zc_down(pos, x));
                        zc_up(read, x) = buffer1(read, ba.sAndH(zc(x) & up(x), index)); // buffer containing ZC position of positive slopes
                        zc_down(read, x) = buffer1(read, ba.sAndH(zc(x) & down(x), index)); // buffer containing ZC position of negative slopes
                        corr = diff(y) / safe_den(xslope(zc_sel, x)); // possible correction for varying slopes
                        xslope(read, x) = buffer1(read, diff(x)); // buffer containing the first derivative of the input signal
                        safe_den(den) = ba.if(den < 0, min(0 - ma.EPSILON, den), max(ma.EPSILON, den));
                    };

            };
    };

in = os.osc(2001);
pos = no.noise * hslider("rand", 0, 0, size, 1) + hslider("pos", 0, 0, size, 1);
pitch = hslider("pitch", 1, -16, 16, .001);
len = hslider("len", 100, 0, 192000, 1); // it should be more or equal to the length of the crossfade
vol1 = hslider("vol1", 0, 0, 1, .001);
vol2 = hslider("vol2", 0, 0, 1, .001);
vol3 = hslider("vol3", 0, 0, 1, .001);

process(x) =    grain1(len, pos, pitch, in) * vol1 , 
                grain1(len, pos, pitch, in) * vol1 
                //grain2(len, pos, pitch, in) * vol2 ,
                //grain3(len, pos, pitch, in) * vol3
    with {
        in = +(r * x) ~ (de.delay(size - 1, size - 1) * (1 - r));
        r = checkbox("rec");
    };
