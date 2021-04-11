import("stdfaust.lib");

zc(x) = x * x' < 0; 
up(x) = diff(x) > 0; 
down(x) = diff(x) < 0; 
diff(x) = x - x'; 

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

ibuffer(r_idx, x) = rwtable(size, .0, W_idx, x, int(ma.modulo(r_idx, size))); 
fbuffer(r_idx, x) = frwtable(5, size, .0, W_idx, x, r_idx);

grains(pos_0, dur_0, ptc_0, x_0) =  loop ~ 
                                    _
    with {
        loop(y) =   grain , 
                    lagrain
            with {
                grain = fbuffer(offset + line, x);
                t = loop ~ 
                    _ 
                    with {
                        loop(reset) = 
                            (fi.pole(1 - reset, 1) >= 
                                ba.sAndH(1 - 1' + reset, dur)) & zc(y);
                    };
                pos = pos_0 * ma.SR;
                dur = dur_0 * ma.SR;
                ptc = ba.sAndH(1 - 1' + t, ba.if(   ptc_0 < 0, 
                                                    min(-1 / 16, ptc_0), 
                                                    max(1 / 16, ptc_0)));
                x = +(r * x_0) ~ 
                    (de.delay(size - 1, size - 1) * (1 - r))
                    with {
                        r = checkbox("[00]Write on buffer");
                    };
                line = fi.pole(1 - t, 1 - t) * ptc; 
                offset = ba.sAndH(t, zc_sel + corr) 
                    with {
                        dir = ma.signum(ptc);
                        zc_sel = ba.if( diff(y) * dir > 0, 
                                        zc_up(pos, x), 
                                        zc_down(pos, x));
                        zc_up(read, x) = 
                            ibuffer(read, ba.sAndH(zc(x) & up(x), W_idx)); 
                        zc_down(read, x) = 
                            ibuffer(read, ba.sAndH(zc(x) & down(x), W_idx)); 
                        corr = y_diff / safe_den(x_diff) + (dir - 1) / 2
                            with {
                                y_diff = diff(y);
                                x_diff = ibuffer(zc_sel, diff(x));
                                safe_den(den) = ba.if(  den < 0,
                                                        min(0 - ma.EPSILON, den),
                                                        max(ma.EPSILON, den));
                            };
                    };
                lagrain = ba.if(idx < L, 
                                lagrangeN(N, x_vals, idx, y_vals), 
                                grain)
                    with {
                        N = 5;
                        L = 16;
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
                                        fbuffer(offset + (L + i) * ptc, x));
                            };
                        idx =   min(L, +(1 - t)) ~ 
                                *(1 - t);
                    };
            };
    };

//in = os.osc(1000);
secs = 5; // buffer size in seconds
SR = 192000;
size = SR * secs; // buffer size in samples
W_idx = ba.period(size); // writing pointer
position = no.noise * hslider("[01]Position jitter", secs, 0, secs, .000001) + hslider("Buffer position (secs)", 0, 0, secs, .000001);
duration = hslider("[02]Grain duration (secs)", .1, 0, secs, .000001); 
pitch = hslider("[03]Grain pitch (factor)", 1, -16, 16, .000001);
vol1 = hslider("[04]Output volume (linear)", 0, 0, 1, .000001);

process(in) = grains(position, duration, pitch, in) : *(vol1) , 
                                                    *(vol1);

// grain1(len, pos, pitch, x) = loop ~ _
//     with {
//         loop(y) = buffer2(offset + line, x)
//             with {
//                 t = loop ~ _ // trigger function; condition: grain dur. passed AND output at a ZC
//                     with {
//                         loop(reset) = (fi.pole(1 - reset, 1) >= ba.sAndH(1 - 1' + reset, len)) & zc(y);
//                     };
//                 pitch_sel = ba.sAndH(t, ba.if(checkbox("ada"), corr, pitch));
//                 line = fi.pole(1 - t, 1 - t) * pitch_sel; // line function: from 0 to at least len - 1
//                 offset = ba.sAndH(t, zc_sel + corr); // grain starting position
//                 //    with {
//                         dir = ma.signum(pitch);
//                         zc_sel = ba.if(diff(y) * dir > 0, zc_up(pos, x), zc_down(pos, x));
//                         //zc_slopes_up(N, read, x) = 
//                         zc_up(read, x) = buffer1(read, ba.sAndH(zc(x) & up(x), index)); // buffer containing ZC position of positive slopes
//                         zc_down(read, x) = buffer1(read, ba.sAndH(zc(x) & down(x), index)); // buffer containing ZC position of negative slopes
//                         corr = y_diff / safe_den(x_diff) + (dir - 1) / 2 : inspect(0, -1000, 1000) // correction for zero-order continuity
//                             with {
//                                 y_diff = diff(y);
//                                 x_diff = buffer1(zc_sel, diff(x)); // buffer containing the first derivative of the input signal
//                                 safe_den(den) = ba.if(  den < 0, 
//                                                         min(0 - ma.EPSILON, den), 
//                                                         max(ma.EPSILON, den));
//                             };
//                   //  };
//                 
//             };
//     };
// 
// grain2(len, pos, pitch, x) = loop ~ _ : ! , _
//     with {
//         loop(y) = buffer2(offset + line, x) , crossfade(N, buffer2(offset + line, x) @ (N), buffer2(offset - (N) + line, x)) 
//             with {
//                 N = hslider("xfade", 8, 8, 1024, 1); // number of crossfading samples at the junction corresponding to the lookahead delay
//                 t = loop ~ _ // trigger function; condition: grain dur. passed AND output at a ZC
//                     with {
//                         loop(reset) = (fi.pole(1 - reset, 1) >= ba.sAndH(1 - 1' + reset, len)) & zc(y);
//                     };
//                 line = fi.pole(1 - t, 1 - t) * pitch; // line function: from 0 to at least len - 1
//                 crossfade(N, x1, x2) = x1 * (1 - line) + x2 * line
//                     with {
//                         line = ((+(1 - t) : min(N - 1)) ~ *(1 - t)) / (N - 1);
//                     };
//                 offset = ba.sAndH(t, zc_sel + corr) // grain starting position
//                     with {
//                         zc_sel = ba.if(diff(y) > 0, zc_up(pos, x), zc_down(pos, x));
//                         zc_up(read, x) = buffer1(read, ba.sAndH(zc(x) & up(x), index)); // buffer containing ZC position of positive slopes
//                         zc_down(read, x) = buffer1(read, ba.sAndH(zc(x) & down(x), index)); // buffer containing ZC position of negative slopes
//                         //corr = y / safe_den(buffer(zc_sel, x));
//                         corr = diff(y) / safe_den(xslope(zc_sel, x)); // possible correction for varying slopes
//                         //corr = 1; // basic correction to avoid zero-order hold at the junction
//                         xslope(read, x) = buffer1(read, diff(x)); // buffer containing the first derivative of the input signal
//                         safe_den(den) = ba.if(den < 0, min(0 - ma.EPSILON, den), max(ma.EPSILON, den));
//                     };
// 
//             };
//     };

