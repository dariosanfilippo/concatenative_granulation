import("stdfaust.lib");

// =============================================================================
//      AUXILIARY FUNCTIONS
// =============================================================================

dl(del, in) = de.delay(size, del, in);
grain(del, in) = de.fdelayltv(4, size, del, in);
index = ba.period(size);
lowpass(cf, x) =    +(x * a0) 
                    ~ *(b1)
with {
    a0 = 1 - b1;
    b1 = exp(-cf * ma.PI);
};

// =============================================================================
//      MATH
// =============================================================================

diff(x) = x - x';
//div(x1, x2) = ba.if(x2 == 0, 0, x1 / x2);
div(x1, x2) = x1 / ba.if(  x2 < 0, 
                           min(ma.EPSILON * -1, x2), 
                           max(ma.EPSILON, x2));
line_reset(rate, reset) =  (rate / ma.SR)' : 
    (+ : *(1 - (reset != 0)))
    ~ _;
wrap(lower, upper, x) = (x - lower) / (upper - lower) : 
    ma.decimal : * (upper - lower) : + (lower);
zc(x) = x * x' < 0;

// =============================================================================
//      BUFFER SIZE (samples)
// =============================================================================

size = 2^20;

// =============================================================================
//      LIVE/LOOPED INPUT FUNCTION
// =============================================================================

input(x) =  +(x * rec)
            ~ de.delay(size - 1, size - 1) * (1 - rec);

// =============================================================================
//      GRANULATOR FUNCTION
// =============================================================================

grains_dl_zc(size) =    loop 
                        ~ _
    with {
        // MAIN LOOP
        loop(out, pitch1, rate1, position1, in) =
            ((ba.sAndH(trigger(out), zc_index(position, in, out)) + 
                shift(trigger(out))) -
                    ba.sAndH(trigger(out), corr(position, in, out)) : 
                        wrap(0, size) ,
            in : grain)
            with {
                // PARAMETERS SETUP
                pitch = ba.sAndH(trigger(out), pitch1);
                rate = abs(rate1);
                position = position1 : wrap(0, size);
                // TRIGGER FUNCTION
                trigger(y) =    loop
                                ~ _
                    with {
                        loop(ready) =   
                            zc(y) ,
                            line_reset(ba.sAndH(1 - 1' + ready, rate),
                                ready) >= 1 : &;
                    };
                // DIRECTION INVERSION
                dir = ma.signum(pitch);
                // READING HEAD FUNCTION
                shift(reset) = 
                    div((1 - pitch), rate) *
                        line_reset(ba.sAndH(reset, rate), reset) ^ 
                            p_mod * ma.SR;
                // ZC POSITION FUNCTION
                zc_index(recall, x, y) = 
                    index - 
                        ba.if(dir * diff(y) >= 0, zc_up, zc_down) : 
                            wrap(0, size)
                    with {
                        zc_up = recall , 
                                ba.sAndH(store, index) : dl
                            with {
                                store = zc(x) ,
                                        (diff(x) > 0) : &;
                            };
                        zc_down =   recall , 
                                    ba.sAndH(store, index) : dl
                            with {
                                store = zc(x) ,
                                        (diff(x) < 0) : &;
                            };
                    };
                // POSITION CORRECTION FUNCTION
                corr(recall, x, y) = div(y_diff, x_diff) + ((dir - 1) / 2)
                    with {
                        y_diff = diff(y);
                        x_diff =    zc_index(recall, x, y) , 
                                    diff(x) : dl;
                    };
            };
    };

// =============================================================================
//      INTERFACE SETTINGS
// =============================================================================

rec = checkbox("[0]Looped/Live");
vol = hslider("[8]Output level", 0, 0, 1, .001) ^ 2;
p = hslider("[2]Pitch factor", -1, -16, 16, .001);
p_mod = pow(16, hslider("[3]Pitch modulation", 0, -1, 1, .001));
r = hslider("[1]Grain rate", 20, 1, 1000, .001);
t = hslider("[4]Time factor", 1, -16, 16, .001);
region = hslider("[5]Buffer region", 0, 0, 1, .001) * size;
asynch(x) = asynch_amount * lowpass(asynch_degree, x);
asynch_amount = hslider("[6]Position self-modulation depth", 0, 0, 1, .001) ^ 
    2 * size * 16;
asynch_degree = hslider("[7]Position self-modulation rate", .5, 0, 1, .001) ^ 4;
pos(x) =    ((+(1 - t) : wrap(0, size)) 
            ~ _) + region + asynch(x);

// =============================================================================
//      MAIN FUNCTION
// =============================================================================

process(x) =    ((p, r, pos, input(x)) : grains_dl_zc(size)) 
                ~ _ <: par(i, 2, *(vol));
