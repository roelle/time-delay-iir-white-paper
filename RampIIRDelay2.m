sysc = tf(1, [1, sqrt(2), 1]);
Ts = 1e-2
[B, A] = butter(2, 0.5*Ts/2);
sysd = tf(B, A, Ts);

tdc = ramp_delay(sysc);
tdd = ramp_delay(sysd);

t = (0:0.01:10)';
f = (0:0.0001:10)';
[Md, Pd] = bode(sysd, f*2*pi);
[Mc, Pc] = bode(sysc, f*2*pi);
M = [squeeze(Mc), squeeze(Md)];
P = [squeeze(Pc), squeeze(Pd)];
fig = 1;

figure(fig)
m = 3;
n = 1;
ii = 0;
a = zeros(m*n, 1);
h = cell(m*n, 1);

ii = ii + 1;
a(ii, 1) = subplot(m, n, ii);
h{ii, 1} = semilogx(f, 20*log10(M));
grid('on');

ii = ii + 1;
a(ii, 1) = subplot(m, n, ii);
h{ii, 1} = semilogx(f, P);
grid('on');

ii = ii + 1;
a(ii, 1) = subplot(m, n, ii);
h{ii, 1} = semilogx(f, -P/360./f, ...
                    10.^[-3, 3], tdc*[1, 1], 'b--', ...
                    10.^[-3, 3], tdd*[1, 1], 'r--');
grid('on');

linkaxes(a, 'x');
set(a, 'XTick', logticks());
set(a(2, 1), 'YTick', -360:45:0);
xlim([0.02, 5]);
ylim(a(1, 1), [-60, 0]);
ylim(a(2, 1), [-180, 0]);
xlabel(a(end, 1), 'Frequency [Hz]');
ylabel(a(1, 1), 'Magnitude [dB]');
ylabel(a(2, 1), 'Phase [degrees]');
ylabel(a(3, 1), 'Time Delay [s]');

s = sym('s');
z = sym('z');
t_sym = poly2sym(sysc.num{1}, s)/poly2sym(sysc.den{1}, s);
d_sym = vpa(poly2sym(sysd.num{1}, s)/poly2sym(sysd.den{1}, z),4);

legend(a(1, 1), {sprintf('$%s$',latex(t_sym))
        sprintf('$%s$',latex(d_sym))}, ...
        'Location', 'SouthWest', ...
        'Interpreter', 'latex')
legend(a(3, 1), h{3, 1}(3:4), ...
       'Continuous', 'Discrete', ...
       'Location', 'SouthWest')
pdfplot(fig, 'RampIIRDelay2', [], [5, 6.5]);

