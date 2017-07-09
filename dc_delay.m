function td = dc_delay(sys)
%DC_DELAY(SYS) calculate the DC (ramp) delay for a transfer function
%   DC_DELAY(SYS) calculates the delay for a ramp signal subject to a
%   unit-gain transfer function. The delay is constant for the transfer
%   function, so the sole argument is the transfer function, sys.
%
%   Usage:
%   	Ts = 1e-2
%       [B, A] = butter(2, 0.5*Ts/2);
%       sysd = tf(B, A, Ts);
%       td = dc_delay(sysd);
%       udot = 10;
%       u = t*udot;
%       t = (0:Ts:10);
%       y = filter(B, A, u); % u is a position at constant velocity
%       u_est = y + td*mean(diff(u)/Ts);
%       subplot(2, 1, 1)
%       plot(t, [u, y, u_est])
%       subplot(2, 1, 2)
%       plot(t, u - u_est)
%
%   See also: TF, DCGAIN 
%
% Copyright 2017 Matt Roelle

assert(isa(sys, 'tf'), 'System must be a transfer function.');
assert(length(sys.den) == 1, 'System must be a single transfer function');

A = sys.den{1};
m = length(A);
B = [zeros(1, m - length(sys.num{1})), sys.num{1}]; % Zero pad B
if sys.Ts > 0
    td = sys.Ts*(0:length(A) - 1)*(B-A)'/sum(A);
else
    td = (A(m-1) - B(m-1))/A(m);
end