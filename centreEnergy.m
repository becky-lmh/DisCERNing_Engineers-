%% Plot energy of collision data
% Prove validity of assumption that Neon gas is at rest

mP= 0.93;           % Units of GeV
mNe = 19;           % Units of GeV
EkP = 7000;         % Units of GeV
EkNe = linspace(0,100,1000);       % Units of GeV
sqrtS = sqrt(mP^2 + mNe^2 + 2*EkP*mP + 2.*EkNe.*mNe);
velocity = sqrt(EkNe*1.6e-19*1e8.*2./(20*1.6e-27));
line = ones(1,1000)*122;
plot(velocity,sqrtS)
hold on
plot(velocity, line)
set(gca,'FontSize',10)
xlabel("Velocity of Neon (m/s)")
ylabel("sqrt(s) (GeV)")
grid on
hold off