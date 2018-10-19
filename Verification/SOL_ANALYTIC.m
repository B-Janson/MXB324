function [] = SOL_ANALYTIC(fig, T, phi_avg, phi_true)

figure(fig)
hold on
plot(T, phi_avg, 'b')
plot(T, phi_true, 'r')
legend('avg', 'true')
title('water content (avg vs analytic)')
xlabel('time (days)')
ylabel('water content (phi)')
hold off
drawnow

end
