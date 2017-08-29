function [t1, J1, t2, J2, t3, J3, profile] = StopTrajGen(a0, v0, d, J, plot_traj)

%% Description
% Based on the initial acceleration (a0), velocity (v0), maximum
% deceleration (d) and a desired Jerk (J), compute a triangular or
% trapezoidal profile in order to stop a motion system.

v_star = d^2/J - a0^2/(2*J);
a02_2J_p = a0^2/(2*J);

if v0 <= v_star && v0 >= a02_2J_p
    a_lim = -sqrt(J*v0 + a0^2/2);
    t1 = -a_lim/J;
    J1 = -J;
    t2 = (a0 - a_lim)/J;
    J2 = J;
    t3 = 0;
    J3 = 0;
    profile = 0;
elseif v0 > v_star
    t1 = d/J;
    J1 = -J;
    t2 = (v0 - v_star)/d;
    J2 = 0;
    t3 = (d + a0)/J;
    J3 = J;
    profile = 1;

end

if plot_traj
    close all
    dt = 0.0001;
    t1 = round(t1/dt)*dt;
    t2 = round(t2/dt)*dt;
    t3 = round(t3/dt)*dt;
    if profile
        t = 0:dt:1.05*(t1 + t2 + t3);
        n = length(t);
        nt1 = t1/dt - 1;
        nt2 = t2/dt;
        nt3 = t3/dt;
        jt = [J1*ones(1,nt3) J2*ones(1,nt2) J3*ones(1,nt1) 0*ones(1,(n-nt1-nt2-nt3))];
    else
        t = 0:dt:1.05*(t1 + t2);
        n = length(t);
        nt1 = round(t1/dt - 1);
        nt2 = round(t2/dt);   
        jt = [J1*ones(1,nt2) J2*ones(1,nt1) 0*ones(1,(n-nt1-nt2))];
    end
    at = zeros(1,n);
    vt = zeros(1,n);
    at(1) = a0;
    vt(1) = v0;
    for i = 2:n
        at(i) = int_tustin(jt(i), jt(i-1), at(i-1), dt);
        vt(i) = int_tustin(at(i), at(i-1), vt(i-1), dt);
    end
    
    figure
    hold on
    plot(t,jt, 'linewidth', 2)
    title('Jerk')
    xlabel('Time [s]')
    ylabel('Jerk [m/s^3]')
    hold off
    
    figure
    hold on
    plot(t,at,'linewidth', 2)
    title('Acceleration')
    xlabel('Time [s]')
    ylabel('Acceleration [m/s^2]')
    hold off
   
    figure
    hold on
    plot(t, vt, 'linewidth', 2)
    title('Velocity')
    xlabel('Time [s]')
    ylabel('Velocity [m/s]')
    hold off   

end
    
    
end

