  % Simulation of an intergrate-and-fire neuron

    % Yair Lahad [205493018]
    % Neriya Mizrahi [315711697]

    % Computation and cognition undergrad - ex2

    clear; close all; clc;
  %% Declare simulation parameters (we will work in units of Volts, Siemens, Farads and Sec)
    gL = 1e-4; % conductance (S/cm^2), 1/R
    C = 1e-6; % capacitance (F/cm^2)
    tau = C/gL; % sec - time constant, equal to RC
    TH_in_mV = 30; % mV
    TH = TH_in_mV/1000; % V
    tau_R = 0.004; % sec - Refractoriness

    dt = 0.0001; % time step for numerical integration
    t_final = 1; % sec, duration of numerical simulation
    n = round(t_final/dt); % number of iterations (steps)
    V = NaN(1,n+1); % initialize array for the voltage trace
    V(1) = 0; % initial condition for the voltage - This is V(0)
    %In matlab we arrays are 1-indexed (V(1) is the first element)
    %This is why we have 10001 columns in I

    t = (0:n)*dt; % time vector (n+1 components)

    %% load external current
    %TODO: After you implement 'the plotting section', load the currents
    %(one by one) and plot the applied current and the neuron's response

     load I_const; % loading I
    % load I_sin;
    % load I_exp;
    % load I_step;

    %% Numerical integration
    RP_flag = 0; %marks the refratory period of the neuron. if the cell is in not the refractory period =0,if it is=1.
    for idx = 1:n %loop running from 1-n, marking the index of each vector calculated below
        VV = V(idx) + (-V(idx)/tau + I(idx)/C)*dt; % voltage value equals the result of the dynamical equation of an RC circuit- Voltage and current (I) per index. 

        if VV >= TH %if the Voltage Value is bigger or equal to the Threshold (30mV)
            VV = 0; % make the Voltage Value 0
            RP_flag = 1; %put the neoron in the Refratory Period  1)
            t_RP_start = t(idx); %mark the time in which the refractory period began
        end
        if RP_flag %if the cell is in the refractory period
            VV = 0; % make the Voltage Value 0
            if (t(idx) - t_RP_start >= tau_R) %if the amount of time passed is either bigger or equal to the duration of the refratory period constant
                RP_flag = 0; % stop refratory period for the neuron
            end
        end

        V(idx+1) = VV; %insert the value calculated above into the vector V
    end

    %% The plotting section
    figure('Color','w');
    MicroA = I*10^(6);
    MiliV = V*10^(3);
    subplot(2,1,1);
    x = t;
    y1 = MicroA;
    plot(x,y1)
    set(gca,'FontSize',16)
    ylabel('Current [\muA]')
    title('current insearted is constant')
    subplot(2,1,2);
    y2 = MiliV;
    plot(x,y2)
    hold on
    %creating dashed (-.) line on the same plot- horizontally from 0-1 on
    %X-axis and at the peak (highest voltage, 30Mv) on the Y-axis:
    subplot(2,1,2);
    xTH = [0,1];
    yTH = [TH_in_mV,TH_in_mV];
    plot(xTH,yTH,'-.r')
    hold on;
    % Note: the simulation computes the voltage of the neuron relative to the
    % resting potential. Thus, the Y-axis starts from zero.
    ylim([0 32]);
    set(gca,'FontSize',16)
    xlabel('Time [S]');
    ylabel('Relative Voltage [mV]')
    
