%% Kinetic model of the mScarlet3/FMN/FMNH2 magnetic field effect.

mSc_init = 0.05; % initial conc. of mScarlet (µM)
FMN_init = 30; % initial conc. of FMN (µM)
FMNH2_init = 10; % initial conc. of FMNH2 (µM)
tspan = 0:0.01:100; % array of timesteps to simulate (seconds)

% Define initial concentrations of all chemical species
% Generally, the initial concentrations of all other tracked states are zero.
% 
% Column 1: mScarlet ground state
% Column 2: mScarlet triplet state
% Column 3: mScarlet dark radical anion
% Column 4: FMN
% Column 5: FMNH2
% Column 6: FMNH radical. 

y0 = [mSc_init, 0, 0, FMN_init, FMNH2_init, 0]; 

% Simulate and store the concentrations at each timestep. 
y = query_model(y0); 

% Plot results for mScarlet3 and mScarlet3 radical. 
figure(1)
hold on
ylim([0 mSc_init])

% Plot the grey bars corresponding to nonzero magnetic field.
for i = 3:2:9 
    plot_patch(10*i, 10*i+10, get(gca, 'ylim'), [.85 .85 .85])
end

plot(tspan, y(:,1)) % Conc. of ground state mScarlet over time
plot(tspan, y(:,3)) % Conc. of mSc radical anion over time 
legend("mSc_0", "mSc^{⚬-}")
xlabel("Time (s)")
ylabel("Concentration (µM)")

%%
function y = query_model(y0, params, tspan)
    % query_model(y0, params) runs the model given y0 (initial concentrations) 
    % and params (optional, model parameters). 
    % 
    % y0: initial concentrations in µM. 
    % [mSc, mSc triplet, mSc radical, FMN, FMNH2, FMNH radical] . 
    % 
    % Returns: 
    % y(t), concentrations over time. 
    %
    % params: model inputs. see below for definitions
    % lis = params(1); % light intensity * cross section [ /s]
    % ksr = params(2); % inverse of excited state lifetime [ /s]
    % ktr = params(3); % inverse of mScarlet triplet lifetime [ /s]
    % kisc = params(4); % mScarlet intersystem crossing rate [ /s]
    % kred = params(5); % rate of reduction of triplet mScarlet by FMNH2 [ /µM /s]
    % kox = params(6); % rate of dark mScarlet reoxidation to ground state [ /µM /s] 
    % kdis = params(7); % rate of free radical FMNH disproportionation [ /µM /s]
    % beta0 = params(8); % SCRP separation fraction at B = 0 [ dimensionless ] 
    % betaB = params(9); % SCRP separation fraction at B = Bsat [ dimensionless ]
    % 
    % beta0 and betaB can also be obtained via simulations, 
    % using the file "SCRP_spin_and_chem_sim.m"

    if ~exist('params', 'var') % default model parameter values
        params = [1652	2.5e8	224	185497.528691897 ... 
        	50.5406 ... 
            0.003675	5000	0.8075, 0.9323];
    end

    if ~exist('tspan', 'var') % default time is 0 to 100s. 
        tspan = 0:0.01:100; % seconds. 
    end
    [~,y] = ode23s(@(t,y) MFE_model_trackedFMN(t, y, params), tspan, y0');

end
%%
function dydt = MFE_model_trackedFMN(t,y,params)

  lis = params(1); % light intensity * cross section [ /s]
  ksr = params(2); % inverse of excited state lifetime [ /s]
  ktr = params(3); % inverse of mScarlet triplet lifetime [ /s]
  kisc = params(4); % mScarlet intersystem crossing rate [ /s]
  kred = params(5); % rate of reduction of triplet mScarlet by FMNH2 [ /µM /s]
  kox = params(6); % rate of dark mScarlet reoxidation to ground state [ /µM /s] 
  kdis = params(7); % rate of free radical FMNH disproportionation [ /µM /s]
  beta0 = params(8); % SCRP separation fraction at B = 0 [ dimensionless ] 
  betaB = params(9); % SCRP separation fraction at B = Bsat [ dimensionless ] 

  dydt = zeros(6,1);
  if t > 30
      if mod(t, 20) < 10 % magnet is off if t mod 20s < 10s
        beta = beta0;
      else
        beta = betaB; % magnet on
      end
  else
        beta = beta0;
  end
  dydt(1) = -(kisc*lis)/(ksr + kisc)*y(1) + (ktr + (1-beta)*kred)*y(5)*y(2) + kox*y(4)*y(3); % mSc3 ground state
  dydt(2) = +(kisc*lis)/(ksr + kisc)*y(1) -(ktr + kred)*y(5)*y(2); % mSc3 triplet
  dydt(3) = beta*kred*y(5)*y(2) - kox*y(4)*y(3); % mSc radical
  dydt(4) = -kox*y(4)*y(3) + kdis*y(6)^2; % FMN
  dydt(5) = -beta*kred*y(5)*y(2) + kdis*y(6)^2; % FMNH2
  dydt(6) = beta*kred*y(5)*y(2) -2*kdis*y(6)^2 + kox*y(4)*y(3); % FMNH radical
end

%%
function plot_patch(start_t,end_t,ylim,color)
    patch([start_t end_t end_t start_t],[min(ylim) min(ylim) max(ylim) max(ylim)], color,'EdgeColor', 'none', 'HandleVisibility', 'off')
    set(gca,'Layer','top')
end