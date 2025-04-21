%% Simulate SCRP spin dynamics with stochastic separation and recombination

nspins = 500; % number of spins to simulate in each triplet sublevel
dt = 0.10; % timestep (nanoseconds) 
nsteps = 500; % number of total timesteps
eta = 1; % ratio of SCRP recombination to separation

Bzero = 0; % zero Tesla
B30 = 0.03; % 30 milliTesla  

[pS_Bzero, beta_zero] = spin_chem_sim(nspins, nsteps, dt, eta, Bzero);
[pS_B30, beta_30] = spin_chem_sim(nspins, nsteps, dt, eta, B30);

%% Display the SCRP separation fraction at 0 and at 30 mT

disp(beta_zero)
disp(beta_30)

%% Plot ensemble-averaged singlet probability in zero and high magnetic field
figure(1)
hold on
tau = (1:nsteps)*dt; 
plot(tau, mean(pS_Bzero, 1), '--', 'Color', "#ff8c19");
plot(tau, mean(pS_B30, 1), '-', 'Color', "#ff8c19");
xlabel("Time (ns)")
ylabel("Ensemble-averaged \newline Singlet Probability")
l = legend("B_{ext} = 0 mT", "B_{ext} = 30 mT");
xlim([0 20])

%% Simulate individual SCRPs for different values of eta and magnetic field

nspins = 200; % number of SCRPs in each initial T0/T+/T- state
nsteps = 50000;
dt = 0.10; % ns 
% bz = (0:30:300)*10^(-4); % Tesla
bz = [0, 3000]*10^(-4); % Tesla, for zero and saturating magnetic field
eta = [0.1, 0.3, 1, 3, 10, 30, 100];

[~, beta_T0born] = spin_chem_sim_analytical(nspins, nsteps, dt, eta, bz, "T0");
[~, beta_Tplusborn] = spin_chem_sim_analytical(nspins, nsteps, dt, eta, bz, "Tpm");
[~, beta_Tminusborn] = spin_chem_sim_analytical(nspins, nsteps, dt, eta, bz, "Tpm");
beta_tot = beta_T0born/3 + beta_Tplusborn/3 + beta_Tminusborn/3; 

%% Plot the results of beta_sat vs. beta_0, as in Fig. 7

figure(1)
scatter(beta_tot(:,1), beta_tot(:,end), 30, log10(eta), 'LineWidth', 1.5)
colormap spring
xlabel("β_0")
ylabel("β_{sat}")
c = colorbar;
c.Label.String = "log_{10}(η)";
c.LineWidth = 1.5;
ylim([0 1])

%% Spin evolution + chemistry via numerical simulation for the spin precession
%  and stochastic simulation of the reaction chemistry.

function [pS, beta] = spin_chem_sim(nspins, nsteps, dt, eta, Bext)
    % nspins: number of spins to simulate in each of the triplet sublevels
    % nsteps: number of simulation timesteps 
    % dt    : timestep duration (nanoseconds) 
    % eta   : ratio of SCRP recombination to separation rate
    % Bext  : external magnetic field (Tesla)

    % Define the Pauli spin matrices
    paulix = [0 1; 1 0];
    pauliy = [0 -1i; 1i 0];
    pauliz = [1 0; 0 -1];
    
    % Define basis vectors terms of the up & down states for spins 1 & 2.
    uu = [1 0 0 0]'; 
    ud = [0 1 0 0]';
    du = [0 0 1 0]';
    dd = [0 0 0 1]';
    
    % Define the singlet and triplet states in terms of the basis vectors
    s = 1/sqrt(2)*(ud - du); % Singlet state 
    t = 1/sqrt(2)*(ud + du); % T0 state
    tp = uu; % T+ state
    tm = dd;  % T- state

    % Define initial spin states with 1/3 in each triplet sublevel
    psi0 = cat(3, repmat(t, 1, 1, nspins), ...
        repmat(tm, 1, 1, nspins), ...
        repmat(tp, 1, 1, nspins));
    nspins_tot = 3*nspins; 
    
    % Define the two-spin operators
    id = eye(2);
    s1x = kron(paulix, id); 
    s1y = kron(pauliy, id);
    s1z = kron(pauliz, id);
    
    s2x = kron(id, paulix); 
    s2y = kron(id, pauliy);
    s2z = kron(id, pauliz);

    s1 = cat(3, s1x, s1y, s1z);
    s2 = cat(3, s2x, s2y, s2z);
    
    % Define the hyperfine field strengths for each spin
    sigma1 = 0.002; % Tesla
    sigma2 = 0.002; % Tesla

    % Define simulation parameters
    gamma = 28; % electron gyromagnetic ratio (GHz/Tesla)    
    tau = (1:nsteps)*dt; % 

    % Define the random static hyperfine fields for each spin
    b1 = sigma1*randn(3, nspins_tot);
    b2 = sigma2*randn(3, nspins_tot);
    
    % Add z-component from an external magnetic field
    b1(3,:) = b1(3,:) + Bext;
    b2(3,:) = b2(3,:) + Bext;
    
    % Define the Hamiltonian
    H = -2*pi*gamma/2*(tensorprod(s1, b1, 3, 1) + tensorprod(s2, b2, 3, 1));

    % Define rates of recombination and separation
    krec = 0.1; % rate of recombination, per nanosecond
    ksep = krec/eta; % rate of separation, per nanosecond
    
    % Initialize matrices
    psi = zeros(4, nsteps, nspins_tot); % Wavefunction at each timestep
    ampS = zeros(nspins_tot, nsteps); % Singlet amplitude at each timestep
    pS = zeros(nspins_tot, nsteps); % Singlet probability at each timestep
    
    % Intialize counts for chemical dynamics
    rcount = 0; % Number of SCRPs which recombine
    scount = 0; % Number of SCRPs which separate
    
    
    for j = 1:nspins_tot
        % Compute the wavefunction via the unitary matrix
        psi(:,:,j) = expmv(-1i*H(:,:,j), psi0(:,j), tau);

        % Compute the singlet amplitude and probability at each timestep
        ampS(j,:) = s'*psi(:,:,j);
        pS(j,:) = ampS(j,:).*conj(ampS(j,:));
         
        % Define the probability of recombination for each timestep,
        % which is proportional to the singlet fraction
        p_rec = pS(j,:)*dt*krec;

        % Define the probability of separation for each timestep
        p_sep = dt*ksep; 

        % Find the first timestep which a generated random number
        % from zero to one is less than the probability of 
        % each event 

        rec_time = find(rand(1, nsteps)<p_rec, 1);
        sep_time = find(rand(1, nsteps)<p_sep, 1);
        
        if sep_time < rec_time
            scount = scount +1;
        elseif rec_time < sep_time
            rcount = rcount + 1;
        % handle case where the separation and recombination
        % events happen in the same timestep
        elseif sep_time == rec_time
            f = randi([0 1]);
            rcount = rcount + f;
            scount = scount + (1-f);
        % handle cases where each event does not happen
        % during the simulation time
        elseif isempty(rec_time) && ~isempty(sep_time)
            scount = scount +1;
        elseif  isempty(sep_time) && ~isempty(rec_time)
            rcount = rcount + 1;            
        end
    end

    % compute beta, the SCRP separation fraction 
    beta= scount/(scount+rcount);

    % display info. for how many of the simulated molecules react
    disp("eta = " + string(eta) + ", Bext = " + string(Bext*1000) +  " mT");
    disp("Total SCRPs reacted = " + string(100*(rcount + scount)/nspins_tot) + "%");
end

%% Spin evolution + chemistry via analytical expressions for the spin precession
%  and stochastic simulation of the reaction chemistry. 

function [pS, beta] = spin_chem_sim_analytical(nspins, nsteps, dt, eta, Bext, initial_state)
    % nspins: number of spins per trial
    % bz: applied magnetic field, T 
    % eta: ratio of recombination rate to separation rate 
    % nsteps: total number of timesteps
    % dt: timestep (ns) 
    % initial_state: "T0" or "Tpm", SCRP initial spin state. 

    gamma = 28; % electron gyromagnetic ratio, GHz/T
    nEta = length(eta);
    nB = length(Bext);
    tau = (1:nsteps)*dt;
    pS = zeros(nsteps,nspins);
    beta = zeros(nEta, nB);
    sigma1 = 0.002; % hyperfine field strength, spin 1 (T)
    sigma2 = 0.002; % hyperfine field strength, spin 2  (T)
    for m = 1:nEta % loop through etas
        for k = 1:nB % loop through Bs
            scount = 0; % initialize count of separating SCRPs
            rcount = 0; % initialize count of recombining SCRPs
            tic
            b1 = sigma1*randn(3, nspins);
            b2 = sigma2*randn(3, nspins);
    
            b1(3,:) = b1(3,:) + Bext(k);
            b2(3,:) = b2(3,:) + Bext(k);
            B1 = squeeze(sum(b1.^2, 1)).^.5;
            B2 = squeeze(sum(b2.^2, 1)).^.5;
        
            for j = 1:nspins
                B1t = 2*pi*gamma*B1(j)*tau/2; 
                B2t = 2*pi*gamma*B2(j)*tau/2;
                if initial_state == "T0"
                    pS(:,j) = B1(j).^(-2).*B2(j).^(-2).*(b1(3,j).*B2(j).*cos(B2t).*sin(B1t)...
                    +(-1).*(B1(j).*b2(3,j).*cos(B1t)+(b1(2,j).*b2(1,j)+(-1) ...
                    .*b1(1,j).*b2(2,j)).*sin(B1t)).*sin(B2t)).^2;
                elseif initial_state == "Tpm"
                     pS(:,j) = (1/4).*B1(j).^(-2).*B2(j).^(-2).*(2.*(b1(1,j).^2+b1(2,j).^2).*B2(j)...
                    .^2.*cos(B2t).^2.*sin(B1t).^2+B1(j).*sin(2.*B1t) ...
                      .*((b1(2,j).*b2(1,j)+(-1).*b1(1,j).*b2(2,j)).*b2(3,j).*((-1)+cos( ...
                      2.*B2t))+(-1).*(b1(1,j).*b2(1,j)+b1(2,j).*b2(2,j)).*B2(j).* ...
                      sin(2.*B2t))+2.*(B1(j).^2.*(b2(1,j).^2+b2(2,j).^2).*cos(B1t ...
                      ).^2.*sin(B2t).^2+sin(B1t).^2.*((b1(3,j).^2.*(b2(1,j) ...
                      .^2+b2(2,j).^2)+(-2).*b1(3,j).*(b1(1,j).*b2(1,j)+b1(2,j).*b2(2,j)) ...
                      .*b2(3,j)+(b1(1,j).^2+b1(2,j).^2).*b2(3,j).^2).*sin(B2t).^2+ ...
                      b1(3,j).*(b1(2,j).*b2(1,j)+(-1).*b1(1,j).*b2(2,j)).*B2(j).*sin(2.* ...
                      B2t))));
                else
                    error("Error: options are T0 and Tpm")
                end
                krec = 0.1; % per nanosecond
                ksep = krec/eta(m); % per nanosecond
        
                for i = 1:size(pS,1)
                    p_rec = pS(i,j)*dt*krec;
                    p_sep = dt*ksep;
                    rec = rand<p_rec;
                    sep = rand<p_sep;
                    if rec+sep == 2
                        f = randi([0 1]);
                        rcount = rcount + f;
                        scount = scount + (1-f);
                        % disp(i)
                        break 
                    elseif rec == 1
                           rcount = rcount + 1;
                           break
                    elseif sep == 1
                            scount = scount +1 ;
                            break
                    end
                end
            end
                beta(m,k) = scount/(scount+rcount);
                disp(strcat("eta = ", string(eta(m)), ", Bext = ", string(Bext(k))));
                disp(strcat("Total SCRPs reacted = ", string(100*(rcount + scount)/nspins), "%"));
                toc
        end
    end
end    
