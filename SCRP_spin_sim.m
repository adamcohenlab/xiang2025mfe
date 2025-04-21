%% Simulate a single spin precessing in an arbitrary magnetic field
% Define the Pauli spin matrices
paulix = [0 1; 1 0];
pauliy = [0 -1i; 1i 0];
pauliz = [1 0; 0 -1];

% Define the up and down spin states. 
u = [1 0]';
d = [0 1]';

% Initialize spin in the up state: 
psi0 = u; 
theta0 = 2*acos(abs(psi0(1,:)));
phi0 = angle(psi0(2,:))-angle(psi0(1,:));

% % Or: pick any angle and phase: 
% theta0 = 0;
% phi0 = 0;
% psi0 = cos(theta0/2)*u + sin(theta0/2)*exp(1i*phi0)*d;

% Provide simulation parameters
sigma1 = 0.002; % hyperfine field strength (Tesla)
gamma = 28; % electron gyromagnetic ratio (GHz/Tesla)
nsteps = 50; % number of timesteps 
dt = 0.15; % simulation timestep (ns)
tau = (1:nsteps)*dt; 

% Initialize magnetic field by drawing each component from a normal dist.
b1 = sigma1*randn(3, 1); 

% Define the Hamiltonian
H = -2*pi*gamma/2*(b1(1)*paulix + b1(2)*pauliy + b1(3)*pauliz);

% Compute the time evolution via the unitary matrix
psi  = expmv(-1i*H, psi0, tau);

% Compute amplitudes and probabilities of being in the up vs. down state
ampU = u'*psi;
ampD = d'*psi;

pU = ampU.*conj(ampU);
pD = ampD.*conj(ampD);

% Plot the state probability vs. time
figure(1)
subplot(1,2,1)
hold on
plot(tau,pU) 
plot(tau, pD)
ylabel("State probability")
xlabel("Time (ns)")
legend("Up", "Down")
ylim([0 1])

% --- Illustrate the spin precession vs. time on the Bloch Sphere --- %

subplot(1,2,2)
hold on

[X, Y, Z] = sphere(40);
a = surf(X, Y, Z, 'FaceColor',[.8 .8 .8],'EdgeColor','none');
alpha(.4);
view(3); daspect([1 1 1]); axis off
lightangle(45,30); lighting phong; 
set(a,'AmbientStrength',.5)
set(a,'SpecularColorReflectance',.8,'SpecularExponent',50);

unitb = 1.5*b1/norm(b1); % Plot the direction of the magnetic field
quiver3(0, 0, 0, unitb(1), unitb(2), unitb(3), 0, 'LineWidth',2) 

% Plot the starting state on the Bloch Sphere
f = 1.05;
x = sin(theta0).*cos(phi0)*f;
y = sin(theta0).*sin(phi0)*f;
z = cos(theta0)*f;
plot3(x,y,z,'bo', 'LineWidth', 3, 'HandleVisibility', 'off')

% Plot the spin state at each timestep on the Bloch Sphere
theta = 2*acos(abs(psi(1,:))); % Find the angle from |0>
phi = angle(psi(2,:))-angle(psi(1,:)); % Find the global phase

x = sin(theta).*cos(phi)*f; % Spherical to cartesian coordinates
y = sin(theta).*sin(phi)*f;
z = cos(theta)*f;

scatter3(x,y,z, 30, tau, "o", 'LineWidth', 1.5);
c = colorbar;
c.LineWidth = 1;
c.Title.String = "Time (ns)";
clim([min(tau), max(tau)]);

%% Simulate a two-spin system in separate arbitrary magnetic fields. 

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

% Initialize the wavefunction in one of the spin states; here, T-
psi0 = tm; 

% Define the two-spin sx, sy, and sz operators 
id = eye(2);
s1x = kron(paulix, id); 
s1y = kron(pauliy, id);
s1z = kron(pauliz, id);

s2x = kron(id, paulix); 
s2y = kron(id, pauliy);
s2z = kron(id, pauliz);

% Provide parameters for the simulation 
sigma1 = 0.002; % spin 1 hyperfine field strength (Tesla)
sigma2 = 0.002; % spin 2 hyperfine field strength (Tesla)
gamma = 28; % electron gyromagnetic ratio (GHz/Tesla)
nsteps = 500; % number of timesteps
dt = 0.05; % timestep (ns)
tau = (1:nsteps)*dt; 

% Initialize magnetic field by drawing each component from a normal dist.
b1 = sigma1*randn(3, 1); % 
b2 = sigma1*randn(3, 1); % 

% % Or define the magnetic fields manually
% b1 = [.01 0 0];
% b2 = [0 0 .01];

% Define the Hamiltonian
H = -2*pi*gamma/2*(b1(1)*s1x + b1(2)*s1y + b1(3)*s1z + ...
    b2(1)*s2x + b2(2)*s2y  + b2(3)*s2z);

% Compute the time evolution via the unitary matrix
psi = expmv(-1i*H, psi0, tau);

% Compute the amplitudes of each of the S, T-, T0, and T+ spin states
ampS = s'*psi;
ampT = t'*psi;
ampTp = tp'*psi;
ampTm = tm'*psi;

% Compute the probabilities of each of the spin states
pS = ampS.*conj(ampS);
pT = ampT.*conj(ampT);
pTp = ampTp.*conj(ampTp);
pTm = ampTm.*conj(ampTm);

% Plot the probabilities vs. time
figure(2)
hold on
plot(tau, pS)
plot(tau, pT)
plot(tau, pTp)
plot(tau, pTm)
% plot(tau, pS+pT+pTp+pTm) % Plot the sum to check normalization 

ylabel("State probability")
xlabel("Time (ns)")
legend("S", "T_0", "T_-", "T_+")