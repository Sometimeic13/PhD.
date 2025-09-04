% for B from 2 Helmholz coils together - W17/75
% you have to choose if in series or in paralell
% Due to heatLim we must calculate with 2 coils together. N is for single coil.
% rules for addition of inductance: https://www.electronics-tutorials.ws/inductor/series-inductors.html, we incorectly ignore additional term
% Expected common B when 2 coils are at distance R is sqrt(2) * B
% for 2 cases: when PheatLimit is and is not relevant

% compare with https://www.accelinstruments.com/Magnetic/Magnetic-field-calculator.html

% TODO:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clf;

rho = 1 / 6e7;
mu = 4*pi*1e-7;

amplifierStronger = 1; % we had 2 variants: 180 and 1000 W, see bellow IampMax and UampMax
frequencyHigher = 0;   % for 2 variants, , see bellow f
applyHeatLimit = 1;    % 1: practical restriction to keep functional temperature regulation, 0: if thermal box can be cooled somehow, e.g. by its compressor

coilDiameter = 0.12;   % inner coil diameter, [m] plus material thickness and winding thickness...10mm
IampMax = [4.5,10];    % Imax of several amlifiers [A], set at r.50
UampMax = [40,100];    % Umax of several amlifiers [V], set at r.50
f = [1000,10000];      % working frequency, [Hz], set at r.48
N0 = 10;               % coil turns in each layer of each of 2 coils
mt = 0.003;            % material thickness, [m]
Ramp = 1.7;            % aditional value of R [Ohm] provided by potentiometer [Ohm], in order to provide amlifier with proper input        
PheatLimit = 2;        % applied if applyHeatLimit = 1; P_joule inside thermal box limited due to box overheating, [W]
connection = 1;        % 1: in series, 0: paralell 

ftsz = 16;

%%%%%%%%%%%%%%%%%%%

phi = 0.0005:0.0001:0.0022; % wire diameter, [m]
N = N0:N0:400;              % coil turns, for each of 2 coils

if amplifierStronger
  UampMax2 = UampMax(2);
else
  UampMax2 = UampMax(1);
end

if amplifierStronger
  IampMax2 = IampMax(2);
else
  IampMax2 = IampMax(1);
end

if frequencyHigher
  fa = f(2);
else
  fa = f(1);
end

r1 = coilDiameter / 2 + mt + phi / 2;    % coil inner radius, diameter for open space 11.5cm + material thickness, [m]

%size(phi)
%size(N)

%%%%%%%%%%%%%%%%%%% 

for i = 1 : length(phi) % wire diameter
    
  gap(i) = 0.0007 / (0.0022 / phi(i));   % gap between wires, [m], more thick higher gap, 0.0007 for phi=0.0022 empiricaly estimated W15/192
  h(i) = N0 * (phi(i) + gap(i));         % coil height  
  
  for j = 1 : length(N) % coil turns
    r2(i,j) = r1(i) + N(j) / (N0 - 1) * phi(i); % each coil outer radius, -1 due to ...
    
    R(i,j) = 0;
    for k = 1 : N(j) / N0 % actual r increases
      R(i,j) = R(i,j) + rho * 2 * pi * (r1(i) + (k - 1) * phi(i)) * N0 / (pi *(phi(i)/2)^2); 
    end
    
    if connection == 1;        % 1: in series, 0: paralell
      R2(i,j) = 2 * R(i,j); % for 2 H. coils in series
    else
      R2(i,j) = R(i,j) / 2; % for 2 H. coils in paralell  
    end
    
    L(i,j) = 31.6 * 1e-6 * N(j)^2 * r1(i)^2 / (6 * r1(i) + 9 * h(i) + 10 * (r2(i,j) - r1(i))); % www.circuits.dk/calculator_multi_layer_aircore.htm
    XL(i,j) = 2 * pi * fa * L(i,j); % for 1 coil
    
    if connection == 1;        % 1: in series, 0: paralell
      XL2(i,j) = 2 * XL(i,j); % for 2 H. coils in series
    else
      XL2(i,j) = XL(i,j) / 2; % for 2 H. coils in paralell
    end
    
%     if R2(i,j) <= Ramp % Ohm nonsense - no disconnection
%       Rprot(i,j) = Ramp - R2(i,j);
%     else
%       Rprot(i,j) = 0;
%     end
    %Rprot(i,j) = 0; % test its role
    %Z(i,j) = sqrt((R2(i,j) + Rprot(i,j))^2 + XL2(i,j)^2);
  
    Z(i,j) = sqrt((R2(i,j) + Ramp)^2 + XL2(i,j)^2);
  
%% W17/82: min not R but Z; Problem - not simple adding, but under the sqrt
%     Z(i,j) = sqrt(R2(i,j)^2 + XL2(i,j)^2);
% 
%     if Z(i,j) <= Ramp % Ohm
%       Rprot(i,j) = Ramp - R2(i,j);
%     else
%       Rprot(i,j) = 0;
%     end
%     %Rprot(i,j) = 0; % test its role
%     
%%    
    
    % I(i,j) = Ua / Z(i,j); % wrong, W15/194

    if applyHeatLimit
      Imax(i,j) = min([IampMax2, UampMax2 / Z(i,j), sqrt(PheatLimit / R(i,j))]);
    else
      Imax(i,j) = min(IampMax2, UampMax2 / Z(i,j));
    end
    
    %Imax(i,j) = 1; % for evaluation what B we get for fixed 1 Amp; no: better use applyHeatLimit=0 and IampMax2=1
    
    %r3(i,j) = (r2(i,j) + r1(i)) / 2; % average
    %B(i,j) = 8 / (5 * sqrt(5)) * mu * Imax(i,j) * N(j) / r3(i,j); % https://en.wikipedia.org/wiki/Helmholtz_coil, B is independent on L and Z

    Bpart(i,j) = 0;
    for k = 1 : N(j) / N0 % more acurate according to actual r increases, this correction makes max 0.5 mt (material thickness) for high N and phi; in region of interest only 0.15
      Bpart(i,j) = Bpart(i,j) + 1 / (r1(i) + (k - 1) * phi(i)); 
    end
    B(i,j) = 8 / (5 * sqrt(5)) * mu * Imax(i,j) * N0 * Bpart(i,j);
    B(i,j) = sqrt(2) * B(i,j); % for 2 Helmholz coils with mutual distance R 
    
    U(i,j) = Imax(i,j) * Z(i,j);
    P(i,j) = U(i,j) * Imax(i,j);

    nLayers(i,j) = N(j) / N0;
    heightLayers(i,j) = nLayers(i,j) * phi(i);

%     if Imax(i,j) <= IampMax2 && U(i,j) <= UampMax2 % for applyHeatLimit = 1 just 2nd condition would be enough
%       Busable(i,j) = B(i,j);
%       Iusable(i,j) = Imax(i,j);
%       Uusable(i,j) = U(i,j);
%       Pusable(i,j) = U(i,j);
%       nLayersUsable(i,j) = nLayers(i,j);
%       heightLayersUsable(i,j) = heightLayers(i,j);
%     else
%       Busable(i,j) = NaN;
%       Iusable(i,j) = NaN;
%       Uusable(i,j) = NaN;
%       Pusable(i,j) = NaN;
%       nLayersUsable(i,j) = NaN;
%       heightLayersUsable(i,j) = NaN;
%     end
    
  end 
end
%size(r2)
%size(Z)

[phim,Nm] = meshgrid(phi,N);
%size(phim)
%size(Nm)

surf(1000*phim',Nm',R2)
xlabel('phi [mm]')
ylabel('N')
zlabel('R [Ohm]')
axis tight

figure
surf(1000*phim',Nm',1000 * L)
xlabel('phi [mm]')
ylabel('N')
zlabel('L [mH]')
axis tight

figure
surf(1000*phim',Nm',XL2)
xlabel('phi [mm]')
ylabel('N')
zlabel('XL [Ohm]')
axis tight

figure
surf(1000*phim',Nm',Z)
xlabel('phi [mm]')
ylabel('N')
zlabel('Z [Ohm]')
axis tight

figure
surf(1000*phim',Nm',Imax)
xlabel('phi [mm]')
ylabel('N')
zlabel('Imax [A]')
axis tight

% figure
% surf(1000*phim',Nm',1000*r3)
% xlabel('phi [mm]')
% ylabel('N')
% zlabel('r3 [mm]')
% axis tight

%colormap hot

figure
surf(1000*phim',Nm',1000*B)
xlabel('\phi [mm]')
ylabel('no. of turns')
%zlabel('Bmax [mT]')
zlabel('B_{max} [mT]')
axis tight
colorbar
set(gca,'FontSize',ftsz)

return

figure
surf(1000*phim',Nm',U)
xlabel('phi [mm]')
ylabel('N')
zlabel('voltage [V]')
axis tight

figure
surf(1000*phim',Nm',P)
xlabel('phi [mm]')
ylabel('N')
zlabel('power [W]')
axis tight

% figure
% surf(1000*phim',Nm',1000*Busable)
% xlabel('phi [mm]')
% ylabel('N')
% zlabel('B usable [mT]')
% axis tight
% 
% figure
% surf(1000*phim',Nm',Iusable)
% xlabel('phi [mm]')
% ylabel('N')
% zlabel('I usable [A]')
% axis tight
% 
% figure
% surf(1000*phim',Nm',Uusable)
% xlabel('phi [mm]')
% ylabel('N')
% zlabel('U usable [V]')
% axis tight
% 
% figure
% surf(1000*phim',Nm',Pusable)
% xlabel('phi [mm]')
% ylabel('N')
% zlabel('P usable [W]')
% axis tight

figure
surf(1000*phim',Nm',nLayers)
xlabel('phi [mm]')
ylabel('N')
zlabel('number of layers')
axis tight

figure
surf(1000*phim',Nm',1000*heightLayers)
xlabel('phi [mm]')
ylabel('N')
zlabel('height of layers [mm]')
axis tight

Bmax = max(max(B));
%find(B == Bmax)
[row,col] = ind2sub(size(B), find(B == Bmax));

fprintf('Ramp[Ohm] phi[mm] N Z[Ohm] Imax[A] B[mT] \n')
[Ramp, 1000 * phi(row), N(col), Z(row,col), Imax(row,col), 1000 * B(row,col)]
fprintf('freq \n')
fa
%max(1000*Busable)
%max(1000*Busable')
