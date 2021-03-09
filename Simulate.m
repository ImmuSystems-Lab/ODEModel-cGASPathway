function Simulate
%Simulate the cGAS/JAKSTAT pathway using an Ordinary Differential Equation
%Model. Provides two options that can be commented out. Option 1 reproduces
%one simulation result from figure 2 in the paper:
%https://doi.org/10.1016/j.jtbi.2018.11.001 and Option 2 plots just the
%dynamic species in the model

% Author: Robert Gregg
% Email: rwg16@pitt.edu
% Website: http://shoemakerlab.pitt.edu/
% Edited: 14-April-2019


%Define some constants for the model
M2nM = 1e9; %converts from M to nM
V = 3e-12; % volume of the cell
Na = 6.02e23; %Avogadro number


%Initialize the 13 states
x0 = zeros(13,1);
%Some states have a non-zero starting value
    %Nonzero states: 'cGAS' 'DNA' 'Sting' 'IRF3'
x0([1:3,5]) = M2nM.*([1e3 1e3 1e3 1e4])./(V*Na);

%Example Parameter set for the ODE model
theta  =[2.68987732928698;4.85048330610477;0.0356149813628553;7.48699230097061;517.405644937774;22328.3852151619;11226.3681590503;0.934069560710524;206.944580408962;10305.4609911003;47639702.9463696;3.84744551149676;13.0064080978012;78.2048219250749;0.0208760041492399;0.00587644557176817;0.00100261806969412;0.0112153351557179;0.00100021202033750;99.9465979604163;15.1436478237262;0.0276031221466779;237539.324901876;61688.2589525162;0.959997274392481;0.346976450528852;12242.8735596321;1.23993808595094;1.51008015416415;0.346976450528852;0.165006100132023;6.92946969842556;0.0177991802039533];

%Mass Balances (total = active + inactive)
cGAStot = (M2nM.*1e3)/(V*Na);
Stingtot = (M2nM.*1e3)/(V*Na);
IRF3tot = (M2nM.*1e4)/(V*Na);

%Time frame (in hours)
t0 = [0,48];

%Simulate the model, pass extra parameters through (theta+mass balances) 
[time,species]=ode15s(@Model,t0,x0,[],theta,cGAStot,Stingtot,IRF3tot);

%%% Option #1: Plot mass balanced species (Figure 2 in paper)
figure()
StatesForData = [species(:,1:3),...
               Stingtot-species(:,3),...
               species(:,4:5),...
               IRF3tot-species(:,5),...
               species(:,6:end)];
AddedStates = {'cGAS' 'DNA' 'Sting' 'Stingc' 'cGAMP' 'IRF3' 'IRF3c' 'IFNbm'... 
          'IFNb' 'STATP2n' 'SOCSm' 'IRF7m' 'TREX1m' 'IRF7Pn' 'TREX1'};
%Create a plot of the simulation
for i=1:length(AddedStates)
   hold on
   subplot(4,4,i)
   plot(time,StatesForData(:,i),'LineWidth',1.5)
   title(AddedStates{i})
   xlabel('Time (hr)')
   xlim(t0)
   xticks(0:12:t0(end))
   ylabel('nM')
end

%%% Option #2: Plot just the species
%Names for the simulated states
figure()
States = {'cGAS' 'DNA' 'Sting' 'cGAMP' 'IRF3' 'IFNbm' 'IFNb' 'STATP2n'...
          'SOCSm' 'IRF7m' 'TREX1m' 'IRF7Pn' 'TREX1'};
%Create a plot of the simulation
for i=1:length(States)
   hold on
   subplot(4,4,i)
   plot(time,species(:,i),'LineWidth',1.5)
   title(States{i})
   xlabel('Time (hr)')
   xlim(t0)
   xticks(0:12:t0(end))
   ylabel('nM')
end

    function dx = Model(t,x,theta,cGAStot,Stingtot,IRF3tot)
        %Unpack all of the States
        cGAS = x(1);
        DNA = x(2);
        Sting = x(3);
        cGAMP = x(4);
        IRF3 = x(5);
        IFNbm = x(6);
        IFNb = x(7);
        STATP2n = x(8);
        SOCSm = x(9);
        IRF7m = x(10);
        TREX1m = x(11);
        IRF7Pn = x(12);
        TREX1 = x(13);

        %Define all the Parameters
        k1f = theta(1);
        k1r = theta(2);
        k3f = theta(3);
        k3r = theta(4);
        k4f = theta(5);
        kcat5 = theta(6);
        Km5 = theta(7);
        k5r = theta(8);
        kcat6 = theta(9);
        Km6 = theta(10);
        kcat7 = theta(11);
        Km7 = theta(12);
        kcat8 = theta(13);
        Km8 = theta(14);
        k8f = theta(15);
        k9f = theta(16);
        k10f1 = theta(17);
        k10f2 = theta(18);
        k11f = theta(19);
        k12f = theta(20);
        k13f = theta(21);
        k6f = theta(22);
        kcat2 = theta(23);
        Km2 = theta(24);
        tau4 = theta(25);
        tau6 = theta(26);
        tau7 = theta(27);
        tau8 = theta(28);
        tau9 = theta(29);
        tau10 = theta(30);
        tau11 = theta(31);
        tau12 = theta(32);
        tau13 = theta(33);
        
        %ODES
        dx = zeros(13,1);
        % cGAS
        dx(1)=-k1f*cGAS*DNA + k1r*(cGAStot - cGAS);
        % vDNA
        dx(2)=-k1f*cGAS*DNA + k1r*(cGAStot - cGAS) - kcat2*TREX1*DNA / (Km2 + DNA);
        % Sting
        dx(3)=-k3f*cGAMP*Sting + k3r*(Stingtot - Sting);
        % cGAMP     
        dx(4)=k4f*(cGAStot - cGAS) - k3f*cGAMP*Sting + k3f*(Stingtot - Sting) - tau4*cGAMP;
        % IRF3         
        dx(5)=-kcat5*IRF3*(Stingtot - Sting) / (Km5 +IRF3) + k5r*(IRF3tot - IRF3);  
        % IFNbm         
        dx(6)=kcat6*(IRF3tot - IRF3) / (Km6 + (IRF3tot - IRF3)) + k6f*IRF7Pn - tau6*IFNbm;
        % IFNb     
        dx(7)=kcat7*IFNbm / (Km7 + IFNbm) - tau7*IFNb;
        % STATP2n
        dx(8)=kcat8*IFNb / (Km8 + IFNb) * 1/(1+k8f*SOCSm) - tau8*STATP2n;
        % SOCSm
        dx(9)=k9f*STATP2n - tau9*SOCSm;
        % IRF7m         
        dx(10)=k10f1*STATP2n + k10f2*IRF7Pn - tau10*IRF7m;
        % TREX1m      
        dx(11)=k11f*STATP2n - tau11*TREX1m;
        % IRF7Pn          
        dx(12)=k12f*IRF7m - tau12*IRF7Pn;
        % TREX1         
        dx(13)=k13f*TREX1m - tau13*TREX1;
    end

end

