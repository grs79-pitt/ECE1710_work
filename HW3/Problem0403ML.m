%Distribution System Modelling and Analysis, Example 4.2
%Written by William Kersting and Robert Kerestes
clear all
clc

%Modified for Problem 4.3
%Loading the phase impedance matrix from Problem 4.2
load('problem0402_zabc.mat') %For Problem 4.3

j = sqrt(-1);

length = 6000; %%Length from
Zabc = zabc*length/5280; %Getting actual impedance
Zl_perphase=25+j*15; %Load Impedance
ZLmatrix=[Zl_perphase,0,0;0,Zl_perphase,0;0,0,Zl_perphase]; %From problem 4.3
Vll=12470; %Source Voltage
Ztot=Zabc+ZLmatrix;
Ztotinv=inv(Ztot); %Inverting matrix
Es = [Vll/sqrt(3);Vll/sqrt(3)*exp(-j*2*pi/3);Vll/sqrt(3)*exp(j*2*pi/3)]; %Source Voltage array

Iabc=Ztotinv*Es;

VLabc = (ZLmatrix*Ztotinv)*Es;

    
[VLabc_mag, VLabc_phase] = rec2pol(VLabc); %Convert Load Voltages in Polar form
[Iabc_mag, Iabc_phase] = rec2pol(Iabc); %Convert Line Currents in Polar form
Vdrop=(abs(Es-VLabc)/Es(1))*100; %Calculate the voltage drop

%Calculating powers
loadpower=VLabc.*conj(Iabc);
realloadpower=real(loadpower);
reactiveloadpower=imag(loadpower);


%Display all data from the script 
 disp(['The Load Voltage at Phase A is ', num2str(VLabc_mag(1)), ' ∠', num2str(VLabc_phase(1)), ' V']);
 disp(['The Load Voltage at Phase B is ', num2str(VLabc_mag(2)), ' ∠', num2str(VLabc_phase(2)), ' V']);
 disp(['The Load Voltage at Phase C is ', num2str(VLabc_mag(3)), ' ∠', num2str(VLabc_phase(3)), ' V', num2str(newline)]);
 
 disp(['The Current at Phase A is ', num2str(Iabc_mag(1)), ' ∠', num2str(Iabc_phase(1)), ' A']);
 disp(['The Current at Phase B is ', num2str(Iabc_mag(2)), ' ∠', num2str(Iabc_phase(2)), ' A']);
 disp(['The Current at Phase C is ', num2str(Iabc_mag(3)), ' ∠', num2str(Iabc_phase(3)), ' A', num2str(newline)]);
 
disp(['The real power at Phase A is, ', num2str(realloadpower(1)/1000), 'kW']);
disp(['The real power at Phase B is, ', num2str(realloadpower(2)/1000), 'kW']);
disp(['The real power at Phase C is, ', num2str(realloadpower(3)/1000), 'kW', num2str(newline)]);

disp(['The reactive power at Phase A is, ', num2str(reactiveloadpower(1)/1000), 'kvar']);
disp(['The reactive power at Phase B is, ', num2str(reactiveloadpower(2)/1000), 'kvar']);
disp(['The reactive power at Phase C is, ', num2str(reactiveloadpower(3)/1000), 'kvar', num2str(newline)]);

disp(['The Voltage Drop at Phase A is, ', num2str(Vdrop(1)), '%']);
disp(['The Voltage Drop at Phase B is, ', num2str(Vdrop(2)), '%']);
disp(['The Voltage Drop at Phase C is, ', num2str(Vdrop(3)), '%']);




