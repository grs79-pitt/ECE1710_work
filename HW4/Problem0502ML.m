%Distribution System Modelling and Analysis, Example 5.2
%Written by William Kersting and Robert Kerestes
clear all
clc

%MODIFIED FOR PROBLEM 5.2
j = sqrt(-1);
f = 60;
omega = 2*pi*f;

%Defining the phase conductor data for the first line
phase1.GMR = 0.00446; %Changed to 1/0 ACSR for Problem 5.2
phase1.resistance = 1.12;  %Changed to 1/0 ACSR for Problem 5.2
phase1.diameter = 0.398;  %Changed to 1/0 ACSR for Problem 5.2
phase1.ncond = 3;

%Defining the neutral conductor data
neutral.GMR = 0.00446;  %Changed to 1/0 ACSR for Problem 5.2
neutral.resistance = 1.12;  %Changed to 1/0 ACSR for Problem 5.2
neutral.diameter = 0.398;  %Changed to 1/0 ACSR for Problem 5.2
neutral.ncond = 1;

ncond = phase1.ncond+neutral.ncond;

%Initializing matrix sizes
r = zeros(ncond,1);
Dshunt = zeros(ncond,ncond);
S  = zeros(ncond,ncond);
zprim = zeros(ncond,ncond);

%Defining the distance vector d

d1 = 0+j*29; d2 = 7 +j*29; d3 = 2.5+j*29;  %Changed for Problem 5.2

d4 = 4+j*25; %Changed for Problem 5.2

d = [d1;d2;d3;d4]; %Changed for Problem 5.2

%Defining the resistance vector r
for i = 1:1:ncond
    
    if i <= phase1.ncond 
        
        r(i) = phase1.resistance;
        
    else 
        
        r(i) = neutral.resistance;
        
    end 
        
end

 
%Calculating the inpedance matrix Dshunt
for i = 1:1:ncond
    
    for k = 1:1:ncond
        
        if i == k && i <= phase1.ncond 
            
            Dshunt(i,k) = phase1.diameter/24;
            
        elseif i == k && i > phase1.ncond
            
            Dshunt(i,k) = neutral.diameter/24;
            
        else
            
            Dshunt(i,k) = abs(d(i) - d(k));
            
        end
        
    end
    
end

%Calculating the image distance matrix S
for i = 1:1:ncond
    
    for k = 1:1:ncond
        
        S(i,k) = abs(d(i)-conj(d(k)));
        
    end
    
end

%Calculating the primitive potential coefficient matrix
for i = 1:1:ncond

    for k = 1:1:ncond
            
            Pprim(i,k) = 11.17689*log(S(i,k)/Dshunt(i,k));  
    end
    
end


% %Partitioning Pprim
for i = 1:1:phase1.ncond

    for k = 1:1:phase1.ncond
        
        Pij(i,k) = Pprim(i,k);
    
    end
    
end

for i = 1:1:phase1.ncond
    
   for k = 1:1:neutral.ncond
    
    Pin(i,k) = Pprim(i,phase1.ncond+k);
    
    end
    
end

for i = 1:1:neutral.ncond
    
   for k = 1:1:phase1.ncond
    
    Pnj(i,k) = Pprim(phase1.ncond+i,k);
    
    end
    
end

for i = 1:1:neutral.ncond
    
   for k = 1:1:neutral.ncond
       
    Pnn(i,k) = Pprim(phase1.ncond+i,phase1.ncond+k);

   end
   
end

%Performing the Kron reduction
Pabc = Pij-Pin*Pnn^-1*Pnj;

%Calculating phase capacitance matrix
Cabc = Pabc^-1;

%Calculating shunt admittance matrix
yabc = j*omega*Cabc;

% %Calculating the partitioned phase impedance matrices
% y11abc = yabc(1:3,1:3);
% y12abc = yabc(1:3,4:6);
% y21abc = yabc(4:6,1:3);
% y22abc = yabc(4:6,4:6);

fprintf('The distance matrix is ')
disp('D = ') 
fprintf('\n')
disp(Dshunt)

fprintf('The phase admittance matrix is: \n\n')
disp('yabc = ')
disp(yabc)

% fprintf('The partitioned phase admittance matrices are \n\n')
% disp('[y11]abc = ') 
% fprintf('\n')
% disp(y11abc)
% 
% disp('[y12]abc = ') 
% fprintf('\n')
% disp(y12abc)
% 
% disp('[y21]abc = ') 
% fprintf('\n')
% disp(y21abc)
% 
% disp('[y22]abc = ') 
% fprintf('\n')
% disp(y22abc)

% disp('The distance matrix in feet for this line is')
% display(D)
% 
% disp('The primitive impedance matrix in ohms/mile for this line is')
% display(zprim)
% 
% disp('In partioned form, the submatrices in ohms/mile are')
% display(zij)
% display(zin)
% display(znj)
% display(znn)
% 
% disp('The "Kron" reduced phase impedance matrix in ohms/mile is ')
% display(zabc)
% 
% 
% 
% disp('The neutral transformation matrix is')
% display(tn)

