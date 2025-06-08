%Distribution System Modelling and Analysis, Example 4.1
%Written by William Kersting and Robert Kerestes

%Modified for Problem 4.5
clear all
clc
j = sqrt(-1);

%Defining the phase conductor data
phase.GMR = 0.0244; %Modified for 336.400 26/7 ACSR
phase.resistance = .306; %Modified for 336.400 26/7 ACSR
phase.ncond = 2; %Two phase

%Defining the neutral conductor data
neutral.GMR = 0.00814; %Modified for 4/0 6/1 ASCR
neutral.resistance = 0.592; %Modified for 4/0 6/1 ASCR
neutral.ncond = 1; %1 neutral

ncond = phase.ncond+neutral.ncond;

%Initializing matrix sizes
r = zeros(ncond,1);
D = zeros(ncond,ncond);
zprim = zeros(ncond,ncond);

%Defining the distance vector d

d1 = 0+j*29; d2 = 7+j*29; d3 = 4+j*25; %Changed for Problem 4.2

d = [d1;d2;d3];

%Defining the resistance vector r
for i = 1:1:ncond
    
    if i <= ncond -1 
        
        r(i) = phase.resistance;
    
    else
        
        r(i) = neutral.resistance;
        
    end 
        
end

%Calculating the inpedance matrix D
for i = 1:1:ncond
    
    for k = 1:1:ncond
        
        if i == k && i <= ncond -1
            
            D(i,k) = phase.GMR;
            
        elseif i == k && i > ncond - 1
            
            D(i,k) = neutral.GMR;
            
        else
            
            D(i,k) = abs(d(i) - d(k));
            
        end
        
    end
    
end

%Calculating the primitive impedance matrix
for i = 1:1:ncond

    for k = 1:1:ncond
       
        if i == k
            
            zprim(i,k) = r(i)+0.0953+j*0.12134*(log(1/D(i,k))+7.93402);
            
        else
            
            zprim(i,k) = 0.0953+j*0.12134*(log(1/D(i,k))+7.93402);
            
        end
          
    end
    
end

%Partitioning Zprim
for i = 1:1:phase.ncond

    for k = 1:1:phase.ncond
        
        zij(i,k) = zprim(i,k);
    
    end
    
end

for i = 1:1:phase.ncond

    for k = 1:1:neutral.ncond
    
        zin(i,k) = zprim(i,k+phase.ncond);

    end
    
end

for i = 1:1:neutral.ncond

    for k = 1:1:phase.ncond
        
        znj(i,k) = zprim(i+phase.ncond,k);
    
    end

end

for i = 1:1:neutral.ncond

    for k = 1:1:neutral.ncond
        
        znn = zprim(i+phase.ncond,k+phase.ncond);

    end

end


%Performing the Kron reduction

%Padding the zeroed elements with zeros depending on the number of
%conductors
zabctemp = zij-zin*znn^-1*znj;
zabc=zeros(3,3);
if (phase.ncond==2)
    zabc(1,1)=zabctemp(1,1);
    zabc(1,3)=zabctemp(1,2);
    zabc(3,1)=zabctemp(2,1);
    zabc(3,3)=zabctemp(2,2);
elseif (phase.ncond==1)
    zabc(2,2)=zabc(1,1);
end

%Calculating the neutral transformation matrix 
tn = -(znn^-1*znj);

% fprintf('The distance matrix is\n')
% fprintf('\n')
% disp('D = ') 
% fprintf('\n')
% disp(D)
% 
% fprintf('The primitive impedance matrix in ohms/mile is\n')
% fprintf('\n')
% disp('[z] = ') 
% fprintf('\n')
% disp(zprim)
% 
% fprintf('In partioned form, the submatrices in ohms/mile are\n')
% fprintf('\n')
% disp('[zij] = ') 
% fprintf('\n')
% disp(zij)
% disp('[zin] = ') 
% fprintf('\n')
% disp(zin)
% disp('[znj] = ') 
% fprintf('\n')
% disp(znj)
% disp('[znn] = ') 
% fprintf('\n')
% disp(znn)

fprintf('The "Kron" reduced phase impedance matrix in ohms/mile is\n ')
fprintf('\n')
disp('[zabc] = ') 
fprintf('\n')
disp(zabc)

% fprintf('The neutral transformation matrix is\n')
% fprintf('\n')
% disp('[tn] = ') 
% fprintf('\n')
% disp(tn)