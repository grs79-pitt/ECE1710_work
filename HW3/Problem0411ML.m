%Distribution System Modelling and Analysis, Example 4.1
%Written by William Kersting and Robert Kerestes

%Modified for Problem 4.2
clear all
clc
j = sqrt(-1);

%Defining the phase conductor data
phase.GMR = 0.0171; %Modified for 1/0 ASCR
phase.resistance = .41; %Modified for 1/0 ASCR
phase.ncond = 3; %Three phase

%Defining the neutral conductor data
neutral.GMR = 0.00208; %Modified for 1/0 ASCR
neutral.resistance = 14.8722; %Modified for 1/0 ASCR
neutral.ncond = 3; %1 neutral
ds=.0641; %inches
dod=1.29; %inches
ncond=phase.ncond+neutral.ncond;


k=13;
R=(dod-ds)/24;
GMR_cn=(neutral.GMR*k*(R^(k-1)))^(1/k);
r_cn=neutral.resistance/k;

%spacings
d1=0.5+j*0; %Location of Phase A Conductor
d2=1+j*0; %Location of Phase B Conductor
d3=0+j*0; %Location of Phase C Conductor

d4=0.5+j*R; %Location of Phase A Concentric Neutral
d5=1+j*R; %Location of Phase B Concentric Neutral
d6=0+j*R; %Location of Phase C Neutral

d = [d1;d2;d3;d4;d5;d6]; %position vector

ncond = phase.ncond+neutral.ncond; %number of total conductors

%Initializing matrix sizes
r = zeros(ncond,1);
D = zeros(ncond,ncond);
zprim = zeros(ncond,ncond);

%Defining the resistance vector r
for i = 1:1:ncond
    
    if i <= ncond -1 
        
        r(i) = phase.resistance;
    
    else
        
        r(i) = r_cn;
        
    end 
        
end

%Calculating the inpedance matrix D
for i = 1:1:ncond
    
    for k = 1:1:ncond
        
        if i == k && i <= ncond/2
            
            D(i,k) = phase.GMR;
            
        elseif i == k && i > ncond/2
            
            D(i,k) = GMR_cn;
            
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
zabc = zij-zin*znn^-1*znj;
zabc_1kFt=zabc*(1000/5280);

%Calculating the neutral transformation matrix 
 tn = -(znn^-1*znj);
% 
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

fprintf('The "Kron" reduced phase impedance matrix in ohms/1k ft is\n ')
fprintf('\n')
disp('[zabc] = ') 
fprintf('\n')
disp(zabc_1kFt)

% fprintf('The neutral transformation matrix is\n')
% fprintf('\n')
% disp('[tn] = ') 
% fprintf('\n')
% disp(tn)