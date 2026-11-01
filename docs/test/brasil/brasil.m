function [p z] = brasil(beta,m)

com = sprintf("/Applications/Anaconda-Navigator.app/Contents/MacOS/python ratapprox.py %f %d", beta, m);
system(com);
p =  csvread(strcat('poles',num2str(m),'.csv'));
z = csvread(strcat('zeros',num2str(m),'.csv'));