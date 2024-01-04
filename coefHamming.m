function w=coefHamming(N,a)
    alfa=0:2*pi/N:2*pi-2*pi/N;
    w=transpose(a-(1-a)*cos(alfa));
end
